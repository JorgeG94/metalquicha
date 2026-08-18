!! Configuration interaction over an occupation-restricted active space
module mqc_ormas_ci
   !! What `mqc_ci` does for a complete active space, for a restricted one.
   !!
   !! The physics is untouched: ORMAS restricts which determinants exist, not
   !! what the Hamiltonian is, so every matrix element here is the same Slater
   !! rule it would be in a CAS. What changes is the bookkeeping around it --
   !! the determinants no longer form a rectangle, so a quantity indexed by one
   !! is a flat vector rather than a matrix, and reaching a determinant means
   !! going through `mqc_ormas_space`.
   !!
   !! The diagonal comes first and on its own. It is the whole Davidson
   !! preconditioner, it needs no excitations, and -- being written by walking
   !! the determinants in storage order while an independent path reaches the
   !! same elements through the address formula -- it checks the addressing
   !! exactly, before anything harder is built on top of it.
   use pic_types, only: dp, int64
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_ormas_space, only: ormas_space_t
   implicit none
   private

   public :: ormas_diagonal

contains

   pure subroutine occupied_orbitals(string, n_orbitals, orbitals, n_found)
      !! The orbitals a string occupies, ascending
      integer(int64), intent(in) :: string
      integer, intent(in) :: n_orbitals
      integer, intent(out) :: orbitals(:)
      integer, intent(out) :: n_found

      integer :: p

      n_found = 0
      do p = 1, n_orbitals
         if (btest(string, p - 1)) then
            n_found = n_found + 1
            orbitals(n_found) = p
         end if
      end do
   end subroutine occupied_orbitals

   subroutine string_terms(strings, n_orbitals, n_electrons, h1e, coulomb, exchange, &
                           one_body, same_spin, screened, occupied)
      !! Everything about a diagonal element that one spin decides alone
      !!
      !! Three things per string: its one-electron sum, its own Coulomb and
      !! exchange, and -- the one worth hoisting -- the vector
      !! `screened(q) = sum over occupied p of (pp|qq)`. Without it the
      !! opposite-spin Coulomb term costs a double loop for every determinant;
      !! with it, a single sum over the other spin's electrons.
      integer(int64), intent(in) :: strings(:)
      integer, intent(in) :: n_orbitals, n_electrons
      real(dp), intent(in) :: h1e(:, :), coulomb(:, :), exchange(:, :)
      real(dp), intent(out) :: one_body(:), same_spin(:)
      real(dp), intent(out) :: screened(:, :)
      integer, intent(out) :: occupied(:, :)

      integer :: s, ip, iq, p, q, n_found

      do s = 1, size(strings)
         call occupied_orbitals(strings(s), n_orbitals, occupied(:, s), n_found)

         one_body(s) = 0.0_dp
         same_spin(s) = 0.0_dp
         screened(:, s) = 0.0_dp

         do ip = 1, n_electrons
            p = occupied(ip, s)
            one_body(s) = one_body(s) + h1e(p, p)
            screened(:, s) = screened(:, s) + coulomb(p, :)
            ! Pairs once each, which is where the usual half goes.
            do iq = ip + 1, n_electrons
               q = occupied(iq, s)
               same_spin(s) = same_spin(s) + coulomb(p, q) - exchange(p, q)
            end do
         end do
      end do
   end subroutine string_terms

   subroutine ormas_diagonal(space, alpha, beta, h1e, eri, diagonal, error)
      !! `<D|H|D>` for every determinant of the space
      !!
      !! For a determinant occupying orbitals `A` with alpha spin and `B` with
      !! beta,
      !!
      !!     sum_{p in A} h_pp + sum_{p in B} h_pp
      !!       + sum_{p<q in A} [(pp|qq) - (pq|qp)]
      !!       + sum_{p<q in B} [(pp|qq) - (pq|qp)]
      !!       + sum_{p in A} sum_{q in B} (pp|qq)
      !!
      !! Coulomb and exchange within a spin, Coulomb only between the two,
      !! because electrons of opposite spin do not exchange. Identical to the
      !! complete-active-space expression -- a restricted space contains fewer
      !! determinants, not different ones.
      !!
      !! Filled by walking the space in storage order: alpha class, then alpha
      !! string, then each compatible beta class, then its strings. That order
      !! is exactly the one the offsets in `mqc_ormas_space` were accumulated
      !! in, so a plain counter tracks the address and the writes run straight
      !! down the array. It also means a disagreement between this counter and
      !! the address formula is a disagreement about the layout itself, which
      !! is why the test asserts they agree at every determinant rather than
      !! only that the energies come out right.
      !!
      !! Takes the raw integrals, not the folded tensor: folding rearranges the
      !! Hamiltonian for `E_pq E_rs` and its diagonal is not the diagonal of
      !! anything physical.
      type(ormas_space_t), intent(in) :: space
      integer(int64), intent(in) :: alpha(:), beta(:)   !! As `ormas_strings` gives them
      real(dp), intent(in) :: h1e(:, :)                 !! (norb, norb), symmetric
      real(dp), intent(in) :: eri(:, :, :, :)           !! chemist (pq|rs)
      real(dp), allocatable, intent(out) :: diagonal(:)  !! (n_determinants)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: coulomb(:, :), exchange(:, :)
      real(dp), allocatable :: one_a(:), same_a(:), screen_a(:, :)
      real(dp), allocatable :: one_b(:), same_b(:), screen_b(:, :)
      integer, allocatable :: occ_a(:, :), occ_b(:, :)
      real(dp) :: cross, alpha_part
      integer :: norb, na, nb, p, q, ga, gb, iq
      integer(int64) :: a, b, at

      if (error%has_error()) return

      norb = space%n_active_orbitals
      if (size(h1e, 1) /= norb .or. size(h1e, 2) /= norb .or. any(shape(eri) /= norb)) then
         call error%set(ERROR_VALIDATION, "the integrals are over "// &
                        to_char(size(h1e, 1))//" orbitals and the partition covers "// &
                        to_char(norb)//".")
         return
      end if

      na = size(alpha)
      nb = size(beta)

      ! The only integrals a diagonal element ever sees.
      allocate (coulomb(norb, norb), exchange(norb, norb))
      do q = 1, norb
         do p = 1, norb
            coulomb(p, q) = eri(p, p, q, q)
            exchange(p, q) = eri(p, q, q, p)
         end do
      end do

      allocate (one_a(na), same_a(na), screen_a(norb, na))
      allocate (one_b(nb), same_b(nb), screen_b(norb, nb))
      allocate (occ_a(max(space%n_alpha, 1), na), occ_b(max(space%n_beta, 1), nb))

      call string_terms(alpha, norb, space%n_alpha, h1e, coulomb, exchange, &
                        one_a, same_a, screen_a, occ_a)
      call string_terms(beta, norb, space%n_beta, h1e, coulomb, exchange, &
                        one_b, same_b, screen_b, occ_b)

      allocate (diagonal(space%n_determinants))

      at = 0_int64
      do ga = 1, space%n_alpha_classes
         do a = space%alpha_offset(ga) + 1_int64, space%alpha_offset(ga + 1)
            alpha_part = one_a(a) + same_a(a)
            do gb = 1, space%n_beta_classes
               if (.not. space%compatible(gb, ga)) cycle
               do b = space%beta_offset(gb) + 1_int64, space%beta_offset(gb + 1)
                  cross = 0.0_dp
                  do iq = 1, space%n_beta
                     cross = cross + screen_a(occ_b(iq, b), a)
                  end do
                  at = at + 1_int64
                  diagonal(at) = alpha_part + one_b(b) + same_b(b) + cross
               end do
            end do
         end do
      end do

      deallocate (coulomb, exchange, one_a, same_a, screen_a, one_b, same_b, screen_b)
      deallocate (occ_a, occ_b)
   end subroutine ormas_diagonal

end module mqc_ormas_ci
