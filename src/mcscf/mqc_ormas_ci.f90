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
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_determinants, only: link_table_t, string_address, build_link_table
   use mqc_ci, only: sigma_vector, absorb_one_electron
   use mqc_davidson, only: sigma_operator_t, davidson_flat
   use pic_blas_interfaces, only: pic_gemm
   use mqc_ormas_space, only: ormas_space_t, ormas_closure_t, closure_address, &
                              build_ormas_closure, ormas_strings
   implicit none
   private

   public :: ormas_diagonal
   public :: full_space_index
   public :: ormas_sigma
   public :: ormas_lowest
   public :: ormas_sigma_direct
   public :: ormas_solve

   type, extends(sigma_operator_t) :: ormas_operator_t
      !! The restricted Hamiltonian, as something Davidson can multiply by
      !!
      !! Holds what `ormas_sigma_direct` needs and nothing else. The vector it
      !! is handed is already flat, which is the reason the method had to be
      !! written that way: a restricted space has no rectangle to reshape into.
      type(ormas_space_t), pointer :: space => null()
      type(ormas_closure_t), pointer :: closure => null()
      real(dp), pointer :: folded(:, :) => null()
      type(link_table_t), pointer :: alpha_table => null(), beta_table => null()
      integer, allocatable :: in_alpha(:), in_beta(:), from_alpha(:), from_beta(:)
   contains
      procedure :: apply => ormas_apply
   end type ormas_operator_t

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

   subroutine full_space_index(space, alpha, beta, in_alpha, in_beta, from_alpha, from_beta)
      !! Where each of the space's strings sits in the unrestricted list
      !!
      !! The restricted space groups its strings by occupation class; the
      !! complete space orders them by bit pattern. Both index the same objects,
      !! so one permutation relates them, and it is wanted because the machinery
      !! being borrowed below -- excitation tables, the sigma build -- is
      !! written against the complete space.
      type(ormas_space_t), intent(in) :: space
      integer(int64), intent(in) :: alpha(:), beta(:)
      integer, allocatable, intent(out) :: in_alpha(:), in_beta(:)
      integer, allocatable, intent(out), optional :: from_alpha(:), from_beta(:)

      integer :: i

      allocate (in_alpha(size(alpha)), in_beta(size(beta)))
      do i = 1, size(alpha)
         in_alpha(i) = string_address(space%n_active_orbitals, space%n_alpha, alpha(i))
      end do
      do i = 1, size(beta)
         in_beta(i) = string_address(space%n_active_orbitals, space%n_beta, beta(i))
      end do

      ! Every string of the complete space is kept, so the map is onto and its
      ! inverse is a permutation rather than a partial one.
      if (present(from_alpha)) then
         allocate (from_alpha(size(alpha)))
         do i = 1, size(alpha)
            from_alpha(in_alpha(i)) = i
         end do
      end if
      if (present(from_beta)) then
         allocate (from_beta(size(beta)))
         do i = 1, size(beta)
            from_beta(in_beta(i)) = i
         end do
      end if
   end subroutine full_space_index

   subroutine ormas_sigma(space, folded, alpha_table, beta_table, in_alpha, in_beta, &
                          ci, sigma, error)
      !! `sigma = P H P c`, with `P` the projection onto the restricted space
      !!
      !! Spread the vector over the complete space, putting zero everywhere the
      !! restriction excludes; apply the complete-space Hamiltonian; keep only
      !! what lands back inside. Because the vector vanishes outside, the sum
      !! that reaches an in-space determinant runs only over in-space
      !! determinants, so this *is* the restricted Hamiltonian -- not an
      !! approximation to it.
      !!
      !! Which makes it correct by construction and, for the moment, the whole
      !! implementation: it inherits a sigma already checked against PySCF
      !! rather than introducing a second one to get wrong.
      !!
      !! What it does not inherit is the saving. The work is that of the
      !! complete space, so this is the route while a space is small enough to
      !! check, and not the route that makes a restricted space worth having.
      !! The intermediate a direct build needs is subtler than it looks: an
      !! excitation may leave the space and a second bring it back, so the
      !! basis it lives on is the space plus one excitation, and truncating it
      !! to the space itself silently drops terms. That is the next piece, and
      !! it has to reproduce this one.
      type(ormas_space_t), intent(in) :: space
      real(dp), intent(in) :: folded(:, :)
      type(link_table_t), intent(in) :: alpha_table, beta_table
      integer, intent(in) :: in_alpha(:), in_beta(:)
      real(dp), intent(in) :: ci(:)                  !! (n_determinants)
      real(dp), intent(out) :: sigma(:)              !! (n_determinants)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: embedded(:, :), acted(:, :)
      integer :: ia, ib, ga, gb
      integer(int64) :: at

      if (error%has_error()) return

      allocate (embedded(alpha_table%n_strings, beta_table%n_strings))
      allocate (acted(alpha_table%n_strings, beta_table%n_strings))
      embedded = 0.0_dp

      do ia = 1, size(in_alpha)
         ga = space%alpha_string_class(ia)
         do ib = 1, size(in_beta)
            gb = space%beta_string_class(ib)
            if (.not. space%compatible(gb, ga)) cycle
            at = determinant_address_local(space, ia, ib)
            embedded(in_alpha(ia), in_beta(ib)) = ci(at)
         end do
      end do

      call sigma_vector(folded, embedded, alpha_table, beta_table, acted, error)
      if (error%has_error()) return

      do ia = 1, size(in_alpha)
         ga = space%alpha_string_class(ia)
         do ib = 1, size(in_beta)
            gb = space%beta_string_class(ib)
            if (.not. space%compatible(gb, ga)) cycle
            at = determinant_address_local(space, ia, ib)
            sigma(at) = acted(in_alpha(ia), in_beta(ib))
         end do
      end do

      deallocate (embedded, acted)
   end subroutine ormas_sigma

   pure function determinant_address_local(space, ia, ib) result(at)
      !! `determinant_address`, spelled out here so this module needs no
      !! circular use of the one it is built on
      type(ormas_space_t), intent(in) :: space
      integer, intent(in) :: ia, ib
      integer(int64) :: at

      integer :: ga, gb

      ga = space%alpha_string_class(ia)
      gb = space%beta_string_class(ib)
      at = space%alpha_base(ia) + space%block_offset(gb, ga) &
           + (int(ib, int64) - space%beta_offset(gb))
   end function determinant_address_local

   subroutine ormas_lowest(space, folded, alpha_table, beta_table, in_alpha, in_beta, &
                           n_roots, energies, error)
      !! The lowest eigenvalues of the restricted Hamiltonian, densely
      !!
      !! Builds the matrix a column at a time by applying `ormas_sigma` to each
      !! unit vector, then diagonalises it. That costs a determinant's worth of
      !! sigma per determinant and is only sane while the space is small -- but
      !! it introduces no matrix elements of its own, so the energies it
      !! returns are exactly what the sigma says they are and nothing has been
      !! written twice.
      type(ormas_space_t), intent(in) :: space
      real(dp), intent(in) :: folded(:, :)
      type(link_table_t), intent(in) :: alpha_table, beta_table
      integer, intent(in) :: in_alpha(:), in_beta(:)
      integer, intent(in) :: n_roots
      real(dp), allocatable, intent(out) :: energies(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: matrix(:, :), unit(:), column(:), values(:)
      integer :: ndet, d, info

      if (error%has_error()) return

      ndet = int(space%n_determinants)
      if (n_roots < 1 .or. n_roots > ndet) then
         call error%set(ERROR_VALIDATION, "asked for "//to_char(n_roots)// &
                        " roots of a space holding "//to_char(ndet)//" determinants.")
         return
      end if

      allocate (matrix(ndet, ndet), unit(ndet), column(ndet), values(ndet))
      do d = 1, ndet
         unit = 0.0_dp
         unit(d) = 1.0_dp
         call ormas_sigma(space, folded, alpha_table, beta_table, in_alpha, in_beta, &
                          unit, column, error)
         if (error%has_error()) return
         matrix(:, d) = column
      end do

      call pic_syev(matrix, values, jobz="N", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "diagonalising the restricted Hamiltonian "// &
                        "failed with LAPACK info "//to_char(info)//".")
         return
      end if

      allocate (energies(n_roots))
      energies = values(1:n_roots)
      deallocate (matrix, unit, column, values)
   end subroutine ormas_lowest

   subroutine ormas_sigma_direct(space, closure, folded, alpha_table, beta_table, &
                                 in_alpha, in_beta, from_alpha, from_beta, ci, sigma, error)
      !! `sigma = H c` without ever visiting the complete space
      !!
      !! The same three passes as the unrestricted build -- apply an excitation
      !! to the vector for every orbital pair at once, contract that against the
      !! Hamiltonian, apply an excitation again -- with the intermediate living
      !! on the closure rather than on the whole rectangle.
      !!
      !! That is the only real difference and it is the whole point. The
      !! intermediate has to be wider than the space, because an excitation can
      !! leave and a second bring it back; it does not have to be as wide as the
      !! complete space, because two excitations from the space is further than
      !! anything here reaches. One excitation is exactly enough, and for a
      !! restricted space that is a far smaller set than the rectangle the
      !! projected build pays for.
      !!
      !! The excitation tables are the unrestricted ones, indexed by bit
      !! pattern; the strings here are grouped by occupation class. `in_*` and
      !! `from_*` carry an index between the two orderings, in each direction.
      type(ormas_space_t), intent(in) :: space
      type(ormas_closure_t), intent(in) :: closure
      real(dp), intent(in) :: folded(:, :)
      type(link_table_t), intent(in) :: alpha_table, beta_table
      integer, intent(in) :: in_alpha(:), in_beta(:)      !! class order -> bit order
      integer, intent(in) :: from_alpha(:), from_beta(:)  !! bit order -> class order
      real(dp), intent(in) :: ci(:)
      real(dp), intent(out) :: sigma(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: gathered(:, :), contracted(:, :)
      integer :: norb, npair, ia, ib, ia2, ib2, ga, gb, row, pair
      real(dp) :: weight

      if (error%has_error()) return

      norb = space%n_active_orbitals
      npair = norb*norb

      allocate (gathered(npair, closure%n_determinants))
      gathered = 0.0_dp

      ! Gather: read the vector where it lives, write the intermediate where the
      ! excitation lands.
      do ia = 1, size(in_alpha)
         ga = space%alpha_string_class(ia)
         do row = 1, alpha_table%n_rows
            pair = alpha_table%cre(row, in_alpha(ia)) &
                   + (alpha_table%des(row, in_alpha(ia)) - 1)*norb
            ia2 = from_alpha(alpha_table%dest(row, in_alpha(ia)))
            weight = real(alpha_table%phase(row, in_alpha(ia)), dp)
            do ib = 1, size(in_beta)
               gb = space%beta_string_class(ib)
               if (.not. space%compatible(gb, ga)) cycle
               gathered(pair, closure_address(closure, space, ia2, ib)) = &
                  gathered(pair, closure_address(closure, space, ia2, ib)) &
                  + weight*ci(determinant_address_local(space, ia, ib))
            end do
         end do
      end do

      do ib = 1, size(in_beta)
         gb = space%beta_string_class(ib)
         do row = 1, beta_table%n_rows
            pair = beta_table%cre(row, in_beta(ib)) &
                   + (beta_table%des(row, in_beta(ib)) - 1)*norb
            ib2 = from_beta(beta_table%dest(row, in_beta(ib)))
            weight = real(beta_table%phase(row, in_beta(ib)), dp)
            do ia = 1, size(in_alpha)
               ga = space%alpha_string_class(ia)
               if (.not. space%compatible(gb, ga)) cycle
               gathered(pair, closure_address(closure, space, ia, ib2)) = &
                  gathered(pair, closure_address(closure, space, ia, ib2)) &
                  + weight*ci(determinant_address_local(space, ia, ib))
            end do
         end do
      end do

      allocate (contracted(npair, closure%n_determinants))
      call pic_gemm(folded, gathered, contracted)

      ! Scatter: read the intermediate where it lives, write the sigma where the
      ! excitation lands -- keeping only what lands back inside the space.
      sigma = 0.0_dp
      do ia = 1, size(in_alpha)
         ga = space%alpha_string_class(ia)
         do row = 1, alpha_table%n_rows
            pair = alpha_table%cre(row, in_alpha(ia)) &
                   + (alpha_table%des(row, in_alpha(ia)) - 1)*norb
            ia2 = from_alpha(alpha_table%dest(row, in_alpha(ia)))
            weight = real(alpha_table%phase(row, in_alpha(ia)), dp)
            do ib = 1, size(in_beta)
               gb = space%beta_string_class(ib)
               if (.not. closure%present(gb, ga)) cycle
               if (.not. space%compatible(gb, space%alpha_string_class(ia2))) cycle
               sigma(determinant_address_local(space, ia2, ib)) = &
                  sigma(determinant_address_local(space, ia2, ib)) &
                  + weight*contracted(pair, closure_address(closure, space, ia, ib))
            end do
         end do
      end do

      do ib = 1, size(in_beta)
         gb = space%beta_string_class(ib)
         do row = 1, beta_table%n_rows
            pair = beta_table%cre(row, in_beta(ib)) &
                   + (beta_table%des(row, in_beta(ib)) - 1)*norb
            ib2 = from_beta(beta_table%dest(row, in_beta(ib)))
            weight = real(beta_table%phase(row, in_beta(ib)), dp)
            do ia = 1, size(in_alpha)
               ga = space%alpha_string_class(ia)
               if (.not. closure%present(gb, ga)) cycle
               if (.not. space%compatible(space%beta_string_class(ib2), ga)) cycle
               sigma(determinant_address_local(space, ia, ib2)) = &
                  sigma(determinant_address_local(space, ia, ib2)) &
                  + weight*contracted(pair, closure_address(closure, space, ia, ib))
            end do
         end do
      end do

      deallocate (gathered, contracted)
   end subroutine ormas_sigma_direct

   subroutine ormas_apply(this, vector, image, error)
      !! One matrix-vector product over the restricted space
      class(ormas_operator_t), intent(inout) :: this
      real(dp), intent(in) :: vector(:)
      real(dp), intent(out) :: image(:)
      type(error_t), intent(inout) :: error

      call ormas_sigma_direct(this%space, this%closure, this%folded, this%alpha_table, &
                              this%beta_table, this%in_alpha, this%in_beta, &
                              this%from_alpha, this%from_beta, vector, image, error)
   end subroutine ormas_apply

   subroutine ormas_solve(space, h1e, eri, n_roots, energies, vectors, error, &
                          tolerance, max_iterations)
      !! The lowest states of a restricted active space, iteratively
      !!
      !! Everything from the partition and the integrals: the strings, the
      !! closure the sigma needs, the excitation tables, the folded Hamiltonian,
      !! the diagonal to precondition with, and then Davidson.
      !!
      !! This is the route that scales. The dense solve alongside it costs a
      !! sigma per determinant and exists to be compared against; this one costs
      !! a sigma per expansion vector, which is tens rather than thousands, and
      !! never forms a matrix.
      type(ormas_space_t), intent(in), target :: space
      real(dp), intent(in) :: h1e(:, :), eri(:, :, :, :)
      integer, intent(in) :: n_roots
      real(dp), allocatable, intent(out) :: energies(:)
      real(dp), allocatable, intent(out) :: vectors(:, :)   !! (n_determinants, n_roots)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations

      type(ormas_operator_t) :: operator
      type(ormas_closure_t), target :: closure
      type(link_table_t), target :: alpha_table, beta_table
      real(dp), allocatable, target :: folded(:, :)
      real(dp), allocatable :: diagonal(:), residuals(:)
      integer(int64), allocatable :: alpha(:), beta(:)
      integer :: iterations_taken, sigma_products
      logical :: converged

      if (error%has_error()) return

      call ormas_strings(space, alpha, beta, error)
      call full_space_index(space, alpha, beta, operator%in_alpha, operator%in_beta, &
                            operator%from_alpha, operator%from_beta)
      call build_ormas_closure(space, closure, error)
      call ormas_diagonal(space, alpha, beta, h1e, eri, diagonal, error)
      call absorb_one_electron(h1e, eri, space%n_alpha + space%n_beta, folded, error)
      call build_link_table(space%n_active_orbitals, space%n_alpha, alpha_table, error)
      call build_link_table(space%n_active_orbitals, space%n_beta, beta_table, error)
      if (error%has_error()) return

      operator%space => space
      operator%closure => closure
      operator%folded => folded
      operator%alpha_table => alpha_table
      operator%beta_table => beta_table

      call davidson_flat(operator, diagonal, n_roots, energies, vectors, residuals, &
                         iterations_taken, sigma_products, converged, error, &
                         tolerance, max_iterations)
      if (error%has_error()) return

      if (.not. converged) then
         call error%set(ERROR_VALIDATION, "the restricted CI did not converge in "// &
                        to_char(iterations_taken)//" iterations; the largest residual "// &
                        "was "//to_char(maxval(residuals))//".")
      end if

      call closure%destroy()
      call alpha_table%destroy()
      call beta_table%destroy()
   end subroutine ormas_solve

end module mqc_ormas_ci
