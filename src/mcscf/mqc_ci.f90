!! The CI Hamiltonian, applied to a vector without ever being formed
module mqc_ci
   !! `sigma = H c` for a complete active space, which is the operation a CI
   !! solver spends all of its time in and the only one that has to be fast.
   !! The Hamiltonian is never built: at CAS(10,10) it would be 32 GB and at
   !! CAS(12,12) it does not bear writing down, while the vector it acts on is
   !! 508 kB and 6.8 MB.
   !!
   !! **One contraction, not three.** The textbook sigma splits into an
   !! alpha-alpha, a beta-beta and an alpha-beta term, each with its own loop
   !! structure and its own sign conventions to get wrong. That is avoidable.
   !! The two-particle part of the Hamiltonian is
   !!
   !!     (pq|rs) E_pr,qs = (pq|rs) (E_pq E_rs - delta_qr E_ps)
   !!
   !! so subtracting the `delta_qr` term from the one-electron Hamiltonian and
   !! then folding *that* into the two-electron tensor -- `absorb_one_electron`
   !! below -- leaves the entire Hamiltonian as a single `E_pq E_rs`. Applying
   !! it is then: apply `E_pq` once, contract with the tensor, apply `E_pq`
   !! again. No case analysis, and the same excitation table for both spins.
   !!
   !! The idea is PySCF's, from `pyscf/fci/direct_spin1.py` and the readable
   !! `fci/fci_slow.py` beside it. The implementation here is its own; what was
   !! taken is the formulation.
   !!
   !! The cost is the intermediate, `n_orbitals^2` by the number of determinants:
   !! 51 MB at CAS(10,10), 983 MB at CAS(12,12), and past that it has to be
   !! blocked over beta strings. Blocking changes nothing structurally and is
   !! deliberately not done yet.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_determinants, only: link_table_t
   implicit none
   private

   public :: absorb_one_electron
   public :: sigma_vector
   public :: ci_hamiltonian
   public :: ci_diagonal
   public :: apply_excitations

contains

   subroutine absorb_one_electron(h1e, eri, n_electrons, folded, error)
      !! Fold the one-electron Hamiltonian into the two-electron tensor
      !!
      !! Returns a tensor `g` such that the whole Hamiltonian is
      !! `sum_pqrs g_pq,rs E_pq E_rs`, so that `sigma_vector` has one term to
      !! evaluate instead of a one-electron term and a two-electron term with
      !! different structures.
      !!
      !! Two steps. First the `delta_qr` correction from rewriting `E_pr,qs`,
      !! which is a one-electron operator and joins `h1e`:
      !!
      !!     f_jk = h_jk - (1/2) sum_i (ji|ik)
      !!
      !! Then `f` is spread over the diagonal of the two-electron tensor, using
      !! `sum_k E_kk = N`: adding `f/N` to `g_kk,rs` and to `g_pq,kk` for every
      !! `k` reproduces `sum f_pq E_pq` exactly, because the number operator
      !! commutes with everything that conserves particles and each of the two
      !! placements contributes half.
      !!
      !! `factor` is folded in at the end and is normally a half, which is where
      !! the `1/2` in front of the two-electron energy goes. `sigma_vector` has
      !! no factors of its own.
      real(dp), intent(in) :: h1e(:, :)          !! (norb, norb), symmetric
      real(dp), intent(in) :: eri(:, :, :, :)    !! (norb^4), chemist (pq|rs)
      integer, intent(in) :: n_electrons         !! Total, both spins
      real(dp), allocatable, intent(out) :: folded(:, :)
         !! (norb^2, norb^2), indexed `(p + (q-1)norb, r + (s-1)norb)`. Already
         !! multiplied by a half.
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: f1e(:, :)
      real(dp) :: correction
      integer :: norb, p, q, r, s, i, k, pq, rs, kk

      if (error%has_error()) return
      norb = size(h1e, 1)
      if (size(h1e, 2) /= norb .or. any(shape(eri) /= norb)) then
         call error%set(ERROR_VALIDATION, "the one- and two-electron integrals "// &
                        "disagree about how many active orbitals there are: h1e is "// &
                        to_char(size(h1e, 1))//" by "//to_char(size(h1e, 2))// &
                        " and the two-electron tensor is "//to_char(size(eri, 1))//".")
         return
      end if
      if (n_electrons <= 0) then
         call error%set(ERROR_VALIDATION, "an active space with "// &
                        to_char(n_electrons)//" electrons has no one-electron "// &
                        "Hamiltonian to fold: the spreading below divides by the "// &
                        "electron count.")
         return
      end if

      allocate (f1e(norb, norb))
      do k = 1, norb
         do p = 1, norb
            correction = 0.0_dp
            do i = 1, norb
               correction = correction + eri(p, i, i, k)
            end do
            f1e(p, k) = (h1e(p, k) - 0.5_dp*correction)/real(n_electrons, dp)
         end do
      end do

      allocate (folded(norb*norb, norb*norb))
      do s = 1, norb
         do r = 1, norb
            rs = r + (s - 1)*norb
            do q = 1, norb
               do p = 1, norb
                  folded(p + (q - 1)*norb, rs) = eri(p, q, r, s)
               end do
            end do
         end do
      end do

      do k = 1, norb
         kk = k + (k - 1)*norb
         do s = 1, norb
            do r = 1, norb
               rs = r + (s - 1)*norb
               folded(kk, rs) = folded(kk, rs) + f1e(r, s)
            end do
         end do
         do q = 1, norb
            do p = 1, norb
               pq = p + (q - 1)*norb
               folded(pq, kk) = folded(pq, kk) + f1e(p, q)
            end do
         end do
      end do

      folded = 0.5_dp*folded
      deallocate (f1e)
   end subroutine absorb_one_electron

   subroutine apply_excitations(ci, alpha, beta, gathered)
      !! `E_pq |c>` for every orbital pair at once
      !!
      !! Returns `(n_orbitals^2, n_determinants)`, column `d` and row
      !! `p + (q-1) n_orb` holding the determinant-`d` component of `E_pq |c>`.
      !!
      !! Shared between the sigma build and the density matrices, which is not
      !! an economy but a correctness argument: both are built entirely out of
      !! this one operation, and it carries the excitation phases. Two copies
      !! would be two places for a sign to go wrong, and the second copy would
      !! be exercised only by the density matrices -- where a sign error does
      !! not change the energy at all and so would survive every test that
      !! looks at one.
      !!
      !! Alpha excitations move the first index of the vector and beta the
      !! second. That is the only place the two spins differ anywhere in the CI.
      real(dp), intent(in) :: ci(:, :)      !! (n_alpha_strings, n_beta_strings)
      type(link_table_t), intent(in) :: alpha, beta
      real(dp), allocatable, intent(out) :: gathered(:, :)

      integer :: norb, na, nb, npair, ndet, source, target_string, row, pair, ia, ib
      real(dp) :: weight

      norb = alpha%n_orbitals
      na = alpha%n_strings
      nb = beta%n_strings
      npair = norb*norb
      ndet = na*nb

      allocate (gathered(npair, ndet))
      gathered = 0.0_dp

      do source = 1, na
         do row = 1, alpha%n_rows
            pair = alpha%cre(row, source) + (alpha%des(row, source) - 1)*norb
            target_string = alpha%dest(row, source)
            weight = real(alpha%phase(row, source), dp)
            do ib = 1, nb
               gathered(pair, target_string + (ib - 1)*na) = &
                  gathered(pair, target_string + (ib - 1)*na) + weight*ci(source, ib)
            end do
         end do
      end do

      do source = 1, nb
         do row = 1, beta%n_rows
            pair = beta%cre(row, source) + (beta%des(row, source) - 1)*norb
            target_string = beta%dest(row, source)
            weight = real(beta%phase(row, source), dp)
            do ia = 1, na
               gathered(pair, ia + (target_string - 1)*na) = &
                  gathered(pair, ia + (target_string - 1)*na) + weight*ci(ia, source)
            end do
         end do
      end do
   end subroutine apply_excitations

   subroutine sigma_vector(folded, ci, alpha, beta, sigma, error)
      !! `sigma = H c`, with `H` the folded tensor acting as `E_pq E_rs`
      !!
      !! Three passes. Apply `E_pq` to the vector for every orbital pair at
      !! once, gathering into an intermediate indexed by pair and determinant;
      !! contract that against the tensor, which is one matrix multiply and
      !! where all the arithmetic is; then apply `E_pq` again, scattering back.
      !!
      !! The gather and the scatter have the same loop, because `E_pq` and its
      !! adjoint are the same table read in opposite directions -- the gather
      !! writes to the destination string and reads the source, the scatter
      !! writes to the destination and reads the source of the *second*
      !! application. Alpha excitations move the first index of the vector and
      !! beta excitations the second, which is the only place the two spins
      !! differ.
      real(dp), intent(in) :: folded(:, :)   !! From `absorb_one_electron`
      real(dp), intent(in) :: ci(:, :)       !! (n_alpha_strings, n_beta_strings)
      type(link_table_t), intent(in) :: alpha, beta
      real(dp), intent(out) :: sigma(:, :)   !! Same shape as `ci`
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: gathered(:, :), contracted(:, :)
      integer :: norb, na, nb, npair, ndet
      integer :: source, target_string, row, pair, ia, ib
      real(dp) :: weight

      if (error%has_error()) return
      norb = alpha%n_orbitals
      na = alpha%n_strings
      nb = beta%n_strings
      npair = norb*norb
      ndet = na*nb

      if (beta%n_orbitals /= norb) then
         call error%set(ERROR_VALIDATION, "the alpha and beta excitation tables "// &
                        "describe active spaces of "//to_char(norb)//" and "// &
                        to_char(beta%n_orbitals)//" orbitals. Both spins occupy the "// &
                        "same orbitals; only the electron counts may differ.")
         return
      end if
      if (size(ci, 1) /= na .or. size(ci, 2) /= nb) then
         call error%set(ERROR_VALIDATION, "the vector is "//to_char(size(ci, 1))// &
                        " by "//to_char(size(ci, 2))//" but the excitation tables "// &
                        "have "//to_char(na)//" alpha and "//to_char(nb)// &
                        " beta strings.")
         return
      end if
      if (size(folded, 1) /= npair .or. size(folded, 2) /= npair) then
         call error%set(ERROR_VALIDATION, "the folded Hamiltonian is "// &
                        to_char(size(folded, 1))//" square but "//to_char(norb)// &
                        " active orbitals need "//to_char(npair)//".")
         return
      end if

      call apply_excitations(ci, alpha, beta, gathered)

      ! The whole two-electron contraction, as one matrix multiply.
      allocate (contracted(npair, ndet))
      call pic_gemm(folded, gathered, contracted)

      sigma = 0.0_dp
      do source = 1, na
         do row = 1, alpha%n_rows
            pair = alpha%cre(row, source) + (alpha%des(row, source) - 1)*norb
            target_string = alpha%dest(row, source)
            weight = real(alpha%phase(row, source), dp)
            do ib = 1, nb
               sigma(target_string, ib) = sigma(target_string, ib) &
                                          + weight*contracted(pair, source + (ib - 1)*na)
            end do
         end do
      end do

      do source = 1, nb
         do row = 1, beta%n_rows
            pair = beta%cre(row, source) + (beta%des(row, source) - 1)*norb
            target_string = beta%dest(row, source)
            weight = real(beta%phase(row, source), dp)
            do ia = 1, na
               sigma(ia, target_string) = sigma(ia, target_string) &
                                          + weight*contracted(pair, ia + (source - 1)*na)
            end do
         end do
      end do

      deallocate (gathered, contracted)
   end subroutine sigma_vector

   subroutine ci_hamiltonian(folded, alpha, beta, matrix, error)
      !! The Hamiltonian as a dense matrix, one column at a time
      !!
      !! Only for spaces small enough that it fits, which past CAS(8,8) is none
      !! of them -- CAS(10,10) would be 32 GB. It exists because it makes the
      !! sigma build testable against a dense eigensolver, separating "is the
      !! Hamiltonian right" from "does the iterative solver find its lowest
      !! eigenvalue", which are different questions and worth failing
      !! separately. It is also what a Davidson preconditioner will want for a
      !! small subspace of determinants.
      !!
      !! Symmetry is not imposed. It comes out symmetric only if the excitation
      !! tables and the folding are consistent, which makes it worth asserting
      !! rather than assuming.
      real(dp), intent(in) :: folded(:, :)
      type(link_table_t), intent(in) :: alpha, beta
      real(dp), allocatable, intent(out) :: matrix(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: unit_vector(:, :), column(:, :)
      integer :: na, nb, ndet, idet, ia, ib

      if (error%has_error()) return
      na = alpha%n_strings
      nb = beta%n_strings
      ndet = na*nb

      allocate (unit_vector(na, nb), column(na, nb), matrix(ndet, ndet))
      do idet = 1, ndet
         ia = mod(idet - 1, na) + 1
         ib = (idet - 1)/na + 1
         unit_vector = 0.0_dp
         unit_vector(ia, ib) = 1.0_dp
         call sigma_vector(folded, unit_vector, alpha, beta, column, error)
         if (error%has_error()) return
         matrix(:, idet) = reshape(column, [ndet])
      end do

      deallocate (unit_vector, column)
   end subroutine ci_hamiltonian

   subroutine ci_diagonal(h1e, eri, alpha, beta, diagonal, error)
      !! The diagonal of the CI Hamiltonian, without touching the rest of it
      !!
      !! `<D|H|D>` for every determinant, which is `n_det` numbers rather than
      !! `n_det^2` and so is affordable where the matrix is not. A Davidson
      !! preconditioner is the reason it is wanted: dividing a residual by
      !! `theta - H_DD` is what turns an expansion that crawls into one that
      !! converges, and it needs nothing but this.
      !!
      !! For a determinant with alpha orbitals `A` and beta orbitals `B`,
      !!
      !!     sum_{p in A} h_pp + sum_{p in B} h_pp
      !!       + (1/2) sum_{p,q in A} [(pp|qq) - (pq|qp)]
      !!       + (1/2) sum_{p,q in B} [(pp|qq) - (pq|qp)]
      !!       + sum_{p in A} sum_{q in B} (pp|qq)
      !!
      !! -- Coulomb and exchange within each spin, Coulomb only between them,
      !! because two electrons of opposite spin do not exchange.
      !!
      !! Takes the raw integrals, not the folded tensor: the folding rearranges
      !! the Hamiltonian into a form convenient for `E_pq E_rs` and its diagonal
      !! is not the diagonal of anything physical.
      real(dp), intent(in) :: h1e(:, :)
      real(dp), intent(in) :: eri(:, :, :, :)
      type(link_table_t), intent(in) :: alpha, beta
      real(dp), allocatable, intent(out) :: diagonal(:, :)   !! (n_alpha, n_beta)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: coulomb(:, :), exchange(:, :)
      real(dp), allocatable :: one_body(:), same_spin(:)
      integer, allocatable :: occupied(:, :)
      real(dp) :: cross
      integer :: norb, na, nb, p, q, ia, ib, ip, iq

      if (error%has_error()) return
      norb = size(h1e, 1)
      na = alpha%n_strings
      nb = beta%n_strings
      if (beta%n_orbitals /= norb .or. alpha%n_orbitals /= norb) then
         call error%set(ERROR_VALIDATION, "the excitation tables and the integrals "// &
                        "disagree about the number of active orbitals.")
         return
      end if

      ! (pp|qq) and (pq|qp), the only integrals a diagonal element sees.
      allocate (coulomb(norb, norb), exchange(norb, norb))
      do q = 1, norb
         do p = 1, norb
            coulomb(p, q) = eri(p, p, q, q)
            exchange(p, q) = eri(p, q, q, p)
         end do
      end do

      ! The occupied orbitals of each string, read back off the diagonal rows of
      ! its excitation table -- there is exactly one per electron and they are
      ! in ascending order, which `link_table_identities` asserts.
      allocate (one_body(max(na, nb)), same_spin(max(na, nb)))

      call string_terms(alpha, coulomb, exchange, h1e, na, one_body, same_spin, occupied)
      allocate (diagonal(na, nb))
      block
         real(dp), allocatable :: one_b(:), same_b(:)
         integer, allocatable :: occupied_b(:, :)
         allocate (one_b(nb), same_b(nb))
         call string_terms(beta, coulomb, exchange, h1e, nb, one_b, same_b, occupied_b)
         do ib = 1, nb
            do ia = 1, na
               cross = 0.0_dp
               do ip = 1, alpha%n_electrons
                  p = occupied(ip, ia)
                  do iq = 1, beta%n_electrons
                     q = occupied_b(iq, ib)
                     cross = cross + coulomb(p, q)
                  end do
               end do
               diagonal(ia, ib) = one_body(ia) + one_b(ib) + same_spin(ia) &
                                  + same_b(ib) + cross
            end do
         end do
         deallocate (one_b, same_b, occupied_b)
      end block

      deallocate (coulomb, exchange, one_body, same_spin, occupied)
   end subroutine ci_diagonal

   subroutine string_terms(table, coulomb, exchange, h1e, n_str, one_body, same_spin, &
                           occupied)
      !! The part of a diagonal element that depends on one spin alone
      type(link_table_t), intent(in) :: table
      real(dp), intent(in) :: coulomb(:, :), exchange(:, :), h1e(:, :)
      integer, intent(in) :: n_str
      real(dp), intent(out) :: one_body(:), same_spin(:)
      integer, allocatable, intent(out) :: occupied(:, :)

      integer :: istr, ip, iq, p, q

      allocate (occupied(max(table%n_electrons, 1), n_str))
      do istr = 1, n_str
         do ip = 1, table%n_electrons
            occupied(ip, istr) = table%des(ip, istr)
         end do

         one_body(istr) = 0.0_dp
         same_spin(istr) = 0.0_dp
         do ip = 1, table%n_electrons
            p = occupied(ip, istr)
            one_body(istr) = one_body(istr) + h1e(p, p)
            do iq = 1, table%n_electrons
               q = occupied(iq, istr)
               same_spin(istr) = same_spin(istr) &
                                 + 0.5_dp*(coulomb(p, q) - exchange(p, q))
            end do
         end do
      end do
   end subroutine string_terms

end module mqc_ci
