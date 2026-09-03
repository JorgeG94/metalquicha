!! The CI Hamiltonian, applied to a vector without ever being formed
module mqc_ci
   !! `sigma = H c` for a complete active space, without ever forming `H`.
   !!
   !! **One contraction, not three.** The two-particle part of the Hamiltonian
   !! is
   !!
   !!     (pq|rs) E_pr,qs = (pq|rs) (E_pq E_rs - delta_qr E_ps)
   !!
   !! so subtracting the `delta_qr` term from the one-electron Hamiltonian and
   !! then folding *that* into the two-electron tensor -- `absorb_one_electron`
   !! below -- leaves the entire Hamiltonian as a single `E_pq E_rs`. Applying
   !! it is then: apply `E_pq` once, contract with the tensor, apply `E_pq`
   !! again. No case analysis, and the same excitation table for both spins.
   !!
   !! The formulation is PySCF's, from `pyscf/fci/direct_spin1.py`; the
   !! implementation is this repository's own.
   !!
   !! The cost is the intermediate, `n_orbitals^2` by the number of
   !! determinants, which is blocked over beta strings.
   use, intrinsic :: iso_fortran_env, only: int64
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
   public :: excitations_block
   public :: beta_strings_per_block

   integer(int64), parameter :: BLOCK_BYTES = 256_int64*1024_int64*1024_int64
      !! Working-set target for one block of the CI intermediate, per buffer
      !! pair. Not a hard limit on the routine: the vector, the diagonal and the
      !! Davidson subspace sit outside it.

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
      !! `k` reproduces `sum f_pq E_pq` exactly, each of the two placements
      !! contributing half.
      !!
      !! The result carries the `1/2` in front of the two-electron energy, so
      !! `sigma_vector` has no factors of its own.
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

   pure function beta_strings_per_block(npair, na, nb) result(width)
      !! How many beta strings one working block may cover
      !!
      !! The intermediate is `(n_orbitals^2, n_determinants)`, or
      !! `norb^2 * na * nb * 8` bytes whole, and the sigma build needs two of
      !! them. Determinants are alpha-major within a beta string, so a range of
      !! beta strings is a contiguous range of columns and the beta string is
      !! the unit blocked on.
      !!
      !! A memory budget, not a fixed count: the cost per beta string grows with
      !! both the orbital count and the alpha string count. A space whose whole
      !! intermediate fits in `BLOCK_BYTES` is not blocked at all.
      integer, intent(in) :: npair, na, nb
      integer :: width

      integer(int64) :: per_string

      ! Two buffers of `npair * na` doubles for each beta string covered: the
      ! gathered intermediate and the contracted one.
      per_string = 2_int64*int(npair, int64)*int(na, int64)*8_int64
      width = int(min(int(nb, int64), max(1_int64, BLOCK_BYTES/per_string)))
   end function beta_strings_per_block

   subroutine excitations_block(ci, alpha, beta, first_beta, last_beta, gathered)
      !! `E_pq |c>` for every orbital pair, over one range of beta strings
      !!
      !! Fills `gathered(npair, na*(last_beta - first_beta + 1))` with the
      !! columns of the full intermediate whose beta index lies in the range.
      !! `apply_excitations` is this over every beta string at once.
      !!
      !! **The two spins block differently.** An alpha excitation leaves the
      !! beta string alone, so it only writes columns whose beta index is
      !! already in the block. A beta excitation *moves* the beta string, so a
      !! contribution landing in this block can come from any source string
      !! anywhere and the whole table is walked with the misses skipped.
      real(dp), intent(in) :: ci(:, :)      !! (n_alpha_strings, n_beta_strings)
      type(link_table_t), intent(in) :: alpha, beta
      integer, intent(in) :: first_beta, last_beta
      real(dp), intent(out) :: gathered(:, :)

      integer :: norb, na, source, target_string, row, pair, ia, ib, column
      integer :: hit, n_hits
      integer, allocatable :: hit_source(:), hit_row(:), hit_target(:)

      norb = alpha%n_orbitals
      na = alpha%n_strings
      gathered = 0.0_dp

      ! Both loops are accumulations, so the parallel axis has to own its output
      ! columns outright; a reduction would need a per-thread copy of the
      ! largest array here, which is the memory the blocking exists to avoid.
      !
      ! An alpha excitation cannot change the beta string, so a column's beta
      ! index is fixed and `ib` owns a whole `na`-wide slice of the output.
      !$omp parallel do default(shared) schedule(static) &
      !$omp   private(ib, source, row, pair, target_string, column)
      do ib = first_beta, last_beta
         do source = 1, na
            do row = 1, alpha%n_rows
               pair = alpha%cre(row, source) + (alpha%des(row, source) - 1)*norb
               target_string = alpha%dest(row, source)
               column = target_string + (ib - first_beta)*na
               gathered(pair, column) = gathered(pair, column) &
                                        + real(alpha%phase(row, source), dp)*ci(source, ib)
            end do
         end do
      end do
      !$omp end parallel do

      ! A beta excitation *moves* the beta string, so the destination is what
      ! picks the column and several sources can land on one. `ib` is therefore
      ! not safe to split on -- but `ia` is, because the column is
      ! `ia + (target - first)*na` and distinct `ia` are distinct columns
      ! whatever the target.
      !
      ! Which sources land inside this block is found once rather than
      ! rediscovered by every `ia`.
      allocate (hit_source(beta%n_strings*beta%n_rows))
      allocate (hit_row(beta%n_strings*beta%n_rows))
      allocate (hit_target(beta%n_strings*beta%n_rows))
      n_hits = 0
      do source = 1, beta%n_strings
         do row = 1, beta%n_rows
            target_string = beta%dest(row, source)
            if (target_string < first_beta .or. target_string > last_beta) cycle
            n_hits = n_hits + 1
            hit_source(n_hits) = source
            hit_row(n_hits) = row
            hit_target(n_hits) = target_string
         end do
      end do

      !$omp parallel do default(shared) schedule(static) &
      !$omp   private(ia, hit, source, row, pair, target_string, column)
      do ia = 1, na
         do hit = 1, n_hits
            source = hit_source(hit)
            row = hit_row(hit)
            target_string = hit_target(hit)
            pair = beta%cre(row, source) + (beta%des(row, source) - 1)*norb
            column = ia + (target_string - first_beta)*na
            gathered(pair, column) = gathered(pair, column) &
                                     + real(beta%phase(row, source), dp)*ci(ia, source)
         end do
      end do
      !$omp end parallel do

      deallocate (hit_source, hit_row, hit_target)
   end subroutine excitations_block

   subroutine apply_excitations(ci, alpha, beta, gathered)
      !! `E_pq |c>` for every orbital pair at once
      !!
      !! Returns `(n_orbitals^2, n_determinants)`, column `d` and row
      !! `p + (q-1) n_orb` holding the determinant-`d` component of `E_pq |c>`.
      !!
      !! Shared with the density matrices, so that the excitation phases are
      !! carried in one place. Alpha excitations move the first index of the
      !! vector and beta the second, which is the only place the two spins
      !! differ anywhere in the CI.
      ! TODO(mqc): `source`, `target_string`, `row`, `pair`, `ia`, `ib` and
      ! `weight` are left over from before the body became one call to
      ! `excitations_block`, and none of them is used.
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
      call excitations_block(ci, alpha, beta, 1, nb, gathered)
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
      !! adjoint are the same table read in opposite directions. Alpha
      !! excitations move the first index of the vector and beta excitations the
      !! second, which is the only place the two spins differ.
      real(dp), intent(in) :: folded(:, :)   !! From `absorb_one_electron`
      real(dp), intent(in) :: ci(:, :)       !! (n_alpha_strings, n_beta_strings)
      type(link_table_t), intent(in) :: alpha, beta
      real(dp), intent(out) :: sigma(:, :)   !! Same shape as `ci`
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: gathered(:, :), contracted(:, :)
      integer :: norb, na, nb, npair, ndet
      integer :: source, target_string, row, pair, ia, ib
      integer :: per_block, first, last, width
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

      ! Blocked over beta strings. The contraction below is column-wise, so it
      ! splits exactly; the scatters accumulate into `sigma`, which is the size
      ! of the vector rather than of the intermediate and stays whole.
      per_block = beta_strings_per_block(npair, na, nb)
      allocate (gathered(npair, na*per_block), contracted(npair, na*per_block))

      sigma = 0.0_dp
      do first = 1, nb, per_block
         last = min(first + per_block - 1, nb)
         width = na*(last - first + 1)

         call excitations_block(ci, alpha, beta, first, last, gathered(:, 1:width))
         call pic_gemm(folded, gathered(:, 1:width), contracted(:, 1:width))

         ! Alpha excitations leave the beta string alone, so the columns read
         ! are those of this block and the range is simply narrowed. Split on
         ! `ib`, which indexes a whole column of `sigma` that the iteration then
         ! owns outright, so the accumulation cannot collide.
         !$omp parallel do default(shared) schedule(static) &
         !$omp   private(ib, source, row, pair, target_string, weight)
         do ib = first, last
            do source = 1, na
               do row = 1, alpha%n_rows
                  pair = alpha%cre(row, source) + (alpha%des(row, source) - 1)*norb
                  target_string = alpha%dest(row, source)
                  weight = real(alpha%phase(row, source), dp)
                  sigma(target_string, ib) = sigma(target_string, ib) &
                                             + weight*contracted(pair, source + (ib - first)*na)
               end do
            end do
         end do
         !$omp end parallel do

         ! A beta excitation reads the column of its *source* string and writes
         ! wherever that string lands, possibly outside the block, which is why
         ! `sigma` accumulates across blocks rather than being filled one block
         ! at a time. Splitting on the source would race, since several sources
         ! can land on one target; splitting on `ia` fixes the row of `sigma`
         ! written and distinct rows are disjoint.
         !$omp parallel do default(shared) schedule(static) &
         !$omp   private(ia, source, row, pair, target_string, weight)
         do ia = 1, na
            do source = first, last
               do row = 1, beta%n_rows
                  pair = beta%cre(row, source) + (beta%des(row, source) - 1)*norb
                  target_string = beta%dest(row, source)
                  weight = real(beta%phase(row, source), dp)
                  sigma(ia, target_string) = sigma(ia, target_string) &
                                             + weight*contracted(pair, ia + (source - first)*na)
               end do
            end do
         end do
         !$omp end parallel do
      end do

      deallocate (gathered, contracted)
   end subroutine sigma_vector

   subroutine ci_hamiltonian(folded, alpha, beta, matrix, error)
      !! The Hamiltonian as a dense matrix, one column at a time
      !!
      !! One sigma build per determinant and an `n_det` square matrix, so only
      !! for spaces small enough to hold both. Symmetry is not imposed: the
      !! matrix comes out symmetric only if the excitation tables and the
      !! folding are consistent, so it is worth asserting rather than assuming.
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
      !! `<D|H|D>` for every determinant: `n_det` numbers rather than `n_det^2`,
      !! and all a Davidson preconditioner needs.
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
      !! Takes the raw integrals, **not** the folded tensor: folding rearranges
      !! the Hamiltonian for `E_pq E_rs` and its diagonal is not the diagonal of
      !! anything physical.
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
