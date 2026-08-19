!! Resolving the energy into atomic and interatomic contributions
module mqc_libcint_ieda
   !! The intrinsic energy decomposition analysis of Del Angel Cruz, Gordon and
   !! Ruedenberg, JACS 147, 42262 (2025). The quasi-atomic orbitals of
   !! `mqc_libcint_quao` give a basis in which every orbital belongs to one
   !! atom; this module sums the energy over that basis and reports what each
   !! atom contributes on its own and what exists only because the atoms are
   !! bonded.
   !!
   !! **This first piece is the kinetic term, and it is not an approximation.**
   !! The total kinetic energy is
   !!
   !!     T = sum_pq gamma_pq T_pq
   !!
   !! over the quasi-atomic basis, and every term in that sum belongs to exactly
   !! one atom (p and q on the same atom) or to exactly one pair of atoms (p and
   !! q on different ones). So the split is a regrouping of a sum, with nothing
   !! discarded and nothing modelled, and the pieces are obliged to add back up.
   !! `kinetic_total` performs that addition and the tests assert it.
   !!
   !! **The factor of two is the trap.** A pair of atoms `{A,B}` collects both
   !! `p in A, q in B` and `p in B, q in A`, which are equal, so the energy of
   !! the pair is twice the one-way sum. `inter` therefore holds the *full*
   !! physical pair energy in both `(A,B)` and `(B,A)`, and `kinetic_total`
   !! halves the matrix rather than summing a triangle. Storing the one-way sum
   !! instead would leave every interatomic number a factor of two small while
   !! still looking entirely plausible, which is why it is stated here and
   !! checked in the tests.
   !!
   !! Related to the kinetic bond order of `kinetic_bond_orders`, but not the
   !! same quantity: that one carries an empirical factor of a tenth to reach
   !! the scale of tabulated bond energies, and the papers say so. Nothing here
   !! is scaled. Where the bond order is a gauge, this is an energy.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_quao, only: quao_result_t
   use mqc_physical_constants, only: HARTREE_TO_KCALMOL
   implicit none
   private

   public :: kinetic_decomposition
   public :: kinetic_total
   public :: print_kinetic_decomposition

   real(dp), parameter :: TO_MILLIHARTREE = 1000.0_dp
   real(dp), parameter :: TOL_SUM_RULE = 1.0e-8_dp
      !! Hartree. How far the regrouped pieces may drift from the total before
      !! the decomposition is refused. The regrouping is exact arithmetic, so
      !! the only thing this tolerates is accumulated rounding; GAMESS applies
      !! the same figure to the equivalent check and aborts on it.
   real(dp), parameter :: TOL_PRINT_PAIR = 1.0e-6_dp
      !! Hartree. Atom pairs contributing less than this are not printed. They
      !! are still counted in the totals, so the printed rows need not add up to
      !! the printed total and the suppressed count says by how much.

contains

   subroutine kinetic_decomposition(quao, interference, n_atoms, intra, inter, error)
      !! Group `gamma_pq T_pq` by the atoms its two orbitals sit on
      !!
      !! `intra(A)` is everything with both orbitals on atom `A`: the kinetic
      !! energy atom `A` would have on its own, plus the intra-atomic
      !! interference its hybridisation in the molecule creates. `inter(A,B)` is
      !! the interatomic kinetic interference of eq (4.26) -- the part of the
      !! kinetic energy that exists only because the two atoms share density,
      !! and in this analysis the dominant source of covalent binding.
      type(quao_result_t), intent(in) :: quao
      real(dp), intent(in) :: interference(:, :)
         !! (n_quao, n_quao). `gamma_pq * T_pq` in hartree, as produced by
         !! `kinetic_bond_orders`. Unscaled -- the kinetic bond order is not
         !! this quantity and must not be passed here.
      integer, intent(in) :: n_atoms
      real(dp), allocatable, intent(out) :: intra(:)
         !! (n_atoms), hartree
      real(dp), allocatable, intent(out) :: inter(:, :)
         !! (n_atoms, n_atoms), hartree. Symmetric, zero on the diagonal, and
         !! each off-diagonal element is the **full** energy of that atom pair.
      type(error_t), intent(inout) :: error

      integer :: n, i, j, a, b
      real(dp) :: residual
      character(len=400) :: message

      if (error%has_error()) return

      n = quao%n_quao
      if (size(interference, 1) /= n .or. size(interference, 2) /= n) then
         call error%set(ERROR_VALIDATION, "the kinetic interference matrix is not "// &
                        "square in the number of quasi-atomic orbitals, so it did "// &
                        "not come from this set of orbitals.")
         return
      end if
      if (n_atoms < 1 .or. maxval(quao%atom_of(1:n)) > n_atoms) then
         call error%set(ERROR_VALIDATION, "a quasi-atomic orbital is assigned to an "// &
                        "atom outside the molecule, so the decomposition would drop "// &
                        "part of the energy silently.")
         return
      end if

      allocate (intra(n_atoms), inter(n_atoms, n_atoms))
      intra = 0.0_dp
      inter = 0.0_dp

      ! Each ordered orbital pair is deposited in **both** elements of its atom
      ! pair, so `inter(A,B)` ends up holding the energy of the whole pair --
      ! `i in A, j in B` and `i in B, j in A` together -- rather than half of
      ! it. That is the number the paper tabulates and the one a reader wants
      ! printed, and it is why `kinetic_total` halves the matrix.
      do i = 1, n
         a = quao%atom_of(i)
         do j = 1, n
            b = quao%atom_of(j)
            if (a == b) then
               intra(a) = intra(a) + interference(i, j)
            else
               inter(a, b) = inter(a, b) + interference(i, j)
               inter(b, a) = inter(b, a) + interference(i, j)
            end if
         end do
      end do

      ! The pieces must add back up to what was handed in. That is the entire
      ! claim this routine makes, it costs one pass over a matrix already in
      ! memory, and it is checked here rather than left to the caller because a
      ! decomposition that quietly loses a term still prints a table that looks
      ! right. It catches a mis-assigned atom, a dropped orbital, and above all
      ! a mishandled factor of two on the interatomic block.
      residual = kinetic_total(intra, inter) - sum(interference)
      if (abs(residual) > TOL_SUM_RULE) then
         write (message, "(a,es12.4,a)") "the kinetic decomposition does not sum "// &
            "back to the kinetic energy it was built from; it is short by ", &
            residual, " hartree. The atom and pair terms are a regrouping of "// &
            "that sum, so any discrepancy means part of the energy was assigned "// &
            "to no atom or to two."
         call error%set(ERROR_VALIDATION, trim(message))
         deallocate (intra, inter)
         return
      end if
   end subroutine kinetic_decomposition

   pure function kinetic_total(intra, inter) result(total)
      !! Add the decomposition back up
      !!
      !! Must reproduce `sum(interference)`, which is the kinetic energy of the
      !! density in the quasi-atomic basis. The only place the halving of
      !! `inter` is written down, so that no caller has to know about it.
      real(dp), intent(in) :: intra(:)
      real(dp), intent(in) :: inter(:, :)
      real(dp) :: total

      total = sum(intra) + 0.5_dp*sum(inter)
   end function kinetic_total

   subroutine print_kinetic_decomposition(intra, inter, element_symbols)
      !! The atom and atom-pair table
      !!
      !! Millihartree rather than kcal/mol as the leading column: these are the
      !! quantities the paper tabulates, and its tables are in millihartree.
      !! kcal/mol follows for comparison with the bond-order table above, which
      !! is on that scale because of the empirical tenth.
      real(dp), intent(in) :: intra(:)
      real(dp), intent(in) :: inter(:, :)
      character(len=*), intent(in) :: element_symbols(:)

      character(len=160) :: line
      character(len=16) :: label_a, label_b
      integer, allocatable :: pair_a(:), pair_b(:), order(:)
      real(dp), allocatable :: magnitude(:)
      integer :: natm, a, b, np, k, printed, suppressed
      real(dp) :: total

      natm = size(intra)

      call logger%info("")
      call logger%info("  kinetic energy decomposition")
      call logger%info("     intra-atomic                    mhartree    kcal/mol")
      do a = 1, natm
         write (label_a, "(a,i0)") trim(adjustl(element_symbols(a)))//" ", a
         write (line, "(4x,a12,18x,f12.3,f12.3)") label_a, &
            intra(a)*TO_MILLIHARTREE, intra(a)*HARTREE_TO_KCALMOL
         call logger%info(trim(line))
      end do

      allocate (pair_a(natm*(natm - 1)/2), pair_b(natm*(natm - 1)/2))
      allocate (magnitude(natm*(natm - 1)/2))
      np = 0
      do a = 1, natm
         do b = a + 1, natm
            np = np + 1
            pair_a(np) = a
            pair_b(np) = b
            magnitude(np) = inter(a, b)
         end do
      end do
      call sort_descending(magnitude(1:np), order)

      call logger%info("")
      call logger%info("     interatomic interference        mhartree    kcal/mol")
      printed = 0
      suppressed = 0
      do k = 1, np
         a = pair_a(order(k))
         b = pair_b(order(k))
         if (abs(inter(a, b)) < TOL_PRINT_PAIR) then
            suppressed = suppressed + 1
            cycle
         end if
         printed = printed + 1
         write (label_a, "(a,i0)") trim(adjustl(element_symbols(a)))//" ", a
         write (label_b, "(a,i0)") trim(adjustl(element_symbols(b)))//" ", b
         write (line, "(4x,a8,a4,a8,6x,f12.3,f12.3)") &
            label_a(1:8), " -- ", label_b(1:8), &
            inter(a, b)*TO_MILLIHARTREE, inter(a, b)*HARTREE_TO_KCALMOL
         call logger%info(trim(line))
      end do
      if (printed == 0) call logger%info("    (none)")
      if (suppressed > 0) then
         write (line, "(a,i0,a)") "    (", suppressed, " pairs below threshold, "// &
            "included in the totals below)"
         call logger%info(trim(line))
      end if

      call logger%info("")
      write (line, "(4x,a24,f16.6,a)") "intra-atomic total", sum(intra), " hartree"
      call logger%info(trim(line))
      write (line, "(4x,a24,f16.6,a)") "interatomic total", 0.5_dp*sum(inter), " hartree"
      call logger%info(trim(line))
      total = kinetic_total(intra, inter)
      write (line, "(4x,a24,f16.6,a)") "kinetic energy", total, " hartree"
      call logger%info(trim(line))

      deallocate (pair_a, pair_b, magnitude, order)
   end subroutine print_kinetic_decomposition

   subroutine sort_descending(values, order)
      !! Indices that put `values` in descending order of magnitude
      real(dp), intent(in) :: values(:)
      integer, allocatable, intent(out) :: order(:)

      integer :: n, i, j, pick, swap
      real(dp) :: best

      n = size(values)
      allocate (order(n))
      do i = 1, n
         order(i) = i
      end do
      do i = 1, n - 1
         pick = i
         best = abs(values(order(i)))
         do j = i + 1, n
            if (abs(values(order(j))) > best) then
               best = abs(values(order(j)))
               pick = j
            end if
         end do
         if (pick /= i) then
            swap = order(i)
            order(i) = order(pick)
            order(pick) = swap
         end if
      end do
   end subroutine sort_descending

end module mqc_libcint_ieda
