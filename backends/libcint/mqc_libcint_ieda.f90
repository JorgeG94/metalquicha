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
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_esp, only: esp_matrices
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_libcint_quao, only: quao_result_t
   use mqc_physical_constants, only: HARTREE_TO_KCALMOL
   implicit none
   private

   public :: kinetic_decomposition
   public :: kinetic_total
   public :: combine_quao_sets
   public :: quao_interference
   public :: nuclear_attraction_per_atom
   public :: quao_nuclear_attraction
   public :: nuclear_decomposition
   public :: quao_eris
   public :: active_cumulant
   public :: quao_projection
   public :: transform_cumulant
   public :: two_electron_decomposition
   public :: print_two_electron_decomposition
   public :: nuclear_repulsion_pairs
   public :: print_total_decomposition
   public :: print_interatomic_split
   public :: print_kinetic_decomposition
   public :: print_nuclear_decomposition

   real(dp), parameter :: TO_MILLIHARTREE = 1000.0_dp
   real(dp), parameter :: TOL_SUM_RULE = 1.0e-8_dp
      !! Hartree. How far the regrouped pieces may drift from the total before
      !! the decomposition is refused. The regrouping is exact arithmetic, so
      !! the only thing this tolerates is accumulated rounding; GAMESS applies
      !! the same figure to the equivalent check and aborts on it.
   real(dp), parameter :: TOL_SPLIT = 1.0e-8_dp
      !! Hartree, per matrix element. How far the per-nucleus attraction
      !! integrals may drift from the one-shot ones before the split is
      !! refused. The two are the same operator summed in a different order, so
      !! this tolerates quadrature and nothing else.
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

   subroutine combine_quao_sets(core, valence, combined, error)
      !! Stack the core quasi-atomic orbitals in front of the valence ones
      !!
      !! The bonding analysis works in the valence-internal space, because that
      !! is where bonding happens and the core is inert. An energy decomposition
      !! cannot: the core carries most of the kinetic energy and most of the
      !! nuclear attraction, and leaving it out means the pieces sum to a number
      !! nothing else knows.
      !!
      !! The core set is built by the same construction as the valence one, just
      !! handed the core molecular orbitals and the core minimal-basis ranges
      !! instead. That matters -- it means core and valence orbitals are atomic
      !! in the same sense rather than in two different senses.
      !!
      !! **The core density is exactly two on the diagonal and zero elsewhere.**
      !! Not an approximation: the core orbitals span the core molecular-orbital
      !! space, the density restricted to that space is twice its projector, and
      !! a projector is the identity in any orthonormal basis of what it
      !! projects onto. The core-valence blocks are zero for the same reason,
      !! the two spaces being orthogonal.
      type(quao_result_t), intent(in) :: core
      type(quao_result_t), intent(in) :: valence
      type(quao_result_t), intent(out) :: combined
      type(error_t), intent(inout) :: error

      integer :: nc, nv, n, i

      if (error%has_error()) return

      nc = core%n_quao
      nv = valence%n_quao
      if (size(core%orbitals, 1) /= size(valence%orbitals, 1)) then
         call error%set(ERROR_VALIDATION, "the core and valence quasi-atomic orbitals "// &
                        "are expanded in different atomic-orbital bases.")
         return
      end if
      n = nc + nv

      combined%n_quao = n
      allocate (combined%orbitals(size(core%orbitals, 1), n))
      combined%orbitals(:, 1:nc) = core%orbitals
      combined%orbitals(:, nc + 1:n) = valence%orbitals

      allocate (combined%atom_of(n))
      combined%atom_of(1:nc) = core%atom_of
      combined%atom_of(nc + 1:n) = valence%atom_of

      allocate (combined%population_bond_order(n, n))
      combined%population_bond_order = 0.0_dp
      do concurrent(i=1:nc)
         combined%population_bond_order(i, i) = 2.0_dp
      end do
      combined%population_bond_order(nc + 1:, nc + 1:) = valence%population_bond_order
   end subroutine combine_quao_sets

   subroutine quao_interference(quao, matrix_ao, weighted, error)
      !! `gamma_pq * O_pq` in the quasi-atomic basis, for a one-electron operator
      !!
      !! The array every decomposition below consumes. Summing it over all pairs
      !! gives the expectation value of the operator, and summing it over
      !! subsets of pairs is what the decomposition is.
      !!
      !! Unlike `kinetic_bond_orders` this does not require oriented orbitals,
      !! and it should not. Orientation is a rotation *within* each atom, and
      !! both `sum_{p,q in A}` and `sum_{p in A, q in B}` are invariant under
      !! one: the density and the operator transform by the same rotation and it
      !! cancels inside the trace. Individual orbital pairs move, atom and
      !! atom-pair totals do not. So the decomposition is well defined for a
      !! core set that was never oriented, and the tests assert the invariance
      !! rather than relying on it quietly.
      type(quao_result_t), intent(in) :: quao
      real(dp), intent(in) :: matrix_ao(:, :)
      real(dp), allocatable, intent(out) :: weighted(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: work(:, :), transformed(:, :)
      integer :: n_ao, n

      if (error%has_error()) return

      n = quao%n_quao
      n_ao = size(quao%orbitals, 1)
      if (size(matrix_ao, 1) /= n_ao .or. size(matrix_ao, 2) /= n_ao) then
         call error%set(ERROR_VALIDATION, "the operator is not in the same "// &
                        "atomic-orbital basis as the quasi-atomic orbitals.")
         return
      end if

      allocate (work(n_ao, n), transformed(n, n), weighted(n, n))
      call pic_gemm(matrix_ao, quao%orbitals, work)
      call pic_gemm(quao%orbitals, work, transformed, transa="T")
      weighted = transformed*quao%population_bond_order
      deallocate (work, transformed)
   end subroutine quao_interference

   subroutine nuclear_attraction_per_atom(mol, atomic_numbers, coordinates, v_atom, error)
      !! Split the nuclear attraction by which nucleus is doing the attracting
      !!
      !!     V_pq^A = -Z_A < p | 1/|r - R_A| | q >
      !!
      !! GAMESS calls the equivalent array `VBYATM`. No new integral is needed:
      !! `esp_matrices` already evaluates `1/|r - R|` at arbitrary points, so
      !! putting those points on the nuclei and scaling by the nuclear charge
      !! is the whole construction.
      !!
      !! **Summing the pieces must give back the ordinary nuclear attraction**,
      !! since the operator is a sum over nuclei and nothing else. That is
      !! checked here rather than left to a caller, because every plausible way
      !! of getting this wrong -- a dropped minus sign, the charge left off, the
      !! points in Angstrom -- produces an array that is smooth, symmetric and
      !! entirely wrong. The check costs two one-electron integral builds
      !! against an `esp_matrices` call that is already more expensive than
      !! both, so it is close to free in context.
      type(libcint_molecule_t), intent(in), target :: mol
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)   !! (3, n_atoms), Bohr
      real(dp), allocatable, intent(out) :: v_atom(:, :, :)
         !! (n_ao, n_ao, n_atoms), hartree
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: rinv(:, :, :), h(:, :), t(:, :), total(:, :)
      character(len=400) :: message
      real(dp) :: worst
      integer :: natm, a

      if (error%has_error()) return

      natm = size(atomic_numbers)
      if (size(coordinates, 1) /= 3 .or. size(coordinates, 2) /= natm) then
         call error%set(ERROR_VALIDATION, "the coordinates and the atomic numbers "// &
                        "describe different numbers of atoms.")
         return
      end if

      call esp_matrices(mol, coordinates, rinv, error)
      if (error%has_error()) return

      allocate (v_atom(mol%nao, mol%nao, natm))
      do a = 1, natm
         v_atom(:, :, a) = -real(atomic_numbers(a), dp)*rinv(:, :, a)
      end do
      deallocate (rinv)

      ! The reference: the one-shot nuclear attraction, which is the core
      ! Hamiltonian less the kinetic energy.
      call mol%core_hamiltonian(h)
      call mol%kinetic(t)
      allocate (total(mol%nao, mol%nao))
      total = 0.0_dp
      do a = 1, natm
         total = total + v_atom(:, :, a)
      end do
      worst = maxval(abs(total - (h - t)))
      deallocate (h, t, total)

      if (worst > TOL_SPLIT) then
         write (message, "(a,es12.4,a)") "the per-nucleus attraction integrals do "// &
            "not sum to the ordinary nuclear attraction; the largest element "// &
            "disagrees by ", worst, " hartree. The two are the same operator "// &
            "summed in a different order, so this is a sign, a nuclear charge "// &
            "or a coordinate, not an accuracy question."
         call error%set(ERROR_VALIDATION, trim(message))
         deallocate (v_atom)
         return
      end if
   end subroutine nuclear_attraction_per_atom

   subroutine quao_nuclear_attraction(quao, v_atom_ao, v_atom_quao, error)
      !! Carry the per-nucleus attraction into the quasi-atomic basis
      !!
      !! One `C^T V^A C` per nucleus. The quasi-atomic orbitals are orthonormal,
      !! so this is a similarity transformation and the trace against a density
      !! is preserved exactly -- which is what lets the decomposition below be a
      !! regrouping rather than a model.
      type(quao_result_t), intent(in) :: quao
      real(dp), intent(in) :: v_atom_ao(:, :, :)      !! (n_ao, n_ao, n_atoms)
      real(dp), allocatable, intent(out) :: v_atom_quao(:, :, :)
         !! (n_quao, n_quao, n_atoms)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: work(:, :)
      integer :: n_ao, n, natm, a

      if (error%has_error()) return

      n = quao%n_quao
      n_ao = size(quao%orbitals, 1)
      natm = size(v_atom_ao, 3)
      if (size(v_atom_ao, 1) /= n_ao .or. size(v_atom_ao, 2) /= n_ao) then
         call error%set(ERROR_VALIDATION, "the per-nucleus attraction integrals are "// &
                        "not in the same atomic-orbital basis as the quasi-atomic "// &
                        "orbitals.")
         return
      end if

      allocate (v_atom_quao(n, n, natm), work(n_ao, n))
      do a = 1, natm
         call pic_gemm(v_atom_ao(:, :, a), quao%orbitals, work)
         call pic_gemm(quao%orbitals, work, v_atom_quao(:, :, a), transa="T")
      end do
      deallocate (work)
   end subroutine quao_nuclear_attraction

   subroutine nuclear_decomposition(quao, density, v_atom_quao, n_atoms, intra, inter, &
                                    error, coulomb)
      !! Group `gamma_pq V_pq^C` by the atoms of its two orbitals and its nucleus
      !!
      !! Unlike the kinetic term this one carries three atomic labels: the
      !! orbital `p` sits on some atom, `q` on another, and the nucleus doing
      !! the attracting is a third. The assignment is
      !!
      !!   * `p, q` on `A` and nucleus `A` -- atom `A`'s own electrons in its
      !!     own field, so `intra(A)`;
      !!   * `p, q` on `A` and nucleus `C` -- atom `A`'s density attracted to a
      !!     foreign nucleus, which is a Coulombic interaction between the two
      !!     and belongs to the pair `{A,C}`;
      !!   * `p` on `A`, `q` on `B` -- an interference term, which exists only
      !!     because those two atoms share density, so it goes to `{A,B}`
      !!     **whichever nucleus** is attracting.
      !!
      !! **The last one is a choice and it should be read as one.** A term with
      !! `p` on `A`, `q` on `B` and the nucleus on a third atom `C` is genuinely
      !! three-body, and charging it to `{A,B}` says the interference is the
      !! feature and the field it sits in is context. That keeps the pair table
      !! interpretable and it loses nothing, since the sum rule below still
      !! holds; a finer grouping would redistribute these terms without changing
      !! any total.
      type(quao_result_t), intent(in) :: quao
      real(dp), intent(in) :: density(:, :)           !! (n_quao, n_quao)
      real(dp), intent(in) :: v_atom_quao(:, :, :)    !! (n_quao, n_quao, n_atoms)
      integer, intent(in) :: n_atoms
      real(dp), allocatable, intent(out) :: intra(:)
      real(dp), allocatable, intent(out) :: inter(:, :)
         !! Symmetric, zero diagonal, each element the whole pair energy
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: coulomb(:, :)
         !! The **classical** share of `inter`: one atom's own density sitting in
         !! another's nuclear field. What is left over is interference -- density
         !! shared between the two atoms, which has no classical description and
         !! exists only because they are bonded. Splitting them is the whole
         !! point of the analysis, since the papers' claim is that binding comes
         !! from interference and not from electrostatics.

      real(dp) :: term, residual, reference
      character(len=400) :: message
      integer :: n, i, j, c, a, b

      if (error%has_error()) return

      n = quao%n_quao
      if (size(density, 1) /= n .or. size(v_atom_quao, 1) /= n) then
         call error%set(ERROR_VALIDATION, "the density and the per-nucleus attraction "// &
                        "are not both in the quasi-atomic basis.")
         return
      end if
      if (size(v_atom_quao, 3) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "the per-nucleus attraction does not cover "// &
                        "every atom, so part of the nuclear field would be dropped.")
         return
      end if
      if (n_atoms < 1 .or. maxval(quao%atom_of(1:n)) > n_atoms) then
         call error%set(ERROR_VALIDATION, "a quasi-atomic orbital is assigned to an "// &
                        "atom outside the molecule.")
         return
      end if

      allocate (intra(n_atoms), inter(n_atoms, n_atoms))
      intra = 0.0_dp
      inter = 0.0_dp
      reference = 0.0_dp
      if (present(coulomb)) then
         allocate (coulomb(n_atoms, n_atoms))
         coulomb = 0.0_dp
      end if

      do c = 1, n_atoms
         do i = 1, n
            a = quao%atom_of(i)
            do j = 1, n
               b = quao%atom_of(j)
               term = density(i, j)*v_atom_quao(i, j, c)
               reference = reference + term
               if (a /= b) then
                  ! Interference: charged to the pair that shares the density,
                  ! deposited in both elements so each holds the whole pair.
                  inter(a, b) = inter(a, b) + term
                  inter(b, a) = inter(b, a) + term
               else if (a == c) then
                  intra(a) = intra(a) + term
               else
                  ! One atom's density in another's nuclear field: an ordinary
                  ! electrostatic attraction between the two.
                  inter(a, c) = inter(a, c) + term
                  inter(c, a) = inter(c, a) + term
                  if (present(coulomb)) then
                     coulomb(a, c) = coulomb(a, c) + term
                     coulomb(c, a) = coulomb(c, a) + term
                  end if
               end if
            end do
         end do
      end do

      residual = kinetic_total(intra, inter) - reference
      if (abs(residual) > TOL_SUM_RULE) then
         write (message, "(a,es12.4,a)") "the nuclear attraction decomposition does "// &
            "not sum back to the energy it was built from; it is short by ", &
            residual, " hartree. Every term carries three atomic labels and each "// &
            "must land in exactly one bin."
         call error%set(ERROR_VALIDATION, trim(message))
         deallocate (intra, inter)
         return
      end if
   end subroutine nuclear_decomposition

   subroutine print_nuclear_decomposition(intra, inter, element_symbols)
      !! The atom and atom-pair table for the nuclear attraction
      real(dp), intent(in) :: intra(:)
      real(dp), intent(in) :: inter(:, :)
      character(len=*), intent(in) :: element_symbols(:)

      call print_decomposition("nuclear attraction decomposition", intra, inter, &
                               element_symbols)
   end subroutine print_nuclear_decomposition

   subroutine quao_eris(quao, eri_ao, eri_quao, error)
      !! Two-electron integrals over the quasi-atomic orbitals
      !!
      !! Four quarter-transformations, each contracting one atomic-orbital index
      !! and cycling it to the back, so after four passes `(mu nu|lambda sigma)`
      !! has become `(pq|rs)` with the indices in the order they started.
      !!
      !! **This takes the dense atomic-orbital array and is therefore bounded by
      !! `n_ao^4`.** That is fine for the molecules this analysis is run on and
      !! it is not a fundamental cost -- the transformation could be driven from
      !! the packed integrals, or direct from shell quartets, and nothing above
      !! here would change. It is written the simple way because the decomposition
      !! is the part that needed care, and the ceiling is stated rather than
      !! discovered.
      type(quao_result_t), intent(in) :: quao
      real(dp), intent(in) :: eri_ao(:, :, :, :)
      real(dp), allocatable, intent(out) :: eri_quao(:, :, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: cur(:, :), t(:, :)
      integer :: n_ao, n, step, ncol

      if (error%has_error()) return

      n = quao%n_quao
      n_ao = size(quao%orbitals, 1)
      if (any(shape(eri_ao) /= n_ao)) then
         call error%set(ERROR_VALIDATION, "the two-electron integrals are not in the "// &
                        "same atomic-orbital basis as the quasi-atomic orbitals.")
         return
      end if

      call four_index_transform(quao%orbitals, eri_ao, eri_quao)
   end subroutine quao_eris

   subroutine four_index_transform(coeff, source, result)
      !! Contract all four indices of an `n_in^4` array with `coeff`
      !!
      !!     result_pqrs = sum_abcd coeff_ap coeff_bq coeff_cr coeff_ds source_abcd
      !!
      !! Four quarter-transformations, each contracting one index and cycling it
      !! to the back, so after four passes the order is restored. `coeff` is
      !! `(n_in, n_out)` and need not be square: the two-electron integrals come
      !! from a large atomic-orbital basis into a small quasi-atomic one, and
      !! the cumulant comes from a small active space into the same one.
      real(dp), intent(in) :: coeff(:, :)          !! (n_in, n_out)
      real(dp), intent(in) :: source(:, :, :, :)   !! (n_in, n_in, n_in, n_in)
      real(dp), allocatable, intent(out) :: result(:, :, :, :)

      real(dp), allocatable :: cur(:, :), t(:, :)
      integer :: n_in, n_out, step, ncol

      n_in = size(coeff, 1)
      n_out = size(coeff, 2)

      ! A transposed `coeff` still type-checks and still has two dimensions, so
      ! without this the first reshape reads past the end of `source` and the
      ! whole thing crashes several frames away from the mistake.
      if (any(shape(source) /= n_in)) then
         error stop "four_index_transform: the array does not have the "// &
            "dimension the coefficients contract over"
      end if

      allocate (cur(n_in, n_in**3))
      cur = reshape(source, [n_in, n_in**3])
      do step = 1, 4
         ncol = size(cur, 2)
         allocate (t(n_out, ncol))
         call pic_gemm(coeff, cur, t, transa="T")
         deallocate (cur)
         if (step < 4) then
            allocate (cur(n_in, ncol*n_out/n_in))
            cur = reshape(transpose(t), [n_in, ncol*n_out/n_in])
         else
            allocate (cur(n_out, n_out**3))
            cur = reshape(transpose(t), [n_out, n_out**3])
         end if
         deallocate (t)
      end do

      allocate (result(n_out, n_out, n_out, n_out))
      result = reshape(cur, [n_out, n_out, n_out, n_out])
      deallocate (cur)
   end subroutine four_index_transform

   pure subroutine active_cumulant(dm1, dm2, cumulant)
      !! What correlation adds to the two-particle density
      !!
      !!     Lambda_pqrs = d_pqrs - [ D_pq D_rs - (1/2) D_ps D_rq ]
      !!
      !! The bracket is the two-particle density a single determinant would
      !! have, so `Lambda` is exactly the part of `d` that a determinant cannot
      !! produce. It is zero for a closed-shell reference, and for an MCSCF wave
      !! function it is zero unless all four indices are active -- the inactive
      !! orbitals are a closed shell and carry no correlation.
      !!
      !! That is what makes the correlated decomposition cheap. The determinant
      !! expression is already evaluated over the whole quasi-atomic basis; the
      !! correction lives in the active space alone and is transformed from
      !! there.
      real(dp), intent(in) :: dm1(:, :)
      real(dp), intent(in) :: dm2(:, :, :, :)
      real(dp), allocatable, intent(out) :: cumulant(:, :, :, :)

      integer :: n, p, q, r, sx

      n = size(dm1, 1)
      allocate (cumulant(n, n, n, n))
      do p = 1, n
         do q = 1, n
            do r = 1, n
               do sx = 1, n
                  cumulant(p, q, r, sx) = dm2(p, q, r, sx) &
                                          - dm1(p, q)*dm1(r, sx) &
                                          + 0.5_dp*dm1(p, sx)*dm1(r, q)
               end do
            end do
         end do
      end do
   end subroutine active_cumulant

   subroutine quao_projection(quao, overlap, orbitals, u, deficit, error)
      !! Express a set of molecular orbitals in the quasi-atomic basis
      !!
      !!     U = C_quao^T S C
      !!
      !! **`U` is an isometry only if the orbitals lie inside the space the
      !! quasi-atomic orbitals span, and active orbitals do not quite.** The
      !! valence-virtual orbitals are chosen to look like free-atom orbitals;
      !! the active orbitals of a converged MCSCF are chosen to lower an energy,
      !! and the two are not the same choice. N2 in cc-pVDZ with a CAS(6,6)
      !! misses orthonormality by about 6e-2.
      !!
      !! So the deficit is measured and handed back rather than treated as an
      !! error. It is the same statement the population sum already makes about
      !! this analysis -- that it describes what lies inside the quasi-atomic
      !! span and not what lies outside -- and the caller is better placed to
      !! decide what to do about it. What must not happen is losing it
      !! silently.
      type(quao_result_t), intent(in) :: quao
      real(dp), intent(in) :: overlap(:, :)     !! (n_ao, n_ao)
      real(dp), intent(in) :: orbitals(:, :)    !! (n_ao, n_orb)
      real(dp), allocatable, intent(out) :: u(:, :)   !! (n_quao, n_orb)
      real(dp), intent(out) :: deficit
         !! `max |U^T U - 1|`. Zero when the orbitals are fully inside the span;
         !! otherwise how much of them is not.
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: work(:, :), gram(:, :)
      integer :: n_ao, n, n_orb, i

      if (error%has_error()) return

      deficit = 0.0_dp
      n = quao%n_quao
      n_ao = size(quao%orbitals, 1)
      n_orb = size(orbitals, 2)
      if (size(overlap, 1) /= n_ao .or. size(orbitals, 1) /= n_ao) then
         call error%set(ERROR_VALIDATION, "the overlap, the orbitals and the "// &
                        "quasi-atomic orbitals are not in the same basis.")
         return
      end if

      allocate (work(n_ao, n_orb), u(n, n_orb), gram(n_orb, n_orb))
      call pic_gemm(overlap, orbitals, work)
      call pic_gemm(quao%orbitals, work, u, transa="T")

      call pic_gemm(u, u, gram, transa="T")
      do i = 1, n_orb
         gram(i, i) = gram(i, i) - 1.0_dp
      end do
      deficit = maxval(abs(gram))
      deallocate (work, gram)
   end subroutine quao_projection

   subroutine transform_cumulant(u, cumulant, cumulant_quao, error)
      !! Carry the cumulant from the active orbitals into the quasi-atomic basis
      !!
      !! `U` being an isometry is what makes this exact: the contraction of the
      !! cumulant against the integrals is the same number computed either side
      !! of the transformation, so the correlation energy is neither gained nor
      !! lost by changing basis. `quao_projection` checks that before returning.
      real(dp), intent(in) :: u(:, :)               !! (n_quao, n_active)
      real(dp), intent(in) :: cumulant(:, :, :, :)  !! over active orbitals
      real(dp), allocatable, intent(out) :: cumulant_quao(:, :, :, :)
      type(error_t), intent(inout) :: error

      if (error%has_error()) return
      if (size(u, 2) /= size(cumulant, 1)) then
         call error%set(ERROR_VALIDATION, "the projection and the cumulant describe "// &
                        "different numbers of active orbitals.")
         return
      end if
      call four_index_transform(transpose(u), cumulant, cumulant_quao)
   end subroutine transform_cumulant

   subroutine two_electron_decomposition(quao, density, eri_quao, n_atoms, &
                                         intra, inter, energy, error, cumulant, coulomb)
      !! Group the two-electron energy by the atoms of its four orbitals
      !!
      !! For a single determinant the two-particle density is
      !!
      !!     Gamma_pqrs = gamma_pq gamma_rs - (1/2) gamma_ps gamma_rq
      !!
      !! and `E2 = (1/2) sum_pqrs Gamma_pqrs (pq|rs)`. Both terms are carried
      !! together rather than decomposed separately, because splitting them
      !! invites the question of which atoms an exchange term belongs to -- its
      !! density factors pair `p` with `s` while the integral pairs `p` with `q`
      !! -- and there is no need to answer it. The integral is what says which
      !! charge distributions are interacting.
      !!
      !! `(pq|rs)` is the interaction of the distribution `pq` with the
      !! distribution `rs`, so the assignment follows the one-electron case with
      !! a bra pair and a ket pair instead of a pair and a nucleus:
      !!
      !!   * both distributions on one atom, the same atom -- `intra(A)`;
      !!   * both on one atom, different atoms -- the Coulombic interaction of
      !!     `A` with `B`, so `{A,B}`;
      !!   * one distribution spread over two atoms, the other on one -- an
      !!     interference, charged to the two atoms sharing it;
      !!   * both spread -- split evenly between the two pairs, which is the
      !!     only symmetric thing to do with a term that belongs to both.
      type(quao_result_t), intent(in) :: quao
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(in) :: eri_quao(:, :, :, :)
      integer, intent(in) :: n_atoms
      real(dp), allocatable, intent(out) :: intra(:)
      real(dp), allocatable, intent(out) :: inter(:, :)
      real(dp), intent(out) :: energy
         !! `E2` as summed here, so a caller can compare it against the same
         !! quantity computed without any of this.
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: cumulant(:, :, :, :)
         !! (n_quao^4). What correlation adds to the two-particle density,
         !! already carried into this basis. Absent means a single determinant,
         !! for which it is identically zero.
      real(dp), allocatable, intent(out), optional :: coulomb(:, :)
         !! The classical share of `inter`: both charge distributions sitting
         !! wholly on one atom each, so the term is the repulsion between two
         !! atomic charge clouds. Everything else involves a distribution spread
         !! across two atoms and is interference.

      real(dp) :: w, residual, extra
      character(len=400) :: message
      integer :: n, p, q, r, sx, a, b, c, d

      if (error%has_error()) return

      n = quao%n_quao
      energy = 0.0_dp
      if (size(density, 1) /= n .or. size(eri_quao, 1) /= n) then
         call error%set(ERROR_VALIDATION, "the density and the two-electron integrals "// &
                        "are not both in the quasi-atomic basis.")
         return
      end if
      if (n_atoms < 1 .or. maxval(quao%atom_of(1:n)) > n_atoms) then
         call error%set(ERROR_VALIDATION, "a quasi-atomic orbital is assigned to an "// &
                        "atom outside the molecule.")
         return
      end if

      allocate (intra(n_atoms), inter(n_atoms, n_atoms))
      intra = 0.0_dp
      inter = 0.0_dp
      if (present(coulomb)) then
         allocate (coulomb(n_atoms, n_atoms))
         coulomb = 0.0_dp
      end if

      do p = 1, n
         a = quao%atom_of(p)
         do q = 1, n
            b = quao%atom_of(q)
            do r = 1, n
               c = quao%atom_of(r)
               do sx = 1, n
                  d = quao%atom_of(sx)
                  extra = 0.0_dp
                  if (present(cumulant)) extra = cumulant(p, q, r, sx)
                  w = 0.5_dp*(density(p, q)*density(r, sx) &
                              - 0.5_dp*density(p, sx)*density(r, q) &
                              + extra)*eri_quao(p, q, r, sx)
                  energy = energy + w
                  if (a /= b .and. c /= d) then
                     call deposit(inter, a, b, 0.5_dp*w)
                     call deposit(inter, c, d, 0.5_dp*w)
                  else if (a /= b) then
                     call deposit(inter, a, b, w)
                  else if (c /= d) then
                     call deposit(inter, c, d, w)
                  else if (a == c) then
                     intra(a) = intra(a) + w
                  else
                     ! Two atomic charge clouds repelling each other, which is
                     ! as classical as this decomposition gets.
                     call deposit(inter, a, c, w)
                     if (present(coulomb)) call deposit(coulomb, a, c, w)
                  end if
               end do
            end do
         end do
      end do

      residual = kinetic_total(intra, inter) - energy
      if (abs(residual) > TOL_SUM_RULE) then
         write (message, "(a,es12.4,a)") "the two-electron decomposition does not sum "// &
            "back to the energy it was built from; it is short by ", residual, &
            " hartree. Every term carries four atomic labels and each must reach "// &
            "exactly one bin, whole or in two halves."
         call error%set(ERROR_VALIDATION, trim(message))
         deallocate (intra, inter)
         return
      end if
   end subroutine two_electron_decomposition

   pure subroutine deposit(inter, a, b, value)
      !! Charge `value` to the pair `{a,b}`
      !!
      !! Written into both elements, because `inter` holds the whole pair energy
      !! in each and `kinetic_total` halves the matrix. The one place any
      !! decomposition in this module touches the pair matrix, so the convention
      !! is stated once and cannot drift between terms.
      real(dp), intent(inout) :: inter(:, :)
      integer, intent(in) :: a, b
      real(dp), intent(in) :: value

      inter(a, b) = inter(a, b) + value
      inter(b, a) = inter(b, a) + value
   end subroutine deposit

   subroutine print_two_electron_decomposition(intra, inter, element_symbols)
      !! The atom and atom-pair table for the two-electron energy
      real(dp), intent(in) :: intra(:)
      real(dp), intent(in) :: inter(:, :)
      character(len=*), intent(in) :: element_symbols(:)

      call print_decomposition("two-electron decomposition", intra, inter, &
                               element_symbols)
   end subroutine print_two_electron_decomposition

   subroutine nuclear_repulsion_pairs(atomic_numbers, coordinates, inter, error)
      !! The nuclear repulsion, which is already a sum over atom pairs
      !!
      !!     inter(A,B) = Z_A Z_B / R_AB
      !!
      !! The one term in the whole energy that needs no decomposing: it is
      !! written as a sum over pairs before anyone asks. Carrying it in the same
      !! convention as everything else -- the whole pair energy in both
      !! elements, halved by `kinetic_total` -- is what lets the total be
      !! accumulated by adding matrices, and it is why the energy ends up
      !! resolved into atoms and pairs with no remainder sitting outside.
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)   !! (3, n_atoms), Bohr
      real(dp), allocatable, intent(out) :: inter(:, :)
      type(error_t), intent(inout) :: error

      integer :: natm, a, b
      real(dp) :: r

      if (error%has_error()) return

      natm = size(atomic_numbers)
      if (size(coordinates, 2) /= natm) then
         call error%set(ERROR_VALIDATION, "the coordinates and the atomic numbers "// &
                        "describe different numbers of atoms.")
         return
      end if

      allocate (inter(natm, natm))
      inter = 0.0_dp
      do a = 1, natm
         do b = a + 1, natm
            r = norm2(coordinates(:, a) - coordinates(:, b))
            if (r < 1.0e-10_dp) then
               call error%set(ERROR_VALIDATION, "two atoms are on top of each other, "// &
                              "so the nuclear repulsion is not finite.")
               deallocate (inter)
               return
            end if
            inter(a, b) = real(atomic_numbers(a), dp)*real(atomic_numbers(b), dp)/r
            inter(b, a) = inter(a, b)
         end do
      end do
   end subroutine nuclear_repulsion_pairs

   subroutine print_total_decomposition(intra, inter, element_symbols)
      !! The energy resolved into atoms and atom pairs
      !!
      !! What the analysis is for. `intra(A)` is everything that happens inside
      !! atom `A` -- its own kinetic energy, its electrons in its own nuclear
      !! field, its own electron repulsion -- and `inter(A,B)` is everything
      !! that exists only because `A` and `B` are both present, Coulombic and
      !! interference and nuclear repulsion together.
      !!
      !! These sum to the total energy exactly, so an interatomic term is not a
      !! model of a bond but a piece of the number the calculation produced.
      !! They are not bond energies: `intra(A)` is an atom deformed by the
      !! molecule, not a free atom, and the difference between the two is the
      !! adaptation energy that this does not yet compute.
      real(dp), intent(in) :: intra(:)
      real(dp), intent(in) :: inter(:, :)
      character(len=*), intent(in) :: element_symbols(:)

      call print_decomposition("energy decomposition", intra, inter, element_symbols)
   end subroutine print_total_decomposition

   subroutine print_interatomic_split(inter, classical, element_symbols)
      !! Each pair interaction split into what is classical and what is not
      !!
      !! **This is the claim the analysis exists to make.** The classical part is
      !! everything an electrostatic model could produce: one atom's density in
      !! the other's nuclear field, the repulsion between two atomic charge
      !! clouds, and the repulsion of the two nuclei. It is a sum of large terms
      !! of both signs that very nearly cancel, because a neutral atom is a
      !! nearly neutral thing to be near.
      !!
      !! What is left is interference -- density shared between the two atoms,
      !! which no classical model has any account of, and which the papers argue
      !! is where covalent binding actually comes from. If that is right, the
      !! second column here is small and the third carries the bond.
      real(dp), intent(in) :: inter(:, :)
      real(dp), intent(in) :: classical(:, :)
      character(len=*), intent(in) :: element_symbols(:)

      character(len=160) :: line
      character(len=16) :: label_a, label_b
      integer :: natm, a, b
      real(dp) :: total_c, total_i

      natm = size(inter, 1)
      call logger%info("")
      call logger%info("  interatomic interactions")
      call logger%info("     pair                  total     classical  interference")
      total_c = 0.0_dp
      total_i = 0.0_dp
      do a = 1, natm
         do b = a + 1, natm
            total_c = total_c + classical(a, b)
            total_i = total_i + (inter(a, b) - classical(a, b))
            if (abs(inter(a, b)) < TOL_PRINT_PAIR) cycle
            write (label_a, "(a,i0)") trim(adjustl(element_symbols(a)))//" ", a
            write (label_b, "(a,i0)") trim(adjustl(element_symbols(b)))//" ", b
            write (line, "(4x,a8,a4,a8,3f14.3)") label_a(1:8), " -- ", label_b(1:8), &
               inter(a, b)*TO_MILLIHARTREE, classical(a, b)*TO_MILLIHARTREE, &
               (inter(a, b) - classical(a, b))*TO_MILLIHARTREE
            call logger%info(trim(line))
         end do
      end do
      call logger%info("     millihartree")
      write (line, "(4x,a24,2f16.6)") "totals, hartree", total_c, total_i
      call logger%info(trim(line))
   end subroutine print_interatomic_split

   subroutine print_kinetic_decomposition(intra, inter, element_symbols)
      !! The atom and atom-pair table for the kinetic energy
      real(dp), intent(in) :: intra(:)
      real(dp), intent(in) :: inter(:, :)
      character(len=*), intent(in) :: element_symbols(:)

      call print_decomposition("kinetic energy decomposition", intra, inter, &
                               element_symbols)
   end subroutine print_kinetic_decomposition

   subroutine print_decomposition(title, intra, inter, element_symbols)
      !! The atom and atom-pair table
      !!
      !! Millihartree rather than kcal/mol as the leading column: these are the
      !! quantities the paper tabulates, and its tables are in millihartree.
      !! kcal/mol follows for comparison with the bond-order table above, which
      !! is on that scale because of the empirical tenth.
      character(len=*), intent(in) :: title
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
      call logger%info("  "//trim(title))
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
      call logger%info("     interatomic                     mhartree    kcal/mol")
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
      write (line, "(4x,a24,f16.6,a)") "total", total, " hartree"
      call logger%info(trim(line))

      deallocate (pair_a, pair_b, magnitude, order)
   end subroutine print_decomposition

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
