!! Resolving the energy into atomic and interatomic contributions
module mqc_libcint_ieda
   !! The intrinsic energy decomposition analysis of Del Angel Cruz, Gordon and
   !! Ruedenberg, JACS 147, 42262 (2025). The quasi-atomic orbitals of
   !! `mqc_libcint_quao` give a basis in which every orbital belongs to one
   !! atom; this module sums the energy over that basis and reports what each
   !! atom contributes on its own and what exists only because the atoms are
   !! bonded.
   !!
   !! Every term of a sum such as `T = sum_pq gamma_pq T_pq` belongs to exactly
   !! one atom or to exactly one pair of atoms, so each decomposition here is a
   !! regrouping and the pieces are obliged to add back up. `kinetic_total`
   !! performs that addition and every routine checks it before returning.
   !!
   !! **The factor of two is the trap.** A pair of atoms `{A,B}` collects both
   !! `p in A, q in B` and `p in B, q in A`, which are equal. `inter` therefore
   !! holds the *full* physical pair energy in both `(A,B)` and `(B,A)`, and
   !! `kinetic_total` halves the matrix rather than summing a triangle.
   !!
   !! Nothing here is scaled, unlike the kinetic bond order of
   !! `kinetic_bond_orders`, which carries an empirical factor of a tenth to
   !! reach the scale of tabulated bond energies.
   use, intrinsic :: iso_fortran_env, only: int64
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_esp, only: esp_matrices
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_aambs, only: aambs_element_counts
   use mqc_determinants, only: generate_strings
   use mqc_libcint_quao, only: quao_result_t
   use mqc_physical_constants, only: HARTREE_TO_KCALMOL
   implicit none
   private

   public :: kinetic_decomposition
   public :: kinetic_total
   public :: project_no_sharing
   public :: combine_quao_sets
   public :: quao_interference
   public :: nuclear_attraction_per_atom
   public :: quao_nuclear_attraction
   public :: nuclear_decomposition
   public :: screened_nucleus_split
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
      !! the only thing this tolerates is accumulated rounding.
   real(dp), parameter :: TOL_SPLIT = 1.0e-8_dp
      !! Hartree, per matrix element. How far the per-nucleus attraction
      !! integrals may drift from the one-shot ones before the split is refused.
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
      !! kinetic energy that exists only because the two atoms share density.
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

      ! Each ordered orbital pair is deposited in both elements of its atom
      ! pair, so `inter(A,B)` holds the energy of the whole pair rather than
      ! half of it.
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

      ! Catches a mis-assigned atom, a dropped orbital, and a mishandled factor
      ! of two on the interatomic block.
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
      !! Add a decomposition back up
      !!
      !! `sum(intra) + 0.5 * sum(inter)`, the only place the halving of `inter`
      !! is written down. Used by every decomposition in this module, not the
      !! kinetic one alone.
      real(dp), intent(in) :: intra(:)
      real(dp), intent(in) :: inter(:, :)
      real(dp) :: total

      total = sum(intra) + 0.5_dp*sum(inter)
   end function kinetic_total

   subroutine project_no_sharing(atom_of, n_atoms, neutral, n_alpha, n_beta, ci, &
                                 recovered, n_kept, error)
      !! Strike out every determinant in which an atom is not neutral
      !!
      !! The no-sharing wave function `Psi-0` of the intrinsic decomposition:
      !! the part of `Psi` in which no charge has moved between atoms. The
      !! energy difference from `Psi` is what charge transfer is worth.
      !!
      !! **`Psi-0` is a projection of `Psi`, not a wave function optimised in
      !! the neutral space.** The amplitudes here are the parent's, with the
      !! rest set to zero and the remainder renormalised; a CI solved inside
      !! the neutral space is a different, lower state.
      !!
      !! **The determinants must be over quasi-atomic orbitals**, since only an
      !! atomic basis can answer how many electrons sit on an atom.
      integer, intent(in) :: atom_of(:)
         !! (n_orbitals) the atom each active orbital belongs to
      integer, intent(in) :: n_atoms
      integer, intent(in) :: neutral(:)
         !! (n_atoms) valence electrons each atom has when neutral
      integer, intent(in) :: n_alpha, n_beta
      real(dp), intent(inout) :: ci(:, :)
         !! (n_alpha_strings, n_beta_strings), overwritten with the projection
      real(dp), intent(out) :: recovered
         !! `|| P Psi ||^2` before renormalising -- a squared norm, not an
         !! overlap
      integer, intent(out) :: n_kept
         !! How many determinants survived
      type(error_t), intent(inout) :: error

      integer(int64), allocatable :: alpha(:), beta(:)
      integer, allocatable :: count_alpha(:, :), count_beta(:, :)
      integer :: n_orb, na, nb, ia, ib, a
      real(dp) :: norm
      logical :: neutral_here

      if (error%has_error()) return

      n_orb = size(atom_of)
      if (n_atoms < 1 .or. maxval(atom_of) > n_atoms .or. size(neutral) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "the orbital-to-atom map and the neutral "// &
                        "electron counts describe different molecules.")
         return
      end if

      call generate_strings(n_orb, n_alpha, alpha, error)
      if (error%has_error()) return
      call generate_strings(n_orb, n_beta, beta, error)
      if (error%has_error()) return
      na = size(alpha)
      nb = size(beta)
      if (size(ci, 1) /= na .or. size(ci, 2) /= nb) then
         call error%set(ERROR_VALIDATION, "the coefficient array is not the shape this "// &
                        "active space produces, so it came from a different space.")
         return
      end if

      ! Per string, not per determinant: the determinant count is the product
      ! of the two string counts.
      allocate (count_alpha(n_atoms, na), count_beta(n_atoms, nb))
      call tally(alpha, atom_of, n_atoms, count_alpha)
      call tally(beta, atom_of, n_atoms, count_beta)

      n_kept = 0
      recovered = 0.0_dp
      do ib = 1, nb
         do ia = 1, na
            neutral_here = .true.
            do a = 1, n_atoms
               if (count_alpha(a, ia) + count_beta(a, ib) /= neutral(a)) then
                  neutral_here = .false.
                  exit
               end if
            end do
            if (neutral_here) then
               n_kept = n_kept + 1
               recovered = recovered + ci(ia, ib)**2
            else
               ci(ia, ib) = 0.0_dp
            end if
         end do
      end do

      if (recovered <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "no determinant leaves every atom neutral, "// &
                        "so there is no no-sharing wave function to project onto. A "// &
                        "closed-shell atom with no valence electrons, or an active "// &
                        "space that is not the whole valence shell, will do this.")
         return
      end if

      norm = 1.0_dp/sqrt(recovered)
      do ib = 1, nb
         do ia = 1, na
            ci(ia, ib) = ci(ia, ib)*norm
         end do
      end do

      deallocate (alpha, beta, count_alpha, count_beta)
   end subroutine project_no_sharing

   pure subroutine tally(strings, atom_of, n_atoms, counts)
      !! How many electrons each string puts on each atom
      integer(int64), intent(in) :: strings(:)
      integer, intent(in) :: atom_of(:)
      integer, intent(in) :: n_atoms
      integer, intent(out) :: counts(:, :)     !! (n_atoms, n_strings)

      integer :: is, p

      counts = 0
      do is = 1, size(strings)
         do p = 1, size(atom_of)
            if (btest(strings(is), p - 1)) then
               counts(atom_of(p), is) = counts(atom_of(p), is) + 1
            end if
         end do
      end do
   end subroutine tally

   subroutine combine_quao_sets(core, valence, combined, error)
      !! Stack the core quasi-atomic orbitals in front of the valence ones
      !!
      !! An energy decomposition needs the core: it carries most of the kinetic
      !! energy and most of the nuclear attraction, where the bonding analysis
      !! can work in the valence-internal space alone.
      !!
      !! **The combined density is exactly two on the core diagonal and zero
      !! elsewhere in the core and core-valence blocks.** Exact, not an
      !! approximation: the core density restricted to the core space is twice
      !! its projector, and the two spaces are orthogonal.
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
      !! Unlike `kinetic_bond_orders` this does not require oriented orbitals.
      !! Orientation is a rotation *within* an atom, so individual orbital pairs
      !! move under one but atom and atom-pair totals do not.
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
      !! Built from `esp_matrices` with the evaluation points on the nuclei.
      !! The pieces must sum back to the ordinary nuclear attraction, which is
      !! checked here and refused on failure: a dropped sign, a missing charge
      !! or coordinates in Angstrom all give an array that looks plausible.
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

      ! The one-shot nuclear attraction: the core Hamiltonian less the kinetic
      ! energy.
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
      !! so the trace against a density is preserved exactly.
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
      !! The last is a choice: a term with `p` on `A`, `q` on `B` and the
      !! nucleus on a third atom is genuinely three-body, and charging it to
      !! `{A,B}` treats the interference as the feature and the field as
      !! context. A finer grouping would redistribute it without changing any
      !! total.
      ! TODO(mqc): the sum-rule failure path deallocates `intra` and `inter` but
      ! not `coulomb`, so an optional argument comes back allocated from a call
      ! that returned an error.
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
         !! another's nuclear field. What is left over is interference --
         !! density shared between the two atoms, which has no classical
         !! description.

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
                  ! Interference: charged to the pair that shares the density.
                  inter(a, b) = inter(a, b) + term
                  inter(b, a) = inter(b, a) + term
               else if (a == c) then
                  intra(a) = intra(a) + term
               else
                  ! One atom's density in another's nuclear field.
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
      !! **Takes the dense atomic-orbital array, so memory is bounded by
      !! `n_ao^4`.** Not a fundamental cost -- the transformation could be
      !! driven from the packed integrals or direct from shell quartets.
      ! TODO(mqc): `cur`, `t`, `step`, `ncol` and `n` are all dead here; the
      ! transformation moved into `four_index_transform` and the declarations
      ! stayed behind.
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
      !! `coeff` is `(n_in, n_out)` and need not be square: the two-electron
      !! integrals come from a large atomic-orbital basis into a small
      !! quasi-atomic one, and the cumulant from a small active space into the
      !! same one.
      real(dp), intent(in) :: coeff(:, :)          !! (n_in, n_out)
      real(dp), intent(in) :: source(:, :, :, :)   !! (n_in, n_in, n_in, n_in)
      real(dp), allocatable, intent(out) :: result(:, :, :, :)

      real(dp), allocatable :: cur(:, :), t(:, :)
      integer :: n_in, n_out, step, ncol

      n_in = size(coeff, 1)
      n_out = size(coeff, 2)

      ! A transposed `coeff` still type-checks, and without this the first
      ! reshape reads past the end of `source`.
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
      !! produce. Zero for a closed-shell reference, and for an MCSCF wave
      !! function zero unless all four indices are active -- the inactive
      !! orbitals are a closed shell and carry no correlation.
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
      !! quasi-atomic orbitals span, and the active orbitals of a converged
      !! MCSCF do not quite.** The deficit is measured and handed back in
      !! `deficit` rather than treated as an error; the caller decides.
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
      !! Exact only where `U` is an isometry, so that the contraction against
      !! the integrals is the same number either side of the transformation.
      !! `quao_projection` measures the departure from that in `deficit`.
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
      !! and `E2 = (1/2) sum_pqrs Gamma_pqrs (pq|rs)`. Coulomb and exchange are
      !! carried together, the assignment following the integral rather than
      !! the density factors.
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
                     ! Two atomic charge clouds repelling each other.
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
      !! in each and `kinetic_total` halves the matrix.
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
      !! Symmetric with a zero diagonal, and in the same convention as every
      !! other pair matrix here -- the whole pair energy in both elements,
      !! halved by `kinetic_total` -- so totals accumulate by adding matrices.
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
      !! `intra(A)` is everything that happens inside atom `A` -- its own
      !! kinetic energy, its electrons in its own nuclear field, its own
      !! electron repulsion -- and `inter(A,B)` is everything that exists only
      !! because `A` and `B` are both present.
      !!
      !! **These are not bond energies.** They sum to the total energy exactly,
      !! but `intra(A)` is an atom deformed by the molecule, not a free atom,
      !! and the adaptation energy between the two is not computed here.
      real(dp), intent(in) :: intra(:)
      real(dp), intent(in) :: inter(:, :)
      character(len=*), intent(in) :: element_symbols(:)

      call print_decomposition("energy decomposition", intra, inter, element_symbols)
   end subroutine print_total_decomposition

   subroutine print_interatomic_split(inter, classical, element_symbols)
      !! Each pair interaction split into what is classical and what is not
      !!
      !! The classical column is everything an electrostatic model could
      !! produce: one atom's density in the other's nuclear field, the repulsion
      !! between two atomic charge clouds, and the repulsion of the two nuclei.
      !! The remainder is interference -- density shared between the two atoms.
      !! Both in millihartree, with the totals in hartree.
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

   subroutine screened_nucleus_split(quao, density, v_atom_quao, n_core_quao, &
                                     atomic_numbers, n_atoms, intra_core, intra_valence, &
                                     error)
      !! Divide an atom's own-nucleus attraction between its core and its valence
      !!
      !! A valence electron does not see the bare nucleus; it sees a nucleus
      !! screened by the core. The papers divide the valence density's attraction
      !! to its own nucleus in the ratio of the charges,
      !!
      !!     V1(Av,Av'/A_core) = (Z_core/Z_A) <Av|-Z_A/r_A|Av'>
      !!     V1(Av,Av'/A_val ) = (Z_val /Z_A) <Av|-Z_A/r_A|Av'>
      !!
      !! and the core density keeps the whole nuclear charge, so the split is
      !! deliberately asymmetric. `Z_core` is twice the chemical-core orbital
      !! count -- 0, 2, 10, 18, ... -- from `aambs_element_counts`.
      !!
      !! **Relabelling, not a screening model.** The two shares sum to the
      !! undivided term because the two charge fractions sum to one, and the
      !! rescaling is applied after the contraction, so it is not
      !! orbital-resolved.
      ! TODO(mqc): unlike every sibling here this validates neither
      ! `size(atomic_numbers)`, `size(v_atom_quao, 3)` nor
      ! `maxval(quao%atom_of)` against `n_atoms`, so an orbital assigned outside
      ! the molecule writes past the end of `intra_core`.
      type(quao_result_t), intent(in) :: quao
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(in) :: v_atom_quao(:, :, :)
      integer, intent(in) :: n_core_quao
         !! How many leading quasi-atomic orbitals are core. `combine_quao_sets`
         !! puts them first, which is the only reason a count suffices.
      integer, intent(in) :: atomic_numbers(:)
      integer, intent(in) :: n_atoms
      real(dp), allocatable, intent(out) :: intra_core(:), intra_valence(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: z_core(:)
      real(dp) :: term, cross, share
      character(len=400) :: message
      integer :: n, i, j, a, core_orbitals, valence_orbitals

      if (error%has_error()) return

      n = quao%n_quao
      if (n_core_quao < 0 .or. n_core_quao > n) then
         call error%set(ERROR_VALIDATION, "the core orbital count is not a prefix of "// &
                        "the quasi-atomic set.")
         return
      end if

      allocate (z_core(n_atoms))
      do a = 1, n_atoms
         call aambs_element_counts(atomic_numbers(a), core_orbitals, valence_orbitals, &
                                   error)
         if (error%has_error()) return
         z_core(a) = 2.0_dp*real(core_orbitals, dp)
      end do

      allocate (intra_core(n_atoms), intra_valence(n_atoms))
      intra_core = 0.0_dp
      intra_valence = 0.0_dp
      cross = 0.0_dp

      do i = 1, n
         a = quao%atom_of(i)
         do j = 1, n
            if (quao%atom_of(j) /= a) cycle
            term = density(i, j)*v_atom_quao(i, j, a)
            if (i <= n_core_quao .and. j <= n_core_quao) then
               intra_core(a) = intra_core(a) + term
            else if (i > n_core_quao .and. j > n_core_quao) then
               share = z_core(a)/real(atomic_numbers(a), dp)
               intra_core(a) = intra_core(a) + share*term
               intra_valence(a) = intra_valence(a) + (1.0_dp - share)*term
            else
               ! Core-valence: zero for a frozen core, since the two spaces are
               ! orthogonal. Accumulated and checked below rather than assumed.
               cross = cross + abs(term)
               intra_core(a) = intra_core(a) + term
            end if
         end do
      end do

      if (cross > TOL_SUM_RULE) then
         write (message, "(a,es12.4,a)") "the density connects core and valence "// &
            "quasi-atomic orbitals, by ", cross, " hartree. The screened-nucleus "// &
            "split assumes it does not, which holds for a frozen core and fails "// &
            "for a correlated one."
         call error%set(ERROR_VALIDATION, trim(message))
         deallocate (intra_core, intra_valence, z_core)
         return
      end if

      deallocate (z_core)
   end subroutine screened_nucleus_split

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
      !! Millihartree in the leading column, kcal/mol beside it, and the three
      !! totals in hartree. Pairs are ordered by descending magnitude.
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
      ! TODO(mqc): a hand-written selection sort, quadratic in the number of
      ! atom pairs and so quartic in the atom count, where `pic_sorting`'s
      ! `sort_index` is already used elsewhere in the project for exactly this.
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
