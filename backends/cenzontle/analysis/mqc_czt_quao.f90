!! Quasi-atomic orbitals: the valence-virtual space
module mqc_czt_quao
   !! The first step of the Ruedenberg construction that needs integrals.
   !!
   !! The valence-virtual orbitals are the chemically meaningful part of a
   !! large-basis virtual space: the antibonding counterparts of the occupied
   !! valence orbitals, recovered by asking which combinations of virtuals look
   !! most like free-atom minimal-basis orbitals. Unlike the ordinary virtual
   !! spectrum they barely depend on the basis.
   !!
   !! Reference: West, Schmidt, Gordon, Ruedenberg, J. Chem. Phys. 139, 234107
   !! (2013), section V.A.1 and Appendix A.1.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cgto, only: molecular_basis_type
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use mqc_aambs, only: aambs_file, aambs_dimensions, aambs_dimensions_t, &
                        aambs_element_counts
   use mqc_czt_integrals, only: czt_molecule_t, mixed_basis_overlap, &
                                atom_ao_blocks
   use mqc_physical_constants, only: HARTREE_TO_KCALMOL
   implicit none
   private

   public :: build_aambs_molecule    !! The free-atom minimal basis on this geometry
   public :: mo_aambs_overlap        !! < MO | orthogonalized AAMBS >
   public :: valence_virtual_orbitals
   public :: vvo_result_t
   public :: aambs_atom_ranges       !! Where each atom's minimal-basis orbitals sit
   public :: quasi_atomic_orbitals
   public :: orient_quasi_atomic_orbitals
   public :: kinetic_bond_orders
   public :: split_localize
   public :: quao_result_t

   real(dp), parameter :: ORIENT_ANGLE_TOL = 1.0e-6_dp
      !! Rotations smaller than this, in radians, are not applied. A skipped
      !! rotation moves a bond order by at most this much times its size,
      !! on the scale of the 1e-6 agreement with reference values the tests ask of
      !! bond orders, and below the four decimals the bond table prints.
   real(dp), parameter :: ORIENT_FLAT = 1.0e-8_dp
      !! Amplitude of a pair's angular dependence below which the pair is left
      !! alone: every angle then gives the same functional, and the angle
      !! atan2 would return is a ratio of round-off.
   integer, parameter :: ORIENT_MAX_SWEEPS = 2000
      !! Default sweep limit for the orientation and the split localization
   real(dp), parameter :: ORIENT_GAIN_CRIT = 1.0e-9_dp
      !! Largest single-rotation increase of eq (5.5) below which an orientation
      !! that reaches its sweep limit is accepted rather than refused.
   real(dp), parameter :: ORIENT_SENSE_SLACK = 1.0e-10_dp
      !! How far eq (5.5) may fall over the whole orientation, as round-off,
      !! before the rotations are judged to have been applied backwards
   real(dp), parameter :: KBO_SCALE = 0.1_dp
      !! Paper II eq (2), and empirical: raw kinetic interference energies in
      !! organic molecules run about an order of magnitude above the bond
      !! energies chemistry quotes, so they are scaled to be comparable.

   real(dp), parameter :: REFINE_ANGLE_TOL = 1.0e-10_dp
      !! Rotations smaller than this, in radians, are not applied, and a sweep
      !! that applies none ends the refinement. Tight because the refinement
      !! fixes the orbitals every later step works in: a residual angle moves
      !! populations and bond orders to first order, while the functional
      !! itself only moves to second.
   real(dp), parameter :: REFINE_FLAT = 1.0e-10_dp
      !! Size of `G_A - G_B`'s traceless part below which every angle gives the
      !! same functional and the pair is left alone
   integer, parameter :: REFINE_MAX_SWEEPS = 5000
      !! The refinement normally settles in a few tens of sweeps. The cap only
      !! stops a run that cannot settle, and is set high enough that one
      !! converging linearly at 0.995 per sweep from an angle of 0.1 still
      !! reaches `REFINE_ANGLE_TOL` inside it.

   type :: vvo_result_t
      !! What the valence-virtual extraction leaves behind
      real(dp), allocatable :: orbitals(:, :)
         !! (n_ao, n_vvo), in the AO basis, orthonormal and orthogonal to the
         !! occupied space because they are combinations of virtuals only
      real(dp), allocatable :: external_orbitals(:, :)
         !! (n_ao, n_virt - n_vvo), the rest of the virtual space. Together with
         !! `orbitals` a complete, orthonormal virtual block, which is what a
         !! calculation running in the valence space rather than analysing it
         !! has to be handed.
      real(dp), allocatable :: singular_values(:)
         !! All of them, descending. The first `n_vvo` are retained; the rest
         !! span the valence-external space.
      integer :: n_vvo = 0
      real(dp) :: smallest_retained = 0.0_dp
      real(dp) :: largest_rejected = 0.0_dp
         !! The gap between these two is the diagnostic: a retained value near
         !! one against a small rejected one. Anything else means the projection
         !! is not finding a clean valence space.
   end type vvo_result_t

   type :: quao_result_t
      !! Quasi-atomic orbitals and the density in their basis
      real(dp), allocatable :: orbitals(:, :)
         !! (n_ao, n_quao), orthonormal, each one belonging to a single atom
      integer, allocatable :: atom_of(:)
         !! Which atom each quasi-atomic orbital sits on
      integer, allocatable :: subspace_of(:)
         !! Which occupation-restricted subspace each one was drawn from, all
         !! ones when the space is unrestricted. No rotation anywhere in this
         !! module may mix two subspaces: a restricted wave function is not
         !! invariant under one that does.
      real(dp), allocatable :: population_bond_order(:, :)
         !! (n_quao, n_quao). Diagonal elements are orbital populations;
         !! interatomic off-diagonal elements are bond orders. Paper I eq (5.2).
      integer :: n_quao = 0
      integer :: sweeps = 0            !! Jacobi sweeps the refinement took
      real(dp), allocatable :: to_valence_internal(:, :)
         !! (n_val, n_quao). Column `i` expands quasi-atomic orbital `i` in the
         !! valence-internal orbitals, so its rows are indexed by those. The
         !! space the split-localization works in.
      logical :: oriented = .false.
      real(dp) :: orientation_sum = 0.0_dp
         !! Paper I eq (5.5) after orientation
      logical :: orientation_stalled = .false.
         !! The orientation reached its sweep limit with angles still above
         !! `ORIENT_ANGLE_TOL` but every rotation gaining under `ORIENT_GAIN_CRIT`,
         !! and was accepted
      real(dp) :: orientation_gain = 0.0_dp
         !! Largest single-rotation increase of eq (5.5) in the final sweep
      real(dp) :: orientation_angle = 0.0_dp
         !! Largest rotation angle in the final sweep, radians; zero when the
         !! orientation settled
      real(dp), allocatable :: character_of(:)
         !! (n_quao) how much of each orbital lies inside its own atom's
         !! free-atom space, as a norm rather than a square, so it is on the
         !! scale the papers quote overlaps on. `atomic_character` is the mean
         !! of the squares of these.
      real(dp) :: atomic_character = 0.0_dp
         !! The maximized functional, Paper II eq (A.5), divided by the number of
         !! orbitals. One means every orbital lies entirely within its own atom's
         !! free-atom space; the deficit is the deformation the molecule imposes.
   end type quao_result_t

contains

   subroutine build_aambs_molecule(atomic_numbers, element_symbols, coordinates, mol, error)
      !! The accurate atomic minimal basis, on this molecule's nuclei
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)   !! (3, natm), Bohr
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: error

      type(molecular_basis_type) :: basis
      character(len=:), allocatable :: path

      if (error%has_error()) return
      call aambs_file(path, error)
      if (error%has_error()) return
      call build_molecular_basis_json(path, element_symbols, basis, error)
      if (error%has_error()) return
      call mol%build(atomic_numbers, coordinates, basis, error)
      call basis%destroy()
   end subroutine build_aambs_molecule

   subroutine inverse_sqrt(matrix, result, error, floor)
      !! `A^(-1/2)` for a symmetric positive definite matrix
      real(dp), intent(in) :: matrix(:, :)
      real(dp), allocatable, intent(out) :: result(:, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: floor

      real(dp), allocatable :: vectors(:, :), values(:), scaled(:, :)
      real(dp) :: small
      integer :: n, i, info

      if (error%has_error()) return
      small = 1.0e-10_dp
      if (present(floor)) small = floor

      n = size(matrix, 1)
      allocate (vectors(n, n), values(n), scaled(n, n), result(n, n))
      vectors = matrix
      call pic_syev(vectors, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "the minimal-basis overlap could not be "// &
                        "diagonalized (info = "//to_char(info)//")")
         return
      end if

      ! A minimal basis on well-separated atoms is far from linearly dependent,
      ! so a small eigenvalue here is a statement about the geometry rather than
      ! a threshold to be tuned.
      if (values(1) < small) then
         call error%set(ERROR_VALIDATION, "the free-atom minimal basis is linearly "// &
                        "dependent on this geometry (smallest overlap eigenvalue "// &
                        to_char(values(1))//"). Two atoms are on top of each other, or "// &
                        "close enough that their free-atom orbitals cannot be told apart.")
         return
      end if

      do i = 1, n
         scaled(:, i) = vectors(:, i)/sqrt(values(i))
      end do
      call pic_gemm(scaled, vectors, result, transb="T")
   end subroutine inverse_sqrt

   subroutine mo_aambs_overlap(orbitals, mixed, aambs_overlap, projection, error, orthogonalize)
      !! `< MO_p | A#a >`, the molecular orbitals against the orthogonalized AAMBS
      !!
      !! **Both bases have to be orthonormal.** Paper I's Appendix A.1 proves
      !! the singular vectors maximize the projection of one space into the
      !! other only under that condition. The molecular orbitals are orthonormal
      !! already; the free-atom orbitals are orthonormal within one atom but not
      !! between atoms, so they are symmetrically orthogonalized across the
      !! molecule first -- Paper I's `|A#a>`.
      real(dp), intent(in) :: orbitals(:, :)        !! C, (n_ao, n_mo)
      real(dp), intent(in) :: mixed(:, :)           !! < AO | AAMBS >, (n_ao, n_mbs)
      real(dp), intent(in) :: aambs_overlap(:, :)   !! < AAMBS | AAMBS >, (n_mbs, n_mbs)
      real(dp), allocatable, intent(out) :: projection(:, :)   !! (n_mo, n_mbs)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: orthogonalize
         !! Default true, giving `|A#a>`. False gives the raw `|A*a>`, which is
         !! what the *per-atom* decompositions want: within one atom the metric
         !! is already the identity, and orthogonalizing across the molecule
         !! would mix in the neighbours the projection is distinguishing from.

      real(dp), allocatable :: half(:, :), work(:, :)
      integer :: n_mo, n_mbs
      logical :: rotate

      if (error%has_error()) return
      rotate = .true.
      if (present(orthogonalize)) rotate = orthogonalize
      n_mo = size(orbitals, 2)
      n_mbs = size(mixed, 2)

      allocate (work(n_mo, n_mbs))
      call pic_gemm(orbitals, mixed, work, transa="T")

      if (.not. rotate) then
         call move_alloc(work, projection)
         return
      end if

      call inverse_sqrt(aambs_overlap, half, error)
      if (error%has_error()) return
      allocate (projection(n_mo, n_mbs))
      call pic_gemm(work, half, projection)
      deallocate (work, half)
   end subroutine mo_aambs_overlap

   subroutine aambs_atom_ranges(atomic_numbers, aambs, core_offset, core_count, &
                                valence_offset, valence_count, error)
      !! Where each atom's core and valence minimal-basis orbitals sit
      !!
      !! Each atom's functions are one contiguous run, and *within* that run the
      !! chemical core comes first. That second fact is the invariant the whole
      !! core/valence split rests on, and it is why the shell ordering in
      !! `aambs.json` is not simply by principal quantum number: scandium
      !! through zinc list 4s before 3d so that the argon core stays a prefix.
      integer, intent(in) :: atomic_numbers(:)
      type(czt_molecule_t), intent(in) :: aambs
      integer, allocatable, intent(out) :: core_offset(:), core_count(:)
      integer, allocatable, intent(out) :: valence_offset(:), valence_count(:)
      type(error_t), intent(inout) :: error

      integer, allocatable :: offsets(:), counts(:)
      integer :: natm, iatom, core, valence

      if (error%has_error()) return
      natm = size(atomic_numbers)
      allocate (offsets(natm), counts(natm))
      call atom_ao_blocks(aambs, offsets, counts)

      allocate (core_offset(natm), core_count(natm))
      allocate (valence_offset(natm), valence_count(natm))
      do iatom = 1, natm
         call aambs_element_counts(atomic_numbers(iatom), core, valence, error)
         if (error%has_error()) return
         if (core + valence /= counts(iatom)) then
            call error%set(ERROR_VALIDATION, "atom "//to_char(iatom)//" has "// &
                           to_char(counts(iatom))//" minimal-basis functions but its "// &
                           "counts say "//to_char(core + valence)//". The shell data "// &
                           "and the orbital counts in aambs.json disagree.")
            return
         end if
         core_offset(iatom) = offsets(iatom)
         core_count(iatom) = core
         valence_offset(iatom) = offsets(iatom) + core
         valence_count(iatom) = valence
      end do
   end subroutine aambs_atom_ranges

   subroutine select_atomic_claims(projection, offset, count, natm, first, last, &
                                   claims, atom_of, subspace_of, error, starved)
      !! Which atom gets which orbital, by competition rather than by decree
      !!
      !! Every atom is offered every subspace. The claims are pooled, ranked by
      !! how much of the subspace each one actually captures, and the strongest
      !! are kept until the subspace is full. **An atom's share is therefore an
      !! output.** It is not read off the free-atom minimal basis, and it can
      !! differ from that atom's free-atom valence count. Only the molecule-wide
      !! sum is checked against the free-atom counts.
      !!
      !! **With one subspace this reproduces the prescribed per-atom rule
      !! exactly.** The claims then sum to the valence dimension by
      !! construction, so the pool fills the space and the ranking selects
      !! everything. The ranking is for a restricted space, where the subspaces
      !! are smaller than the whole and the atoms genuinely compete.
      !!
      !! Each claim's coefficients are confined to the rows of its own subspace,
      !! so two claims from different subspaces have disjoint support and are
      !! orthogonal whatever else happens to them.
      real(dp), intent(in) :: projection(:, :)   !! (n_val, n_mbs)
      integer, intent(in) :: offset(:), count(:)
      integer, intent(in) :: natm
      integer, intent(in) :: first(:), last(:)   !! (n_subspaces) row ranges
      real(dp), allocatable, intent(out) :: claims(:, :)     !! (n_val, n_val)
      integer, allocatable, intent(out) :: atom_of(:), subspace_of(:)
      type(error_t), intent(inout) :: error
      logical, intent(out), optional :: starved
         !! Set when the partition leaves some atom with no orbital, instead of
         !! raising an error. A property of the partition rather than a fault:
         !! occupation-restricted subspaces are usually grouped by orbital
         !! energy, so one strong atom can win every slot. Absent, the same
         !! situation raises an error.

      real(dp), allocatable :: pool(:, :), weight(:), sigma(:, :), block(:, :)
      real(dp), allocatable :: piece(:)
      integer, allocatable :: pool_atom(:), pool_sub(:), taken(:), share(:)
      integer :: n_val, n_sub, k, iatom, m, width, i, j, n_cand, filled, best, kept
      real(dp) :: top

      if (error%has_error()) return
      n_val = size(projection, 1)
      n_sub = size(first)

      n_cand = 0
      do k = 1, n_sub
         width = last(k) - first(k) + 1
         do iatom = 1, natm
            n_cand = n_cand + min(count(iatom), width)
         end do
      end do

      allocate (pool(n_val, n_cand), weight(n_cand))
      allocate (pool_atom(n_cand), pool_sub(n_cand))
      pool = 0.0_dp
      filled = 0
      do k = 1, n_sub
         width = last(k) - first(k) + 1
         do iatom = 1, natm
            m = min(count(iatom), width)
            if (m == 0 .or. count(iatom) == 0) cycle
            allocate (sigma(width, count(iatom)))
            sigma = projection(first(k):last(k), &
                               offset(iatom) + 1:offset(iatom) + count(iatom))
            call atom_adapted_block(sigma, block, error, weights=piece, n_wanted=m)
            deallocate (sigma)
            if (error%has_error()) return
            do i = 1, m
               filled = filled + 1
               ! Confined to this subspace's rows, zero elsewhere.
               pool(first(k):last(k), filled) = block(:, i)
               weight(filled) = piece(i)
               pool_atom(filled) = iatom
               pool_sub(filled) = k
            end do
            deallocate (block, piece)
         end do
      end do

      ! Keep the strongest claims on each subspace until it is full.
      allocate (taken(n_cand), share(natm))
      taken = 0
      share = 0
      do k = 1, n_sub
         width = last(k) - first(k) + 1
         do kept = 1, width
            best = 0
            top = -1.0_dp
            do i = 1, n_cand
               if (taken(i) /= 0 .or. pool_sub(i) /= k) cycle
               if (weight(i) > top) then
                  top = weight(i)
                  best = i
               end if
            end do
            if (best == 0) then
               call error%set(ERROR_VALIDATION, "subspace "//to_char(k)//" has "// &
                              to_char(width)//" orbitals but the atoms between them "// &
                              "claim fewer, so it cannot be filled.")
               return
            end if
            taken(best) = 1
            share(pool_atom(best)) = share(pool_atom(best)) + 1
         end do
      end do

      ! Every atom that had something to claim must end up with something: an
      ! atom with no quasi-atomic orbital has no population, no bond order and
      ! no place in the decomposition. Atoms with nothing to claim are skipped,
      ! since this set is built for the core orbitals as well as the valence
      ! ones and hydrogen has no core.
      if (present(starved)) starved = .false.
      do iatom = 1, natm
         if (count(iatom) == 0) cycle
         if (share(iatom) < 1) then
            if (present(starved)) then
               starved = .true.
               return
            end if
            call error%set(ERROR_VALIDATION, "atom "//to_char(iatom)//" was left "// &
                           "with no quasi-atomic orbital, so it has no population "// &
                           "and no place in the analysis. The subspaces are too "// &
                           "small, or too few, to give every atom one.")
            return
         end if
      end do

      ! Laid out subspace by subspace and, inside each, atom by atom, so with
      ! one subspace the ordering is atom-major.
      allocate (claims(n_val, n_val), atom_of(n_val), subspace_of(n_val))
      j = 0
      do k = 1, n_sub
         do iatom = 1, natm
            do i = 1, n_cand
               if (taken(i) == 0) cycle
               if (pool_sub(i) /= k .or. pool_atom(i) /= iatom) cycle
               j = j + 1
               claims(:, j) = pool(:, i)
               atom_of(j) = iatom
               subspace_of(j) = k
            end do
         end do
      end do
      if (j /= n_val) then
         call error%set(ERROR_VALIDATION, "the subspaces between them hold "// &
                        to_char(j)//" orbitals and the valence space has "// &
                        to_char(n_val)//".")
         return
      end if

      deallocate (pool, weight, pool_atom, pool_sub, taken, share)
   end subroutine select_atomic_claims

   subroutine atom_adapted_block(sigma, block, error, weights, n_wanted)
      !! The orbitals of one subspace that best match one atom's free-atom set
      !!
      !! `sigma` is the overlap of an orthonormal set of molecular orbitals with
      !! one atom's free-atom orbitals, so its left singular vectors are the
      !! combinations of those orbitals lying closest to that atom -- Paper I
      !! section V.B.1. Taken through the small Gram matrix `sigma^T sigma`,
      !! whose order is that atom's minimal-basis dimension.
      real(dp), intent(in) :: sigma(:, :)     !! (n_rows, m_atom)
      real(dp), allocatable, intent(out) :: block(:, :)   !! (n_rows, n_wanted)
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: weights(:)
         !! The squared projections of the vectors returned, descending. What an
         !! atom's claim on a subspace is ranked by when several atoms compete.
      integer, intent(in), optional :: n_wanted
         !! How many vectors to return, default all of them. Fewer when the
         !! subspace is smaller than the atom's minimal basis and there are not
         !! that many dimensions to give.

      real(dp), allocatable :: gram(:, :), values(:), work(:, :)
      integer :: m, i, info, taken

      if (error%has_error()) return
      m = size(sigma, 2)
      taken = m
      if (present(n_wanted)) taken = min(n_wanted, m)
      if (taken < 1) then
         allocate (block(size(sigma, 1), 0))
         if (present(weights)) allocate (weights(0))
         return
      end if
      allocate (gram(m, m), values(m))
      call pic_gemm(sigma, sigma, gram, transa="T")

      ! Negated so the largest projections come first, matching the ordering the
      ! selection downstream assumes.
      gram = -gram
      call pic_syev(gram, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "an atomic projection could not be "// &
                        "diagonalized (info = "//to_char(info)//")")
         return
      end if
      values = -values

      ! Only the vectors actually being taken have to be representable: where
      ! the subspace is smaller than the atom's minimal basis the trailing
      ! projections are legitimately zero.
      if (values(taken) < 1.0e-10_dp) then
         call error%set(ERROR_VALIDATION, "one of this atom's free-atom orbitals has "// &
                        "no counterpart in the molecular orbital space (projection "// &
                        to_char(values(taken))//"). The orbital basis cannot represent "// &
                        "the minimal basis, which usually means it is far too small.")
         return
      end if

      allocate (work(size(sigma, 1), m), block(size(sigma, 1), taken))
      call pic_gemm(sigma, gram, work)
      do i = 1, taken
         block(:, i) = work(:, i)/sqrt(values(i))
      end do
      if (present(weights)) weights = values(1:taken)
      deallocate (gram, values, work)
   end subroutine atom_adapted_block

   subroutine refine_atomic_character(coefficients, projection, atom_of, subspace_of, &
                                      offset, count, sweeps, functional, error, &
                                      character_of)
      !! Rotate between atoms so each orbital sits as fully as possible on its own
      !!
      !! Paper II's Appendix, eqs (A.1)-(A.11), which replaces the plain
      !! symmetric orthogonalization of Paper I. It maximizes
      !!
      !!     P = sum_A sum_{a on A} sum_{alpha on A} <A*alpha | Aa>^2
      !!
      !! over orthogonal transformations, by 2x2 rotations between orbitals on
      !! *different* atoms. A rotation inside one atom leaves P unchanged, so
      !! those pairs are never visited.
      !!
      !! For orbitals `i` on atom A and `j` on atom B, let `G_A` and `G_B` be the
      !! 2x2 Gram matrices of their projections onto A's and B's free-atom
      !! orbitals. After rotating to `(c i + s j, -s i + c j)` the pair's share
      !! of P depends on the unit vector `e = (c, s)` as
      !!
      !!     P_ij(e) = e^T G_A e + (I - e e^T) : G_B = tr G_B + e^T (G_A - G_B) e
      !!
      !! so the best rotation turns `i` onto the leading eigenvector of
      !! `M = G_A - G_B`, which is `2 theta = atan2(M_12, (M_11 - M_22)/2)`.
      real(dp), intent(inout) :: coefficients(:, :)   !! (n_rows, n_quao)
      real(dp), intent(in) :: projection(:, :)        !! (n_rows, n_mbs), raw AAMBS
      integer, intent(in) :: atom_of(:)
      integer, intent(in) :: subspace_of(:)
         !! Rotations across two subspaces are skipped: a restricted wave
         !! function is not invariant under them, and one would destroy the
         !! disjoint support that makes the transformation block-diagonal. With
         !! a single subspace the guard never fires.
      integer, intent(in) :: offset(:), count(:)      !! AAMBS range per atom
      integer, intent(out) :: sweeps
      real(dp), intent(out) :: functional
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: character_of(:)
         !! Per orbital, the norm of its projection onto its own atom's
         !! free-atom space -- the quantity the papers quote per QUAO

      real(dp), allocatable :: overlap(:, :), own_share(:)
      integer, allocatable :: pairs(:, :)
      real(dp) :: gram_own(3), gram_other(3), mismatch(3)
      real(dp) :: half_gap, radius, theta
      integer :: n_orb, sweep, ipair, first, second
      logical :: rotated

      if (error%has_error()) return
      n_orb = size(coefficients, 2)
      call rotatable_pairs(atom_of, subspace_of, .false., pairs)
      allocate (overlap(n_orb, size(projection, 2)))

      sweeps = 0
      rotated = .true.
      do sweep = 1, REFINE_MAX_SWEEPS
         ! overlap(a, alpha) = < orbital a | free-atom orbital alpha >, taken
         ! afresh at the top of each sweep and kept current within it by the
         ! same rotations the coefficients receive.
         call pic_gemm(coefficients, projection, overlap, transa="T")
         rotated = .false.
         do ipair = 1, size(pairs, 2)
            first = pairs(1, ipair)
            second = pairs(2, ipair)
            ! Packed as (11, 22, 12).
            gram_own = pair_gram(overlap, first, second, offset(atom_of(first)), &
                                 count(atom_of(first)))
            gram_other = pair_gram(overlap, first, second, offset(atom_of(second)), &
                                   count(atom_of(second)))
            mismatch = gram_own - gram_other
            half_gap = 0.5_dp*(mismatch(1) - mismatch(2))
            radius = hypot(half_gap, mismatch(3))
            ! A vanishing radius means M is a multiple of the identity and every
            ! angle gives the same P: two orbitals that are symmetry images of
            ! each other land here, and atan2 of two round-off values would
            ! otherwise rotate them at random.
            if (radius < REFINE_FLAT) cycle
            theta = 0.5_dp*atan2(mismatch(3), half_gap)
            if (abs(theta) < REFINE_ANGLE_TOL) cycle
            call rotate_columns(coefficients, first, second, cos(theta), sin(theta))
            call rotate_rows(overlap, first, second, cos(theta), sin(theta))
            rotated = .true.
         end do
         sweeps = sweep
         if (.not. rotated) exit
      end do

      if (rotated) then
         call error%set(ERROR_VALIDATION, "the quasi-atomic refinement did not settle "// &
                        "in "//to_char(REFINE_MAX_SWEEPS)//" sweeps.")
         return
      end if

      ! Each orbital's squared projection onto its own atom's free-atom space.
      call pic_gemm(coefficients, projection, overlap, transa="T")
      allocate (own_share(n_orb))
      do first = 1, n_orb
         associate (lo => offset(atom_of(first)) + 1, &
                    hi => offset(atom_of(first)) + count(atom_of(first)))
            own_share(first) = sum(overlap(first, lo:hi)**2)
         end associate
      end do
      functional = sum(own_share)/real(n_orb, dp)
      if (present(character_of)) character_of = sqrt(own_share)
   end subroutine refine_atomic_character

   pure function pair_gram(overlap, first, second, start, length) result(gram)
      !! Gram matrix of two orbitals' projections onto one atom's free-atom block
      !!
      !! Returned packed as `(G_11, G_22, G_12)`. The block is the `length`
      !! columns of `overlap` after column `start`.
      real(dp), intent(in) :: overlap(:, :)   !! (n_orb, n_mbs)
      integer, intent(in) :: first, second, start, length
      real(dp) :: gram(3)

      associate (u => overlap(first, start + 1:start + length), &
                 v => overlap(second, start + 1:start + length))
         gram(1) = dot_product(u, u)
         gram(2) = dot_product(v, v)
         gram(3) = dot_product(u, v)
      end associate
   end function pair_gram

   subroutine rotatable_pairs(atom_of, subspace_of, same_atom, pairs)
      !! The orbital pairs a Jacobi sweep visits, in the order it visits them
      !!
      !! Every `(a, b)` with `a < b` in the same subspace, and either on the
      !! same atom (`same_atom`) or on different atoms (not `same_atom`).
      !! Fixed for the whole sweep sequence, since neither map changes.
      integer, intent(in) :: atom_of(:), subspace_of(:)
      logical, intent(in) :: same_atom
      integer, allocatable, intent(out) :: pairs(:, :)   !! (2, n_pairs)

      integer, allocatable :: buffer(:, :)
      integer :: n_orb, a, b, n_pairs

      n_orb = size(atom_of)
      allocate (buffer(2, max(1, n_orb*(n_orb - 1)/2)))
      n_pairs = 0
      do a = 1, n_orb
         do b = a + 1, n_orb
            if (subspace_of(a) /= subspace_of(b)) cycle
            if ((atom_of(a) == atom_of(b)) .neqv. same_atom) cycle
            n_pairs = n_pairs + 1
            buffer(:, n_pairs) = [a, b]
         end do
      end do
      pairs = buffer(:, 1:n_pairs)
   end subroutine rotatable_pairs

   pure subroutine rotate_columns(matrix, first, second, c, s)
      !! Columns `first, second` <- `(c first + s second, -s first + c second)`
      real(dp), intent(inout) :: matrix(:, :)
      integer, intent(in) :: first, second
      real(dp), intent(in) :: c, s

      real(dp) :: old_first(size(matrix, 1))

      old_first = matrix(:, first)
      matrix(:, first) = c*old_first + s*matrix(:, second)
      matrix(:, second) = -s*old_first + c*matrix(:, second)
   end subroutine rotate_columns

   pure subroutine rotate_rows(matrix, first, second, c, s)
      !! Rows `first, second` <- `(c first + s second, -s first + c second)`
      real(dp), intent(inout) :: matrix(:, :)
      integer, intent(in) :: first, second
      real(dp), intent(in) :: c, s

      real(dp) :: old_first(size(matrix, 2))

      old_first = matrix(first, :)
      matrix(first, :) = c*old_first + s*matrix(second, :)
      matrix(second, :) = -s*old_first + c*matrix(second, :)
   end subroutine rotate_rows

   subroutine valence_virtual_orbitals(orbitals, projection, dims, result, error)
      !! Recover the chemically meaningful part of the virtual space
      !!
      !! The number kept is fixed by counting, not by a threshold:
      !! `n_vvo = n_mbs - n_occ` exactly. The singular values are reported so a
      !! caller can check there is a gap where the cut falls, but they do not
      !! decide where it falls.
      !!
      !! Solved as the eigenproblem of `Sigma Sigma^T` rather than by an SVD:
      !! Appendix A.1 shows these have the same left singular vectors, and the
      !! external block of the rectangular problem is full of degenerate zero
      !! singular values whose vectors are arbitrary.
      real(dp), intent(in) :: orbitals(:, :)     !! C, (n_ao, n_mo)
      real(dp), intent(in) :: projection(:, :)   !! < MO | A#a >, (n_mo, n_mbs)
      type(aambs_dimensions_t), intent(in) :: dims
      type(vvo_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: sigma(:, :), b(:, :), values(:), rotated(:, :)
      real(dp), allocatable :: c_virt(:, :)
      integer :: n_ao, n_mo, n_virt, n_occ, i, info

      if (error%has_error()) return
      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_occ = dims%n_occupied
      n_virt = n_mo - n_occ

      if (dims%n_vvo > n_virt) then
         call error%set(ERROR_VALIDATION, "the minimal basis asks for "// &
                        to_char(dims%n_vvo)//" valence-virtual orbitals but there are "// &
                        "only "//to_char(n_virt)//" virtuals. The orbital basis is "// &
                        "smaller than the minimal basis, which cannot describe the "// &
                        "valence space at all.")
         return
      end if

      ! The virtual block of the projection, and its Gram matrix.
      allocate (sigma(n_virt, size(projection, 2)))
      sigma = projection(n_occ + 1:n_mo, :)
      allocate (b(n_virt, n_virt), values(n_virt))
      call pic_gemm(sigma, sigma, b, transb="T")

      ! Negated so the largest projections come first once the sign is undone,
      ! which keeps the retained block at the front where every downstream index
      ! expects it.
      b = -b
      call pic_syev(b, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "the valence-virtual eigenproblem failed "// &
                        "(info = "//to_char(info)//")")
         return
      end if
      values = -values

      allocate (result%singular_values(n_virt))
      do i = 1, n_virt
         ! Negative only by rounding, at eigenvalues that are zero to begin with.
         result%singular_values(i) = sqrt(max(values(i), 0.0_dp))
      end do

      result%n_vvo = dims%n_vvo
      result%smallest_retained = result%singular_values(max(dims%n_vvo, 1))
      if (dims%n_vvo < n_virt) then
         result%largest_rejected = result%singular_values(dims%n_vvo + 1)
      else
         result%largest_rejected = 0.0_dp
      end if

      allocate (c_virt(n_ao, n_virt), rotated(n_ao, dims%n_vvo))
      c_virt = orbitals(:, n_occ + 1:n_mo)
      call pic_gemm(c_virt, b(:, 1:dims%n_vvo), rotated)
      call move_alloc(rotated, result%orbitals)

      ! The same rotation applied to the columns the count rejected: the
      ! valence-external space, the other half of a complete virtual block.
      if (dims%n_vvo < n_virt) then
         allocate (rotated(n_ao, n_virt - dims%n_vvo))
         call pic_gemm(c_virt, b(:, dims%n_vvo + 1:n_virt), rotated)
         call move_alloc(rotated, result%external_orbitals)
      else
         allocate (result%external_orbitals(n_ao, 0))
      end if

      deallocate (sigma, b, values, c_virt)
   end subroutine valence_virtual_orbitals

   subroutine quasi_atomic_orbitals(atomic_numbers, valence_internal, n_occupied_valence, &
                                    mixed, offset, count, result, error, &
                                    valence_density, subspaces, starved)
      !! Quasi-atomic orbitals for the valence-internal space
      !!
      !! Each atom claims the combinations of valence-internal orbitals that
      !! look most like its own free-atom orbitals. Those claims are orthonormal
      !! within an atom and not between atoms, so the collection is
      !! symmetrically orthogonalized and then refined by the Paper II sweep,
      !! which recovers the atomic character the orthogonalization costs.
      !!
      !! The result is a basis for the same space the molecular orbitals span,
      !! in which every orbital belongs to one atom. The density in that basis
      !! is the population-bond-order matrix: its diagonal counts electrons on
      !! an orbital, and its interatomic off-diagonal elements are bond orders.
      !!
      !! **Bond-order signs are meaningless on their own.** Nothing here fixes a
      !! phase -- the eigenvectors come back however LAPACK produced them -- and
      !! the papers say so explicitly: a negative bond order routinely describes
      !! a bonding interaction. Compare magnitudes, or use the kinetic bond
      !! order, in which the two phase factors cancel.
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: valence_internal(:, :)
         !! (n_ao, n_val): the occupied valence orbitals followed by the
         !! valence-virtual ones, orthonormal
      integer, intent(in) :: n_occupied_valence
         !! How many of those columns are occupied. The rest are empty. Used to
         !! build the closed-shell density when `valence_density` is absent, and
         !! for nothing else.
      real(dp), intent(in) :: mixed(:, :)      !! < AO | AAMBS >, (n_ao, n_mbs)
      integer, intent(in) :: offset(:), count(:)   !! Valence AAMBS range per atom
      type(quao_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      logical, intent(out), optional :: starved
         !! Reports a partition that cannot give every atom an orbital, rather
         !! than failing on it. Only meaningful beside `subspaces`.
      integer, intent(in), optional :: subspaces(:)
         !! The valence-internal orbital each occupation-restricted subspace
         !! starts at, ascending, as `ormas_space_t%first_orbital` gives them.
         !! Absent means one subspace covering everything, which reproduces the
         !! prescribed per-atom assignment exactly.
      real(dp), intent(in), optional :: valence_density(:, :)
         !! (n_val, n_val), the one-particle density matrix in the
         !! valence-internal orbital basis. Absent means the closed-shell
         !! density -- two electrons in each of the first `n_occupied_valence`
         !! orbitals, none in the rest -- which is what a Hartree-Fock or
         !! Kohn-Sham reference gives and is the only case Paper I treats.
         !!
         !! **This is the argument a correlated wave function needs.** An MCSCF
         !! density is not idempotent -- its active orbitals carry fractional
         !! occupations -- so the closed-shell shortcut does not apply to it.
         !! Nothing else in the construction cares; the quasi-atomic orbitals
         !! never see an occupation.
         !!
         !! **It has to be in the basis of the `valence_internal` passed with
         !! it, and that basis is not reproducible.** Where the molecule has
         !! degenerate valence orbitals a threaded SCF returns a different basis
         !! of the degenerate space from one run to the next, and a density
         !! built against another run's orbitals describes a different set of
         !! orbitals in the same numbers with nothing here able to detect it.
         !! Build both from one wave function and pass them together.

      real(dp), allocatable :: projection(:, :)
      real(dp), allocatable :: claims(:, :), gram(:, :), half(:, :), orthogonal(:, :)
      real(dp), allocatable :: density(:, :), work(:, :)
      integer, allocatable :: first(:), last(:)
      integer :: n_ao, n_val, natm, i, k

      if (error%has_error()) return
      n_ao = size(valence_internal, 1)
      n_val = size(valence_internal, 2)
      natm = size(atomic_numbers)

      if (sum(count) /= n_val) then
         call error%set(ERROR_VALIDATION, "the atoms contribute "//to_char(sum(count))// &
                        " valence minimal-basis orbitals but the valence-internal "// &
                        "space has "//to_char(n_val)//". These are the same number by "// &
                        "construction, so they disagreeing means the space was "// &
                        "assembled from the wrong columns.")
         return
      end if

      ! Raw free-atom orbitals, not the orthogonalized ones: the decomposition
      ! below is per atom, and within one atom the metric is already the
      ! identity.
      allocate (projection(n_val, size(mixed, 2)))
      call pic_gemm(valence_internal, mixed, projection, transa="T")

      ! One subspace unless told otherwise, which is the unrestricted case.
      if (present(subspaces)) then
         allocate (first(size(subspaces)), last(size(subspaces)))
         first = subspaces
         do k = 1, size(subspaces) - 1
            last(k) = subspaces(k + 1) - 1
         end do
         last(size(subspaces)) = n_val
         if (first(1) /= 1 .or. any(last < first)) then
            call error%set(ERROR_VALIDATION, "the subspaces do not partition the "// &
                           "valence space: they must start at orbital one, ascend, "// &
                           "and none may be empty.")
            return
         end if
      else
         allocate (first(1), last(1))
         first(1) = 1
         last(1) = n_val
      end if

      call select_atomic_claims(projection, offset, count, natm, first, last, &
                                claims, result%atom_of, result%subspace_of, error, &
                                starved)
      if (error%has_error()) return
      if (present(starved)) then
         if (starved) return
      end if

      ! Orthonormal within an atom, not between atoms. Symmetric
      ! orthogonalization moves every orbital as little as possible, which is
      ! what keeps them atomic.
      allocate (gram(n_val, n_val), orthogonal(n_val, n_val))
      call pic_gemm(claims, claims, gram, transa="T")
      call inverse_sqrt(gram, half, error)
      if (error%has_error()) return
      call pic_gemm(claims, half, orthogonal)

      call refine_atomic_character(orthogonal, projection, result%atom_of, &
                                   result%subspace_of, offset, count, result%sweeps, &
                                   result%atomic_character, error, &
                                   character_of=result%character_of)
      if (error%has_error()) return

      allocate (result%orbitals(n_ao, n_val))
      call pic_gemm(valence_internal, orthogonal, result%orbitals)
      result%to_valence_internal = orthogonal
      result%n_quao = n_val

      ! The density in the quasi-atomic basis, Paper I eq (5.2), where U is the
      ! transformation from valence-internal orbitals to quasi-atomic ones:
      !
      !     P = U^T D U
      !
      ! U is orthogonal -- symmetric orthogonalization followed by Jacobi
      ! rotations -- so the trace survives and the quasi-atomic populations sum
      ! to the electron count for any D whatsoever, idempotent or not.
      allocate (density(n_val, n_val))
      if (present(valence_density)) then
         if (size(valence_density, 1) /= n_val .or. size(valence_density, 2) /= n_val) then
            call error%set(ERROR_VALIDATION, "the supplied density is "// &
                           to_char(size(valence_density, 1))//" by "// &
                           to_char(size(valence_density, 2))//" but the "// &
                           "valence-internal space has "//to_char(n_val)//" orbitals. "// &
                           "It has to be expressed in that basis, in that order.")
            return
         end if
         density = valence_density
      else
         density = 0.0_dp
         do i = 1, n_occupied_valence
            density(i, i) = 2.0_dp
         end do
      end if

      allocate (work(n_val, n_val), result%population_bond_order(n_val, n_val))
      call pic_gemm(density, orthogonal, work)
      call pic_gemm(orthogonal, work, result%population_bond_order, transa="T")

      deallocate (projection, claims, gram, half, orthogonal, density, work)
   end subroutine quasi_atomic_orbitals

   subroutine orient_quasi_atomic_orbitals(quao, error, max_sweeps)
      !! Rotate within each atom until the bonding pattern appears
      !!
      !! Paper I eq (5.5). The sum of squares of each interatomic block of the
      !! population-bond-order matrix is invariant under rotations inside an
      !! atom -- eq (5.4) -- so no intra-atomic rotation can change how much
      !! bonding there is between two atoms. What it can change is how that
      !! total is distributed, and maximizing the sum of *fourth* powers
      !! concentrates it: a few large bond orders and many near-zero ones,
      !! rather than the same total smeared across every pair. Hybridization is
      !! an output of that rather than an assumption.
      !!
      !! Each 2x2 rotation is closed form. For orbitals `a` and `b` on one atom,
      !! write each partner `k` on another atom as the complex number
      !! `z_k = P_ak + i P_bk`. Rotating `(a, b)` by `theta` turns `z_k` into
      !! `z_k exp(-i theta)`, and since `x^4 + y^4 = (3 |z|^4 + Re z^4) / 4`, the
      !! pair's share of eq (5.5) is
      !!
      !!     sum_k (3 |z_k|^4)/4 + Re( W exp(-4 i theta) )/4,   W = sum_k z_k^4
      !!
      !! The first term is fixed, and the second is largest at
      !! `theta = arg(W)/4`, where it has risen by `(|W| - Re W)/4`.
      ! TODO(mqc): overwrites `quao%sweeps`, which `quasi_atomic_orbitals` set
      ! from the refinement and `quao_result_t` documents as the refinement's
      ! count. After orientation the field silently means something else.
      type(quao_result_t), intent(inout) :: quao
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_sweeps
         !! At least one; `ORIENT_MAX_SWEEPS` when absent

      real(dp), allocatable :: p(:, :), rotation(:, :)
      integer, allocatable :: pairs(:, :)
      logical, allocatable :: on_other_atom(:)
      complex(dp) :: quartic
      real(dp) :: amplitude, theta, initial_sum, final_sum
      real(dp) :: best_gain, widest_angle
      character(len=16) :: gain_text
      integer :: n_orb, ipair, a, b, sweep, sweep_limit
      logical :: rotated

      if (error%has_error()) return
      sweep_limit = ORIENT_MAX_SWEEPS
      if (present(max_sweeps)) sweep_limit = max_sweeps
      if (sweep_limit < 1) then
         call error%set(ERROR_VALIDATION, "the orientation needs at least one sweep, "// &
                        "and was allowed "//to_char(sweep_limit)//".")
         return
      end if

      n_orb = quao%n_quao
      p = quao%population_bond_order
      allocate (rotation(n_orb, n_orb), on_other_atom(n_orb))
      rotation = 0.0_dp
      do a = 1, n_orb
         rotation(a, a) = 1.0_dp
      end do
      ! Same atom, same subspace: a rotation between atoms would reassign an
      ! orbital, and one across subspaces would break a restricted wave
      ! function's invariance.
      call rotatable_pairs(quao%atom_of, quao%subspace_of, .true., pairs)

      initial_sum = orientation_sum(p, quao%atom_of)
      quao%sweeps = 0
      rotated = .true.
      best_gain = 0.0_dp
      widest_angle = 0.0_dp
      do sweep = 1, sweep_limit
         rotated = .false.
         best_gain = 0.0_dp
         widest_angle = 0.0_dp
         do ipair = 1, size(pairs, 2)
            a = pairs(1, ipair)
            b = pairs(2, ipair)
            on_other_atom = quao%atom_of /= quao%atom_of(a)
            quartic = sum(cmplx(p(a, :), p(b, :), kind=dp)**4, mask=on_other_atom)
            amplitude = 0.25_dp*abs(quartic)
            ! The pair's share does not depend on the angle at all.
            if (amplitude < ORIENT_FLAT) cycle
            theta = 0.25_dp*atan2(aimag(quartic), real(quartic, dp))
            if (abs(theta) < ORIENT_ANGLE_TOL) cycle
            best_gain = max(best_gain, amplitude - 0.25_dp*real(quartic, dp))
            widest_angle = max(widest_angle, abs(theta))
            call rotate_rows(p, a, b, cos(theta), sin(theta))
            call rotate_columns(p, a, b, cos(theta), sin(theta))
            call rotate_columns(rotation, a, b, cos(theta), sin(theta))
            rotated = .true.
         end do
         quao%sweeps = sweep
         if (.not. rotated) exit
      end do
      quao%orientation_gain = best_gain
      quao%orientation_angle = widest_angle

      ! Out of sweeps with rotations still above the angle tolerance. Two
      ! orbitals on one atom coupled through a bond to two on another -- a
      ! carbonyl C and O -- can creep towards their optimum by well under a
      ! tenth of a percent per sweep, so the angle stays above tolerance long
      ! after eq (5.5) has stopped moving. Accepted, and flagged, when the
      ! largest gain of the last sweep is negligible; refused otherwise.
      if (rotated) then
         if (best_gain >= ORIENT_GAIN_CRIT) then
            write (gain_text, "(es8.1)") best_gain
            call error%set(ERROR_VALIDATION, "the orientation did not settle in "// &
                           to_char(sweep_limit)//" sweeps; the last still "// &
                           "raised the functional by "//trim(adjustl(gain_text))// &
                           ". Raise properties.bonding_analysis.orientation_max_sweeps.")
            return
         end if
         quao%orientation_stalled = .true.
      end if

      ! Every accepted rotation raises eq (5.5), so a lower final value can
      ! only mean the rotations were applied in the opposite sense to the one
      ! the closed form solved for.
      final_sum = orientation_sum(p, quao%atom_of)
      if (final_sum < initial_sum - ORIENT_SENSE_SLACK) then
         call error%set(ERROR_VALIDATION, "the orientation reduced the functional it "// &
                        "is supposed to maximize, which means the rotations are being "// &
                        "applied with the wrong sense.")
         return
      end if
      quao%orientation_sum = final_sum

      call fix_phases(p, quao%atom_of, rotation)

      block
         real(dp), allocatable :: rotated_orbitals(:, :)
         allocate (rotated_orbitals(size(quao%orbitals, 1), n_orb))
         call pic_gemm(quao%orbitals, rotation, rotated_orbitals)
         call move_alloc(rotated_orbitals, quao%orbitals)
         allocate (rotated_orbitals(size(quao%to_valence_internal, 1), n_orb))
         call pic_gemm(quao%to_valence_internal, rotation, rotated_orbitals)
         call move_alloc(rotated_orbitals, quao%to_valence_internal)
      end block
      quao%population_bond_order = p
      quao%oriented = .true.
   end subroutine orient_quasi_atomic_orbitals

   pure function orientation_sum(p, atom_of) result(total)
      !! Paper I eq (5.5), the quantity the orientation maximizes
      real(dp), intent(in) :: p(:, :)
      integer, intent(in) :: atom_of(:)
      real(dp) :: total
      integer :: i, j

      total = 0.0_dp
      do i = 1, size(atom_of)
         do j = i + 1, size(atom_of)
            if (atom_of(i) == atom_of(j)) cycle
            total = total + p(i, j)**4
         end do
      end do
   end function orientation_sum

   subroutine fix_phases(p, atom_of, rotation)
      !! Choose signs so the strongest interatomic coupling of each orbital is positive
      !!
      !! Nothing up to here fixes a phase -- eigenvectors come back however
      !! LAPACK produced them -- so bond-order signs are arbitrary. This makes
      !! them reproducible rather than meaningful: a convention, not physics.
      !! The kinetic bond order does not need it, since both of its factors
      !! change sign together. Orbitals are visited in index order and each
      !! sees the signs earlier ones were given.
      real(dp), intent(inout) :: p(:, :)
      integer, intent(in) :: atom_of(:)
      real(dp), intent(inout) :: rotation(:, :)

      integer :: orb, strongest

      do orb = 1, size(atom_of)
         ! First of equals on a tie; zero when the orbital has no partner atom.
         strongest = maxloc(abs(p(orb, :)), dim=1, mask=atom_of /= atom_of(orb))
         if (strongest == 0) cycle
         if (.not. p(orb, strongest) < 0.0_dp) cycle
         p(orb, :) = -p(orb, :)
         p(:, orb) = -p(:, orb)
         rotation(:, orb) = -rotation(:, orb)
      end do
   end subroutine fix_phases

   subroutine kinetic_bond_orders(quao, kinetic_ao, kbo, error, interference)
      !! Kinetic bond orders, kcal/mol, negative for a bonding interaction
      !!
      !!     k_{Aa,Bb} = 0.1 * < Aa | -1/2 nabla^2 | Bb > * p_{Aa,Bb}
      !!
      !! Paper II eq (2). For orthonormal quasi-atomic orbitals the kinetic
      !! interference energy is the product of the bond order and the kinetic
      !! integral; the factor of 0.1 is empirical, bringing the numbers onto the
      !! scale of tabulated bond energies.
      !!
      !! **The sign is meaningful where the bond order's is not.** A phase flip
      !! changes both factors, so it cancels in the product. Negative means
      !! bonding.
      !!
      !! **Two quantities, and they get confused.** Paper II eq (1) is the
      !! kinetic interference energy `t = p * T` in hartree; eq (2) is the
      !! kinetic bond order `k = 0.1 t`, in kcal/mol. Both are returned so
      !! neither has to be reconstructed from the other.
      type(quao_result_t), intent(in) :: quao
      real(dp), intent(in) :: kinetic_ao(:, :)    !! (n_ao, n_ao)
      real(dp), allocatable, intent(out) :: kbo(:, :)   !! (n_quao, n_quao)
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: interference(:, :)
         !! The unscaled kinetic interference energy `p * T`, in hartree --
         !! Paper II eq (1) rather than eq (2). Its sum over all orbital pairs,
         !! core included, is the total kinetic energy, and it is what
         !! `mqc_czt_ieda` decomposes.

      real(dp), allocatable :: work(:, :), t(:, :)
      integer :: n

      if (error%has_error()) return
      if (.not. quao%oriented) then
         call error%set(ERROR_VALIDATION, "kinetic bond orders are defined for the "// &
                        "oriented quasi-atomic orbitals. Before orientation the "// &
                        "interatomic couplings are distributed arbitrarily within each "// &
                        "atom, so the individual numbers would not mean anything.")
         return
      end if

      n = quao%n_quao
      allocate (work(size(kinetic_ao, 1), n), t(n, n), kbo(n, n))
      call pic_gemm(kinetic_ao, quao%orbitals, work)
      call pic_gemm(quao%orbitals, work, t, transa="T")

      if (present(interference)) interference = t*quao%population_bond_order
      kbo = KBO_SCALE*t*quao%population_bond_order*HARTREE_TO_KCALMOL
      deallocate (work, t)
   end subroutine kinetic_bond_orders

   subroutine split_localize(quao, n_occupied_valence, valence_internal, localized, error)
      !! Localize the occupied and the empty valence orbitals separately
      !!
      !! Paper I section V.C: *"each localized occupied molecular orbital as
      !! well as each localized valence-virtual orbital shall cover as few of
      !! the oriented quasi-atomic orbitals as possible"*. Done in the two
      !! blocks independently -- which is what "split" means, and what makes the
      !! result different from ordinary localization.
      !!
      !! The occupied space alone gives bonds and lone pairs; the empty valence
      !! space alone gives the antibonding partner of each, on the same two
      !! atoms, from a calculation never told which atoms are bonded.
      !!
      !! The criterion is a fourth power again, and for the same reason as the
      !! orientation: the sum of squares of each column is fixed at one by
      !! orthonormality, so only higher moments can distinguish a spread-out
      !! orbital from a concentrated one. The closed form is Paper I
      !! eqs (A17)-(A24), **with the two corrigenda from Paper II** -- in (A22a)
      !! `2 P_1112` should read `2 P_1122`, and in `P_c` `-6 P_1112` should read
      !! `-6 P_1122`. Both are implemented in their corrected form; using the
      !! printed ones rotates by a wrong angle that still converges, to a
      !! different and worse answer.
      type(quao_result_t), intent(in) :: quao
      integer, intent(in) :: n_occupied_valence
      real(dp), intent(in) :: valence_internal(:, :)   !! (n_ao, n_val)
      real(dp), allocatable, intent(out) :: localized(:, :)
         !! (n_ao, n_val): the localized occupied orbitals, then the localized
         !! empty valence ones, in that order
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: block_t(:, :), coefficients(:, :), transform(:, :)
      integer :: n_val, n_empty

      if (error%has_error()) return
      if (.not. quao%oriented) then
         call error%set(ERROR_VALIDATION, "split localization is defined against the "// &
                        "*oriented* quasi-atomic orbitals: the criterion counts how "// &
                        "many of them each orbital covers, and before orientation they "// &
                        "are an arbitrary basis of each atom's space.")
         return
      end if

      n_val = quao%n_quao
      n_empty = n_val - n_occupied_valence
      allocate (localized(size(valence_internal, 1), n_val))
      allocate (transform(n_val, n_val))
      transform = 0.0_dp

      ! `to_valence_internal(k, a)` expands quasi-atomic orbital `a` in the
      ! valence-internal orbitals, so its transpose gives < QUAO | orbital >,
      ! which is the overlap the criterion is built from.
      if (n_occupied_valence > 0) then
         coefficients = transpose(quao%to_valence_internal(1:n_occupied_valence, :))
         call fourth_power_localize(coefficients, block_t, error)
         if (error%has_error()) return
         transform(1:n_occupied_valence, 1:n_occupied_valence) = block_t
         deallocate (coefficients, block_t)
      end if

      if (n_empty > 0) then
         coefficients = transpose(quao%to_valence_internal(n_occupied_valence + 1:, :))
         call fourth_power_localize(coefficients, block_t, error)
         if (error%has_error()) return
         transform(n_occupied_valence + 1:, n_occupied_valence + 1:) = block_t
         deallocate (coefficients, block_t)
      end if

      call pic_gemm(valence_internal, transform, localized)
      deallocate (transform)
   end subroutine split_localize

   subroutine fourth_power_localize(coefficients, transform, error)
      !! Maximize `sum_n sum_alpha (C T)_{alpha n}^4` over orthogonal `T`
      !!
      !! Paper I Appendix A.2. Each 2x2 rotation is closed form; the quarter
      !! angle is what makes it a fourth-power criterion rather than a
      !! second-power one.
      ! TODO(mqc): converging on sweep `ORIENT_MAX_SWEEPS` exactly is reported
      ! as a failure to converge, because the test after the loop is on the
      ! sweep counter rather than on whether the last sweep rotated anything.
      real(dp), intent(in) :: coefficients(:, :)   !! (n_quao, n_orb)
      real(dp), allocatable, intent(out) :: transform(:, :)   !! (n_orb, n_orb)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: r(:, :), ci(:), cj(:)
      real(dp) :: p1111, p2222, p1122, p1112, p2221, pc, ps, q, gamma, c, s
      integer :: n, i, j, sweep
      logical :: moved

      if (error%has_error()) return
      n = size(coefficients, 2)
      allocate (transform(n, n), r(size(coefficients, 1), n))
      allocate (ci(size(coefficients, 1)), cj(size(coefficients, 1)))
      transform = 0.0_dp
      do i = 1, n
         transform(i, i) = 1.0_dp
      end do
      r = coefficients
      if (n < 2) return

      do sweep = 1, ORIENT_MAX_SWEEPS
         moved = .false.
         do i = 1, n
            do j = i + 1, n
               p1111 = sum(r(:, i)**4)
               p2222 = sum(r(:, j)**4)
               p1122 = sum((r(:, i)*r(:, j))**2)
               p1112 = sum(r(:, i)**3*r(:, j))
               p2221 = sum(r(:, j)**3*r(:, i))

               ! Paper II's corrigenda are in these two lines.
               pc = 0.25_dp*(p1111 + p2222 - 6.0_dp*p1122)
               ps = p1112 - p2221

               q = sqrt(pc*pc + ps*ps)
               if (q < ORIENT_FLAT) cycle
               ! Restricted to (-pi/4, pi/4), which atan2 already gives once
               ! quartered: the functional has period pi/2 in this angle.
               gamma = 0.25_dp*atan2(ps, pc)
               if (abs(gamma) < ORIENT_ANGLE_TOL) cycle

               c = cos(gamma)
               s = sin(gamma)
               ci = c*r(:, i) + s*r(:, j)
               cj = -s*r(:, i) + c*r(:, j)
               r(:, i) = ci
               r(:, j) = cj
               ci = c*transform(:, i) + s*transform(:, j)
               cj = -s*transform(:, i) + c*transform(:, j)
               transform(:, i) = ci
               transform(:, j) = cj
               moved = .true.
            end do
         end do
         if (.not. moved) exit
      end do

      if (sweep >= ORIENT_MAX_SWEEPS) then
         call error%set(ERROR_VALIDATION, "the split localization did not settle in "// &
                        to_char(ORIENT_MAX_SWEEPS)//" sweeps.")
         return
      end if
      deallocate (r, ci, cj)
   end subroutine fourth_power_localize

end module mqc_czt_quao
