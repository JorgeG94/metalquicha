!! Quasi-atomic orbitals: the valence-virtual space
module mqc_libcint_quao
   !! The first step of the Ruedenberg construction that needs integrals.
   !!
   !! The virtual space of a large-basis SCF is mostly polarization and diffuse
   !! functions with no chemical content -- you cannot point at an antibonding
   !! orbital in a cc-pVTZ virtual space, and even the *symmetry* of the lowest
   !! virtual changes with the basis. The valence-virtual orbitals are the part
   !! of it that does have chemical content: the antibonding counterparts of the
   !! occupied valence orbitals, recovered by asking which combinations of
   !! virtuals look most like free-atom minimal-basis orbitals.
   !!
   !! What makes them worth having is that the answer barely depends on the
   !! basis. Paper I reports valence-virtual orbital energies that move less
   !! between 6-31G and aug-cc-pVQZ than the occupied ones do, while the
   !! ordinary virtual spectrum changes beyond recognition.
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
   use mqc_libcint_integrals, only: libcint_molecule_t, mixed_basis_overlap, &
                                    atom_ao_blocks
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

   real(dp), parameter :: ORIENT_CRIT = 1.0e-6_dp
      !! GAMESS's `CVGLOC` default, tested against the rotation angle.
   real(dp), parameter :: ORIENT_FLOOR = 1.0e-8_dp
      !! Below this the functional does not depend on the angle; skip the pair.
   integer, parameter :: ORIENT_MAX_SWEEPS = 2000
   real(dp), parameter :: KBO_SCALE = 0.1_dp
      !! Paper II eq (2). Empirical, and admitted to be: raw kinetic
      !! interference energies in organic molecules run about an order of
      !! magnitude above the bond energies chemistry quotes, so they are scaled
      !! to be comparable with them. Paper II says the treatment of resonance
      !! integrals will be revisited.
   real(dp), parameter :: HARTREE_TO_KCAL = 627.5094740631_dp

   real(dp), parameter :: JACOBI_CRIT = 1.0e-10_dp
      !! GAMESS's `CRIT` in `LOCAL_MODELORB_DRIV`. Both the rotation magnitude
      !! and the angle are tested against it.
   integer, parameter :: JACOBI_MAX_SWEEPS = 10000
      !! GAMESS's `MXITER`. Neither paper states a criterion or a cap, so the
      !! reference implementation's are used rather than invented.

   type :: vvo_result_t
      !! What the valence-virtual extraction leaves behind
      real(dp), allocatable :: orbitals(:, :)
         !! (n_ao, n_vvo), in the AO basis, orthonormal and orthogonal to the
         !! occupied space because they are combinations of virtuals only
      real(dp), allocatable :: external_orbitals(:, :)
         !! (n_ao, n_virt - n_vvo), the rest of the virtual space. Kept because
         !! together with `orbitals` they are a complete, orthonormal virtual
         !! block, which is what anything wanting to *run* in the valence space
         !! rather than analyse it needs -- a CASSCF has to be handed every
         !! orbital, not only the interesting ones.
      real(dp), allocatable :: singular_values(:)
         !! All of them, descending. The first `n_vvo` are retained; the rest
         !! span the valence-external space.
      integer :: n_vvo = 0
      real(dp) :: smallest_retained = 0.0_dp
      real(dp) :: largest_rejected = 0.0_dp
         !! The gap between these two is the diagnostic. Paper I Table I reports
         !! 0.99999 against 0.105-0.272 across eight molecules; anything else
         !! means the projection is not finding a clean valence space.
   end type vvo_result_t

   type :: quao_result_t
      !! Quasi-atomic orbitals and the density in their basis
      real(dp), allocatable :: orbitals(:, :)
         !! (n_ao, n_quao), orthonormal, each one belonging to a single atom
      integer, allocatable :: atom_of(:)
         !! Which atom each quasi-atomic orbital sits on
      integer, allocatable :: subspace_of(:)
         !! Which occupation-restricted subspace each one was drawn from. All
         !! ones when the space is unrestricted, which is the case that says why
         !! it exists: where there is more than one subspace, no rotation may
         !! mix two of them, because a restricted wave function is not invariant
         !! under one that does.
      real(dp), allocatable :: population_bond_order(:, :)
         !! (n_quao, n_quao). Diagonal elements are orbital populations;
         !! interatomic off-diagonal elements are bond orders. Paper I eq (5.2).
      integer :: n_quao = 0
      integer :: sweeps = 0            !! Jacobi sweeps the refinement took
      real(dp), allocatable :: to_valence_internal(:, :)
         !! (n_val, n_quao). Column `i` expands quasi-atomic orbital `i` in the
         !! valence-internal orbitals, so its rows are indexed by those. Kept
         !! because the split-localization works in this space rather than in
         !! the atomic-orbital one.
      logical :: oriented = .false.
      real(dp) :: orientation_sum = 0.0_dp
         !! Paper I eq (5.5) after orientation
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
      type(libcint_molecule_t), intent(out) :: mol
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
      ! so a small eigenvalue here is a statement about the geometry -- atoms on
      ! top of one another -- rather than a threshold to be tuned.
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
      !! **Both bases have to be orthonormal**, and that is not a detail of
      !! taste. Paper I's Appendix A.1 proves the singular vectors maximize the
      !! projection of one space into the other only under that condition; with
      !! a non-orthogonal basis the largest singular value picks out the
      !! direction with the largest *coefficient*, which is not the same
      !! statement and not the one the method needs. The molecular orbitals are
      !! orthonormal already; the free-atom orbitals are orthonormal within one
      !! atom but not between atoms, so they are symmetrically orthogonalized
      !! across the molecule first -- Paper I's `|A#a>`.
      !!
      !! Paper I notes an earlier formulation of this step "failed to account
      !! for the necessary orthogonalities". This is that account.
      real(dp), intent(in) :: orbitals(:, :)        !! C, (n_ao, n_mo)
      real(dp), intent(in) :: mixed(:, :)           !! < AO | AAMBS >, (n_ao, n_mbs)
      real(dp), intent(in) :: aambs_overlap(:, :)   !! < AAMBS | AAMBS >, (n_mbs, n_mbs)
      real(dp), allocatable, intent(out) :: projection(:, :)   !! (n_mo, n_mbs)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: orthogonalize
         !! Default true, giving `|A#a>`. False gives the raw `|A*a>`, which is
         !! what the *per-atom* decompositions want: within one atom the
         !! free-atom orbitals are already orthonormal, so the metric there is
         !! the identity and orthogonalizing across the molecule would mix in
         !! the neighbours the projection is trying to distinguish from.

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
      !! The extraction asserts it for every element; this consumes it.
      integer, intent(in) :: atomic_numbers(:)
      type(libcint_molecule_t), intent(in) :: aambs
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
                                   claims, atom_of, subspace_of, error)
      !! Which atom gets which orbital, by competition rather than by decree
      !!
      !! Every atom is offered every subspace. The claims are pooled, ranked by
      !! how much of the subspace each one actually captures, and the strongest
      !! are kept until the subspace is full. **An atom's share is therefore an
      !! output.** It is not read off the free-atom minimal basis, and it can
      !! differ from that atom's free-atom valence count.
      !!
      !! That is GAMESS's rule and it is deliberate there -- `LOCAL_PPASVD`
      !! offers every atom every subspace under the comment "NO KEEPER LINE HERE
      !! BECAUSE WE WANT TO LOOK (NOT USE) ALL POSSIBLE SVD VALUES", pools into
      !! one flat list and sorts by (subspace, singular value), and Aaron West's
      !! changelog calls the result auto-assignment. Only the molecule-wide sum
      !! is checked against the free-atom counts; no atom's individual share is.
      !!
      !! **With one subspace this is exactly the old prescribed rule**, and that
      !! is worth knowing rather than hoping. Each atom offers as many claims as
      !! it has minimal-basis valence orbitals, and those counts sum to the
      !! valence dimension by construction -- the check above enforces it -- so
      !! the pool exactly fills the space and the ranking selects everything.
      !! What the ranking is *for* is a restricted space, where the subspaces
      !! are smaller than the whole and the atoms genuinely compete.
      !!
      !! Each claim's coefficients are confined to the rows of its own subspace,
      !! which is what makes the eventual transformation block-diagonal by
      !! construction rather than by assertion: two claims from different
      !! subspaces have disjoint support and are orthogonal whatever else
      !! happens to them.
      real(dp), intent(in) :: projection(:, :)   !! (n_val, n_mbs)
      integer, intent(in) :: offset(:), count(:)
      integer, intent(in) :: natm
      integer, intent(in) :: first(:), last(:)   !! (n_subspaces) row ranges
      real(dp), allocatable, intent(out) :: claims(:, :)     !! (n_val, n_val)
      integer, allocatable, intent(out) :: atom_of(:), subspace_of(:)
      type(error_t), intent(inout) :: error

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

      ! Keep the strongest claims on each subspace until it is full. Selection
      ! rather than a full sort: the pool is a few dozen columns and the winners
      ! are wanted in a particular order anyway.
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

      ! Every atom must end up with something. GAMESS asserts the same and says
      ! what to do about it; an atom with no quasi-atomic orbital has no
      ! population, no bond order and no place in the decomposition.
      do iatom = 1, natm
         if (share(iatom) < 1) then
            call error%set(ERROR_VALIDATION, "atom "//to_char(iatom)//" was left "// &
                           "with no quasi-atomic orbital, so it has no population "// &
                           "and no place in the analysis. The subspaces are too "// &
                           "small, or too few, to give every atom one.")
            return
         end if
      end do

      ! Laid out subspace by subspace and, inside each, atom by atom -- so with
      ! one subspace the ordering is atom-major exactly as before. GAMESS sorts
      ! its survivors by atom label for the same reason.
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
      !! section V.B.1.
      !!
      !! Taken through the *small* Gram matrix `sigma^T sigma`, which is one
      !! atom's minimal-basis dimension square: five for carbon, one for
      !! hydrogen. The left vectors follow as `sigma V / s`. The alternative is
      !! an SVD of a matrix whose row count is the whole valence space, for
      !! exactly the same answer.
      real(dp), intent(in) :: sigma(:, :)     !! (n_rows, m_atom)
      real(dp), allocatable, intent(out) :: block(:, :)   !! (n_rows, n_wanted)
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: weights(:)
         !! The squared projections of the vectors returned, descending. These
         !! are what an atom's claim on a subspace is *ranked* by when several
         !! atoms compete for it, so they are an output rather than an internal.
      integer, intent(in), optional :: n_wanted
         !! How many vectors to return, default all of them. Fewer is asked for
         !! when the subspace is smaller than the atom's minimal basis, where
         !! the atom cannot be given a vector per free-atom orbital because
         !! there are not that many dimensions to give.

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

      ! Only the vectors actually being taken have to be representable. When
      ! the subspace is smaller than the atom's minimal basis the trailing
      ! projections are legitimately zero -- there is no room for them -- and
      ! refusing on that would refuse every restricted space.
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
                                      offset, count, sweeps, functional, error)
      !! Rotate between atoms so each orbital sits as fully as possible on its own
      !!
      !! Paper II's Appendix, eqs (A.1)-(A.11), which replaces the plain
      !! symmetric orthogonalization of Paper I and applies retroactively to the
      !! Hartree-Fock case. It maximizes
      !!
      !!     P = sum_A sum_{a on A} sum_{alpha on A} <A*alpha | Aa>^2
      !!
      !! over orthogonal transformations. Only rotations between orbitals on
      !! *different* atoms change it -- the functional is blind to intra-atomic
      !! mixing, which is what the orientation step later exploits -- so the
      !! sweep skips same-atom pairs entirely.
      !!
      !! Each rotation is closed form. With
      !!
      !!     D = (Q^A_ii - Q^A_jj + Q^B_jj - Q^B_ii) / 2
      !!     F =  Q^A_ij - Q^B_ij
      !!
      !! the optimum is theta = atan2(F, D) / 2. GAMESS computes exactly this in
      !! `LOCAL_MODELORB_JACOBI`, down to the halved angle and the sign flip on
      !! the second atom's contribution, and its thresholds are the ones used
      !! here because neither paper states any.
      real(dp), intent(inout) :: coefficients(:, :)   !! (n_rows, n_quao)
      real(dp), intent(in) :: projection(:, :)        !! (n_rows, n_mbs), raw AAMBS
      integer, intent(in) :: atom_of(:)
      integer, intent(in) :: subspace_of(:)
         !! Rotations across two subspaces are skipped. A restricted wave
         !! function is not invariant under them -- mixing two active orbitals
         !! stops being redundant the moment the space is incomplete -- and one
         !! would also destroy the disjoint support that makes the
         !! transformation block-diagonal to begin with. With a single subspace
         !! the guard never fires.
      integer, intent(in) :: offset(:), count(:)      !! AAMBS range per atom
      integer, intent(out) :: sweeps
      real(dp), intent(out) :: functional
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: amb(:, :), ci(:), cj(:)
      real(dp), allocatable :: ri(:), rj(:)
      real(dp) :: d, f, radius, theta, a, b, qii, qjj, qij
      integer :: n, i, j, k, lo, hi, sweep
      logical :: moved

      if (error%has_error()) return
      n = size(coefficients, 2)
      allocate (amb(n, size(projection, 2)), ci(size(coefficients, 1)))
      allocate (cj(size(coefficients, 1)))
      allocate (ri(size(projection, 2)), rj(size(projection, 2)))

      sweeps = 0
      do sweep = 1, JACOBI_MAX_SWEEPS
         moved = .false.
         ! `< QUAO | A*alpha >` for the current orbitals. Rebuilt each sweep
         ! rather than updated in place: two rows change per rotation, and the
         ! bookkeeping to track that is more error-prone than the multiply.
         call pic_gemm(coefficients, projection, amb, transa="T")

         do i = 1, n
            do j = i + 1, n
               if (atom_of(i) == atom_of(j)) cycle
               if (subspace_of(i) /= subspace_of(j)) cycle   ! invisible to the functional

               d = 0.0_dp
               f = 0.0_dp
               lo = offset(atom_of(i)) + 1
               hi = offset(atom_of(i)) + count(atom_of(i))
               qii = sum(amb(i, lo:hi)**2)
               qjj = sum(amb(j, lo:hi)**2)
               qij = sum(amb(i, lo:hi)*amb(j, lo:hi))
               d = d + 0.5_dp*(qii - qjj)
               f = f + qij

               lo = offset(atom_of(j)) + 1
               hi = offset(atom_of(j)) + count(atom_of(j))
               qii = sum(amb(i, lo:hi)**2)
               qjj = sum(amb(j, lo:hi)**2)
               qij = sum(amb(i, lo:hi)*amb(j, lo:hi))
               d = d + 0.5_dp*(qjj - qii)
               f = f - qij

               radius = sqrt(d*d + f*f)
               if (radius < JACOBI_CRIT) cycle
               ! P is independent of the angle here; any rotation is as good as
               ! any other, so leave the pair alone. Neither paper mentions this
               ! case; it is reachable whenever two orbitals are symmetric
               ! partners.
               theta = 0.5_dp*atan2(f, d)
               if (abs(theta) < JACOBI_CRIT) cycle

               a = cos(theta)
               b = sin(theta)
               ci = coefficients(:, i)
               cj = coefficients(:, j)
               coefficients(:, i) = a*ci + b*cj
               coefficients(:, j) = -b*ci + a*cj
               ! Both rows from the *old* values. Updating in place would feed
               ! the new row i into row j and silently rotate by the wrong angle.
               ri = amb(i, :)
               rj = amb(j, :)
               amb(i, :) = a*ri + b*rj
               amb(j, :) = -b*ri + a*rj
               moved = .true.
            end do
         end do

         sweeps = sweep
         if (.not. moved) exit
      end do

      if (sweeps >= JACOBI_MAX_SWEEPS) then
         call error%set(ERROR_VALIDATION, "the quasi-atomic refinement did not settle "// &
                        "in "//to_char(JACOBI_MAX_SWEEPS)//" sweeps.")
         return
      end if

      call pic_gemm(coefficients, projection, amb, transa="T")
      functional = 0.0_dp
      do i = 1, n
         k = atom_of(i)
         lo = offset(k) + 1
         hi = offset(k) + count(k)
         functional = functional + sum(amb(i, lo:hi)**2)
      end do
      functional = functional/real(n, dp)
      deallocate (amb, ci, cj, ri, rj)
   end subroutine refine_atomic_character

   subroutine valence_virtual_orbitals(orbitals, projection, dims, result, error)
      !! Recover the chemically meaningful part of the virtual space
      !!
      !! The number kept is fixed by counting, not by a threshold:
      !! `n_vvo = n_mbs - n_occ` exactly. The singular values are reported so a
      !! caller can check there is a gap where the cut falls, but they do not
      !! decide where it falls -- a construction that selected on magnitude
      !! would return a different number of orbitals for the same molecule in a
      !! different basis, which is the property the method exists to avoid.
      !!
      !! Solved as the eigenproblem of `Sigma Sigma^T` rather than by an SVD.
      !! Appendix A.1 shows these have the same left singular vectors, and the
      !! external block of the rectangular problem is full of degenerate zero
      !! singular values whose vectors are arbitrary -- GAMESS takes the
      !! eigenvalue route for exactly that reason, and this follows it.
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

      ! Negated so the eigenvalues come back ascending in magnitude-reversed
      ! order, i.e. the largest projections first once the sign is undone. The
      ! same device GAMESS uses, and it keeps the retained block at the front
      ! where every downstream index expects it.
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

      ! The same rotation applied to the columns the count rejected. They are
      ! the valence-external space and are wanted for nothing except being the
      ! other half of a complete virtual block.
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
                                    valence_density, subspaces)
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
      integer, intent(in), optional :: subspaces(:)
         !! The valence-internal orbital each occupation-restricted subspace
         !! starts at, ascending, as `ormas_space_t%first_orbital` gives them.
         !! Absent means one subspace covering everything, which is the
         !! unrestricted case and reproduces the prescribed per-atom assignment
         !! exactly -- see `select_atomic_claims` for why that is a theorem
         !! rather than a hope.
      real(dp), intent(in), optional :: valence_density(:, :)
         !! (n_val, n_val), the one-particle density matrix in the
         !! valence-internal orbital basis. Absent means the closed-shell
         !! density -- two electrons in each of the first `n_occupied_valence`
         !! orbitals, none in the rest -- which is what a Hartree-Fock or
         !! Kohn-Sham reference gives and is the only case Paper I treats.
         !!
         !! **This is the argument a correlated wave function needs.** An MCSCF
         !! density is not idempotent: its active orbitals carry fractional
         !! occupations, so it cannot be written as a projector onto some set of
         !! filled orbitals and the closed-shell shortcut below does not apply
         !! to it. Nothing else in the construction cares -- the quasi-atomic
         !! orbitals themselves come from overlaps with free-atom orbitals and
         !! never see an occupation -- so this one argument is the whole
         !! difference between analysing a reference and analysing a correlated
         !! wave function.
         !!
         !! **It has to be in the basis of the `valence_internal` passed with
         !! it, and that basis is not reproducible.** Where the molecule has
         !! degenerate valence orbitals a threaded SCF returns a different basis
         !! of the degenerate space from one run to the next -- measured at 1.4
         !! in the transformation matrix for water between two identical calls,
         !! while the quasi-atomic orbitals and the populations that come out
         !! agree to 4e-13. So a density built against orbitals from a *different*
         !! run of the SCF describes a different set of orbitals in the same
         !! numbers, and nothing here can detect that. Build both from one wave
         !! function and pass them together.

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
      ! below is per atom, and within one atom the metric is already the identity.
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
                                claims, result%atom_of, result%subspace_of, error)
      if (error%has_error()) return

      ! Orthonormal within an atom, not between atoms. Symmetric
      ! orthogonalization is the choice that moves every orbital as little as
      ! possible, which is what keeps them atomic.
      allocate (gram(n_val, n_val), orthogonal(n_val, n_val))
      call pic_gemm(claims, claims, gram, transa="T")
      call inverse_sqrt(gram, half, error)
      if (error%has_error()) return
      call pic_gemm(claims, half, orthogonal)

      call refine_atomic_character(orthogonal, projection, result%atom_of, &
                                   result%subspace_of, offset, count, result%sweeps, &
                                   result%atomic_character, error)
      if (error%has_error()) return

      allocate (result%orbitals(n_ao, n_val))
      call pic_gemm(valence_internal, orthogonal, result%orbitals)
      result%to_valence_internal = orthogonal
      result%n_quao = n_val

      ! The density in the quasi-atomic basis, Paper I eq (5.2): a change of
      ! basis of the one-particle density matrix by U, the transformation from
      ! valence-internal orbitals to quasi-atomic ones.
      !
      !     P = U^T D U
      !
      ! For a closed shell D is 2 on the occupied diagonal and zero elsewhere,
      ! and this collapses to 2 U_occ U_occ^T -- which is what the code used to
      ! compute directly. It is written out in full instead because the general
      ! form is the one a correlated density needs and two routes to the same
      ! matrix is one route too many: the closed-shell shortcut would be exercised
      ! by every test and the general path by almost none.
      !
      ! U is orthogonal -- symmetric orthogonalization followed by Jacobi
      ! rotations -- so the trace survives, and the quasi-atomic populations sum
      ! to the electron count for any D whatsoever. That is the invariant worth
      ! knowing here, because it is the one that does not depend on idempotency.
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

   subroutine orient_quasi_atomic_orbitals(quao, error)
      !! Rotate within each atom until the bonding pattern appears
      !!
      !! Paper I eq (5.5). The sum of squares of each interatomic block of the
      !! population-bond-order matrix is invariant under rotations inside an
      !! atom -- eq (5.4) -- so no intra-atomic rotation can change how much
      !! bonding there is between two atoms. What it can change is how that
      !! total is distributed, and maximizing the sum of *fourth* powers
      !! concentrates it: a few large bond orders and many near-zero ones,
      !! rather than the same total smeared across every pair.
      !!
      !! That is what turns a basis into a picture. Before this step an O-H bond
      !! of 0.93 might sit as 0.66 and 0.63 across two oxygen orbitals; after it,
      !! one oxygen orbital points at the hydrogen and the others do not.
      !! Hybridization is an output here rather than an assumption -- nothing
      !! told it to build sp3 orbitals.
      !!
      !! The rotation angle is closed form and quartered, because the functional
      !! is quartic where the atomic-character one was quadratic. GAMESS's
      !! `ORIEN` computes the same R2 and R3; its thresholds are used here.
      type(quao_result_t), intent(inout) :: quao
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: p(:, :), rotation(:, :), pi(:), pj(:)
      real(dp) :: r2, r3, q, theta, c, s, before, after
      integer :: n, i, j, k, sweep
      logical :: moved

      if (error%has_error()) return
      n = quao%n_quao
      allocate (p(n, n), rotation(n, n), pi(n), pj(n))
      p = quao%population_bond_order
      rotation = 0.0_dp
      do i = 1, n
         rotation(i, i) = 1.0_dp
      end do

      before = orientation_sum(p, quao%atom_of)
      quao%sweeps = 0
      do sweep = 1, ORIENT_MAX_SWEEPS
         moved = .false.
         do i = 1, n
            do j = i + 1, n
               ! Only *within* an atom. A rotation between atoms would change
               ! which atom an orbital belongs to, which is the one thing the
               ! previous stage established.
               if (quao%atom_of(i) /= quao%atom_of(j)) cycle
               ! And not across two subspaces, for the reason the refinement
               ! above gives: a restricted wave function is not invariant under
               ! a rotation that mixes them.
               if (quao%subspace_of(i) /= quao%subspace_of(j)) cycle

               r2 = 0.0_dp
               r3 = 0.0_dp
               do k = 1, n
                  if (quao%atom_of(k) == quao%atom_of(i)) cycle
                  r2 = r2 + p(i, k)**4 - 6.0_dp*p(i, k)**2*p(j, k)**2 + p(j, k)**4
                  r3 = r3 + (p(i, k)**2 - p(j, k)**2)*p(i, k)*p(j, k)
               end do
               r2 = 0.25_dp*r2

               q = sqrt(r2*r2 + r3*r3)
               if (q < ORIENT_FLOOR) cycle
               theta = 0.25_dp*atan2(r3/q, r2/q)
               if (abs(theta) < ORIENT_CRIT) cycle

               c = cos(theta)
               s = sin(theta)
               pi = c*p(i, :) + s*p(j, :)
               pj = -s*p(i, :) + c*p(j, :)
               p(i, :) = pi
               p(j, :) = pj
               pi = c*p(:, i) + s*p(:, j)
               pj = -s*p(:, i) + c*p(:, j)
               p(:, i) = pi
               p(:, j) = pj
               pi = c*rotation(:, i) + s*rotation(:, j)
               pj = -s*rotation(:, i) + c*rotation(:, j)
               rotation(:, i) = pi
               rotation(:, j) = pj
               moved = .true.
            end do
         end do
         quao%sweeps = sweep
         if (.not. moved) exit
      end do

      if (quao%sweeps >= ORIENT_MAX_SWEEPS) then
         call error%set(ERROR_VALIDATION, "the orientation did not settle in "// &
                        to_char(ORIENT_MAX_SWEEPS)//" sweeps.")
         return
      end if

      ! The functional is being maximized, so it cannot come out lower than it
      ! started. It would if the rotation sense disagreed with the angle the
      ! closed form solves for, which is the sort of sign error that otherwise
      ! produces a plausible-looking picture of the wrong thing.
      after = orientation_sum(p, quao%atom_of)
      if (after < before - 1.0e-10_dp) then
         call error%set(ERROR_VALIDATION, "the orientation reduced the functional it "// &
                        "is supposed to maximize, which means the rotations are being "// &
                        "applied with the wrong sense.")
         return
      end if
      quao%orientation_sum = after

      call fix_phases(p, quao%atom_of, rotation)

      block
         real(dp), allocatable :: rotated(:, :)
         allocate (rotated(size(quao%orbitals, 1), n))
         call pic_gemm(quao%orbitals, rotation, rotated)
         call move_alloc(rotated, quao%orbitals)
         allocate (rotated(size(quao%to_valence_internal, 1), n))
         call pic_gemm(quao%to_valence_internal, rotation, rotated)
         call move_alloc(rotated, quao%to_valence_internal)
      end block
      quao%population_bond_order = p
      quao%oriented = .true.
      deallocate (p, rotation, pi, pj)
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
      !! change sign together.
      real(dp), intent(inout) :: p(:, :)
      integer, intent(in) :: atom_of(:)
      real(dp), intent(inout) :: rotation(:, :)

      integer :: n, i, j, best
      real(dp) :: largest

      n = size(atom_of)
      do i = 1, n
         best = 0
         largest = 0.0_dp
         do j = 1, n
            if (atom_of(j) == atom_of(i)) cycle
            if (abs(p(i, j)) > largest) then
               largest = abs(p(i, j))
               best = j
            end if
         end do
         if (best == 0) cycle
         if (p(i, best) >= 0.0_dp) cycle
         p(i, :) = -p(i, :)
         p(:, i) = -p(:, i)
         rotation(:, i) = -rotation(:, i)
      end do
   end subroutine fix_phases

   subroutine kinetic_bond_orders(quao, kinetic_ao, kbo, error, interference)
      !! Kinetic bond orders, kcal/mol, negative for a bonding interaction
      !!
      !!     k_{Aa,Bb} = 0.1 * < Aa | -1/2 nabla^2 | Bb > * p_{Aa,Bb}
      !!
      !! Paper II eq (2). The reasoning is that the covalent binding in an ab
      !! initio wave function comes from kinetic interference between atomic
      !! orbitals, and for orthonormal quasi-atomic orbitals that energy is
      !! simply the product of the bond order and the kinetic integral.
      !!
      !! **The sign is meaningful where the bond order's is not.** A phase flip
      !! changes both factors, so it cancels in the product -- which is exactly
      !! why Paper II introduced this quantity. Negative means bonding, in every
      !! case the papers report where the interaction is manifestly so.
      !!
      !! The factor of 0.1 is empirical and the papers say so. It brings the
      !! numbers onto the scale of tabulated bond energies; it is not derived.
      !!
      !! **Two quantities, and they get confused.** Paper II eq (1) is the
      !! kinetic interference energy `t = p * T` in hartree; eq (2) is the
      !! kinetic bond order `k = 0.1 t`, quoted in kcal/mol. GAMESS prints the
      !! *former* under the heading KEI-BO, so a value of -0.578 there is
      !! -36.3 kcal/mol here. Both are returned so neither has to be
      !! reconstructed from the other.
      type(quao_result_t), intent(in) :: quao
      real(dp), intent(in) :: kinetic_ao(:, :)    !! (n_ao, n_ao)
      real(dp), allocatable, intent(out) :: kbo(:, :)   !! (n_quao, n_quao)
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: interference(:, :)
         !! The unscaled kinetic interference energy `p * T`, in hartree --
         !! Paper II eq (1) rather than eq (2). This is what GAMESS prints, as
         !! its "novel oriented density matrix" and the KEI-BO column, so it is
         !! the quantity to compare against directly. Its sum over all orbital
         !! pairs, core included, is the total kinetic energy, which GAMESS
         !! prints either side of the transformation as a check that the
         !! rotation was orthogonal.

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
      kbo = KBO_SCALE*t*quao%population_bond_order*HARTREE_TO_KCAL
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
      !! Localizing the occupied space alone gives bonds and lone pairs.
      !! Localizing the empty valence space alone gives the antibonding partner
      !! of each. Together they pair up, and that pairing is the point: a
      !! bonding orbital and its antibonding counterpart, on the same two atoms,
      !! from a calculation that was never told which atoms are bonded.
      !!
      !! Paper I notes that conventional localization instead produces
      !! "rabbit-ear" lone pairs on oxygen, where this yields one of sigma type
      !! and one of p type.
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
               if (q < ORIENT_FLOOR) cycle
               ! Restricted to (-pi/4, pi/4), which atan2 already gives once
               ! quartered: the functional has period pi/2 in this angle.
               gamma = 0.25_dp*atan2(ps, pc)
               if (abs(gamma) < ORIENT_CRIT) cycle

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

end module mqc_libcint_quao
