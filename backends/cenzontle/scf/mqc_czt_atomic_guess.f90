!! Free-atom initial guesses for the CPU backend
module mqc_czt_atomic_guess
   !! Builds a molecular starting density by superposing converged free atoms.
   !!
   !! Each distinct element is solved once as an isolated atom -- unrestricted,
   !! at its Hund's-rule ground-state multiplicity, in exactly the basis
   !! functions it contributes to the molecule -- and its density is dropped into
   !! the diagonal block belonging to each atom of that element. The guess is
   !! block-diagonal by construction.
   !!
   !! Two flavours:
   !!
   !! * **SAC** keeps the free atom's own alpha and beta densities, so an
   !!   open-shell atom's guess arrives with its spatial and spin symmetry
   !!   already broken -- what a radical wants, and what a closed shell then has
   !!   to undo.
   !!
   !! * **SAD** spherically averages the total atomic density and halves it, so
   !!   both spins start from the same rotation-invariant block. The
   !!   closed-shell default.
   !!
   !! Atomic solutions are cached for the lifetime of the process, keyed by the
   !! basis data they were converged in; under MPI each rank keeps its own.
   !!
   !! **The cache is module state and is not thread-safe.** Two threads missing
   !! on the same element would both increment `n_cached` and write the same
   !! slot, corrupting a guess rather than failing.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_atomic_guess_common, only: hund_multiplicity
   use mqc_string_utils, only: int_to_text
   use mqc_czt_integrals, only: czt_molecule_t, atom_ao_blocks, subshell_layout
   use mqc_czt_rhf, only: SCF_GUESS_PROJ, rhf_result_t, run_czt_uhf, &
                          SCF_GUESS_CORE, SCF_GUESS_GWH, SCF_GUESS_SAC, SCF_GUESS_SAD
   implicit none
   private

   public :: build_atomic_guess     !! Molecular guess densities from free-atom solutions
   public :: build_restricted_guess  !! Name + density a restricted SCF starts from
   public :: parse_guess_name       !! Deck spelling to SCF_GUESS_*
   public :: guess_display_name     !! SCF_GUESS_* back to a name, for reporting
   public :: hund_multiplicity      !! Ground-state multiplicity of a free atom
   public :: spherical_average      !! Rotationally invariant part of an atomic density
   public :: free_atom_energies     !! Converged energy of each atom on its own
   public :: clear_atomic_cache     !! Drop cached atomic solutions

   integer, parameter :: MAX_CACHED = 64
      !! How many distinct (element, basis) atomic solutions are kept at once,
      !! one per element in the calculation.

   integer, parameter :: ATOMIC_MAX_ITER = 200
   real(dp), parameter :: ATOMIC_ENERGY_TOL = 1.0e-8_dp
   real(dp), parameter :: ATOMIC_DENSITY_TOL = 1.0e-6_dp
      !! The free-atom SCF's convergence settings: looser than a production SCF
      !! because the answer is a guess, with a high iteration cap because a
      !! first-row transition metal at Hund multiplicity is slow.

   integer, parameter :: ATOMIC_INCORE_MAX = 50
      !! Above this many functions the free atom is solved directly rather than
      !! in core. One atom's n^4 tensor is 50 MB at 50 functions, 800 MB at 100.

   type :: atomic_solution_t
      !! One converged free-atom calculation, with both guesses precomputed
      integer :: atomic_number = 0
      integer :: n_ao = 0
      real(dp) :: energy = 0.0_dp
         !! The converged free-atom energy, hartree. Solved in exactly the basis
         !! functions this atom contributes to the parent molecule, so the
         !! difference from an atom-in-molecule energy carries no basis-set
         !! superposition error.
      real(dp), allocatable :: d_alpha(:, :), d_beta(:, :)  !! As converged, for SAC
      real(dp), allocatable :: d_half(:, :)                 !! Averaged and halved, for SAD
      logical :: average_exact = .true.
         !! False when a Cartesian shell of l >= 2 was left unaveraged. See
         !! `spherical_average`.
      ! TODO(mqc): `average_exact` is assigned in `solve_free_atom` and read
      ! nowhere, so the one case where the SAD guess is not spherically
      ! symmetric is recorded and never reported.

      ! The basis this was solved in, kept verbatim: a basis set *name* would
      ! not identify it, since the same name reaches different atoms as
      ! different shells once general contractions have been merged.
      integer, allocatable :: bas(:, :)
      real(dp), allocatable :: env(:)
      logical :: cartesian = .false.
      logical :: valid = .false.
   end type atomic_solution_t

   type(atomic_solution_t), save :: cache(MAX_CACHED)
   integer, save :: n_cached = 0

contains

   subroutine parse_guess_name(text, kind, error)
      !! Deck spelling to a guess constant, refusing anything unrecognised
      !!
      !! The schema validates that the *key* exists; nothing but this validates
      !! the value. `"auto"` and `""` resolve to SAD, this backend's default.
      character(len=*), intent(in) :: text
      integer, intent(out) :: kind
      type(error_t), intent(inout) :: error

      select case (trim(adjustl(text)))
      case ("auto", "")
         kind = SCF_GUESS_SAD
      case ("core")
         kind = SCF_GUESS_CORE
      case ("gwh")
         kind = SCF_GUESS_GWH
      case ("sac")
         kind = SCF_GUESS_SAC
      case ("sad")
         kind = SCF_GUESS_SAD
      case ("basis_set_projection", "projection")
         kind = SCF_GUESS_PROJ
      case default
         kind = SCF_GUESS_SAD
         call error%set(ERROR_VALIDATION, "unknown initial guess '"//trim(adjustl(text))// &
                        "'; the CPU backend has core, gwh, sac, sad and "// &
                        "basis_set_projection")
      end select
   end subroutine parse_guess_name

   subroutine build_restricted_guess(mol, guess_name, guess_kind, guess_total, error)
      !! Resolve a deck guess name and build the density a restricted SCF starts from
      !!
      !! For SAD and SAC, the free-atom spin channels summed into the one total a
      !! closed shell takes. An atomic guess that will not build warns and falls
      !! back to GWH rather than failing the run.
      !!
      !! Core and GWH need nothing here -- the SCF builds them from H and S -- so
      !! `guess_total` comes back unallocated and only the kind is set.
      !! Basis-set projection needs a sub-SCF ladder the caller drives and is
      !! refused.
      type(czt_molecule_t), intent(in) :: mol
      character(len=*), intent(in) :: guess_name
      integer, intent(out) :: guess_kind
      real(dp), allocatable, intent(out) :: guess_total(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: guess_a(:, :), guess_b(:, :)
      type(error_t) :: guess_error

      call parse_guess_name(guess_name, guess_kind, error)
      if (error%has_error()) return

      if (guess_kind == SCF_GUESS_PROJ) then
         call error%set(ERROR_VALIDATION, "the basis-set-projection guess needs a sub-SCF "// &
                        "ladder and is not available here; choose sad, gwh or core")
         return
      end if

      if (guess_kind == SCF_GUESS_SAC .or. guess_kind == SCF_GUESS_SAD) then
         call build_atomic_guess(mol, guess_kind, guess_a, guess_b, guess_error)
         if (guess_error%has_error()) then
            call logger%warning("initial guess: "//guess_error%get_message()// &
                                " -- falling back to gwh")
            guess_kind = SCF_GUESS_GWH
         else
            guess_total = guess_a + guess_b
         end if
      end if
   end subroutine build_restricted_guess

   pure function guess_display_name(kind) result(name)
      !! A guess constant back to its deck spelling, so a run can report itself
      integer, intent(in) :: kind
      character(len=:), allocatable :: name

      select case (kind)
      case (SCF_GUESS_CORE)
         name = "core"
      case (SCF_GUESS_GWH)
         name = "gwh"
      case (SCF_GUESS_SAC)
         name = "sac"
      case (SCF_GUESS_PROJ)
         name = "basis_set_projection"
      case (SCF_GUESS_SAD)
         name = "sad"
      case default
         name = "unknown"
      end select
   end function guess_display_name

   subroutine spherical_average(atom_mol, density, averaged, exact)
      !! The rotationally invariant part of a one-atom density
      !!
      !! A spherically symmetric density has one form and only one:
      !!
      !!    D[(n l m), (n' l' m')] = delta_ll' delta_mm' d[n n' l]
      !!
      !! so averaging over all rotations amounts to keeping the same-l same-m
      !! entries, averaging them over the 2l+1 components, and zeroing
      !! everything else. Same-l blocks from *different* contraction columns
      !! survive -- 1s with 2s is a legitimate radial off-diagonal -- which is
      !! why this walks columns rather than shells.
      !!
      !! Exact for spherical harmonics at any l, and for Cartesian s and p. Not
      !! exact for Cartesian d and above, where six Cartesian d functions are
      !! five d plus an s and averaging over them would mix angular momenta:
      !! those columns are left as converged and `exact` comes back false.
      type(czt_molecule_t), intent(in) :: atom_mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(out) :: averaged(:, :)
      logical, intent(out) :: exact

      integer, allocatable :: ang(:), first(:), ncomp(:)
      integer :: n_sub, a, b, k
      real(dp) :: mean

      call subshell_layout(atom_mol, ang, first, ncomp, n_sub)

      averaged = 0.0_dp
      exact = .true.

      do b = 1, n_sub
         do a = 1, n_sub
            ! Different angular momenta cannot mix in a spherical density --
            ! in a basis of spherical harmonics. A CARTESIAN shell of l>=2 is
            ! not pure l: (xx+yy+zz) is s-like, so a spherically symmetric
            ! density genuinely has s-d entries, and zeroing them throws real
            ! density away. Measured on oxygen: Tr(D S) came back 8.014 rather
            ! than 8, and cc-pVTZ, which adds f, 7.933.
            if (ang(a) /= ang(b)) then
               if (atom_mol%cartesian .and. max(ang(a), ang(b)) >= 2) then
                  exact = .false.
                  averaged(first(a) + 1:first(a) + ncomp(a), &
                           first(b) + 1:first(b) + ncomp(b)) = &
                     density(first(a) + 1:first(a) + ncomp(a), &
                             first(b) + 1:first(b) + ncomp(b))
               end if
               cycle
            end if

            if (atom_mol%cartesian .and. ang(a) >= 2) then
               exact = .false.
               averaged(first(a) + 1:first(a) + ncomp(a), first(b) + 1:first(b) + ncomp(b)) = &
                  density(first(a) + 1:first(a) + ncomp(a), first(b) + 1:first(b) + ncomp(b))
               cycle
            end if

            mean = 0.0_dp
            do k = 1, ncomp(a)
               mean = mean + density(first(a) + k, first(b) + k)
            end do
            mean = mean/real(ncomp(a), dp)
            do k = 1, ncomp(a)
               averaged(first(a) + k, first(b) + k) = mean
            end do
         end do
      end do
   end subroutine spherical_average

   subroutine build_atomic_guess(mol, kind, d_alpha, d_beta, error)
      !! Molecular guess densities, one atomic block at a time
      !!
      !! The blocks are neutral free atoms, so the guess carries sum(Z)
      !! electrons whatever the molecule's own charge is; the Fock matrix it
      !! builds is then occupied with the right number.
      type(czt_molecule_t), intent(in) :: mol
      integer, intent(in) :: kind    !! SCF_GUESS_SAC or SCF_GUESS_SAD
      real(dp), allocatable, intent(out) :: d_alpha(:, :), d_beta(:, :)
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: atom_mol
      integer, allocatable :: offsets(:), counts(:)
      integer :: iatom, idx, z, lo, hi

      if (kind /= SCF_GUESS_SAC .and. kind /= SCF_GUESS_SAD) then
         call error%set(ERROR_VALIDATION, "atomic guess: asked for a guess that needs no "// &
                        "free-atom solutions")
         return
      end if

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)

      allocate (d_alpha(mol%nao, mol%nao), d_beta(mol%nao, mol%nao))
      d_alpha = 0.0_dp
      d_beta = 0.0_dp

      do iatom = 1, mol%natm
         z = nint(mol%charges(iatom))
         if (z <= 0) cycle   ! a ghost centre carries no electrons to superpose

         ! Sliced before the cache is consulted, because the slice is the key.
         call mol%atom_subset(iatom, atom_mol, error)
         if (error%has_error()) return

         idx = cached_index(z, atom_mol)
         if (idx == 0) then
            if (n_cached >= MAX_CACHED) then
               ! Refused rather than evicted: a run that has seen 64 distinct
               ! elements will see them again.
               call error%set(ERROR_VALIDATION, "atomic guess: solution cache exhausted")
               call atom_mol%destroy()
               return
            end if
            n_cached = n_cached + 1
            call solve_free_atom(atom_mol, z, cache(n_cached), error)
            if (error%has_error()) then
               ! Left invalid so a retry re-solves rather than reading a
               ! half-filled entry.
               cache(n_cached)%valid = .false.
               n_cached = n_cached - 1
               call atom_mol%destroy()
               return
            end if
            idx = n_cached
         end if

         if (cache(idx)%n_ao /= counts(iatom)) then
            call error%set(ERROR_VALIDATION, "atomic guess: the free atom's basis function "// &
                           "count does not match the block it belongs to")
            call atom_mol%destroy()
            return
         end if

         lo = offsets(iatom) + 1
         hi = offsets(iatom) + counts(iatom)
         if (kind == SCF_GUESS_SAC) then
            d_alpha(lo:hi, lo:hi) = cache(idx)%d_alpha
            d_beta(lo:hi, lo:hi) = cache(idx)%d_beta
         else
            d_alpha(lo:hi, lo:hi) = cache(idx)%d_half
            d_beta(lo:hi, lo:hi) = cache(idx)%d_half
         end if

         call atom_mol%destroy()
      end do

      deallocate (offsets, counts)
   end subroutine build_atomic_guess

   subroutine solve_free_atom(atom_mol, atomic_number, solution, error)
      !! Converge one isolated atom and precompute both guesses from it
      !!
      !! The atomic SCF starts from GWH, never from an atomic guess, which is
      !! what keeps this from recursing into itself.
      type(czt_molecule_t), intent(in) :: atom_mol
      integer, intent(in) :: atomic_number
      type(atomic_solution_t), intent(inout) :: solution
      type(error_t), intent(inout) :: error

      type(rhf_result_t) :: scf
      real(dp), allocatable :: total(:, :), averaged(:, :)
      integer :: multiplicity
      logical :: exact

      multiplicity = hund_multiplicity(atomic_number)

      ! mqc: scf-subset -- a per-element atomic SCF that exists only to build the
      ! SAD guess. Its thresholds are the ATOMIC_* constants above rather than
      ! the molecule's, deliberately: the answer here is a starting density.
      call run_czt_uhf(atom_mol, atomic_number, multiplicity, ATOMIC_MAX_ITER, &
                       ATOMIC_ENERGY_TOL, ATOMIC_DENSITY_TOL, .false., scf, error, &
                       diis_vectors=8, in_core=(atom_mol%nao <= ATOMIC_INCORE_MAX), &
                       guess=SCF_GUESS_GWH)
      if (error%has_error()) then
         call error%add_context("atomic guess: free atom Z="//int_to_text(atomic_number))
         return
      end if
      if (.not. scf%converged) then
         call error%set(ERROR_VALIDATION, "atomic guess: the free atom Z="// &
                        int_to_text(atomic_number)//" did not converge in "// &
                        int_to_text(ATOMIC_MAX_ITER)//" cycles")
         return
      end if

      allocate (total(atom_mol%nao, atom_mol%nao), averaged(atom_mol%nao, atom_mol%nao))
      total = scf%density + scf%density_beta
      call spherical_average(atom_mol, total, averaged, exact)

      solution%atomic_number = atomic_number
      solution%n_ao = atom_mol%nao
      solution%energy = scf%energy
      solution%d_alpha = scf%density
      solution%d_beta = scf%density_beta
      solution%d_half = 0.5_dp*averaged
      solution%average_exact = exact
      solution%bas = atom_mol%bas
      solution%env = atom_mol%env
      solution%cartesian = atom_mol%cartesian
      solution%valid = .true.

      deallocate (total, averaged)
   end subroutine solve_free_atom

   function cached_index(atomic_number, atom_mol) result(idx)
      !! Position of a cached solution for this atom in this basis, or 0
      !!
      !! The basis data is compared element by element rather than hashed.
      integer, intent(in) :: atomic_number
      type(czt_molecule_t), intent(in) :: atom_mol
      integer :: idx, i

      idx = 0
      do i = 1, n_cached
         if (.not. cache(i)%valid) cycle
         if (cache(i)%atomic_number /= atomic_number) cycle
         if (cache(i)%cartesian .neqv. atom_mol%cartesian) cycle
         if (size(cache(i)%bas, 2) /= size(atom_mol%bas, 2)) cycle
         if (size(cache(i)%env) /= size(atom_mol%env)) cycle
         if (any(cache(i)%bas /= atom_mol%bas)) cycle
         if (any(cache(i)%env /= atom_mol%env)) cycle
         idx = i
         return
      end do
   end function cached_index

   subroutine free_atom_energies(mol, energies, error)
      !! The energy each atom would have on its own, in its own basis
      !!
      !! One unrestricted SCF per *distinct* element at its Hund's-rule ground
      !! state, sharing the cache the SAD and SAC guesses fill.
      !!
      !! **Each atom is solved in exactly the basis functions it contributes to
      !! the molecule**, so subtracting these from atom-in-molecule energies
      !! leaves no basis-set superposition error to correct. A ghost centre
      !! carries no electrons and gets zero.
      ! TODO(mqc): the cache miss-and-solve block below is a verbatim copy of
      ! the one in `build_atomic_guess`, so a change to the eviction or failure
      ! policy has to be made twice or the two drift apart.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), allocatable, intent(out) :: energies(:)   !! (n_atoms), hartree
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: atom_mol
      integer :: iatom, idx, z

      if (error%has_error()) return

      allocate (energies(mol%natm))
      energies = 0.0_dp

      do iatom = 1, mol%natm
         z = nint(mol%charges(iatom))
         if (z <= 0) cycle

         call mol%atom_subset(iatom, atom_mol, error)
         if (error%has_error()) return

         idx = cached_index(z, atom_mol)
         if (idx == 0) then
            if (n_cached >= MAX_CACHED) then
               call error%set(ERROR_VALIDATION, "free atom energies: solution cache "// &
                              "exhausted")
               call atom_mol%destroy()
               return
            end if
            n_cached = n_cached + 1
            call solve_free_atom(atom_mol, z, cache(n_cached), error)
            if (error%has_error()) then
               cache(n_cached)%valid = .false.
               n_cached = n_cached - 1
               call atom_mol%destroy()
               return
            end if
            idx = n_cached
         end if

         energies(iatom) = cache(idx)%energy
         call atom_mol%destroy()
      end do
   end subroutine free_atom_energies

   subroutine clear_atomic_cache()
      !! Drop every cached atomic solution
      !!
      !! For the tests, which count how many free-atom SCFs a call performs; a
      !! run needs it only if it wants the memory back, since the cache keys
      !! itself on the basis data and a changed basis simply misses.
      integer :: i

      do i = 1, MAX_CACHED
         if (allocated(cache(i)%d_alpha)) deallocate (cache(i)%d_alpha)
         if (allocated(cache(i)%d_beta)) deallocate (cache(i)%d_beta)
         if (allocated(cache(i)%d_half)) deallocate (cache(i)%d_half)
         if (allocated(cache(i)%bas)) deallocate (cache(i)%bas)
         if (allocated(cache(i)%env)) deallocate (cache(i)%env)
         cache(i)%valid = .false.
         cache(i)%atomic_number = 0
         cache(i)%n_ao = 0
      end do
      n_cached = 0
   end subroutine clear_atomic_cache

end module mqc_czt_atomic_guess
