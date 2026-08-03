!! Superposition-of-atomic-coefficients (SAC) initial guess
module mqc_cuest_atomic_guess
   !! Builds a molecular starting guess from converged isolated-atom solutions.
   !!
   !! Each distinct element is solved once as a free atom -- unrestricted, at
   !! its Hund's-rule ground-state multiplicity, in the same orbital basis,
   !! auxiliary basis, functional and grid as the parent calculation -- and the
   !! occupied coefficients are dropped into the molecular AO block belonging
   !! to each atom of that element.
   !!
   !! Coefficients rather than densities, because cuEST's density-fitted
   !! exchange only accepts occupied MO coefficients. A density guess would
   !! need a Fock build that skips exchange on the first iteration.
   !!
   !! Why it matters here: a core-Hamiltonian guess gets the sigma/pi ordering
   !! wrong in radicals, and since a wrong occupation can be perfectly
   !! self-consistent, the SCF then converges tidily onto an excited state. OH
   !! lands ~4 eV high that way, which is its A2Sigma+ <- X2Pi gap almost
   !! exactly.
   !!
   !! Atomic solutions are cached for the lifetime of the process, keyed by the
   !! settings they were converged under. That matters more than it looks: a
   !! Hessian is 6N SCFs on the same elements, and a fragmented run is
   !! thousands.
   use, intrinsic :: iso_c_binding, only: c_int64_t
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cgto, only: molecular_basis_type
   use mqc_cuest_context, only: cuest_context_t
   use mqc_cuest_integrals, only: cuest_system_t
   use mqc_cuest_scf, only: scf_result_t, run_uks_scf, SCF_GUESS_GWH
   implicit none
   private

   public :: build_sac_guess       !! Molecular guess coefficients from atomic solutions
   public :: clear_atomic_cache    !! Drop cached atomic solutions
   public :: hund_multiplicity     !! Ground-state multiplicity of a free atom
   public :: atom_ao_counts        !! AO functions contributed by each atom

   integer, parameter :: MAX_CACHED = 64

   type :: atomic_solution_t
      !! One converged free-atom calculation
      integer :: atomic_number = 0
      integer :: n_ao = 0
      integer :: n_alpha = 0, n_beta = 0
      real(dp), allocatable :: c_alpha(:, :), c_beta(:, :)
      logical :: valid = .false.
   end type atomic_solution_t

   type(atomic_solution_t), save :: cache(MAX_CACHED)
   integer, save :: n_cached = 0
   character(len=64), save :: cache_basis = ''
   integer, save :: cache_functional = -99
   logical, save :: cache_spherical = .true.

contains

   pure function hund_multiplicity(atomic_number) result(multiplicity)
      !! Ground-state spin multiplicity of a neutral free atom
      !!
      !! Aufbau filling in Madelung order, with Hund's first rule inside each
      !! subshell: electrons occupy distinct orbitals with parallel spin before
      !! pairing. Good enough for a guess; it gets the transition metals
      !! nominally wrong where the real ground state defies Madelung (Cr, Cu),
      !! which costs iterations rather than correctness.
      integer, intent(in) :: atomic_number
      integer :: multiplicity

      ! Subshell capacities in Madelung (n+l, then n) order: 1s 2s 2p 3s 3p 4s
      ! 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s
      integer, parameter :: CAPACITY(16) = [2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2]
      integer :: remaining, i, in_shell, degeneracy, unpaired

      remaining = atomic_number
      unpaired = 0
      do i = 1, size(CAPACITY)
         if (remaining <= 0) exit
         degeneracy = CAPACITY(i)/2
         in_shell = min(remaining, CAPACITY(i))
         remaining = remaining - in_shell
         ! Hund: singly occupy all orbitals first, then pair up.
         if (in_shell <= degeneracy) then
            unpaired = in_shell
         else
            unpaired = 2*degeneracy - in_shell
         end if
      end do
      multiplicity = unpaired + 1
   end function hund_multiplicity

   pure subroutine atom_ao_counts(mol_basis, use_spherical, counts)
      !! Number of AO basis functions contributed by each atom
      !!
      !! cuEST builds its AO basis from the per-atom shell list in atom order,
      !! so molecular AO indices are grouped by atom in the same order. That is
      !! what lets an atomic block be copied into a contiguous row range.
      type(molecular_basis_type), intent(in) :: mol_basis
      logical, intent(in) :: use_spherical
      integer, intent(out) :: counts(:)

      integer :: iatom, ishell, l

      counts = 0
      do iatom = 1, mol_basis%nelements
         do ishell = 1, mol_basis%elements(iatom)%nshells
            l = mol_basis%elements(iatom)%shells(ishell)%ang_mom
            if (use_spherical) then
               counts(iatom) = counts(iatom) + (2*l + 1)
            else
               counts(iatom) = counts(iatom) + (l + 1)*(l + 2)/2
            end if
         end do
      end do
   end subroutine atom_ao_counts

   subroutine clear_atomic_cache()
      !! Drop every cached atomic solution
      integer :: i

      do i = 1, n_cached
         if (allocated(cache(i)%c_alpha)) deallocate (cache(i)%c_alpha)
         if (allocated(cache(i)%c_beta)) deallocate (cache(i)%c_beta)
         cache(i)%valid = .false.
      end do
      n_cached = 0
      cache_basis = ''
      cache_functional = -99
   end subroutine clear_atomic_cache

   subroutine invalidate_if_settings_changed(basis_name, functional_id, use_spherical)
      !! An atomic solution is only reusable under the settings it was made with
      character(len=*), intent(in) :: basis_name
      integer, intent(in) :: functional_id
      logical, intent(in) :: use_spherical

      if (trim(cache_basis) /= trim(basis_name) .or. cache_functional /= functional_id &
          .or. (cache_spherical .neqv. use_spherical)) then
         call clear_atomic_cache()
         cache_basis = basis_name
         cache_functional = functional_id
         cache_spherical = use_spherical
      end if
   end subroutine invalidate_if_settings_changed

   function cached_index(atomic_number) result(idx)
      !! Position of a cached solution, or 0
      integer, intent(in) :: atomic_number
      integer :: idx, i

      idx = 0
      do i = 1, n_cached
         if (cache(i)%valid .and. cache(i)%atomic_number == atomic_number) then
            idx = i
            return
         end if
      end do
   end function cached_index

   subroutine solve_free_atom(context, atomic_number, mol_basis, aux_basis, use_spherical, &
                              functional_id, n_radial, n_angular, solution, error)
      !! Converge one isolated atom and keep its occupied coefficients
      !!
      !! The atomic SCF uses the GWH guess, never SAC, so there is no recursion
      !! back into this routine.
      type(cuest_context_t), intent(inout) :: context
      integer, intent(in) :: atomic_number
      type(molecular_basis_type), intent(in) :: mol_basis   !! One-atom basis
      type(molecular_basis_type), intent(in) :: aux_basis   !! One-atom auxiliary
      logical, intent(in) :: use_spherical
      integer, intent(in) :: functional_id, n_radial, n_angular
      type(atomic_solution_t), intent(out) :: solution
      type(error_t), intent(inout) :: error

      type(cuest_system_t) :: system
      type(scf_result_t) :: scf
      integer :: numbers(1), multiplicity, n_alpha, n_beta
      real(dp) :: coordinates(3, 1)

      numbers(1) = atomic_number
      coordinates = 0.0_dp
      multiplicity = hund_multiplicity(atomic_number)
      n_beta = (atomic_number - (multiplicity - 1))/2
      n_alpha = n_beta + multiplicity - 1

      call system%create(context, numbers, coordinates, mol_basis, aux_basis, &
                         use_spherical, n_alpha, functional_id, n_radial, n_angular, &
                         error, n_occ_beta=n_beta)
      if (.not. error%has_error()) then
         call run_uks_scf(system, numbers, coordinates, atomic_number, multiplicity, &
                          200, 1.0e-8_dp, 1.0e-6_dp, .true., 8, .false., scf, error, &
                          guess=SCF_GUESS_GWH)
      end if
      call system%destroy()

      if (error%has_error()) then
         call error%add_context("atomic guess: free atom Z="//trim(int_to_string(atomic_number)))
         return
      end if

      solution%atomic_number = atomic_number
      solution%n_ao = size(scf%occupied, 1)
      solution%n_alpha = scf%n_occupied
      solution%n_beta = scf%n_occupied_beta
      solution%c_alpha = scf%occupied
      solution%c_beta = scf%occupied_beta
      solution%valid = .true.
   end subroutine solve_free_atom

   pure function int_to_string(value) result(text)
      !! Minimal integer formatting, to keep error messages readable
      integer, intent(in) :: value
      character(len=12) :: text

      write (text, '(I0)') value
   end function int_to_string

   subroutine build_sac_guess(context, atomic_numbers, mol_basis, aux_basis, use_spherical, &
                              functional_id, n_radial, n_angular, atom_bases, atom_aux_bases, &
                              guess_alpha, guess_beta, error)
      !! Assemble molecular guess coefficients from free-atom solutions
      !!
      !! Column count is the sum of the atomic occupations, which generally
      !! differs from the molecule's own alpha/beta counts -- these coefficients
      !! are only ever used to build a guess Fock matrix, never occupied
      !! directly.
      type(cuest_context_t), intent(inout) :: context
      integer, intent(in) :: atomic_numbers(:)
      type(molecular_basis_type), intent(in) :: mol_basis, aux_basis
      logical, intent(in) :: use_spherical
      integer, intent(in) :: functional_id, n_radial, n_angular
      type(molecular_basis_type), intent(in) :: atom_bases(:)      !! One-atom orbital basis per atom
      type(molecular_basis_type), intent(in) :: atom_aux_bases(:)  !! One-atom auxiliary per atom
      real(dp), allocatable, intent(out) :: guess_alpha(:, :), guess_beta(:, :)
      type(error_t), intent(inout) :: error

      integer, allocatable :: counts(:), offsets(:)
      integer :: n_atoms, iatom, idx, n_ao, col_a, col_b, total_a, total_b, i

      n_atoms = size(atomic_numbers)
      allocate (counts(n_atoms), offsets(n_atoms))
      call atom_ao_counts(mol_basis, use_spherical, counts)

      offsets(1) = 0
      do iatom = 2, n_atoms
         offsets(iatom) = offsets(iatom - 1) + counts(iatom - 1)
      end do
      n_ao = sum(counts)

      ! ---- converge each distinct element once ------------------------------
      call invalidate_if_settings_changed(cache_basis, functional_id, use_spherical)
      total_a = 0
      total_b = 0
      do iatom = 1, n_atoms
         idx = cached_index(atomic_numbers(iatom))
         if (idx == 0) then
            if (n_cached >= MAX_CACHED) then
               call error%set(ERROR_VALIDATION, "atomic guess: cache exhausted")
               return
            end if
            n_cached = n_cached + 1
            call solve_free_atom(context, atomic_numbers(iatom), atom_bases(iatom), &
                                 atom_aux_bases(iatom), use_spherical, functional_id, &
                                 n_radial, n_angular, cache(n_cached), error)
            if (error%has_error()) return
            idx = n_cached
         end if
         total_a = total_a + cache(idx)%n_alpha
         total_b = total_b + cache(idx)%n_beta
      end do

      ! ---- drop each atomic block into its own rows -------------------------
      allocate (guess_alpha(n_ao, max(total_a, 1)), guess_beta(n_ao, max(total_b, 1)))
      guess_alpha = 0.0_dp
      guess_beta = 0.0_dp

      col_a = 0
      col_b = 0
      do iatom = 1, n_atoms
         idx = cached_index(atomic_numbers(iatom))
         if (idx == 0) cycle
         if (cache(idx)%n_ao /= counts(iatom)) then
            call error%set(ERROR_VALIDATION, &
                           "atomic guess: free-atom AO count does not match its molecular block")
            return
         end if
         do i = 1, cache(idx)%n_alpha
            col_a = col_a + 1
            guess_alpha(offsets(iatom) + 1:offsets(iatom) + counts(iatom), col_a) = &
               cache(idx)%c_alpha(:, i)
         end do
         do i = 1, cache(idx)%n_beta
            col_b = col_b + 1
            guess_beta(offsets(iatom) + 1:offsets(iatom) + counts(iatom), col_b) = &
               cache(idx)%c_beta(:, i)
         end do
      end do

      deallocate (counts, offsets)
   end subroutine build_sac_guess

end module mqc_cuest_atomic_guess
