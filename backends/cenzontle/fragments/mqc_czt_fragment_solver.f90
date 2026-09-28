!! One fragment or n-mer, solved for whatever HF/MP2 method the deck asked for
module mqc_czt_fragment_solver
   !! The one fragment-solver interface FMO, EE-MBE and EFMO call: a fragment
   !! or n-mer, an optional embedding operator and Fock projector in, an energy
   !! and a density out. See `mqc_docs/source/developer_fragment_solver.rst`.
   !!
   !! Restricted Hartree-Fock, and MP2 or RI-MP2 on top of it, are the whole of
   !! what runs here; anything else is refused before this module is reached
   !! (`mqc_fragment_capabilities`), and arriving here regardless is an
   !! internal-consistency error rather than a silent Hartree-Fock answer.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_scf_types, only: scf_numerics_t
   use mqc_fock_projector, only: fock_projector_t
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mp2, only: mp2_result_t, run_czt_mp2, run_czt_ri_mp2
   use mqc_elements, only: core_orbital_count
   implicit none
   private

   public :: fragment_request_t
   public :: fragment_outcome_t
   public :: solve_fragment_method

   real(dp), parameter :: RETRY_LEVEL_SHIFT = 0.5_dp
      !! Hartree, the least shift the retry of `retry_level_shift` uses

   type :: fragment_request_t
      !! How one fragment or n-mer is to be solved
      real(dp), allocatable :: h_extra(:, :)
         !! The embedding operator `u`, over `mol`'s AOs. Unallocated is no
         !! embedding.
      type(fock_projector_t), allocatable :: projector
         !! A detached bond's frozen orbitals. Allocated only when the group
         !! holds one.
      integer, allocatable :: guess
         !! One of `SCF_GUESS_*`. Unallocated takes the backend's own default.
      real(dp), allocatable :: guess_density(:, :)
      type(scf_numerics_t) :: drive
         !! How the SCF is driven -- the `scf=` argument of `run_czt_rhf`.
      integer :: max_iter = 100
      real(dp) :: energy_tol = 1.0e-9_dp
      real(dp) :: density_tol = 1.0e-7_dp
      real(dp), allocatable :: grad_tol
         !! Passed as `grad_tol=` only when allocated.
      logical :: retry_level_shift = .false.
         !! One retry at `max(drive%level_shift, RETRY_LEVEL_SHIFT)` if the
         !! first attempt does not converge; the n-mer retry.
      logical :: verbose = .false.
   end type fragment_request_t

   type :: fragment_outcome_t
      !! What a fragment or n-mer solve leaves behind
      real(dp) :: energy = 0.0_dp
         !! `E`: the reference plus correlation, including `Tr(D u)`.
      real(dp) :: internal = 0.0_dp
         !! `E' = E - Tr(D_esp u)`; equals `energy` when there is no `h_extra`.
      real(dp) :: reference = 0.0_dp
         !! The SCF energy, including `Tr(D u)`.
      real(dp) :: correlation = 0.0_dp
      real(dp), allocatable :: density(:, :)
         !! `D_esp`: the reference SCF's total density.
      type(rhf_result_t) :: scf
         !! The reference SCF itself -- iterations, commutator, orbitals.
      logical :: converged = .false.
   end type fragment_outcome_t

contains

   subroutine solve_fragment_method(method, mol, nelec, real_z, request, outcome, error, &
                                    reference, aux, label)
      !! Solve one fragment or n-mer for `method`: Hartree-Fock, or MP2 or
      !! RI-MP2 on it
      !!
      !! An unconverged SCF is not an error here: `outcome%converged` says so
      !! and the caller decides.
      type(cuest_scf_settings_t), intent(in) :: method
         !! What runs: `functional == "" .and. .not. run_cc` is Hartree-Fock,
         !! plus MP2 or RI-MP2 when `run_mp2` is set. Anything else is an
         !! internal-consistency error.
      type(czt_molecule_t), intent(in) :: mol
      integer, intent(in) :: nelec
      integer, intent(in) :: real_z(:)
         !! Atomic numbers of the real (non-ghost) atoms, for the frozen-core
         !! count MP2 uses. Unread on a Hartree-Fock request.
      type(fragment_request_t), intent(in) :: request
      type(fragment_outcome_t), intent(out) :: outcome
      type(error_t), intent(inout) :: error
      type(rhf_result_t), intent(in), optional :: reference
         !! An SCF already converged on `mol`. Present skips the SCF and adds
         !! only correlation; absent runs `run_czt_rhf` here.
      type(czt_molecule_t), intent(in), optional :: aux
         !! Fitting basis for RI-MP2, built by the caller. Required only when
         !! `method%run_mp2 .and. method%corr_density_fitting`.
      character(len=*), intent(in), optional :: label
         !! Names the term in the retry's log line. "fragment" when absent.

      type(scf_numerics_t) :: drive
      character(len=:), allocatable :: retry_label
      integer :: attempt, max_attempt

      if (len_trim(method%functional) > 0 .or. method%run_cc) then
         call error%set(ERROR_VALIDATION, "fragment solver: only Hartree-Fock and an "// &
                        "MP2 correction on top of it run here; this request should "// &
                        "have been refused by fragment_refusal before it reached the "// &
                        "backend.")
         return
      end if

      if (present(reference)) then
         outcome%scf = reference
      else
         retry_label = "fragment"
         if (present(label)) retry_label = label
         drive = request%drive
         max_attempt = 1
         if (request%retry_level_shift) max_attempt = 2
         do attempt = 1, max_attempt
            ! An unallocated component of `request` arrives absent.
            call run_czt_rhf(mol, nelec, request%max_iter, request%energy_tol, &
                             request%density_tol, request%verbose, outcome%scf, error, &
                             scf=drive, guess=request%guess, &
                             guess_density=request%guess_density, h_extra=request%h_extra, &
                             projector=request%projector, grad_tol=request%grad_tol)
            if (error%has_error()) exit
            if (outcome%scf%converged .or. attempt == max_attempt) exit
            call logger%verbose("  fmo: "//retry_label//" did not converge; retrying "// &
                                "with a level shift")
            drive%level_shift = max(drive%level_shift, RETRY_LEVEL_SHIFT)
         end do
         if (error%has_error()) return
      end if

      outcome%reference = outcome%scf%energy
      outcome%density = outcome%scf%density
      outcome%converged = outcome%scf%converged

      outcome%correlation = 0.0_dp
      if (method%run_mp2) then
         call fragment_correlation(method, mol, nelec, real_z, outcome%scf, aux, &
                                   outcome%correlation, error)
         if (error%has_error()) return
      end if
      outcome%energy = outcome%reference + outcome%correlation

      if (allocated(request%h_extra)) then
         outcome%internal = outcome%reference - sum(outcome%scf%density*request%h_extra)
         outcome%internal = outcome%internal + outcome%correlation
      else
         outcome%internal = outcome%energy
      end if
   end subroutine solve_fragment_method

   subroutine fragment_correlation(method, mol, nelec, real_z, scf, aux, correlation, error)
      !! MP2 or RI-MP2 on `scf`'s orbitals, EFMO's rule for the frozen core
      type(cuest_scf_settings_t), intent(in) :: method
      type(czt_molecule_t), intent(in) :: mol
      integer, intent(in) :: nelec
      integer, intent(in) :: real_z(:)
      type(rhf_result_t), intent(in) :: scf
      type(czt_molecule_t), intent(in), optional :: aux
      real(dp), intent(out) :: correlation
      type(error_t), intent(inout) :: error

      type(mp2_result_t) :: mp2
      integer :: frozen

      correlation = 0.0_dp
      frozen = method%n_frozen_core
      if (frozen < 0) frozen = core_orbital_count(real_z)
      if (.not. method%freeze_core) frozen = 0

      if (method%corr_density_fitting) then
         if (.not. present(aux)) then
            call error%set(ERROR_VALIDATION, "fragment solver: a fitted correlation "// &
                           "needs an auxiliary basis. Set model.aux_basis, or ask for "// &
                           "'mp2' rather than 'ri-mp2'.")
            return
         end if
         call run_czt_ri_mp2(mol, aux, scf%orbitals, scf%orbital_energies, nelec/2, &
                             scf%energy, mp2, error, n_frozen=frozen)
      else
         call run_czt_mp2(mol, scf%orbitals, scf%orbital_energies, nelec/2, scf%energy, &
                          mp2, error, n_frozen=frozen)
      end if
      if (error%has_error()) return

      ! One and one unless `method` scaled them, so a plain MP2 request is
      ! `mp2%same_spin + mp2%opposite_spin` to the bit: multiplying by exactly
      ! 1.0_dp changes nothing.
      correlation = method%scs_ss*mp2%same_spin + method%scs_os*mp2%opposite_spin
   end subroutine fragment_correlation

end module mqc_czt_fragment_solver
