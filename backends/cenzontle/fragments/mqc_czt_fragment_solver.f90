!! One fragment or n-mer, solved for whatever HF/DFT/MP2 method the deck asked for
module mqc_czt_fragment_solver
   !! The one fragment-solver interface FMO, EE-MBE and EFMO call: a fragment
   !! or n-mer, an optional embedding operator and Fock projector in, an energy
   !! and a density out. See `mqc_docs/source/developer_fragment_solver.rst`.
   !!
   !! Restricted Hartree-Fock, restricted Kohn-Sham, and MP2 or RI-MP2 on top
   !! of Hartree-Fock, are the whole of what runs here; anything else is
   !! refused before this module is reached (`mqc_fragment_capabilities`), and
   !! arriving here regardless is an internal-consistency error rather than a
   !! silent Hartree-Fock answer.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_scf_types, only: scf_numerics_t
   use mqc_fock_projector, only: fock_projector_t
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_czt_mp2, only: mp2_result_t, run_czt_mp2, run_czt_ri_mp2
   use mqc_elements, only: core_orbital_count
   use mqc_dispersion_apply, only: dispersion_apply
   implicit none
   private

   public :: fragment_request_t
   public :: fragment_outcome_t
   public :: solve_fragment_method
   public :: fragment_xc_context

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
      character(len=16) :: dispersion = "none"
         !! Empirical dispersion to add to this group's energy: `"none"`, or
         !! `"d3bj"` or `"d4"` as `keywords.dft.dispersion` spells them. It
         !! is evaluated on the group's own real atoms, with the damping
         !! parameters of `method%functional` (`"hf"` when that is empty), and
         !! is what GAMESS's `DFTDSM` does for each fragment and n-mer.
      real(dp), allocatable :: dispersion_xyz(:, :)
         !! (3, size(real_z)) Bohr: the real atoms `dispersion` is evaluated
         !! on, in the order of `real_z`. Never a ghost centre. Required when
         !! `dispersion` is not `"none"`.
      real(dp) :: dispersion_charge = 0.0_dp
         !! The group's total charge as `dispersion` sees it. Only D4 reads it.
   end type fragment_request_t

   type :: fragment_outcome_t
      !! What a fragment or n-mer solve leaves behind
      real(dp) :: energy = 0.0_dp
         !! `E`: the reference plus correlation and dispersion, including
         !! `Tr(D u)`.
      real(dp) :: internal = 0.0_dp
         !! `E' = E - Tr(D_esp u)`; equals `energy` when there is no `h_extra`.
         !! The dispersion does not depend on the density, so it is in both.
      real(dp) :: reference = 0.0_dp
         !! The SCF energy, including `Tr(D u)`.
      real(dp) :: correlation = 0.0_dp
      real(dp) :: dispersion = 0.0_dp
         !! The empirical dispersion inside `energy` and `internal`, from
         !! `fragment_request_t%dispersion`. Zero when none was asked for.
      real(dp), allocatable :: density(:, :)
         !! `D_esp`: the reference SCF's total density.
      type(rhf_result_t) :: scf
         !! The reference SCF itself -- iterations, commutator, orbitals.
      logical :: converged = .false.
   end type fragment_outcome_t

contains

   subroutine solve_fragment_method(method, mol, nelec, real_z, request, outcome, error, &
                                    reference, aux, label)
      !! Solve one fragment or n-mer for `method`: Hartree-Fock or Kohn-Sham,
      !! plus MP2 or RI-MP2 on top of a Hartree-Fock reference
      !!
      !! An unconverged SCF is not an error here: `outcome%converged` says so
      !! and the caller decides.
      !!
      !! For either reference, `outcome%internal` is `outcome%energy` less
      !! `Tr(D u)` with the SCF's own density.
      !!
      !! With `request%dispersion` set, the dispersion of the group's own
      !! real atoms is added once, to `outcome%energy` and `outcome%internal`,
      !! and reported apart in `outcome%dispersion`.
      type(cuest_scf_settings_t), intent(in) :: method
         !! What runs: an empty `functional` is Hartree-Fock and a non-empty
         !! one is restricted Kohn-Sham, plus MP2 or RI-MP2 on a Hartree-Fock
         !! reference when `run_mp2` is set. `run_cc` is an
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
      type(xc_context_t) :: xc
      character(len=:), allocatable :: retry_label
      integer :: attempt, max_attempt
      logical :: kohn_sham

      if (method%run_cc) then
         call error%set(ERROR_VALIDATION, "fragment solver: coupled-cluster fragments "// &
                        "are not implemented; this request should have been refused "// &
                        "by fragment_refusal before it reached the backend.")
         return
      end if

      kohn_sham = len_trim(method%functional) > 0

      ! Before the SCF, as the unfragmented path does it: it costs
      ! microseconds and depends on the nuclei alone, so a functional with no
      ! damping parameters is refused before a fragment is solved.
      if (trim(request%dispersion) /= "none") then
         call fragment_dispersion(method, real_z, request, outcome%dispersion, error)
         if (error%has_error()) return
      end if

      if (present(reference)) then
         outcome%scf = reference
      else
         retry_label = "fragment"
         if (present(label)) retry_label = label
         drive = request%drive
         max_attempt = 1
         if (request%retry_level_shift) max_attempt = 2

         if (kohn_sham) then
            call fragment_xc_context(method, mol, xc, error)
            if (error%has_error()) return
         end if

         do attempt = 1, max_attempt
            ! An unallocated component of `request` arrives absent. `xc`
            ! arrives absent for Hartree-Fock too: it is passed only in the
            ! Kohn-Sham branch, so a Hartree-Fock fragment calls `run_czt_rhf`
            ! exactly as it did before this branch existed.
            if (kohn_sham) then
               call run_czt_rhf(mol, nelec, request%max_iter, request%energy_tol, &
                                request%density_tol, request%verbose, outcome%scf, error, &
                                scf=drive, guess=request%guess, &
                                guess_density=request%guess_density, xc=xc, &
                                h_extra=request%h_extra, projector=request%projector, &
                                grad_tol=request%grad_tol)
            else
               call run_czt_rhf(mol, nelec, request%max_iter, request%energy_tol, &
                                request%density_tol, request%verbose, outcome%scf, error, &
                                scf=drive, guess=request%guess, &
                                guess_density=request%guess_density, h_extra=request%h_extra, &
                                projector=request%projector, grad_tol=request%grad_tol)
            end if
            if (error%has_error()) exit
            if (outcome%scf%converged .or. attempt == max_attempt) exit
            call logger%verbose("  fmo: "//retry_label//" did not converge; retrying "// &
                                "with a level shift")
            drive%level_shift = max(drive%level_shift, RETRY_LEVEL_SHIFT)
         end do
         if (kohn_sham) call xc%destroy()
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

      ! Added last, and only when asked for, so a run without dispersion is
      ! arithmetically what it was.
      if (trim(request%dispersion) /= "none") then
         outcome%energy = outcome%energy + outcome%dispersion
         outcome%internal = outcome%internal + outcome%dispersion
      end if
   end subroutine solve_fragment_method

   subroutine fragment_dispersion(method, real_z, request, energy, error)
      !! The empirical dispersion of one group's real atoms, in Hartree
      !!
      !! The ordinary correction of the whole group as if it were a molecule of
      !! its own, at the damping parameters of the deck's functional.
      type(cuest_scf_settings_t), intent(in) :: method
      integer, intent(in) :: real_z(:)
      type(fragment_request_t), intent(in) :: request
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      type(error_t) :: derr
      character(len=:), allocatable :: functional

      energy = 0.0_dp
      if (.not. allocated(request%dispersion_xyz)) then
         call error%set(ERROR_VALIDATION, "fragment solver: dispersion '"// &
                        trim(request%dispersion)//"' was requested without the "// &
                        "coordinates of the atoms it is evaluated on.")
         return
      end if
      if (size(request%dispersion_xyz, 2) /= size(real_z)) then
         call error%set(ERROR_VALIDATION, "fragment solver: dispersion was given "// &
                        "coordinates for a different number of atoms than real_z.")
         return
      end if
      functional = "hf"
      if (len_trim(method%functional) > 0) functional = trim(method%functional)

      call dispersion_apply(trim(request%dispersion), functional, request%dispersion_charge, &
                            real_z, request%dispersion_xyz, energy, error=derr)
      if (derr%has_error()) then
         call error%set(ERROR_VALIDATION, "fragment solver: "//derr%get_message())
      end if
   end subroutine fragment_dispersion

   subroutine fragment_xc_context(method, mol, xc, error)
      !! The Kohn-Sham functional `method` names, on `mol`'s grid
      !!
      !! Built with the functional, grid, non-local, screening and blocking
      !! settings of `method`, so the energy of any density evaluated with it
      !! is the one the SCF of a fragment solved for `method` minimises. A
      !! double hybrid is refused. The caller destroys `xc`.
      type(cuest_scf_settings_t), intent(in) :: method
      type(czt_molecule_t), intent(in) :: mol
      type(xc_context_t), intent(inout) :: xc
      type(error_t), intent(inout) :: error

      if (.not. xc_available()) then
         call error%set(ERROR_VALIDATION, "a functional was requested ('"// &
                        trim(method%functional)//"') but this build has no "// &
                        "libxc: configure with -DMQC_ENABLE_LIBXC=ON")
         return
      end if
      call xc_context_create(mol, trim(method%functional), xc, error, &
                             level=method%grid_level, polarized=.false., &
                             nlc_level=method%nlc_grid_level, &
                             screen_tol=method%screening_tolerance, &
                             point_block=method%block_size, &
                             n_radial=method%radial_points, &
                             n_angular=method%angular_points)
      if (error%has_error()) return
      if (xc%pt2_fraction /= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "fragment solver: model.functional '"// &
                        trim(method%functional)//"' is a double hybrid, whose "// &
                        "perturbative correlation this solver does not add; this "// &
                        "request should have been refused by fragment_refusal.")
         call xc%destroy()
      end if
   end subroutine fragment_xc_context

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
