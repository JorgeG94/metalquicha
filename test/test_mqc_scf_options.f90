!! The shared SCF options reach the backend's settings
module test_mqc_scf_options
   !! One concept used to be maintained as three hand-copied instances -- the
   !! options type, the `configure_*` that filled it, and the
   !! `settings%x = options%x` block in each method. Missing any one of them
   !! failed silently: `keywords.scf.allow_crap_scf` was read from the deck,
   !! passed the schema, reached the config, and was dropped before the
   !! backend, so a Kohn-Sham run hard-stopped while the deck said to keep
   !! going. `scf_options_t`, `configure_scf` and `apply_scf_settings` removed
   !! the three copies; this is what keeps them removed.
   !!
   !! **It tests the class, not the symptom.** A deck that fails to converge
   !! would exercise `allow_crap_scf` and nothing else, needs an SCF to run,
   !! and has no stable energy to compare. `apply_scf_settings` is a pure
   !! function of its inputs, so every shared field can be checked at once in
   !! milliseconds with no molecule: give each a value nothing defaults to,
   !! hand it over, and require it to arrive.
   !!
   !! **Adding a field to `scf_options_t` means adding it here.** A field that
   !! nothing asserts is a field that can go missing again -- this catches a
   !! copy that was dropped, not a copy that was never written, and the two
   !! look identical from the backend's side.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_method_config, only: scf_options_t, properties_config_t
   use mqc_cuest_iface, only: cuest_scf_settings_t, apply_scf_settings, &
                              apply_properties_settings
   use mqc_scf_types, only: guess_step_t
   implicit none
   private

   public :: collect_mqc_scf_options_tests

contains

   subroutine collect_mqc_scf_options_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("every_shared_field_reaches_the_settings", every_field_arrives), &
                  new_unittest("an_unallocated_guess_ladder_is_not_copied", ladder_optional), &
                  new_unittest("every_property_reaches_the_settings", every_property_arrives), &
                  new_unittest("unset_allocatable_properties_are_not_copied", properties_optional) &
                  ]
   end subroutine collect_mqc_scf_options_tests

   !> An options object with nothing left at its default, so a field that fails
   !> to arrive shows up as the default it should no longer hold.
   subroutine populate(options)
      type(scf_options_t), intent(out) :: options

      options%basis_set = "cc-pvdz"
      options%ecp_set = "def2-ecp"
      options%aux_basis_set = "cc-pvdz-jkfit"
      options%aux_basis_named = .true.
      options%cartesian = .true.
      options%spherical = .false.
      options%density_fitting = .true.
      options%freeze_core = .true.
      options%n_frozen_core = 3
      options%unrestricted = .true.
      options%guess = "gwh"
      options%max_iter = 217
      options%allow_crap_scf = .true.
      options%energy_tol = 1.5e-9_dp
      options%density_tol = 2.5e-7_dp
      options%grad_tol = 3.5e-6_dp
      options%level_shift = 0.35_dp
      options%linear_dependence = 1.0e-5_dp
      options%use_diis = .false.
      options%diis_size = 11
      ! Set to the non-default, so a copy that never happens is caught. This
      ! one was added after the deck value reached `mqc_config_t`, was carried
      ! through every settings type, and then was not passed at the SCF call
      ! site -- a flag that parsed, validated and did nothing.
      options%incremental_fock = .false.
      options%accelerator = "ediis"
      options%verbose = .true.
      options%device_rank = 2
      options%pcm%enabled = .true.
      options%pcm%dielectric = 32.7_dp
      allocate (options%guess_steps(2))
      options%guess_steps(1)%basis = "sto-3g"
      options%guess_steps(1)%maxiter = 7
      options%guess_steps(2)%basis = "6-31g"
      options%guess_steps(2)%tolerance = 1.0e-4_dp
   end subroutine populate

   subroutine every_field_arrives(error)
      type(error_type), allocatable, intent(out) :: error

      type(scf_options_t) :: options
      type(cuest_scf_settings_t) :: settings

      call populate(options)
      call apply_scf_settings(settings, options)

      call check(error, settings%basis_set == options%basis_set, "basis_set")
      if (allocated(error)) return
      call check(error, settings%ecp_set == options%ecp_set, "ecp_set")
      if (allocated(error)) return
      call check(error, settings%aux_basis_set == options%aux_basis_set, "aux_basis_set")
      if (allocated(error)) return
      call check(error, settings%aux_basis_named .eqv. options%aux_basis_named, "aux_basis_named")
      if (allocated(error)) return
      call check(error, settings%cartesian .eqv. options%cartesian, "cartesian")
      if (allocated(error)) return
      call check(error, settings%spherical .eqv. options%spherical, "spherical")
      if (allocated(error)) return
      call check(error, settings%density_fitting .eqv. options%density_fitting, "density_fitting")
      if (allocated(error)) return
      call check(error, settings%freeze_core .eqv. options%freeze_core, "freeze_core")
      if (allocated(error)) return
      call check(error, settings%n_frozen_core == options%n_frozen_core, "n_frozen_core")
      if (allocated(error)) return
      call check(error, settings%unrestricted .eqv. options%unrestricted, "unrestricted")
      if (allocated(error)) return
      call check(error, settings%guess == options%guess, "guess")
      if (allocated(error)) return
      call check(error, settings%max_iter == options%max_iter, "max_iter")
      if (allocated(error)) return
      ! The field the whole refactor exists for: honoured on one reference and
      ! silently dropped on the other, for as long as this block was written
      ! out by hand in each method.
      call check(error, settings%allow_crap_scf .eqv. options%allow_crap_scf, "allow_crap_scf")
      if (allocated(error)) return
      call check(error, settings%energy_tol == options%energy_tol, "energy_tol")
      if (allocated(error)) return
      call check(error, settings%density_tol == options%density_tol, "density_tol")
      if (allocated(error)) return
      call check(error, settings%grad_tol == options%grad_tol, "grad_tol")
      if (allocated(error)) return
      call check(error, settings%level_shift == options%level_shift, "level_shift")
      if (allocated(error)) return
      call check(error, settings%linear_dependence == options%linear_dependence, "linear_dependence")
      if (allocated(error)) return
      call check(error, settings%use_diis .eqv. options%use_diis, "use_diis")
      if (allocated(error)) return
      call check(error, settings%diis_size == options%diis_size, "diis_size")
      if (allocated(error)) return
      call check(error, settings%incremental_fock .eqv. options%incremental_fock, &
                 "incremental_fock")
      if (allocated(error)) return
      call check(error, settings%accelerator == options%accelerator, "accelerator")
      if (allocated(error)) return
      call check(error, settings%verbose .eqv. options%verbose, "verbose")
      if (allocated(error)) return
      call check(error, settings%device_rank == options%device_rank, "device_rank")
      if (allocated(error)) return
      call check(error, settings%pcm%enabled .eqv. options%pcm%enabled, "pcm%enabled")
      if (allocated(error)) return
      call check(error, settings%pcm%dielectric == options%pcm%dielectric, "pcm%dielectric")
      if (allocated(error)) return

      call check(error, allocated(settings%guess_steps), "guess_steps not allocated")
      if (allocated(error)) return
      call check(error, size(settings%guess_steps) == 2, "guess_steps size")
      if (allocated(error)) return
      call check(error, settings%guess_steps(1)%basis == "sto-3g", "guess_steps(1)%basis")
      if (allocated(error)) return
      call check(error, settings%guess_steps(2)%maxiter == options%guess_steps(2)%maxiter, &
                 "guess_steps(2)%maxiter")
   end subroutine every_field_arrives

   !> The ladder is optional, and copying an unallocated one would be a crash
   !> rather than a wrong answer -- worth pinning separately because it is the
   !> one field in the block guarded by a condition.
   subroutine ladder_optional(error)
      type(error_type), allocatable, intent(out) :: error

      type(scf_options_t) :: options
      type(cuest_scf_settings_t) :: settings

      options%max_iter = 42
      call apply_scf_settings(settings, options)
      call check(error, settings%max_iter == 42, "an unallocated ladder stopped the copy")
      if (allocated(error)) return
      call check(error,.not. allocated(settings%guess_steps), &
                 "a ladder appeared from nowhere")
   end subroutine ladder_optional

   !> The properties block had diverged three ways before it was shared: HF
   !> unpacked all eight, Kohn-Sham four and MCSCF six. Same populate-and-assert
   !> as the SCF fields, for the same reason.
   subroutine every_property_arrives(error)
      type(error_type), allocatable, intent(out) :: error

      type(properties_config_t) :: properties
      type(cuest_scf_settings_t) :: settings

      properties%bonding_analysis = "quao"
      properties%bonding_threshold = 0.125_dp
      properties%bonding_energy = .true.
      properties%bonding_no_sharing = .true.
      properties%bonding_no_sharing_ci = "project"
      properties%bonding_restrict_localization = .true.
      properties%fukui_population = "mulliken"
      properties%charges_scheme = "chelpg"
      ! The ion SCF group. None of the three fields this replaced was ever
      ! asserted here, despite `scf_options_t` promising that every field is --
      ! so all three could have stopped crossing and nothing would have said
      ! so. Setting each to something no default would produce.
      properties%fukui_scf%max_iter = 173
      properties%fukui_scf%energy_tol = 3.5e-7_dp
      properties%fukui_scf%density_tol = 2.5e-5_dp
      properties%fukui_scf%level_shift = 0.375_dp
      properties%fukui_scf%diis_size = 11
      properties%fukui_scf%accelerator = "ediis"
      properties%fukui_scf%incremental_fock = .false.
      properties%fukui_scf%inherit_scf = .false.

      call apply_properties_settings(settings, properties)

      call check(error, settings%bonding_analysis == properties%bonding_analysis, &
                 "bonding_analysis")
      if (allocated(error)) return
      call check(error, settings%bonding_threshold == properties%bonding_threshold, &
                 "bonding_threshold")
      if (allocated(error)) return
      ! The one that guards a refusal rather than a feature: Kohn-Sham never
      ! set it, so the check that stops a bonding-energy decomposition running
      ! under a continuum could not fire.
      call check(error, settings%bonding_energy .eqv. properties%bonding_energy, &
                 "bonding_energy")
      if (allocated(error)) return
      call check(error, settings%bonding_no_sharing .eqv. properties%bonding_no_sharing, &
                 "bonding_no_sharing")
      if (allocated(error)) return
      call check(error, settings%bonding_no_sharing_ci == properties%bonding_no_sharing_ci, &
                 "bonding_no_sharing_ci")
      if (allocated(error)) return
      call check(error, settings%bonding_restrict_localization .eqv. &
                 properties%bonding_restrict_localization, "bonding_restrict_localization")
      if (allocated(error)) return
      call check(error, allocated(settings%fukui_population), "fukui_population unset")
      if (allocated(error)) return
      call check(error, settings%fukui_population == "mulliken", "fukui_population")
      if (allocated(error)) return
      call check(error, allocated(settings%charges_scheme), "charges_scheme unset")
      if (allocated(error)) return
      call check(error, settings%charges_scheme == "chelpg", "charges_scheme")
      if (allocated(error)) return
      call check(error, settings%fukui_scf%max_iter == 173, "fukui_scf%max_iter")
      if (allocated(error)) return
      call check(error, settings%fukui_scf%energy_tol == 3.5e-7_dp, "fukui_scf%energy_tol")
      if (allocated(error)) return
      call check(error, settings%fukui_scf%density_tol == 2.5e-5_dp, "fukui_scf%density_tol")
      if (allocated(error)) return
      call check(error, settings%fukui_scf%level_shift == 0.375_dp, "fukui_scf%level_shift")
      if (allocated(error)) return
      call check(error, settings%fukui_scf%diis_size == 11, "fukui_scf%diis_size")
      if (allocated(error)) return
      call check(error, settings%fukui_scf%accelerator == "ediis", "fukui_scf%accelerator")
      if (allocated(error)) return
      call check(error,.not. settings%fukui_scf%incremental_fock, &
                 "fukui_scf%incremental_fock")
      if (allocated(error)) return
      call check(error,.not. settings%fukui_scf%inherit_scf, "fukui_scf%inherit_scf")
   end subroutine every_property_arrives

   !> Two of the eight are allocatable, and an unset one must not overwrite
   !> whatever the backend already holds.
   subroutine properties_optional(error)
      type(error_type), allocatable, intent(out) :: error

      type(properties_config_t) :: properties
      type(cuest_scf_settings_t) :: settings

      properties%bonding_analysis = "quao"
      call apply_properties_settings(settings, properties)
      call check(error, settings%bonding_analysis == "quao", &
                 "an unset allocatable stopped the copy")
      if (allocated(error)) return
      call check(error,.not. allocated(settings%fukui_population), &
                 "fukui_population appeared from nowhere")
      if (allocated(error)) return
      call check(error,.not. allocated(settings%charges_scheme), &
                 "charges_scheme appeared from nowhere")
   end subroutine properties_optional

end module test_mqc_scf_options

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_scf_options, only: collect_mqc_scf_options_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_scf_options", collect_mqc_scf_options_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
