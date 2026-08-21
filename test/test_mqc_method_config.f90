!! The method configuration: its defaults, and what it says it is set to
module test_mqc_method_config
   !! `method_config_t` is what every backend is handed, and it is assembled by
   !! the reader, the adapter and the C interface separately. Two things about
   !! it are worth pinning and neither is checked anywhere else:
   !!
   !!   * **`reset` restores the documented defaults.** The type declares its
   !!     defaults on the components and `reset` restates them, which is a
   !!     second copy: a default changed in one place and not the other means a
   !!     freshly built config and a reset one describe different calculations,
   !!     and every number that comes out is plausible either way. The test is
   !!     the comparison the code cannot make of itself -- a default config
   !!     against a modified-then-reset one, field by field.
   !!   * **The solvation summary reports what was configured.** `has_solvation`
   !!     decides whether a solvent reaches tblite at all, and the info lines
   !!     are the only place a run says which model it used. A gas-phase energy
   !!     reported as solvated is the failure this prevents, and it looks
   !!     exactly like a solvated energy.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_method_config, only: method_config_t
   use mqc_method_types, only: METHOD_TYPE_UNKNOWN, METHOD_TYPE_DFT
   implicit none
   private

   public :: collect_mqc_method_config_tests

contains

   subroutine collect_mqc_method_config_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("reset_restores_every_default", test_reset), &
                  new_unittest("no_solvent_means_no_solvation", test_no_solvation), &
                  new_unittest("configure_defaults_the_model_to_alpb", test_default_model), &
                  new_unittest("cpcm_reports_its_grid_and_scale", test_cpcm_lines), &
                  new_unittest("alpb_reports_its_optional_terms", test_alpb_lines) &
                  ]
   end subroutine collect_mqc_method_config_tests

   subroutine test_reset(error)
      type(error_type), allocatable, intent(out) :: error

      type(method_config_t) :: fresh, touched

      ! Move everything a run might: the method, the basis, the SCF, the DFT
      ! block, the correlation block, the active space, and an allocatable
      ! inside the MCSCF settings, which is the one `reset` has to deallocate
      ! rather than merely overwrite.
      touched%method_type = METHOD_TYPE_DFT
      touched%basis_set = "cc-pvqz"
      touched%verbose = .true.
      touched%use_spherical = .false.
      touched%scf%max_iter = 7
      touched%scf%energy_convergence = 1.0e-3_dp
      touched%scf%use_diis = .false.
      touched%dft%functional = "pbe0"
      touched%dft%radial_points = 5
      touched%corr%freeze_core = .false.
      touched%corr%use_scs = .true.
      touched%cc%include_triples = .true.
      touched%cc%max_iter = 3
      touched%mcscf%n_active_electrons = 8
      touched%mcscf%optimize_orbitals = .false.
      allocate (touched%mcscf%ormas_subspaces(2))
      touched%mcscf%ormas_subspaces = [1, 3]
      touched%xtb%solvent = "water"
      touched%xtb%cpcm_nang = 590

      call touched%reset()

      call check(error, touched%method_type == fresh%method_type .and. &
                 touched%method_type == METHOD_TYPE_UNKNOWN, "method_type")
      if (allocated(error)) return
      call check(error, trim(touched%basis_set) == trim(fresh%basis_set), "basis_set")
      if (allocated(error)) return
      call check(error, touched%verbose .eqv. fresh%verbose, "verbose")
      if (allocated(error)) return
      call check(error, touched%use_spherical .eqv. fresh%use_spherical, "use_spherical")
      if (allocated(error)) return
      call check(error, touched%scf%max_iter == fresh%scf%max_iter, "scf%max_iter")
      if (allocated(error)) return
      call check(error, abs(touched%scf%energy_convergence - fresh%scf%energy_convergence) &
                 < tiny(1.0_dp), "scf%energy_convergence")
      if (allocated(error)) return
      call check(error, touched%scf%use_diis .eqv. fresh%scf%use_diis, "scf%use_diis")
      if (allocated(error)) return
      call check(error, trim(touched%dft%functional) == trim(fresh%dft%functional), &
                 "dft%functional")
      if (allocated(error)) return
      call check(error, touched%dft%radial_points == fresh%dft%radial_points, &
                 "dft%radial_points")
      if (allocated(error)) return
      call check(error, touched%corr%freeze_core .eqv. fresh%corr%freeze_core, &
                 "corr%freeze_core")
      if (allocated(error)) return
      call check(error, touched%corr%use_scs .eqv. fresh%corr%use_scs, "corr%use_scs")
      if (allocated(error)) return
      call check(error, touched%cc%include_triples .eqv. fresh%cc%include_triples, &
                 "cc%include_triples")
      if (allocated(error)) return
      call check(error, touched%cc%max_iter == fresh%cc%max_iter, "cc%max_iter")
      if (allocated(error)) return
      call check(error, touched%mcscf%n_active_electrons == fresh%mcscf%n_active_electrons, &
                 "mcscf%n_active_electrons")
      if (allocated(error)) return
      call check(error, touched%mcscf%optimize_orbitals .eqv. fresh%mcscf%optimize_orbitals, &
                 "mcscf%optimize_orbitals")
      if (allocated(error)) return
      ! The allocatable: left allocated, a later run inherits an active-space
      ! partition nobody asked for and the CI it solves is not the one the deck
      ! describes.
      call check(error,.not. allocated(touched%mcscf%ormas_subspaces), &
                 "reset left an ORMAS partition behind")
      if (allocated(error)) return
      call check(error, len_trim(touched%xtb%solvent) == 0, "xtb%solvent")
      if (allocated(error)) return
      call check(error, touched%xtb%cpcm_nang == fresh%xtb%cpcm_nang, "xtb%cpcm_nang")
   end subroutine test_reset

   subroutine test_no_solvation(error)
      type(error_type), allocatable, intent(out) :: error

      type(method_config_t) :: config
      character(len=128) :: lines(4)
      integer :: n

      call check(error,.not. config%xtb%has_solvation(), &
                 "a fresh config claims to be solvated")
      if (allocated(error)) return

      call config%xtb%get_solvation_info(lines, n)
      call check(error, n == 0, "a gas-phase run reported solvation settings")
      if (allocated(error)) return

      ! A dielectric on its own is a solvation request too -- CPCM takes one
      ! without a named solvent, and reading only the name would silently drop
      ! the continuum.
      config%xtb%dielectric = 78.39_dp
      call check(error, config%xtb%has_solvation(), &
                 "a dielectric alone did not count as solvation")
   end subroutine test_no_solvation

   subroutine test_default_model(error)
      type(error_type), allocatable, intent(out) :: error

      type(method_config_t) :: config

      ! Naming a solvent and no model is the common case; it must resolve to a
      ! model rather than reaching tblite as an empty string.
      call config%xtb%configure(use_cds=.true., use_shift=.true., dielectric=-1.0_dp, &
                                cpcm_nang=110, cpcm_rscale=1.0_dp, solvent="water")
      call check(error, config%xtb%has_solvation(), "the solvent did not take")
      if (allocated(error)) return
      call check(error, trim(config%xtb%solvation_model) == "alpb", &
                 "an unnamed solvation model did not default to ALPB")
      if (allocated(error)) return

      ! And an explicit one is kept.
      call config%xtb%configure(use_cds=.true., use_shift=.true., dielectric=-1.0_dp, &
                                cpcm_nang=110, cpcm_rscale=1.0_dp, solvent="water", &
                                solvation_model="gbsa")
      call check(error, trim(config%xtb%solvation_model) == "gbsa", &
                 "an explicit solvation model was overwritten")
   end subroutine test_default_model

   subroutine test_cpcm_lines(error)
      type(error_type), allocatable, intent(out) :: error

      type(method_config_t) :: config
      character(len=128) :: lines(4)
      integer :: n

      call config%xtb%configure(use_cds=.true., use_shift=.true., dielectric=-1.0_dp, &
                                cpcm_nang=590, cpcm_rscale=1.2_dp, solvent="water", &
                                solvation_model="cpcm")

      call config%xtb%get_solvation_info(lines, n)
      call check(error, n == 3, "CPCM reported something other than three lines")
      if (allocated(error)) return
      call check(error, index(lines(1), "cpcm") > 0, "the model is not named")
      if (allocated(error)) return
      ! The grid and the radii scale change the cavity and therefore the
      ! energy, so a run that does not print them cannot be reproduced.
      call check(error, index(lines(2), "590") > 0, "the CPCM grid size is not reported")
      if (allocated(error)) return
      call check(error, index(lines(3), "1.2") > 0, "the radii scale is not reported")
   end subroutine test_cpcm_lines

   subroutine test_alpb_lines(error)
      type(error_type), allocatable, intent(out) :: error

      type(method_config_t) :: config
      character(len=128) :: lines(4)
      integer :: n, bare

      call config%xtb%configure(use_cds=.false., use_shift=.false., dielectric=-1.0_dp, &
                                cpcm_nang=110, cpcm_rscale=1.0_dp, solvent="water", &
                                solvation_model="alpb")
      call config%xtb%get_solvation_info(lines, bare)
      call check(error, bare == 1, "ALPB with no extras reported more than its model")
      if (allocated(error)) return

      ! Both extras add a line each, and both change the energy: a solvation
      ! energy with the CDS terms and one without are different numbers.
      call config%xtb%configure(use_cds=.true., use_shift=.true., dielectric=-1.0_dp, &
                                cpcm_nang=110, cpcm_rscale=1.0_dp, solvent="water", &
                                solvation_model="alpb")
      call config%xtb%get_solvation_info(lines, n)
      call check(error, n == bare + 2, &
                 "the CDS and state-shift terms are not reported")
      if (allocated(error)) return
      call check(error, index(lines(2), "CDS") > 0, "the CDS line is missing")
      if (allocated(error)) return
      call check(error, index(lines(3), "shift") > 0 .or. index(lines(3), "Shift") > 0, &
                 "the solution-state shift line is missing")
   end subroutine test_alpb_lines

end module test_mqc_method_config

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_method_config, only: collect_mqc_method_config_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_method_config", collect_mqc_method_config_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
