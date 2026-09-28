!! `fragment_refusal`'s table, case by case
module test_mqc_fragment_capabilities
   !! One call per combination, no SCF: `fragment_capabilities` and
   !! `fragment_refusal` are pure functions of a `method_config_t` and a
   !! `fragment_needs_t`, so every branch of the table in
   !! `mqc_docs/source/developer_fragment_solver.rst` is milliseconds.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_fragment_capabilities, only: fragment_refusal, fragment_needs_t, &
                                        FRAGMENT_SCHEME_FMO, FRAGMENT_SCHEME_EE_MBE, &
                                        FRAGMENT_SCHEME_EFMO
   use mqc_method_config, only: method_config_t
   use mqc_method_types, only: METHOD_TYPE_HF, METHOD_TYPE_DFT, METHOD_TYPE_MP2, &
                               METHOD_TYPE_CCSD
   use mqc_calc_types, only: CALC_TYPE_GRADIENT, CALC_TYPE_HESSIAN
   implicit none
   private

   public :: collect_fragment_capabilities

contains

   subroutine collect_fragment_capabilities(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fmo_hf_energy_ok", test_fmo_hf_ok), &
                  new_unittest("fmo_hf_pieda_ok", test_fmo_hf_pieda_ok), &
                  new_unittest("efmo_mp2_pieda_refused", test_efmo_mp2_pieda_refused), &
                  new_unittest("fmo_dft_ok", test_fmo_dft_ok), &
                  new_unittest("fmo_double_hybrid_refused", test_fmo_double_hybrid_refused), &
                  new_unittest("fmo_dft_cut_ok", test_fmo_dft_cut_ok), &
                  new_unittest("eembe_dft_ok", test_eembe_dft_ok), &
                  new_unittest("fmo_dft_dispersion_refused", test_fmo_dft_dispersion_refused), &
                  new_unittest("fmo_dft_pieda_refused", test_fmo_dft_pieda_refused), &
                  new_unittest("fmo_mp2_refused", test_fmo_mp2_refused), &
                  new_unittest("fmo_hf_gradient_refused", test_fmo_gradient_refused), &
                  new_unittest("fmo_hf_hessian_refused", test_fmo_hessian_refused), &
                  new_unittest("fmo_hf_unrestricted_refused", test_fmo_unrestricted_refused), &
                  new_unittest("eembe_hf_gradient_refused", test_eembe_gradient_refused), &
                  new_unittest("efmo_hf_ok", test_efmo_hf_ok), &
                  new_unittest("efmo_hf_cut_ok", test_efmo_hf_cut_ok), &
                  new_unittest("efmo_mp2_ok", test_efmo_mp2_ok), &
                  new_unittest("efmo_ri_mp2_ok", test_efmo_ri_mp2_ok), &
                  new_unittest("efmo_scs_mp2_refused", test_efmo_scs_refused), &
                  new_unittest("efmo_mp2_cut_refused", test_efmo_mp2_cut_refused), &
                  new_unittest("efmo_dft_refused", test_efmo_dft_refused), &
                  new_unittest("efmo_ccsd_refused", test_efmo_ccsd_refused), &
                  new_unittest("efmo_hf_gradient_refused", test_efmo_gradient_refused), &
                  new_unittest("efmo_hf_unrestricted_refused", test_efmo_unrestricted_refused) &
                  ]
   end subroutine collect_fragment_capabilities

   subroutine test_fmo_hf_ok(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_HF
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_FMO, config, needs)), 0, &
                 "FMO+HF energy was refused")
   end subroutine test_fmo_hf_ok

   subroutine test_fmo_hf_pieda_ok(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_HF
      needs%pieda = .true.
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_FMO, config, needs)), 0, &
                 "FMO+HF with pieda was refused")
   end subroutine test_fmo_hf_pieda_ok

   subroutine test_efmo_mp2_pieda_refused(error)
      !! PIEDA is a property of the method: MP2 has no decomposition yet
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs
      character(len=:), allocatable :: why

      config%method_type = METHOD_TYPE_MP2
      needs%pieda = .true.
      why = fragment_refusal(FRAGMENT_SCHEME_EFMO, config, needs)
      call check(error, len(why) > 0, "EFMO+MP2 with pieda was not refused")
      if (allocated(error)) return
      call check(error, index(why, "keywords.fragmentation.pieda") > 0, &
                 "the refusal does not name the key to change")
   end subroutine test_efmo_mp2_pieda_refused

   subroutine test_fmo_dft_ok(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_DFT
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_FMO, config, needs)), 0, &
                 "FMO+DFT energy was refused")
   end subroutine test_fmo_dft_ok

   subroutine test_fmo_double_hybrid_refused(error)
      !! A double hybrid's PT2 part is not added under FMO, so it is refused
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs
      character(len=:), allocatable :: why

      config%method_type = METHOD_TYPE_DFT
      config%dft%functional = "b2plyp"
      why = fragment_refusal(FRAGMENT_SCHEME_FMO, config, needs)
      call check(error, len(why) > 0, "FMO+B2PLYP was not refused")
      if (allocated(error)) return
      call check(error, index(why, "model.functional") > 0, &
                 "the refusal does not name the key to change")
   end subroutine test_fmo_double_hybrid_refused

   subroutine test_fmo_dft_cut_ok(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_DFT
      needs%cut = .true.
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_FMO, config, needs)), 0, &
                 "FMO+DFT with a detached bond was refused")
   end subroutine test_fmo_dft_cut_ok

   subroutine test_eembe_dft_ok(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_DFT
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_EE_MBE, config, needs)), 0, &
                 "EE-MBE+DFT was refused")
   end subroutine test_eembe_dft_ok

   subroutine test_fmo_dft_dispersion_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs
      character(len=:), allocatable :: why

      config%method_type = METHOD_TYPE_DFT
      needs%dispersion = .true.
      why = fragment_refusal(FRAGMENT_SCHEME_FMO, config, needs)
      call check(error, len(why) > 0, "FMO+DFT with dispersion was not refused")
      if (allocated(error)) return
      call check(error, index(why, "keywords.dft.dispersion") > 0, &
                 "the refusal does not name the key to change")
   end subroutine test_fmo_dft_dispersion_refused

   subroutine test_fmo_dft_pieda_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs
      character(len=:), allocatable :: why

      config%method_type = METHOD_TYPE_DFT
      needs%pieda = .true.
      why = fragment_refusal(FRAGMENT_SCHEME_FMO, config, needs)
      call check(error, len(why) > 0, "FMO+DFT with pieda was not refused")
      if (allocated(error)) return
      call check(error, index(why, "keywords.fragmentation.pieda") > 0, &
                 "the refusal does not name the key to change")
   end subroutine test_fmo_dft_pieda_refused

   subroutine test_fmo_mp2_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_MP2
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_FMO, config, needs)) > 0, &
                 "FMO+MP2 was not refused")
   end subroutine test_fmo_mp2_refused

   subroutine test_fmo_gradient_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs
      character(len=:), allocatable :: why

      config%method_type = METHOD_TYPE_HF
      needs%calc_type = CALC_TYPE_GRADIENT
      why = fragment_refusal(FRAGMENT_SCHEME_FMO, config, needs)
      call check(error, len(why) > 0, "FMO+HF gradient was not refused")
      if (allocated(error)) return
      call check(error, index(why, "gradient") > 0, &
                 "the refusal does not name the driver that was asked for")
   end subroutine test_fmo_gradient_refused

   subroutine test_fmo_hessian_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs
      character(len=:), allocatable :: why

      config%method_type = METHOD_TYPE_HF
      needs%calc_type = CALC_TYPE_HESSIAN
      why = fragment_refusal(FRAGMENT_SCHEME_FMO, config, needs)
      call check(error, len(why) > 0, "FMO+HF hessian was not refused")
      if (allocated(error)) return
      call check(error, index(why, "hessian") > 0, &
                 "the refusal does not name the driver that was asked for")
   end subroutine test_fmo_hessian_refused

   subroutine test_fmo_unrestricted_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs
      character(len=:), allocatable :: why

      config%method_type = METHOD_TYPE_HF
      needs%unrestricted = .true.
      why = fragment_refusal(FRAGMENT_SCHEME_FMO, config, needs)
      call check(error, len(why) > 0, "FMO+HF unrestricted was not refused")
      if (allocated(error)) return
      call check(error, index(why, "model.unrestricted") > 0, &
                 "the refusal does not name the key to change")
   end subroutine test_fmo_unrestricted_refused

   subroutine test_eembe_gradient_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_HF
      needs%calc_type = CALC_TYPE_GRADIENT
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_EE_MBE, config, needs)) > 0, &
                 "EE-MBE+HF gradient was not refused")
   end subroutine test_eembe_gradient_refused

   subroutine test_efmo_hf_ok(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_HF
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_EFMO, config, needs)), 0, &
                 "EFMO+HF was refused")
   end subroutine test_efmo_hf_ok

   subroutine test_efmo_hf_cut_ok(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_HF
      needs%cut = .true.
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_EFMO, config, needs)), 0, &
                 "EFMO+HF with a cut was refused")
   end subroutine test_efmo_hf_cut_ok

   subroutine test_efmo_mp2_ok(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_MP2
      config%corr%use_df = .false.
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_EFMO, config, needs)), 0, &
                 "EFMO+MP2 was refused")
   end subroutine test_efmo_mp2_ok

   subroutine test_efmo_ri_mp2_ok(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_MP2
      config%corr%use_df = .true.
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_EFMO, config, needs)), 0, &
                 "EFMO+RI-MP2 was refused")
   end subroutine test_efmo_ri_mp2_ok

   subroutine test_efmo_scs_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_MP2
      config%corr%use_scs = .true.
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_EFMO, config, needs)) > 0, &
                 "EFMO+SCS-MP2 was not refused")
   end subroutine test_efmo_scs_refused

   subroutine test_efmo_mp2_cut_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_MP2
      needs%cut = .true.
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_EFMO, config, needs)) > 0, &
                 "EFMO+MP2 with a cut was not refused")
   end subroutine test_efmo_mp2_cut_refused

   subroutine test_efmo_dft_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_DFT
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_EFMO, config, needs)) > 0, &
                 "EFMO+DFT was not refused")
   end subroutine test_efmo_dft_refused

   subroutine test_efmo_ccsd_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_CCSD
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_EFMO, config, needs)) > 0, &
                 "EFMO+CCSD was not refused")
   end subroutine test_efmo_ccsd_refused

   subroutine test_efmo_gradient_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs
      character(len=:), allocatable :: why

      config%method_type = METHOD_TYPE_HF
      needs%calc_type = CALC_TYPE_GRADIENT
      why = fragment_refusal(FRAGMENT_SCHEME_EFMO, config, needs)
      call check(error, len(why) > 0, "EFMO+HF gradient was not refused")
      if (allocated(error)) return
      call check(error, index(why, "gradient") > 0, &
                 "the refusal does not name the driver that was asked for")
   end subroutine test_efmo_gradient_refused

   subroutine test_efmo_unrestricted_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(method_config_t) :: config
      type(fragment_needs_t) :: needs

      config%method_type = METHOD_TYPE_HF
      needs%unrestricted = .true.
      call check(error, len(fragment_refusal(FRAGMENT_SCHEME_EFMO, config, needs)) > 0, &
                 "EFMO+HF unrestricted was not refused")
   end subroutine test_efmo_unrestricted_refused

end module test_mqc_fragment_capabilities

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fragment_capabilities, only: collect_fragment_capabilities
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fragment_capabilities", collect_fragment_capabilities)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
