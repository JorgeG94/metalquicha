module test_mqc_optimizer_types
   !! The optimizer's vocabulary and defaults
   !!
   !! Deliberately tests `mqc_optimizer_types` and not the DL-FIND bridge:
   !! this module is compiled whether or not `MQC_ENABLE_DLFIND` is on, so
   !! these run in CI, where the optimizer backend is off by default. What is
   !! checked here is what a deck can say and what it means -- the part a
   !! misspelling has to be caught by.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_optimizer_types, only: optimizer_settings_t, &
                                  coordinates_from_string, coordinates_to_string, &
                                  algorithm_from_string, algorithm_to_string, &
                                  OPT_COORDS_UNKNOWN, OPT_COORDS_CARTESIAN, &
                                  OPT_COORDS_HDLC, OPT_COORDS_DLC, &
                                  OPT_ALGO_UNKNOWN, OPT_ALGO_SD, OPT_ALGO_CG, &
                                  OPT_ALGO_LBFGS, OPT_ALGO_PRFO
   use mqc_calc_types, only: calc_type_from_string, calc_type_to_string, CALC_TYPE_OPTIMIZE
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_optimizer_types_tests

contains

   subroutine collect_mqc_optimizer_types_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("coordinates_from_string", test_coordinates_from_string), &
                  new_unittest("coordinates_unknown", test_coordinates_unknown), &
                  new_unittest("coordinates_roundtrip", test_coordinates_roundtrip), &
                  new_unittest("algorithm_from_string", test_algorithm_from_string), &
                  new_unittest("algorithm_unknown", test_algorithm_unknown), &
                  new_unittest("algorithm_roundtrip", test_algorithm_roundtrip), &
                  new_unittest("case_insensitive", test_case_insensitive), &
                  new_unittest("settings_defaults", test_settings_defaults), &
                  new_unittest("optimize_driver_parses", test_optimize_driver_parses) &
                  ]
   end subroutine collect_mqc_optimizer_types_tests

   subroutine test_coordinates_from_string(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, coordinates_from_string("cartesian") == OPT_COORDS_CARTESIAN, &
                 "cartesian should map to OPT_COORDS_CARTESIAN")
      if (allocated(error)) return

      call check(error, coordinates_from_string("xyz") == OPT_COORDS_CARTESIAN, &
                 "xyz is a spelling of cartesian")
      if (allocated(error)) return

      call check(error, coordinates_from_string("hdlc") == OPT_COORDS_HDLC, &
                 "hdlc should map to OPT_COORDS_HDLC")
      if (allocated(error)) return

      call check(error, coordinates_from_string("dlc") == OPT_COORDS_DLC, &
                 "dlc should map to OPT_COORDS_DLC")
   end subroutine test_coordinates_from_string

   subroutine test_coordinates_unknown(error)
      !! An unrecognised spelling must be recognisably unrecognised
      !!
      !! This is the case that matters: falling back to Cartesians for a deck
      !! that asked for something else would optimize the right molecule in the
      !! wrong coordinate system and report success.
      type(error_type), allocatable, intent(out) :: error

      call check(error, coordinates_from_string("hdcl") == OPT_COORDS_UNKNOWN, &
                 "a typo must not silently become a valid coordinate system")
      if (allocated(error)) return

      call check(error, coordinates_from_string("") == OPT_COORDS_UNKNOWN, &
                 "an empty string is not a coordinate system")
   end subroutine test_coordinates_unknown

   subroutine test_coordinates_roundtrip(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, coordinates_from_string(coordinates_to_string(OPT_COORDS_HDLC)) &
                 == OPT_COORDS_HDLC, "hdlc should survive a round trip")
      if (allocated(error)) return

      call check(error, coordinates_from_string(coordinates_to_string(OPT_COORDS_DLC)) &
                 == OPT_COORDS_DLC, "dlc should survive a round trip")
      if (allocated(error)) return

      call check(error, coordinates_from_string(coordinates_to_string(OPT_COORDS_CARTESIAN)) &
                 == OPT_COORDS_CARTESIAN, "cartesian should survive a round trip")
   end subroutine test_coordinates_roundtrip

   subroutine test_algorithm_from_string(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, algorithm_from_string("lbfgs") == OPT_ALGO_LBFGS, &
                 "lbfgs should map to OPT_ALGO_LBFGS")
      if (allocated(error)) return

      call check(error, algorithm_from_string("l-bfgs") == OPT_ALGO_LBFGS, &
                 "l-bfgs is a spelling of lbfgs")
      if (allocated(error)) return

      call check(error, algorithm_from_string("cg") == OPT_ALGO_CG, &
                 "cg should map to OPT_ALGO_CG")
      if (allocated(error)) return

      call check(error, algorithm_from_string("sd") == OPT_ALGO_SD, &
                 "sd should map to OPT_ALGO_SD")
      if (allocated(error)) return

      call check(error, algorithm_from_string("prfo") == OPT_ALGO_PRFO, &
                 "prfo should map to OPT_ALGO_PRFO")
   end subroutine test_algorithm_from_string

   subroutine test_algorithm_unknown(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, algorithm_from_string("bfgs") == OPT_ALGO_UNKNOWN, &
                 "bfgs is not lbfgs and must not be taken for it")
      if (allocated(error)) return

      call check(error, algorithm_from_string("newton") == OPT_ALGO_UNKNOWN, &
                 "an algorithm that is not offered must be refused")
   end subroutine test_algorithm_unknown

   subroutine test_algorithm_roundtrip(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, algorithm_from_string(algorithm_to_string(OPT_ALGO_LBFGS)) &
                 == OPT_ALGO_LBFGS, "lbfgs should survive a round trip")
      if (allocated(error)) return

      call check(error, algorithm_from_string(algorithm_to_string(OPT_ALGO_CG)) &
                 == OPT_ALGO_CG, "cg should survive a round trip")
      if (allocated(error)) return

      call check(error, algorithm_from_string(algorithm_to_string(OPT_ALGO_SD)) &
                 == OPT_ALGO_SD, "sd should survive a round trip")
      if (allocated(error)) return

      call check(error, algorithm_from_string(algorithm_to_string(OPT_ALGO_PRFO)) &
                 == OPT_ALGO_PRFO, "prfo should survive a round trip")
   end subroutine test_algorithm_roundtrip

   subroutine test_case_insensitive(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, coordinates_from_string("HDLC") == OPT_COORDS_HDLC, &
                 "the coordinate system should not depend on case")
      if (allocated(error)) return

      call check(error, algorithm_from_string("LBFGS") == OPT_ALGO_LBFGS, &
                 "the algorithm should not depend on case")
      if (allocated(error)) return

      call check(error, coordinates_from_string("  Cartesian  ") == OPT_COORDS_CARTESIAN, &
                 "surrounding blanks should not matter")
   end subroutine test_case_insensitive

   subroutine test_settings_defaults(error)
      !! A deck that says nothing should still describe a runnable optimization
      type(error_type), allocatable, intent(out) :: error

      type(optimizer_settings_t) :: settings

      call check(error, settings%coordinates == OPT_COORDS_CARTESIAN, &
                 "the default coordinate system should be cartesian")
      if (allocated(error)) return

      call check(error, settings%algorithm == OPT_ALGO_LBFGS, &
                 "the default algorithm should be lbfgs")
      if (allocated(error)) return

      call check(error, settings%max_steps > 0, &
                 "a default step cap of zero would refuse every optimization")
      if (allocated(error)) return

      call check(error, settings%gradient_tolerance > 0.0_dp, &
                 "the gradient tolerance must be a real threshold, not a sentinel")
      if (allocated(error)) return

      ! Negative is the "leave the engine's default alone" sentinel, so these
      ! must stay negative or the engine's own choice is silently overridden
      ! with a number nobody picked.
      call check(error, settings%max_step < 0.0_dp, &
                 "max_step should default to the engine's own choice")
      if (allocated(error)) return

      call check(error, settings%lbfgs_memory < 0, &
                 "lbfgs_memory should default to the engine's own choice")
   end subroutine test_settings_defaults

   subroutine test_optimize_driver_parses(error)
      !! `"driver": "Optimize"` has to reach the optimizer
      type(error_type), allocatable, intent(out) :: error

      call check(error, calc_type_from_string("Optimize") == CALC_TYPE_OPTIMIZE, &
                 "Optimize should map to CALC_TYPE_OPTIMIZE")
      if (allocated(error)) return

      call check(error, calc_type_from_string("optimize") == CALC_TYPE_OPTIMIZE, &
                 "the driver should not depend on case")
      if (allocated(error)) return

      call check(error, calc_type_from_string("opt") == CALC_TYPE_OPTIMIZE, &
                 "opt is a spelling of optimize")
      if (allocated(error)) return

      call check(error, calc_type_to_string(CALC_TYPE_OPTIMIZE) == "optimize", &
                 "CALC_TYPE_OPTIMIZE should name itself")
   end subroutine test_optimize_driver_parses

end module test_mqc_optimizer_types

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_optimizer_types, only: collect_mqc_optimizer_types_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_optimizer_types", collect_mqc_optimizer_types_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
