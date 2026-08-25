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
                                  OPT_ALGO_LBFGS, OPT_ALGO_PRFO, &
                                  OPT_ALGO_CG_AUTO, OPT_ALGO_NR, OPT_ALGO_DAMPED, &
                                  algorithm_needs_hessian, &
                                  hessian_update_from_string, &
                                  OPT_HESSIAN_UPDATE_ENGINE, OPT_HESSIAN_UPDATE_NONE, &
                                  OPT_HESSIAN_UPDATE_POWELL, OPT_HESSIAN_UPDATE_BOFILL, &
                                  constraint_from_string, constraint_atom_count, &
                                  target_from_string, target_to_string, &
                                  algorithm_finds_saddle, &
                                  OPT_TARGET_MINIMUM, OPT_TARGET_SADDLE, &
                                  neb_ends_from_string, neb_ends_to_string, &
                                  NEB_ENDS_FROZEN, NEB_ENDS_PERPENDICULAR, &
                                  NEB_ENDS_FREE, MIN_NEB_IMAGES
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
                  new_unittest("optimize_driver_parses", test_optimize_driver_parses), &
                  new_unittest("new_algorithms_parse", test_new_algorithms), &
                  new_unittest("hessian_update_parses", test_hessian_update), &
                  new_unittest("constraint_atom_counts", test_constraint_atoms), &
                  new_unittest("target_parses", test_target_parses), &
                  new_unittest("target_roundtrip", test_target_roundtrip), &
                  new_unittest("target_defaults_to_minimum", test_target_default), &
                  new_unittest("only_saddle_capable_algorithms", test_saddle_capable), &
                  new_unittest("neb_endpoints_parse", test_neb_ends), &
                  new_unittest("neb_defaults_are_off", test_neb_defaults) &
                  ]
   end subroutine collect_mqc_optimizer_types_tests

   subroutine test_neb_ends(error)
      !! How a deck spells the endpoint treatment
      type(error_type), allocatable, intent(out) :: error

      call check(error, neb_ends_from_string("frozen") == NEB_ENDS_FROZEN)
      if (allocated(error)) return
      call check(error, neb_ends_from_string("fixed") == NEB_ENDS_FROZEN)
      if (allocated(error)) return
      call check(error, neb_ends_from_string("perpendicular") == NEB_ENDS_PERPENDICULAR)
      if (allocated(error)) return
      call check(error, neb_ends_from_string("free") == NEB_ENDS_FREE)
      if (allocated(error)) return
      call check(error, neb_ends_from_string("FROZEN") == NEB_ENDS_FROZEN, &
                 "endpoint spellings should be case-insensitive like their siblings")
      if (allocated(error)) return
      call check(error, neb_ends_from_string("banana") < NEB_ENDS_FROZEN, &
                 "an unrecognised endpoint treatment must not fall back to a valid one")
      if (allocated(error)) return
      call check(error, neb_ends_from_string(neb_ends_to_string(NEB_ENDS_PERPENDICULAR)) == &
                 NEB_ENDS_PERPENDICULAR, "the endpoint treatment should round-trip")
   end subroutine test_neb_ends

   subroutine test_neb_defaults(error)
      !! A path is opt-in, and its floor is real
      !!
      !! The default matters more than it looks: every optimization in every
      !! existing deck runs through these settings, and a non-zero image count
      !! or an allocated endpoint would turn one of them into a band.
      type(error_type), allocatable, intent(out) :: error
      type(optimizer_settings_t) :: settings

      call check(error,.not. allocated(settings%endpoint), &
                 "an optimization is a single structure unless a deck says otherwise")
      if (allocated(error)) return
      call check(error, settings%n_images == 0, &
                 "zero images means the engine default, not a band")
      if (allocated(error)) return
      call check(error, settings%neb_spring < 0.0_dp, &
                 "a negative spring constant means the engine default")
      if (allocated(error)) return
      call check(error, settings%neb_ends == NEB_ENDS_FROZEN, &
                 "endpoints supplied by a user are usually already minima")
      if (allocated(error)) return
      ! Two images are the two endpoints. DL-FIND refuses fewer than three from
      ! inside the engine; this is the same floor, named where a deck is read.
      call check(error, MIN_NEB_IMAGES == 3)
   end subroutine test_neb_defaults

   subroutine test_target_parses(error)
      !! Every spelling a deck is likely to use, and one it is not
      type(error_type), allocatable, intent(out) :: error

      call check(error, target_from_string("minimum") == OPT_TARGET_MINIMUM)
      if (allocated(error)) return
      call check(error, target_from_string("min") == OPT_TARGET_MINIMUM)
      if (allocated(error)) return
      call check(error, target_from_string("saddle") == OPT_TARGET_SADDLE)
      if (allocated(error)) return
      call check(error, target_from_string("ts") == OPT_TARGET_SADDLE)
      if (allocated(error)) return
      call check(error, target_from_string("transition_state") == OPT_TARGET_SADDLE)
      if (allocated(error)) return
      call check(error, target_from_string("SADDLE") == OPT_TARGET_SADDLE, &
                 "the target should be case-insensitive like its siblings")
      if (allocated(error)) return
      ! Negative rather than a default. A misspelled target that quietly became
      ! "minimum" would run a completely different calculation from the one
      ! asked for and report success.
      call check(error, target_from_string("banana") < OPT_TARGET_MINIMUM, &
                 "an unrecognised target must not fall back to a valid one")
   end subroutine test_target_parses

   subroutine test_target_roundtrip(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, target_from_string(target_to_string(OPT_TARGET_MINIMUM)) == &
                 OPT_TARGET_MINIMUM)
      if (allocated(error)) return
      call check(error, target_from_string(target_to_string(OPT_TARGET_SADDLE)) == &
                 OPT_TARGET_SADDLE)
   end subroutine test_target_roundtrip

   subroutine test_target_default(error)
      !! An optimization nobody said anything about is a minimisation
      type(error_type), allocatable, intent(out) :: error
      type(optimizer_settings_t) :: settings

      call check(error, settings%target == OPT_TARGET_MINIMUM, &
                 "the default target must stay 'minimum'; every existing deck relies on it")
   end subroutine test_target_default

   subroutine test_saddle_capable(error)
      !! Which algorithms can be asked for a saddle
      !!
      !! The downhill methods are the point of this test. Asking steepest
      !! descent for a transition state is a deck that cannot be satisfied, and
      !! it has to be refused rather than run to a minimum and reported as one.
      type(error_type), allocatable, intent(out) :: error

      call check(error, algorithm_finds_saddle(OPT_ALGO_PRFO), &
                 "P-RFO is the saddle optimizer")
      if (allocated(error)) return
      call check(error, algorithm_finds_saddle(OPT_ALGO_NR), &
                 "Newton-Raphson seeks a stationary point of any curvature")
      if (allocated(error)) return
      call check(error,.not. algorithm_finds_saddle(OPT_ALGO_LBFGS))
      if (allocated(error)) return
      call check(error,.not. algorithm_finds_saddle(OPT_ALGO_SD))
      if (allocated(error)) return
      call check(error,.not. algorithm_finds_saddle(OPT_ALGO_CG))
      if (allocated(error)) return
      call check(error,.not. algorithm_finds_saddle(OPT_ALGO_DAMPED), &
                 "damped dynamics crosses barriers; it does not stop on them")
   end subroutine test_saddle_capable

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
      if (allocated(error)) return

      ! Off unless asked for. The check costs a whole Hessian, and a deck that
      ! did not ask for one should not find it in the bill.
      call check(error,.not. settings%hess_end, &
                 "hess_end should be off until a deck asks for it")
   end subroutine test_settings_defaults

   subroutine test_new_algorithms(error)
      !! The algorithms that were mapped but unreachable, and which need curvature
      !!
      !! `cg` and `cg-auto` are different DL-FIND algorithms and must not
      !! collapse onto one another -- they are the two restart policies for the
      !! same method, and a deck asking for one should not silently get the
      !! other.
      type(error_type), allocatable, intent(out) :: error

      call check(error, algorithm_from_string("prfo") == OPT_ALGO_PRFO, &
                 "prfo should parse")
      if (allocated(error)) return
      call check(error, algorithm_from_string("nr") == OPT_ALGO_NR, &
                 "nr should parse")
      if (allocated(error)) return
      call check(error, algorithm_from_string("newton-raphson") == OPT_ALGO_NR, &
                 "and its long spelling")
      if (allocated(error)) return
      call check(error, algorithm_from_string("damped") == OPT_ALGO_DAMPED, &
                 "damped should parse")
      if (allocated(error)) return

      call check(error, algorithm_from_string("cg-auto") == OPT_ALGO_CG_AUTO, &
                 "cg-auto should parse")
      if (allocated(error)) return
      call check(error, algorithm_from_string("cg") /= algorithm_from_string("cg-auto"), &
                 "the two conjugate gradients are different algorithms")
      if (allocated(error)) return

      ! Only these two hold a Hessian, and that is what decides whether the
      ! driver offers one at all.
      call check(error, algorithm_needs_hessian(OPT_ALGO_PRFO), &
                 "p-rfo needs curvature")
      if (allocated(error)) return
      call check(error, algorithm_needs_hessian(OPT_ALGO_NR), &
                 "so does newton-raphson")
      if (allocated(error)) return
      call check(error,.not. algorithm_needs_hessian(OPT_ALGO_LBFGS), &
                 "l-bfgs approximates its own and must not be handed one")
      if (allocated(error)) return
      call check(error,.not. algorithm_needs_hessian(OPT_ALGO_DAMPED), &
                 "damped dynamics integrates motion and holds no Hessian")
   end subroutine test_new_algorithms

   subroutine test_hessian_update(error)
      !! The update schemes, and that an unknown one is distinguishable
      !!
      !! "auto" and a misspelling must not both come back as the engine
      !! default: the first is a request to leave DL-FIND alone and the second
      !! is a mistake that should be reported.
      type(error_type), allocatable, intent(out) :: error

      call check(error, hessian_update_from_string("bofill") == OPT_HESSIAN_UPDATE_BOFILL, &
                 "bofill should parse")
      if (allocated(error)) return
      call check(error, hessian_update_from_string("powell") == OPT_HESSIAN_UPDATE_POWELL, &
                 "powell should parse")
      if (allocated(error)) return
      call check(error, hessian_update_from_string("none") == OPT_HESSIAN_UPDATE_NONE, &
                 "none should parse")
      if (allocated(error)) return
      call check(error, hessian_update_from_string("auto") == OPT_HESSIAN_UPDATE_ENGINE, &
                 "auto should mean the engine's own choice")
      if (allocated(error)) return
      call check(error, hessian_update_from_string("bofil") < OPT_HESSIAN_UPDATE_ENGINE, &
                 "a misspelling must be distinguishable from auto")
   end subroutine test_hessian_update

   subroutine test_constraint_atoms(error)
      !! Each constraint type over the number of atoms it is measured on
      !!
      !! This count is the only check on a constraint that can catch a wrong
      !! one: three atom indices under "torsion" is a well-formed list of
      !! integers and a meaningless constraint, and DL-FIND would read the
      !! fourth slot as atom zero.
      type(error_type), allocatable, intent(out) :: error

      call check(error, constraint_atom_count(constraint_from_string("bond")) == 2, &
                 "a bond is two atoms")
      if (allocated(error)) return
      call check(error, constraint_atom_count(constraint_from_string("angle")) == 3, &
                 "an angle is three")
      if (allocated(error)) return
      call check(error, constraint_atom_count(constraint_from_string("torsion")) == 4, &
                 "a torsion is four")
      if (allocated(error)) return
      call check(error, constraint_atom_count(constraint_from_string("cartesian")) == 1, &
                 "a cartesian constraint is one atom")
      if (allocated(error)) return
      call check(error, constraint_atom_count(constraint_from_string("dihedral")) == 4, &
                 "dihedral is a spelling of torsion")
      if (allocated(error)) return

      ! Zero is how an unknown type is caught, so it has to stay zero.
      call check(error, constraint_atom_count(constraint_from_string("wiggle")) == 0, &
                 "an unknown constraint type names no atoms")
   end subroutine test_constraint_atoms

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
