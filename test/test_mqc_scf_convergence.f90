!! The rule that stops an SCF, tested where it lives rather than through a run
module test_mqc_scf_convergence
   !! `scf_convergence_t` is a pure decision over three numbers, so it can be
   !! tested as one. That matters more than the convenience: the interesting
   !! cases are the ones an SCF cannot easily be made to produce on demand --
   !! an iterate whose energy has stopped moving while its commutator has not,
   !! which is what an interpolating accelerator does and what an energy-only
   !! test gets wrong.
   !!
   !! Driving that through a real SCF means finding a molecule, a basis and an
   !! accelerator that stall together, and the test then fails for a different
   !! reason the day any of the three changes. Here the stalled iterate is just
   !! a pair of numbers.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_scf_convergence, only: scf_convergence_t, parse_convergence_metric, &
                                  convergence_metric_name, CONV_METRIC_STANDARD, &
                                  CONV_METRIC_ENERGY, CONV_METRIC_COMMUTATOR, &
                                  CONV_METRIC_DENSITY, CONV_METRIC_UNKNOWN
   implicit none
   private

   public :: collect_scf_convergence

contains

   subroutine collect_scf_convergence(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("standard_is_the_default_and_needs_both", test_standard), &
                  new_unittest("standard_derives_the_commutator_bound", test_derived_bound), &
                  new_unittest("a_stated_gradient_tolerance_wins", test_stated_bound), &
                  new_unittest("single_metrics_read_only_their_own_number", test_single), &
                  new_unittest("energy_alone_accepts_a_stalled_iterate", test_stall), &
                  new_unittest("metric_spellings_and_aliases", test_parse), &
                  new_unittest("metric_names_round_trip", test_names) &
                  ]
   end subroutine collect_scf_convergence

   subroutine test_standard(error)
      !! The default is both gates, which is what every loop tested before
      !! this type existed. If this moves, introducing the type moved a number.
      type(error_type), allocatable, intent(out) :: error

      type(scf_convergence_t) :: conv

      call check(error, conv%metric == CONV_METRIC_STANDARD, "default metric")
      if (allocated(error)) return
      call check(error, conv%tolerance == 1.0e-8_dp, "default tolerance")
      if (allocated(error)) return

      ! Both met.
      call check(error, conv%is_converged(1.0e-9_dp, 1.0e-7_dp, 1.0e-5_dp), &
                 "both gates met must converge")
      if (allocated(error)) return
      ! Energy met, commutator not: sqrt(1e-8) is 1e-4.
      call check(error,.not. conv%is_converged(1.0e-9_dp, 1.0e-7_dp, 1.0e-3_dp), &
                 "a loose commutator must not converge under standard")
      if (allocated(error)) return
      ! Commutator met, energy not.
      call check(error,.not. conv%is_converged(1.0e-6_dp, 1.0e-7_dp, 1.0e-5_dp), &
                 "a loose energy must not converge under standard")
   end subroutine test_standard

   subroutine test_derived_bound(error)
      !! `sqrt(tolerance)`, and it is the number the table prints
      type(error_type), allocatable, intent(out) :: error

      type(scf_convergence_t) :: conv

      conv%tolerance = 1.0e-6_dp
      call check(error, abs(conv%commutator_bound() - 1.0e-3_dp) < 1.0e-15_dp, &
                 "sqrt(1e-6) should be 1e-3")
      if (allocated(error)) return
      call check(error, index(conv%describe(), "FDS-SDF") > 0, &
                 "the description must name the commutator: "//conv%describe())
   end subroutine test_derived_bound

   subroutine test_stated_bound(error)
      !! A stated `gradient_tolerance` replaces the derivation
      !!
      !! This is how a density consumer asks for what it needs. Deriving 1e-8
      !! through `sqrt` would need a tolerance of 1e-16, below what a molecular
      !! energy resolves, so the two demands genuinely take two numbers.
      type(error_type), allocatable, intent(out) :: error

      type(scf_convergence_t) :: conv

      conv%tolerance = 1.0e-6_dp
      conv%gradient_tolerance = 1.0e-8_dp
      call check(error, conv%commutator_bound() == 1.0e-8_dp, "a stated bound must win")
      if (allocated(error)) return
      call check(error,.not. conv%is_converged(1.0e-9_dp, 1.0e-9_dp, 1.0e-5_dp), &
                 "1e-5 must not pass a stated bound of 1e-8")
      if (allocated(error)) return
      call check(error, conv%is_converged(1.0e-9_dp, 1.0e-9_dp, 1.0e-9_dp), &
                 "1e-9 must pass a stated bound of 1e-8")
   end subroutine test_stated_bound

   subroutine test_single(error)
      !! Naming a metric bounds that quantity and no other
      !!
      !! Each case hands in two numbers that are wildly out of bounds and one
      !! that is met, and requires the answer to follow the one that was named.
      type(error_type), allocatable, intent(out) :: error

      type(scf_convergence_t) :: conv

      conv%tolerance = 1.0e-6_dp

      conv%metric = CONV_METRIC_ENERGY
      call check(error, conv%is_converged(1.0e-7_dp, 1.0_dp, 1.0_dp), &
                 "energy metric must read only dE")
      if (allocated(error)) return
      call check(error,.not. conv%is_converged(1.0e-5_dp, 0.0_dp, 0.0_dp), &
                 "energy metric must fail on dE")
      if (allocated(error)) return

      conv%metric = CONV_METRIC_COMMUTATOR
      call check(error, conv%is_converged(1.0_dp, 1.0_dp, 1.0e-7_dp), &
                 "commutator metric must read only gnorm")
      if (allocated(error)) return
      call check(error, conv%commutator_bound() == 1.0e-6_dp, &
                 "under the commutator metric the tolerance IS the bound, unsquared")
      if (allocated(error)) return

      conv%metric = CONV_METRIC_DENSITY
      call check(error, conv%is_converged(1.0_dp, 1.0e-7_dp, 1.0_dp), &
                 "density metric must read only drms")
   end subroutine test_single

   subroutine test_stall(error)
      !! **The case that makes the metric choice matter.**
      !!
      !! An interpolating accelerator can hold the energy still while the
      !! commutator stays large: EDIIS on water/6-31G sits at `dE` 1e-13 with
      !! the commutator at 1.1e-2, and stopping there lands 5.8e-5 hartree from
      !! the answer. Those are the measured numbers, used here as data.
      !!
      !! Energy alone accepts that iterate. The default does not. Anyone
      !! selecting the energy metric is choosing this, which is why it is a
      !! choice and not a default.
      type(error_type), allocatable, intent(out) :: error

      type(scf_convergence_t) :: conv
      real(dp), parameter :: STALLED_DE = 5.684e-13_dp
      real(dp), parameter :: STALLED_GNORM = 1.144e-2_dp

      conv%tolerance = 1.0e-9_dp

      conv%metric = CONV_METRIC_ENERGY
      call check(error, conv%is_converged(STALLED_DE, 1.0e-11_dp, STALLED_GNORM), &
                 "the energy metric accepts a stalled iterate -- that is what it means")
      if (allocated(error)) return

      conv%metric = CONV_METRIC_STANDARD
      call check(error,.not. conv%is_converged(STALLED_DE, 1.0e-11_dp, STALLED_GNORM), &
                 "the default must REFUSE a stalled iterate")
      if (allocated(error)) return

      conv%metric = CONV_METRIC_COMMUTATOR
      call check(error,.not. conv%is_converged(STALLED_DE, 1.0e-11_dp, STALLED_GNORM), &
                 "the commutator metric must refuse a stalled iterate")
   end subroutine test_stall

   subroutine test_parse(error)
      !! The spellings a deck may use, including the two aliases
      type(error_type), allocatable, intent(out) :: error

      integer :: metric
      logical :: ok

      call parse_convergence_metric("energy", metric, ok)
      call check(error, ok .and. metric == CONV_METRIC_ENERGY, "energy")
      if (allocated(error)) return

      ! All three name FDS-SDF: it is what DIIS extrapolates against and it is
      ! the orbital gradient. Refusing two of the three would be pedantry.
      call parse_convergence_metric("commutator", metric, ok)
      call check(error, ok .and. metric == CONV_METRIC_COMMUTATOR, "commutator")
      if (allocated(error)) return
      call parse_convergence_metric("diis", metric, ok)
      call check(error, ok .and. metric == CONV_METRIC_COMMUTATOR, "diis alias")
      if (allocated(error)) return
      call parse_convergence_metric("gradient", metric, ok)
      call check(error, ok .and. metric == CONV_METRIC_COMMUTATOR, "gradient alias")
      if (allocated(error)) return

      call parse_convergence_metric("DENSITY", metric, ok)
      call check(error, ok .and. metric == CONV_METRIC_DENSITY, "case is not significant")
      if (allocated(error)) return

      call parse_convergence_metric("wibble", metric, ok)
      call check(error,.not. ok, "an unknown spelling must be refused")
      if (allocated(error)) return
      call check(error, metric == CONV_METRIC_UNKNOWN, "and report UNKNOWN")
   end subroutine test_parse

   subroutine test_names(error)
      !! Every kind has a canonical spelling that parses back to itself
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: kinds(4) = [CONV_METRIC_STANDARD, CONV_METRIC_ENERGY, &
                                        CONV_METRIC_COMMUTATOR, CONV_METRIC_DENSITY]
      integer :: i, back
      logical :: ok

      do i = 1, size(kinds)
         call parse_convergence_metric(convergence_metric_name(kinds(i)), back, ok)
         call check(error, ok, "the canonical name must parse: "// &
                    convergence_metric_name(kinds(i)))
         if (allocated(error)) return
         call check(error, back == kinds(i), "and round trip: "// &
                    convergence_metric_name(kinds(i)))
         if (allocated(error)) return
      end do
   end subroutine test_names

end module test_mqc_scf_convergence

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_scf_convergence, only: collect_scf_convergence
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_scf_convergence", collect_scf_convergence)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
