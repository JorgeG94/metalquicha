!! The assembled EFP2 interaction energy
module test_mqc_efp_energy
   !! Every term of this was validated on its own against GAMESS, on this exact
   !! dimer. What is left for the assembly to get wrong is not the physics but the
   !! bookkeeping: a term dropped, a term counted twice, a pair loop that visits
   !! `(a,b)` and `(b,a)`, dispersion totalled as the undamped `E6` instead of the
   !! damped `E6+E7+E8`. So this checks the breakdown term by term against the
   !! numbers already in `EFP_PLAN.md`, and then that the total is their sum.
   !!
   !! Checking the total alone would be much weaker: the terms differ in sign, and
   !! two compensating mistakes in a sum of five numbers is not a remote possibility.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_efp_potential, only: efp_potential_t, make_efp_potential, &
                                write_efp_potential
   use mqc_efp_read, only: efp_fragment_t, read_efp_potential
   use mqc_efp_energy, only: efp_energy_t, efp_interaction_energy
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_efp_energy_tests

   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp

contains

   subroutine collect_mqc_efp_energy_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("efp_energy_breakdown", test_breakdown), &
                  new_unittest("efp_energy_single_fragment", test_single) &
                  ]
   end subroutine collect_mqc_efp_energy_tests

   subroutine water_fragment(frag, err)
      type(efp_fragment_t), intent(out) :: frag
      type(error_t), intent(inout) :: err

      type(efp_potential_t) :: pot
      real(dp) :: c(3, 3)
      integer :: z(3)
      character(len=2) :: symbols(3)
      character(len=*), parameter :: path = "test_efp_energy.efp"
      integer :: unit, stat

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.00000000000000_dp, 0.00000000009155_dp, 0.10077199490609_dp, &
                   0.00000000000000_dp, 0.77250895271063_dp, -0.46780199741728_dp, &
                   0.00000000000000_dp, -0.77250895280218_dp, -0.46780199748881_dp], &
                  [3, 3])*ANG
      call make_efp_potential(z, symbols, c, "6-31g*", "WATER", pot, err)
      if (err%has_error()) return
      call write_efp_potential(pot, path, err)
      if (err%has_error()) return
      call read_efp_potential(path, frag, err)
      call pot%destroy()
      open (newunit=unit, file=path, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete")
   end subroutine water_fragment

   subroutine test_breakdown(error)
      !! All five terms of the water dimer, against the numbers each was pinned on
      type(error_type), allocatable, intent(out) :: error

      type(efp_fragment_t) :: frags(2)
      type(error_t) :: err
      type(efp_energy_t) :: e
      real(dp) :: shifts(3, 2)

      call water_fragment(frags(1), err)
      call water_fragment(frags(2), err)
      call check(error,.not. err%has_error(), "building the fragments failed")
      if (allocated(error)) return

      shifts = 0.0_dp
      shifts(1, 2) = 3.0_dp*ANG
      e = efp_interaction_energy(frags, shifts, err)
      call check(error,.not. err%has_error(), "the interaction energy failed")
      if (allocated(error)) return

      call check(error, e%electrostatics, 0.004959639_dp, thr=1.0e-8_dp, &
                 message="electrostatics")
      if (allocated(error)) return
      call check(error, e%polarization, -0.000218123_dp, thr=1.0e-8_dp, &
                 message="polarization")
      if (allocated(error)) return
      call check(error, e%exchange_repulsion, -0.001172851_dp, thr=1.0e-8_dp, &
                 message="exchange repulsion")
      if (allocated(error)) return
      ! Damped, and the three dispersion orders separately -- E7 is positive where
      ! the other two are negative, so a sign slip in the sum would not hide.
      call check(error, e%dispersion_e6, -0.0005163981_dp, thr=1.0e-9_dp, &
                 message="dispersion E6")
      if (allocated(error)) return
      call check(error, e%dispersion_e7, 0.0000598618_dp, thr=1.0e-9_dp, &
                 message="dispersion E7")
      if (allocated(error)) return
      call check(error, e%dispersion_e8, -0.0001143337_dp, thr=1.0e-9_dp, &
                 message="dispersion E8")
      if (allocated(error)) return
      call check(error, e%dispersion, -0.000570870_dp, thr=1.0e-8_dp, &
                 message="dispersion total")
      if (allocated(error)) return
      call check(error, e%charge_transfer, 0.000081783_dp, thr=1.0e-8_dp, &
                 message="charge transfer")
      if (allocated(error)) return

      ! And that the total is the sum of the parts rather than anything else.
      call check(error, e%total, e%electrostatics + e%polarization &
                 + e%exchange_repulsion + e%dispersion + e%charge_transfer, &
                 thr=1.0e-14_dp, message="the total is not the sum of its terms")
   end subroutine test_breakdown

   subroutine test_single(error)
      !! One fragment interacts with nothing
      !!
      !! Worth asserting because every term is a loop over pairs and an empty loop is
      !! easy to write as a special case that is never reached. A single fragment is
      !! also what a caller building up a cluster starts from.
      type(error_type), allocatable, intent(out) :: error

      type(efp_fragment_t) :: frags(1)
      type(error_t) :: err
      type(efp_energy_t) :: e
      real(dp) :: shifts(3, 1)

      call water_fragment(frags(1), err)
      call check(error,.not. err%has_error(), "building the fragment failed")
      if (allocated(error)) return

      shifts = 0.0_dp
      e = efp_interaction_energy(frags, shifts, err)
      call check(error,.not. err%has_error(), "the interaction energy failed")
      if (allocated(error)) return
      call check(error, e%total, 0.0_dp, thr=0.0_dp, &
                 message="a lone fragment should have no interaction energy")
   end subroutine test_single

end module test_mqc_efp_energy

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_efp_energy, only: collect_mqc_efp_energy_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_efp_energy", collect_mqc_efp_energy_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
