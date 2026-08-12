module test_mqc_pcm_radii
   !! Pins the cavity radii, because nothing downstream can check them.
   !!
   !! A continuum model's solvation energy scales roughly as 1/R, so a radius that
   !! is wrong by a few percent gives a solvation energy wrong by a few percent and
   !! there is no identity anywhere that objects -- unlike a density, whose
   !! electron count gives it away. The table is therefore data we are responsible
   !! for, and this test is what makes a mistyped entry fail rather than converge.
   !!
   !! What is checked: the handful of values a chemist recognises on sight, the
   !! ordering relations that a transposed pair would break, the unit conversion,
   !! and that an element outside the table is refused rather than guessed.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_pcm_radii, only: vdw_radius, cavity_radius, DEFAULT_RADII_SCALE, &
                            MAX_PCM_ELEMENT
   implicit none
   private
   public :: collect_mqc_pcm_radii_tests

   real(dp), parameter :: BOHR_PER_ANGSTROM = 1.8897261254578281_dp
   real(dp), parameter :: TOL = 1.0e-8_dp

contains

   subroutine collect_mqc_pcm_radii_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("known_radii_are_the_published_values", test_known), &
                  new_unittest("radii_are_returned_in_bohr", test_units), &
                  new_unittest("periodic_trends_hold", test_trends), &
                  new_unittest("scaling_is_applied", test_scaling), &
                  new_unittest("every_element_in_range_has_one", test_complete), &
                  new_unittest("an_element_outside_the_table_is_refused", test_refused), &
                  new_unittest("a_nonpositive_scale_is_refused", test_bad_scale) &
                  ]
   end subroutine collect_mqc_pcm_radii_tests

   subroutine test_known(error)
      !! The values themselves, in Angstrom, against Bondi and Mantina
      !!
      !! Spelled out rather than looped: this is the test that catches a typo, and
      !! it can only do that if the expected number is written here independently
      !! of the table it checks.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: r

      call vdw_radius(1, r, err)
      call check(error, abs(r/BOHR_PER_ANGSTROM - 1.20_dp) < TOL, "H is 1.20 A")
      if (allocated(error)) return
      call vdw_radius(6, r, err)
      call check(error, abs(r/BOHR_PER_ANGSTROM - 1.70_dp) < TOL, "C is 1.70 A")
      if (allocated(error)) return
      call vdw_radius(7, r, err)
      call check(error, abs(r/BOHR_PER_ANGSTROM - 1.55_dp) < TOL, "N is 1.55 A")
      if (allocated(error)) return
      call vdw_radius(8, r, err)
      call check(error, abs(r/BOHR_PER_ANGSTROM - 1.52_dp) < TOL, "O is 1.52 A")
      if (allocated(error)) return
      call vdw_radius(9, r, err)
      call check(error, abs(r/BOHR_PER_ANGSTROM - 1.47_dp) < TOL, "F is 1.47 A")
      if (allocated(error)) return
      call vdw_radius(16, r, err)
      call check(error, abs(r/BOHR_PER_ANGSTROM - 1.80_dp) < TOL, "S is 1.80 A")
      if (allocated(error)) return
      call vdw_radius(17, r, err)
      call check(error, abs(r/BOHR_PER_ANGSTROM - 1.75_dp) < TOL, "Cl is 1.75 A")
   end subroutine test_known

   subroutine test_units(error)
      !! Bohr out, Angstrom in the table
      !!
      !! Carbon's 1.70 Angstrom is 3.21 Bohr. A table returned unconverted would
      !! give every cavity roughly half its radius, which is a large error in the
      !! solvation energy and a small-looking one in the source.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: r

      call vdw_radius(6, r, err)
      call check(error,.not. err%has_error(), "carbon must have a radius")
      if (allocated(error)) return
      call check(error, r > 3.2_dp .and. r < 3.22_dp, &
                 "carbon's radius must be about 3.21 Bohr, not 1.70 -- an "// &
                 "unconverted table halves every cavity")
   end subroutine test_units

   subroutine test_trends(error)
      !! The orderings a transposed pair of entries would break
      !!
      !! Cheap, and independent of the exact values: within a period the radius
      !! falls from left to right as the nuclear charge pulls the shell in, and
      !! down a group it grows as a shell is added. Sodium against magnesium and
      !! oxygen against sulfur are the two clearest cases.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: r_c, r_n, r_o, r_f, r_s, r_na, r_mg

      call vdw_radius(6, r_c, err)
      call vdw_radius(7, r_n, err)
      call vdw_radius(8, r_o, err)
      call vdw_radius(9, r_f, err)
      call vdw_radius(16, r_s, err)
      call vdw_radius(11, r_na, err)
      call vdw_radius(12, r_mg, err)
      call check(error,.not. err%has_error(), "all seven must have radii")
      if (allocated(error)) return

      call check(error, r_c > r_n .and. r_n > r_o .and. r_o > r_f, &
                 "the second period must shrink from carbon to fluorine")
      if (allocated(error)) return
      call check(error, r_s > r_o, "sulfur must be larger than oxygen")
      if (allocated(error)) return
      call check(error, r_na > r_mg, "sodium must be larger than magnesium")
   end subroutine test_trends

   subroutine test_scaling(error)
      !! A cavity radius is the van der Waals one times the scale
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: bare, scaled

      call vdw_radius(8, bare, err)
      call cavity_radius(8, DEFAULT_RADII_SCALE, scaled, err)
      call check(error,.not. err%has_error(), "oxygen must scale")
      if (allocated(error)) return
      call check(error, abs(scaled - DEFAULT_RADII_SCALE*bare) < TOL, &
                 "the cavity radius must be the scaled van der Waals radius")
      if (allocated(error)) return
      call check(error, abs(DEFAULT_RADII_SCALE - 1.2_dp) < TOL, &
                 "the default scale is 1.2; changing it changes every solvation "// &
                 "energy this program has ever printed, so it should not move quietly")
   end subroutine test_scaling

   subroutine test_complete(error)
      !! No gaps: every element the table claims to cover has a positive radius
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: r
      integer :: z

      do z = 1, MAX_PCM_ELEMENT
         call err%clear()
         call vdw_radius(z, r, err)
         call check(error,.not. err%has_error(), "every element in range must resolve")
         if (allocated(error)) return
         call check(error, r > 0.0_dp, "every radius must be positive")
         if (allocated(error)) return
      end do
   end subroutine test_complete

   subroutine test_refused(error)
      !! Outside the table is an error, not an extrapolation
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: r

      call vdw_radius(MAX_PCM_ELEMENT + 1, r, err)
      call check(error, err%has_error(), &
                 "an element past the table must be refused rather than guessed")
      if (allocated(error)) return
      call err%clear()
      call vdw_radius(0, r, err)
      call check(error, err%has_error(), "atomic number zero must be refused")
   end subroutine test_refused

   subroutine test_bad_scale(error)
      !! A zero or negative scale would collapse or invert the cavity
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: r

      call cavity_radius(8, 0.0_dp, r, err)
      call check(error, err%has_error(), "a zero scale must be refused")
      if (allocated(error)) return
      call err%clear()
      call cavity_radius(8, -1.2_dp, r, err)
      call check(error, err%has_error(), "a negative scale must be refused")
   end subroutine test_bad_scale

end module test_mqc_pcm_radii

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_pcm_radii, only: collect_mqc_pcm_radii_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_pcm_radii", collect_mqc_pcm_radii_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
