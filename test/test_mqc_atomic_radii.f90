!! Unit tests for the consolidated atomic radii tables
module test_mqc_atomic_radii
   !! These tables were gathered into one module from four places. The point of
   !! these tests is that the gathering moved no digits and merged no two tables
   !! that disagree, so the spot values below are quoted from the original
   !! sources rather than read back out of the module.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_atomic_radii, only: covalent_radius_cordero, covalent_radius_emsley, &
                               vdw_radius_bondi, vdw_radius_geodesic, &
                               MAX_Z_CORDERO, MAX_Z_EMSLEY, MAX_Z_BONDI, &
                               MAX_Z_GEODESIC, GEODESIC_RADIUS_DEFAULT
   implicit none
   private

   public :: collect_mqc_atomic_radii

   real(dp), parameter :: TOL = 1.0e-12_dp

contains

   subroutine collect_mqc_atomic_radii(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("radii_cordero_matches_the_published_values", test_cordero), &
                  new_unittest("radii_emsley_matches_the_gamess_values", test_emsley), &
                  new_unittest("radii_bondi_matches_the_published_values", test_bondi), &
                  new_unittest("radii_geodesic_matches_the_gamess_values", test_geodesic), &
                  new_unittest("radii_the_tables_stay_distinct", test_distinct), &
                  new_unittest("radii_refuse_out_of_range_elements", test_out_of_range), &
                  new_unittest("radii_are_positive_across_each_table", test_positive) &
                  ]
   end subroutine collect_mqc_atomic_radii

   subroutine test_cordero(error)
      !! Cordero et al., Dalton Trans. 2008, 2832 -- Angstrom
      type(error_type), allocatable, intent(out) :: error

      call check(error, covalent_radius_cordero(1), 0.31_dp, thr=TOL)
      if (allocated(error)) return
      call check(error, covalent_radius_cordero(6), 0.76_dp, thr=TOL)
      if (allocated(error)) return
      call check(error, covalent_radius_cordero(8), 0.66_dp, thr=TOL)
      if (allocated(error)) return
      ! Curium, the last element Cordero tabulates.
      call check(error, covalent_radius_cordero(MAX_Z_CORDERO), 1.69_dp, thr=TOL)
   end subroutine test_cordero

   subroutine test_emsley(error)
      !! Emsley's table as GAMESS uses it -- deliberately not Cordero's
      type(error_type), allocatable, intent(out) :: error

      call check(error, covalent_radius_emsley(1), 0.32_dp, thr=TOL)
      if (allocated(error)) return
      call check(error, covalent_radius_emsley(6), 0.77_dp, thr=TOL)
      if (allocated(error)) return
      call check(error, covalent_radius_emsley(8), 0.73_dp, thr=TOL)
      if (allocated(error)) return
      ! Krypton, where the table stops.
      call check(error, covalent_radius_emsley(MAX_Z_EMSLEY), 1.12_dp, thr=TOL)
   end subroutine test_emsley

   subroutine test_bondi(error)
      !! Bondi (1964), with Mantina (2009) for six light metals
      type(error_type), allocatable, intent(out) :: error

      call check(error, vdw_radius_bondi(1), 1.20_dp, thr=TOL)
      if (allocated(error)) return
      call check(error, vdw_radius_bondi(6), 1.70_dp, thr=TOL)
      if (allocated(error)) return
      call check(error, vdw_radius_bondi(8), 1.52_dp, thr=TOL)
      if (allocated(error)) return
      ! Sodium is one of the Mantina values.
      call check(error, vdw_radius_bondi(11), 2.27_dp, thr=TOL)
      if (allocated(error)) return
      call check(error, vdw_radius_bondi(MAX_Z_BONDI), 1.88_dp, thr=TOL)
   end subroutine test_bondi

   subroutine test_geodesic(error)
      !! The screening grid's radii, including the substituted default
      type(error_type), allocatable, intent(out) :: error

      call check(error, vdw_radius_geodesic(1), 1.20_dp, thr=TOL)
      if (allocated(error)) return
      call check(error, vdw_radius_geodesic(6), 1.50_dp, thr=TOL)
      if (allocated(error)) return
      call check(error, vdw_radius_geodesic(MAX_Z_GEODESIC), 1.80_dp, thr=TOL)
      if (allocated(error)) return
      ! Helium sits inside the table's range but has no entry, so it takes the
      ! default -- this one substitutes rather than refuses, because the grid is
      ! a fitting aid and a slightly wrong sphere moves points, not an energy.
      call check(error, vdw_radius_geodesic(2), GEODESIC_RADIUS_DEFAULT, thr=TOL)
      if (allocated(error)) return
      call check(error, vdw_radius_geodesic(50), GEODESIC_RADIUS_DEFAULT, thr=TOL)
   end subroutine test_geodesic

   subroutine test_distinct(error)
      !! The reason these are four accessors and not one
      !!
      !! If a later cleanup ever collapses them into a single table, this fails.
      !! Carbon is the clearest case: four parametrisations, four numbers.
      type(error_type), allocatable, intent(out) :: error

      call check(error, abs(covalent_radius_cordero(6) - covalent_radius_emsley(6)) > 1.0e-6_dp, &
                 "Cordero and Emsley must stay separate tables")
      if (allocated(error)) return
      call check(error, abs(vdw_radius_bondi(6) - vdw_radius_geodesic(6)) > 1.0e-6_dp, &
                 "Bondi and geodesic must stay separate tables")
      if (allocated(error)) return
      ! A covalent radius is much smaller than a van der Waals one; if these
      ! ever cross, two tables have been swapped somewhere.
      call check(error, covalent_radius_cordero(6) < vdw_radius_bondi(6), &
                 "a covalent radius must be smaller than a van der Waals radius")
   end subroutine test_distinct

   subroutine test_out_of_range(error)
      !! Zero is a refusal, so a caller cannot invent a bond from a default
      type(error_type), allocatable, intent(out) :: error

      call check(error, covalent_radius_cordero(0), 0.0_dp, thr=TOL)
      if (allocated(error)) return
      call check(error, covalent_radius_cordero(MAX_Z_CORDERO + 1), 0.0_dp, thr=TOL)
      if (allocated(error)) return
      call check(error, covalent_radius_emsley(MAX_Z_EMSLEY + 1), 0.0_dp, thr=TOL)
      if (allocated(error)) return
      call check(error, vdw_radius_bondi(MAX_Z_BONDI + 1), 0.0_dp, thr=TOL)
      if (allocated(error)) return
      call check(error, vdw_radius_bondi(-1), 0.0_dp, thr=TOL)
   end subroutine test_out_of_range

   subroutine test_positive(error)
      !! Every tabulated entry is a real length
      type(error_type), allocatable, intent(out) :: error
      integer :: z

      do z = 1, MAX_Z_CORDERO
         call check(error, covalent_radius_cordero(z) > 0.0_dp, "Cordero entry must be positive")
         if (allocated(error)) return
      end do
      do z = 1, MAX_Z_EMSLEY
         call check(error, covalent_radius_emsley(z) > 0.0_dp, "Emsley entry must be positive")
         if (allocated(error)) return
      end do
      do z = 1, MAX_Z_BONDI
         call check(error, vdw_radius_bondi(z) > 0.0_dp, "Bondi entry must be positive")
         if (allocated(error)) return
      end do
      ! The geodesic table has deliberate zero entries, but the accessor must
      ! never hand one back.
      do z = 1, MAX_Z_GEODESIC
         call check(error, vdw_radius_geodesic(z) > 0.0_dp, "geodesic lookup must be positive")
         if (allocated(error)) return
      end do
   end subroutine test_positive

end module test_mqc_atomic_radii

program tester_mqc_atomic_radii
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_atomic_radii, only: collect_mqc_atomic_radii
   implicit none

   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_atomic_radii", collect_mqc_atomic_radii)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_atomic_radii
