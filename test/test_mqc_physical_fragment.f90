module test_mqc_physical_fragment
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_physical_fragment, only: to_angstrom, to_bohr, initialize_system_geometry, &
                                    build_fragment_from_indices, fragment_centroid, &
                                    fragment_center_of_mass, distance_between_points, &
                                    distance_between_fragments, minimal_distance_between_fragments, &
                                    system_geometry_t, physical_fragment_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_physical_fragment_tests

   ! Tolerance for floating point comparisons
   real(dp), parameter :: tol = 1.0e-8_dp

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_physical_fragment_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("unit_conversions_angstrom_to_bohr", test_angstrom_to_bohr), &
                  new_unittest("unit_conversions_bohr_to_angstrom", test_bohr_to_angstrom), &
                  new_unittest("unit_conversions_round_trip", test_unit_conversions_round_trip), &
                  new_unittest("initialize_system_geometry_valid", test_initialize_system_geometry), &
                  new_unittest("initialize_system_geometry_invalid", test_initialize_system_invalid), &
                  new_unittest("build_fragment_single_monomer", test_build_fragment_single), &
                  new_unittest("build_fragment_multiple_monomers", test_build_fragment_multiple), &
                  new_unittest("fragment_centroid", test_fragment_centroid), &
                  new_unittest("fragment_center_of_mass", test_fragment_center_of_mass), &
                  new_unittest("distance_between_points", test_distance_points), &
                  new_unittest("distance_between_fragments_centroid", test_distance_fragments_centroid), &
                  new_unittest("distance_between_fragments_com", test_distance_fragments_com), &
                  new_unittest("minimal_distance_between_fragments", test_minimal_distance_fragments), &
                  new_unittest("fragment_destroy", test_fragment_destroy_proc), &
                  new_unittest("system_destroy", test_system_destroy_proc) &
                  ]
   end subroutine collect_mqc_physical_fragment_tests

   subroutine test_angstrom_to_bohr(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: bohr_val

      ! 1 Angstrom should be approximately 1.889726 Bohr
      bohr_val = to_bohr(1.0_dp)
      call check(error, abs(bohr_val - 1.889726_dp) < 1.0e-5_dp, &
                 "1 Angstrom should convert to approximately 1.889726 Bohr")
   end subroutine test_angstrom_to_bohr

   subroutine test_bohr_to_angstrom(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: angstrom_val

      ! 1 Bohr should be approximately 0.529177 Angstrom
      angstrom_val = to_angstrom(1.0_dp)
      call check(error, abs(angstrom_val - 0.529177_dp) < 1.0e-5_dp, &
                 "1 Bohr should convert to approximately 0.529177 Angstrom")
   end subroutine test_bohr_to_angstrom

   subroutine test_unit_conversions_round_trip(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: original, converted

      original = 5.0_dp
      converted = to_angstrom(to_bohr(original))
      call check(error, abs(converted - original) < tol, &
                 "Round trip conversion should preserve value")
   end subroutine test_unit_conversions_round_trip

   subroutine test_initialize_system_geometry(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      integer :: stat
      character(len=:), allocatable :: errmsg

      ! Create test files
      call create_test_water_trimer()

      ! Initialize system geometry
      call initialize_system_geometry("test_water_trimer.xyz", "test_water_monomer.xyz", &
                                      sys_geom, stat, errmsg)

      if (stat /= 0) then
         call check(error, .false., "Failed to initialize system: "//errmsg)
         call cleanup_test_files()
         return
      end if

      ! Check that system was initialized correctly
      call check(error, sys_geom%n_monomers, 3, "Should have 3 monomers")
      if (allocated(error)) then
         call sys_geom%destroy()
         call cleanup_test_files()
         return
      end if

      call check(error, sys_geom%atoms_per_monomer, 3, "Should have 3 atoms per monomer")
      if (allocated(error)) then
         call sys_geom%destroy()
         call cleanup_test_files()
         return
      end if

      call check(error, sys_geom%total_atoms, 9, "Should have 9 total atoms")
      if (allocated(error)) then
         call sys_geom%destroy()
         call cleanup_test_files()
         return
      end if

      ! Check that coordinates are allocated
      call check(error, allocated(sys_geom%coordinates), &
                 "Coordinates should be allocated")
      if (allocated(error)) then
         call sys_geom%destroy()
         call cleanup_test_files()
         return
      end if

      call check(error, allocated(sys_geom%element_numbers), &
                 "Element numbers should be allocated")

      call sys_geom%destroy()
      call cleanup_test_files()
   end subroutine test_initialize_system_geometry

   subroutine test_initialize_system_invalid(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      integer :: stat
      character(len=:), allocatable :: errmsg

      ! Create test files with mismatched sizes
      call create_test_water_trimer()
      call create_mismatched_monomer()

      ! Initialize system geometry - should fail
      call initialize_system_geometry("test_water_trimer.xyz", "test_mismatched.xyz", &
                                      sys_geom, stat, errmsg)

      call check(error, stat /= 0, "Should fail with mismatched monomer size")

      call cleanup_test_files()
      call delete_file("test_mismatched.xyz")
   end subroutine test_initialize_system_invalid

   subroutine test_build_fragment_single(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      type(physical_fragment_t) :: fragment
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_water_trimer()

      call initialize_system_geometry("test_water_trimer.xyz", "test_water_monomer.xyz", &
                                      sys_geom, stat, errmsg)

      if (stat /= 0) then
         call check(error, .false., "Failed to initialize system")
         call cleanup_test_files()
         return
      end if

      ! Build fragment from first monomer only
      call build_fragment_from_indices(sys_geom, [1], fragment)

      call check(error, fragment%n_atoms, 3, "Fragment should have 3 atoms")
      if (allocated(error)) then
         call fragment%destroy()
         call sys_geom%destroy()
         call cleanup_test_files()
         return
      end if

      ! Check first atom is oxygen (element 8)
      call check(error, fragment%element_numbers(1), 8, &
                 "First atom should be oxygen")

      call fragment%destroy()
      call sys_geom%destroy()
      call cleanup_test_files()
   end subroutine test_build_fragment_single

   subroutine test_build_fragment_multiple(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      type(physical_fragment_t) :: fragment
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_water_trimer()

      call initialize_system_geometry("test_water_trimer.xyz", "test_water_monomer.xyz", &
                                      sys_geom, stat, errmsg)

      if (stat /= 0) then
         call check(error, .false., "Failed to initialize system")
         call cleanup_test_files()
         return
      end if

      ! Build fragment from monomers 1 and 3
      call build_fragment_from_indices(sys_geom, [1, 3], fragment)

      call check(error, fragment%n_atoms, 6, "Fragment should have 6 atoms (2 waters)")
      if (allocated(error)) then
         call fragment%destroy()
         call sys_geom%destroy()
         call cleanup_test_files()
         return
      end if

      ! Check that we have 2 oxygen atoms
      call check(error, fragment%element_numbers(1), 8, &
                 "First atom should be oxygen")
      if (allocated(error)) then
         call fragment%destroy()
         call sys_geom%destroy()
         call cleanup_test_files()
         return
      end if

      call check(error, fragment%element_numbers(4), 8, &
                 "Fourth atom should be oxygen (from second water)")

      call fragment%destroy()
      call sys_geom%destroy()
      call cleanup_test_files()
   end subroutine test_build_fragment_multiple

   subroutine test_fragment_centroid(error)
      type(error_type), allocatable, intent(out) :: error
      type(physical_fragment_t) :: fragment
      real(dp) :: centroid(3)

      ! Create a simple fragment with known coordinates
      fragment%n_atoms = 3
      allocate (fragment%element_numbers(3))
      allocate (fragment%coordinates(3, 3))

      fragment%element_numbers = [1, 1, 1]  ! All hydrogen for simplicity
      ! Create a triangle at (0,0,0), (3,0,0), (0,3,0)
      fragment%coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      fragment%coordinates(:, 2) = [3.0_dp, 0.0_dp, 0.0_dp]
      fragment%coordinates(:, 3) = [0.0_dp, 3.0_dp, 0.0_dp]

      centroid = fragment_centroid(fragment)

      ! Centroid should be at (1, 1, 0)
      call check(error, abs(centroid(1) - 1.0_dp) < tol, &
                 "Centroid x-coordinate should be 1.0")
      if (allocated(error)) then
         call fragment%destroy()
         return
      end if

      call check(error, abs(centroid(2) - 1.0_dp) < tol, &
                 "Centroid y-coordinate should be 1.0")
      if (allocated(error)) then
         call fragment%destroy()
         return
      end if

      call check(error, abs(centroid(3)) < tol, &
                 "Centroid z-coordinate should be 0.0")

      call fragment%destroy()
   end subroutine test_fragment_centroid

   subroutine test_fragment_center_of_mass(error)
      type(error_type), allocatable, intent(out) :: error
      type(physical_fragment_t) :: fragment
      real(dp) :: com(3)

      ! Create a water molecule with known coordinates
      fragment%n_atoms = 3
      allocate (fragment%element_numbers(3))
      allocate (fragment%coordinates(3, 3))

      fragment%element_numbers = [8, 1, 1]  ! O, H, H
      ! Simple linear arrangement for easy testing
      fragment%coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]  ! O at origin
      fragment%coordinates(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]  ! H at x=1
      fragment%coordinates(:, 3) = [-1.0_dp, 0.0_dp, 0.0_dp]  ! H at x=-1

      com = fragment_center_of_mass(fragment)

      ! COM should be very close to origin since O is much heavier than H
      ! and positioned at origin, with H atoms symmetrically placed
      call check(error, abs(com(1)) < 0.2_dp, &
                 "COM x-coordinate should be near 0 (heavy O at origin)")
      if (allocated(error)) then
         call fragment%destroy()
         return
      end if

      call check(error, abs(com(2)) < tol .and. abs(com(3)) < tol, &
                 "COM y and z coordinates should be 0")

      call fragment%destroy()
   end subroutine test_fragment_center_of_mass

   subroutine test_distance_points(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: point1(3), point2(3), dist

      point1 = [0.0_dp, 0.0_dp, 0.0_dp]
      point2 = [3.0_dp, 4.0_dp, 0.0_dp]

      dist = distance_between_points(point1, point2)

      ! Distance should be 5.0 (3-4-5 triangle)
      call check(error, abs(dist - 5.0_dp) < tol, &
                 "Distance should be 5.0 for 3-4-5 triangle")
   end subroutine test_distance_points

   subroutine test_distance_fragments_centroid(error)
      type(error_type), allocatable, intent(out) :: error
      type(physical_fragment_t) :: frag1, frag2
      real(dp) :: dist

      ! Create two simple fragments
      frag1%n_atoms = 2
      allocate (frag1%element_numbers(2))
      allocate (frag1%coordinates(3, 2))
      frag1%element_numbers = [1, 1]
      frag1%coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      frag1%coordinates(:, 2) = [2.0_dp, 0.0_dp, 0.0_dp]
      ! Centroid of frag1 is at (1, 0, 0)

      frag2%n_atoms = 2
      allocate (frag2%element_numbers(2))
      allocate (frag2%coordinates(3, 2))
      frag2%element_numbers = [1, 1]
      frag2%coordinates(:, 1) = [4.0_dp, 0.0_dp, 0.0_dp]
      frag2%coordinates(:, 2) = [6.0_dp, 0.0_dp, 0.0_dp]
      ! Centroid of frag2 is at (5, 0, 0)

      dist = distance_between_fragments(frag1, frag2, .false.)

      ! Distance between centroids should be 4.0
      call check(error, abs(dist - 4.0_dp) < tol, &
                 "Distance between fragment centroids should be 4.0")

      call frag1%destroy()
      call frag2%destroy()
   end subroutine test_distance_fragments_centroid

   subroutine test_distance_fragments_com(error)
      type(error_type), allocatable, intent(out) :: error
      type(physical_fragment_t) :: frag1, frag2
      real(dp) :: dist

      ! Create two fragments with different masses
      frag1%n_atoms = 2
      allocate (frag1%element_numbers(2))
      allocate (frag1%coordinates(3, 2))
      frag1%element_numbers = [8, 1]  ! O and H
      frag1%coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      frag1%coordinates(:, 2) = [2.0_dp, 0.0_dp, 0.0_dp]

      frag2%n_atoms = 2
      allocate (frag2%element_numbers(2))
      allocate (frag2%coordinates(3, 2))
      frag2%element_numbers = [8, 1]  ! O and H
      frag2%coordinates(:, 1) = [10.0_dp, 0.0_dp, 0.0_dp]
      frag2%coordinates(:, 2) = [12.0_dp, 0.0_dp, 0.0_dp]

      dist = distance_between_fragments(frag1, frag2, .true.)

      ! Distance should be positive
      call check(error, dist > 0.0_dp, &
                 "Distance between fragments should be positive")

      call frag1%destroy()
      call frag2%destroy()
   end subroutine test_distance_fragments_com

   subroutine test_minimal_distance_fragments(error)
      type(error_type), allocatable, intent(out) :: error
      type(physical_fragment_t) :: frag1, frag2
      real(dp) :: min_dist

      ! Create two fragments with known positions
      ! Fragment 1: atoms at (0,0,0) and (1,0,0)
      frag1%n_atoms = 2
      allocate (frag1%element_numbers(2))
      allocate (frag1%coordinates(3, 2))
      frag1%element_numbers = [1, 1]
      frag1%coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      frag1%coordinates(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]

      ! Fragment 2: atoms at (5,0,0) and (6,0,0)
      frag2%n_atoms = 2
      allocate (frag2%element_numbers(2))
      allocate (frag2%coordinates(3, 2))
      frag2%element_numbers = [1, 1]
      frag2%coordinates(:, 1) = [5.0_dp, 0.0_dp, 0.0_dp]
      frag2%coordinates(:, 2) = [6.0_dp, 0.0_dp, 0.0_dp]

      ! Minimal distance should be between (1,0,0) and (5,0,0) = 4.0
      min_dist = minimal_distance_between_fragments(frag1, frag2)

      call check(error, abs(min_dist - 4.0_dp) < tol, &
                 "Minimal distance should be 4.0")

      call frag1%destroy()
      call frag2%destroy()

      if (allocated(error)) return

      ! Test with fragments at different positions
      ! Fragment 1: triangle with vertices at (0,0,0), (2,0,0), (1,2,0)
      frag1%n_atoms = 3
      allocate (frag1%element_numbers(3))
      allocate (frag1%coordinates(3, 3))
      frag1%element_numbers = [1, 1, 1]
      frag1%coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      frag1%coordinates(:, 2) = [2.0_dp, 0.0_dp, 0.0_dp]
      frag1%coordinates(:, 3) = [1.0_dp, 2.0_dp, 0.0_dp]

      ! Fragment 2: single atom at (1.5, 3.0, 0.0)
      frag2%n_atoms = 1
      allocate (frag2%element_numbers(1))
      allocate (frag2%coordinates(3, 1))
      frag2%element_numbers = [1]
      frag2%coordinates(:, 1) = [1.5_dp, 3.0_dp, 0.0_dp]

      ! Minimal distance should be from (1,2,0) to (1.5,3.0,0)
      ! Distance = sqrt((1.5-1)^2 + (3-2)^2 + 0^2) = sqrt(0.25 + 1) = sqrt(1.25)
      min_dist = minimal_distance_between_fragments(frag1, frag2)

      call check(error, abs(min_dist - sqrt(1.25_dp)) < tol, &
                 "Minimal distance should be sqrt(1.25)")

      call frag1%destroy()
      call frag2%destroy()
   end subroutine test_minimal_distance_fragments

   subroutine test_fragment_destroy_proc(error)
      type(error_type), allocatable, intent(out) :: error
      type(physical_fragment_t) :: fragment

      fragment%n_atoms = 3
      allocate (fragment%element_numbers(3))
      allocate (fragment%coordinates(3, 3))

      call fragment%destroy()

      call check(error,.not. allocated(fragment%element_numbers), &
                 "Element numbers should be deallocated")
      if (allocated(error)) return

      call check(error,.not. allocated(fragment%coordinates), &
                 "Coordinates should be deallocated")
      if (allocated(error)) return

      call check(error, fragment%n_atoms == 0, &
                 "n_atoms should be reset to 0")
   end subroutine test_fragment_destroy_proc

   subroutine test_system_destroy_proc(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom

      sys_geom%n_monomers = 5
      sys_geom%atoms_per_monomer = 3
      sys_geom%total_atoms = 15
      allocate (sys_geom%element_numbers(15))
      allocate (sys_geom%coordinates(3, 15))

      call sys_geom%destroy()

      call check(error,.not. allocated(sys_geom%element_numbers), &
                 "Element numbers should be deallocated")
      if (allocated(error)) return

      call check(error,.not. allocated(sys_geom%coordinates), &
                 "Coordinates should be deallocated")
      if (allocated(error)) return

      call check(error, sys_geom%n_monomers == 0, &
                 "n_monomers should be reset to 0")
      if (allocated(error)) return

      call check(error, sys_geom%atoms_per_monomer == 0, &
                 "atoms_per_monomer should be reset to 0")
      if (allocated(error)) return

      call check(error, sys_geom%total_atoms == 0, &
                 "total_atoms should be reset to 0")
   end subroutine test_system_destroy_proc

   ! Helper subroutines for creating test files

   subroutine create_test_water_monomer()
      integer :: unit

      open (newunit=unit, file="test_water_monomer.xyz", status="replace", action="write")
      write (unit, '(a)') "3"
      write (unit, '(a)') "Water monomer"
      write (unit, '(a)') "O  0.0  0.0  0.0"
      write (unit, '(a)') "H  0.757  0.586  0.0"
      write (unit, '(a)') "H -0.757  0.586  0.0"
      close (unit)
   end subroutine create_test_water_monomer

   subroutine create_test_water_trimer()
      integer :: unit

      call create_test_water_monomer()

      open (newunit=unit, file="test_water_trimer.xyz", status="replace", action="write")
      write (unit, '(a)') "9"
      write (unit, '(a)') "Water trimer"
      ! First water
      write (unit, '(a)') "O  0.0  0.0  0.0"
      write (unit, '(a)') "H  0.757  0.586  0.0"
      write (unit, '(a)') "H -0.757  0.586  0.0"
      ! Second water
      write (unit, '(a)') "O  3.0  0.0  0.0"
      write (unit, '(a)') "H  3.757  0.586  0.0"
      write (unit, '(a)') "H  2.243  0.586  0.0"
      ! Third water
      write (unit, '(a)') "O  1.5  2.5  0.0"
      write (unit, '(a)') "H  2.257  3.086  0.0"
      write (unit, '(a)') "H  0.743  3.086  0.0"
      close (unit)
   end subroutine create_test_water_trimer

   subroutine create_mismatched_monomer()
      integer :: unit

      open (newunit=unit, file="test_mismatched.xyz", status="replace", action="write")
      write (unit, '(a)') "4"
      write (unit, '(a)') "Mismatched monomer"
      write (unit, '(a)') "O  0.0  0.0  0.0"
      write (unit, '(a)') "H  0.757  0.586  0.0"
      write (unit, '(a)') "H -0.757  0.586  0.0"
      write (unit, '(a)') "H  0.0  -0.5  0.0"
      close (unit)
   end subroutine create_mismatched_monomer

   subroutine cleanup_test_files()
      call delete_file("test_water_monomer.xyz")
      call delete_file("test_water_trimer.xyz")
   end subroutine cleanup_test_files

   subroutine delete_file(filename)
      character(len=*), intent(in) :: filename
      integer :: unit, stat

      open (newunit=unit, file=filename, status="old", iostat=stat)
      if (stat == 0) then
         close (unit, status="delete")
      end if
   end subroutine delete_file

end module test_mqc_physical_fragment

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_physical_fragment, only: collect_mqc_physical_fragment_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_physical_fragment", collect_mqc_physical_fragment_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester
