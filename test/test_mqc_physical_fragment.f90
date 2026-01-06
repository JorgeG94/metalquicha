module test_mqc_physical_fragment
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_physical_fragment, only: to_angstrom, to_bohr, initialize_system_geometry, &
                                    build_fragment_from_indices, &
                                    system_geometry_t, physical_fragment_t
   use mqc_config_parser, only: bond_t
   use mqc_error, only: error_t
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
                  new_unittest("fragment_destroy", test_fragment_destroy_proc), &
                  new_unittest("system_destroy", test_system_destroy_proc), &
                  new_unittest("h_cap_single_broken_bond", test_h_cap_single_broken), &
                  new_unittest("h_cap_no_broken_bonds", test_h_cap_no_broken), &
                  new_unittest("h_cap_dimer_intact_internal", test_h_cap_dimer_intact), &
                  new_unittest("h_cap_positions_and_elements", test_h_cap_positions) &
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
      type(error_t) :: parse_error

      ! Create test files
      call create_test_water_trimer()

      ! Initialize system geometry
      call initialize_system_geometry("test_water_trimer.xyz", "test_water_monomer.xyz", &
                                      sys_geom, parse_error)

      if (parse_error%has_error()) then
         call check(error, .false., "Failed to initialize system: "//parse_error%get_message())
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
      type(error_t) :: parse_error

      ! Create test files with mismatched sizes
      call create_test_water_trimer()
      call create_mismatched_monomer()

      ! Initialize system geometry - should fail
      call initialize_system_geometry("test_water_trimer.xyz", "test_mismatched.xyz", &
                                      sys_geom, parse_error)

      call check(error, parse_error%has_error(), "Should fail with mismatched monomer size")

      call cleanup_test_files()
      call delete_file("test_mismatched.xyz")
   end subroutine test_initialize_system_invalid

   subroutine test_build_fragment_single(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      type(physical_fragment_t) :: fragment
      type(error_t) :: parse_error
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_water_trimer()

      call initialize_system_geometry("test_water_trimer.xyz", "test_water_monomer.xyz", &
                                      sys_geom, parse_error)

      if (parse_error%has_error()) then
         call check(error, .false., "Failed to initialize system")
         call cleanup_test_files()
         return
      end if

      ! Build fragment from first monomer only
      call build_fragment_from_indices(sys_geom, [1], fragment, parse_error)

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
      type(error_t) :: parse_error
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_water_trimer()

      call initialize_system_geometry("test_water_trimer.xyz", "test_water_monomer.xyz", &
                                      sys_geom, parse_error)

      if (parse_error%has_error()) then
         call check(error, .false., "Failed to initialize system")
         call cleanup_test_files()
         return
      end if

      ! Build fragment from monomers 1 and 3
      call build_fragment_from_indices(sys_geom, [1, 3], fragment, parse_error)

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

   ! ============================================================================
   ! Hydrogen capping tests
   ! ============================================================================

   subroutine test_h_cap_single_broken(error)
      !! Test single monomer with one broken bond → should add 1 H-cap
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      type(physical_fragment_t) :: fragment
      type(bond_t), allocatable :: bonds(:)
      type(error_t) :: parse_error
      integer :: stat
      character(len=:), allocatable :: errmsg

      ! Create test system: 2 monomers, 3 atoms each
      call create_test_water_trimer()
      call initialize_system_geometry("test_water_trimer.xyz", "test_water_monomer.xyz", &
                                      sys_geom, parse_error)

      if (parse_error%has_error()) then
         call check(error, .false., "Failed to initialize system: "//parse_error%get_message())
         call cleanup_test_files()
         return
      end if

      ! Create one broken bond between monomer 1 (atom 0) and monomer 2 (atom 3)
      ! Using 0-based indexing as in the bond_t structure
      allocate (bonds(1))
      bonds(1)%atom_i = 0  ! First atom of monomer 1
      bonds(1)%atom_j = 3  ! First atom of monomer 2
      bonds(1)%order = 1
      bonds(1)%is_broken = .true.

      ! Build fragment from monomer 1 only
      call build_fragment_from_indices(sys_geom, [1], fragment, parse_error, bonds)

      ! Should have 3 real atoms + 1 H-cap = 4 total
      call check(error, fragment%n_atoms, 4, &
                 "Fragment with 1 broken bond should have 4 atoms (3 real + 1 cap)")
      if (allocated(error)) return

      ! Check n_caps field
      call check(error, fragment%n_caps, 1, &
                 "Should have exactly 1 H-cap")
      if (allocated(error)) return

      ! Check that last atom is hydrogen (element 1)
      call check(error, fragment%element_numbers(4), 1, &
                 "Last atom should be hydrogen (H-cap)")
      if (allocated(error)) return

      ! Check cap_replaces_atom is set correctly (should be atom 3, 0-indexed)
      call check(error, allocated(fragment%cap_replaces_atom), &
                 "cap_replaces_atom should be allocated")
      if (allocated(error)) return

      call check(error, fragment%cap_replaces_atom(1), 3, &
                 "H-cap should replace atom 3")

      call fragment%destroy()
      call sys_geom%destroy()
      call cleanup_test_files()
      if (allocated(bonds)) deallocate (bonds)
   end subroutine test_h_cap_single_broken

   subroutine test_h_cap_no_broken(error)
      !! Test monomer with no broken bonds → should add 0 H-caps
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      type(physical_fragment_t) :: fragment
      type(bond_t), allocatable :: bonds(:)
      type(error_t) :: parse_error
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_water_trimer()
      call initialize_system_geometry("test_water_trimer.xyz", "test_water_monomer.xyz", &
                                      sys_geom, parse_error)

      if (parse_error%has_error()) then
         call check(error, .false., "Failed to initialize system: "//parse_error%get_message())
         call cleanup_test_files()
         return
      end if

      ! Create bond but mark it as NOT broken
      allocate (bonds(1))
      bonds(1)%atom_i = 0
      bonds(1)%atom_j = 3
      bonds(1)%order = 1
      bonds(1)%is_broken = .false.  ! Not broken!

      ! Build fragment from monomer 1 only
      call build_fragment_from_indices(sys_geom, [1], fragment, parse_error, bonds)

      ! Should have only 3 atoms (no caps)
      call check(error, fragment%n_atoms, 3, &
                 "Fragment with no broken bonds should have 3 atoms (no caps)")
      if (allocated(error)) return

      ! Check n_caps field
      call check(error, fragment%n_caps, 0, &
                 "Should have 0 H-caps")

      call fragment%destroy()
      call sys_geom%destroy()
      call cleanup_test_files()
      if (allocated(bonds)) deallocate (bonds)
   end subroutine test_h_cap_no_broken

   subroutine test_h_cap_dimer_intact(error)
      !! Test dimer with intact internal bond → should add 0 H-caps between monomers
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      type(physical_fragment_t) :: fragment
      type(bond_t), allocatable :: bonds(:)
      type(error_t) :: parse_error
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_water_trimer()
      call initialize_system_geometry("test_water_trimer.xyz", "test_water_monomer.xyz", &
                                      sys_geom, parse_error)

      if (parse_error%has_error()) then
         call check(error, .false., "Failed to initialize system: "//parse_error%get_message())
         call cleanup_test_files()
         return
      end if

      ! Create broken bond between atoms 0 and 3 (monomers 1-2)
      ! But we're building a fragment with BOTH monomers
      allocate (bonds(1))
      bonds(1)%atom_i = 0  ! Monomer 1
      bonds(1)%atom_j = 3  ! Monomer 2
      bonds(1)%order = 1
      bonds(1)%is_broken = .true.

      ! Build fragment from monomers 1 AND 2 (dimer)
      call build_fragment_from_indices(sys_geom, [1, 2], fragment, parse_error, bonds)

      ! Should have 6 atoms (2 waters) and NO caps (bond is internal to fragment)
      call check(error, fragment%n_atoms, 6, &
                 "Dimer with internal bond should have 6 atoms (no caps)")
      if (allocated(error)) return

      call check(error, fragment%n_caps, 0, &
                 "Dimer with intact internal bond should have 0 caps")

      call fragment%destroy()
      call sys_geom%destroy()
      call cleanup_test_files()
      if (allocated(bonds)) deallocate (bonds)
   end subroutine test_h_cap_dimer_intact

   subroutine test_h_cap_positions(error)
      !! Test that H-cap positions match replaced atom positions
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      type(physical_fragment_t) :: fragment
      type(bond_t), allocatable :: bonds(:)
      type(error_t) :: parse_error
      integer :: stat, i
      character(len=:), allocatable :: errmsg
      real(dp) :: replaced_atom_coords(3), cap_coords(3)

      call create_test_water_trimer()
      call initialize_system_geometry("test_water_trimer.xyz", "test_water_monomer.xyz", &
                                      sys_geom, parse_error)

      if (parse_error%has_error()) then
         call check(error, .false., "Failed to initialize system: "//parse_error%get_message())
         call cleanup_test_files()
         return
      end if

      ! Create broken bond between atom 0 (monomer 1) and atom 3 (monomer 2)
      allocate (bonds(1))
      bonds(1)%atom_i = 0
      bonds(1)%atom_j = 3  ! This is 0-indexed
      bonds(1)%order = 1
      bonds(1)%is_broken = .true.

      ! Get coordinates of atom 3 before building fragment
      ! atom_j is 0-indexed (3), so in Fortran array it's at index 4
      replaced_atom_coords = sys_geom%coordinates(:, 4)

      ! Build fragment from monomer 1 only
      call build_fragment_from_indices(sys_geom, [1], fragment, parse_error, bonds)

      ! H-cap should be at position 4 (last atom)
      cap_coords = fragment%coordinates(:, 4)

      ! Check that H-cap coordinates match replaced atom coordinates
      do i = 1, 3
         call check(error, abs(cap_coords(i) - replaced_atom_coords(i)) < tol, &
                    "H-cap coordinate should match replaced atom")
         if (allocated(error)) return
      end do

      ! Also verify element number is hydrogen
      call check(error, fragment%element_numbers(4), 1, &
                 "H-cap element should be hydrogen (1)")

      call fragment%destroy()
      call sys_geom%destroy()
      call cleanup_test_files()
      if (allocated(bonds)) deallocate (bonds)
   end subroutine test_h_cap_positions

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

      open (newunit=unit, file=filename, status="old", action="read", iostat=stat)
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
