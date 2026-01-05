module test_mqc_fragment_from_atoms
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_physical_fragment, only: physical_fragment_t, system_geometry_t, build_fragment_from_atom_list
   use mqc_config_parser, only: bond_t
   use pic_types, only: dp
   use mqc_error, only: error_t
   implicit none
   private
   public :: collect_mqc_fragment_from_atoms_tests

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_fragment_from_atoms_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("simple_intersection_no_caps", test_simple_intersection_no_caps), &
                  new_unittest("intersection_with_broken_bonds", test_intersection_with_broken_bonds), &
                  new_unittest("intersection_charge_neutral", test_intersection_charge_neutral) &
                  ]
   end subroutine collect_mqc_fragment_from_atoms_tests

   subroutine test_simple_intersection_no_caps(error)
      !! Test building intersection fragment without broken bonds
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: my_error
      type(system_geometry_t) :: sys_geom
      type(physical_fragment_t) :: fragment
      integer :: atom_indices(3)

      ! Create simple system with 6 atoms (e.g., two overlapping fragments)
      sys_geom%total_atoms = 6
      sys_geom%n_monomers = 2
      allocate (sys_geom%element_numbers(6))
      allocate (sys_geom%coordinates(3, 6))

      ! Atoms: H, C, N, C, N, H (0-indexed: 0,1,2,3,4,5)
      sys_geom%element_numbers = [1, 6, 7, 6, 7, 1]
      sys_geom%coordinates = reshape([ &
                                     0.0_dp, 0.0_dp, 0.0_dp, &  ! H (0)
                                     1.0_dp, 0.0_dp, 0.0_dp, &  ! C (1)
                                     2.0_dp, 0.0_dp, 0.0_dp, &  ! N (2)
                                     3.0_dp, 0.0_dp, 0.0_dp, &  ! C (3)
                                     4.0_dp, 0.0_dp, 0.0_dp, &  ! N (4)
                                     5.0_dp, 0.0_dp, 0.0_dp &   ! H (5)
                                     ], [3, 6])

      ! Build intersection fragment from atoms [1, 2, 3] (C-N-C)
      atom_indices = [1, 2, 3]
      call build_fragment_from_atom_list(sys_geom, atom_indices, 3, fragment, my_error)

      ! Check basic properties
      call check(error, fragment%n_atoms, 3, "Fragment should have 3 atoms")
      if (allocated(error)) return

      call check(error, fragment%n_caps, 0, "Fragment should have 0 caps (no bonds provided)")
      if (allocated(error)) return

      call check(error, fragment%element_numbers(1), 6, "First atom should be C")
      if (allocated(error)) return

      call check(error, fragment%element_numbers(2), 7, "Second atom should be N")
      if (allocated(error)) return

      call check(error, fragment%element_numbers(3), 6, "Third atom should be C")
      if (allocated(error)) return

      ! Check intersection is NEUTRAL
      call check(error, fragment%charge, 0, "Intersection fragment should have charge=0")
      if (allocated(error)) return
      call check(error, fragment%multiplicity, 1, "Intersection fragment should have multiplicity=1")
      if (allocated(error)) return

      ! Check electron count (C=6, N=7, C=6, total=19 electrons for neutral)
      call check(error, fragment%nelec, 19, "Neutral C-N-C should have 19 electrons")
      if (allocated(error)) return

      call fragment%destroy()
      call sys_geom%destroy()

   end subroutine test_simple_intersection_no_caps

   subroutine test_intersection_with_broken_bonds(error)
      !! Test building intersection fragment with hydrogen caps for broken bonds
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: my_error
      type(system_geometry_t) :: sys_geom
      type(physical_fragment_t) :: fragment
      type(bond_t), allocatable :: bonds(:)
      integer :: atom_indices(3)

      ! Create system with 6 atoms: H-C-N-C-N-H
      sys_geom%total_atoms = 6
      sys_geom%n_monomers = 2
      allocate (sys_geom%element_numbers(6))
      allocate (sys_geom%coordinates(3, 6))

      sys_geom%element_numbers = [1, 6, 7, 6, 7, 1]
      sys_geom%coordinates = reshape([ &
                                     0.0_dp, 0.0_dp, 0.0_dp, &  ! H (0)
                                     1.0_dp, 0.0_dp, 0.0_dp, &  ! C (1)
                                     2.0_dp, 0.0_dp, 0.0_dp, &  ! N (2)
                                     3.0_dp, 0.0_dp, 0.0_dp, &  ! C (3)
                                     4.0_dp, 0.0_dp, 0.0_dp, &  ! N (4)
                                     5.0_dp, 0.0_dp, 0.0_dp &   ! H (5)
                                     ], [3, 6])

      ! Define bonds (0-indexed): 0-1, 1-2, 2-3, 3-4, 4-5
      ! Intersection will be [1, 2, 3], so bonds to 0 and 4 are broken
      allocate (bonds(5))
      bonds(1) = bond_t(atom_i=0, atom_j=1, is_broken=.true.)   ! Broken: 0 not in intersection
      bonds(2) = bond_t(atom_i=1, atom_j=2, is_broken=.false.)  ! Preserved: both in intersection
      bonds(3) = bond_t(atom_i=2, atom_j=3, is_broken=.false.)  ! Preserved: both in intersection
      bonds(4) = bond_t(atom_i=3, atom_j=4, is_broken=.true.)   ! Broken: 4 not in intersection
      bonds(5) = bond_t(atom_i=4, atom_j=5, is_broken=.false.)  ! Not involved

      ! Build intersection fragment from atoms [1, 2, 3] (C-N-C)
      atom_indices = [1, 2, 3]
      call build_fragment_from_atom_list(sys_geom, atom_indices, 3, fragment, my_error, bonds)

      ! Check we added 2 caps (for bonds 0-1 and 3-4)
      call check(error, fragment%n_caps, 2, "Should have 2 hydrogen caps")
      if (allocated(error)) return

      call check(error, fragment%n_atoms, 5, "Fragment should have 3 + 2 = 5 atoms")
      if (allocated(error)) return

      ! First 3 atoms should be C, N, C
      call check(error, fragment%element_numbers(1), 6, "Atom 1 should be C")
      if (allocated(error)) return
      call check(error, fragment%element_numbers(2), 7, "Atom 2 should be N")
      if (allocated(error)) return

      call check(error, fragment%element_numbers(3), 6, "Atom 3 should be C")
      if (allocated(error)) return

      ! Last 2 atoms should be H (caps)
      call check(error, fragment%element_numbers(4), 1, "Cap 1 should be H")
      if (allocated(error)) return

      call check(error, fragment%element_numbers(5), 1, "Cap 2 should be H")
      if (allocated(error)) return

      ! Check cap_replaces_atom
      call check(error, allocated(fragment%cap_replaces_atom), .true., &
                 "cap_replaces_atom should be allocated")
      if (allocated(error)) return

      ! Electron count: C(6) + N(7) + C(6) + H(1) + H(1) = 21 electrons
      call check(error, fragment%nelec, 21, "C-N-C with 2 H-caps should have 21 electrons")
      if (allocated(error)) return

      call fragment%destroy()
      call sys_geom%destroy()
      deallocate (bonds)

   end subroutine test_intersection_with_broken_bonds

   subroutine test_intersection_charge_neutral(error)
      !! Test that intersection fragments are always neutral, even if parent fragments are charged
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: my_error
      type(system_geometry_t) :: sys_geom
      type(physical_fragment_t) :: fragment
      integer :: atom_indices(2)

      ! Create system where fragments might be charged
      sys_geom%total_atoms = 4
      sys_geom%n_monomers = 2
      allocate (sys_geom%element_numbers(4))
      allocate (sys_geom%coordinates(3, 4))
      allocate (sys_geom%fragment_charges(2))
      allocate (sys_geom%fragment_multiplicities(2))

      sys_geom%element_numbers = [7, 1, 1, 1]  ! N-H-H-H (like NH3)
      sys_geom%coordinates = reshape([ &
                                     0.0_dp, 0.0_dp, 0.0_dp, &
                                     1.0_dp, 0.0_dp, 0.0_dp, &
                                     0.0_dp, 1.0_dp, 0.0_dp, &
                                     0.0_dp, 0.0_dp, 1.0_dp &
                                     ], [3, 4])

      ! Fragments might have charges (but intersections should not)
      sys_geom%fragment_charges = [1, 1]  ! Both charged +1
      sys_geom%fragment_multiplicities = [1, 1]
      sys_geom%charge = 2
      sys_geom%multiplicity = 1

      ! Build intersection from atoms [1, 2] (two hydrogens)
      atom_indices = [1, 2]
      call build_fragment_from_atom_list(sys_geom, atom_indices, 2, fragment, my_error)

      ! Intersection should ALWAYS be neutral
      call check(error, fragment%charge, 0, &
                 "Intersection should be neutral even if parent fragments are charged")
      if (allocated(error)) return

      call check(error, fragment%multiplicity, 1, &
                 "Intersection should have multiplicity=1")
      if (allocated(error)) return
      ! H + H = 2 electrons (neutral)
      call check(error, fragment%nelec, 2, "Two neutral H atoms should have 2 electrons")
      if (allocated(error)) return

      call fragment%destroy()
      call sys_geom%destroy()

   end subroutine test_intersection_charge_neutral

end module test_mqc_fragment_from_atoms

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fragment_from_atoms, only: collect_mqc_fragment_from_atoms_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_fragment_from_atoms", collect_mqc_fragment_from_atoms_tests) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester
