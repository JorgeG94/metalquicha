module test_mqc_gmbe_intersection
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_frag_utils, only: find_fragment_intersection, generate_intersections, &
                             compute_polymer_atoms, generate_polymer_intersections
   use mqc_physical_fragment, only: system_geometry_t
   use pic_types, only: default_int, dp
   implicit none
   private
   public :: collect_mqc_gmbe_intersection_tests

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_gmbe_intersection_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("no_intersection", test_no_intersection), &
                  new_unittest("single_atom_intersection", test_single_atom_intersection), &
                  new_unittest("multiple_atom_intersection", test_multiple_atom_intersection), &
                  new_unittest("complete_overlap", test_complete_overlap), &
                  new_unittest("partial_overlap_beginning", test_partial_overlap_beginning), &
                  new_unittest("partial_overlap_end", test_partial_overlap_end), &
                  new_unittest("generate_no_overlaps", test_generate_no_overlaps), &
                  new_unittest("generate_single_overlap", test_generate_single_overlap), &
                  new_unittest("generate_multiple_overlaps", test_generate_multiple_overlaps), &
                  new_unittest("compute_polymer_atoms_dimer", test_compute_polymer_atoms_dimer), &
                  new_unittest("compute_polymer_atoms_overlapping", test_compute_polymer_atoms_overlapping), &
                  new_unittest("generate_polymer_intersections_dimers", test_generate_polymer_intersections_dimers) &
                  ]
   end subroutine collect_mqc_gmbe_intersection_tests

   subroutine test_no_intersection(error)
      !! Test two fragments with no shared atoms
      type(error_type), allocatable, intent(out) :: error
      integer :: frag1(4), frag2(4)
      integer, allocatable :: intersection(:)
      integer :: n_intersect
      logical :: has_intersection

      ! Fragment 1: atoms [0, 1, 2, 3]
      frag1 = [0, 1, 2, 3]

      ! Fragment 2: atoms [4, 5, 6, 7] (no overlap)
      frag2 = [4, 5, 6, 7]

      has_intersection = find_fragment_intersection(frag1, size(frag1), &
                                                    frag2, size(frag2), &
                                                    intersection, n_intersect)

      call check(error, has_intersection, .false., &
                 "Should have no intersection between non-overlapping fragments")
      if (allocated(error)) return

      call check(error, n_intersect, 0, &
                 "Intersection count should be 0")
      if (allocated(error)) return

   end subroutine test_no_intersection

   subroutine test_single_atom_intersection(error)
      !! Test two fragments with one shared atom
      type(error_type), allocatable, intent(out) :: error
      integer :: frag1(4), frag2(4)
      integer, allocatable :: intersection(:)
      integer :: n_intersect
      logical :: has_intersection

      ! Fragment 1: atoms [0, 1, 2, 3]
      frag1 = [0, 1, 2, 3]

      ! Fragment 2: atoms [3, 4, 5, 6] (shares atom 3)
      frag2 = [3, 4, 5, 6]

      has_intersection = find_fragment_intersection(frag1, size(frag1), &
                                                    frag2, size(frag2), &
                                                    intersection, n_intersect)

      call check(error, has_intersection, .true., &
                 "Should have intersection with one shared atom")
      if (allocated(error)) return

      call check(error, n_intersect, 1, &
                 "Intersection should contain 1 atom")
      if (allocated(error)) return

      call check(error, intersection(1), 3, &
                 "Shared atom should be atom 3")
      if (allocated(error)) return

      deallocate (intersection)

   end subroutine test_single_atom_intersection

   subroutine test_multiple_atom_intersection(error)
      !! Test two fragments with multiple shared atoms (polypeptide-like)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nfrag = 6
      integer :: frag1(nfrag), frag2(nfrag)
      integer, allocatable :: intersection(:)
      integer :: n_intersect
      logical :: has_intersection

      ! Fragment 1: atoms [0, 1, 2, 3, 4, 5]
      frag1 = [0, 1, 2, 3, 4, 5]

      ! Fragment 2: atoms [3, 4, 5, 6, 7, 8] (shares atoms 3, 4, 5)
      frag2 = [3, 4, 5, 6, 7, 8]

      has_intersection = find_fragment_intersection(frag1, size(frag1), &
                                                    frag2, size(frag2), &
                                                    intersection, n_intersect)

      call check(error, has_intersection, .true., &
                 "Should have intersection with three shared atoms")
      if (allocated(error)) return

      call check(error, n_intersect, 3, &
                 "Intersection should contain 3 atoms")
      if (allocated(error)) return

      ! Check that all three shared atoms are present
      call check(error, any(intersection(1:n_intersect) == 3), .true., &
                 "Intersection should contain atom 3")
      if (allocated(error)) then
         deallocate (intersection)
         return
      end if

      call check(error, any(intersection(1:n_intersect) == 4), .true., &
                 "Intersection should contain atom 4")
      if (allocated(error)) then
         deallocate (intersection)
         return
      end if

      call check(error, any(intersection(1:n_intersect) == 5), .true., &
                 "Intersection should contain atom 5")
      if (allocated(error)) then
         deallocate (intersection)
         return
      end if

      deallocate (intersection)

   end subroutine test_multiple_atom_intersection

   subroutine test_complete_overlap(error)
      !! Test two identical fragments (complete overlap)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nfrag = 4
      integer :: frag1(nfrag), frag2(nfrag)
      integer, allocatable :: intersection(:)
      integer :: n_intersect
      logical :: has_intersection

      ! Both fragments: atoms [0, 1, 2, 3]
      frag1 = [0, 1, 2, 3]
      frag2 = [0, 1, 2, 3]

      has_intersection = find_fragment_intersection(frag1, size(frag1), &
                                                    frag2, size(frag2), &
                                                    intersection, n_intersect)

      call check(error, has_intersection, .true., &
                 "Should have intersection for identical fragments")
      if (allocated(error)) return

      call check(error, n_intersect, 4, &
                 "Complete overlap should have 4 atoms")
      if (allocated(error)) return

      deallocate (intersection)

   end subroutine test_complete_overlap

   subroutine test_partial_overlap_beginning(error)
      !! Test fragments that overlap at the beginning of frag2
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nfrag = 5
      integer :: frag1(nfrag), frag2(nfrag)
      integer, allocatable :: intersection(:)
      integer :: n_intersect
      logical :: has_intersection

      ! Fragment 1: atoms [5, 6, 7, 8, 9]
      frag1 = [5, 6, 7, 8, 9]

      ! Fragment 2: atoms [8, 9, 10, 11, 12] (shares 8, 9 at beginning)
      frag2 = [8, 9, 10, 11, 12]

      has_intersection = find_fragment_intersection(frag1, size(frag1), &
                                                    frag2, size(frag2), &
                                                    intersection, n_intersect)

      call check(error, has_intersection, .true., &
                 "Should have intersection")
      if (allocated(error)) return

      call check(error, n_intersect, 2, &
                 "Should have 2 shared atoms")
      if (allocated(error)) return

      deallocate (intersection)

   end subroutine test_partial_overlap_beginning

   subroutine test_partial_overlap_end(error)
      !! Test fragments that overlap at the end of frag1
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nfrag = 5
      integer :: frag1(nfrag), frag2(nfrag)
      integer, allocatable :: intersection(:)
      integer :: n_intersect
      logical :: has_intersection

      ! Fragment 1: atoms [0, 1, 2, 3, 4]
      frag1 = [0, 1, 2, 3, 4]

      ! Fragment 2: atoms [3, 4, 5, 6, 7] (shares 3, 4 at end of frag1)
      frag2 = [3, 4, 5, 6, 7]

      has_intersection = find_fragment_intersection(frag1, size(frag1), &
                                                    frag2, size(frag2), &
                                                    intersection, n_intersect)

      call check(error, has_intersection, .true., &
                 "Should have intersection")
      if (allocated(error)) return

      call check(error, n_intersect, 2, &
                 "Should have 2 shared atoms")
      if (allocated(error)) return

      call check(error, any(intersection(1:n_intersect) == 3), .true., &
                 "Should contain atom 3")
      if (allocated(error)) then
         deallocate (intersection)
         return
      end if

      call check(error, any(intersection(1:n_intersect) == 4), .true., &
                 "Should contain atom 4")
      if (allocated(error)) then
         deallocate (intersection)
         return
      end if

      deallocate (intersection)

   end subroutine test_partial_overlap_end

   subroutine test_generate_no_overlaps(error)
      !! Test generate_intersections with non-overlapping fragments
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      integer, allocatable :: monomers(:), polymers(:, :)
      integer, allocatable :: intersections(:, :), intersection_sets(:, :), intersection_levels(:)
      integer :: n_intersections

      ! Create a simple system with 2 non-overlapping fragments
      sys_geom%n_monomers = 2
      sys_geom%total_atoms = 8

      ! Allocate fragment info
      allocate (sys_geom%fragment_sizes(2))
      allocate (sys_geom%fragment_atoms(4, 2))

      ! Fragment 1: atoms [0, 1, 2, 3]
      sys_geom%fragment_sizes(1) = 4
      sys_geom%fragment_atoms(1:4, 1) = [0, 1, 2, 3]

      ! Fragment 2: atoms [4, 5, 6, 7]
      sys_geom%fragment_sizes(2) = 4
      sys_geom%fragment_atoms(1:4, 2) = [4, 5, 6, 7]

      ! Generate intersections
      allocate (monomers(2))
      monomers = [1, 2]
      allocate (polymers(2, 1))

      call generate_intersections(sys_geom, monomers, polymers, 2, 999, &
                                  intersections, intersection_sets, intersection_levels, n_intersections)

      ! Should have 0 intersections
      call check(error, n_intersections, 0, &
                 "Non-overlapping fragments should produce 0 intersections")
      if (allocated(error)) return

      ! Clean up
      deallocate (monomers, polymers)
      call sys_geom%destroy()

   end subroutine test_generate_no_overlaps

   subroutine test_generate_single_overlap(error)
      !! Test generate_intersections with two overlapping fragments
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      integer, allocatable :: monomers(:), polymers(:, :)
      integer, allocatable :: intersections(:, :), intersection_sets(:, :), intersection_levels(:)
      integer :: n_intersections

      ! Create a system with 2 overlapping fragments (polypeptide-like)
      sys_geom%n_monomers = 2
      sys_geom%total_atoms = 9

      ! Allocate fragment info
      allocate (sys_geom%fragment_sizes(2))
      allocate (sys_geom%fragment_atoms(6, 2))

      ! Fragment 1: atoms [0, 1, 2, 3, 4, 5]
      sys_geom%fragment_sizes(1) = 6
      sys_geom%fragment_atoms(1:6, 1) = [0, 1, 2, 3, 4, 5]

      ! Fragment 2: atoms [3, 4, 5, 6, 7, 8] (shares 3, 4, 5)
      sys_geom%fragment_sizes(2) = 6
      sys_geom%fragment_atoms(1:6, 2) = [3, 4, 5, 6, 7, 8]

      ! Generate intersections
      allocate (monomers(2))
      monomers = [1, 2]
      allocate (polymers(2, 1))

      call generate_intersections(sys_geom, monomers, polymers, 2, 999, &
                                  intersections, intersection_sets, intersection_levels, n_intersections)

      ! Should have 1 intersection
      call check(error, n_intersections, 1, &
                 "Two overlapping fragments should produce 1 intersection")
      if (allocated(error)) goto 100

      ! Check intersection level (should be 2 for pairwise)
      call check(error, intersection_levels(1), 2, &
                 "Intersection level should be 2 for pairwise")
      if (allocated(error)) goto 100

      ! Check intersection sets
      call check(error, intersection_sets(1, 1), 1, &
                 "First fragment in set should be 1")
      if (allocated(error)) goto 100

      call check(error, intersection_sets(2, 1), 2, &
                 "Second fragment in set should be 2")
      if (allocated(error)) goto 100

      ! Check intersection atoms (should be 3, 4, 5)
      ! Order doesn't matter, just check count
      call check(error, size(intersections, 1) >= 3, .true., &
                 "Intersection should have space for at least 3 atoms")
      if (allocated(error)) goto 100

100   continue
      ! Clean up
      deallocate (monomers, polymers)
      if (allocated(intersections)) deallocate (intersections)
      if (allocated(intersection_sets)) deallocate (intersection_sets)
      if (allocated(intersection_levels)) deallocate (intersection_levels)
      call sys_geom%destroy()

   end subroutine test_generate_single_overlap

   subroutine test_generate_multiple_overlaps(error)
      !! Test generate_intersections with three fragments, multiple overlaps
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      integer, allocatable :: monomers(:), polymers(:, :)
      integer, allocatable :: intersections(:, :), intersection_sets(:, :), intersection_levels(:)
      integer :: n_intersections

      ! Create a system with 3 fragments in a chain (like tripeptide)
      ! frag1: [0,1,2,3,4,5], frag2: [3,4,5,6,7,8], frag3: [6,7,8,9,10,11]
      sys_geom%n_monomers = 3
      sys_geom%total_atoms = 12

      ! Allocate fragment info
      allocate (sys_geom%fragment_sizes(3))
      allocate (sys_geom%fragment_atoms(6, 3))

      ! Fragment 1: atoms [0, 1, 2, 3, 4, 5]
      sys_geom%fragment_sizes(1) = 6
      sys_geom%fragment_atoms(1:6, 1) = [0, 1, 2, 3, 4, 5]

      ! Fragment 2: atoms [3, 4, 5, 6, 7, 8] (overlaps with 1 and 3)
      sys_geom%fragment_sizes(2) = 6
      sys_geom%fragment_atoms(1:6, 2) = [3, 4, 5, 6, 7, 8]

      ! Fragment 3: atoms [6, 7, 8, 9, 10, 11] (overlaps with 2)
      sys_geom%fragment_sizes(3) = 6
      sys_geom%fragment_atoms(1:6, 3) = [6, 7, 8, 9, 10, 11]

      ! Generate intersections
      allocate (monomers(3))
      monomers = [1, 2, 3]
      allocate (polymers(3, 1))

      call generate_intersections(sys_geom, monomers, polymers, 3, 999, &
                                  intersections, intersection_sets, intersection_levels, n_intersections)

      ! Should have 2 intersections: (1,2) and (2,3)
      ! Note: (1,3) should have no intersection
      call check(error, n_intersections, 2, &
                 "Three fragments in chain should produce 2 intersections")
      if (allocated(error)) goto 200

200   continue
      ! Clean up
      deallocate (monomers, polymers)
      if (allocated(intersections)) deallocate (intersections)
      if (allocated(intersection_sets)) deallocate (intersection_sets)
      if (allocated(intersection_levels)) deallocate (intersection_levels)
      call sys_geom%destroy()

   end subroutine test_generate_multiple_overlaps

   subroutine test_compute_polymer_atoms_dimer(error)
      !! Test compute_polymer_atoms with non-overlapping fragments forming a dimer
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      integer :: polymer(2)
      integer, allocatable :: atom_list(:)
      integer :: n_atoms

      ! Create system with 2 non-overlapping fragments
      sys_geom%n_monomers = 2
      sys_geom%total_atoms = 8

      allocate (sys_geom%fragment_sizes(2))
      allocate (sys_geom%fragment_atoms(4, 2))

      ! Fragment 1: atoms [0, 1, 2, 3]
      sys_geom%fragment_sizes(1) = 4
      sys_geom%fragment_atoms(1:4, 1) = [0, 1, 2, 3]

      ! Fragment 2: atoms [4, 5, 6, 7]
      sys_geom%fragment_sizes(2) = 4
      sys_geom%fragment_atoms(1:4, 2) = [4, 5, 6, 7]

      ! Dimer: (1, 2)
      polymer = [1, 2]

      call compute_polymer_atoms(sys_geom, polymer, 2, atom_list, n_atoms)

      ! Should have 8 atoms (union of both fragments)
      call check(error, n_atoms, 8, &
                 "Dimer of non-overlapping fragments should have 8 atoms")
      if (allocated(error)) goto 100

      ! Check all atoms are present
      call check(error, any(atom_list == 0), .true., "Should contain atom 0")
      if (allocated(error)) goto 100
      call check(error, any(atom_list == 7), .true., "Should contain atom 7")
      if (allocated(error)) goto 100

100   continue
      if (allocated(atom_list)) deallocate (atom_list)
      call sys_geom%destroy()

   end subroutine test_compute_polymer_atoms_dimer

   subroutine test_compute_polymer_atoms_overlapping(error)
      !! Test compute_polymer_atoms with overlapping fragments
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      integer :: polymer(2)
      integer, allocatable :: atom_list(:)
      integer :: n_atoms

      ! Create system with 2 overlapping fragments
      sys_geom%n_monomers = 2
      sys_geom%total_atoms = 9

      allocate (sys_geom%fragment_sizes(2))
      allocate (sys_geom%fragment_atoms(6, 2))

      ! Fragment 1: atoms [0, 1, 2, 3, 4, 5]
      sys_geom%fragment_sizes(1) = 6
      sys_geom%fragment_atoms(1:6, 1) = [0, 1, 2, 3, 4, 5]

      ! Fragment 2: atoms [3, 4, 5, 6, 7, 8] (shares 3, 4, 5)
      sys_geom%fragment_sizes(2) = 6
      sys_geom%fragment_atoms(1:6, 2) = [3, 4, 5, 6, 7, 8]

      ! Dimer: (1, 2)
      polymer = [1, 2]

      call compute_polymer_atoms(sys_geom, polymer, 2, atom_list, n_atoms)

      ! Should have 9 unique atoms (6 + 6 - 3 overlap = 9)
      call check(error, n_atoms, 9, &
                 "Dimer of overlapping fragments should have 9 unique atoms")
      if (allocated(error)) goto 100

      ! Verify boundary atoms
      call check(error, any(atom_list == 0), .true., "Should contain atom 0")
      if (allocated(error)) goto 100
      call check(error, any(atom_list == 8), .true., "Should contain atom 8")
      if (allocated(error)) goto 100

100   continue
      if (allocated(atom_list)) deallocate (atom_list)
      call sys_geom%destroy()

   end subroutine test_compute_polymer_atoms_overlapping

   subroutine test_generate_polymer_intersections_dimers(error)
      !! Test generate_polymer_intersections with overlapping dimers
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys_geom
      integer, allocatable :: polymers(:, :)
      integer, allocatable :: intersections(:, :), intersection_sets(:, :), intersection_levels(:)
      integer :: n_intersections

      ! Create system with 3 overlapping fragments (tripeptide-like)
      sys_geom%n_monomers = 3
      sys_geom%total_atoms = 12

      allocate (sys_geom%fragment_sizes(3))
      allocate (sys_geom%fragment_atoms(6, 3))

      ! Fragment 1: atoms [0, 1, 2, 3, 4, 5]
      sys_geom%fragment_sizes(1) = 6
      sys_geom%fragment_atoms(1:6, 1) = [0, 1, 2, 3, 4, 5]

      ! Fragment 2: atoms [3, 4, 5, 6, 7, 8]
      sys_geom%fragment_sizes(2) = 6
      sys_geom%fragment_atoms(1:6, 2) = [3, 4, 5, 6, 7, 8]

      ! Fragment 3: atoms [6, 7, 8, 9, 10, 11]
      sys_geom%fragment_sizes(3) = 6
      sys_geom%fragment_atoms(1:6, 3) = [6, 7, 8, 9, 10, 11]

      ! Generate dimers: (1,2), (1,3), (2,3)
      allocate (polymers(3, 2))
      polymers(1, :) = [1, 2]
      polymers(2, :) = [1, 3]
      polymers(3, :) = [2, 3]

      ! Generate intersections between dimers
      call generate_polymer_intersections(sys_geom, polymers, 3, 2, &
                                          intersections, intersection_sets, intersection_levels, n_intersections)

      ! Dimer(1,2) = {0,1,2,3,4,5,6,7,8}
      ! Dimer(1,3) = {0,1,2,3,4,5,6,7,8,9,10,11}
      ! Dimer(2,3) = {3,4,5,6,7,8,9,10,11}
      ! Intersection (1,2)∩(1,3) = {0,1,2,3,4,5,6,7,8}
      ! Intersection (1,2)∩(2,3) = {3,4,5,6,7,8}
      ! Intersection (2,3)∩(1,3) = {3,4,5,6,7,8,9,10,11}
      ! Should have 3 pairwise intersections

      call check(error, n_intersections >= 3, .true., &
                 "Should have at least 3 intersections between 3 dimers")
      if (allocated(error)) goto 100

100   continue
      deallocate (polymers)
      if (allocated(intersections)) deallocate (intersections)
      if (allocated(intersection_sets)) deallocate (intersection_sets)
      if (allocated(intersection_levels)) deallocate (intersection_levels)
      call sys_geom%destroy()

   end subroutine test_generate_polymer_intersections_dimers

end module test_mqc_gmbe_intersection

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_gmbe_intersection, only: collect_mqc_gmbe_intersection_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_gmbe_intersection", collect_mqc_gmbe_intersection_tests) &
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
