module test_mqc_interaction_bonding
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp, int64
   use mqc_error, only: error_t
   use mqc_result_types, only: calculation_result_t
   use mqc_json_output_types, only: interaction_bonding_term_t
   use mqc_quao_rows, only: QUAO_ROW_BOND, QUAO_ROW_DELOCALIZATION, QUAO_TYPE_SIGMA, &
                            QUAO_TYPE_LONE_PAIR
   use mqc_physical_fragment, only: system_geometry_t, bond_t
   use mqc_interaction_bonding, only: collect_interaction_bonding
   implicit none
   private
   public :: collect_mqc_interaction_bonding_tests

   integer, parameter :: REFERENCE = 3
      !! Monomer 3, atoms 4 and 5 (0-based), is the reference throughout

contains

   subroutine collect_mqc_interaction_bonding_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("only_crossing_rows_survive", test_crossing_rows), &
                  new_unittest("terms_without_rows_are_skipped", test_term_selection), &
                  new_unittest("atom_pairs_sorted_and_thresholded", test_atom_pairs) &
                  ]
   end subroutine collect_mqc_interaction_bonding_tests

   subroutine three_monomers(sys_geom)
      !! Six carbons on a line, two per monomer, with the bond between atoms 1
      !! and 2 (0-based) cut by the partition
      !!
      !! A term of monomers 1 and 3 therefore carries one cap, on atom 1's side,
      !! as local atom 5 after the four real atoms.
      type(system_geometry_t), intent(out) :: sys_geom

      integer :: a

      sys_geom%total_atoms = 6
      sys_geom%n_monomers = 3
      sys_geom%atoms_per_monomer = 2
      sys_geom%charge = 0
      sys_geom%multiplicity = 1
      allocate (sys_geom%element_numbers(6), sys_geom%coordinates(3, 6))
      sys_geom%element_numbers = 6
      sys_geom%coordinates = 0.0_dp
      do a = 1, 6
         sys_geom%coordinates(1, a) = 2.9_dp*real(a - 1, dp)
      end do
      allocate (sys_geom%bonds(1))
      sys_geom%bonds(1) = bond_t(atom_i=1, atom_j=2, order=1, is_broken=.true.)
   end subroutine three_monomers

   subroutine term_rows(result)
      !! Four orbital pairs over the 5 local atoms of the term (1, 3)
      !!
      !! Local atoms 1-2 are monomer 1, 3-4 the reference, 5 the cap.
      !!   row 1  bond 1-2, both on monomer 1                 -> not crossing
      !!   row 2  3 donates into 1, 1's partner the cap        -> kept
      !!   row 3  2 into the cap, both on monomer 1's side     -> not crossing
      !!   row 4  4 into the cap: crosses, but reaches a cap   -> omitted
      type(calculation_result_t), intent(inout) :: result

      call result%quao_rows%unpack( &
         [QUAO_ROW_BOND, 3, 4, 1, 2, QUAO_TYPE_SIGMA, QUAO_TYPE_SIGMA, 1, 1, 2, 1, &
          QUAO_ROW_DELOCALIZATION, 7, 5, 3, 1, QUAO_TYPE_LONE_PAIR, QUAO_TYPE_SIGMA, 1, 1, 0, 5, &
          QUAO_ROW_DELOCALIZATION, 6, 9, 2, 5, QUAO_TYPE_SIGMA, QUAO_TYPE_SIGMA, 1, 0, 1, 2, &
          QUAO_ROW_DELOCALIZATION, 8, 9, 4, 5, QUAO_TYPE_LONE_PAIR, QUAO_TYPE_SIGMA, 1, 0, 0, 2], &
         [1.0_dp, 1.0_dp, 0.95_dp, -60.0_dp, &
          1.9_dp, 1.0_dp, 0.10_dp, -2.5_dp, &
          1.0_dp, 1.0_dp, 0.20_dp, -4.0_dp, &
          1.9_dp, 1.0_dp, 0.30_dp, -5.0_dp])
      allocate (result%quao_rows%atom_bond_index(5, 5), &
                result%quao_rows%atom_kinetic_bond_order(5, 5))
      result%quao_rows%atom_bond_index = 0.0_dp
      result%quao_rows%atom_kinetic_bond_order = 0.0_dp
      ! Reference atoms are local 3 and 4. Local 3 against 1 and 2 above the
      ! threshold, 4 against 2 below it, 4 against the cap large but a cap.
      result%quao_rows%atom_bond_index(3, 1) = 0.02_dp
      result%quao_rows%atom_kinetic_bond_order(3, 1) = -1.5_dp
      result%quao_rows%atom_bond_index(3, 2) = 0.04_dp
      result%quao_rows%atom_kinetic_bond_order(3, 2) = -3.0_dp
      result%quao_rows%atom_kinetic_bond_order(4, 2) = -0.5_dp
      result%quao_rows%atom_kinetic_bond_order(4, 5) = -9.0_dp
      result%quao_rows%threshold = 1.0_dp
      result%has_quao_rows = .true.
   end subroutine term_rows

   subroutine collected(terms, err)
      !! Four terms: (1,3) with rows, (2,3) without, (3) alone, (1,2) without
      !! the reference. Only the first should come back.
      type(interaction_bonding_term_t), allocatable, intent(out) :: terms(:)
      type(error_t), intent(out) :: err

      type(system_geometry_t) :: sys_geom
      type(calculation_result_t) :: results(4)
      integer :: polymers(4, 2)

      call three_monomers(sys_geom)
      polymers(1, :) = [1, 3]
      polymers(2, :) = [2, 3]
      polymers(3, :) = [3, 0]
      polymers(4, :) = [1, 2]
      call term_rows(results(1))
      call term_rows(results(3))
      call term_rows(results(4))
      call collect_interaction_bonding(polymers, 4_int64, REFERENCE, results, sys_geom, &
                                       terms, err)
   end subroutine collected

   subroutine test_term_selection(error)
      !! A term is reported only when it is real, holds the reference and
      !! another monomer, and its calculation returned rows
      type(error_type), allocatable, intent(out) :: error
      type(interaction_bonding_term_t), allocatable :: terms(:)
      type(error_t) :: err

      call collected(terms, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, allocated(terms), "the term (1, 3) should have been reported")
      if (allocated(error)) return
      call check(error, size(terms), 1, "only the term (1, 3) qualifies")
      if (allocated(error)) return
      call check(error, terms(1)%term, 1, "and it is row 1 of the term list")
      if (allocated(error)) return
      call check(error, all(terms(1)%monomers == [1, 3]), "with its own monomers")
   end subroutine test_term_selection

   subroutine test_crossing_rows(error)
      !! Rows with both ends on one side are dropped silently; a crossing row
      !! reaching a cap is dropped and counted; the survivor is renumbered
      type(error_type), allocatable, intent(out) :: error
      type(interaction_bonding_term_t), allocatable :: terms(:)
      type(error_t) :: err

      call collected(terms, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      associate (t => terms(1))
         call check(error, t%rows%n, 1, "one orbital pair crosses cleanly")
         if (allocated(error)) return
         call check(error, t%omitted_rows, 1, "one crossing pair reached the cap")
         if (allocated(error)) return
         ! Local 3 is the reference's first atom, system atom 5 (1-based);
         ! local 1 is system atom 1.
         call check(error, all(t%rows%atom(:, 1) == [5, 1]), "atoms renumbered to the system")
         if (allocated(error)) return
         call check(error, all(t%monomer_of_atom(:, 1) == [REFERENCE, 1]), &
                    "donor on the reference, acceptor on monomer 1")
         if (allocated(error)) return
         call check(error, t%rows%orbital(1, 1), 7, "orbital numbers carried through")
         if (allocated(error)) return
         call check(error, t%rows%partner_atom(1, 1), 0, "a lone pair has no partner")
         if (allocated(error)) return
         call check(error, t%rows%partner_atom(2, 1), 0, &
                    "a partner that is a cap is not a system atom, so none is named")
         if (allocated(error)) return
         call check(error, t%rows%kinetic_bond_order(1), -2.5_dp, "values carried through")
      end associate
   end subroutine test_crossing_rows

   subroutine test_atom_pairs(error)
      !! Reference atom against environment atom, at the threshold, strongest
      !! bonding first, never against a cap
      type(error_type), allocatable, intent(out) :: error
      type(interaction_bonding_term_t), allocatable :: terms(:)
      type(error_t) :: err

      call collected(terms, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      associate (t => terms(1))
         call check(error, size(t%pair_kinetic_bond_order), 2, &
                    "two pairs reach 1 kcal/mol; the cap's -9 is not a pair")
         if (allocated(error)) return
         call check(error, all(t%pair_atoms(:, 1) == [5, 2]), "the -3.0 pair first")
         if (allocated(error)) return
         call check(error, t%pair_kinetic_bond_order(1), -3.0_dp, "strongest first")
         if (allocated(error)) return
         call check(error, t%pair_bond_index(1), 0.04_dp, "its bond index")
         if (allocated(error)) return
         call check(error, all(t%pair_atoms(:, 2) == [5, 1]), "then the -1.5 pair")
      end associate
   end subroutine test_atom_pairs

end module test_mqc_interaction_bonding

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_interaction_bonding, only: collect_mqc_interaction_bonding_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_interaction_bonding", collect_mqc_interaction_bonding_tests) &
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
