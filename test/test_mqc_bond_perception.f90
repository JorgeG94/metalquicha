module test_mqc_bond_perception
   !! Tests distance-based bond perception on geometries with a known answer.
   !!
   !! The two mistakes worth guarding against pull in opposite directions. Too
   !! generous a tolerance calls a hydrogen bond covalent, which for a water
   !! cluster fuses the monomers into one molecule and produces caps nobody
   !! wanted. Too tight misses a stretched bond in an unrelaxed structure and
   !! leaves a real cut undeclared. The water dimer below pins the first and
   !! the carbon chain the second.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_bond_perception, only: perceive_bonds, missing_broken_bonds, &
                                  auto_monomers, DEFAULT_BOND_TOLERANCE
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: system_geometry_t, bond_t, to_bohr
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_bond_perception_tests

contains

   subroutine collect_mqc_bond_perception_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("water_dimer_hydrogen_bond_is_not_covalent", test_water_dimer), &
                  new_unittest("carbon_chain_finds_every_backbone_bond", test_chain), &
                  new_unittest("cuts_are_marked_from_the_partition", test_cuts_marked), &
                  new_unittest("audit_finds_an_undeclared_cut", test_audit_finds), &
                  new_unittest("audit_is_quiet_when_all_cuts_declared", test_audit_quiet), &
                  new_unittest("declared_either_way_round_counts", test_audit_order), &
                  new_unittest("element_without_a_radius_bonds_to_nothing", test_no_radius), &
                  new_unittest("auto_monomers_splits_a_cluster", test_auto_cluster), &
                  new_unittest("auto_monomers_refuses_one_molecule", test_auto_refuses) &
                  ]
   end subroutine collect_mqc_bond_perception_tests

   subroutine test_water_dimer(error)
      !! Four real O-H bonds; the ~1.9 Angstrom hydrogen bond is not one
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(bond_t), allocatable :: bonds(:)
      integer :: n_bonds

      call water_dimer(sys)
      call perceive_bonds(sys, bonds, n_bonds)

      call check(error, n_bonds, 4, &
                 "a water dimer has four covalent bonds; a fifth means the "// &
                 "hydrogen bond was perceived as covalent")
      call sys%destroy()
   end subroutine test_water_dimer

   subroutine test_chain(error)
      !! Four carbons at 1.54 Angstrom: three bonds, and no 1-3 contacts
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(bond_t), allocatable :: bonds(:)
      integer :: n_bonds

      call carbon_chain(sys)
      call perceive_bonds(sys, bonds, n_bonds)

      call check(error, n_bonds, 3, &
                 "a four-carbon chain has three bonds; more means next-nearest "// &
                 "neighbours at 3.08 Angstrom were bonded")
      if (allocated(error)) return
      call check(error, bonds(1)%atom_i, 0, "bonds are reported 0-based")
      if (allocated(error)) return
      call check(error, bonds(1)%atom_j, 1)
      if (allocated(error)) return
      call check(error, bonds(1)%order, 1, "perception reports order 1 always")

      call sys%destroy()
   end subroutine test_chain

   subroutine test_cuts_marked(error)
      !! One carbon per monomer, so every backbone bond is a cut
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(bond_t), allocatable :: bonds(:)
      integer :: n_bonds, ibond

      call carbon_chain(sys)
      call partition_one_atom_each(sys, 4)
      call perceive_bonds(sys, bonds, n_bonds)

      call check(error, n_bonds, 3)
      if (allocated(error)) return
      do ibond = 1, n_bonds
         call check(error, bonds(ibond)%is_broken, &
                    "every bond crosses a monomer boundary and should be marked broken")
         if (allocated(error)) return
      end do

      ! The same geometry as one monomer cuts nothing.
      call partition_all_together(sys, 4)
      call perceive_bonds(sys, bonds, n_bonds)
      do ibond = 1, n_bonds
         call check(error,.not. bonds(ibond)%is_broken, &
                    "with one monomer no bond is cut")
         if (allocated(error)) return
      end do

      call sys%destroy()
   end subroutine test_cuts_marked

   subroutine test_audit_finds(error)
      !! The case the declaration check cannot see: a cut nobody mentioned
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(bond_t) :: declared(2)
      integer, allocatable :: missing_i(:), missing_j(:)
      integer :: n_missing

      call carbon_chain(sys)
      call partition_one_atom_each(sys, 4)

      ! Two of the three cuts declared; 2-3 left out entirely.
      declared(1)%atom_i = 0; declared(1)%atom_j = 1
      declared(1)%order = 1; declared(1)%is_broken = .true.
      declared(2)%atom_i = 1; declared(2)%atom_j = 2
      declared(2)%order = 1; declared(2)%is_broken = .true.

      call missing_broken_bonds(sys, declared, 2, missing_i, missing_j, n_missing)

      call check(error, n_missing, 1, "one cut was left undeclared")
      if (allocated(error)) return
      call check(error, missing_i(1), 2, "the missing cut starts at atom 2")
      if (allocated(error)) return
      call check(error, missing_j(1), 3, "and ends at atom 3")

      call sys%destroy()
   end subroutine test_audit_finds

   subroutine test_audit_quiet(error)
      !! A complete declaration produces no complaints
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(bond_t) :: declared(3)
      integer, allocatable :: missing_i(:), missing_j(:)
      integer :: n_missing, ibond

      call carbon_chain(sys)
      call partition_one_atom_each(sys, 4)
      do ibond = 1, 3
         declared(ibond)%atom_i = ibond - 1
         declared(ibond)%atom_j = ibond
         declared(ibond)%order = 1
         declared(ibond)%is_broken = .true.
      end do

      call missing_broken_bonds(sys, declared, 3, missing_i, missing_j, n_missing)
      call check(error, n_missing, 0, "nothing should be missing")

      call sys%destroy()
   end subroutine test_audit_quiet

   subroutine test_audit_order(error)
      !! A bond declared j-i is the same bond as i-j
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(bond_t) :: declared(3)
      integer, allocatable :: missing_i(:), missing_j(:)
      integer :: n_missing, ibond

      call carbon_chain(sys)
      call partition_one_atom_each(sys, 4)
      do ibond = 1, 3
         ! Reversed relative to how perception reports them.
         declared(ibond)%atom_i = ibond
         declared(ibond)%atom_j = ibond - 1
         declared(ibond)%order = 1
         declared(ibond)%is_broken = .true.
      end do

      call missing_broken_bonds(sys, declared, 3, missing_i, missing_j, n_missing)
      call check(error, n_missing, 0, &
                 "atom order within a bond should not matter")

      call sys%destroy()
   end subroutine test_audit_order

   subroutine test_no_radius(error)
      !! An element past the tabulated set bonds to nothing rather than guessing
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(bond_t), allocatable :: bonds(:)
      integer :: n_bonds

      call sys%destroy()
      sys%total_atoms = 2
      sys%n_monomers = 0
      allocate (sys%element_numbers(2))
      ! Oganesson: real element, no consensus covalent radius.
      sys%element_numbers = [118, 118]
      allocate (sys%coordinates(3, 2))
      sys%coordinates = 0.0_dp
      sys%coordinates(1, 2) = to_bohr(1.5_dp)

      call perceive_bonds(sys, bonds, n_bonds)
      call check(error, n_bonds, 0, &
                 "without a tabulated radius no bond should be invented")

      call sys%destroy()
   end subroutine test_no_radius

   subroutine test_auto_cluster(error)
      !! A hydrogen-bonded dimer comes apart into its two molecules
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(error_t) :: err

      call water_dimer(sys)
      call auto_monomers(sys, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      call check(error, sys%n_monomers, 2, "two waters, two monomers")
      if (allocated(error)) return
      call check(error, sys%atoms_per_monomer, 3, "each of three atoms")
      if (allocated(error)) return
      call check(error, all(sys%fragment_sizes == 3), "sizes agree")
      if (allocated(error)) return
      ! Atoms 0,1,2 are the first water; the partition should say so.
      call check(error, all(sys%fragment_atoms(1:3, 1) == [0, 1, 2]), &
                 "the first monomer holds the first three atoms")
      if (allocated(error)) return
      call check(error, all(sys%fragment_atoms(1:3, 2) == [3, 4, 5]), &
                 "the second holds the rest")
      if (allocated(error)) return
      call check(error, all(sys%fragment_charges == 0), "neutral by default")

      call sys%destroy()
   end subroutine test_auto_cluster

   subroutine test_auto_refuses(error)
      !! One connected molecule is exactly where the user has to choose
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(error_t) :: err

      call carbon_chain(sys)
      call auto_monomers(sys, err)

      call check(error, err%has_error(), &
                 "a single covalent molecule has no automatic partition")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "explicitly") > 0, &
                 "and the message should say to supply monomers explicitly")

      call sys%destroy()
   end subroutine test_auto_refuses

   ! ---- geometries ----------------------------------------------------------

   subroutine water_dimer(sys)
      !! A near-linear dimer: O-H ~0.96, and the O...H contact ~1.9 Angstrom
      type(system_geometry_t), intent(out) :: sys

      real(dp) :: xyz(3, 6)

      xyz(:, 1) = [0.000_dp, 0.000_dp, 0.000_dp]   ! O
      xyz(:, 2) = [0.960_dp, 0.000_dp, 0.000_dp]   ! H, donor
      xyz(:, 3) = [-0.240_dp, 0.930_dp, 0.000_dp]  ! H
      xyz(:, 4) = [2.860_dp, 0.000_dp, 0.000_dp]   ! O, acceptor at 1.90 from H
      xyz(:, 5) = [3.170_dp, 0.900_dp, 0.000_dp]   ! H
      xyz(:, 6) = [3.170_dp, -0.450_dp, 0.780_dp]  ! H

      sys%total_atoms = 6
      sys%n_monomers = 0
      allocate (sys%element_numbers(6))
      sys%element_numbers = [8, 1, 1, 8, 1, 1]
      allocate (sys%coordinates(3, 6))
      sys%coordinates = to_bohr(xyz)
   end subroutine water_dimer

   subroutine carbon_chain(sys)
      !! Four carbons in a line at 1.54 Angstrom
      type(system_geometry_t), intent(out) :: sys

      integer :: iatom

      sys%total_atoms = 4
      sys%n_monomers = 0
      allocate (sys%element_numbers(4))
      sys%element_numbers = 6
      allocate (sys%coordinates(3, 4))
      sys%coordinates = 0.0_dp
      do iatom = 1, 4
         sys%coordinates(1, iatom) = to_bohr(1.54_dp*real(iatom - 1, dp))
      end do
   end subroutine carbon_chain

   subroutine partition_one_atom_each(sys, n)
      type(system_geometry_t), intent(inout) :: sys
      integer, intent(in) :: n
      integer :: imon

      if (allocated(sys%fragment_sizes)) deallocate (sys%fragment_sizes)
      if (allocated(sys%fragment_atoms)) deallocate (sys%fragment_atoms)
      allocate (sys%fragment_sizes(n), sys%fragment_atoms(1, n))
      sys%n_monomers = n
      do imon = 1, n
         sys%fragment_sizes(imon) = 1
         sys%fragment_atoms(1, imon) = imon - 1
      end do
   end subroutine partition_one_atom_each

   subroutine partition_all_together(sys, n)
      type(system_geometry_t), intent(inout) :: sys
      integer, intent(in) :: n
      integer :: iatom

      if (allocated(sys%fragment_sizes)) deallocate (sys%fragment_sizes)
      if (allocated(sys%fragment_atoms)) deallocate (sys%fragment_atoms)
      allocate (sys%fragment_sizes(1), sys%fragment_atoms(n, 1))
      sys%n_monomers = 1
      sys%fragment_sizes(1) = n
      do iatom = 1, n
         sys%fragment_atoms(iatom, 1) = iatom - 1
      end do
   end subroutine partition_all_together

end module test_mqc_bond_perception

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_bond_perception, only: collect_mqc_bond_perception_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_bond_perception", collect_mqc_bond_perception_tests)]
   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
