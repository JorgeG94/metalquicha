module test_mqc_validate
   !! Tests the semantic validator: the checks that stop a calculation which
   !! would converge to a confident wrong number.
   !!
   !! None of these failures crashes -- that is the whole point of the module,
   !! and the reason it is worth testing carefully. A partition that leaves an
   !! atom out, fragment charges that do not sum, a screened term list missing
   !! the subsets its assembly needs: each yields a plausible answer that is
   !! wrong. The suite drives a hand-built `system_geometry_t` and `fraglist_t`
   !! through both modes -- strict, which refuses, and permissive, which only
   !! warns -- and checks each refusal fires exactly when it should.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_validate, only: validate_system, validate_terms
   use mqc_physical_fragment, only: system_geometry_t, to_bohr
   use mqc_fraglist, only: fraglist_t
   use mqc_error, only: error_t
   use pic_types, only: dp, default_int, int64
   implicit none
   private
   public :: collect_mqc_validate_tests

contains

   subroutine collect_mqc_validate_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("valid_system_passes", test_system_ok), &
                  new_unittest("no_atoms_is_refused", test_no_atoms), &
                  new_unittest("bad_multiplicity_is_refused", test_bad_multiplicity), &
                  new_unittest("no_monomers_is_refused", test_no_monomers), &
                  new_unittest("orphan_atom_is_refused", test_orphan_atom), &
                  new_unittest("charge_mismatch_is_refused", test_charge_mismatch), &
                  new_unittest("permissive_mode_only_warns", test_permissive), &
                  new_unittest("bond_audit_flags_undeclared_cut", test_bond_audit), &
                  new_unittest("empty_term_list_is_refused", test_terms_empty), &
                  new_unittest("closed_term_list_passes", test_terms_ok), &
                  new_unittest("out_of_range_monomer_is_refused", test_terms_range), &
                  new_unittest("repeated_monomer_is_refused", test_terms_repeat), &
                  new_unittest("unclosed_term_list_is_refused", test_terms_closure) &
                  ]
   end subroutine collect_mqc_validate_tests

   ! ---- validate_system ----------------------------------------------------

   subroutine test_system_ok(error)
      !! A well-formed two-water cluster passes strict validation
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: geom
      type(error_t) :: err

      geom = two_water_geom()
      call validate_system(geom, strict=.true., error=err)
      call check(error,.not. err%has_error(), &
                 "a valid system should pass: "//err%get_message())

      call geom%destroy()
   end subroutine test_system_ok

   subroutine test_no_atoms(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: geom
      type(error_t) :: err

      geom = two_water_geom()
      geom%total_atoms = 0
      call validate_system(geom, strict=.true., error=err)
      call check(error, err%has_error(), "a system with no atoms should be refused")

      call geom%destroy()
   end subroutine test_no_atoms

   subroutine test_bad_multiplicity(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: geom
      type(error_t) :: err

      geom = two_water_geom()
      geom%multiplicity = 0
      call validate_system(geom, strict=.true., error=err)
      call check(error, err%has_error(), "multiplicity below 1 should be refused")

      call geom%destroy()
   end subroutine test_bad_multiplicity

   subroutine test_no_monomers(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: geom
      type(error_t) :: err

      geom = two_water_geom()
      geom%n_monomers = 0
      call validate_system(geom, strict=.true., error=err)
      call check(error, err%has_error(), "nothing to fragment should be refused")

      call geom%destroy()
   end subroutine test_no_monomers

   subroutine test_orphan_atom(error)
      !! An atom in no monomer is in no fragment, so its electrons vanish
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: geom
      type(error_t) :: err

      geom = two_water_geom()
      ! Drop atom 5 from its monomer: the second monomer now covers only 3,4.
      geom%fragment_sizes(2) = 2
      geom%fragment_atoms(3, 2) = 0
      call validate_system(geom, strict=.true., error=err)
      call check(error, err%has_error(), "an atom owned by no monomer should be refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "no monomer") > 0, &
                 "the message should name the orphaned atom, got: "//err%get_message())

      call geom%destroy()
   end subroutine test_orphan_atom

   subroutine test_charge_mismatch(error)
      !! Fragment charges that do not sum to the molecule's describe other atoms
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: geom
      type(error_t) :: err

      geom = two_water_geom()
      geom%fragment_charges = [1, 0]   ! sums to 1, but the system is neutral
      call validate_system(geom, strict=.true., error=err)
      call check(error, err%has_error(), "a charge that does not sum should be refused")

      call geom%destroy()
   end subroutine test_charge_mismatch

   subroutine test_permissive(error)
      !! The same fault only warns when strict is off
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: geom
      type(error_t) :: err

      geom = two_water_geom()
      geom%fragment_charges = [1, 0]
      call validate_system(geom, strict=.false., error=err)
      call check(error,.not. err%has_error(), &
                 "permissive mode reports rather than refusing")

      call geom%destroy()
   end subroutine test_permissive

   subroutine test_bond_audit(error)
      !! The audit catches a cut the geometry implies but the bond list omits
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: geom
      type(error_t) :: err

      ! Two oxygens 1.2 A apart, one per monomer, with an empty bond list.
      geom%total_atoms = 2
      geom%n_monomers = 2
      geom%atoms_per_monomer = 1
      geom%charge = 0
      geom%multiplicity = 1
      allocate (geom%element_numbers(2), source=[8, 8])
      allocate (geom%coordinates(3, 2))
      geom%coordinates = 0.0_dp
      geom%coordinates(1, 2) = to_bohr(1.2_dp)
      allocate (geom%fragment_sizes(2), source=[1, 1])
      allocate (geom%fragment_atoms(1, 2))
      geom%fragment_atoms(1, 1) = 0
      geom%fragment_atoms(1, 2) = 1
      allocate (geom%fragment_charges(2), source=[0, 0])
      allocate (geom%fragment_multiplicities(2), source=[1, 1])
      allocate (geom%bonds(0))

      ! Without the audit the empty bond list has nothing to check against.
      call validate_system(geom, strict=.true., error=err)
      call check(error,.not. err%has_error(), &
                 "with no audit an empty bond list passes: "//err%get_message())
      if (allocated(error)) then
         call geom%destroy(); return
      end if

      err = error_t()
      call validate_system(geom, strict=.true., error=err, check_bonds=.true.)
      call check(error, err%has_error(), &
                 "the audit should catch the undeclared cut")

      call geom%destroy()
   end subroutine test_bond_audit

   ! ---- validate_terms -----------------------------------------------------

   subroutine test_terms_empty(error)
      !! An empty list has nothing to compute
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: geom
      type(fraglist_t) :: terms
      type(error_t) :: err

      geom = monomer_only_geom(4)
      call validate_terms(terms, geom, strict=.true., error=err)
      call check(error, err%has_error(), "an empty term list should be refused")

      call geom%destroy()
   end subroutine test_terms_empty

   subroutine test_terms_ok(error)
      !! A generated-then-closed list is assemblable and passes
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: geom
      type(fraglist_t) :: terms
      type(error_t) :: err

      geom = monomer_only_geom(4)
      call terms%create(4_default_int, 2_default_int, err)
      call terms%close_subsets(err)   ! pulls in the monomers the dimers need
      call validate_terms(terms, geom, strict=.true., error=err)
      call check(error,.not. err%has_error(), &
                 "a closed list should pass: "//err%get_message())

      call terms%destroy()
      call geom%destroy()
   end subroutine test_terms_ok

   subroutine test_terms_range(error)
      !! A term naming a monomer the system does not have is refused
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: geom
      type(fraglist_t) :: terms
      type(error_t) :: err
      integer(default_int) :: mine(1, 2)

      geom = monomer_only_geom(4)
      mine(1, :) = [9, 1]   ! monomer 9 does not exist in a 4-monomer system
      call terms%replace(mine, 1_int64, 2_default_int, err)
      call validate_terms(terms, geom, strict=.true., error=err)
      call check(error, err%has_error(), "an out-of-range monomer should be refused")

      call terms%destroy()
      call geom%destroy()
   end subroutine test_terms_range

   subroutine test_terms_repeat(error)
      !! A term naming the same monomer twice is refused
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: geom
      type(fraglist_t) :: terms
      type(error_t) :: err
      integer(default_int) :: mine(1, 2)

      geom = monomer_only_geom(4)
      mine(1, :) = [2, 2]
      call terms%replace(mine, 1_int64, 2_default_int, err)
      call validate_terms(terms, geom, strict=.true., error=err)
      call check(error, err%has_error(), "a repeated monomer should be refused")

      call terms%destroy()
      call geom%destroy()
   end subroutine test_terms_repeat

   subroutine test_terms_closure(error)
      !! A lone trimer cannot be assembled without its subsets
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: geom
      type(fraglist_t) :: terms
      type(error_t) :: err
      integer(default_int) :: mine(1, 3)

      geom = monomer_only_geom(4)
      mine(1, :) = [1, 2, 3]
      call terms%replace(mine, 1_int64, 3_default_int, err)
      call validate_terms(terms, geom, strict=.true., error=err)
      call check(error, err%has_error(), "an unclosed list should be refused")
      if (allocated(error)) then
         call terms%destroy(); call geom%destroy(); return
      end if
      call check(error, index(err%get_message(), "close_subsets") > 0, &
                 "the message should point at the fix, got: "//err%get_message())

      call terms%destroy()
      call geom%destroy()
   end subroutine test_terms_closure

   ! ---- fixtures -----------------------------------------------------------

   function two_water_geom() result(geom)
      !! A neutral two-water cluster, each molecule a monomer, atoms 0-based
      type(system_geometry_t) :: geom

      geom%total_atoms = 6
      geom%n_monomers = 2
      geom%atoms_per_monomer = 3
      geom%charge = 0
      geom%multiplicity = 1
      allocate (geom%element_numbers(6), source=[8, 1, 1, 8, 1, 1])
      allocate (geom%coordinates(3, 6))
      geom%coordinates = 0.0_dp
      allocate (geom%fragment_sizes(2), source=[3, 3])
      allocate (geom%fragment_atoms(3, 2))
      geom%fragment_atoms(:, 1) = [0, 1, 2]
      geom%fragment_atoms(:, 2) = [3, 4, 5]
      allocate (geom%fragment_charges(2), source=[0, 0])
      allocate (geom%fragment_multiplicities(2), source=[1, 1])
   end function two_water_geom

   function monomer_only_geom(n_monomers) result(geom)
      !! The bare partition validate_terms reads -- only the monomer count matters
      integer, intent(in) :: n_monomers
      type(system_geometry_t) :: geom

      geom%total_atoms = n_monomers
      geom%n_monomers = n_monomers
      geom%atoms_per_monomer = 1
      geom%charge = 0
      geom%multiplicity = 1
   end function monomer_only_geom

end module test_mqc_validate

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_validate, only: collect_mqc_validate_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_validate", collect_mqc_validate_tests) &
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
