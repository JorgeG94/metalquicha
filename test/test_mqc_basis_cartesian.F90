#ifdef MQC_HAVE_BASIS_SETS
module test_mqc_basis_cartesian
   !! Pins which angular form the reader takes each basis set to be in.
   !!
   !! Basis Set Exchange records a `function_type` on every shell and the
   !! reader used to ignore it, so every basis was built spherical. 6-31G*
   !! specifies Cartesian d functions and BSE marks them `gto_cartesian`: a
   !! Cartesian d shell carries six functions where a spherical one carries
   !! five, so a deck asking for 6-31G* was computed in a smaller basis of the
   !! same name, 1.4 mHartree out on water with nothing said about it.
   !!
   !! What these check is the decision, not the integrals -- which convention
   !! the reader concludes, and which combinations it refuses. The routing that
   !! decision drives is checked in `test_mqc_libcint_cartesian`, where it
   !! shows up as a basis function count.
   !!
   !! The distinction that matters here: s and p shells carry the same number
   !! of functions either way, so a shell at or below p is evidence of nothing.
   !! Hydrogen in 6-31G* has only s shells and so does not disagree with the
   !! oxygen next to it -- it simply has no opinion.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_json_basis_reader, only: read_json_basis_element, build_molecular_basis_json
   use mqc_cgto, only: atomic_basis_type, molecular_basis_type, &
                       ANGULAR_FORM_UNSET, ANGULAR_FORM_SPHERICAL, ANGULAR_FORM_CARTESIAN
   use mqc_error, only: error_t
   implicit none
   private
   public :: collect_mqc_basis_cartesian_tests

   character(len=*), parameter :: POPLE_STAR = "../basis_sets/6-31g_st_.json"
   character(len=*), parameter :: POPLE_PLAIN = "../basis_sets/6-31g.json"
   character(len=*), parameter :: CC_PVDZ = "../basis_sets/cc-pvdz.json"
   character(len=*), parameter :: DEF2_SVP = "../basis_sets/def2-svp.json"
   character(len=*), parameter :: JKFIT = "../basis_sets/def2-universal-jkfit.json"

contains

   subroutine collect_mqc_basis_cartesian_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("pople_star_oxygen_is_cartesian", test_oxygen_cartesian), &
                  new_unittest("pople_star_hydrogen_has_no_opinion", test_hydrogen_unset), &
                  new_unittest("correlation_consistent_is_spherical", test_ccpvdz_spherical), &
                  new_unittest("karlsruhe_is_spherical", test_def2_spherical), &
                  new_unittest("jkfit_is_spherical", test_jkfit_spherical), &
                  new_unittest("water_in_pople_star_is_cartesian", test_water_molecule), &
                  new_unittest("water_in_plain_pople_is_not", test_water_plain_pople), &
                  new_unittest("water_in_ccpvdz_is_spherical", test_water_ccpvdz), &
                  new_unittest("mixed_element_is_refused", test_mixed_element_refused), &
                  new_unittest("mixed_molecule_is_refused", test_mixed_molecule_refused) &
                  ]
   end subroutine collect_mqc_basis_cartesian_tests

   subroutine read_element(path, symbol, form, ok, error)
      !! Parse one element, failing the test if the file will not load
      character(len=*), intent(in) :: path, symbol
      integer, intent(out) :: form
      logical, intent(out) :: ok
      type(error_type), allocatable, intent(out) :: error

      type(atomic_basis_type) :: basis
      type(error_t) :: parse_error

      form = ANGULAR_FORM_UNSET
      call read_json_basis_element(path, symbol, basis, parse_error)
      ok = .not. parse_error%has_error()
      call check(error, ok, trim(path)//" must parse (run cmake to extract it): "// &
                 trim(parse_error%get_message()))
      if (.not. ok) return
      form = basis%angular_form
      call basis%destroy()
   end subroutine read_element

   subroutine test_oxygen_cartesian(error)
      !! 6-31G* oxygen: the d shell is `gto_cartesian`, so oxygen is 6d
      type(error_type), allocatable, intent(out) :: error
      integer :: form
      logical :: ok

      call read_element(POPLE_STAR, "O", form, ok, error)
      if (allocated(error) .or. .not. ok) return

      call check(error, form, ANGULAR_FORM_CARTESIAN, &
                 "6-31G* oxygen carries a Cartesian d shell")
   end subroutine test_oxygen_cartesian

   subroutine test_hydrogen_unset(error)
      !! 6-31G* hydrogen has only s shells, which look the same either way
      !!
      !! This is the case that makes the difference between a per-element flag
      !! that works and one that does not. Calling hydrogen "spherical" here
      !! would put it in conflict with the oxygen it is bonded to, and water in
      !! 6-31G* would be refused as a mixed basis rather than computed.
      type(error_type), allocatable, intent(out) :: error
      integer :: form
      logical :: ok

      call read_element(POPLE_STAR, "H", form, ok, error)
      if (allocated(error) .or. .not. ok) return

      call check(error, form, ANGULAR_FORM_UNSET, &
                 "6-31G* hydrogen has no shell above s and so votes for neither form")
   end subroutine test_hydrogen_unset

   subroutine test_ccpvdz_spherical(error)
      !! cc-pVDZ oxygen: `gto_spherical` d, which is 5d
      type(error_type), allocatable, intent(out) :: error
      integer :: form
      logical :: ok

      call read_element(CC_PVDZ, "O", form, ok, error)
      if (allocated(error) .or. .not. ok) return

      call check(error, form, ANGULAR_FORM_SPHERICAL, "cc-pVDZ oxygen is spherical")
   end subroutine test_ccpvdz_spherical

   subroutine test_def2_spherical(error)
      !! def2-SVP oxygen is spherical -- the set the CPU validation suite runs
      type(error_type), allocatable, intent(out) :: error
      integer :: form
      logical :: ok

      call read_element(DEF2_SVP, "O", form, ok, error)
      if (allocated(error) .or. .not. ok) return

      call check(error, form, ANGULAR_FORM_SPHERICAL, "def2-SVP oxygen is spherical")
   end subroutine test_def2_spherical

   subroutine test_jkfit_spherical(error)
      !! Every fitting set here is spherical, which is why a Cartesian orbital
      !! basis cannot be fitted with one of them
      type(error_type), allocatable, intent(out) :: error
      integer :: form
      logical :: ok

      call read_element(JKFIT, "O", form, ok, error)
      if (allocated(error) .or. .not. ok) return

      call check(error, form, ANGULAR_FORM_SPHERICAL, &
                 "def2-universal-JKFIT oxygen is spherical")
   end subroutine test_jkfit_spherical

   subroutine test_water_molecule(error)
      !! Water in 6-31G*: oxygen decides, hydrogen abstains, the basis is 6d
      type(error_type), allocatable, intent(out) :: error
      type(molecular_basis_type) :: basis
      type(error_t) :: parse_error
      character(len=2) :: symbols(3)

      symbols = ["O ", "H ", "H "]
      call build_molecular_basis_json(POPLE_STAR, symbols, basis, parse_error)
      call check(error,.not. parse_error%has_error(), &
                 "water must build in 6-31G*: "//trim(parse_error%get_message()))
      if (allocated(error)) return

      call check(error, basis%angular_form, ANGULAR_FORM_CARTESIAN, &
                 "water in 6-31G* is a Cartesian basis")
      if (allocated(error)) return
      call check(error, basis%is_cartesian(), "is_cartesian() must agree")
      call basis%destroy()
   end subroutine test_water_molecule

   subroutine test_water_plain_pople(error)
      !! Plain 6-31G on water has no shell above p, so nothing is decided
      !!
      !! The file does mark Cartesian d shells -- for potassium upwards. None
      !! of them belongs to an atom in this molecule, and the convention has to
      !! be read from the elements actually present or every deck using a set
      !! with a transition metal in the file would turn Cartesian.
      type(error_type), allocatable, intent(out) :: error
      type(molecular_basis_type) :: basis
      type(error_t) :: parse_error
      character(len=2) :: symbols(3)

      symbols = ["O ", "H ", "H "]
      call build_molecular_basis_json(POPLE_PLAIN, symbols, basis, parse_error)
      call check(error,.not. parse_error%has_error(), &
                 "water must build in 6-31G: "//trim(parse_error%get_message()))
      if (allocated(error)) return

      call check(error, basis%angular_form, ANGULAR_FORM_UNSET, &
                 "water in 6-31G reaches only p, so neither form is implied")
      if (allocated(error)) return
      call check(error,.not. basis%is_cartesian(), &
                 "an undecided basis takes the spherical entry points")
      call basis%destroy()
   end subroutine test_water_plain_pople

   subroutine test_water_ccpvdz(error)
      !! The spherical case has to stay spherical; this is what did not change
      type(error_type), allocatable, intent(out) :: error
      type(molecular_basis_type) :: basis
      type(error_t) :: parse_error
      character(len=2) :: symbols(3)

      symbols = ["O ", "H ", "H "]
      call build_molecular_basis_json(CC_PVDZ, symbols, basis, parse_error)
      call check(error,.not. parse_error%has_error(), &
                 "water must build in cc-pVDZ: "//trim(parse_error%get_message()))
      if (allocated(error)) return

      call check(error,.not. basis%is_cartesian(), "water in cc-pVDZ is spherical")
      call basis%destroy()
   end subroutine test_water_ccpvdz

   subroutine test_mixed_element_refused(error)
      !! 6-31G* scandium marks its d shell Cartesian and its f shell spherical
      !!
      !! libcint picks the form per *call*, so there is no way to honour both
      !! inside one molecule. Choosing one silently is the bug this was all
      !! fixed for, only quieter, so the element is refused by name.
      type(error_type), allocatable, intent(out) :: error
      type(atomic_basis_type) :: basis
      type(error_t) :: parse_error

      call read_json_basis_element(POPLE_STAR, "Sc", basis, parse_error)

      call check(error, parse_error%has_error(), &
                 "6-31G* scandium mixes a Cartesian d with a spherical f and must be refused")
      if (allocated(error)) return
      call check(error, index(parse_error%get_message(), "Cartesian") > 0, &
                 "the refusal has to say what it refused: "//trim(parse_error%get_message()))
      call basis%destroy()
   end subroutine test_mixed_element_refused

   subroutine test_mixed_molecule_refused(error)
      !! Two elements of one molecule disagreeing is refused as well
      !!
      !! Oxygen in 6-31G* is Cartesian; the same file gives argon a Cartesian d
      !! and scandium a spherical f, so a molecule reaching both is a basis
      !! that cannot be built in one form. The refusal comes from the element
      !! reader here, which is the first place the disagreement is visible.
      type(error_type), allocatable, intent(out) :: error
      type(molecular_basis_type) :: basis
      type(error_t) :: parse_error
      character(len=2) :: symbols(2)

      symbols = ["O ", "Sc"]
      call build_molecular_basis_json(POPLE_STAR, symbols, basis, parse_error)

      call check(error, parse_error%has_error(), &
                 "a molecule whose elements disagree about the angular form must be refused")
      call basis%destroy()
   end subroutine test_mixed_molecule_refused

end module test_mqc_basis_cartesian

program tester_mqc_basis_cartesian
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_basis_cartesian, only: collect_mqc_basis_cartesian_tests
   implicit none

   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [ &
                new_testsuite("mqc_basis_cartesian", collect_mqc_basis_cartesian_tests) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_basis_cartesian

#else
program tester_mqc_basis_cartesian
   !! The basis files this reads are unpacked from the Basis Set Exchange
   !! bundle by cmake/modules/ExtractBasisSets.cmake, which only the CMake
   !! build runs. They are gitignored, so under fpm there is nothing on disk
   !! and the test would fail on a missing file rather than on the reader.
   implicit none
   write (*, "(A)") "# mqc_basis_cartesian: skipped (basis sets not extracted)"
end program tester_mqc_basis_cartesian
#endif
