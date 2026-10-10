!! The literature references libxc records for a functional
module test_mqc_xc_references
   !! `xc_functional_references` reads the papers to cite out of libxc's own
   !! reference table rather than out of a list kept here. These checks pin the
   !! three things that can go wrong reading it: walking off the end of a
   !! component's list (libxc's `number` comes back -1 after the last), losing
   !! a component of a functional built from several, and listing a paper twice
   !! when two components share it -- PBE exchange and PBE correlation both cite
   !! Perdew, Burke and Ernzerhof (1996).
   !!
   !! Every reference is printed, so a run shows what a user will be told to
   !! cite. A build without libxc has nothing to check.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available, xc_reference_t, &
                         xc_functional_references
   implicit none
   private

   public :: collect_mqc_xc_references_tests

   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.190459_dp, &
                           0.0_dp, 1.459853_dp, -0.884001_dp, &
                           0.0_dp, -1.459853_dp, -0.884001_dp], [3, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]

contains

   subroutine collect_mqc_xc_references_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("pbe_cites_its_paper_once", test_pbe), &
                  new_unittest("b3lyp_cites_stephens", test_b3lyp), &
                  new_unittest("every_reference_has_a_citation", test_well_formed) &
                  ]
   end subroutine collect_mqc_xc_references_tests

   subroutine references_of(functional, refs, error)
      !! The references for `functional` on water, printed
      character(len=*), intent(in) :: functional
      type(xc_reference_t), allocatable, intent(out) :: refs(:)
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(xc_context_t) :: xc
      type(error_t) :: err
      integer :: k

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return
      call xc_context_create(mol, functional, xc, err, level=1)
      call check(error,.not. err%has_error(), "the context for "//functional//" should build")
      if (allocated(error)) return

      call xc_functional_references(xc, refs)
      print "(a, a, a, i0, a)", "    ", functional, ": ", size(refs), " reference(s)"
      do k = 1, size(refs)
         print "(a, a, a, a, a, a)", "      [", refs(k)%functional, "] ", refs(k)%citation, &
            "  doi:", refs(k)%doi
      end do
   end subroutine references_of

   integer function count_containing(refs, text) result(n)
      !! How many citations contain `text`
      type(xc_reference_t), intent(in) :: refs(:)
      character(len=*), intent(in) :: text
      integer :: k
      n = 0
      do k = 1, size(refs)
         if (index(refs(k)%citation, text) > 0) n = n + 1
      end do
   end function count_containing

   subroutine test_pbe(error)
      type(error_type), allocatable, intent(out) :: error
      type(xc_reference_t), allocatable :: refs(:)

      if (.not. xc_available()) return
      call references_of("pbe", refs, error)
      if (allocated(error)) return
      call check(error, count_containing(refs, "3865") == 1, &
                 "Perdew, Burke and Ernzerhof (1996) is listed exactly once, though both "// &
                 "components cite it")
   end subroutine test_pbe

   subroutine test_b3lyp(error)
      type(error_type), allocatable, intent(out) :: error
      type(xc_reference_t), allocatable :: refs(:)

      if (.not. xc_available()) return
      call references_of("b3lyp", refs, error)
      if (allocated(error)) return
      call check(error, count_containing(refs, "Stephens") >= 1, &
                 "B3LYP's references include Stephens et al. (1994)")
   end subroutine test_b3lyp

   subroutine test_well_formed(error)
      type(error_type), allocatable, intent(out) :: error
      type(xc_reference_t), allocatable :: refs(:)
      integer :: k

      if (.not. xc_available()) return
      ! A double hybrid, so the context holds more than one libxc component.
      call references_of("b2plyp", refs, error)
      if (allocated(error)) return
      call check(error, size(refs) >= 1, "a functional has at least one reference")
      if (allocated(error)) return
      do k = 1, size(refs)
         call check(error, len(refs(k)%citation) > 0 .and. len(refs(k)%functional) > 0, &
                    "every reference names its component and has a citation")
         if (allocated(error)) return
      end do
   end subroutine test_well_formed

end module test_mqc_xc_references

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_xc_references, only: collect_mqc_xc_references_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_xc_references", collect_mqc_xc_references_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
