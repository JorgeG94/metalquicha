module test_mqc_xc_spec
   !! Pins the double-hybrid compositions, which are the only functional
   !! definitions this repository owns.
   !!
   !! Everything else is libxc's and needs no test here -- a wrong name fails at
   !! resolution with libxc's own message. Double hybrids are different: libxc
   !! 7.1.2 carries none of them, so their coefficients live in `mqc_xc_spec` and
   !! nothing outside this file would notice a mistyped digit. A wrong weight gives
   !! an energy of entirely the right magnitude.
   !!
   !! The assertion that does the work is the defining relation of a
   !! Grimme-type double hybrid: the semilocal exchange and correlation weights are
   !! the complements of the exact-exchange and perturbative fractions. That is a
   !! property of the construction rather than of any one functional, so it holds
   !! for every row and catches a transposed or mistyped pair without needing a
   !! reference energy.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_xc_spec, only: xc_spec_t, xc_spec_from_name
   implicit none
   private
   public :: collect_mqc_xc_spec_tests

contains

   subroutine collect_mqc_xc_spec_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("double_hybrid_weights_are_complements", test_complements), &
                  new_unittest("plain_functionals_pass_through_to_libxc", test_passthrough), &
                  new_unittest("published_fractions_are_what_we_think", test_published), &
                  new_unittest("mpw2plyp_matches_its_paper", test_mpw2plyp), &
                  new_unittest("an_empty_name_is_refused", test_empty) &
                  ]
   end subroutine collect_mqc_xc_spec_tests

   subroutine test_complements(error)
      !! For every double hybrid: semilocal weights = 1 - the fractions
      type(error_type), allocatable, intent(out) :: error
      type(xc_spec_t) :: spec
      type(error_t) :: err
      character(len=16), parameter :: NAMES(3) = ["b2plyp    ", "b2gp-plyp ", "mpw2plyp  "]
      integer :: i

      do i = 1, size(NAMES)
         call xc_spec_from_name(trim(NAMES(i)), spec, err)
         call check(error,.not. err%has_error(), err%get_message())
         if (allocated(error)) return
         call check(error, spec%is_double_hybrid(), trim(NAMES(i))//" should be a double hybrid")
         if (allocated(error)) return
         call check(error, spec%n_components == 2, trim(NAMES(i))//" should have two components")
         if (allocated(error)) return
         ! Exchange component complements the exact-exchange fraction.
         call check(error, abs(spec%component(1)%weight - (1.0_dp - spec%exx_fraction)) < 1.0e-12_dp, &
                    trim(NAMES(i))//": exchange weight is not 1 - exx_fraction")
         if (allocated(error)) return
         ! Correlation component complements the perturbative fraction.
         call check(error, abs(spec%component(2)%weight - (1.0_dp - spec%pt2_fraction)) < 1.0e-12_dp, &
                    trim(NAMES(i))//": correlation weight is not 1 - pt2_fraction")
         if (allocated(error)) return
         ! Both fractions are genuine fractions. A transposed pair would still
         ! satisfy the complements above, so their ranges are checked too.
         call check(error, spec%exx_fraction > 0.0_dp .and. spec%exx_fraction < 1.0_dp, &
                    trim(NAMES(i))//": exx_fraction is not in (0,1)")
         if (allocated(error)) return
         call check(error, spec%pt2_fraction > 0.0_dp .and. spec%pt2_fraction < 1.0_dp, &
                    trim(NAMES(i))//": pt2_fraction is not in (0,1)")
         if (allocated(error)) return
         ! Every double hybrid in view takes more exact exchange than MP2.
         call check(error, spec%exx_fraction > spec%pt2_fraction, &
                    trim(NAMES(i))//": exx and pt2 fractions look transposed")
         if (allocated(error)) return
      end do
   end subroutine test_complements

   subroutine test_passthrough(error)
      !! A functional libxc owns arrives as one component and claims nothing
      type(error_type), allocatable, intent(out) :: error
      type(xc_spec_t) :: spec
      type(error_t) :: err

      call xc_spec_from_name("PBE", spec, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, spec%from_libxc, "pbe should be left to libxc")
      if (allocated(error)) return
      call check(error, spec%n_components == 1, "pbe should be one component")
      if (allocated(error)) return
      call check(error, spec%name == "pbe", "the name should be lowercased")
      if (allocated(error)) return
      call check(error,.not. spec%is_double_hybrid(), "pbe is not a double hybrid")
      if (allocated(error)) return
      ! Crucially zero: a plain hybrid's fraction is libxc's to report, and a
      ! number here would be a second source of truth that could disagree.
      call check(error,.not. spec%needs_exact_exchange(), &
                 "a libxc-owned functional must not state its own exx fraction")
      if (allocated(error)) return
      call xc_spec_from_name("b3lyp", spec, err)
      call check(error, spec%from_libxc .and. .not. spec%needs_exact_exchange(), &
                 "b3lyp is a hybrid but libxc owns its fraction")
   end subroutine test_passthrough

   subroutine test_published(error)
      !! The fractions, against the papers, one row at a time
      !!
      !! Spelled out separately from the structural test because a composition can
      !! be internally consistent and still be the wrong functional.
      type(error_type), allocatable, intent(out) :: error
      type(xc_spec_t) :: spec
      type(error_t) :: err

      ! Cross-checked against an independent implementation as well as the paper:
      ! pyscf-forge's `pyscf.dh.DFDH` reports cx and c_os/c_ss for each of these,
      ! and all three rows agree exactly with the values below. That matters more
      ! than the citation does -- a transcription error survives a citation and does
      ! not survive another code reporting the same numbers.
      !
      ! Grimme, JCP 124, 034108 (2006); DFDH cx=0.53, c_os=c_ss=0.27
      call xc_spec_from_name("b2plyp", spec, err)
      call check(error, abs(spec%exx_fraction - 0.53_dp) < 1.0e-12_dp, "B2PLYP a_x is 0.53")
      if (allocated(error)) return
      call check(error, abs(spec%pt2_fraction - 0.27_dp) < 1.0e-12_dp, "B2PLYP a_c is 0.27")
      if (allocated(error)) return
      call check(error, spec%component(1)%name == "gga_x_b88", "B2PLYP exchange is B88")
      if (allocated(error)) return
      call check(error, spec%component(2)%name == "gga_c_lyp", "B2PLYP correlation is LYP")
      if (allocated(error)) return

      ! Karton et al., JPCA 112, 12868 (2008); DFDH cx=0.65, c_os=c_ss=0.36
      call xc_spec_from_name("b2gp-plyp", spec, err)
      call check(error, abs(spec%exx_fraction - 0.65_dp) < 1.0e-12_dp, "B2GP-PLYP a_x is 0.65")
      if (allocated(error)) return
      call check(error, abs(spec%pt2_fraction - 0.36_dp) < 1.0e-12_dp, "B2GP-PLYP a_c is 0.36")
   end subroutine test_published

   subroutine test_mpw2plyp(error)
      !! mPW2-PLYP, the third row, against DFDH's reported composition
      type(error_type), allocatable, intent(out) :: error
      type(xc_spec_t) :: spec
      type(error_t) :: err

      ! Schwabe and Grimme, PCCP 8, 4398 (2006); DFDH cx=0.55, c_os=c_ss=0.25
      call xc_spec_from_name("mpw2plyp", spec, err)
      call check(error, abs(spec%exx_fraction - 0.55_dp) < 1.0e-12_dp, "mPW2-PLYP a_x is 0.55")
      if (allocated(error)) return
      call check(error, abs(spec%pt2_fraction - 0.25_dp) < 1.0e-12_dp, "mPW2-PLYP a_c is 0.25")
      if (allocated(error)) return
      call check(error, spec%component(1)%name == "gga_x_mpw91", "mPW2-PLYP exchange is mPW91")
   end subroutine test_mpw2plyp

   subroutine test_empty(error)
      !! No functional named is an error, not a default
      type(error_type), allocatable, intent(out) :: error
      type(xc_spec_t) :: spec
      type(error_t) :: err

      call xc_spec_from_name("   ", spec, err)
      call check(error, err%has_error(), "an empty functional name must be refused")
   end subroutine test_empty

end module test_mqc_xc_spec

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_xc_spec, only: collect_mqc_xc_spec_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'
   stat = 0
   testsuites = [new_testsuite("mqc_xc_spec", collect_mqc_xc_spec_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
