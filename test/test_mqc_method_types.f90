module test_mqc_method_types
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_method_types, only: method_type_from_string, method_type_to_string, &
                               METHOD_TYPE_GFN1, METHOD_TYPE_GFN2, METHOD_TYPE_HF, METHOD_TYPE_DFT, &
                               METHOD_TYPE_MCSCF, METHOD_TYPE_MP2, METHOD_TYPE_MP2_F12, METHOD_TYPE_CCSD, &
                               METHOD_TYPE_CCSD_T, METHOD_TYPE_CCSD_F12, METHOD_TYPE_CCSD_T_F12, METHOD_TYPE_UNKNOWN
   use pic_types, only: int32
   implicit none
   private
   public :: collect_mqc_method_types_tests

contains

   subroutine collect_mqc_method_types_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("gfn1_from_string", test_gfn1_from_string), &
                  new_unittest("gfn2_from_string", test_gfn2_from_string), &
                  new_unittest("hf_from_string", test_hf_from_string), &
                  new_unittest("dft_from_string", test_dft_from_string), &
                  new_unittest("mcscf_from_string", test_mcscf_from_string), &
                  new_unittest("mp2_from_string", test_mp2_from_string), &
                  new_unittest("ccsd_from_string", test_ccsd_from_string), &
                  new_unittest("unknown_from_string", test_unknown_from_string), &
                  new_unittest("case_insensitive", test_case_insensitive), &
                  new_unittest("to_string_roundtrip", test_to_string_roundtrip) &
                  ]
   end subroutine collect_mqc_method_types_tests

   subroutine test_gfn1_from_string(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, method_type_from_string("gfn1") == METHOD_TYPE_GFN1, &
                 "gfn1 should map to METHOD_TYPE_GFN1")
      if (allocated(error)) return

      call check(error, method_type_from_string("gfn1-xtb") == METHOD_TYPE_GFN1, &
                 "gfn1-xtb should map to METHOD_TYPE_GFN1")
   end subroutine test_gfn1_from_string

   subroutine test_gfn2_from_string(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, method_type_from_string("gfn2") == METHOD_TYPE_GFN2, &
                 "gfn2 should map to METHOD_TYPE_GFN2")
      if (allocated(error)) return

      call check(error, method_type_from_string("gfn2-xtb") == METHOD_TYPE_GFN2, &
                 "gfn2-xtb should map to METHOD_TYPE_GFN2")
   end subroutine test_gfn2_from_string

   subroutine test_hf_from_string(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, method_type_from_string("hf") == METHOD_TYPE_HF, &
                 "hf should map to METHOD_TYPE_HF")
      if (allocated(error)) return

      call check(error, method_type_from_string("rhf") == METHOD_TYPE_HF, &
                 "rhf should map to METHOD_TYPE_HF")
      if (allocated(error)) return

      call check(error, method_type_from_string("uhf") == METHOD_TYPE_HF, &
                 "uhf should map to METHOD_TYPE_HF")
      if (allocated(error)) return

      call check(error, method_type_from_string("hartree-fock") == METHOD_TYPE_HF, &
                 "hartree-fock should map to METHOD_TYPE_HF")
   end subroutine test_hf_from_string

   subroutine test_dft_from_string(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, method_type_from_string("dft") == METHOD_TYPE_DFT, &
                 "dft should map to METHOD_TYPE_DFT")
      if (allocated(error)) return

      call check(error, method_type_from_string("ks") == METHOD_TYPE_DFT, &
                 "ks should map to METHOD_TYPE_DFT")
      if (allocated(error)) return

      call check(error, method_type_from_string("kohn-sham") == METHOD_TYPE_DFT, &
                 "kohn-sham should map to METHOD_TYPE_DFT")
   end subroutine test_dft_from_string

   subroutine test_mcscf_from_string(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, method_type_from_string("mcscf") == METHOD_TYPE_MCSCF, &
                 "mcscf should map to METHOD_TYPE_MCSCF")
      if (allocated(error)) return

      call check(error, method_type_from_string("casscf") == METHOD_TYPE_MCSCF, &
                 "casscf should map to METHOD_TYPE_MCSCF")
      if (allocated(error)) return

      call check(error, method_type_from_string("casci") == METHOD_TYPE_MCSCF, &
                 "casci should map to METHOD_TYPE_MCSCF")
   end subroutine test_mcscf_from_string

   subroutine test_mp2_from_string(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, method_type_from_string("mp2") == METHOD_TYPE_MP2, &
                 "mp2 should map to METHOD_TYPE_MP2")
      if (allocated(error)) return

      call check(error, method_type_from_string("ri-mp2") == METHOD_TYPE_MP2, &
                 "ri-mp2 should map to METHOD_TYPE_MP2")
      if (allocated(error)) return

      call check(error, method_type_from_string("mp2-f12") == METHOD_TYPE_MP2_F12, &
                 "mp2-f12 should map to METHOD_TYPE_MP2_F12")
   end subroutine test_mp2_from_string

   subroutine test_ccsd_from_string(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, method_type_from_string("ccsd") == METHOD_TYPE_CCSD, &
                 "ccsd should map to METHOD_TYPE_CCSD")
      if (allocated(error)) return

      call check(error, method_type_from_string("ccsd(t)") == METHOD_TYPE_CCSD_T, &
                 "ccsd(t) should map to METHOD_TYPE_CCSD_T")
      if (allocated(error)) return

      call check(error, method_type_from_string("ccsd-f12") == METHOD_TYPE_CCSD_F12, &
                 "ccsd-f12 should map to METHOD_TYPE_CCSD_F12")
      if (allocated(error)) return

      call check(error, method_type_from_string("ccsd(t)-f12") == METHOD_TYPE_CCSD_T_F12, &
                 "ccsd(t)-f12 should map to METHOD_TYPE_CCSD_T_F12")
   end subroutine test_ccsd_from_string

   subroutine test_unknown_from_string(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, method_type_from_string("garbage") == METHOD_TYPE_UNKNOWN, &
                 "unknown string should map to METHOD_TYPE_UNKNOWN")
      if (allocated(error)) return

      call check(error, method_type_from_string("") == METHOD_TYPE_UNKNOWN, &
                 "empty string should map to METHOD_TYPE_UNKNOWN")
   end subroutine test_unknown_from_string

   subroutine test_case_insensitive(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, method_type_from_string("GFN2") == METHOD_TYPE_GFN2, &
                 "GFN2 (uppercase) should map to METHOD_TYPE_GFN2")
      if (allocated(error)) return

      call check(error, method_type_from_string("Gfn2") == METHOD_TYPE_GFN2, &
                 "Gfn2 (mixed case) should map to METHOD_TYPE_GFN2")
      if (allocated(error)) return

      call check(error, method_type_from_string("CCSD(T)") == METHOD_TYPE_CCSD_T, &
                 "CCSD(T) (uppercase) should map to METHOD_TYPE_CCSD_T")
   end subroutine test_case_insensitive

   subroutine test_to_string_roundtrip(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, method_type_to_string(METHOD_TYPE_GFN1) == "gfn1", &
                 "METHOD_TYPE_GFN1 should convert to 'gfn1'")
      if (allocated(error)) return

      call check(error, method_type_to_string(METHOD_TYPE_GFN2) == "gfn2", &
                 "METHOD_TYPE_GFN2 should convert to 'gfn2'")
      if (allocated(error)) return

      call check(error, method_type_to_string(METHOD_TYPE_HF) == "hf", &
                 "METHOD_TYPE_HF should convert to 'hf'")
      if (allocated(error)) return

      call check(error, method_type_to_string(METHOD_TYPE_DFT) == "dft", &
                 "METHOD_TYPE_DFT should convert to 'dft'")
      if (allocated(error)) return

      call check(error, method_type_to_string(METHOD_TYPE_MCSCF) == "mcscf", &
                 "METHOD_TYPE_MCSCF should convert to 'mcscf'")
      if (allocated(error)) return

      call check(error, method_type_to_string(METHOD_TYPE_CCSD_T) == "ccsd(t)", &
                 "METHOD_TYPE_CCSD_T should convert to 'ccsd(t)'")
      if (allocated(error)) return

      call check(error, method_type_to_string(METHOD_TYPE_UNKNOWN) == "unknown", &
                 "METHOD_TYPE_UNKNOWN should convert to 'unknown'")
   end subroutine test_to_string_roundtrip

end module test_mqc_method_types

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_method_types, only: collect_mqc_method_types_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_method_types", collect_mqc_method_types_tests) &
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
