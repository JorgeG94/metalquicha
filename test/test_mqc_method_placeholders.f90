!! What the method objects promise, now that most of them are not placeholders
!!
!! Hartree-Fock and Kohn-Sham are real calculations on this path: analytic
!! gradients and finite-difference Hessians both exist for them. Only MCSCF is
!! still a stub, and its tests below are the only ones that assert stub
!! behaviour. The file keeps its name because the target does.
!!
!! **The contract every test here shares** is not that a particular number comes
!! back. It is that exactly one of "produced a result" and "reported an error"
!! is true. Which one depends on the build and the machine -- these fragments
!! carry no basis set, and the test environment has no `basis_sets/` in its
!! working directory -- so asserting success would make the suite a test of the
!! environment. Asserting the exclusive-or tests the thing that must hold
!! everywhere: nothing is ever invented, and nothing that succeeded is ever
!! reported as failed.
module test_mqc_method_placeholders
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_method_hf, only: hf_method_t, hf_options_t
   use mqc_method_dft, only: dft_method_t, dft_options_t
   use mqc_method_mcscf, only: mcscf_method_t, mcscf_options_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_elements, only: element_symbol_to_number
   use pic_types, only: dp
   use pic_test_helpers, only: is_equal
   implicit none
   private
   public :: collect_mqc_method_placeholders_tests

contains

   subroutine collect_mqc_method_placeholders_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("hf_energy", test_hf_energy), &
                  new_unittest("hf_gradient", test_hf_gradient), &
                  new_unittest("hf_hessian", test_hf_hessian), &
                  new_unittest("dft_energy", test_dft_energy), &
                  new_unittest("dft_gradient", test_dft_gradient), &
                  new_unittest("dft_hessian", test_dft_hessian), &
                  new_unittest("mcscf_energy_no_active_space", test_mcscf_no_active), &
                  new_unittest("mcscf_energy_with_active_space", test_mcscf_with_active), &
                  new_unittest("mcscf_gradient", test_mcscf_gradient) &
                  ]
   end subroutine collect_mqc_method_placeholders_tests

   subroutine create_test_fragment(fragment)
      type(physical_fragment_t), intent(out) :: fragment

      fragment%n_atoms = 3
      allocate (fragment%element_numbers(3))
      allocate (fragment%coordinates(3, 3))

      fragment%element_numbers(1) = element_symbol_to_number("O")
      fragment%element_numbers(2) = element_symbol_to_number("H")
      fragment%element_numbers(3) = element_symbol_to_number("H")

      fragment%coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.22_dp]
      fragment%coordinates(:, 2) = [0.0_dp, 1.43_dp, -0.88_dp]
      fragment%coordinates(:, 3) = [0.0_dp, -1.43_dp, -0.88_dp]

      fragment%charge = 0
      fragment%multiplicity = 1
      fragment%nelec = 10
   end subroutine create_test_fragment

   subroutine test_hf_energy(error)
      !! HF is no longer a placeholder: it either produces a real energy or
      !! reports an error, and must never invent a number.
      !!
      !! What it does here depends on the build and the machine. Without the
      !! cuEST backend it reports "no integral backend"; with it, this test
      !! environment has no basis_sets/ directory in the working directory and
      !! no supported GPU, so it reports that instead. Either way the contract
      !! under test is the same: no silent fake energy.
      type(error_type), allocatable, intent(out) :: error
      type(hf_method_t) :: method
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call create_test_fragment(fragment)
      call method%calc_energy(fragment, result)

      call check(error, result%has_error .neqv. result%has_energy, &
                 "HF must either produce an energy or report an error, not both or neither")
      if (allocated(error)) return

      if (result%has_error) then
         call check(error, len_trim(result%error%get_message()) > 0, &
                    "A failed HF calculation must carry a diagnostic message")
      end if

      call fragment%destroy()
   end subroutine test_hf_energy

   subroutine test_hf_gradient(error)
      type(error_type), allocatable, intent(out) :: error
      type(hf_method_t) :: method
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call create_test_fragment(fragment)
      call method%calc_gradient(fragment, result)

      ! The analytic HF gradient is implemented. What this asserts is the
      ! contract around it rather than its value, which
      ! `validation/check_scf_gradient` covers against PySCF and finite
      ! differences: one of a gradient and an error, never both, never neither,
      ! and never a zero gradient passed off as real.
      call check(error, result%has_error .neqv. result%has_gradient, &
                 "HF must either produce a gradient or report an error, not both "// &
                 "or neither")
      if (allocated(error)) return

      if (result%has_error) then
         call check(error, len_trim(result%error%get_message()) > 0, &
                    "A failed HF gradient must carry a diagnostic message")
      else
         call check_gradient_shape(error, result, fragment%n_atoms, "HF")
      end if

      call fragment%destroy()
   end subroutine test_hf_gradient

   subroutine test_hf_hessian(error)
      type(error_type), allocatable, intent(out) :: error
      type(hf_method_t) :: method
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call create_test_fragment(fragment)
      call method%calc_hessian(fragment, result)

      ! Built by finite differences of the analytic gradient, so it exists
      ! wherever that does.
      call check(error, result%has_error .neqv. result%has_hessian, &
                 "HF must either produce a Hessian or report an error, not both "// &
                 "or neither")
      if (allocated(error)) return

      if (result%has_error) then
         call check(error, len_trim(result%error%get_message()) > 0, &
                    "A failed HF Hessian must carry a diagnostic message")
      else
         call check_hessian_shape(error, result, fragment%n_atoms, "HF")
      end if

      call fragment%destroy()
   end subroutine test_hf_hessian

   subroutine test_dft_energy(error)
      type(error_type), allocatable, intent(out) :: error
      type(dft_method_t) :: method
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call create_test_fragment(fragment)

      method%options%verbose = .true.
      method%options%use_density_fitting = .true.

      call method%calc_energy(fragment, result)

      ! DFT is a real calculation now: it either produces an energy or reports
      ! why it could not. Which of the two happens depends on the build and on
      ! whether a supported GPU and a basis_sets/ directory are present, but
      ! inventing a number is never acceptable.
      call check(error, result%has_error .neqv. result%has_energy, &
                 "DFT must either produce an energy or report an error, not both or neither")
      if (allocated(error)) return

      if (result%has_error) then
         call check(error, len_trim(result%error%get_message()) > 0, &
                    "A failed DFT calculation must carry a diagnostic message")
         if (allocated(error)) return
      end if

      ! Dispersion is not implemented for the cuEST backend, so asking for it
      ! must fail loudly rather than quietly return an undispersed energy.
      method%options%use_dispersion = .true.
      call method%calc_energy(fragment, result)
      call check(error, result%has_error, &
                 "Requesting unimplemented dispersion must report an error")

      call fragment%destroy()
   end subroutine test_dft_energy

   subroutine test_dft_gradient(error)
      type(error_type), allocatable, intent(out) :: error
      type(dft_method_t) :: method
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call create_test_fragment(fragment)
      method%options%verbose = .true.

      call method%calc_gradient(fragment, result)

      ! Implemented through hybrid GGA. Meta-GGAs and range-separated hybrids
      ! are refused rather than approximated, which is an error and so is
      ! covered by the same exclusive-or.
      call check(error, result%has_error .neqv. result%has_gradient, &
                 "DFT must either produce a gradient or report an error, not both "// &
                 "or neither")
      if (allocated(error)) return

      if (result%has_error) then
         call check(error, len_trim(result%error%get_message()) > 0, &
                    "A failed DFT gradient must carry a diagnostic message")
      else
         call check_gradient_shape(error, result, fragment%n_atoms, "DFT")
      end if

      call fragment%destroy()
   end subroutine test_dft_gradient

   subroutine test_dft_hessian(error)
      type(error_type), allocatable, intent(out) :: error
      type(dft_method_t) :: method
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call create_test_fragment(fragment)
      method%options%verbose = .true.

      call method%calc_hessian(fragment, result)

      call check(error, result%has_error .neqv. result%has_hessian, &
                 "DFT must either produce a Hessian or report an error, not both "// &
                 "or neither")
      if (allocated(error)) return

      if (result%has_error) then
         call check(error, len_trim(result%error%get_message()) > 0, &
                    "A failed DFT Hessian must carry a diagnostic message")
      else
         call check_hessian_shape(error, result, fragment%n_atoms, "DFT")
      end if

      call fragment%destroy()
   end subroutine test_dft_hessian

   subroutine check_gradient_shape(error, result, n_atoms, label)
      !! A gradient that claims to exist has to be one
      !!
      !! Shape and finiteness rather than value. An unallocated array behind a
      !! true `has_gradient`, or a NaN in it, is the failure this catches -- and
      !! it is exactly what a zero-initialised placeholder would *not* have been
      !! caught doing by the old assertion, which only asked whether the flag
      !! was false.
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), intent(in) :: result
      integer, intent(in) :: n_atoms
      character(len=*), intent(in) :: label

      call check(error, allocated(result%gradient), &
                 label//" claims a gradient but did not allocate one")
      if (allocated(error)) return

      call check(error, size(result%gradient, 1) == 3 .and. &
                 size(result%gradient, 2) == n_atoms, &
                 label//" gradient is not (3, n_atoms)")
      if (allocated(error)) return

      call check(error, all(ieee_is_finite(result%gradient)), &
                 label//" gradient contains a non-finite element")
   end subroutine check_gradient_shape

   subroutine check_hessian_shape(error, result, n_atoms, label)
      !! The same, for the second derivatives
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), intent(in) :: result
      integer, intent(in) :: n_atoms
      character(len=*), intent(in) :: label

      call check(error, allocated(result%hessian), &
                 label//" claims a Hessian but did not allocate one")
      if (allocated(error)) return

      call check(error, size(result%hessian, 1) == 3*n_atoms .and. &
                 size(result%hessian, 2) == 3*n_atoms, &
                 label//" Hessian is not (3 n_atoms, 3 n_atoms)")
      if (allocated(error)) return

      call check(error, all(ieee_is_finite(result%hessian)), &
                 label//" Hessian contains a non-finite element")
   end subroutine check_hessian_shape

   subroutine test_mcscf_no_active(error)
      type(error_type), allocatable, intent(out) :: error
      type(mcscf_method_t) :: method
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call create_test_fragment(fragment)

      method%options%n_active_electrons = 0
      method%options%n_active_orbitals = 0

      call method%calc_energy(fragment, result)

      call check(error, result%has_error, &
                 "MCSCF without active space should set has_error")

      call fragment%destroy()
   end subroutine test_mcscf_no_active

   subroutine test_mcscf_with_active(error)
      type(error_type), allocatable, intent(out) :: error
      type(mcscf_method_t) :: method
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call create_test_fragment(fragment)

      method%options%verbose = .true.
      method%options%n_active_electrons = 4
      method%options%n_active_orbitals = 4
      method%options%n_states = 2
      method%options%use_pt2 = .true.

      call method%calc_energy(fragment, result)

      call check(error,.not. result%has_error, &
                 "MCSCF with active space should not have error")
      if (allocated(error)) return

      call check(error, result%has_energy, "MCSCF energy should set has_energy")

      call fragment%destroy()
   end subroutine test_mcscf_with_active

   subroutine test_mcscf_gradient(error)
      type(error_type), allocatable, intent(out) :: error
      type(mcscf_method_t) :: method
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call create_test_fragment(fragment)

      method%options%verbose = .true.
      method%options%n_active_electrons = 4
      method%options%n_active_orbitals = 4

      call method%calc_gradient(fragment, result)

      call check(error,.not. result%has_error, &
                 "MCSCF gradient should not have error")
      if (allocated(error)) return

      call check(error, result%has_gradient, "MCSCF gradient should set has_gradient")
      if (allocated(error)) return

      call check(error, allocated(result%gradient), "MCSCF gradient should be allocated")

      call fragment%destroy()
   end subroutine test_mcscf_gradient

end module test_mqc_method_placeholders

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_method_placeholders, only: collect_mqc_method_placeholders_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_method_placeholders", collect_mqc_method_placeholders_tests) &
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
