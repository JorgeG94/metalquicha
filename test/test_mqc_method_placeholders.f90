module test_mqc_method_placeholders
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

      ! Analytic HF gradients are not implemented yet. The contract is that
      ! this is reported, not that a zero gradient is handed back as if real.
      call check(error, result%has_error, &
                 "HF gradients are unimplemented and must report an error")
      if (allocated(error)) return

      call check(error,.not. result%has_gradient, &
                 "An unimplemented HF gradient must not claim to have one")

      call fragment%destroy()
   end subroutine test_hf_gradient

   subroutine test_hf_hessian(error)
      type(error_type), allocatable, intent(out) :: error
      type(hf_method_t) :: method
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call create_test_fragment(fragment)
      call method%calc_hessian(fragment, result)

      call check(error,.not. result%has_hessian, &
                 "HF Hessian placeholder should report not implemented")

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

      ! Analytic DFT gradients are not implemented; the contract is that this
      ! is reported rather than a zero gradient being passed off as real.
      call check(error, result%has_error, &
                 "DFT gradients are unimplemented and must report an error")
      if (allocated(error)) return

      call check(error,.not. result%has_gradient, &
                 "An unimplemented DFT gradient must not claim to have one")

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

      ! Analytic DFT Hessians are not implemented; report, never fabricate.
      call check(error, result%has_error, &
                 "DFT Hessians are unimplemented and must report an error")
      if (allocated(error)) return

      call check(error,.not. result%has_hessian, &
                 "An unimplemented DFT Hessian must not claim to have one")

      call fragment%destroy()
   end subroutine test_dft_hessian

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
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester
