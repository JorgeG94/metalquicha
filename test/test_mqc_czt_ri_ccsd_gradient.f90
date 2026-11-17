module test_mqc_czt_ri_ccsd_gradient
   !! The RI-CCSD gradient against its numpy reference.
   !!
   !! `tools/cpu_validation/ri_ccsd_gradient.py` is the same gradient written
   !! out in numpy over PySCF's integrals, amplitudes and densities, and is
   !! itself checked against central differences of the RI-CCSD energy to
   !! 1.5e-7 hartree/bohr. The numbers below are its output at this geometry.
   !! It shares nothing with the Fortran beyond the equations: not the Lambda
   !! solver, not the density code, not the integral library, not the
   !! coupled-perturbed solver.
   !!
   !! The reference differentiates the energy exactly, so the cases that are
   !! held to it ask for the exact response operator. The default fits that
   !! operator, which moves the gradient by the fitting error, and one case
   !! bounds how far.
   !!
   !! Two things are checked with no reference: the forces on the atoms add up
   !! to nothing, and a molecule in a plane feels no force out of it.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_rcc, only: rcc_result_t, rcc_state_t, run_czt_rccsd
   use mqc_czt_rcc_lambda, only: rcc_lambda_t, run_czt_rccsd_lambda
   use mqc_czt_rcc_triples_terms, only: rcc_triples_terms_t, rccsd_triples_terms
   use mqc_czt_ri_ccsd_gradient, only: czt_ri_ccsd_gradient
   implicit none
   private
   public :: collect_mqc_czt_ri_ccsd_gradient_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

   ! Hartree/bohr, (component, atom), from the numpy reference.
   real(dp), parameter :: G_STO3G(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, -0.1151777821_dp, &
                           0.0_dp, -0.0409816319_dp, 0.0575888911_dp, &
                           0.0_dp, 0.0409816317_dp, 0.0575888910_dp], [3, 3])
   real(dp), parameter :: G_PVDZ_F1(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, -0.0162144691_dp, &
                           0.0_dp, 0.0021651591_dp, 0.0081072346_dp, &
                           0.0_dp, -0.0021651592_dp, 0.0081072345_dp], [3, 3])

   ! The same with the perturbative triples.
   real(dp), parameter :: G_T_STO3G(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, -0.1151431065_dp, &
                           0.0_dp, -0.0410635589_dp, 0.0575715533_dp, &
                           0.0_dp, 0.0410635587_dp, 0.0575715532_dp], [3, 3])
   real(dp), parameter :: G_T_PVDZ_F1(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, -0.0181309958_dp, &
                           0.0_dp, 0.0013577995_dp, 0.0090654979_dp, &
                           0.0_dp, -0.0013577996_dp, 0.0090654979_dp], [3, 3])

contains

   subroutine collect_mqc_czt_ri_ccsd_gradient_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("gradient_sto3g_matches_reference", test_sto3g), &
                  new_unittest("gradient_cc_pvdz_frozen_core_matches_reference", test_pvdz), &
                  new_unittest("fitted_response_moves_it_by_the_fitting_error", test_fitted), &
                  new_unittest("triples_gradient_sto3g_matches_reference", test_t_sto3g), &
                  new_unittest("triples_gradient_cc_pvdz_frozen_core", test_t_pvdz) &
                  ]
   end subroutine collect_mqc_czt_ri_ccsd_gradient_tests

   subroutine water(mol, err, basis)
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      character(len=*), intent(in) :: basis
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.10077199490609_dp*ANG, &
                   0.0_dp, 0.77250895271063_dp*ANG, -0.46780199741728_dp*ANG, &
                   0.0_dp, -0.77250895280218_dp*ANG, -0.46780199748881_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
   end subroutine water

   subroutine gradient_of(basis, frozen, fit_response, gradient, error, with_triples)
      !! SCF, CCSD, Lambda and the gradient, for water
      character(len=*), intent(in) :: basis
      integer, intent(in) :: frozen
      logical, intent(in) :: fit_response
      logical, intent(in), optional :: with_triples   !! The gradient of CCSD(T)
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(rcc_result_t) :: cc
      type(rcc_state_t) :: state
      type(rcc_lambda_t) :: lambda
      type(rcc_triples_terms_t) :: terms
      type(error_t) :: err
      logical :: triples

      triples = .false.
      if (present(with_triples)) triples = with_triples
      call water(mol, err, basis)
      call water(aux, err, "cc-pvdz-rifit")
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return
      call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the SCF should converge")
      if (allocated(error)) return
      call run_czt_rccsd(mol, scf%orbitals, scf%orbital_energies, 5, frozen, 100, &
                         1.0e-13_dp, .false., .false., cc, err, aux=aux, state=state)
      call check(error,.not. err%has_error() .and. cc%converged, "CCSD should converge")
      if (allocated(error)) return
      if (triples) then
         call rccsd_triples_terms(state, terms)
         call run_czt_rccsd_lambda(state, 200, 1.0e-11_dp, .false., lambda, err, &
                                   triples=terms)
      else
         call run_czt_rccsd_lambda(state, 100, 1.0e-11_dp, .false., lambda, err)
      end if
      call check(error,.not. err%has_error() .and. lambda%converged, &
                 "the Lambda equations should converge")
      if (allocated(error)) return
      if (triples) then
         call czt_ri_ccsd_gradient(mol, aux, scf%orbitals, scf%orbital_energies, 5, &
                                   frozen, state, lambda, gradient, err, &
                                   fit_response=fit_response, triples=terms)
      else
         call czt_ri_ccsd_gradient(mol, aux, scf%orbitals, scf%orbital_energies, 5, &
                                   frozen, state, lambda, gradient, err, &
                                   fit_response=fit_response)
      end if
      call check(error,.not. err%has_error(), "the gradient failed: "//err%get_message())
      call mol%destroy()
      call aux%destroy()
   end subroutine gradient_of

   subroutine hold_to(gradient, reference, tolerance, label, error)
      real(dp), intent(in) :: gradient(:, :), reference(:, :)
      real(dp), intent(in) :: tolerance
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error

      call check(error, maxval(abs(sum(gradient, dim=2))) < 1.0e-9_dp, &
                 "the forces do not add up to nothing, "//label)
      if (allocated(error)) return
      call check(error, maxval(abs(gradient(1, :))) < 1.0e-9_dp, &
                 "there is a force out of the molecule's plane, "//label)
      if (allocated(error)) return
      call check(error, maxval(abs(gradient - reference)) < tolerance, &
                 "the gradient is "//real_to_text(maxval(abs(gradient - reference)))// &
                 " hartree/bohr from the reference, "//label)
   end subroutine hold_to

   subroutine test_sto3g(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: gradient(:, :)

      call gradient_of("sto-3g", 0, .false., gradient, error)
      if (allocated(error)) return
      call hold_to(gradient, G_STO3G, 5.0e-8_dp, "STO-3G", error)
   end subroutine test_sto3g

   subroutine test_pvdz(error)
      !! A d shell, and a frozen core, whose rotations into the active occupied
      !! space the amplitudes never see and the relaxation has to supply
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: gradient(:, :)

      call gradient_of("cc-pvdz", 1, .false., gradient, error)
      if (allocated(error)) return
      call hold_to(gradient, G_PVDZ_F1, 5.0e-8_dp, "cc-pVDZ with a frozen core", error)
   end subroutine test_pvdz

   subroutine test_fitted(error)
      !! The default, which fits the response operator
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: gradient(:, :)

      call gradient_of("cc-pvdz", 1, .true., gradient, error)
      if (allocated(error)) return
      call hold_to(gradient, G_PVDZ_F1, 1.0e-5_dp, "cc-pVDZ, fitted response", error)
   end subroutine test_fitted

   subroutine test_t_sto3g(error)
      !! RI-CCSD(T). The correction is not invariant to a rotation among the
      !! occupied orbitals or among the virtuals, as CCSD is, so this also
      !! exercises the relaxation of those two blocks, which CCSD leaves at
      !! what the amplitudes gave.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: gradient(:, :)

      call gradient_of("sto-3g", 0, .false., gradient, error, with_triples=.true.)
      if (allocated(error)) return
      call hold_to(gradient, G_T_STO3G, 5.0e-8_dp, "STO-3G with (T)", error)
   end subroutine test_t_sto3g

   subroutine test_t_pvdz(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: gradient(:, :)

      call gradient_of("cc-pvdz", 1, .false., gradient, error, with_triples=.true.)
      if (allocated(error)) return
      call hold_to(gradient, G_T_PVDZ_F1, 5.0e-8_dp, "cc-pVDZ with a frozen core and (T)", &
                   error)
   end subroutine test_t_pvdz

   function real_to_text(value) result(text)
      !! A number in a failure message, so the message says how badly
      real(dp), intent(in) :: value
      character(len=:), allocatable :: text

      character(len=32) :: buffer

      write (buffer, "(es16.8)") value
      text = trim(adjustl(buffer))
   end function real_to_text

end module test_mqc_czt_ri_ccsd_gradient

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_ri_ccsd_gradient, only: collect_mqc_czt_ri_ccsd_gradient_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_ri_ccsd_gradient", &
                               collect_mqc_czt_ri_ccsd_gradient_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
