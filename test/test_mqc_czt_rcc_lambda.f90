module test_mqc_czt_rcc_lambda
   !! The Lambda amplitudes of spin-adapted CCSD against PySCF's.
   !!
   !! An amplitude cannot be compared element by element: the sign of every
   !! orbital is arbitrary, and so is the sign of any element with an odd number
   !! of indices on one orbital. What is compared is what does not depend on
   !! that -- the squared norms of `l1` and `l2`, and their overlaps with `t1`
   !! and `t2`. Four numbers per case, each summed over every element, so an
   !! index transposed in one contraction moves all of them.
   !!
   !! The cases are the ones the energy is tested on next door: conventional and
   !! fitted integrals, a minimal basis and one with a d shell, and a frozen
   !! core. The fitted references are `pyscf.cc.dfccsd` and the conventional
   !! ones `pyscf.cc.ccsd`, both through `solve_lambda`, on the geometry below.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_rcc, only: rcc_result_t, rcc_state_t, run_czt_rccsd
   use mqc_czt_rcc_lambda, only: rcc_lambda_t, run_czt_rccsd_lambda
   implicit none
   private
   public :: collect_mqc_czt_rcc_lambda_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: TOL = 2.0e-8_dp

contains

   subroutine collect_mqc_czt_rcc_lambda_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("lambda_conventional_sto3g", test_conventional_sto3g), &
                  new_unittest("lambda_conventional_cc_pvdz", test_conventional_pvdz), &
                  new_unittest("lambda_conventional_frozen_core", test_conventional_frozen), &
                  new_unittest("lambda_fitted_sto3g", test_fitted_sto3g), &
                  new_unittest("lambda_fitted_cc_pvdz", test_fitted_pvdz), &
                  new_unittest("lambda_fitted_frozen_core", test_fitted_frozen), &
                  new_unittest("lambda_needs_converged_amplitudes", test_refusal) &
                  ]
   end subroutine collect_mqc_czt_rcc_lambda_tests

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

   subroutine compare(basis, frozen, fitted, reference, error)
      !! Solve CCSD and Lambda, and hold the four invariants to the reference
      character(len=*), intent(in) :: basis
      integer, intent(in) :: frozen
      logical, intent(in) :: fitted
      real(dp), intent(in) :: reference(4)
         !! `l1.l1`, `l2.l2`, `l1.t1`, `l2.t2`, each summed over every element
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(rcc_result_t) :: cc
      type(rcc_state_t) :: state
      type(rcc_lambda_t) :: lambda
      type(error_t) :: err
      real(dp) :: got(4)
      character(len=*), parameter :: NAMES(4) = ["l1.l1", "l2.l2", "l1.t1", "l2.t2"]
      integer :: k

      call water(mol, err, basis)
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return
      call run_czt_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the SCF should converge")
      if (allocated(error)) return

      if (fitted) then
         call water(aux, err, "cc-pvdz-rifit")
         call check(error,.not. err%has_error(), "the auxiliary basis should build")
         if (allocated(error)) return
         call run_czt_rccsd(mol, scf%orbitals, scf%orbital_energies, 5, frozen, 80, &
                            1.0e-12_dp, .false., .false., cc, err, aux=aux, state=state)
         call aux%destroy()
      else
         call run_czt_rccsd(mol, scf%orbitals, scf%orbital_energies, 5, frozen, 80, &
                            1.0e-12_dp, .false., .false., cc, err, state=state)
      end if
      call check(error,.not. err%has_error() .and. cc%converged, &
                 "CCSD should converge: "//err%get_message())
      if (allocated(error)) return

      call run_czt_rccsd_lambda(state, 80, 1.0e-10_dp, .false., lambda, err)
      call check(error,.not. err%has_error() .and. lambda%converged, &
                 "the Lambda equations should converge: "//err%get_message())
      if (allocated(error)) return

      got(1) = sum(lambda%l1**2)
      got(2) = sum(lambda%l2**2)
      got(3) = sum(lambda%l1*state%t1)
      got(4) = sum(lambda%l2*state%t2)
      do k = 1, 4
         call check(error, abs(got(k) - reference(k)) < TOL, &
                    NAMES(k)//" differs from PySCF's by "//real_to_text(got(k) - reference(k))// &
                    " in "//basis)
         if (allocated(error)) return
      end do
      call mol%destroy()
   end subroutine compare

   subroutine test_conventional_sto3g(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("sto-3g", 0, .false., [0.000153229609_dp, 0.025195947997_dp, &
                                          0.000174769023_dp, 0.025558304220_dp], error)
   end subroutine test_conventional_sto3g

   subroutine test_conventional_pvdz(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("cc-pvdz", 0, .false., [0.000175689636_dp, 0.045162392572_dp, &
                                           0.000217088181_dp, 0.045879832593_dp], error)
   end subroutine test_conventional_pvdz

   subroutine test_conventional_frozen(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("cc-pvdz", 1, .false., [0.000177106454_dp, 0.045179218288_dp, &
                                           0.000218612900_dp, 0.045896578227_dp], error)
   end subroutine test_conventional_frozen

   subroutine test_fitted_sto3g(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("sto-3g", 0, .true., [0.000153304717_dp, 0.025199653032_dp, &
                                         0.000174857379_dp, 0.025562120259_dp], error)
   end subroutine test_fitted_sto3g

   subroutine test_fitted_pvdz(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("cc-pvdz", 0, .true., [0.000175720390_dp, 0.045242818226_dp, &
                                          0.000217297541_dp, 0.045961856005_dp], error)
   end subroutine test_fitted_pvdz

   subroutine test_fitted_frozen(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("cc-pvdz", 1, .true., [0.000177135588_dp, 0.045259748875_dp, &
                                          0.000218823172_dp, 0.045978707779_dp], error)
   end subroutine test_fitted_frozen

   subroutine test_refusal(error)
      !! A state no CCSD calculation filled is refused, not dereferenced
      type(error_type), allocatable, intent(out) :: error

      type(rcc_state_t) :: empty
      type(rcc_lambda_t) :: lambda
      type(error_t) :: err

      call run_czt_rccsd_lambda(empty, 10, 1.0e-8_dp, .false., lambda, err)
      call check(error, err%has_error(), "Lambda ran without amplitudes to start from")
   end subroutine test_refusal

   function real_to_text(value) result(text)
      !! A number in a failure message, so the message says how badly
      real(dp), intent(in) :: value
      character(len=:), allocatable :: text

      character(len=32) :: buffer

      write (buffer, "(es16.8)") value
      text = trim(adjustl(buffer))
   end function real_to_text

end module test_mqc_czt_rcc_lambda

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_rcc_lambda, only: collect_mqc_czt_rcc_lambda_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_rcc_lambda", collect_mqc_czt_rcc_lambda_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
