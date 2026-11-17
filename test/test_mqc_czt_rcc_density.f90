module test_mqc_czt_rcc_density
   !! The CCSD one-particle density against PySCF's.
   !!
   !! A density matrix element changes sign with the sign of either orbital, so
   !! it is compared by what does not: the trace of the occupied block, the sum
   !! of the squares of every element, and the trace against the orbital
   !! energies, which is the one-particle part of the correlation energy. Two
   !! exact properties are held as well, to rounding and with no reference: the
   !! matrix is symmetric, and it has no trace, because the reference
   !! determinant's electrons are not in it and correlation makes none.
   !!
   !! References are `pyscf.cc.ccsd_rdm._make_rdm1` with `with_mf=False`, over
   !! `ccsd` for conventional integrals and `dfccsd` for fitted ones.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_rcc, only: rcc_result_t, rcc_state_t, run_czt_rccsd
   use mqc_czt_rcc_lambda, only: rcc_lambda_t, run_czt_rccsd_lambda
   use mqc_czt_rcc_density, only: rccsd_one_particle_density
   implicit none
   private
   public :: collect_mqc_czt_rcc_density_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: TOL = 2.0e-8_dp

contains

   subroutine collect_mqc_czt_rcc_density_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("density_conventional_sto3g", test_conventional_sto3g), &
                  new_unittest("density_conventional_cc_pvdz", test_conventional_pvdz), &
                  new_unittest("density_conventional_frozen_core", test_conventional_frozen), &
                  new_unittest("density_fitted_frozen_core", test_fitted_frozen) &
                  ]
   end subroutine collect_mqc_czt_rcc_density_tests

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
      real(dp), intent(in) :: reference(3)
         !! Trace of the occupied block, sum of squares, trace with the orbital
         !! energies
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(rcc_result_t) :: cc
      type(rcc_state_t) :: state
      type(rcc_lambda_t) :: lambda
      type(error_t) :: err
      real(dp), allocatable :: gamma(:, :)
      real(dp) :: got(3)
      character(len=*), parameter :: NAMES(3) = ["occupied trace ", "sum of squares ", &
                                                 "energy trace   "]
      integer :: k, p, no

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

      call rccsd_one_particle_density(state, lambda, gamma)
      no = state%no
      call check(error, maxval(abs(gamma - transpose(gamma))) < 1.0e-13_dp, &
                 "the density is not symmetric in "//basis)
      if (allocated(error)) return
      got = 0.0_dp
      do p = 1, no
         got(1) = got(1) + gamma(p, p)
         got(3) = got(3) + gamma(p, p)*state%eps_o(p)
      end do
      call check(error, abs(got(1) + sum([(gamma(p, p), p=no + 1, size(gamma, 1))])) &
                 < 1.0e-12_dp, "the correlation density has a trace in "//basis)
      if (allocated(error)) return
      do p = 1, state%nv
         got(3) = got(3) + gamma(no + p, no + p)*state%eps_v(p)
      end do
      got(2) = sum(gamma**2)
      do k = 1, 3
         call check(error, abs(got(k) - reference(k)) < TOL, &
                    trim(NAMES(k))//" differs from PySCF's by "// &
                    real_to_text(got(k) - reference(k))//" in "//basis)
         if (allocated(error)) return
      end do
      call mol%destroy()
   end subroutine compare

   subroutine test_conventional_sto3g(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("sto-3g", 0, .false., [-0.052472153713_dp, 0.003875182014_dp, &
                                          0.069674205689_dp], error)
   end subroutine test_conventional_sto3g

   subroutine test_conventional_pvdz(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("cc-pvdz", 0, .false., [-0.111391546180_dp, 0.006905712254_dp, &
                                           0.225793627516_dp], error)
   end subroutine test_conventional_pvdz

   subroutine test_conventional_frozen(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("cc-pvdz", 1, .false., [-0.111402441778_dp, 0.006923403222_dp, &
                                           0.223986675289_dp], error)
   end subroutine test_conventional_frozen

   subroutine test_fitted_frozen(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("cc-pvdz", 1, .true., [-0.111514889082_dp, 0.006932680085_dp, &
                                          0.224326411618_dp], error)
   end subroutine test_fitted_frozen

   function real_to_text(value) result(text)
      !! A number in a failure message, so the message says how badly
      real(dp), intent(in) :: value
      character(len=:), allocatable :: text

      character(len=32) :: buffer

      write (buffer, "(es16.8)") value
      text = trim(adjustl(buffer))
   end function real_to_text

end module test_mqc_czt_rcc_density

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_rcc_density, only: collect_mqc_czt_rcc_density_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_rcc_density", collect_mqc_czt_rcc_density_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
