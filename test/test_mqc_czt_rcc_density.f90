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
   use mqc_czt_rcc_triples_terms, only: rcc_triples_terms_t, rccsd_triples_terms
   use mqc_czt_rcc_density, only: rccsd_one_particle_density, &
                                  rccsd_fitted_two_particle_density, rcc_active_fit
   implicit none
   private
   public :: collect_mqc_czt_rcc_density_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: TOL = 2.0e-8_dp
   integer, parameter :: N_TRIPLES_TERMS = 10
      !! Numbers the triples' terms are compared by
   integer, parameter :: N_TRIPLES_REFERENCES = N_TRIPLES_TERMS + 1
      !! Those, and the correction itself

contains

   subroutine collect_mqc_czt_rcc_density_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("density_conventional_sto3g", test_conventional_sto3g), &
                  new_unittest("density_conventional_cc_pvdz", test_conventional_pvdz), &
                  new_unittest("density_conventional_frozen_core", test_conventional_frozen), &
                  new_unittest("density_fitted_frozen_core", test_fitted_frozen), &
                  new_unittest("two_particle_density_sto3g", test_two_particle_sto3g), &
                  new_unittest("two_particle_density_cc_pvdz", test_two_particle_pvdz), &
                  new_unittest("two_particle_density_frozen_core", test_two_particle_frozen), &
                  new_unittest("two_particle_density_needs_a_fit", test_needs_fit), &
                  new_unittest("triples_terms_sto3g", test_triples_sto3g), &
                  new_unittest("triples_terms_cc_pvdz_frozen_core", test_triples_pvdz) &
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

   subroutine compare_two_particle(basis, frozen, reference, error)
      !! The fitted two-particle density, by two things that need no element
      !!
      !! **The energy.** The correlation energy is the one-particle density
      !! against the Fock matrix plus half the two-particle density against the
      !! integrals, and the integrals are the fit times itself. So the two
      !! densities built here must give back the energy the amplitudes
      !! converged to. Nothing external enters, and every block of the density
      !! contributes: a sign or a factor wrong in one of them moves the sum.
      !!
      !! **The norm.** The energy weights each block by its integrals and would
      !! miss an error that the integrals happen to cancel, so the sum of
      !! squares is compared with PySCF's as well. It does not depend on the
      !! sign of an orbital or on which factor of the metric the fit was built
      !! with.
      character(len=*), intent(in) :: basis
      integer, intent(in) :: frozen
      real(dp), intent(in) :: reference   !! Sum of squares of the contracted density
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(rcc_result_t) :: cc
      type(rcc_state_t) :: state
      type(rcc_lambda_t) :: lambda
      type(error_t) :: err
      real(dp), allocatable :: gamma(:, :), g3(:, :, :), b(:, :, :)
      real(dp) :: e_one, e_two
      integer :: p

      call water(mol, err, basis)
      call water(aux, err, "cc-pvdz-rifit")
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return
      call run_czt_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the SCF should converge")
      if (allocated(error)) return
      call run_czt_rccsd(mol, scf%orbitals, scf%orbital_energies, 5, frozen, 80, &
                         1.0e-12_dp, .false., .false., cc, err, aux=aux, state=state)
      call check(error,.not. err%has_error() .and. cc%converged, "CCSD should converge")
      if (allocated(error)) return
      call run_czt_rccsd_lambda(state, 80, 1.0e-10_dp, .false., lambda, err)
      call check(error,.not. err%has_error() .and. lambda%converged, &
                 "the Lambda equations should converge")
      if (allocated(error)) return

      call rccsd_one_particle_density(state, lambda, gamma)
      call rcc_active_fit(state, b, err)
      call rccsd_fitted_two_particle_density(state, lambda, b, g3, err)
      call check(error,.not. err%has_error(), "the density failed: "//err%get_message())
      if (allocated(error)) return

      e_one = 0.0_dp
      do p = 1, state%no
         e_one = e_one + gamma(p, p)*state%eps_o(p)
      end do
      do p = 1, state%nv
         e_one = e_one + gamma(state%no + p, state%no + p)*state%eps_v(p)
      end do
      e_two = 0.5_dp*sum(g3*b)
      call check(error, abs(e_one + e_two - cc%e_correlation) < 1.0e-9_dp, &
                 "the densities give back a correlation energy off by "// &
                 real_to_text(e_one + e_two - cc%e_correlation)//" in "//basis)
      if (allocated(error)) return
      call check(error, maxval(abs(g3 - reshape(g3, shape(g3), order=[2, 1, 3]))) &
                 < 1.0e-13_dp, "the contracted density is not symmetric in "//basis)
      if (allocated(error)) return
      call check(error, abs(sum(g3**2) - reference) < TOL, &
                 "the sum of squares differs from PySCF's by "// &
                 real_to_text(sum(g3**2) - reference)//" in "//basis)
      call mol%destroy()
      call aux%destroy()
   end subroutine compare_two_particle

   subroutine test_two_particle_sto3g(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_two_particle("sto-3g", 0, 0.037657971213_dp, error)
   end subroutine test_two_particle_sto3g

   subroutine test_two_particle_pvdz(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_two_particle("cc-pvdz", 0, 0.175031155204_dp, error)
   end subroutine test_two_particle_pvdz

   subroutine test_two_particle_frozen(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_two_particle("cc-pvdz", 1, 0.174157632852_dp, error)
   end subroutine test_two_particle_frozen

   subroutine test_needs_fit(error)
      !! Conventional integrals leave no three-index tensor, and that is said
      type(error_type), allocatable, intent(out) :: error

      type(rcc_state_t) :: state
      type(rcc_lambda_t) :: lambda
      type(error_t) :: err
      real(dp), allocatable :: g3(:, :, :), b(:, :, :)

      call rcc_active_fit(state, b, err)
      call check(error, err%has_error(), "a fit was handed on where none was built")
      if (allocated(error)) return
      allocate (b(1, 1, 1), source=0.0_dp)
      call rccsd_fitted_two_particle_density(state, lambda, b, g3, err)
      call check(error, err%has_error(), "a density was built with no amplitudes")
   end subroutine test_needs_fit

   subroutine compare_triples(basis, frozen, reference, error)
      !! What (T) adds to the Lambda equations and to the densities
      !!
      !! **Term by term**, against `pyscf.cc.ccsd_t_lambda_slow` and
      !! `ccsd_t_rdm_slow`: the two Lambda sources by their norms and their
      !! overlaps with the amplitudes, the one-particle diagonal by its sum and
      !! its sum of squares in each space, and the three two-particle blocks by
      !! their norms. None depends on the sign of an orbital.
      !!
      !! **Together**, with no reference. The correction's denominators are
      !! orbital energies and its numerators are integrals, so the energy is
      !! homogeneous in the pair, and the one- and two-particle densities
      !! against the Fock matrix and the integrals must give back the CCSD
      !! energy *and* the triples correction. The diagonal the correction adds
      !! to the one-particle density must also have no trace.
      character(len=*), intent(in) :: basis
      integer, intent(in) :: frozen
      real(dp), intent(in) :: reference(N_TRIPLES_REFERENCES)
         !! Sums of squares of the `l1` and `l2` sources; their overlaps with
         !! `t1` and `t2`; the sum of the occupied diagonal, negated; the sums
         !! of squares of the occupied and virtual diagonals; and the sums of
         !! squares of the additions to `dovov`, `dooov` and `dovvv`. The last
         !! is the correction itself.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(rcc_result_t) :: cc
      type(rcc_state_t) :: state
      type(rcc_lambda_t) :: lambda
      type(rcc_triples_terms_t) :: terms
      type(error_t) :: err
      real(dp), allocatable :: gamma(:, :), g3(:, :, :), b(:, :, :)
      real(dp) :: got(N_TRIPLES_TERMS)
      real(dp) :: e_one, e_two
      character(len=*), parameter :: NAMES(N_TRIPLES_TERMS) = [character(len=24) :: &
                                           "l1 source norm", "l2 source norm", "l1 source with t1", "l2 source with t2", &
                                             "occupied diagonal sum", "occupied diagonal norm", "virtual diagonal norm", &
                                                               "dovov addition", "dooov addition", "dovvv addition"]
      integer :: k, p

      call water(mol, err, basis)
      call water(aux, err, "cc-pvdz-rifit")
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return
      call run_czt_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the SCF should converge")
      if (allocated(error)) return
      call run_czt_rccsd(mol, scf%orbitals, scf%orbital_energies, 5, frozen, 80, &
                         1.0e-12_dp, .true., .false., cc, err, aux=aux, state=state)
      call check(error,.not. err%has_error() .and. cc%converged, "CCSD(T) should converge")
      if (allocated(error)) return
      call check(error, abs(cc%e_triples - reference(N_TRIPLES_REFERENCES)) < TOL, &
                 "the triples correction itself differs from PySCF's in "//basis)
      if (allocated(error)) return

      call rccsd_triples_terms(state, terms)
      got(1) = sum(terms%l1_source**2)
      got(2) = sum(terms%l2_source**2)
      got(3) = sum(terms%l1_source*state%t1)
      got(4) = sum(terms%l2_source*state%t2)
      got(5) = -sum(terms%gamma_occ)
      got(6) = sum(terms%gamma_occ**2)
      got(7) = sum(terms%gamma_vir**2)
      got(8) = sum(terms%dovov**2)
      got(9) = sum(terms%dooov**2)
      got(10) = sum(terms%dovvv**2)
      do k = 1, N_TRIPLES_TERMS
         call check(error, abs(got(k) - reference(k)) < TOL*max(1.0_dp, abs(reference(k))) &
                    .and. abs(got(k) - reference(k)) < 1.0e-6_dp*abs(reference(k)), &
                    trim(NAMES(k))//" differs from PySCF's by "// &
                    real_to_text(got(k) - reference(k))//" in "//basis)
         if (allocated(error)) return
      end do
      call check(error, abs(sum(terms%gamma_occ) + sum(terms%gamma_vir)) < 1.0e-13_dp, &
                 "the triples' one-particle diagonal has a trace in "//basis)
      if (allocated(error)) return

      call run_czt_rccsd_lambda(state, 200, 1.0e-10_dp, .false., lambda, err, triples=terms)
      call check(error,.not. err%has_error() .and. lambda%converged, &
                 "the Lambda equations should converge with their triples sources")
      if (allocated(error)) return
      call rccsd_one_particle_density(state, lambda, gamma, terms)
      call rcc_active_fit(state, b, err)
      call rccsd_fitted_two_particle_density(state, lambda, b, g3, err, terms)
      call check(error,.not. err%has_error(), "the density failed: "//err%get_message())
      if (allocated(error)) return
      e_one = 0.0_dp
      do p = 1, state%no
         e_one = e_one + gamma(p, p)*state%eps_o(p)
      end do
      do p = 1, state%nv
         e_one = e_one + gamma(state%no + p, state%no + p)*state%eps_v(p)
      end do
      e_two = 0.5_dp*sum(g3*b)
      call check(error, abs(e_one + e_two - cc%e_correlation) < 1.0e-9_dp, &
                 "the CCSD(T) densities give back a correlation energy off by "// &
                 real_to_text(e_one + e_two - cc%e_correlation)//" in "//basis)
      call mol%destroy()
      call aux%destroy()
   end subroutine compare_triples

   subroutine test_triples_sto3g(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_triples("sto-3g", 0, &
                           [1.450961230632e-07_dp, 8.797389984237e-07_dp, &
                            -4.686920165586e-06_dp, 5.681792254221e-05_dp, &
                            5.591504728020e-05_dp, 7.977139787946e-10_dp, &
                            1.587261986887e-09_dp, 1.578094046696e-09_dp, &
                            6.252371391667e-07_dp, 7.236460319319e-07_dp, &
                            -0.000072751439_dp], error)
   end subroutine test_triples_sto3g

   subroutine test_triples_pvdz(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_triples("cc-pvdz", 1, &
                           [2.854075656811e-06_dp, 9.766911826857e-05_dp, &
                            -2.346618685550e-05_dp, 1.502037729958e-03_dp, &
                            1.447772647282e-03_dp, 5.794903305181e-07_dp, &
                            1.520217469078e-07_dp, 3.342912848745e-08_dp, &
                            9.163956572795e-06_dp, 3.499270962604e-05_dp, &
                            -0.003019775465_dp], error)
   end subroutine test_triples_pvdz

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
