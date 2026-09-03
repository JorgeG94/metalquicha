!! The fixed-density second-derivative skeleton, pinned to what already works
module test_mqc_mp2_hessian_skeleton
   !! Unit 1.3 of the MP2 Hessian ladder (`mp2-hessian-phased-plan.md`): the
   !! per-atom-pair scalar `s = Gamma_eff . eri^{XY} + D_rel . h^{XY} +
   !! I . S^{XY}`, checked here without any external reference.
   !!
   !! **Why these two checks and not a pycc comparison.** The cross-code gate
   !! (pycc's `skel_s`, all four oracle configurations, term by term) was run
   !! when this unit landed and its residuals are recorded in the commit; it is
   !! deliberately not stored as a test, for the reason `test_mqc_hess_ints`'
   !! header gives -- a stored comparison against an external library pins our
   !! layout to that library's conventions rather than to anything true. What
   !! stays in-tree is what must not drift silently:
   !!
   !! * the SCF reference skeleton this sweep deposits equals the one
   !!   `rhf_hessian` assembles through `partial_hessian` -- the duplicated
   !!   deposit rules are the same rules, statement for statement;
   !! * the correlation block is symmetric in the atom pair,
   !!   `s(a, b, A, B) = s(b, a, B, A)`, which the assembly does not impose:
   !!   the two elements come from different derivative placements on
   !!   different integrals, so their agreement is earned, not built in.
   !!
   !! One thread throughout: the reference comparison is at 4e-14, where the
   !! unordered OpenMP merges in `hess_2e_contract` scatter more than the gate.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mp2_gradient, only: czt_mp2_gradient, build_amplitudes, &
                                   build_effective_2pdm_ao
   use mqc_czt_mp2, only: transform_ovov
   use mqc_czt_mp2_hessian, only: mp2_skeleton_hessian
   use mqc_czt_hessian, only: response_hessian, rhf_hessian, &
                              nuclear_repulsion_hessian
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   implicit none
   private

   public :: collect_mqc_mp2_hessian_skeleton_tests

   !! pycc's pinned case: water/6-31G, bohr, frame pinned -- the geometry every
   !! unit of the ladder gates on (`test_mqc_mp2_hessian_fd` holds the same).
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.000000_dp, 0.000000_dp, &
                           0.0_dp, 0.000000_dp, 1.814137_dp, &
                           0.0_dp, 1.756000_dp, -0.454300_dp], [3, 3])
   integer, parameter :: WATER_NELEC = 10

contains

   subroutine collect_mqc_mp2_hessian_skeleton_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("reference_block_matches_rhf_hessian", reference_matches), &
                  new_unittest("correlation_block_is_pair_symmetric", pair_symmetric) &
                  ]
   end subroutine collect_mqc_mp2_hessian_skeleton_tests

   !! Everything the skeleton needs, from one SCF and one gradient call: the
   !! effective two-particle density rebuilt the way the gradient builds its
   !! own pieces, and both skeleton blocks.
   subroutine skeleton_at(n_frozen, mol, scf, hess_corr, hess_ref, err)
      integer, intent(in) :: n_frozen
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      real(dp), allocatable, intent(out) :: hess_corr(:, :, :, :)
      real(dp), allocatable, intent(out) :: hess_ref(:, :, :, :)
      type(error_t), intent(inout) :: err

      real(dp), allocatable :: gradient(:, :), dm1mo(:, :), w_ao(:, :)
      real(dp), allocatable :: eri_packed(:, :), ovov(:, :, :, :), t2(:, :, :, :)
      real(dp), allocatable :: gamma_eff(:, :, :, :)
      integer :: n_ao, n_mo, n_o

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "6-31g", mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, .false., &
                       scf, err)
      if (err%has_error()) return

      call czt_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, &
                            WATER_NELEC/2, gradient, err, n_frozen=n_frozen, &
                            relaxed_density_mo=dm1mo, energy_weighted_ao=w_ao)
      if (err%has_error()) return

      n_ao = mol%nao
      n_mo = size(scf%orbitals, 2)
      n_o = WATER_NELEC/2
      call mol%eris_packed(eri_packed)
      call transform_ovov(eri_packed, scf%orbitals, n_frozen, n_o, n_ao, n_mo, ovov)
      call build_amplitudes(ovov, scf%orbital_energies, n_frozen, n_o - n_frozen, &
                            n_mo - n_o, n_o, t2)
      call build_effective_2pdm_ao(t2, dm1mo, scf%orbitals, n_ao, n_mo, n_o, &
                                   n_frozen, gamma_eff)

      call mp2_skeleton_hessian(mol, gamma_eff, dm1mo, w_ao, scf%orbitals, &
                                scf%orbital_energies, n_o, hess_corr, hess_ref, err)
   end subroutine skeleton_at

   !! Gate 1.3b: the reference skeleton this sweep deposits, plus the nuclear
   !! repulsion and the response `rhf_hessian` itself would add, equals
   !! `rhf_hessian` -- so the second-derivative integrals can be generated once
   !! for both densities without the five duplicated deposit lines diverging.
   !!
   !! 1e-13 and not 1e-16: the two-electron walks differ (`hess_2e_contract`
   !! restricts the ket pair and doubles, this sweep runs every ordering), and
   !! two summation orders round differently. Measured 3.8e-14 when this
   !! landed. The one-electron deposits, by contrast, are statement-for-
   !! statement the same and contribute nothing at this tolerance.
   subroutine reference_matches(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: hess_corr(:, :, :, :), hess_ref(:, :, :, :)
      real(dp), allocatable :: nuc(:, :, :, :), resp(:, :, :, :), full(:, :, :, :)
      integer :: threads

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call skeleton_at(0, mol, scf, hess_corr, hess_ref, err)
      call check(error,.not. err%has_error(), &
                 "the skeleton did not evaluate: "//err%get_message())
      if (allocated(error)) then
         call omp_set_num_threads(threads)
         return
      end if

      call nuclear_repulsion_hessian(WATER_Z, mol%coords, nuc, err)
      call response_hessian(mol, scf%density, scf%orbitals, scf%orbital_energies, &
                            WATER_NELEC/2, resp, err)
      call rhf_hessian(mol, WATER_Z, scf%density, scf%orbitals, &
                       scf%orbital_energies, WATER_NELEC/2, full, err)
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), &
                 "the reference Hessian did not evaluate: "//err%get_message())
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      write (*, "(a, es10.3)") "        max |skeleton ref - rhf_hessian| = ", &
         maxval(abs(nuc + hess_ref + resp - full))
      call check(error, maxval(abs(nuc + hess_ref + resp - full)) < 1.0e-13_dp, &
                 "the reference skeleton diverged from rhf_hessian")
      call mol%destroy()
   end subroutine reference_matches

   !! `d2E/dA dB = d2E/dB dA` with the transposed Cartesian pair. The two
   !! elements are assembled from different derivative placements -- `ipvip1`
   !! on one ordering against `ipvip1` on another, `ip1ip2` against its
   !! ket-swapped evaluation -- so this checks the deposit map, not a symmetry
   !! the loop enforces. Measured ~1e-15 all-electron and frozen-core; the
   !! tolerance leaves room for integral rounding, not for a wrong deposit.
   subroutine pair_symmetric(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: hess_corr(:, :, :, :), hess_ref(:, :, :, :)
      real(dp) :: worst
      integer :: threads, ia, ja, a, b

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      ! Frozen core, so the pair map is checked on the fuller assembly: the
      ! core-active coupling in the relaxed density and the projector both
      ! reach every term here.
      call skeleton_at(1, mol, scf, hess_corr, hess_ref, err)
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), &
                 "the skeleton did not evaluate: "//err%get_message())
      if (allocated(error)) return

      worst = 0.0_dp
      do ja = 1, size(hess_corr, 4)
         do ia = 1, size(hess_corr, 3)
            do b = 1, 3
               do a = 1, 3
                  worst = max(worst, abs(hess_corr(a, b, ia, ja) &
                                         - hess_corr(b, a, ja, ia)))
               end do
            end do
         end do
      end do
      write (*, "(a, es10.3)") "        max |s(a,b,A,B) - s(b,a,B,A)| = ", worst
      call check(error, worst < 1.0e-12_dp, &
                 "the correlation skeleton is not symmetric in the atom pair")
      call mol%destroy()
   end subroutine pair_symmetric

end module test_mqc_mp2_hessian_skeleton

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mp2_hessian_skeleton, only: collect_mqc_mp2_hessian_skeleton_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_mp2_hessian_skeleton", &
                               collect_mqc_mp2_hessian_skeleton_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
