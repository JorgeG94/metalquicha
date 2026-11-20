!! Cartesian polarizability derivatives, against PySCF and against themselves
module test_mqc_czt_raman
   !! The pinned numbers are PySCF's coupled-perturbed polarizability of
   !! water in 6-31G, differenced over the same displacement by a script that
   !! shares nothing with this program: three coordinates on three different
   !! atoms, six tensor components each.
   !!
   !! Beside them, what needs no other program: moving every atom together
   !! does not change a polarizability, so the derivatives with respect to
   !! one direction must sum to zero over the atoms, for each component. They
   !! do so to the error of the central differences and no better, which
   !! makes this the one place that error is measured.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_raman, only: polarizability_derivatives, polarizability_derivatives_field, &
                            polarizability_derivatives_analytic
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_scf_convergence, only: scf_convergence_t
   use mqc_diis, only: ACCEL_DIIS
   implicit none
   private

   public :: collect_mqc_czt_raman_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   ! PySCF, water/6-31G. Components xx, yy, zz, xy, xz, yz; Bohr^3 and Bohr^2.
   real(dp), parameter :: ALPHA_REF(4) = &
                          [5.2420415435_dp, 1.3940576432_dp, 5.8188035990_dp, -1.0842638457_dp]
      !! xx, yy, zz and xz of the polarizability itself
   real(dp), parameter :: D_O_Z(6) = [-0.0629687918_dp, 0.3506691657_dp, -7.9166878450_dp, &
                                      0.0_dp, -1.3439554052_dp, 0.0_dp]
   real(dp), parameter :: D_H1_Z(6) = [0.3723992802_dp, -0.4427933520_dp, 8.5097308792_dp, &
                                       0.0_dp, -0.5980592160_dp, 0.0_dp]
   real(dp), parameter :: D_H2_X(6) = [7.8657959733_dp, -0.4342201375_dp, 1.0614212215_dp, &
                                       0.0_dp, -2.0744888267_dp, 0.0_dp]
   real(dp), parameter :: FIELD_TOL = 1.0e-4_dp
      !! The two routes' truncation errors, 3e-5 each, with room
   real(dp), parameter :: TOL = 2.0e-6_dp
      !! Two independent response solves at 1e-9, divided by a step of 5e-3.
   ! The same three columns with the step's error removed: PySCF's
   ! polarizability differenced at 1e-2 and 5e-3 Bohr and extrapolated to zero.
   real(dp), parameter :: X_O_Z(6) = [-0.0629616879_dp, 0.3506727893_dp, -7.9166845431_dp, &
                                      0.0_dp, -1.3439660108_dp, 0.0_dp]
   real(dp), parameter :: X_H1_Z(6) = [0.3723987273_dp, -0.4427983667_dp, 8.5097162996_dp, &
                                       0.0_dp, -0.5980602745_dp, 0.0_dp]
   real(dp), parameter :: X_H2_X(6) = [7.8657824311_dp, -0.4342248504_dp, 1.0614212240_dp, &
                                       0.0_dp, -2.0744928878_dp, 0.0_dp]
   real(dp), parameter :: ANALYTIC_TOL = 1.0e-6_dp
      !! The extrapolated reference's own remaining error; 3.5e-7 measured
   real(dp), parameter :: SUM_RULE_TOL = 1.0e-7_dp
      !! What the analytic route leaves of an exact zero

contains

   subroutine collect_mqc_czt_raman_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("water_derivatives_match_pyscf", test_water), &
                  new_unittest("finite_field_route_agrees", test_field), &
                  new_unittest("analytic_route_agrees", test_analytic) &
                  ]
   end subroutine collect_mqc_czt_raman_tests

   subroutine test_water(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(cuest_scf_settings_t) :: settings
      type(error_t) :: err
      type(scf_convergence_t) :: conv
      real(dp), allocatable :: dalpha(:, :)
      real(dp) :: c(3, 3), alpha(3, 3), drift
      integer :: comp, k

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
      call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err, &
                       grad_tol=1.0e-10_dp)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "water must converge: "//err%get_full_trace())
      if (allocated(error)) return

      settings%basis_set = "6-31g"
      settings%max_iter = 200
      call polarizability_derivatives(settings, [8, 1, 1], ["O ", "H ", "H "], c, 10, &
                                      scf%density, ACCEL_DIIS, 8, conv, alpha, dalpha, err)
      call check(error,.not. err%has_error(), "derivatives: "//err%get_full_trace())
      if (allocated(error)) return

      write (*, "(4x,a,4f14.8)") "alpha xx yy zz xz ", alpha(1, 1), alpha(2, 2), &
         alpha(3, 3), alpha(1, 3)
      write (*, "(4x,a,6f12.7)") "d/dz(O)  ", dalpha(:, 3)
      write (*, "(4x,a,6f12.7)") "d/dz(H1) ", dalpha(:, 6)
      write (*, "(4x,a,6f12.7)") "d/dx(H2) ", dalpha(:, 7)

      call check(error, abs(alpha(1, 1) - ALPHA_REF(1)) < 1.0e-7_dp &
                 .and. abs(alpha(2, 2) - ALPHA_REF(2)) < 1.0e-7_dp &
                 .and. abs(alpha(3, 3) - ALPHA_REF(3)) < 1.0e-7_dp &
                 .and. abs(alpha(1, 3) - ALPHA_REF(4)) < 1.0e-7_dp, &
                 "the polarizability does not match PySCF's")
      if (allocated(error)) return
      call check(error, maxval(abs(dalpha(:, 3) - D_O_Z)) < TOL, &
                 "the derivative along z of oxygen does not match PySCF's")
      if (allocated(error)) return
      call check(error, maxval(abs(dalpha(:, 6) - D_H1_Z)) < TOL, &
                 "the derivative along z of the first hydrogen does not match PySCF's")
      if (allocated(error)) return
      call check(error, maxval(abs(dalpha(:, 7) - D_H2_X)) < TOL, &
                 "the derivative along x of the second hydrogen does not match PySCF's")
      if (allocated(error)) return

      ! Translation: for each direction and each component, the three atoms'
      ! derivatives cancel.
      drift = 0.0_dp
      do comp = 1, 3
         do k = 1, 6
            drift = max(drift, abs(dalpha(k, comp) + dalpha(k, 3 + comp) + dalpha(k, 6 + comp)))
         end do
      end do
      write (*, "(4x,a,es10.2)") "translational sum rule ", drift
      ! Not to `TOL`: each derivative is a central difference, and three of
      ! them sum to the translation's derivative only up to the term
      ! quadratic in the step. That term is what this measures -- 2.8e-5
      ! here -- so the bound is on the differences' own error, which the
      ! comparison above cannot see because the reference shares the step.
      call check(error, drift < 1.0e-4_dp, &
                 "the polarizability changed when the whole molecule was moved")
      call mol%destroy()
   end subroutine test_water

   subroutine test_field(error)
      !! The finite-field route against the same PySCF numbers and against
      !! the displacement route.
      !!
      !! The two routes share nothing past the reference SCF: one differences
      !! a polarizability over nuclear positions, the other differences a
      !! gradient twice over field strengths. They are different finite
      !! differences of the same third derivative, so they agree to the larger
      !! of their two truncation errors and no better, and the bound here is
      !! that error measured, not a tolerance chosen.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(cuest_scf_settings_t) :: settings
      type(scf_convergence_t) :: conv
      type(error_t) :: err
      real(dp), allocatable :: by_field(:, :), by_moving(:, :)
      real(dp) :: c(3, 3), alpha(3, 3), alpha_moving(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
      call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err, &
                       grad_tol=1.0e-10_dp)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "water must converge: "//err%get_full_trace())
      if (allocated(error)) return

      settings%basis_set = "6-31g"
      settings%max_iter = 200
      call polarizability_derivatives_field(settings, mol, 10, scf%density, ACCEL_DIIS, 8, &
                                            conv, alpha, by_field, err)
      call check(error,.not. err%has_error(), "field route: "//err%get_full_trace())
      if (allocated(error)) return
      call polarizability_derivatives(settings, [8, 1, 1], ["O ", "H ", "H "], c, 10, &
                                      scf%density, ACCEL_DIIS, 8, conv, alpha_moving, &
                                      by_moving, err)
      call check(error,.not. err%has_error(), "displacement route: "//err%get_full_trace())
      if (allocated(error)) return

      write (*, "(4x,a,6f12.7)") "d/dz(O)  by field ", by_field(:, 3)
      write (*, "(4x,a,6f12.7)") "d/dz(H1) by field ", by_field(:, 6)
      write (*, "(4x,a,6f12.7)") "d/dx(H2) by field ", by_field(:, 7)
      write (*, "(4x,a,es10.2)") "largest difference between the two routes ", &
         maxval(abs(by_field - by_moving))
      write (*, "(4x,a,es10.2)") "largest difference from PySCF            ", &
         max(maxval(abs(by_field(:, 3) - D_O_Z)), maxval(abs(by_field(:, 6) - D_H1_Z)), &
             maxval(abs(by_field(:, 7) - D_H2_X)))

      call check(error, maxval(abs(alpha - alpha_moving)) < 1.0e-8_dp, &
                 "the two routes report different polarizabilities")
      if (allocated(error)) return
      call check(error, maxval(abs(by_field - by_moving)) < FIELD_TOL, &
                 "the finite-field derivatives do not match the displacement ones")
      if (allocated(error)) return
      call check(error, maxval(abs(by_field(:, 3) - D_O_Z)) < FIELD_TOL &
                 .and. maxval(abs(by_field(:, 6) - D_H1_Z)) < FIELD_TOL &
                 .and. maxval(abs(by_field(:, 7) - D_H2_X)) < FIELD_TOL, &
                 "the finite-field derivatives do not match PySCF's")
      call mol%destroy()
   end subroutine test_field

   subroutine test_analytic(error)
      !! The analytic derivatives, which difference nothing.
      !!
      !! Two checks of different kinds. PySCF's differenced polarizabilities
      !! carry a truncation error of 3e-5, more than this route's own, so the
      !! reference is those differences at two steps extrapolated to zero.
      !! The translational sum rule needs no reference: moving every atom
      !! together moves nothing, so each component's derivatives sum to zero
      !! over the atoms, and for a finite difference that sum is its
      !! truncation error while here it is the response solves' residual.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: d(:, :)
      real(dp) :: c(3, 3), alpha(3, 3), drift
      integer :: comp

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
      call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err, &
                       grad_tol=1.0e-10_dp)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "water must converge: "//err%get_full_trace())
      if (allocated(error)) return

      call polarizability_derivatives_analytic(mol, scf%orbitals, scf%orbital_energies, &
                                               scf%n_occupied, alpha, d, err)
      call check(error,.not. err%has_error(), "analytic route: "//err%get_full_trace())
      if (allocated(error)) return

      drift = 0.0_dp
      do comp = 1, 3
         drift = max(drift, maxval(abs(d(:, comp) + d(:, 3 + comp) + d(:, 6 + comp))))
      end do
      write (*, "(4x,a,6f12.7)") "d/dz(O)  analytic ", d(:, 3)
      write (*, "(4x,a,6f12.7)") "d/dz(H1) analytic ", d(:, 6)
      write (*, "(4x,a,6f12.7)") "d/dx(H2) analytic ", d(:, 7)
      write (*, "(4x,a,es10.2)") "worst difference from PySCF, extrapolated:  ", &
         max(maxval(abs(d(:, 3) - X_O_Z)), maxval(abs(d(:, 6) - X_H1_Z)), &
             maxval(abs(d(:, 7) - X_H2_X)))
      write (*, "(4x,a,es10.2)") "translational sum rule:                    ", drift

      call check(error, abs(alpha(1, 1) - ALPHA_REF(1)) < TOL &
                 .and. abs(alpha(2, 2) - ALPHA_REF(2)) < TOL &
                 .and. abs(alpha(3, 3) - ALPHA_REF(3)) < TOL &
                 .and. abs(alpha(1, 3) - ALPHA_REF(4)) < TOL, &
                 "the polarizability does not match PySCF's")
      if (allocated(error)) return
      call check(error, maxval(abs(d(:, 3) - X_O_Z)) < ANALYTIC_TOL &
                 .and. maxval(abs(d(:, 6) - X_H1_Z)) < ANALYTIC_TOL &
                 .and. maxval(abs(d(:, 7) - X_H2_X)) < ANALYTIC_TOL, &
                 "the analytic derivatives do not match PySCF's")
      if (allocated(error)) return
      call check(error, drift < SUM_RULE_TOL, &
                 "the analytic derivatives do not sum to zero over the atoms")
      if (allocated(error)) return
      ! Symmetry forbids these outright, and with nothing differenced they
      ! are zero rather than a truncation error.
      call check(error, abs(d(4, 3)) < SUM_RULE_TOL .and. abs(d(6, 3)) < SUM_RULE_TOL, &
                 "a component symmetry forbids is not zero")
   end subroutine test_analytic

end module test_mqc_czt_raman

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_raman, only: collect_mqc_czt_raman_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_raman", collect_mqc_czt_raman_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
