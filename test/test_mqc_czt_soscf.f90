!! The second-order SCF: its conventions, its step, and where it lands
module test_mqc_czt_soscf
   !! Four questions, and the first two are the ones everything else rests on.
   !!
   !! **Is the gradient the derivative of the energy?** `soscf_gradient` returns
   !! `4 F_ai`, and both the factor and the sign are conventions -- `C exp(kappa)`
   !! and `C exp(-kappa)` both appear in the literature, and a factor can hide
   !! because it cancels out of a Newton step. So the energy is differentiated
   !! numerically along a rotation and compared with `dot(g, v)`. Nothing is
   !! asserted; the finite difference says.
   !!
   !! **Is the Hessian the second derivative?** Same molecule, same direction,
   !! the second difference of the energy against
   !! `HESSIAN_SCALE * dot(v, (A+B) v)`. This is the check that catches the
   !! factor of four, which the gradient test cannot: it cancels between
   !! gradient and Hessian in the step and shows up only in the curvature.
   !!
   !! **Does the Krylov solve solve the Newton equations?** The step comes out
   !! of a subspace, so it is only an approximate solution by construction. The
   !! residual `H k + g` is formed through the operator and checked against the
   !! tolerance the solver was asked for.
   !!
   !! **Does it converge to the same place as DIIS?** A second-order SCF that
   !! reaches a different energy than DIIS has not converged faster, it has
   !! converged somewhere else.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t
   use mqc_scf_common, only: build_density_closed_shell
   use mqc_scf_convergence, only: scf_convergence_t, CONV_METRIC_COMMUTATOR
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, build_fock
   use mqc_czt_hessian, only: nuclear_response_t
   use mqc_czt_ov_hessian, only: ov_hessian_t, build_scf_ov_hessian
   use mqc_czt_soscf, only: soscf_gradient, soscf_semicanonicalize, &
                            soscf_newton_step, HESSIAN_SCALE
   use mqc_orbital_rotation, only: rotation_matrix, MAX_ROTATION
   implicit none
   private

   public :: collect_mqc_czt_soscf_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: STEP = 1.0e-3_dp
      !! Finite-difference displacement in the rotation angle, in radians.
      !!
      !! Central differences, so the truncation error is `O(STEP^2)` -- about
      !! 1e-6 on a first derivative. The second derivative divides the same
      !! rounding by `STEP^2`, which for an energy of order 100 hartree leaves
      !! about 1e-8, so both checks have four orders of headroom over their
      !! tolerances and neither is at the limit of the arithmetic.

contains

   subroutine collect_mqc_czt_soscf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("soscf_gradient_is_the_energy_derivative", &
                               test_gradient_by_finite_difference), &
                  new_unittest("soscf_hessian_is_the_second_derivative", &
                               test_curvature_by_finite_difference), &
                  new_unittest("soscf_semicanonicalization_moves_no_energy", &
                               test_semicanonicalization), &
                  new_unittest("soscf_step_solves_the_newton_equations", &
                               test_newton_residual), &
                  new_unittest("soscf_lands_where_diis_lands", test_same_energy_as_diis), &
                  new_unittest("soscf_refuses_a_continuum_solvent", test_refuses_pcm) &
                  ]
   end subroutine collect_mqc_czt_soscf_tests

   ! ---- the molecule and an energy at arbitrary orbitals ------------------

   subroutine water(mol, scf, err, second_order)
      !! Water in STO-3G, the geometry the stability and CPHF tests use
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err
      logical, intent(in), optional :: second_order

      real(dp) :: c(3, 3)
      type(scf_convergence_t) :: conv

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "sto-3g", mol, err)
      if (err%has_error()) return
      ! On the commutator alone, and tightly: this is what both paths are held
      ! to, so the energies they reach are comparable at the same threshold
      ! rather than at whatever `sqrt(energy_tol)` happened to derive.
      conv%metric = CONV_METRIC_COMMUTATOR
      conv%tolerance = 1.0e-9_dp
      call run_czt_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, &
                       in_core=.true., convergence=conv, second_order=second_order)
   end subroutine water

   subroutine energy_at(h, eri, coeff, n_occ, energy)
      !! The Hartree-Fock electronic energy of the determinant these orbitals build
      !!
      !! Without the nuclear repulsion, which is a constant and cancels out of
      !! every difference taken below.
      real(dp), intent(in) :: h(:, :)
      real(dp), allocatable, intent(in) :: eri(:, :, :, :)
      real(dp), intent(in) :: coeff(:, :)
      integer, intent(in) :: n_occ
      real(dp), intent(out) :: energy

      real(dp), allocatable :: density(:, :), fock(:, :)

      allocate (density(size(h, 1), size(h, 1)), fock(size(h, 1), size(h, 1)))
      call build_density_closed_shell(coeff, n_occ, density)
      call build_fock(h, eri, density, fock)
      energy = 0.5_dp*sum(density*(h + fock))
      deallocate (density, fock)
   end subroutine energy_at

   subroutine rotated(coeff, kappa, scale, moved)
      !! `C exp(scale * kappa)`
      real(dp), intent(in) :: coeff(:, :), kappa(:, :)
      real(dp), intent(in) :: scale
      real(dp), allocatable, intent(out) :: moved(:, :)

      real(dp), allocatable :: rotation(:, :)

      call rotation_matrix(scale*kappa, rotation)
      allocate (moved(size(coeff, 1), size(coeff, 2)))
      call pic_gemm(coeff, rotation, moved, beta=0.0_dp)
      deallocate (rotation)
   end subroutine rotated

   subroutine direction(n_mo, n_occ, vector, kappa)
      !! A reproducible occupied-virtual direction, flat and as a `kappa`
      !!
      !! Deterministic rather than random, so a failure fails the same way
      !! twice, and spread over every rotation rather than a single one, so a
      !! per-index mistake in the flattening cannot pass.
      integer, intent(in) :: n_mo, n_occ
      real(dp), allocatable, intent(out) :: vector(:)
      real(dp), allocatable, intent(out) :: kappa(:, :)

      integer :: n_vir, i, a
      real(dp) :: norm

      n_vir = n_mo - n_occ
      allocate (vector(n_vir*n_occ), kappa(n_mo, n_mo))
      kappa = 0.0_dp
      do i = 1, n_occ
         do a = 1, n_vir
            vector(a + (i - 1)*n_vir) = sin(0.7_dp*real(a, dp) + 1.3_dp*real(i, dp))
         end do
      end do
      norm = sqrt(dot_product(vector, vector))
      vector = vector/norm
      do i = 1, n_occ
         do a = 1, n_vir
            kappa(n_occ + a, i) = vector(a + (i - 1)*n_vir)
            kappa(i, n_occ + a) = -vector(a + (i - 1)*n_vir)
         end do
      end do
   end subroutine direction

   ! ---- the cases ---------------------------------------------------------

   subroutine test_gradient_by_finite_difference(error)
      !! `dot(g, v)` against `dE/dt` along `C exp(t v)`
      !!
      !! Taken at orbitals that have been *displaced* off the solution, so the
      !! gradient is not zero and the comparison has something in it. At the
      !! stationary point both sides are zero and any sign would pass.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: eri(:, :, :, :), h(:, :), coeff(:, :)
      real(dp), allocatable :: vector(:), kappa(:, :), offset_kappa(:, :)
      real(dp), allocatable :: fock(:, :), fock_mo(:, :), work(:, :), gradient(:)
      real(dp), allocatable :: density(:, :), plus(:, :), minus(:, :)
      real(dp) :: e_plus, e_minus, numeric, analytic
      integer :: n_ao, n_mo, n_occ

      call water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF failed")
      if (allocated(error)) return

      n_ao = mol%nao
      n_mo = size(scf%orbitals, 2)
      n_occ = scf%n_occupied
      call mol%core_hamiltonian(h)
      call mol%eris(eri)

      ! Off the solution by a tenth of a radian along a different direction
      ! than the one differentiated, so nothing is accidentally stationary.
      allocate (offset_kappa(n_mo, n_mo))
      offset_kappa = 0.0_dp
      offset_kappa(n_occ + 1, 1) = 0.1_dp
      offset_kappa(1, n_occ + 1) = -0.1_dp
      offset_kappa(n_occ + 2, 2) = -0.07_dp
      offset_kappa(2, n_occ + 2) = 0.07_dp
      call rotated(scf%orbitals, offset_kappa, 1.0_dp, coeff)

      ! The analytic gradient at those orbitals, 4 F_ai in their own basis.
      allocate (density(n_ao, n_ao), fock(n_ao, n_ao), work(n_ao, n_mo), fock_mo(n_mo, n_mo))
      call build_density_closed_shell(coeff, n_occ, density)
      call build_fock(h, eri, density, fock)
      call pic_gemm(fock, coeff, work, beta=0.0_dp)
      call pic_gemm(coeff, work, fock_mo, transa="T", beta=0.0_dp)
      call soscf_gradient(fock_mo, n_occ, gradient)

      call direction(n_mo, n_occ, vector, kappa)
      analytic = dot_product(gradient, vector)

      call rotated(coeff, kappa, STEP, plus)
      call rotated(coeff, kappa, -STEP, minus)
      call energy_at(h, eri, plus, n_occ, e_plus)
      call energy_at(h, eri, minus, n_occ, e_minus)
      numeric = (e_plus - e_minus)/(2.0_dp*STEP)

      call check(error, abs(numeric) > 1.0e-3_dp, "the displaced orbitals turned out "// &
                 "to be stationary along this direction, so the case tests nothing; "// &
                 "numeric derivative "//real_to_text(numeric))
      if (allocated(error)) return
      call check(error, abs(numeric - analytic) < 1.0e-5_dp*max(1.0_dp, abs(numeric)), &
                 "the orbital gradient is not the derivative of the energy: analytic "// &
                 real_to_text(analytic)//", finite difference "//real_to_text(numeric))

      deallocate (offset_kappa)
      call mol%destroy()
   end subroutine test_gradient_by_finite_difference

   subroutine test_curvature_by_finite_difference(error)
      !! `HESSIAN_SCALE * dot(v, (A+B) v)` against `d2E/dt2`
      !!
      !! At the converged orbitals, which is where `(A+B)` is the Hessian: the
      !! operator's first term is the orbital-energy gaps, and those are only
      !! `delta_ij F_ab - delta_ab F_ij` when the Fock matrix is diagonal in
      !! each block. This is the case that pins the factor of four.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(nuclear_response_t), target :: response
      type(ov_hessian_t) :: hessian
      real(dp), allocatable :: eri(:, :, :, :), h(:, :)
      real(dp), allocatable :: vector(:), kappa(:, :), image(:)
      real(dp), allocatable :: plus(:, :), minus(:, :)
      real(dp) :: e_plus, e_minus, e_zero, numeric, analytic
      integer :: n_mo, n_occ

      call water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF failed")
      if (allocated(error)) return

      n_mo = size(scf%orbitals, 2)
      n_occ = scf%n_occupied
      call mol%core_hamiltonian(h)
      call mol%eris(eri)

      call build_scf_ov_hessian(mol, scf%orbitals, scf%orbital_energies, n_occ, &
                                response, hessian, err)
      call check(error,.not. err%has_error(), "the Hessian would not build: "// &
                 err%get_message())
      if (allocated(error)) return

      call direction(n_mo, n_occ, vector, kappa)
      allocate (image(size(vector)))
      call hessian%apply(vector, image)
      call check(error,.not. hessian%error%has_error(), "applying the Hessian failed: "// &
                 hessian%error%get_message())
      if (allocated(error)) return
      analytic = HESSIAN_SCALE*dot_product(vector, image)

      call energy_at(h, eri, scf%orbitals, n_occ, e_zero)
      call rotated(scf%orbitals, kappa, STEP, plus)
      call rotated(scf%orbitals, kappa, -STEP, minus)
      call energy_at(h, eri, plus, n_occ, e_plus)
      call energy_at(h, eri, minus, n_occ, e_minus)
      numeric = (e_plus - 2.0_dp*e_zero + e_minus)/STEP**2

      call check(error, abs(numeric) > 1.0e-2_dp, "the curvature along this direction "// &
                 "is too small for the case to test anything; "//real_to_text(numeric))
      if (allocated(error)) return
      ! A quarter of the true curvature is what leaving `HESSIAN_SCALE` out
      ! would give, and this tolerance is four orders inside that.
      call check(error, abs(numeric - analytic) < 1.0e-5_dp*max(1.0_dp, abs(numeric)), &
                 "the electronic Hessian is not the second derivative of the energy: "// &
                 "analytic "//real_to_text(analytic)//", finite difference "// &
                 real_to_text(numeric))

      call mol%destroy()
   end subroutine test_curvature_by_finite_difference

   subroutine test_semicanonicalization(error)
      !! Rotating within the blocks changes the representation, not the state
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: eri(:, :, :, :), h(:, :), coeff(:, :)
      real(dp), allocatable :: mixing(:, :), fock(:, :), density(:, :)
      real(dp), allocatable :: fock_mo(:, :), energies(:)
      real(dp) :: before, after, worst
      integer :: n_ao, n_mo, n_occ, i, j

      call water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF failed")
      if (allocated(error)) return

      n_ao = mol%nao
      n_mo = size(scf%orbitals, 2)
      n_occ = scf%n_occupied
      call mol%core_hamiltonian(h)
      call mol%eris(eri)

      ! Mix the occupied orbitals among themselves, which a closed-shell
      ! determinant cannot notice, and the virtuals among themselves, which
      ! nothing can.
      allocate (mixing(n_mo, n_mo))
      mixing = 0.0_dp
      do i = 1, n_occ - 1
         mixing(i + 1, i) = 0.3_dp
         mixing(i, i + 1) = -0.3_dp
      end do
      do i = n_occ + 1, n_mo - 1
         mixing(i + 1, i) = 0.4_dp
         mixing(i, i + 1) = -0.4_dp
      end do
      call rotated(scf%orbitals, mixing, 1.0_dp, coeff)

      call energy_at(h, eri, coeff, n_occ, before)
      allocate (density(n_ao, n_ao), fock(n_ao, n_ao))
      call build_density_closed_shell(coeff, n_occ, density)
      call build_fock(h, eri, density, fock)
      call soscf_semicanonicalize(fock, coeff, n_occ, fock_mo, energies, err)
      call check(error,.not. err%has_error(), "semicanonicalization failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call energy_at(h, eri, coeff, n_occ, after)

      call check(error, abs(after - before) < 1.0e-12_dp, "semicanonicalization moved "// &
                 "the energy by "//real_to_text(after - before)//", so it is a step "// &
                 "and not a change of basis")
      if (allocated(error)) return

      ! And the blocks really are diagonal now, which is the whole point.
      worst = 0.0_dp
      do j = 1, n_occ
         do i = 1, n_occ
            if (i /= j) worst = max(worst, abs(fock_mo(i, j)))
         end do
      end do
      do j = n_occ + 1, n_mo
         do i = n_occ + 1, n_mo
            if (i /= j) worst = max(worst, abs(fock_mo(i, j)))
         end do
      end do
      call check(error, worst < 1.0e-10_dp, "the Fock matrix is not block diagonal "// &
                 "after semicanonicalization; worst element "//real_to_text(worst))
      if (allocated(error)) return
      call check(error, all(abs(energies - [(fock_mo(i, i), i=1, n_mo)]) < 1.0e-14_dp), &
                 "the orbital energies are not the diagonal they are meant to be")

      call mol%destroy()
   end subroutine test_semicanonicalization

   subroutine test_newton_residual(error)
      !! `H k + g` is as small as the subspace solve promised
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(nuclear_response_t), target :: response
      type(ov_hessian_t) :: hessian
      real(dp), allocatable :: eri(:, :, :, :), h(:, :), coeff(:, :)
      real(dp), allocatable :: offset_kappa(:, :), fock(:, :), fock_mo(:, :)
      real(dp), allocatable :: energies(:), gradient(:), kappa(:, :)
      real(dp), allocatable :: flat(:), image(:), density(:, :)
      real(dp) :: lowest, predicted, residual, gnorm
      integer :: n_ao, n_mo, n_occ, n_vir, products, i, a

      call water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF failed")
      if (allocated(error)) return

      n_ao = mol%nao
      n_mo = size(scf%orbitals, 2)
      n_occ = scf%n_occupied
      n_vir = n_mo - n_occ
      call mol%core_hamiltonian(h)
      call mol%eris(eri)

      ! Displaced, so there is a gradient for the Newton equations to solve
      ! against, and then semicanonicalised, which is what the SCF does before
      ! it differentiates.
      allocate (offset_kappa(n_mo, n_mo))
      offset_kappa = 0.0_dp
      offset_kappa(n_occ + 1, 1) = 0.05_dp
      offset_kappa(1, n_occ + 1) = -0.05_dp
      call rotated(scf%orbitals, offset_kappa, 1.0_dp, coeff)

      allocate (density(n_ao, n_ao), fock(n_ao, n_ao))
      call build_density_closed_shell(coeff, n_occ, density)
      call build_fock(h, eri, density, fock)
      call soscf_semicanonicalize(fock, coeff, n_occ, fock_mo, energies, err)
      call check(error,.not. err%has_error(), "semicanonicalization failed")
      if (allocated(error)) return
      call soscf_gradient(fock_mo, n_occ, gradient)
      gnorm = sqrt(dot_product(gradient, gradient))
      call check(error, gnorm > 1.0e-4_dp, "there is no gradient to solve against")
      if (allocated(error)) return

      call build_scf_ov_hessian(mol, coeff, energies, n_occ, response, hessian, err)
      call check(error,.not. err%has_error(), "the Hessian would not build: "// &
                 err%get_message())
      if (allocated(error)) return

      call soscf_newton_step(hessian, gradient, MAX_ROTATION, kappa, lowest, &
                             predicted, products, err)
      call check(error,.not. err%has_error(), "the Newton step failed: "//err%get_message())
      if (allocated(error)) return
      call check(error, products > 0, "the Newton step reported no Hessian-vector "// &
                 "products, so it never reached the operator")
      if (allocated(error)) return
      call check(error, lowest > 0.0_dp, "water in STO-3G near its solution should "// &
                 "have positive curvature; the step reported "//real_to_text(lowest))
      if (allocated(error)) return
      call check(error, predicted > 0.0_dp, "a Newton step from a non-zero gradient "// &
                 "should predict a decrease, not "//real_to_text(predicted))
      if (allocated(error)) return

      ! `kappa` back to a flat vector, and the residual of the unshifted Newton
      ! equations through the operator. No level shift is active here -- the
      ! curvature was just checked positive and well above `MIN_CURVATURE` --
      ! so this is the equation the solver was actually asked to solve.
      allocate (flat(n_vir*n_occ), image(n_vir*n_occ))
      do i = 1, n_occ
         do a = 1, n_vir
            flat(a + (i - 1)*n_vir) = kappa(n_occ + a, i)
         end do
      end do
      call hessian%apply(flat, image)
      call check(error,.not. hessian%error%has_error(), "applying the Hessian failed")
      if (allocated(error)) return
      image = HESSIAN_SCALE*image + gradient
      residual = sqrt(dot_product(image, image))

      call check(error, residual < 1.0e-2_dp*gnorm, "the Newton step does not solve "// &
                 "the Newton equations: residual "//real_to_text(residual)// &
                 " against a gradient of "//real_to_text(gnorm))
      if (allocated(error)) return

      ! And the step descends: to first order the energy change is dot(g, k).
      call check(error, dot_product(gradient, flat) < 0.0_dp, "the Newton step points "// &
                 "uphill, which is the sign convention being wrong")

      call mol%destroy()
   end subroutine test_newton_residual

   subroutine test_same_energy_as_diis(error)
      !! Two ways to the same stationary point
      !!
      !! Both SCFs stop on the same commutator threshold, so the energies are
      !! comparable at the same tightness. A disagreement here is two different
      !! solutions, not a slower method.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol_diis, mol_soscf
      type(rhf_result_t) :: diis, soscf
      type(error_t) :: err
      real(dp) :: difference

      call water(mol_diis, diis, err)
      call check(error,.not. err%has_error() .and. diis%converged, "the DIIS SCF failed")
      if (allocated(error)) return

      call water(mol_soscf, soscf, err, second_order=.true.)
      call check(error,.not. err%has_error(), "the second-order SCF errored: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, soscf%converged, "the second-order SCF did not converge")
      if (allocated(error)) return
      call check(error, soscf%second_order_started_at > 0, "the second-order SCF never "// &
                 "handed over from DIIS, so this case only ran DIIS twice")
      if (allocated(error)) return
      call check(error, soscf%second_order_iterations > 0, "the second-order phase "// &
                 "took no Newton steps")
      if (allocated(error)) return
      call check(error, soscf%second_order_hessian_products > 0, "the Newton steps "// &
                 "spent no Hessian-vector products")
      if (allocated(error)) return

      ! Well inside the threshold both were converged to. The energy error of
      ! an SCF goes as the square of the commutator, so at 1e-9 there is no
      ! room for a genuine difference above 1e-10.
      difference = abs(soscf%energy - diis%energy)
      call check(error, difference < 1.0e-10_dp, "the second-order SCF converged to "// &
                 real_to_text(soscf%energy)//" and DIIS to "//real_to_text(diis%energy)// &
                 ", a difference of "//real_to_text(difference)//" -- two different "// &
                 "solutions, not two routes to one")
      if (allocated(error)) return

      ! The curvature the last step saw, reported whether or not it was needed.
      call check(error, soscf%has_second_order_curvature, "the second-order phase did "// &
                 "not report a curvature")
      if (allocated(error)) return
      call check(error, soscf%second_order_curvature > 0.0_dp, "a converged closed "// &
                 "shell should not report negative curvature; it said "// &
                 real_to_text(soscf%second_order_curvature))

      call mol_diis%destroy()
      call mol_soscf%destroy()
   end subroutine test_same_energy_as_diis

   subroutine test_refuses_pcm(error)
      !! What it cannot do it declines, naming what is missing
      !!
      !! A continuum solvent adds a term to the energy whose curvature is not
      !! in the orbital-rotation Hessian, so a Newton step would be taken on
      !! the wrong surface. Refused rather than approximated.
      type(error_type), allocatable, intent(out) :: error

      use_pcm: block
         use mqc_czt_pcm, only: pcm_context_t
         type(czt_molecule_t) :: mol
         type(rhf_result_t) :: scf
         type(error_t) :: err
         type(pcm_context_t) :: pcm
         real(dp) :: c(3, 3)

         c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                      0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                      0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
         call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "sto-3g", mol, err)
         call check(error,.not. err%has_error(), "the molecule would not build")
         if (allocated(error)) return

         pcm%enabled = .true.
         call run_czt_rhf(mol, 10, 20, 1.0e-8_dp, 1.0e-6_dp, .false., scf, err, &
                          in_core=.true., pcm=pcm, second_order=.true.)
         call check(error, err%has_error(), "a second-order SCF in a continuum "// &
                    "solvent should be refused")
         if (allocated(error)) return
         call check(error, index(err%get_message(), "second_order") > 0, &
                    "the refusal should name the keyword; it said: "//err%get_message())
         call mol%destroy()
      end block use_pcm
   end subroutine test_refuses_pcm

   function real_to_text(value) result(text)
      !! A number in a failure message, so the message says how badly
      real(dp), intent(in) :: value
      character(len=:), allocatable :: text

      character(len=32) :: buffer

      write (buffer, "(es16.8)") value
      text = trim(adjustl(buffer))
   end function real_to_text

end module test_mqc_czt_soscf

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_soscf, only: collect_mqc_czt_soscf_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_soscf", collect_mqc_czt_soscf_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
