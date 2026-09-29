!! The matrix-free SA-CASSCF Hessian-vector product, a genuinely
!! multi-dimensional CI space
module test_mqc_sa_hessian_multistate
   !! `test_mqc_sa_hessian.f90`'s LiH/STO-3G SA-2-CAS(2,2) case leaves only
   !! one non-redundant CI direction per state (three singlet determinants,
   !! minus the two projected reference states), so its gates B/C/D barely
   !! exercise the CI-CI and orbital-CI blocks. This file uses LiH/6-31G,
   !! CAS(4,4), SA-3: 11 basis functions, no frozen core (4 active electrons
   !! is every electron LiH has, so `n_inactive = 0` and the "active" space
   !! includes the Li 1s-like orbital -- not a physically sensible active
   !! space, but a fine one for exercising the linear algebra), 6 alpha and
   !! 6 beta strings, 36 determinants, 21 of them singlet-symmetric, minus
   !! the 3 reference states leaves 18 genuine non-redundant CI directions
   !! per state.
   !!
   !! Also covers unequal weights (0.7, 0.3): `sa_hessian_apply`'s own
   !! docstring states the redundancy projection is exact only for equal
   !! weights, so this checks what the operator actually does in that case
   !! -- not that PySCF's `NotImplementedError` is reproduced (phase 4's
   !! problem), but that the *projected* Hessian still matches a finite
   !! difference of the *projected* gradient, consistently gauge-fixed on
   !! both sides exactly as the equal-weight gate does.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mcscf, only: run_czt_casscf, casscf_result_t, rotation_matrix
   use mqc_czt_sa_hessian, only: sa_hessian_t, build_sa_hessian, destroy_sa_hessian, &
                                 sa_hessian_n_param, sa_gradient, sa_hessian_apply, &
                                 project_ci_block
   implicit none
   private

   public :: collect_mqc_sa_hessian_multistate_tests

   real(dp), parameter :: LIH(3, 2) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 0.0_dp, 3.0139241961656_dp], [3, 2])
   integer, parameter :: LIH_Z(2) = [3, 1]
   character(len=2), parameter :: LIH_SYM(2) = ["Li", "H "]

contains

   subroutine collect_mqc_sa_hessian_multistate_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("hvp_symmetric_cas44_sa3", test_symmetry_sa3), &
                  new_unittest("hvp_against_fd_cas44_sa3", test_finite_difference_sa3), &
                  new_unittest("hvp_against_fd_cas44_unequal_weights", &
                               test_finite_difference_unequal_weights) &
                  ]
   end subroutine collect_mqc_sa_hessian_multistate_tests

   subroutine lih_631g_reference(mol, orbitals, err)
      !! LiH/6-31G, converged RHF orbitals (11 basis functions)
      type(czt_molecule_t), intent(out) :: mol
      real(dp), allocatable, intent(out) :: orbitals(:, :)
      type(error_t), intent(inout) :: err

      type(rhf_result_t) :: scf

      call build_czt_molecule(LIH_Z, LIH_SYM, LIH, "6-31g", mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, 4, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      if (err%has_error()) return
      if (.not. scf%converged) then
         call err%set(ERROR_VALIDATION, "the LiH/6-31G RHF reference did not converge")
         return
      end if
      orbitals = scf%orbitals
   end subroutine lih_631g_reference

   subroutine converged_cas44(mol, orbitals, weights, result, err, ok)
      !! CAS(4,4), no frozen core, `size(weights)` singlet states
      type(czt_molecule_t), intent(out) :: mol
      real(dp), allocatable, intent(out) :: orbitals(:, :)
      real(dp), intent(in) :: weights(:)
      type(casscf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      ok = .false.
      call lih_631g_reference(mol, orbitals, err)
      if (err%has_error()) return
      call run_czt_casscf(mol, orbitals, 0, 4, 2, 2, result, err, max_iterations=500, &
                          gradient_tol=1.0e-10_dp, n_states=size(weights), weights=weights)
      if (err%has_error()) return
      ok = result%converged
   end subroutine converged_cas44

   function flat_kappa_direction(n_rot, phase) result(v)
      integer, intent(in) :: n_rot
      real(dp), intent(in) :: phase
      real(dp) :: v(n_rot)

      integer :: l

      do l = 1, n_rot
         v(l) = sin(0.31_dp*real(l, dp) + phase)
      end do
   end function flat_kappa_direction

   subroutine ci_direction(na, nb, phase, v)
      integer, intent(in) :: na, nb
      real(dp), intent(in) :: phase
      real(dp), intent(out) :: v(na, nb)

      integer :: ia, ib

      do ib = 1, nb
         do ia = 1, na
            v(ia, ib) = sin(0.41_dp*real(ia, dp) + 0.23_dp*real(ib, dp) + phase)
         end do
      end do
   end subroutine ci_direction

   subroutine test_symmetry_sa3(error)
      !! |<y,Hx> - <x,Hy>| on the CAS(4,4) SA-3 system, random directions
      !! mixing the orbital and all three CI blocks
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), allocatable :: x(:, :), y(:, :), hx(:, :), hy(:, :)
      real(dp), allocatable :: kdir(:), vci(:, :)
      real(dp) :: lhs, rhs, scale
      real(dp), parameter :: WEIGHTS(3) = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]
      integer :: np, n_rot, na, nb, j
      logical :: ok

      call converged_cas44(mol, orbitals, WEIGHTS, result, err, ok)
      call check(error, ok, "CAS(4,4) SA-3 should converge")
      if (allocated(error)) return
      call build_sa_hessian(mol, result%orbitals, 0, 4, 2, 2, result%ci_vectors, &
                            result%energies, WEIGHTS, state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      n_rot = state%n_rot
      na = state%alpha%n_strings
      nb = state%beta%n_strings
      np = sa_hessian_n_param(state)
      allocate (x(np, 1), y(np, 1), hx(np, 1), hy(np, 1), vci(na, nb))

      kdir = flat_kappa_direction(n_rot, 0.15_dp)
      x(1:n_rot, 1) = kdir
      do j = 1, 3
         call ci_direction(na, nb, 0.3_dp + 0.1_dp*real(j, dp), vci)
         call project_ci_block(state, vci)
         vci = 0.5_dp*(vci + transpose(vci))
         x(n_rot + (j - 1)*state%n_det + 1:n_rot + j*state%n_det, 1) = reshape(vci, [na*nb])
      end do

      kdir = flat_kappa_direction(n_rot, 1.7_dp)
      y(1:n_rot, 1) = kdir
      do j = 1, 3
         call ci_direction(na, nb, 1.2_dp + 0.2_dp*real(j, dp), vci)
         call project_ci_block(state, vci)
         vci = 0.5_dp*(vci + transpose(vci))
         y(n_rot + (j - 1)*state%n_det + 1:n_rot + j*state%n_det, 1) = reshape(vci, [na*nb])
      end do

      call sa_hessian_apply(state, x, hx, err)
      call check(error,.not. err%has_error(), "sa_hessian_apply(x) should not error")
      if (allocated(error)) return
      call sa_hessian_apply(state, y, hy, err)
      call check(error,.not. err%has_error(), "sa_hessian_apply(y) should not error")
      if (allocated(error)) return

      lhs = sum(y(:, 1)*hx(:, 1))
      rhs = sum(x(:, 1)*hy(:, 1))
      scale = max(1.0_dp, abs(lhs), abs(rhs))
      call check(error, abs(lhs - rhs) < 1.0e-9_dp*scale, "<y,Hx> should equal <x,Hy>")

      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_symmetry_sa3

   subroutine fd_gate(state, mol, weights, error, label_suffix)
      !! The shared FD-gate body: orbital-only, CI-only and mixed directions,
      !! each CI block ( one per state) checked and reported separately
      type(sa_hessian_t), intent(in) :: state
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: weights(:)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label_suffix

      type(error_t) :: err
      real(dp), allocatable :: kappa_dir_flat(:), kappa_dir(:, :)
      real(dp), allocatable :: xci_dir(:, :, :)
      real(dp), allocatable :: x(:, :), hx(:, :)
      real(dp), allocatable :: g_plus_orb(:, :), g_minus_orb(:, :)
      real(dp), allocatable :: g_plus_ci(:, :, :), g_minus_ci(:, :, :)
      real(dp), allocatable :: fd_orb(:), fd_ci_j(:, :)
      real(dp) :: h_step, denom
      integer :: np, n_rot, na, nb, j, kase, n_states
      character(len=200) :: msg

      n_states = size(weights)
      n_rot = state%n_rot
      na = state%alpha%n_strings
      nb = state%beta%n_strings
      np = sa_hessian_n_param(state)
      ! Smaller than the LiH/STO-3G SA-2 case's 1e-4: this system's CI-CI
      ! block runs several Hartree, not several hundredths, so the same 1e-4
      ! step leaves an absolute central-difference truncation error of order
      ! 1e-6 -- fine in absolute terms but ~2e-7 relative to a ~4-9 Hartree
      ! quantity, over the 1e-7 gate below for no physical reason. Confirmed
      ! by re-running at 1e-4 (~1.7e-7 relative) and here (~5e-8 relative):
      ! the error tracks h^2, so it is truncation, not a wrong formula.
      h_step = 2.0e-5_dp

      allocate (kappa_dir(state%n_mo, state%n_mo), xci_dir(na, nb, n_states))
      allocate (x(np, 1), hx(np, 1))
      allocate (g_plus_ci(na, nb, n_states), g_minus_ci(na, nb, n_states))
      allocate (fd_ci_j(na, nb))

      do kase = 1, 3
         kappa_dir_flat = flat_kappa_direction(n_rot, 0.6_dp)
         kappa_dir = 0.0_dp
         if (kase /= 2) then
            block
               integer :: l
               do l = 1, n_rot
                  kappa_dir(state%rows(l), state%cols(l)) = kappa_dir_flat(l)
                  kappa_dir(state%cols(l), state%rows(l)) = -kappa_dir_flat(l)
               end do
            end block
         else
            kappa_dir_flat = 0.0_dp
         end if

         xci_dir = 0.0_dp
         if (kase /= 1) then
            do j = 1, n_states
               call ci_direction(na, nb, 0.8_dp + 0.2_dp*real(j, dp), xci_dir(:, :, j))
               call project_ci_block(state, xci_dir(:, :, j))
               xci_dir(:, :, j) = 0.5_dp*(xci_dir(:, :, j) + transpose(xci_dir(:, :, j)))
            end do
         end if

         x = 0.0_dp
         x(1:n_rot, 1) = kappa_dir_flat
         do j = 1, n_states
            x(n_rot + (j - 1)*state%n_det + 1:n_rot + j*state%n_det, 1) = &
               reshape(xci_dir(:, :, j), [na*nb])
         end do
         call sa_hessian_apply(state, x, hx, err)
         call check(error,.not. err%has_error(), "sa_hessian_apply should not error")
         if (allocated(error)) exit

         call sa_gradient(state, h_step*kappa_dir, h_step*xci_dir, g_plus_orb, g_plus_ci, err)
         call check(error,.not. err%has_error(), "sa_gradient(+h) should not error")
         if (allocated(error)) exit
         call sa_gradient(state, -h_step*kappa_dir, -h_step*xci_dir, g_minus_orb, g_minus_ci, err)
         call check(error,.not. err%has_error(), "sa_gradient(-h) should not error")
         if (allocated(error)) exit

         allocate (fd_orb(n_rot))
         block
            integer :: l
            do l = 1, n_rot
               fd_orb(l) = (g_plus_orb(state%rows(l), state%cols(l)) &
                            - g_minus_orb(state%rows(l), state%cols(l)))/(2.0_dp*h_step)
            end do
         end block
         denom = max(1.0_dp, maxval(abs(hx(1:n_rot, 1))))
         write (msg, "(a,i0,a,a)") "kase ", kase, ": the orbital output should match the "// &
            "finite difference", label_suffix
         call check(error, maxval(abs(hx(1:n_rot, 1) - fd_orb)) < 1.0e-7_dp*denom, trim(msg))
         deallocate (fd_orb)
         if (allocated(error)) exit

         do j = 1, n_states
            fd_ci_j = (g_plus_ci(:, :, j) - g_minus_ci(:, :, j))/(2.0_dp*h_step)
            ! Gauge-consistent comparison: see `sa_hessian_apply`'s own
            ! docstring and `test_mqc_sa_hessian.f90`'s finite-difference
            ! gate for why the FD side is projected too.
            call project_ci_block(state, fd_ci_j)
            denom = max(1.0_dp, maxval(abs(hx(n_rot + (j - 1)*state%n_det + 1: &
                                              n_rot + j*state%n_det, 1))))
            write (msg, "(a,i0,a,i0,a,a)") "kase ", kase, " state ", j, ": the CI output "// &
               "should match the finite difference", label_suffix
            call check(error, maxval(abs(reshape(hx(n_rot + (j - 1)*state%n_det + 1: &
                                                    n_rot + j*state%n_det, 1), [na, nb]) &
                                         - fd_ci_j)) < 1.0e-7_dp*denom, trim(msg))
            if (allocated(error)) exit
         end do
         if (allocated(error)) exit
      end do
   end subroutine fd_gate

   subroutine test_finite_difference_sa3(error)
      !! Gate B on CAS(4,4) SA-3: orbital-only, CI-only and mixed directions,
      !! every one of the three CI blocks checked and reported separately
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), parameter :: WEIGHTS(3) = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]
      logical :: ok

      call converged_cas44(mol, orbitals, WEIGHTS, result, err, ok)
      call check(error, ok, "CAS(4,4) SA-3 should converge")
      if (allocated(error)) return
      call build_sa_hessian(mol, result%orbitals, 0, 4, 2, 2, result%ci_vectors, &
                            result%energies, WEIGHTS, state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      call fd_gate(state, mol, WEIGHTS, error, " (equal weights)")

      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_finite_difference_sa3

   subroutine test_finite_difference_unequal_weights(error)
      !! Gate B with weights (0.7, 0.3): the in-space state-mixing rotation
      !! is no longer redundant (curvature `(w_J-w_K)(E_K-E_J)`), but the
      !! *projected* operator this module builds should still match a
      !! *projected* finite difference, consistently gauge-fixed on both
      !! sides -- see `sa_hessian_apply`'s docstring.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), parameter :: WEIGHTS(2) = [0.7_dp, 0.3_dp]
      logical :: ok

      call converged_cas44(mol, orbitals, WEIGHTS, result, err, ok)
      call check(error, ok, "CAS(4,4) SA-2 (0.7, 0.3) should converge")
      if (allocated(error)) return
      call build_sa_hessian(mol, result%orbitals, 0, 4, 2, 2, result%ci_vectors, &
                            result%energies, WEIGHTS, state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      call fd_gate(state, mol, WEIGHTS, error, " (weights 0.7/0.3)")

      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_finite_difference_unequal_weights

end module test_mqc_sa_hessian_multistate

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_sa_hessian_multistate, only: collect_mqc_sa_hessian_multistate_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_sa_hessian_multistate", &
                               collect_mqc_sa_hessian_multistate_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
