!! Unit tests for the SCF linear algebra shared by the libcint and cuEST paths
module test_mqc_scf_common
   !! These routines used to exist twice, once per backend. Now that one copy
   !! serves both, its properties are worth pinning: an orthogonaliser that
   !! actually orthogonalises, a density that is idempotent under the overlap,
   !! and a spin contamination that is zero for a closed shell.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_scf_common, only: build_orthogonalizer, build_density_closed_shell, &
                             build_density_spin, spin_contamination, &
                             LINEAR_DEPENDENCE_TOL
   implicit none
   private

   public :: collect_mqc_scf_common

   real(dp), parameter :: TOL = 1.0e-10_dp

contains

   subroutine collect_mqc_scf_common(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("scf_orthogonalizer_diagonalizes_the_overlap", test_orth_identity), &
                  new_unittest("scf_orthogonalizer_keeps_every_mode_when_well_conditioned", test_orth_full_rank), &
                  new_unittest("scf_orthogonalizer_drops_a_null_mode", test_orth_drops_null), &
                  new_unittest("scf_orthogonalizer_refuses_a_singular_overlap", test_orth_singular), &
                  new_unittest("scf_closed_shell_density_is_idempotent", test_density_idempotent), &
                  new_unittest("scf_closed_shell_density_traces_to_the_electron_count", test_density_trace), &
                  new_unittest("scf_spin_density_is_half_the_closed_shell_one", test_density_spin), &
                  new_unittest("scf_density_of_no_electrons_is_zero", test_density_empty), &
                  new_unittest("scf_spin_contamination_vanishes_for_a_closed_shell", test_s2_closed), &
                  new_unittest("scf_spin_contamination_is_exact_for_a_lone_electron", test_s2_doublet) &
                  ]
   end subroutine collect_mqc_scf_common

   pure subroutine model_overlap(n, off_diagonal, overlap)
      !! A symmetric positive-definite overlap: 1 on the diagonal, `s` off it
      integer, intent(in) :: n
      real(dp), intent(in) :: off_diagonal
      real(dp), allocatable, intent(out) :: overlap(:, :)
      integer :: i, j

      allocate (overlap(n, n))
      do j = 1, n
         do i = 1, n
            if (i == j) then
               overlap(i, j) = 1.0_dp
            else
               overlap(i, j) = off_diagonal
            end if
         end do
      end do
   end subroutine model_overlap

   subroutine test_orth_identity(error)
      !! X^T S X = I, which is the entire contract of an orthogonaliser
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: overlap(:, :), x(:, :), product(:, :)
      type(error_t) :: err
      integer :: n_mo, i, j

      call model_overlap(4, 0.25_dp, overlap)
      call build_orthogonalizer(overlap, x, n_mo, err)
      call check(error,.not. err%has_error(), "a well-conditioned overlap must succeed")
      if (allocated(error)) return

      product = matmul(transpose(x), matmul(overlap, x))
      do j = 1, n_mo
         do i = 1, n_mo
            if (i == j) then
               call check(error, product(i, j), 1.0_dp, thr=TOL)
            else
               call check(error, product(i, j), 0.0_dp, thr=TOL)
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine test_orth_identity

   subroutine test_orth_full_rank(error)
      !! Nothing is discarded when no eigenvalue is near zero
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: overlap(:, :), x(:, :)
      type(error_t) :: err
      integer :: n_mo

      call model_overlap(5, 0.1_dp, overlap)
      call build_orthogonalizer(overlap, x, n_mo, err)
      call check(error,.not. err%has_error(), "a well-conditioned overlap must succeed")
      if (allocated(error)) return
      call check(error, n_mo, 5)
      if (allocated(error)) return
      call check(error, size(x, 2), 5)
   end subroutine test_orth_full_rank

   subroutine test_orth_drops_null(error)
      !! A duplicated basis function is a null mode, and must be dropped
      !!
      !! Two identical rows/columns make the overlap exactly rank-deficient,
      !! which is the near-linear-dependence case the canonical orthogonaliser
      !! exists to survive.
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: overlap(:, :), x(:, :)
      type(error_t) :: err
      integer :: n_mo

      allocate (overlap(3, 3))
      overlap = reshape([1.0_dp, 0.0_dp, 0.0_dp, &
                         0.0_dp, 1.0_dp, 1.0_dp, &
                         0.0_dp, 1.0_dp, 1.0_dp], [3, 3])
      call build_orthogonalizer(overlap, x, n_mo, err)
      call check(error,.not. err%has_error(), "a rank-deficient overlap is survivable")
      if (allocated(error)) return
      call check(error, n_mo, 2)
   end subroutine test_orth_drops_null

   subroutine test_orth_singular(error)
      !! An all-zero overlap has no surviving modes and must be reported
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: overlap(:, :), x(:, :)
      type(error_t) :: err
      integer :: n_mo

      allocate (overlap(3, 3))
      overlap = 0.0_dp
      call build_orthogonalizer(overlap, x, n_mo, err)
      call check(error, err%has_error(), "a singular overlap must be refused, not returned")
   end subroutine test_orth_singular

   subroutine test_density_idempotent(error)
      !! D S D = 2 D for the closed-shell density, the standard identity
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: overlap(:, :), x(:, :), density(:, :), check_matrix(:, :)
      type(error_t) :: err
      integer :: n_mo, i, j

      call model_overlap(4, 0.2_dp, overlap)
      call build_orthogonalizer(overlap, x, n_mo, err)
      call check(error,.not. err%has_error(), "orthogonaliser must succeed")
      if (allocated(error)) return

      ! Orthonormal columns of X are legitimate orbitals for this overlap.
      allocate (density(4, 4))
      call build_density_closed_shell(x, 2, density)

      check_matrix = matmul(density, matmul(overlap, density))
      do j = 1, 4
         do i = 1, 4
            call check(error, check_matrix(i, j), 2.0_dp*density(i, j), thr=1.0e-9_dp)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_density_idempotent

   subroutine test_density_trace(error)
      !! tr(D S) is the electron count: two per occupied orbital
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: overlap(:, :), x(:, :), density(:, :), ds(:, :)
      type(error_t) :: err
      integer :: n_mo, i
      real(dp) :: trace

      call model_overlap(4, 0.15_dp, overlap)
      call build_orthogonalizer(overlap, x, n_mo, err)
      call check(error,.not. err%has_error(), "orthogonaliser must succeed")
      if (allocated(error)) return

      allocate (density(4, 4))
      call build_density_closed_shell(x, 3, density)
      ds = matmul(density, overlap)
      trace = 0.0_dp
      do i = 1, 4
         trace = trace + ds(i, i)
      end do
      call check(error, trace, 6.0_dp, thr=1.0e-9_dp)
   end subroutine test_density_trace

   subroutine test_density_spin(error)
      !! The spin density carries one electron per orbital, not two
      !!
      !! The two routines are deliberately separate so that a spin density is
      !! never built with the closed-shell factor; this is that guarantee.
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: coeff(:, :), closed(:, :), spin(:, :)
      integer :: i, j

      allocate (coeff(3, 3), closed(3, 3), spin(3, 3))
      coeff = reshape([1.0_dp, 0.0_dp, 0.0_dp, &
                       0.0_dp, 1.0_dp, 0.0_dp, &
                       0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
      call build_density_closed_shell(coeff, 2, closed)
      call build_density_spin(coeff, 2, spin)

      do j = 1, 3
         do i = 1, 3
            call check(error, closed(i, j), 2.0_dp*spin(i, j), thr=TOL)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_density_spin

   subroutine test_density_empty(error)
      !! No occupied orbitals means a zero density, not stale memory
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: coeff(:, :), density(:, :)

      allocate (coeff(3, 3), density(3, 3))
      coeff = 1.0_dp
      density = 99.0_dp
      call build_density_closed_shell(coeff, 0, density)
      call check(error, maxval(abs(density)), 0.0_dp, thr=TOL)
   end subroutine test_density_empty

   subroutine test_s2_closed(error)
      !! Identical alpha and beta orbitals: <S^2> = 0 exactly
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: coeff(:, :), overlap(:, :)
      real(dp) :: s2

      allocate (coeff(3, 3))
      coeff = reshape([1.0_dp, 0.0_dp, 0.0_dp, &
                       0.0_dp, 1.0_dp, 0.0_dp, &
                       0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
      call model_overlap(3, 0.0_dp, overlap)

      s2 = spin_contamination(coeff, coeff, overlap, 2, 2)
      call check(error, s2, 0.0_dp, thr=TOL)
   end subroutine test_s2_closed

   subroutine test_s2_doublet(error)
      !! One unpaired electron in an orthogonal set: <S^2> = 0.75 exactly
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: coeff(:, :), overlap(:, :)
      real(dp) :: s2

      allocate (coeff(3, 3))
      coeff = reshape([1.0_dp, 0.0_dp, 0.0_dp, &
                       0.0_dp, 1.0_dp, 0.0_dp, &
                       0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
      call model_overlap(3, 0.0_dp, overlap)

      ! Two alpha, one beta: S_z = 1/2, and the beta orbital is one of the
      ! alpha ones, so the overlap sum removes the contamination entirely.
      s2 = spin_contamination(coeff, coeff, overlap, 2, 1)
      call check(error, s2, 0.75_dp, thr=TOL)
      if (allocated(error)) return

      ! With no beta electrons at all the sum is empty and <S^2> = S_z(S_z+1).
      s2 = spin_contamination(coeff, coeff, overlap, 2, 0)
      call check(error, s2, 2.0_dp, thr=TOL)
   end subroutine test_s2_doublet

end module test_mqc_scf_common

program tester_mqc_scf_common
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_scf_common, only: collect_mqc_scf_common
   implicit none

   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_scf_common", collect_mqc_scf_common)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_scf_common
