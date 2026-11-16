!! The Coulomb-only density fit against the dense fitted tensor
module test_mqc_czt_jfit
   !! `build_df_jfit` keeps the raw `(mu nu|P)` for the kept shell pairs and the
   !! metric's factor; `build_df_tensor` forms `B = (mn|Q) J^-1/2` over the full
   !! square. They are two evaluations of one fitted Coulomb matrix, so they
   !! have to agree to rounding, and an SCF built on either has to land on the
   !! same energy.
   !!
   !! Three questions:
   !!
   !! **Is J the same, with the screen removing nothing?** A compact cluster,
   !! compared against `B (B^T D)`.
   !!
   !! **And with it removing pairs?** Two waters far enough apart that the
   !! screen drops the cross pairs. The screened packed fit is compared with the
   !! screened dense one, and the unscreened packed fit with a dense reference
   !! built from the metric's inverse square root and no screen at all.
   !!
   !! **Does an RKS SCF land on the same energy?** The switch that forces the
   !! dense tensor is a testing one and is released again afterwards.
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, build_df_tensor, &
                                df_jfit_t, build_df_jfit, df_jfit_coulomb, &
                                three_centre, two_centre, metric_inverse_sqrt
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, set_df_jfit_enabled
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available
   implicit none
   private

   public :: collect_mqc_czt_jfit_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: WATER(3, 3) = reshape([ &
                                                0.0_dp, 0.0_dp, 0.1173_dp, &
                                                0.0_dp, 0.7572_dp, -0.4692_dp, &
                                                0.0_dp, -0.7572_dp, -0.4692_dp], [3, 3])*ANG

contains

   subroutine collect_mqc_czt_jfit_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("jfit_matches_dense_in_a_cluster", test_cluster), &
                  new_unittest("jfit_matches_dense_when_the_screen_bites", test_screened), &
                  new_unittest("rks_energy_is_the_same_on_both_paths", test_rks_energy) &
                  ]
   end subroutine collect_mqc_czt_jfit_tests

   subroutine place(n_water, shifts, orb, aux, error)
      !! `n_water` waters at `shifts(:, i)` bohr from the first, in def2-SVP and its J fit set
      integer, intent(in) :: n_water
      real(dp), intent(in) :: shifts(3, n_water)
      type(czt_molecule_t), intent(out) :: orb, aux
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: c(:, :)
      integer, allocatable :: z(:)
      character(len=2), allocatable :: sym(:)
      integer :: w, a

      allocate (c(3, 3*n_water), z(3*n_water), sym(3*n_water))
      do w = 1, n_water
         do a = 1, 3
            c(:, 3*(w - 1) + a) = WATER(:, a) + shifts(:, w)
            z(3*(w - 1) + a) = merge(8, 1, a == 1)
            sym(3*(w - 1) + a) = merge("O ", "H ", a == 1)
         end do
      end do
      call build_czt_molecule(z, sym, c, "def2-svp", orb, error)
      if (error%has_error()) return
      call build_czt_molecule(z, sym, c, "def2-universal-jfit", aux, error)
   end subroutine place

   subroutine random_density(n, d)
      !! A symmetric matrix of order one, fixed seed so a failure repeats
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: d(:, :)

      integer :: i, j
      real(dp) :: x

      allocate (d(n, n))
      x = 0.37_dp
      do j = 1, n
         do i = 1, j
            x = mod(x*9301.0_dp + 49297.0_dp, 233280.0_dp)
            d(i, j) = x/233280.0_dp - 0.5_dp
            d(j, i) = d(i, j)
         end do
      end do
   end subroutine random_density

   function relative_difference(a, b) result(r)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp) :: r

      r = maxval(abs(a - b))/max(maxval(abs(b)), tiny(1.0_dp))
   end function relative_difference

   subroutine dense_j(b, d, j_mat)
      !! `B (B^T D)`, the fitted Coulomb matrix from the dense tensor
      real(dp), intent(in) :: b(:, :), d(:, :)
      real(dp), allocatable, intent(out) :: j_mat(:, :)

      real(dp), allocatable :: c(:), flat(:)
      integer :: n

      n = size(d, 1)
      flat = reshape(d, [n*n])
      allocate (c(size(b, 2)))
      c = matmul(flat, b)
      flat = matmul(b, c)
      j_mat = reshape(flat, [n, n])
   end subroutine dense_j

   subroutine test_cluster(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: orb, aux
      type(df_jfit_t) :: jfit
      type(error_t) :: err
      real(dp), allocatable :: b(:, :), d(:, :), j_dense(:, :), j_fit(:, :)
      real(dp) :: shifts(3, 3)
      real(dp) :: diff

      shifts = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                        5.5_dp, 0.0_dp, 0.8_dp, &
                        0.6_dp, 5.7_dp, -0.9_dp], [3, 3])
      call place(3, shifts, orb, aux, err)
      call check(error,.not. err%has_error(), "the bases were not built")
      if (allocated(error)) return

      call build_df_tensor(orb, aux, b, err)
      call check(error,.not. err%has_error(), "the dense tensor was not built")
      if (allocated(error)) return
      call build_df_jfit(orb, aux, jfit, err)
      call check(error,.not. err%has_error(), "the packed fit was not built")
      if (allocated(error)) return

      call random_density(orb%nao, d)
      call dense_j(b, d, j_dense)
      allocate (j_fit(orb%nao, orb%nao))
      call df_jfit_coulomb(jfit, d, j_fit)
      diff = relative_difference(j_fit, j_dense)
      call check(error, diff < 1.0e-12_dp, "J from the packed fit differs from B (B^T D) by "// &
                 real_to_text(diff))
      if (allocated(error)) return
      call check(error, maxval(abs(j_fit - transpose(j_fit))) == 0.0_dp, &
                 "J is not exactly symmetric")
      if (allocated(error)) return
      call check(error, jfit%bytes() < 8.0_dp*real(orb%nao, dp)**2*real(aux%nao, dp), &
                 "the packed fit is not smaller than the dense tensor")
      call orb%destroy()
      call aux%destroy()
   end subroutine test_cluster

   subroutine test_screened(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: orb, aux
      type(df_jfit_t) :: jfit, jfit_all
      type(error_t) :: err
      real(dp), allocatable :: b(:, :), d(:, :), j_dense(:, :), j_fit(:, :)
      real(dp), allocatable :: metric(:, :), half(:, :), three(:, :), c(:), flat(:)
      real(dp) :: shifts(3, 2)
      real(dp) :: diff
      integer :: n

      shifts = reshape([0.0_dp, 0.0_dp, 0.0_dp, 90.0_dp, 0.0_dp, 0.0_dp], [3, 2])
      call place(2, shifts, orb, aux, err)
      call check(error,.not. err%has_error(), "the bases were not built")
      if (allocated(error)) return
      n = orb%nao

      call build_df_tensor(orb, aux, b, err)
      call build_df_jfit(orb, aux, jfit, err)
      call build_df_jfit(orb, aux, jfit_all, err, screen=.false.)
      call check(error,.not. err%has_error(), "a fit was not built")
      if (allocated(error)) return

      call check(error, jfit%nrow < jfit_all%nrow, "the screen removed no function pairs")
      if (allocated(error)) return
      call check(error, jfit_all%nrow == n*(n + 1)/2, "the unscreened fit is not every pair")
      if (allocated(error)) return

      call random_density(n, d)
      allocate (j_fit(n, n))

      call dense_j(b, d, j_dense)
      call df_jfit_coulomb(jfit, d, j_fit)
      diff = relative_difference(j_fit, j_dense)
      call check(error, diff < 1.0e-12_dp, "screened: J differs from B (B^T D) by "// &
                 real_to_text(diff))
      if (allocated(error)) return

      ! No screen anywhere, and a factor of the metric that shares nothing with
      ! the Cholesky: J = T J^-1/2 J^-1/2 T^T D.
      call two_centre(aux, metric)
      call metric_inverse_sqrt(metric, half, err)
      call three_centre(orb, aux, three)
      flat = reshape(d, [n*n])
      c = matmul(half, matmul(half, matmul(flat, three)))
      flat = matmul(three, c)
      j_dense = reshape(flat, [n, n])
      call df_jfit_coulomb(jfit_all, d, j_fit)
      diff = relative_difference(j_fit, j_dense)
      call check(error, diff < 1.0e-10_dp, "unscreened: J differs from the unscreened reference by "// &
                 real_to_text(diff))
      call orb%destroy()
      call aux%destroy()
   end subroutine test_screened

   subroutine test_rks_energy(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: orb, aux
      type(xc_context_t) :: ctx_new, ctx_old
      type(rhf_result_t) :: scf_new, scf_old
      type(error_t) :: err
      real(dp) :: shifts(3, 2)
      real(dp) :: diff

      if (.not. xc_available()) then
         call skip_test(error, "this build has no libxc")
         return
      end if
      shifts = reshape([0.0_dp, 0.0_dp, 0.0_dp, 5.6_dp, 0.0_dp, 0.8_dp], [3, 2])
      call place(2, shifts, orb, aux, err)
      call check(error,.not. err%has_error(), "the bases were not built")
      if (allocated(error)) return

      call xc_context_create(orb, "pbe", ctx_new, err, level=1)
      call check(error,.not. err%has_error(), "the functional was not built")
      if (allocated(error)) return
      call run_czt_rhf(orb, 20, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf_new, err, aux=aux, xc=ctx_new)
      call check(error,.not. err%has_error(), "the packed-fit SCF failed")
      if (allocated(error)) return

      call set_df_jfit_enabled(.false.)
      call xc_context_create(orb, "pbe", ctx_old, err, level=1)
      call run_czt_rhf(orb, 20, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf_old, err, aux=aux, xc=ctx_old)
      call set_df_jfit_enabled(.true.)
      call check(error,.not. err%has_error(), "the dense-tensor SCF failed")
      if (allocated(error)) return

      call check(error, scf_new%converged .and. scf_old%converged, "an SCF did not converge: "// &
                 real_to_text(scf_new%energy)//" after "//real_to_text(real(scf_new%iterations, dp))// &
                 " iterations against "//real_to_text(scf_old%energy)//" after "// &
                 real_to_text(real(scf_old%iterations, dp)))
      if (allocated(error)) return
      diff = abs(scf_new%energy - scf_old%energy)
      call check(error, diff < 1.0e-9_dp, "RKS energies differ by "//real_to_text(diff)// &
                 " hartree between the packed and the dense fit")
      call orb%destroy()
      call aux%destroy()
   end subroutine test_rks_energy

   function real_to_text(value) result(text)
      !! A number in a failure message, so the message says how badly
      real(dp), intent(in) :: value
      character(len=:), allocatable :: text

      character(len=32) :: buffer

      write (buffer, "(es16.8)") value
      text = trim(adjustl(buffer))
   end function real_to_text

end module test_mqc_czt_jfit

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_jfit, only: collect_mqc_czt_jfit_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_jfit", collect_mqc_czt_jfit_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
