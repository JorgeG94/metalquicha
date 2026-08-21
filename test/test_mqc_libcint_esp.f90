!! The potential integral and its gradient, checked against each other
module test_mqc_libcint_esp
   !! `drinv_matrices` is the gradient with respect to the centre of what
   !! `esp_matrices` computes at that centre, so the two are checkable against each
   !! other with nothing external -- and they need to be, because libcint's
   !! `nabla-rinv` could be the gradient with respect to `r` or with respect to `C`
   !! and those differ by exactly a minus.
   !!
   !! **This test is what fixes the convention for every caller.** Charge transfer's
   !! dipole rank is `-sum_C d_C . grad_C <mu|1/|r-C||nu>`; with the sign the other way
   !! it is a plausible-looking energy of the wrong sign, and the only reference that
   !! would catch it is a GAMESS number one rung further up a ladder. Pinning it here
   !! instead means the rank above can be written without a sign to guess at.
   !!
   !! The reference is a central difference of `esp_matrices`, which is itself checked
   !! against PySCF elementwise. A central difference is second-order accurate, so at
   !! `h = 1e-4` Bohr the comparison is good to about `1e-8` -- far inside the gap
   !! between a right answer and a sign error, which is the thing being tested.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_esp, only: esp_matrices, drinv_matrices, ddrinv_matrices
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_libcint_esp_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_libcint_esp_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("esp_drinv_is_the_centre_gradient", test_drinv), &
                  new_unittest("esp_ddrinv_is_the_second_derivative", test_ddrinv) &
                  ]
   end subroutine collect_mqc_libcint_esp_tests

   subroutine test_drinv(error)
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: grad(:, :, :, :)
      real(dp), allocatable :: plus(:, :, :), minus(:, :, :)
      real(dp) :: c(3, 3), probe(3, 2), shifted(3, 2)
      real(dp), parameter :: H = 1.0e-4_dp
      real(dp) :: numerical, worst, scale
      integer :: z(3), x, p, i, j
      character(len=2) :: symbols(3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_libcint_molecule(z, symbols, c, "6-31g*", mol, err)
      call check(error,.not. err%has_error(), "building the molecule failed: "//err%get_full_trace())
      if (allocated(error)) return

      ! Two centres off the nuclei: on an atom the integral is still finite but its
      ! gradient is where the integrand is worst behaved, and a point that happens to
      ! sit on a symmetry element can hide a wrong component in a zero.
      probe = reshape([0.7_dp, 0.4_dp, 1.3_dp, &
                       -1.1_dp, 0.9_dp, -0.6_dp], [3, 2])

      call drinv_matrices(mol, probe, grad, err)
      call check(error,.not. err%has_error(), "drinv failed: "//err%get_full_trace())
      if (allocated(error)) return

      worst = 0.0_dp
      scale = 0.0_dp
      do p = 1, 2
         do x = 1, 3
            shifted = probe
            shifted(x, p) = probe(x, p) + H
            call esp_matrices(mol, shifted, plus, err)
            shifted(x, p) = probe(x, p) - H
            call esp_matrices(mol, shifted, minus, err)
            call check(error,.not. err%has_error(), "esp failed: "//err%get_full_trace())
            if (allocated(error)) return
            do j = 1, mol%nao
               do i = 1, mol%nao
                  numerical = (plus(i, j, p) - minus(i, j, p))/(2.0_dp*H)
                  worst = max(worst, abs(grad(i, j, x, p) - numerical))
                  scale = max(scale, abs(numerical))
               end do
            end do
            deallocate (plus, minus)
         end do
      end do

      ! The gradient is order 1 here, so this is a relative agreement of about 1e-8 --
      ! central-difference truncation, not integral error. A sign error would show up
      ! as `worst` of order `2*scale` instead.
      call check(error, worst < 1.0e-7_dp, &
                 "drinv is not the centre-gradient of the potential integral")
      if (allocated(error)) return
      ! Guards the guard: if both were zero the comparison above would pass for the
      ! wrong reason, and a centre far from the molecule would make that nearly true.
      call check(error, scale > 0.1_dp, "the probe gradients are too small to test")

      call mol%destroy()
      deallocate (grad)
   end subroutine test_drinv

   subroutine test_ddrinv(error)
      !! The second centre-derivative, against a difference of the first
      !!
      !! Differencing `drinv_matrices` rather than `esp_matrices` twice keeps this a
      !! first-order difference of an already-pinned quantity: two nested central
      !! differences of the potential would carry `h^-2` amplification and force a
      !! looser bound, which is exactly the regime where a wrong assembly could hide.
      !!
      !! What this catches is the assembly, since `ddrinv_matrices` is four terms across
      !! two libcint integrals and their transposes. It is deliberately checked on both
      !! diagonal and off-diagonal `(a,b)`: a mistake in the transposed placement leaves
      !! the diagonal components right and the off-diagonal ones wrong.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: second(:, :, :, :, :)
      real(dp), allocatable :: plus(:, :, :, :), minus(:, :, :, :)
      real(dp) :: c(3, 3), probe(3, 1), shifted(3, 1)
      real(dp), parameter :: H = 1.0e-5_dp
      real(dp) :: numerical, worst, scale, asym
      integer :: z(3), a, b, i, j
      character(len=2) :: symbols(3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_libcint_molecule(z, symbols, c, "6-31g*", mol, err)
      call check(error,.not. err%has_error(), "building the molecule failed: "//err%get_full_trace())
      if (allocated(error)) return

      ! **On a nucleus, deliberately.** Every centre an EFP potential supplies is one of
      ! its own atoms, so the coincident case is the only one production ever uses --
      ! and it is the one where the integrand is worst behaved and where libcint's own
      ! screening is most likely to take a different branch. A probe floating in empty
      ! space tests a regime no caller is in.
      probe = reshape(c(:, 1), [3, 1])
      call ddrinv_matrices(mol, probe, second, err)
      call check(error,.not. err%has_error(), "ddrinv failed: "//err%get_full_trace())
      if (allocated(error)) return

      worst = 0.0_dp
      scale = 0.0_dp
      asym = 0.0_dp
      do b = 1, 3
         shifted(:, 1) = probe(:, 1)
         shifted(b, 1) = probe(b, 1) + H
         call drinv_matrices(mol, shifted, plus, err)
         shifted(b, 1) = probe(b, 1) - H
         call drinv_matrices(mol, shifted, minus, err)
         call check(error,.not. err%has_error(), "drinv failed: "//err%get_full_trace())
         if (allocated(error)) return
         do a = 1, 3
            do j = 1, mol%nao
               do i = 1, mol%nao
                  numerical = (plus(i, j, a, 1) - minus(i, j, a, 1))/(2.0_dp*H)
                  worst = max(worst, abs(second(i, j, a, b, 1) - numerical))
                  scale = max(scale, abs(numerical))
                  asym = max(asym, abs(second(i, j, a, b, 1) - second(i, j, b, a, 1)))
               end do
            end do
         end do
         deallocate (plus, minus)
      end do

      ! Relative, and loose, on purpose. On a nucleus the second derivative is ~6e+02
      ! and the third -- which sets the truncation error -- is larger still, so a
      ! central difference converges to this at `h^2` from a long way off: the measured
      ! gap is 3.3e+01 at h=1e-2, 3.9e-01 at 1e-3, 3.9e-03 at 1e-4 and 3.9e-05 at 1e-5,
      ! a clean second-order sequence. Tightening this bound would only be measuring the
      ! difference formula. What it still separates by orders of magnitude is a wrong
      ! assembly, which does not converge at all.
      call check(error, worst < 1.0e-4_dp*scale, &
                 "ddrinv is not the second centre-derivative of the potential integral")
      if (allocated(error)) return
      call check(error, scale > 0.1_dp, "the probe second derivatives are too small")
      if (allocated(error)) return
      ! Partial derivatives commute, so this must hold exactly rather than numerically.
      ! It is the one property that would survive a consistent-but-wrong assembly, so it
      ! is checked separately from the comparison above rather than instead of it.
      call check(error, asym < 1.0e-12_dp, "ddrinv is not symmetric in its two indices")

      call mol%destroy()
      deallocate (second)
   end subroutine test_ddrinv

end module test_mqc_libcint_esp

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_libcint_esp, only: collect_mqc_libcint_esp_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_libcint_esp", collect_mqc_libcint_esp_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
