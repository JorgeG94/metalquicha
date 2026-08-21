!! Unit tests for the Treutler-Ahlrichs radial quadrature
module test_mqc_dft_radial
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_dft_radial, only: treutler_xi, bragg_radius, treutler_ahlrichs_radial, &
                             radial_volume_weights
   implicit none
   private

   public :: collect_mqc_dft_radial

   real(dp), parameter :: PI = 3.14159265358979323846_dp

contains

   subroutine collect_mqc_dft_radial(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("radial_xi_matches_the_published_table", test_xi), &
                  new_unittest("radial_bragg_radii_cover_the_table", test_bragg), &
                  new_unittest("radial_nodes_ascend_and_are_positive", test_nodes), &
                  new_unittest("radial_weights_are_positive", test_weights_positive), &
                  new_unittest("radial_reaches_further_with_more_points", test_reach), &
                  new_unittest("radial_integrates_a_gaussian_exactly", test_gaussian), &
                  new_unittest("radial_tail_converges_as_the_mesh_grows", test_tail), &
                  new_unittest("radial_volume_weight_is_four_pi_r_squared_dr", test_volume), &
                  new_unittest("radial_refuses_an_empty_mesh", test_empty) &
                  ]
   end subroutine collect_mqc_dft_radial

   subroutine test_xi(error)
      !! H to Kr are the values from the original paper; heavier from Psi4
      type(error_type), allocatable, intent(out) :: error

      call check(error, treutler_xi(1), 0.8_dp, thr=1.0e-12_dp)
      if (allocated(error)) return
      call check(error, treutler_xi(36), 0.9_dp, thr=1.0e-12_dp)
      if (allocated(error)) return
      call check(error, treutler_xi(79), 1.5_dp, thr=1.0e-12_dp)
      if (allocated(error)) return
      ! Index 0 is a ghost atom, and beyond the table the fallback is 1.0.
      call check(error, treutler_xi(0), 1.0_dp, thr=1.0e-12_dp)
      if (allocated(error)) return
      call check(error, treutler_xi(500), 1.0_dp, thr=1.0e-12_dp)
   end subroutine test_xi

   subroutine test_bragg(error)
      !! The table runs to Z=130, so nothing real falls back
      !!
      !! This is the property the branch exists for: the previous helper
      !! returned 1.0 beyond krypton, which silently gave every heavy element
      !! an unoptimised mesh.
      type(error_type), allocatable, intent(out) :: error
      integer :: z

      do z = 1, 118
         call check(error, bragg_radius(z) > 0.0_dp, "every element needs a radius")
         if (allocated(error)) return
      end do

      ! Past krypton must not be the 1.0 fallback any more.
      call check(error, abs(bragg_radius(79) - 1.0_dp) > 1.0e-6_dp, &
                 "gold must have a real radius, not the fallback")
      if (allocated(error)) return
      call check(error, abs(bragg_radius(92) - 1.0_dp) > 1.0e-6_dp, &
                 "uranium must have a real radius, not the fallback")
   end subroutine test_bragg

   subroutine test_nodes(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:), dr(:)
      type(error_t) :: err
      integer :: i, n

      n = 40
      call treutler_ahlrichs_radial(n, 8, r, dr, err)
      call check(error,.not. err%has_error(), "mesh must build: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, size(r), n)
      if (allocated(error)) return

      call check(error, r(1) > 0.0_dp, "the innermost node must be positive")
      if (allocated(error)) return
      do i = 2, n
         call check(error, r(i) > r(i - 1), "nodes must ascend")
         if (allocated(error)) return
      end do
   end subroutine test_nodes

   subroutine test_weights_positive(error)
      !! A negative radial weight would make a density integrate negative
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:), dr(:)
      type(error_t) :: err
      integer :: z

      do z = 1, 92, 13
         call treutler_ahlrichs_radial(50, z, r, dr, err)
         call check(error, all(dr > 0.0_dp), "mapping weights must be positive")
         deallocate (r, dr)
         if (allocated(error)) return
      end do
   end subroutine test_weights_positive

   subroutine test_reach(error)
      !! More points means further out, not just denser
      !!
      !! The outer reach is what limits a diffuse density, so it is worth
      !! asserting that it grows rather than assuming it.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r1(:), d1(:), r2(:), d2(:)
      type(error_t) :: err

      call treutler_ahlrichs_radial(50, 1, r1, d1, err)
      call treutler_ahlrichs_radial(200, 1, r2, d2, err)
      call check(error, maxval(r2) > maxval(r1), "a finer mesh must reach further")
   end subroutine test_reach

   subroutine test_gaussian(error)
      !! exp(-r^2) is contained inside the mesh, so this probes the Jacobian
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:), dr(:), w(:)
      type(error_t) :: err
      real(dp) :: total

      call treutler_ahlrichs_radial(200, 1, r, dr, err)
      allocate (w(size(r)))
      call radial_volume_weights(r, dr, w)
      total = sum(w*exp(-r*r))
      call check(error, abs(total - PI**1.5_dp) < 1.0e-12_dp, &
                 "a contained Gaussian must integrate exactly")
   end subroutine test_gaussian

   subroutine test_tail(error)
      !! exp(-r) is not contained, so its error is the truncated tail
      !!
      !! At 200 points the mesh reaches about 17 Bohr, where exp(-r) is still
      !! 4e-8. That is a property of the mapping rather than a defect -- PySCF's
      !! mesh gives the same -- so the assertion is that refining fixes it.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:), dr(:), w(:)
      type(error_t) :: err
      real(dp) :: exact, coarse, fine

      exact = 8.0_dp*PI

      call treutler_ahlrichs_radial(100, 1, r, dr, err)
      allocate (w(size(r)))
      call radial_volume_weights(r, dr, w)
      coarse = abs(sum(w*exp(-r)) - exact)/exact
      deallocate (r, dr, w)

      call treutler_ahlrichs_radial(800, 1, r, dr, err)
      allocate (w(size(r)))
      call radial_volume_weights(r, dr, w)
      fine = abs(sum(w*exp(-r)) - exact)/exact

      call check(error, fine < coarse, "refining must reduce the tail error")
      if (allocated(error)) return
      call check(error, fine < 1.0e-8_dp, "800 points must get exp(-r) to 1e-8")
   end subroutine test_tail

   subroutine test_volume(error)
      !! radial_volume_weights applies 4*pi*r^2 and nothing else
      !!
      !! Worth pinning: the molecular grid supplies the 4*pi itself because the
      !! Lebedev weights sum to one, and applying it in both places is a silent
      !! factor of 4*pi on every integral.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:), dr(:), w(:)
      type(error_t) :: err
      integer :: k

      call treutler_ahlrichs_radial(20, 6, r, dr, err)
      allocate (w(size(r)))
      call radial_volume_weights(r, dr, w)
      do k = 1, size(r)
         call check(error, w(k), 4.0_dp*PI*r(k)*r(k)*dr(k), thr=1.0e-14_dp)
         if (allocated(error)) return
      end do
   end subroutine test_volume

   subroutine test_empty(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:), dr(:)
      type(error_t) :: err

      call treutler_ahlrichs_radial(0, 1, r, dr, err)
      call check(error, err%has_error(), "a mesh of no points must be refused")
   end subroutine test_empty

end module test_mqc_dft_radial

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_dft_radial, only: collect_mqc_dft_radial
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_dft_radial", collect_mqc_dft_radial)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
