!! Unit tests for the assembled molecular DFT grid
module test_mqc_dft_grid
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_dft_grid, only: dft_grid_t, build_dft_grid, grid_level_radial, &
                           grid_level_angular, MIN_GRID_LEVEL, MAX_GRID_LEVEL, &
                           DEFAULT_GRID_LEVEL
   use mqc_dft_prune, only: PRUNE_NONE
   implicit none
   private

   public :: collect_mqc_dft_grid

   integer, parameter :: N_DIM = 3
   real(dp), parameter :: PI = 3.14159265358979323846_dp

contains

   subroutine collect_mqc_dft_grid(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("grid_level_tables_grow_with_level", test_level_tables), &
                  new_unittest("grid_level_is_clamped_not_rejected", test_clamp), &
                  new_unittest("grid_point_count_matches_the_tables", test_counts), &
                  new_unittest("grid_integrates_a_gaussian_sum", test_gaussians), &
                  new_unittest("grid_integrates_cusped_functions", test_slaters), &
                  new_unittest("grid_improves_with_level", test_convergence), &
                  new_unittest("grid_overrides_apply_to_every_atom", test_override), &
                  new_unittest("grid_rejects_a_half_given_override", test_half_override), &
                  new_unittest("grid_rejects_a_mismatched_system", test_mismatch) &
                  ]
   end subroutine collect_mqc_dft_grid

   subroutine water(coords, z)
      real(dp), allocatable, intent(out) :: coords(:, :)
      integer, allocatable, intent(out) :: z(:)
      allocate (coords(N_DIM, 3), z(3))
      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, -1.4308_dp, 1.1078_dp, &
                        0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3])
      z = [8, 1, 1]
   end subroutine water

   subroutine test_level_tables(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: level

      do level = MIN_GRID_LEVEL + 1, MAX_GRID_LEVEL
         call check(error, grid_level_radial(8, level) >= grid_level_radial(8, level - 1), &
                    "radial count must not shrink with level")
         if (allocated(error)) return
         call check(error, grid_level_angular(8, level) >= grid_level_angular(8, level - 1), &
                    "angular count must not shrink with level")
         if (allocated(error)) return
      end do

      ! Heavier elements get at least as much as lighter ones.
      call check(error, grid_level_radial(79, DEFAULT_GRID_LEVEL) >= &
                 grid_level_radial(1, DEFAULT_GRID_LEVEL), &
                 "a heavy element needs at least as many radial shells")
   end subroutine test_level_tables

   subroutine test_clamp(error)
      !! Out-of-range levels saturate rather than erroring
      type(error_type), allocatable, intent(out) :: error

      call check(error, grid_level_radial(8, -5), grid_level_radial(8, MIN_GRID_LEVEL))
      if (allocated(error)) return
      call check(error, grid_level_radial(8, 999), grid_level_radial(8, MAX_GRID_LEVEL))
   end subroutine test_clamp

   subroutine test_counts(error)
      !! Unpruned, the count is exactly sum over atoms of n_radial * n_angular
      !!
      !! Says PRUNE_NONE explicitly since pruning became the default: the
      !! identity being checked is the unpruned one, and inheriting whatever
      !! the default happens to be would make this test silently change meaning
      !! the next time that default moves.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: coords(:, :)
      integer, allocatable :: z(:)
      type(dft_grid_t) :: grid
      type(error_t) :: err
      integer :: expected, ia

      call water(coords, z)
      call build_dft_grid(coords, z, grid, err, level=1, prune=PRUNE_NONE)
      call check(error,.not. err%has_error(), "grid must build")
      if (allocated(error)) return

      expected = 0
      do ia = 1, size(z)
         expected = expected + grid_level_radial(z(ia), 1)*grid_level_angular(z(ia), 1)
      end do
      call check(error, grid%n_points, expected)
      if (allocated(error)) return
      call check(error, size(grid%weights), expected)
      if (allocated(error)) return
      call check(error, all(grid%atom >= 1 .and. grid%atom <= 3), &
                 "every point must belong to a real atom")
      call grid%destroy()
   end subroutine test_counts

   subroutine test_gaussians(error)
      !! A sum of normalised Gaussians integrates to the number of atoms
      !!
      !! Smooth, so any error is in the assembly rather than the mesh: a
      !! missing r^2, a 4*pi applied twice or not at all, the partition
      !! multiplied into the wrong points. Each of those returns a smooth wrong
      !! number instead of failing.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: coords(:, :)
      integer, allocatable :: z(:)
      type(dft_grid_t) :: grid
      type(error_t) :: err
      real(dp) :: total

      call water(coords, z)
      call build_dft_grid(coords, z, grid, err, level=3)
      total = integrate_gaussians(grid, coords, 1.3_dp)
      ! Relative, and calibrated: level 3 is the production default rather than
      ! a converged grid, and it reaches 4.3e-9 here. An absolute 1e-8 bound
      ! looks equivalent and is not -- three atoms make it 1.3e-8.
      call check(error, abs(total - 3.0_dp)/3.0_dp < 1.0e-8_dp, &
                 "a sum of Gaussians must integrate to the atom count")
      call grid%destroy()
   end subroutine test_gaussians

   subroutine test_slaters(error)
      !! Cusped at every nucleus, which is what the radial mesh is built for
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: coords(:, :)
      integer, allocatable :: z(:)
      type(dft_grid_t) :: grid
      type(error_t) :: err
      real(dp) :: total

      call water(coords, z)
      call build_dft_grid(coords, z, grid, err, level=3)
      total = integrate_slaters(grid, coords, 1.7_dp)
      call check(error, abs(total - 3.0_dp)/3.0_dp < 1.0e-8_dp, &
                 "a sum of Slaters must integrate to the atom count")
      call grid%destroy()
   end subroutine test_slaters

   subroutine test_convergence(error)
      !! Refining must improve it
      !!
      !! A grid that is merely close at one level could be wrong in a way that
      !! cancels. Only unpruned: the NWChem scheme fixes its inner two zones at
      !! 50 and 86 points whatever the level, so a pruned grid saturates rather
      !! than converging, and asserting otherwise would assert something the
      !! scheme does not promise.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: coords(:, :)
      integer, allocatable :: z(:)
      type(dft_grid_t) :: grid
      type(error_t) :: err
      real(dp) :: prev, now
      integer :: level

      call water(coords, z)
      prev = huge(1.0_dp)
      do level = 1, 4
         call build_dft_grid(coords, z, grid, err, level=level)
         now = abs(integrate_gaussians(grid, coords, 1.3_dp) - 3.0_dp)
         call check(error, now < prev, "refining the level must reduce the error")
         call grid%destroy()
         if (allocated(error)) return
         prev = now
      end do
   end subroutine test_convergence

   subroutine test_override(error)
      !! An explicit size applies to every atom, ignoring the level tables
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: coords(:, :)
      integer, allocatable :: z(:)
      type(dft_grid_t) :: grid
      type(error_t) :: err

      call water(coords, z)
      call build_dft_grid(coords, z, grid, err, n_radial=20, n_angular=50, &
                          prune=PRUNE_NONE)
      call check(error,.not. err%has_error(), "an override must be accepted")
      if (allocated(error)) return
      call check(error, grid%n_points, 3*20*50)
      call grid%destroy()
   end subroutine test_override

   subroutine test_half_override(error)
      !! One of the two is refused rather than half-applied
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: coords(:, :)
      integer, allocatable :: z(:)
      type(dft_grid_t) :: grid
      type(error_t) :: err

      call water(coords, z)
      call build_dft_grid(coords, z, grid, err, n_radial=20)
      call check(error, err%has_error(), "n_radial without n_angular must be refused")
   end subroutine test_half_override

   subroutine test_mismatch(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: coords(N_DIM, 2)
      type(dft_grid_t) :: grid
      type(error_t) :: err

      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 2])
      call build_dft_grid(coords, [8], grid, err)
      call check(error, err%has_error(), "atomic_numbers must match the coordinates")
   end subroutine test_mismatch

   function integrate_gaussians(grid, coords, alpha) result(total)
      type(dft_grid_t), intent(in) :: grid
      real(dp), intent(in) :: coords(:, :)
      real(dp), intent(in) :: alpha
      real(dp) :: total, f, d2, norm
      integer :: k, ia

      norm = (alpha/PI)**1.5_dp
      total = 0.0_dp
      do k = 1, grid%n_points
         f = 0.0_dp
         do ia = 1, size(coords, 2)
            d2 = sum((grid%coords(:, k) - coords(:, ia))**2)
            f = f + norm*exp(-alpha*d2)
         end do
         total = total + grid%weights(k)*f
      end do
   end function integrate_gaussians

   function integrate_slaters(grid, coords, zeta) result(total)
      type(dft_grid_t), intent(in) :: grid
      real(dp), intent(in) :: coords(:, :)
      real(dp), intent(in) :: zeta
      real(dp) :: total, f, d, norm
      integer :: k, ia

      norm = zeta**3/PI
      total = 0.0_dp
      do k = 1, grid%n_points
         f = 0.0_dp
         do ia = 1, size(coords, 2)
            d = norm2(grid%coords(:, k) - coords(:, ia))
            f = f + norm*exp(-2.0_dp*zeta*d)
         end do
         total = total + grid%weights(k)*f
      end do
   end function integrate_slaters

end module test_mqc_dft_grid

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_dft_grid, only: collect_mqc_dft_grid
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_dft_grid", collect_mqc_dft_grid)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
