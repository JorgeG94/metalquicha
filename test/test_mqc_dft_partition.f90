!! Unit tests for the Becke partition of space into atomic cells
module test_mqc_dft_partition
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: iso_fortran_env, only: int64
   use pic_types, only: dp
   use mqc_error, only: error_t
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   use mqc_dft_partition, only: becke_partition_weights, becke_partition_derivatives, partition_scheme_name, &
                                PARTITION_BECKE, PARTITION_STRATMANN, &
                                ADJUST_NONE, ADJUST_BECKE, ADJUST_TREUTLER, &
                                becke_cutoff_derivative, stratmann_cutoff_derivative, &
                                becke_cutoff_second_derivative, &
                                stratmann_cutoff_second_derivative
   implicit none
   private

   public :: collect_mqc_dft_partition

   integer, parameter :: N_DIM = 3

contains

   subroutine collect_mqc_dft_partition(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("partition_weights_sum_to_one", test_unity), &
                  new_unittest("partition_of_one_atom_is_everything", test_single_atom), &
                  new_unittest("partition_splits_identical_atoms_evenly", test_symmetry), &
                  new_unittest("partition_favours_the_nearer_atom", test_nearer), &
                  new_unittest("partition_size_adjustment_favours_the_larger", test_adjust), &
                  new_unittest("partition_stratmann_saturates_far_away", test_stratmann_cutoff), &
                  new_unittest("partition_scheme_names", test_names), &
                  new_unittest("partition_rejects_bad_input", test_bad_input), &
                  new_unittest("cutoff_second_derivatives", test_cutoff_second), &
                  new_unittest("partition_becke_screening_bounded", test_screening_becke), &
                  new_unittest("partition_stratmann_screening_exact", test_screening_stratmann), &
                  new_unittest("partition_screening_small_systems", test_screening_small), &
                  new_unittest("partition_thread_count_independent", test_thread_count), &
                  new_unittest("partition_derivative_screening", test_derivative_screening) &
                  ]
   end subroutine collect_mqc_dft_partition

   subroutine water(coords, z)
      real(dp), allocatable, intent(out) :: coords(:, :)
      integer, allocatable, intent(out) :: z(:)
      allocate (coords(N_DIM, 3), z(3))
      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, -1.4308_dp, 1.1078_dp, &
                        0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3])
      z = [8, 1, 1]
   end subroutine water

   subroutine test_cutoff_second_helper(error, becke, x, tag)
      !! One cutoff, one point: the analytic second derivative against a central
      !! difference of the analytic first
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: becke
      real(dp), intent(in) :: x
      character(len=*), intent(in) :: tag

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp) :: fd, exact

      if (becke) then
         fd = (becke_cutoff_derivative(x + H) - becke_cutoff_derivative(x - H))/(2.0_dp*H)
         exact = becke_cutoff_second_derivative(x)
      else
         fd = (stratmann_cutoff_derivative(x + H) - stratmann_cutoff_derivative(x - H))/(2.0_dp*H)
         exact = stratmann_cutoff_second_derivative(x)
      end if
      call check(error, exact, fd, thr=1.0e-5_dp, more=tag)
   end subroutine test_cutoff_second_helper

   subroutine test_cutoff_second(error)
      !! Second derivatives of both cutoffs, against differences of the first
      !!
      !! These are the innermost scalars of a partition-weight Hessian: every
      !! term above them carries a chain factor built from these, so a sign
      !! error here is invisible until a whole Hessian is wrong by a few
      !! percent. Differencing the *analytic first derivative* rather than the
      !! value is deliberate -- it tests the new code against code that is
      !! already validated against finite differences itself, so a disagreement
      !! localises to the second derivative rather than to either.
      type(error_type), allocatable, intent(out) :: error

      call test_cutoff_second_helper(error, .true., 0.3_dp, "becke at nu = 0.3")
      if (allocated(error)) return
      call test_cutoff_second_helper(error, .true., -0.7_dp, "becke at nu = -0.7")
      if (allocated(error)) return
      call test_cutoff_second_helper(error, .true., 0.0_dp, "becke at nu = 0")
      if (allocated(error)) return
      call test_cutoff_second_helper(error, .false., 0.2_dp, "stratmann at mu = 0.2")
      if (allocated(error)) return
      call test_cutoff_second_helper(error, .false., -0.4_dp, "stratmann at mu = -0.4")
      if (allocated(error)) return
      ! Outside the support both derivatives are identically zero, and staying
      ! zero there is what keeps the far field from contributing noise.
      call check(error, stratmann_cutoff_second_derivative(1.0_dp) == 0.0_dp, &
                 "the Stratmann cutoff must be flat outside its support")
   end subroutine test_cutoff_second

   subroutine test_unity(error)
      !! The defining property: the cells partition space, so summing the
      !! weight of every atom at any point must give exactly one.
      !!
      !! Strong for the price. It fails if the pair loop is asymmetric, if a
      !! size-adjustment shift has the wrong sign, or if an atom is dropped
      !! from a product -- none of which a single-atom weight would reveal.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: coords(:, :), probes(:, :), w(:), unity(:)
      integer, allocatable :: z(:), owner(:)
      type(error_t) :: err
      integer :: schemes(2), adjusts(3)
      integer :: s, a, ia, i, j, k, n
      real(dp) :: worst

      call water(coords, z)
      schemes = [PARTITION_BECKE, PARTITION_STRATMANN]
      adjusts = [ADJUST_NONE, ADJUST_BECKE, ADJUST_TREUTLER]

      allocate (probes(N_DIM, 125))
      n = 0
      do i = 0, 4
         do j = 0, 4
            do k = 0, 4
               n = n + 1
               probes(:, n) = [-4.0_dp + 2.0_dp*i, -4.0_dp + 2.0_dp*j, -4.0_dp + 2.0_dp*k]
            end do
         end do
      end do
      allocate (w(n), owner(n), unity(n))

      worst = 0.0_dp
      do s = 1, 2
         do a = 1, 3
            unity = 0.0_dp
            do ia = 1, 3
               owner = ia
               call becke_partition_weights(probes, coords, z, owner, &
                                            schemes(s), adjusts(a), w, err)
               call check(error,.not. err%has_error(), "partition must build: "//err%get_full_trace())
               if (allocated(error)) return
               unity = unity + w
            end do
            worst = max(worst, maxval(abs(unity - 1.0_dp)))
         end do
      end do
      call check(error, worst < 1.0e-13_dp, "weights must sum to one everywhere")
   end subroutine test_unity

   subroutine test_single_atom(error)
      !! One atom owns all of space; the pair loop has nothing to do
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: coords(N_DIM, 1), probes(N_DIM, 3), w(3)
      integer :: z(1), owner(3)
      type(error_t) :: err

      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 1])
      z = [8]
      probes = reshape([0.1_dp, 0.0_dp, 0.0_dp, &
                        5.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, 0.0_dp, -20.0_dp], [N_DIM, 3])
      owner = 1
      call becke_partition_weights(probes, coords, z, owner, PARTITION_BECKE, &
                                   ADJUST_BECKE, w, err)
      call check(error,.not. err%has_error(), "one atom must be allowed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, all(abs(w - 1.0_dp) < 1.0e-15_dp), "one atom owns everything")
   end subroutine test_single_atom

   subroutine test_symmetry(error)
      !! Two identical atoms split their midpoint exactly in half
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: coords(N_DIM, 2), probes(N_DIM, 1), w(1)
      integer :: z(2), owner(1)
      type(error_t) :: err

      coords = reshape([-2.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 2])
      z = [8, 8]
      probes = reshape([0.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 1])

      owner = 1
      call becke_partition_weights(probes, coords, z, owner, PARTITION_BECKE, &
                                   ADJUST_BECKE, w, err)
      call check(error, abs(w(1) - 0.5_dp) < 1.0e-14_dp, &
                 "identical atoms must split the midpoint evenly")
   end subroutine test_symmetry

   subroutine test_nearer(error)
      !! A point beside one nucleus belongs almost entirely to it
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: coords(N_DIM, 2), probes(N_DIM, 1), w(1)
      integer :: z(2), owner(1)
      type(error_t) :: err

      coords = reshape([-2.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 2])
      z = [8, 8]
      probes = reshape([-1.95_dp, 0.0_dp, 0.0_dp], [N_DIM, 1])

      owner = 1
      call becke_partition_weights(probes, coords, z, owner, PARTITION_BECKE, &
                                   ADJUST_BECKE, w, err)
      call check(error, w(1) > 0.99_dp, "a point at a nucleus belongs to that atom")
   end subroutine test_nearer

   subroutine test_adjust(error)
      !! The size adjustment moves the boundary toward the smaller atom
      !!
      !! With no adjustment the midpoint of any pair splits evenly whatever the
      !! elements are. With it, oxygen takes more than half from hydrogen --
      !! which is the whole point of the correction, and its sign is the easiest
      !! thing to get backwards.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: coords(N_DIM, 2), probes(N_DIM, 1), w(1)
      integer :: z(2), owner(1)
      type(error_t) :: err
      real(dp) :: plain, adjusted

      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, 4.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 2])
      z = [8, 1]
      probes = reshape([2.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 1])
      owner = 1

      call becke_partition_weights(probes, coords, z, owner, PARTITION_BECKE, &
                                   ADJUST_NONE, w, err)
      plain = w(1)
      call becke_partition_weights(probes, coords, z, owner, PARTITION_BECKE, &
                                   ADJUST_BECKE, w, err)
      adjusted = w(1)

      call check(error, abs(plain - 0.5_dp) < 1.0e-14_dp, &
                 "without adjustment the midpoint splits evenly")
      if (allocated(error)) return
      call check(error, adjusted > plain, &
                 "oxygen must take more than half from hydrogen")
   end subroutine test_adjust

   subroutine test_stratmann_cutoff(error)
      !! Stratmann reaches exactly 0 and 1, which is what makes it screenable
      !!
      !! Becke's cutoff is smooth and never quite reaches either, so a point
      !! deep inside one cell has a tiny but non-zero weight for the other.
      !! Stratmann's is exactly 1, and that difference is the reason it exists.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: coords(N_DIM, 2), probes(N_DIM, 1), w(1)
      integer :: z(2), owner(1)
      type(error_t) :: err

      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, 30.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 2])
      z = [8, 8]
      probes = reshape([0.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 1])
      owner = 1

      call becke_partition_weights(probes, coords, z, owner, PARTITION_STRATMANN, &
                                   ADJUST_NONE, w, err)
      call check(error, w(1) == 1.0_dp, "Stratmann must saturate to exactly one")
   end subroutine test_stratmann_cutoff

   subroutine test_names(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, partition_scheme_name(PARTITION_BECKE), "Becke")
      if (allocated(error)) return
      call check(error, partition_scheme_name(PARTITION_STRATMANN), "Stratmann")
      if (allocated(error)) return
      call check(error, partition_scheme_name(-1), "unknown")
   end subroutine test_names

   subroutine test_bad_input(error)
      !! Mismatched arrays and out-of-range owners are refused
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: coords(N_DIM, 2), probes(N_DIM, 2), w(2)
      integer :: owner(2)
      type(error_t) :: err

      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 2])
      probes = coords
      owner = 1

      call becke_partition_weights(probes, coords, [8], owner, PARTITION_BECKE, &
                                   ADJUST_BECKE, w, err)
      call check(error, err%has_error(), "a short atomic_numbers must be refused")
      if (allocated(error)) return

      call err%clear()
      owner = [1, 99]
      call becke_partition_weights(probes, coords, [8, 8], owner, PARTITION_BECKE, &
                                   ADJUST_BECKE, w, err)
      call check(error, err%has_error(), "an owner past the atom list must be refused")
      if (allocated(error)) return

      call err%clear()
      owner = 1
      call becke_partition_weights(probes, coords, [8, 8], owner, -5, &
                                   ADJUST_BECKE, w, err)
      call check(error, err%has_error(), "an unknown scheme must be refused")
   end subroutine test_bad_input

   subroutine make_cluster(n_water, coords, z)
      !! `n_water` slightly distorted waters on a cubic lattice, 5.6 bohr apart
      integer, intent(in) :: n_water
      real(dp), allocatable, intent(out) :: coords(:, :)
      integer, allocatable, intent(out) :: z(:)

      real(dp), allocatable :: w0(:, :)
      integer, allocatable :: z0(:)
      integer :: iw, ia, side, ix, iy, iz, seed
      real(dp) :: offset(N_DIM)

      call water(w0, z0)
      allocate (coords(N_DIM, 3*n_water), z(3*n_water))
      side = 1
      do while (side**3 < n_water)
         side = side + 1
      end do

      seed = 12345
      iw = 0
      do iz = 0, side - 1
         do iy = 0, side - 1
            do ix = 0, side - 1
               iw = iw + 1
               if (iw > n_water) exit
               offset = 5.6_dp*real([ix, iy, iz], dp)
               do ia = 1, 3
                  coords(:, 3*(iw - 1) + ia) = w0(:, ia) + offset &
                                               + 0.3_dp*[next_random(seed), next_random(seed), next_random(seed)]
                  z(3*(iw - 1) + ia) = z0(ia)
               end do
            end do
         end do
      end do
   end subroutine make_cluster

   function next_random(seed) result(x)
      !! A reproducible number in (-0.5, 0.5), from a linear congruential generator
      integer, intent(inout) :: seed
      real(dp) :: x

      seed = int(mod(int(seed, kind=int64)*1103515245_int64 + 12345_int64, 2147483648_int64))
      x = real(seed, dp)/2147483648.0_dp - 0.5_dp
   end function next_random

   subroutine make_points(coords, per_atom, points, owner)
      !! `per_atom` points about each atom, from the nucleus out to 9 bohr
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: per_atom
      real(dp), allocatable, intent(out) :: points(:, :)
      integer, allocatable, intent(out) :: owner(:)

      integer :: ia, ip, n, seed
      real(dp) :: r
      real(dp) :: dir(N_DIM)

      n = size(coords, 2)
      allocate (points(N_DIM, n*per_atom), owner(n*per_atom))
      seed = 777
      do ia = 1, n
         do ip = 1, per_atom
            dir = [next_random(seed), next_random(seed), next_random(seed)]
            dir = dir/norm2(dir)
            r = 0.05_dp + 9.0_dp*(0.5_dp + next_random(seed))**2
            points(:, (ia - 1)*per_atom + ip) = coords(:, ia) + r*dir
            owner((ia - 1)*per_atom + ip) = ia
         end do
      end do
   end subroutine make_points

   subroutine compare_screening(error, n_water, scheme, adjust, thr, tag)
      !! Screened against unscreened weights on one cluster
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: n_water, scheme, adjust
      real(dp), intent(in) :: thr
      character(len=*), intent(in) :: tag

      real(dp), allocatable :: coords(:, :), points(:, :), fast(:), ref(:)
      integer, allocatable :: z(:), owner(:)
      type(error_t) :: err

      call make_cluster(n_water, coords, z)
      call make_points(coords, 12, points, owner)
      allocate (fast(size(owner)), ref(size(owner)))

      call becke_partition_weights(points, coords, z, owner, scheme, adjust, fast, err)
      call becke_partition_weights(points, coords, z, owner, scheme, adjust, ref, err, screening=.false.)
      call check(error,.not. err%has_error(), "partition reported an error: "//tag)
      if (allocated(error)) return

      ! Not a vacuous comparison: the weights have to span the range, both
      ! points that belong to their atom and points that do not.
      call check(error, maxval(ref) > 0.99_dp .and. minval(ref) < 0.5_dp, &
                 tag//": the probe points do not exercise the partition")
      if (allocated(error)) return

      call check(error, maxval(abs(fast - ref)), 0.0_dp, thr=thr, more=tag)
   end subroutine compare_screening

   subroutine test_screening_becke(error)
      !! The bounded Becke screening stays within 1e-14 of the full evaluation
      !!
      !! `n_atoms * SCREEN_TOL` is the stated bound, and it is 1e-14 up to a few
      !! thousand atoms; the clusters here are 48, 60 and 96 atoms, all above `MIN_BOUNDED_ATOMS`, under all
      !! three size adjustments.
      type(error_type), allocatable, intent(out) :: error

      call compare_screening(error, 16, PARTITION_BECKE, ADJUST_TREUTLER, 1.0e-14_dp, "becke 48 atoms")
      if (allocated(error)) return
      call compare_screening(error, 20, PARTITION_BECKE, ADJUST_TREUTLER, 1.0e-14_dp, "becke 60 atoms")
      if (allocated(error)) return
      call compare_screening(error, 20, PARTITION_BECKE, ADJUST_BECKE, 1.0e-14_dp, "becke 60 atoms, Bragg")
      if (allocated(error)) return
      call compare_screening(error, 20, PARTITION_BECKE, ADJUST_NONE, 1.0e-14_dp, "becke 60 atoms, no adjust")
      if (allocated(error)) return
      call compare_screening(error, 32, PARTITION_BECKE, ADJUST_TREUTLER, 1.0e-14_dp, "becke 96 atoms")
   end subroutine test_screening_becke

   subroutine test_screening_stratmann(error)
      !! Stratmann screening skips only exact zeros, so nothing moves at all
      type(error_type), allocatable, intent(out) :: error

      call compare_screening(error, 4, PARTITION_STRATMANN, ADJUST_TREUTLER, 0.0_dp, "stratmann 12 atoms")
      if (allocated(error)) return
      call compare_screening(error, 12, PARTITION_STRATMANN, ADJUST_TREUTLER, 0.0_dp, "stratmann 36 atoms")
      if (allocated(error)) return
      call compare_screening(error, 12, PARTITION_STRATMANN, ADJUST_NONE, 0.0_dp, "stratmann 36 atoms, no adjust")
      if (allocated(error)) return
      call compare_screening(error, 32, PARTITION_STRATMANN, ADJUST_BECKE, 0.0_dp, "stratmann 96 atoms")
   end subroutine test_screening_stratmann

   subroutine test_screening_small(error)
      !! Water and a water dimer, below the size where screening starts
      type(error_type), allocatable, intent(out) :: error

      call compare_screening(error, 1, PARTITION_BECKE, ADJUST_TREUTLER, 0.0_dp, "becke water")
      if (allocated(error)) return
      call compare_screening(error, 2, PARTITION_STRATMANN, ADJUST_TREUTLER, 0.0_dp, "stratmann dimer")
      if (allocated(error)) return
      call compare_screening(error, 2, PARTITION_BECKE, ADJUST_BECKE, 0.0_dp, "becke dimer")
   end subroutine test_screening_small

   subroutine test_thread_count(error)
      !! One thread and several give the same weights and derivatives, bit for bit
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: coords(:, :), points(:, :), w1(:), wn(:), d1(:, :, :), dn(:, :, :)
      integer, allocatable :: z(:), owner(:)
      type(error_t) :: err
      integer :: saved, scheme

      saved = 1
!$    saved = omp_get_max_threads()
      if (saved < 2) return

      call make_cluster(12, coords, z)
      call make_points(coords, 6, points, owner)
      allocate (w1(size(owner)), wn(size(owner)))
      allocate (d1(N_DIM, size(z), size(owner)), dn(N_DIM, size(z), size(owner)))

      do scheme = PARTITION_BECKE, PARTITION_STRATMANN
!$       call omp_set_num_threads(1)
         call becke_partition_weights(points, coords, z, owner, scheme, ADJUST_TREUTLER, w1, err)
         call becke_partition_derivatives(points, coords, z, owner, scheme, ADJUST_TREUTLER, d1, err)
!$       call omp_set_num_threads(saved)
         call becke_partition_weights(points, coords, z, owner, scheme, ADJUST_TREUTLER, wn, err)
         call becke_partition_derivatives(points, coords, z, owner, scheme, ADJUST_TREUTLER, dn, err)
         call check(error,.not. err%has_error(), "partition reported an error")
         if (allocated(error)) return
         call check(error, maxval(abs(w1 - wn)), 0.0_dp, thr=1.0e-15_dp, more="threaded weights")
         if (allocated(error)) return
         call check(error, maxval(abs(d1 - dn)), 0.0_dp, thr=1.0e-15_dp, more="threaded derivatives")
         if (allocated(error)) return
      end do
   end subroutine test_thread_count

   subroutine test_derivative_screening(error)
      !! Stratmann derivative screening is exact, and both schemes match differences
      !!
      !! The finite differences are of the unscreened weights, so this checks the
      !! restructured derivative loop against the definition of the quantity
      !! rather than against itself.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-5_dp
      real(dp), allocatable :: coords(:, :), points(:, :), shifted(:, :), fast(:, :, :), ref(:, :, :)
      real(dp), allocatable :: wp(:), wm(:)
      integer, allocatable :: z(:), owner(:)
      type(error_t) :: err
      integer :: scheme, ia, ic, ip
      real(dp) :: worst, fd

      call make_cluster(4, coords, z)
      call make_points(coords, 3, points, owner)
      allocate (fast(N_DIM, size(z), size(owner)), ref(N_DIM, size(z), size(owner)))
      allocate (wp(size(owner)), wm(size(owner)), shifted(N_DIM, size(z)))

      do scheme = PARTITION_BECKE, PARTITION_STRATMANN
         call becke_partition_derivatives(points, coords, z, owner, scheme, ADJUST_TREUTLER, fast, err)
         call becke_partition_derivatives(points, coords, z, owner, scheme, ADJUST_TREUTLER, ref, err, &
                                          screening=.false.)
         call check(error,.not. err%has_error(), "partition reported an error")
         if (allocated(error)) return
         call check(error, maxval(abs(fast - ref)), 0.0_dp, thr=0.0_dp, more="screened derivatives")
         if (allocated(error)) return

         worst = 0.0_dp
         do ia = 1, size(z)
            do ic = 1, N_DIM
               shifted = coords
               shifted(ic, ia) = coords(ic, ia) + H
               call becke_partition_weights(points, shifted, z, owner, scheme, ADJUST_TREUTLER, wp, err, &
                                            screening=.false.)
               shifted(ic, ia) = coords(ic, ia) - H
               call becke_partition_weights(points, shifted, z, owner, scheme, ADJUST_TREUTLER, wm, err, &
                                            screening=.false.)
               do ip = 1, size(owner)
                  fd = (wp(ip) - wm(ip))/(2.0_dp*H)
                  worst = max(worst, abs(fd - ref(ic, ia, ip)))
               end do
            end do
         end do
         call check(error, worst, 0.0_dp, thr=1.0e-7_dp, more="derivative against differences")
         if (allocated(error)) return
      end do
   end subroutine test_derivative_screening

end module test_mqc_dft_partition

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_dft_partition, only: collect_mqc_dft_partition
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_dft_partition", collect_mqc_dft_partition)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
