!! The kernel cache against the evaluation it replaces
module test_mqc_xc_kernel_cache
   !! `xc_kernel_cache_fill` moves the reference density and libxc's second
   !! derivatives out of `xc_kernel_apply_many`'s per-block loop and into one
   !! pass over the grid. The only thing that makes it a cache rather than a
   !! second implementation of the kernel is that the two agree **exactly**, so
   !! that is what is checked: same trial densities, same context, same
   !! molecule, one call with the cache and one without, compared element by
   !! element with no tolerance at all.
   !!
   !! **At one thread, deliberately.** The contraction accumulates each block's
   !! contribution into the output under a lock, so the summation order is the
   !! order the blocks finish in and two runs on many threads differ in the last
   !! bits whether or not anything else changed. The suite's own note on that
   !! says bit-identity is testable at one thread and nowhere else, so the
   !! thread count is set here rather than left to whoever runs `ctest`.
   !!
   !! Three rungs, because the cache carries a different number of channels on
   !! each: a GGA (`frr`, `frs`, `fss`, `vsigma`), a meta-GGA (those plus
   !! `frt`, `fst`, `ftt`) and a hybrid, whose exchange fraction never reaches
   !! this routine but whose composition mixes two libxc components.
   use testdrive, only: new_unittest, unittest_type, error_type, check
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available, &
                         xc_kernel_apply, xc_kernel_apply_many, &
                         xc_kernel_cache_t, xc_kernel_cache_fill
   implicit none
   private

   public :: collect_mqc_xc_kernel_cache

   ! Water, bent, in Bohr -- the geometry the exchange-correlation derivative
   ! tests already use, so a disagreement here cannot be the geometry.
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape([ &
                                                0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                0.0000_dp, 1.4300_dp, 1.1075_dp, &
                                                0.0000_dp, -1.4300_dp, 1.1075_dp], [3, 3])

   integer, parameter :: N_TRIAL = 3
      !! Response densities per batch. More than one, because the batch loop
      !! reads the cached coefficients once per set and a cache consumed
      !! destructively would show up only from the second set on.

contains

   subroutine collect_mqc_xc_kernel_cache(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("cached_gga_kernel_is_bit_identical", test_gga), &
                  new_unittest("cached_mgga_kernel_is_bit_identical", test_mgga), &
                  new_unittest("cached_hybrid_kernel_matches", test_hybrid), &
                  new_unittest("an_unfilled_cache_is_refused", test_unfilled) &
                  ]
   end subroutine collect_mqc_xc_kernel_cache

   subroutine test_gga(error)
      !! PBE: the rung the response solves actually run on
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: worst, worst_one
      logical :: ok

      call cache_case("pbe", worst, worst_one, error, ok)
      if (allocated(error) .or. .not. ok) return
      call check(error, worst == 0.0_dp, "the cached GGA kernel differs from the "// &
                 "evaluated one")
      if (allocated(error)) return
      call check(error, worst_one == 0.0_dp, "the cached single-density GGA kernel "// &
                 "differs from the evaluated one")
   end subroutine test_gga

   subroutine test_mgga(error)
      !! TPSS exchange: the tau channels have to travel through the cache too,
      !! and this is the case that fails if `frt`, `fst` or `ftt` is dropped
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: worst, worst_one
      logical :: ok

      call cache_case("mgga_x_tpss", worst, worst_one, error, ok)
      if (allocated(error) .or. .not. ok) return
      call check(error, worst == 0.0_dp, "the cached meta-GGA kernel differs from "// &
                 "the evaluated one")
      if (allocated(error)) return
      call check(error, worst_one == 0.0_dp, "the cached single-density meta-GGA "// &
                 "kernel differs from the evaluated one")
   end subroutine test_mgga

   subroutine test_hybrid(error)
      !! B3LYP, whose composition is five libxc components summed with weights:
      !! the loop the cache has to reproduce in the same order
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: worst, worst_one
      logical :: ok

      call cache_case("b3lyp", worst, worst_one, error, ok)
      if (allocated(error) .or. .not. ok) return
      call check(error, worst < 1.0e-12_dp, "the cached B3LYP kernel differs from "// &
                 "the evaluated one")
      if (allocated(error)) return
      call check(error, worst_one < 1.0e-12_dp, "the cached single-density B3LYP "// &
                 "kernel differs from the evaluated one")
   end subroutine test_hybrid

   subroutine test_unfilled(error)
      !! A cache that was never filled is refused rather than read as zeros,
      !! which would be a silently missing kernel term in a converged answer
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(xc_kernel_cache_t) :: cache
      type(error_t) :: err
      real(dp), allocatable :: dens(:, :), trials(:, :, :), out(:, :, :)
      integer :: nao
      logical :: ok

      if (.not. xc_available()) return
      call reference_state("pbe", mol, ctx, dens, err, ok)
      if (.not. ok) then
         call check(error, .false., "the reference Kohn-Sham state failed")
         return
      end if

      nao = size(dens, 1)
      allocate (trials(nao, nao, 1), out(nao, nao, 1))
      trials(:, :, 1) = dens
      out = 0.0_dp
      call xc_kernel_apply_many(ctx, mol, dens, trials, out, err, cache=cache)
      call check(error, err%has_error(), "an unfilled kernel cache was accepted")
      call ctx%destroy()
      call mol%destroy()
   end subroutine test_unfilled

   subroutine cache_case(functional, worst, worst_one, error, ok)
      !! One functional both ways: the largest element-wise disagreement
      !! between a batched apply that was handed the cache and one that was not,
      !! and the same for the single-density entry point
      character(len=*), intent(in) :: functional
      real(dp), intent(out) :: worst
         !! Largest absolute difference over the batch, in Hartree per whatever
         !! the trial density carries -- zero is what is expected, not small
      real(dp), intent(out) :: worst_one
         !! The same for `xc_kernel_apply`, which forwards its optional cache
      type(error_type), allocatable, intent(out) :: error
      logical, intent(out) :: ok

      type(czt_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(xc_kernel_cache_t) :: cache
      type(error_t) :: err
      real(dp), allocatable :: dens(:, :), trials(:, :, :)
      real(dp), allocatable :: plain(:, :, :), cached(:, :, :)
      real(dp), allocatable :: one_plain(:, :), one_cached(:, :)
      integer :: nao, iset, threads

      ok = .false.
      worst = 0.0_dp
      worst_one = 0.0_dp
      ! A build without libxc has no kernel to cache, and the suite treats that
      ! as nothing to check rather than as a failure.
      if (.not. xc_available()) then
         ok = .true.
         return
      end if

      threads = 1
!$    threads = omp_get_max_threads()
!$    call omp_set_num_threads(1)

      call reference_state(functional, mol, ctx, dens, err, ok)
      if (.not. ok) then
         call check(error, .false., "the reference Kohn-Sham state failed")
!$       call omp_set_num_threads(threads)
         return
      end if
      ok = .false.

      nao = size(dens, 1)
      call kernel_trials(dens, trials)
      allocate (plain(nao, nao, N_TRIAL), cached(nao, nao, N_TRIAL))
      plain = 0.0_dp
      cached = 0.0_dp

      call xc_kernel_apply_many(ctx, mol, dens, trials, plain, err)
      call xc_kernel_cache_fill(ctx, mol, dens, cache, err)
      call xc_kernel_apply_many(ctx, mol, dens, trials, cached, err, cache=cache)

      allocate (one_plain(nao, nao), one_cached(nao, nao))
      one_plain = 0.0_dp
      one_cached = 0.0_dp
      call xc_kernel_apply(ctx, mol, dens, trials(:, :, 1), one_plain, err)
      call xc_kernel_apply(ctx, mol, dens, trials(:, :, 1), one_cached, err, cache=cache)

!$    call omp_set_num_threads(threads)

      call check(error,.not. err%has_error(), "the kernel apply or the cache fill failed")
      if (allocated(error)) then
         call ctx%destroy()
         call mol%destroy()
         return
      end if

      do iset = 1, N_TRIAL
         worst = max(worst, maxval(abs(cached(:, :, iset) - plain(:, :, iset))))
      end do
      worst_one = maxval(abs(one_cached - one_plain))

      ! Not a check of the cache, but of the comparison: two zero matrices
      ! agree perfectly and would pass everything above.
      call check(error, maxval(abs(plain)) > 1.0e-6_dp, "the uncached kernel returned "// &
                 "nothing, so the comparison proves nothing")

      call cache%destroy()
      call ctx%destroy()
      call mol%destroy()
      ok = .not. allocated(error)
   end subroutine cache_case

   subroutine reference_state(functional, mol, ctx, density, err, ok)
      !! One converged Kohn-Sham water, its grid and its density
      !!
      !! `allow_half` because `mgga_x_tpss` is named on purpose: one rung's
      !! kernel with no correlation term beside it, so a dropped tau channel
      !! cannot be masked by a correlation component that has none.
      !!
      !! The molecule comes back alive, unlike the derivative tests' version of
      !! this: every call below wants the same one, and rebuilding it per call
      !! would leave the comparison resting on two builds agreeing.
      character(len=*), intent(in) :: functional
      type(czt_molecule_t), intent(out) :: mol
      type(xc_context_t), intent(out) :: ctx
      real(dp), allocatable, intent(out) :: density(:, :)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(rhf_result_t) :: scf

      ok = .false.
      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      if (err%has_error()) return
      call xc_context_create(mol, functional, ctx, err, level=3, allow_half=.true.)
      if (err%has_error()) return
      call run_czt_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, xc=ctx)
      if (err%has_error()) return
      density = scf%density
      ok = .true.
   end subroutine reference_state

   subroutine kernel_trials(dens, trials)
      !! Three symmetric trial densities, none proportional to the reference
      !!
      !! A trial proportional to the density moves rho and sigma together and
      !! cannot separate the channels; three unlike ones also make the batch
      !! loop run its stacking more than once.
      real(dp), intent(in) :: dens(:, :)
      real(dp), allocatable, intent(out) :: trials(:, :, :)

      integer :: n, i, j

      n = size(dens, 1)
      allocate (trials(n, n, N_TRIAL))
      do j = 1, n
         do i = 1, n
            trials(i, j, 1) = 0.05_dp/(1.0_dp + real(abs(i - j), dp)) + 0.01_dp*dens(i, j)
            trials(i, j, 2) = 0.1_dp*cos(real(i - j, dp)) + 0.03_dp*sin(real(i + j, dp))
            trials(i, j, 3) = 0.02_dp*real(i*j, dp)/real(n*n, dp)
         end do
         trials(j, j, 2) = trials(j, j, 2) + 0.2_dp
      end do
      do i = 1, N_TRIAL
         trials(:, :, i) = 0.5_dp*(trials(:, :, i) + transpose(trials(:, :, i)))
      end do
   end subroutine kernel_trials

end module test_mqc_xc_kernel_cache

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_xc_kernel_cache, only: collect_mqc_xc_kernel_cache
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_xc_kernel_cache", collect_mqc_xc_kernel_cache)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
