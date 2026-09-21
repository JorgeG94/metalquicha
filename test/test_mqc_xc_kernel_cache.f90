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
   use mqc_memory, only: set_memory_budget
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available, &
                         xc_kernel_apply, xc_kernel_apply_many, &
                         xc_kernel_cache_t, xc_kernel_cache_fill, &
                         xc_kernel_apply_uks_many, xc_kernel_cache_uks_t, &
                         xc_kernel_cache_uks_fill
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

   real(dp), parameter :: TOL_POLARISED = 1.0e-13_dp
      !! The polarised kernel against the two restricted ones it contains.
      !!
      !! Not zero, unlike the cache comparisons: the restricted side asks
      !! libxc for the unpolarised functional's second derivative and this
      !! side asks for the polarised one at `rho_a = rho_b`, which are two
      !! evaluations of the same analytic quantity and agree to the last few
      !! bits rather than to all of them.

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
                  new_unittest("cached_hybrid_kernel_is_bit_identical", test_hybrid), &
                  new_unittest("an_unfilled_cache_is_refused", test_unfilled), &
                  new_unittest("a_cache_from_another_rung_is_refused", test_wrong_rung), &
                  new_unittest("an_over_budget_cache_is_declined", test_over_budget), &
                  new_unittest("cached_triplet_kernel_is_bit_identical", test_triplet), &
                  new_unittest("a_singlet_only_cache_refuses_a_triplet", test_triplet_refused), &
                  new_unittest("the_polarised_gga_kernel_is_the_two_restricted_ones", &
                               test_uks_pbe), &
                  new_unittest("the_polarised_hybrid_kernel_is_the_two_restricted_ones", &
                               test_uks_b3lyp), &
                  new_unittest("the_polarised_lda_kernel_is_the_two_restricted_ones", &
                               test_uks_lda) &
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
      call check(error, worst == 0.0_dp, "the cached B3LYP kernel differs from "// &
                 "the evaluated one")
      if (allocated(error)) return
      call check(error, worst_one == 0.0_dp, "the cached single-density B3LYP "// &
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

   subroutine test_triplet(error)
      !! The triplet kernel through the cache is the triplet kernel without it
      !!
      !! Two paths reach the polarised evaluation: `xc_kernel_cache_fill` makes
      !! it once over the whole grid, and `kernel_apply_batch` makes it per
      !! block when there is no cache. A response solve only ever takes the
      !! first, so without this the second is code nothing runs -- and it is
      !! the one the first is supposed to be a cache *of*.
      !!
      !! The second check is not about caching at all. A triplet kernel that
      !! was quietly the singlet one would agree with itself through both
      !! paths and pass everything above; what says the polarised evaluation
      !! happened is that the two differ, which for PBE they do by three
      !! orders more than this floor.
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: worst, manifold_gap
      logical :: ok

      ! Redundant now that the helpers report a libxc-less build as `.not. ok`,
      ! and kept because it states at the point of use what the whole file
      ! depends on: no assertion below runs on a build with no kernel to
      ! evaluate. Reading those zeros as measurements is what turned three of
      ! the eleven CI jobs red -- the three built with MQC_ENABLE_LIBXC=OFF.
      if (.not. xc_available()) return

      call triplet_case("pbe", worst, manifold_gap, error, ok)
      if (allocated(error) .or. .not. ok) return
      call check(error, worst == 0.0_dp, "the cached triplet kernel differs from "// &
                 "the evaluated one")
      if (allocated(error)) return
      call check(error, manifold_gap > 1.0e-6_dp, "the triplet kernel is indistinguishable "// &
                 "from the singlet one, so the polarised evaluation did not happen")
   end subroutine test_triplet

   subroutine test_triplet_refused(error)
      !! A cache filled for singlets only cannot serve a triplet contraction
      !!
      !! It is the right shape and the right grid, and read as though it were
      !! a triplet cache it would contract the singlet coefficients and return
      !! a converged spectrum of the wrong manifold.
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
      call xc_kernel_cache_fill(ctx, mol, dens, cache, err)
      call xc_kernel_apply_many(ctx, mol, dens, trials, out, err, cache=cache, &
                                triplet=.true.)
      call check(error, err%has_error(), "a singlet-only cache was read as a triplet one")
      call cache%destroy()
      call ctx%destroy()
      call mol%destroy()
   end subroutine test_triplet_refused

   subroutine test_over_budget(error)
      !! A fill that will not fit declines rather than allocating over the
      !! budget: no error, an unfilled cache, and a consumer left on the path
      !! it would have taken had nobody asked for a cache at all. Forced by
      !! fixing the budget at a byte, which is what a deck's `system.memory_gb`
      !! reaches the same code through.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(xc_kernel_cache_t) :: cache
      type(error_t) :: err
      real(dp), allocatable :: dens(:, :)
      logical :: ok

      if (.not. xc_available()) return
      call reference_state("pbe", mol, ctx, dens, err, ok)
      if (.not. ok) then
         call check(error, .false., "the reference Kohn-Sham state failed")
         return
      end if

      call set_memory_budget(1.0e-9_dp)     ! one byte
      call xc_kernel_cache_fill(ctx, mol, dens, cache, err)
      call set_memory_budget(-1.0_dp)       ! back to the machine

      call check(error,.not. err%has_error(), "a declined kernel cache raised an "// &
                 "error, where its consumers expect to be left on the uncached path")
      if (.not. allocated(error)) then
         call check(error,.not. cache%filled, "a kernel cache was filled over its "// &
                    "budget")
      end if
      call cache%destroy()
      call ctx%destroy()
      call mol%destroy()
   end subroutine test_over_budget

   subroutine test_wrong_rung(error)
      !! A cache filled against a GGA and handed to a meta-GGA contraction is
      !! three channels short of the answer, and the shapes do not say so: the
      !! two contexts are the same grid over the same molecule, so every array
      !! is exactly the right length. Refused on the recorded rung, not on a
      !! size.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol_gga, mol_mgga
      type(xc_context_t) :: ctx_gga, ctx_mgga
      type(xc_kernel_cache_t) :: cache
      type(error_t) :: err
      real(dp), allocatable :: dens_gga(:, :), dens_mgga(:, :)
      real(dp), allocatable :: trials(:, :, :), out(:, :, :)
      integer :: nao
      logical :: ok

      if (.not. xc_available()) return

      call reference_state("pbe", mol_gga, ctx_gga, dens_gga, err, ok)
      if (.not. ok) then
         call check(error, .false., "the reference GGA Kohn-Sham state failed")
         return
      end if
      call reference_state("mgga_x_tpss", mol_mgga, ctx_mgga, dens_mgga, err, ok)
      if (.not. ok) then
         call check(error, .false., "the reference meta-GGA Kohn-Sham state failed")
         call ctx_gga%destroy()
         call mol_gga%destroy()
         return
      end if

      call xc_kernel_cache_fill(ctx_gga, mol_gga, dens_gga, cache, err)
      call check(error,.not. err%has_error(), "the GGA cache fill failed")
      if (allocated(error)) then
         call cleanup()
         return
      end if
      call check(error, cache%filled, "the GGA cache was declined, so the refusal "// &
                 "below would fire for the wrong reason")
      if (allocated(error)) then
         call cleanup()
         return
      end if

      ! The same grid, so the cache's arrays are the length this context wants.
      call check(error, cache%n_points == ctx_mgga%grid%n_points, "the two contexts "// &
                 "do not share a grid, so a length check alone would catch this")
      if (allocated(error)) then
         call cleanup()
         return
      end if

      nao = size(dens_mgga, 1)
      allocate (trials(nao, nao, 1), out(nao, nao, 1))
      trials(:, :, 1) = dens_mgga
      out = 0.0_dp
      call xc_kernel_apply_many(ctx_mgga, mol_mgga, dens_mgga, trials, out, err, cache=cache)
      call check(error, err%has_error(), "a GGA kernel cache was accepted by a "// &
                 "meta-GGA contraction")
      call cleanup()

   contains

      subroutine cleanup()
         call cache%destroy()
         call ctx_gga%destroy()
         call mol_gga%destroy()
         call ctx_mgga%destroy()
         call mol_mgga%destroy()
      end subroutine cleanup
   end subroutine test_wrong_rung

   subroutine triplet_case(functional, worst, manifold_gap, error, ok)
      !! The triplet kernel both ways, and how far it is from the singlet one
      character(len=*), intent(in) :: functional
      real(dp), intent(out) :: worst
         !! Largest absolute difference between the cached and the uncached
         !! triplet contraction. Zero is what is expected, not small.
      real(dp), intent(out) :: manifold_gap
         !! Largest absolute difference between the triplet contraction and
         !! the singlet one, which says the polarised evaluation ran
      type(error_type), allocatable, intent(out) :: error
      logical, intent(out) :: ok
         !! `.true.` means nothing failed, **not** that the outputs are
         !! meaningful: a build without libxc returns `.true.` with both
         !! differences left at zero. Guard on `xc_available()` before
         !! asserting anything about them.

      type(czt_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(xc_kernel_cache_t) :: cache
      type(error_t) :: err
      real(dp), allocatable :: dens(:, :), trials(:, :, :)
      real(dp), allocatable :: plain(:, :, :), cached(:, :, :), singlet(:, :, :)
      integer :: nao, iset, threads

      ok = .false.
      worst = 0.0_dp
      manifold_gap = 0.0_dp
      if (.not. xc_available()) return

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
      allocate (plain(nao, nao, N_TRIAL), cached(nao, nao, N_TRIAL), &
                singlet(nao, nao, N_TRIAL))
      plain = 0.0_dp
      cached = 0.0_dp
      singlet = 0.0_dp

      call xc_kernel_apply_many(ctx, mol, dens, trials, plain, err, triplet=.true.)
      call xc_kernel_cache_fill(ctx, mol, dens, cache, err, triplet=.true.)
      call xc_kernel_apply_many(ctx, mol, dens, trials, cached, err, cache=cache, &
                                triplet=.true.)
      call xc_kernel_apply_many(ctx, mol, dens, trials, singlet, err, cache=cache)

!$    call omp_set_num_threads(threads)

      call check(error,.not. err%has_error(), "the triplet kernel apply or its cache "// &
                 "fill failed: "//err%get_message())
      if (allocated(error)) then
         call ctx%destroy()
         call mol%destroy()
         return
      end if

      do iset = 1, N_TRIAL
         worst = max(worst, maxval(abs(cached(:, :, iset) - plain(:, :, iset))))
         manifold_gap = max(manifold_gap, &
                            maxval(abs(cached(:, :, iset) - singlet(:, :, iset))))
      end do

      call cache%destroy()
      call ctx%destroy()
      call mol%destroy()
      ok = .not. allocated(error)
   end subroutine triplet_case

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
      ! as nothing to check rather than as a failure. `ok` stays false: it says
      ! the outputs are measurements, and here there was nothing to measure.
      ! The zeros happen to pass this helper's assertions, which are all
      ! "the difference is zero" -- that made the skip look harmless here and
      ! cost three CI jobs where the assertion was instead "the gap is large".
      if (.not. xc_available()) return

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

   subroutine test_uks_pbe(error)
      !! PBE: the rung an unrestricted response solve actually runs on
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: singlet, triplet, cached, spread
      logical :: ok

      call uks_case("pbe", singlet, triplet, cached, spread, error, ok)
      if (allocated(error) .or. .not. ok) return
      call check(error, singlet < TOL_POLARISED, "the polarised kernel on an equal "// &
                 "spin perturbation is not the restricted singlet kernel")
      if (allocated(error)) return
      call check(error, triplet < TOL_POLARISED, "the polarised kernel on an opposite "// &
                 "spin perturbation is not the restricted triplet kernel")
      if (allocated(error)) return
      call check(error, cached == 0.0_dp, "the cached polarised kernel differs from "// &
                 "the evaluated one")
      if (allocated(error)) return
      ! The two restricted kernels have to be far apart for the two checks
      ! above to mean anything: agreeing with both would otherwise be one
      ! statement, not two.
      call check(error, spread > 1.0e-3_dp, "the singlet and triplet kernels are "// &
                 "too close together for this comparison to say anything")
   end subroutine test_uks_pbe

   subroutine test_uks_b3lyp(error)
      !! B3LYP: five libxc components, summed with weights, twice over
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: singlet, triplet, cached, spread
      logical :: ok

      call uks_case("b3lyp", singlet, triplet, cached, spread, error, ok)
      if (allocated(error) .or. .not. ok) return
      call check(error, singlet < TOL_POLARISED, "the polarised B3LYP kernel is not "// &
                 "the restricted singlet one on an equal spin perturbation")
      if (allocated(error)) return
      call check(error, triplet < TOL_POLARISED, "the polarised B3LYP kernel is not "// &
                 "the restricted triplet one on an opposite spin perturbation")
      if (allocated(error)) return
      call check(error, cached == 0.0_dp, "the cached polarised B3LYP kernel differs "// &
                 "from the evaluated one")
      if (allocated(error)) return
      call check(error, spread > 1.0e-3_dp, "the singlet and triplet B3LYP kernels "// &
                 "are too close together for this comparison to say anything")
   end subroutine test_uks_b3lyp

   subroutine test_uks_lda(error)
      !! LDA: only `v2rho2` survives, so this is the rung where a wrong
      !! `v2rhosigma` or `v2sigma2` index cannot hide the answer
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: singlet, triplet, cached, spread
      logical :: ok

      call uks_case("svwn", singlet, triplet, cached, spread, error, ok)
      if (allocated(error) .or. .not. ok) return
      call check(error, singlet < TOL_POLARISED, "the polarised LDA kernel is not "// &
                 "the restricted singlet one")
      if (allocated(error)) return
      call check(error, triplet < TOL_POLARISED, "the polarised LDA kernel is not "// &
                 "the restricted triplet one")
      if (allocated(error)) return
      call check(error, cached == 0.0_dp, "the cached polarised LDA kernel differs "// &
                 "from the evaluated one")
      if (allocated(error)) return
      call check(error, spread > 1.0e-3_dp, "the singlet and triplet LDA kernels are "// &
                 "too close together for this comparison to say anything")
   end subroutine test_uks_lda

   subroutine uks_case(functional, singlet, triplet, cached, spread, error, ok)
      !! The polarised kernel against the two restricted ones it contains
      !!
      !! On a **closed shell**, where all three are defined and related
      !! exactly. Perturb the two spins together, `drho_a = drho_b = drho/2`,
      !! and the polarised kernel's alpha response is the restricted singlet
      !! one; perturb them oppositely and it is the triplet one. Those are the
      !! two combinations the restricted code computes directly, out of the
      !! unpolarised functional and out of `triplet_component`; this routine
      !! computes them out of the full spin-resolved second derivative, and
      !! nothing but the physics makes the three agree.
      !!
      !! The disagreement that survives is libxc's own: the restricted side
      !! evaluates the unpolarised functional and this side the polarised one
      !! at `rho_a = rho_b`. The same number by two routes, to `TOL_POLARISED`.
      character(len=*), intent(in) :: functional
      real(dp), intent(out) :: singlet
         !! Largest disagreement with the restricted singlet kernel
      real(dp), intent(out) :: triplet
         !! The same against the restricted triplet kernel
      real(dp), intent(out) :: cached
         !! Cached against uncached, which is expected to be zero
      real(dp), intent(out) :: spread
         !! How far apart the two restricted kernels are, so the two
         !! comparisons above are two statements and not one
      type(error_type), allocatable, intent(out) :: error
      logical, intent(out) :: ok
         !! Whether the four outputs above are measurements rather than the
         !! zeros they start as. False both when libxc is absent and when the
         !! setup below failed; the callers return on it either way.

      type(czt_molecule_t) :: mol
      type(xc_context_t) :: ctx, pol
      type(xc_kernel_cache_uks_t) :: cache
      type(error_t) :: err
      real(dp), allocatable :: dens(:, :), trials(:, :, :), half(:, :)
      real(dp), allocatable :: sing_r(:, :, :), trip_r(:, :, :)
      real(dp), allocatable :: pa(:, :, :), pb(:, :, :), qa(:, :, :), qb(:, :, :)
      real(dp), allocatable :: ca(:, :, :), cb(:, :, :)
      real(dp), allocatable :: plus_a(:, :, :), plus_b(:, :, :)
      real(dp), allocatable :: minus_a(:, :, :), minus_b(:, :, :)
      integer :: nao, iset, threads

      ok = .false.
      singlet = 0.0_dp
      triplet = 0.0_dp
      cached = 0.0_dp
      spread = 0.0_dp

      ! Saying `.true.` here told the callers that four zeros were a result.
      ! They return on `.not. ok` and were right to; what they got instead was
      ! a `spread` of zero, which fails "the two kernels are far enough apart
      ! for this comparison to say anything" on every job built without libxc.
      ! The flag means the outputs are real, so with nothing to measure it is
      ! false -- the same answer as a setup that failed, and the same handling.
      if (.not. xc_available()) return

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
      ! The same functional and the same grid level, spin-polarised: libxc
      ! fixes the channel at initialisation, so the unrestricted kernel needs
      ! its own handles rather than the restricted context's.
      call xc_context_create(mol, functional, pol, err, level=3, allow_half=.true., &
                             polarized=.true.)
      if (err%has_error()) then
         call check(error, .false., "the polarised context failed: "//err%get_message())
!$       call omp_set_num_threads(threads)
         call ctx%destroy()
         call mol%destroy()
         return
      end if

      nao = size(dens, 1)
      call kernel_trials(dens, trials)
      half = 0.5_dp*dens
      allocate (sing_r(nao, nao, N_TRIAL), trip_r(nao, nao, N_TRIAL))
      allocate (plus_a(nao, nao, N_TRIAL), plus_b(nao, nao, N_TRIAL))
      allocate (minus_a(nao, nao, N_TRIAL), minus_b(nao, nao, N_TRIAL))
      allocate (pa(nao, nao, N_TRIAL), pb(nao, nao, N_TRIAL))
      allocate (qa(nao, nao, N_TRIAL), qb(nao, nao, N_TRIAL))
      allocate (ca(nao, nao, N_TRIAL), cb(nao, nao, N_TRIAL))
      sing_r = 0.0_dp
      trip_r = 0.0_dp
      plus_a = 0.0_dp
      plus_b = 0.0_dp
      minus_a = 0.0_dp
      minus_b = 0.0_dp
      ca = 0.0_dp
      cb = 0.0_dp
      ! `drho_a = drho_b = drho/2` for the singlet combination, and
      ! `drho_a = -drho_b = drho/2` for the triplet one.
      pa = 0.5_dp*trials
      pb = pa
      qa = pa
      qb = -pa

      call xc_kernel_apply_many(ctx, mol, dens, trials, sing_r, err)
      call xc_kernel_apply_many(ctx, mol, dens, trials, trip_r, err, triplet=.true.)
      call xc_kernel_apply_uks_many(pol, mol, half, half, pa, pb, plus_a, plus_b, err)
      call xc_kernel_apply_uks_many(pol, mol, half, half, qa, qb, minus_a, minus_b, err)
      call xc_kernel_cache_uks_fill(pol, mol, half, half, cache, err)
      call xc_kernel_apply_uks_many(pol, mol, half, half, pa, pb, ca, cb, err, &
                                    cache=cache)

!$    call omp_set_num_threads(threads)

      call check(error,.not. err%has_error(), "a kernel apply or cache fill failed: "// &
                 err%get_message())
      if (allocated(error)) then
         call pol%destroy()
         call ctx%destroy()
         call mol%destroy()
         return
      end if

      do iset = 1, N_TRIAL
         singlet = max(singlet, maxval(abs(plus_a(:, :, iset) - sing_r(:, :, iset))))
         ! And the beta half, which a routine that filled only alpha would
         ! leave at zero.
         singlet = max(singlet, maxval(abs(plus_b(:, :, iset) - sing_r(:, :, iset))))
         triplet = max(triplet, maxval(abs(minus_a(:, :, iset) - trip_r(:, :, iset))))
         triplet = max(triplet, maxval(abs(minus_b(:, :, iset) + trip_r(:, :, iset))))
         cached = max(cached, maxval(abs(ca(:, :, iset) - plus_a(:, :, iset))))
         cached = max(cached, maxval(abs(cb(:, :, iset) - plus_b(:, :, iset))))
         spread = max(spread, maxval(abs(sing_r(:, :, iset) - trip_r(:, :, iset))))
      end do

      call cache%destroy()
      call pol%destroy()
      call ctx%destroy()
      call mol%destroy()
      ok = .not. allocated(error)
   end subroutine uks_case

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
