module mqc_libcint_direct
   !! Direct Fock construction: integrals consumed as produced, never stored
   !!
   !! The in-core path in `mqc_libcint_integrals` holds all n^4 integrals and
   !! contracts them each iteration. That is 800 MB at a hundred basis functions
   !! and 65 GB at three hundred, so it is a reference implementation and not a
   !! method. This builds the Fock matrix without ever forming the tensor.
   !!
   !! Three things make it work, following Huang, Sherrill and Chow
   !! [JCP 152, 024122 (2020)], Algorithm 1:
   !!
   !! **Permutational symmetry.** (MN|PQ) is invariant under swapping M and N,
   !! swapping P and Q, and swapping the pairs -- eightfold. Only quartets with
   !! M >= N, P >= Q and (MN) >= (PQ) are computed, and each contributes to six
   !! Fock blocks rather than one.
   !!
   !! **Schwarz screening.** |(MN|PQ)| <= Q_MN Q_PQ with Q_MN = sqrt(max|(MN|MN)|),
   !! by Cauchy-Schwarz. Quartets whose bound is below threshold are skipped
   !! without being computed. This is what turns the formal n^4 into roughly n^2
   !! for extended systems, and it is the reason recomputing integrals every
   !! iteration is cheaper than streaming them from memory.
   !!
   !! **libcint's optimizer.** CINTOpt precomputes per-shell-pair data that is
   !! reused across quartet calls. It is created once per Fock build and costs
   !! nothing to use.
   !!
   !! The degeneracy factors deserve care. A quartet where M == N still has both
   !! (mu nu|..) and (nu mu|..) inside its own block -- the block enumerates all
   !! function pairs -- so the permutation is already covered and must not be
   !! counted again. That is why the factors test shell equality, not function
   !! equality, and why getting it wrong produces a Fock matrix that is wrong by
   !! a factor of two on its diagonal blocks only. `check_direct` compares
   !! against the in-core build elementwise, which catches exactly that.
   use pic_types, only: dp, int64, int_index
   use pic_sorting, only: sort_index
   use pic_timer, only: timer_type
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim, &
                                    two_electron_block, two_electron_optimizer, &
                                    eri_shell_table_t, eri_shell_table, &
                                    eri_schwarz_collapse
   use libcint_fortran, only: libcint_del_optimizer, LIBCINT_PTR_RANGE_OMEGA
   implicit none
   private

   public :: schwarz_bounds
   public :: shell_density_max
   public :: build_fock_direct
   public :: build_fock_direct_many
   public :: build_fock_direct_nosym
   public :: build_fock_direct_uhf
   public :: direct_stats_t
   public :: DEFAULT_SCREEN_TOL

   !> Quartets whose Schwarz bound falls below this are skipped.
   !>
   !> The GTFock paper uses 1e-11 on shell quartets. That is tighter than the
   !> 1e-10 often seen and costs little, since the count of surviving quartets
   !> is insensitive to the threshold over a decade or so.
   real(dp), parameter :: DEFAULT_SCREEN_TOL = 1.0e-11_dp

   type :: direct_stats_t
      !! What a Fock build did, so screening can be reported rather than assumed
      integer(int64) :: quartets_total = 0     !! Unique quartets before screening
      integer(int64) :: quartets_computed = 0  !! Quartets actually handed to libcint
      integer(int64) :: quartets_screened = 0  !! Skipped, either test
      integer(int64) :: screened_schwarz = 0
         !! Skipped because the integral itself is negligible. Depends on the
         !! basis and the geometry only, so this number is the same at every
         !! iteration of an SCF.
      real(dp) :: thread_imbalance = 1.0_dp
         !! Slowest thread's time on the quartet loop divided by the mean.
         !!
         !! One is perfect balance. This is what a work-ordering change has to
         !! move, and it is measurable where wall time is not: it is a ratio
         !! within a single run, so the contention that swamps a before/after
         !! comparison on a shared node largely divides out.
      integer(int64) :: screened_density = 0
         !! Skipped because the *contribution* is negligible although the
         !! integral is not -- the extra reach the density weighting buys.
         !!
         !! Counted separately because it is the only honest way to measure that
         !! weighting on a shared machine: run-to-run wall time on a contended
         !! node varies by more than the effect, but these counts are exactly
         !! reproducible.
   contains
      procedure :: screened_fraction => stats_screened_fraction
   end type direct_stats_t

contains

   pure function stats_screened_fraction(this) result(fraction)
      !! Fraction of unique quartets that screening removed
      class(direct_stats_t), intent(in) :: this
      real(dp) :: fraction

      if (this%quartets_total <= 0) then
         fraction = 0.0_dp
      else
         fraction = real(this%quartets_screened, dp)/real(this%quartets_total, dp)
      end if
   end function stats_screened_fraction

   subroutine schwarz_bounds(mol, bounds, error)
      !! Q_MN = sqrt(max |(MN|MN)|) for every shell pair
      !!
      !! Costs n_shells^2 quartet evaluations, done once per geometry -- the
      !! bounds depend only on the basis and the positions, not the density, so
      !! a whole SCF reuses one set.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), allocatable, intent(out) :: bounds(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: buf(:)
      integer :: shls(4)
      integer :: ish, jsh, di, dj, ret, block_max
      real(dp) :: largest

      block_max = 1
      do ish = 1, mol%nbas
         block_max = max(block_max, mol%shell_offset(ish + 1) - mol%shell_offset(ish))
      end do

      allocate (bounds(mol%nbas, mol%nbas))
      allocate (buf(block_max**4))
      bounds = 0.0_dp

      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         do jsh = 1, ish
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            shls = [ish - 1, jsh - 1, ish - 1, jsh - 1]
            ret = two_electron_block(mol%cartesian, buf, shls, mol%atm, mol%natm, &
                                     mol%bas, mol%nbas, mol%env)
            if (ret == 0) then
               largest = 0.0_dp
            else
               ! (MN|MN) is positive semidefinite on its diagonal, but take the
               ! absolute value anyway: the bound is what matters, not the sign.
               largest = maxval(abs(buf(1:(di*dj)**2)))
            end if
            bounds(ish, jsh) = sqrt(largest)
            bounds(jsh, ish) = bounds(ish, jsh)
         end do
      end do

      if (all(bounds <= 0.0_dp)) then
         call error%set(ERROR_VALIDATION, "Schwarz bounds are all zero; the basis is empty or broken")
      end if

      deallocate (buf)
   end subroutine schwarz_bounds

   integer function omp_threads() result(n)
      !! Threads that ran the last parallel region, or one without OpenMP
!$    use omp_lib, only: omp_get_max_threads
      n = 1
!$    n = omp_get_max_threads()
   end function omp_threads

   subroutine pair_work_order(pair_i, pair_j, dims, order)
      !! Task indices sorted by descending cost
      !!
      !! Task `ij` runs `kl = 1..ij`, so its cost is the block size of pair `ij`
      !! times the block sizes of every pair at or before it -- growing roughly
      !! linearly in `ij`, but not monotonically, because a pair of two d shells
      !! is thirty-six function pairs where two s shells are one.
      !!
      !! `schedule(dynamic)` hands tasks out in order, so ascending order
      !! dispatches the most expensive task last and finishes with every other
      !! thread idle behind it. Reversing fixes the trend but not the variation
      !! within it; sorting on the actual cost fixes both, and `pic_sorting`'s
      !! `sort_index` is an introsort that returns the permutation rather than
      !! the sorted values, which is exactly the question being asked.
      !!
      !! The cost is the unscreened function-quartet count. Screening will remove
      !! some of it unevenly, so this is an estimate -- but it is the estimate
      !! available before the work is done, and the ordering only has to be
      !! roughly right for the tail to disappear.
      integer, intent(in) :: pair_i(:), pair_j(:), dims(:)
      integer, allocatable, intent(out) :: order(:)

      integer(int64), allocatable :: cost(:)
      integer(int_index), allocatable :: perm(:)
      integer(int64) :: cumulative, block
      integer :: k, npair

      npair = size(pair_i)
      allocate (cost(npair), perm(npair), order(npair))

      cumulative = 0_int64
      do k = 1, npair
         block = int(dims(pair_i(k)), int64)*int(dims(pair_j(k)), int64)
         cumulative = cumulative + block
         cost(k) = block*cumulative
      end do

      call sort_index(cost, perm, reverse=.true.)
      order = int(perm)
   end subroutine pair_work_order

   subroutine shell_density_max(mol, density, dsh)
      !! The largest |D| in each shell-pair block
      !!
      !! The Schwarz bound says how big an integral can be. It says nothing about
      !! whether that integral matters, and by the middle of an SCF most of them
      !! do not: a quartet multiplies six density elements, and if all six are
      !! negligible the quartet contributes nothing however large it is.
      !!
      !! Blocked to shells rather than kept per function because the screening
      !! decision is made per shell quartet -- one number per pair is all the test
      !! can use, and it has to be the largest in the block or the bound would not
      !! be one.
      !!
      !! Rebuilt every Fock build, unlike the Schwarz bounds: this one is a
      !! property of the density and the density is what changes. It costs one
      !! pass over an n-by-n matrix against a quartet loop that is the rest of the
      !! iteration, so the cost does not show up.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: dsh(:, :)

      integer, allocatable :: offs(:), dims(:)
      integer :: ish

      allocate (offs(mol%nbas), dims(mol%nbas))
      do ish = 1, mol%nbas
         offs(ish) = mol%shell_offset(ish)
         dims(ish) = shell_dim(mol%cartesian, ish - 1, mol%bas)
      end do
      call block_density_max(density, mol%nbas, offs, dims, dsh)
   end subroutine shell_density_max

   subroutine block_density_max(density, nbas, offs, dims, dsh)
      !! `shell_density_max` over an explicit shell blocking
      !!
      !! Split out so the Fock builders can block the density to whatever shell
      !! table their quartet loop runs over -- the fused-sp view when the
      !! molecule carries one -- while `shell_density_max` keeps answering per
      !! split shell for any caller blocked that way. The SCF gradient, once
      !! the reason this note said "never fused", now blocks against the view
      !! it loops over, inline in `two_electron_deriv`.
      real(dp), intent(in) :: density(:, :)
      integer, intent(in) :: nbas
      integer, intent(in) :: offs(:)   !! First AO of each shell, 0-based
      integer, intent(in) :: dims(:)   !! Functions per shell
      real(dp), allocatable, intent(out) :: dsh(:, :)

      integer :: ish, jsh, oi, oj, di, dj

      allocate (dsh(nbas, nbas))
      dsh = 0.0_dp
      do ish = 1, nbas
         oi = offs(ish)
         di = dims(ish)
         do jsh = 1, ish
            oj = offs(jsh)
            dj = dims(jsh)
            dsh(ish, jsh) = maxval(abs(density(oi + 1:oi + di, oj + 1:oj + dj)))
            dsh(jsh, ish) = dsh(ish, jsh)
         end do
      end do
   end subroutine block_density_max

   pure function density_weight(dsh, s1, s2, s3, s4, jf, kq, deg) result(denmax)
      !! Largest density contribution this quartet can make, per unit integral
      !!
      !! The six elements are the six the Fock updates read, and the weights are
      !! the ones they are multiplied by, so `bound * denmax` is an upper bound on
      !! the largest change any one Fock element can see from this quartet. That
      !! is what the screening test should compare against a tolerance -- the
      !! Schwarz bound alone compares the integral, which is not the quantity
      !! anyone cares about being small.
      !!
      !! `deg` belongs in here. The permutational factor multiplies the
      !! contribution, so leaving it out would make the bound smaller than the
      !! thing it bounds and the test could discard a quartet that mattered.
      !! Including it costs a little screening and buys rigour.
      real(dp), intent(in) :: dsh(:, :)
      integer, intent(in) :: s1, s2, s3, s4
      real(dp), intent(in) :: jf, kq, deg
      real(dp) :: denmax

      denmax = deg*max(jf*dsh(s1, s2), jf*dsh(s3, s4), &
                       kq*dsh(s1, s3), kq*dsh(s1, s4), &
                       kq*dsh(s2, s3), kq*dsh(s2, s4))
   end function density_weight

   pure function pair_degeneracy(s1, s2, s3, s4) result(deg)
      !! The eightfold permutational weight this quartet carries
      !!
      !! Shell equality, not function equality: a block with s1 == s2 already
      !! contains both orderings of its function pair, so the permutation is
      !! covered and must not be counted again. Getting it wrong produces a Fock
      !! matrix wrong by a factor of two on its diagonal blocks only.
      integer, intent(in) :: s1, s2, s3, s4
      real(dp) :: deg

      deg = 1.0_dp
      if (s1 /= s2) deg = deg*2.0_dp
      if (s3 /= s4) deg = deg*2.0_dp
      if (.not. (s1 == s3 .and. s2 == s4)) deg = deg*2.0_dp
   end function pair_degeneracy

   subroutine build_fock_direct(mol, h, density, bounds, fock, stats, error, screen_tol, &
                                k_scale, j_scale, omega, density_screen)
      !! F = H + J - K/2, without forming the integral tensor
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h(:, :)          !! Core Hamiltonian
      real(dp), intent(in) :: density(:, :)    !! D = 2 C_occ C_occ^T
      real(dp), intent(in) :: bounds(:, :)     !! From `schwarz_bounds`
      real(dp), intent(out) :: fock(:, :)
      type(direct_stats_t), intent(out) :: stats
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: screen_tol
      real(dp), intent(in), optional :: k_scale
         !! Fraction of exact exchange to keep. One is Hartree-Fock and the
         !! default; zero is pure density-functional exchange; a hybrid is between.
         !! Present so a Kohn-Sham build needs no second routine and cannot drift
         !! out of step with this one.
      real(dp), intent(in), optional :: j_scale
         !! Fraction of the Coulomb term, default one. Zero is what a
         !! range-separated hybrid's *second* pass wants: it needs a long-range
         !! exchange matrix and nothing else, and asking for it here rather than in
         !! an exchange-only routine avoids a second copy of the quartet loop.
      real(dp), intent(in), optional :: omega
         !! Range-separation parameter. Zero, the default, is the full Coulomb
         !! kernel. Positive gives the long-range erf-attenuated one and negative
         !! the short-range complement -- libcint switches on this through
         !! `env(PTR_RANGE_OMEGA)`, so the same entry points and the same optimizer
         !! serve both and no new integral code is needed.
         !!
         !! The Schwarz bounds passed in are the full-kernel ones and stay valid
         !! here without being rebuilt: erf(omega r)/r <= 1/r pointwise, so an
         !! attenuated quartet is never larger than the bound screening it, and
         !! the screening can only become conservative rather than wrong.
      logical, intent(in), optional :: density_screen
         !! Weight the Schwarz bound by the density it multiplies before screening
         !! (default true). An SCF wants it: a quartet whose six density elements
         !! are all negligible contributes nothing however large its integral, so
         !! dropping it is free accuracy. A CPHF response build must pass false.
         !! Its `density` is a Krylov trial vector the solver drives towards zero,
         !! so a screen keyed on the density's magnitude tightens as the solve
         !! proceeds -- the operator stops being the same linear map from one
         !! matvec to the next and the iteration cannot converge. False recovers
         !! the plain Schwarz screen, which depends on the basis alone and leaves
         !! the operator fixed. Bit-for-bit equal to `build_fock_direct_many` with
         !! one density, which is the invariant the response path relies on.

      real(dp), allocatable :: buf(:), g(:, :), g_local(:, :), d_half(:, :)
      real(dp), allocatable :: dsh(:, :)
      real(dp), allocatable :: env_local(:)
      real(dp), allocatable :: bq(:, :)
      type(eri_shell_table_t) :: tab
      real(dp) :: jf
      type(c_ptr) :: opt
      integer :: s1, s2, s3, s4
      integer :: d1, d2, d3, d4, o1, o2, o3, o4
      integer :: shls(4)
      integer :: f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, block_max, n
      integer :: ij, kl, npair, ipair
      integer, allocatable :: pair_i(:), pair_j(:), dims(:), offs(:), order(:)
      integer :: itask
      type(timer_type) :: thread_clock
      real(dp) :: t_thread, t_thread_max, t_thread_sum
      integer(int64) :: n_total, n_computed, n_screened, n_schwarz, n_density
      real(dp) :: schwarz
      real(dp) :: tol, deg, value, scaled
      real(dp) :: kq
      logical :: weight_density

      n = mol%nao
      if (size(h, 1) /= n .or. size(density, 1) /= n .or. size(fock, 1) /= n) then
         call error%set(ERROR_VALIDATION, "direct Fock: matrix dimensions do not match the basis")
         return
      end if

      tol = DEFAULT_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol

      weight_density = .true.
      if (present(density_screen)) weight_density = density_screen

      ! The shells the quartet loop runs over: the fused-sp view when the
      ! molecule carries one, its split shells otherwise. libfint's `int2e` is
      ! the one driver that reads a fused L shell, and it is exactly what this
      ! loop calls, so the view is safe here and nowhere looser -- and since
      ! libfint learned to carry a derivative's tensor component through an L
      ! shell, the SCF gradient's int2e_ip1 loop takes the view too. The
      ! Schwarz bounds arrive per split shell either way, straight from
      ! `schwarz_bounds`, and are collapsed to match. Dimensions
      ! and offsets are copied out up front: both are needed inside the
      ! parallel region, and looking them up there would mean every thread
      ! reaching into `bas` for something that does not change.
      call eri_shell_table(mol, tab)
      call eri_schwarz_collapse(mol, bounds, bq)
      dims = tab%dims
      offs = tab%offs(1:tab%nbas)
      block_max = tab%block_max

      ! The quartet loop, flattened onto one index so it can be handed out.
      !
      ! The nested form -- s3 up to s1, and s4 up to s2 only when s3 equals s1 --
      ! is exactly "every shell pair (s3,s4) at or before (s1,s2)" in the
      ! canonical pair ordering, so enumerating pairs and taking kl <= ij covers
      ! the same quartets once each. Flattening is what makes the outer loop
      ! divisible; the triangular nest is not.
      npair = tab%nbas*(tab%nbas + 1)/2
      allocate (pair_i(npair), pair_j(npair))
      ipair = 0
      do s1 = 1, tab%nbas
         do s2 = 1, s1
            ipair = ipair + 1
            pair_i(ipair) = s1
            pair_j(ipair) = s2
         end do
      end do
      call pair_work_order(pair_i, pair_j, dims, order)

      allocate (g(n, n), d_half(n, n))
      g = 0.0_dp

      ! The six-update form below is written for D = C_occ C_occ^T, the density
      ! without its factor of two. `build_density` produces D = 2 C_occ C_occ^T,
      ! so halve it here. Skipping this makes both J and K exactly twice too
      ! large -- an error that still converges, to a badly wrong energy.
      d_half = 0.5_dp*density

      ! Built from `d_half` rather than `density` because `d_half` is what the six
      ! updates below actually multiply -- a bound has to be on the quantity that
      ! is used, not on a factor of two away from it. Skipped when the caller has
      ! opted out: `dsh` is then never read, so a zero-size placeholder keeps the
      ! `shared` clause legal without paying for the pass over the density.
      if (weight_density) then
         call block_density_max(d_half, tab%nbas, offs, dims, dsh)
      else
         allocate (dsh(0, 0))
      end if

      ! The four exchange updates below carry a quarter each, so folding the
      ! exchange fraction in here scales all four at once at no per-quartet cost.
      kq = 0.25_dp
      if (present(k_scale)) kq = 0.25_dp*k_scale
      jf = 1.0_dp
      if (present(j_scale)) jf = j_scale

      ! A copy, because the range-separation parameter lives in `env` and the
      ! molecule is read-only here. Cheap -- env is a few thousand doubles -- and it
      ! keeps a range-separated build from mutating state a caller shares.
      !
      ! `+ 1` because the `env` pointer constants are libcint's own 0-based
      ! offsets and are *not* converted by the Fortran interface -- unlike the
      ! `atm`/`bas` column constants, which are. Getting this wrong is silent:
      ! slot 8 is `PTR_RINV_ZETA`, which a plain two-electron integral ignores,
      ! so the attenuated build returns full-range exchange, the two passes sum
      ! back to unscaled K, and the SCF converges several Hartree out.
      env_local = tab%env
      if (present(omega)) env_local(LIBCINT_PTR_RANGE_OMEGA + 1) = omega

      ! Created once and reused for every quartet in this build.
      opt = c_null_ptr
      call two_electron_optimizer(mol%cartesian, opt, mol%atm, mol%natm, tab%bas, &
                                  tab%nbas, env_local)

      n_total = 0_int64
      n_computed = 0_int64
      n_screened = 0_int64
      n_schwarz = 0_int64
      n_density = 0_int64
      t_thread_max = 0.0_dp
      t_thread_sum = 0.0_dp

      ! Threaded over bra pairs. libcint carries no mutable state across calls
      ! -- the 2e path has no static globals, and `opt` is written once here and
      ! only read inside -- so the integrals themselves need nothing but a
      ! private buffer each.
      !
      ! The accumulator is the part that does need care. The six updates below
      ! scatter into `g` at positions that depend on all four shells, so two
      ! threads holding different quartets can land on the same element. Each
      ! thread therefore fills its own copy and adds it in once at the end.
      ! Atomics on the innermost statement would be correct too and far slower:
      ! six of them per integral, on the hottest line in the program.
      !
      ! The cost of that is one n*n array per thread. At forty threads it is 18
      ! MB for the 237-function case here, and it grows as n^2 -- worth watching
      ! if this ever meets a few thousand functions, where a blocked scheme
      ! along the lines of GTFock's would be the answer.
      !
      ! `schedule(dynamic)`: pair ij does ij quartets, so the last chunk is
      ! thousands of times the first, and a static split would leave most
      ! threads idle waiting for the tail.
      !$omp parallel default(none) &
      !$omp    shared(mol, tab, bq, dsh, d_half, g, dims, offs, pair_i, pair_j, order, npair, tol, opt, n, &
      !$omp           block_max, kq, jf, env_local, weight_density) &
      !$omp    private(itask, ij, kl, s1, s2, s3, s4, d1, d2, d3, d4, o1, o2, o3, o4, &
      !$omp            shls, f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, deg, value, scaled, schwarz, &
      !$omp            buf, g_local, thread_clock, t_thread) &
      !$omp    reduction(+:n_total, n_computed, n_screened, n_schwarz, n_density) &
      !$omp    reduction(max:t_thread_max) reduction(+:t_thread_sum)
      allocate (buf(block_max**4))
      allocate (g_local(n, n))
      g_local = 0.0_dp

      call thread_clock%start()
      !$omp do schedule(dynamic)
      do itask = 1, npair
         ij = order(itask)
         s1 = pair_i(ij)
         s2 = pair_j(ij)
         d1 = dims(s1)
         o1 = offs(s1)
         d2 = dims(s2)
         o2 = offs(s2)

         do kl = 1, ij
            s3 = pair_i(kl)
            s4 = pair_j(kl)
            d3 = dims(s3)
            o3 = offs(s3)
            d4 = dims(s4)
            o4 = offs(s4)

            n_total = n_total + 1_int64

            ! One criterion, not two. The quantity that must be small is the
            ! *contribution*, so that is what the test compares; screening on the
            ! bare Schwarz bound as well would discard quartets whose weight
            ! exceeds one -- `deg` reaches eight -- and those are exactly the ones
            ! that matter most. The two counters only attribute the decision:
            ! whether the Schwarz bound alone would have been enough to reach it.
            deg = pair_degeneracy(s1, s2, s3, s4)
            schwarz = bq(s1, s2)*bq(s3, s4)
            if (weight_density) then
               if (schwarz*density_weight(dsh, s1, s2, s3, s4, jf, kq, deg) < tol) then
                  n_screened = n_screened + 1_int64
                  if (schwarz < tol) then
                     n_schwarz = n_schwarz + 1_int64
                  else
                     n_density = n_density + 1_int64
                  end if
                  cycle
               end if
            else
               if (schwarz < tol) then
                  n_screened = n_screened + 1_int64
                  n_schwarz = n_schwarz + 1_int64
                  cycle
               end if
            end if

            shls = [s1 - 1, s2 - 1, s3 - 1, s4 - 1]
            ret = two_electron_block(mol%cartesian, buf, shls, mol%atm, mol%natm, &
                                     tab%bas, tab%nbas, env_local, opt)
            if (ret == 0) then
               n_screened = n_screened + 1_int64
               cycle
            end if
            n_computed = n_computed + 1_int64

            do f4 = 1, d4
               b4 = o4 + f4
               do f3 = 1, d3
                  b3 = o3 + f3
                  do f2 = 1, d2
                     b2 = o2 + f2
                     do f1 = 1, d1
                        b1 = o1 + f1

                        idx = f1 + (f2 - 1)*d1 + (f3 - 1)*d1*d2 + (f4 - 1)*d1*d2*d3
                        value = buf(idx)
                        scaled = value*deg

                        ! Two Coulomb and four exchange contributions. g is
                        ! not symmetric as it stands; symmetrising at the
                        ! end is what makes these six updates equivalent to
                        ! the full sum over all eight permutations.
                        g_local(b1, b2) = g_local(b1, b2) + jf*d_half(b3, b4)*scaled
                        g_local(b3, b4) = g_local(b3, b4) + jf*d_half(b1, b2)*scaled
                        g_local(b1, b3) = g_local(b1, b3) - kq*d_half(b2, b4)*scaled
                        g_local(b2, b4) = g_local(b2, b4) - kq*d_half(b1, b3)*scaled
                        g_local(b1, b4) = g_local(b1, b4) - kq*d_half(b2, b3)*scaled
                        g_local(b2, b3) = g_local(b2, b3) - kq*d_half(b1, b4)*scaled
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end do nowait
      ! `nowait`, so what is timed is this thread's share of the loop and not
      ! the wait for the slowest one -- the wait is the thing being measured.
      call thread_clock%stop()
      t_thread = thread_clock%get_elapsed_time()
      t_thread_max = max(t_thread_max, t_thread)
      t_thread_sum = t_thread_sum + t_thread

      !$omp critical(mqc_direct_fock_accumulate)
      g = g + g_local
      !$omp end critical(mqc_direct_fock_accumulate)

      deallocate (buf, g_local)
      !$omp end parallel

      stats%quartets_total = n_total
      stats%quartets_computed = n_computed
      stats%quartets_screened = n_screened
      stats%screened_schwarz = n_schwarz
      stats%screened_density = n_density
      if (t_thread_sum > 0.0_dp) then
         stats%thread_imbalance = t_thread_max/(t_thread_sum/real(omp_threads(), dp))
      end if

      call libcint_del_optimizer(opt)

      fock = h + 0.5_dp*(g + transpose(g))

      deallocate (g, d_half, dims, offs, pair_i, pair_j, env_local)
   end subroutine build_fock_direct

   subroutine build_fock_direct_many(mol, h, densities, bounds, focks, stats, error, &
                                     screen_tol, k_scale)
      !! F = H + J - K/2 for many densities, over one pass of the integrals
      !!
      !! **Why this exists.** In a direct scheme the integral evaluation dominates
      !! and the contractions against it are nearly free, so computing a quartet
      !! once and contracting it against every density in hand is the difference
      !! between one integral pass and N of them. The coupled-perturbed equations
      !! for the dynamic polarizabilities need roughly a hundred right-hand sides
      !! -- nine perturbations times twelve imaginary frequencies -- and the matvec
      !! for each is a Fock build on a different response density, so without this
      !! the frequency loop pays for the integrals a hundred times over. DIIS trial
      !! densities and the distributed multipoles want the same amortization.
      !!
      !! **Screening is shared, deliberately.** The Schwarz bound depends on the
      !! basis and not on any density, so every set sees exactly the same quartets
      !! and no set can be screened differently from its neighbours. That keeps a
      !! batch bit-for-bit equal to the same densities passed one at a time, which
      !! is what makes the single-density wrapper above safe.
      !!
      !! **The win saturates around fourfold, and the accumulator is why.** Measured
      !! on methanol in cc-pVTZ, 116 functions, against the same densities passed one
      !! at a time: 3.9x at 6 densities, 4.3x at 12, and 3.8x at 24 -- so it stops
      !! improving and then reverses. The reason is that this is not
      !! integral-dominated once several sets are in flight. The six updates below
      !! scatter into an `n^2 * n_set` accumulator per thread, which is memory-bound
      !! and grows with the batch, while the integral it reuses does not; at forty
      !! threads and 24 densities that accumulator is already 100 MB and stops
      !! fitting in cache.
      !!
      !! So batch in **chunks of about a dozen** rather than passing a hundred
      !! right-hand sides at once. A hundred solved twelve at a time is nine integral
      !! passes instead of a hundred, which is the fourfold that is actually
      !! available -- not the hundredfold that counting integral passes alone would
      !! suggest.
      !!
      !! **Symmetric densities only.** The six updates below, with their degeneracy
      !! factors, assume `D` is symmetric in two separate places: the factor of two
      !! for `s1 /= s2` stands in for the `mu <-> nu` permutation without adding
      !! `D(nu,mu)`, and likewise for `s3 /= s4`. Hand this an antisymmetric density
      !! and those permutations double where they should cancel, so the Coulomb term
      !! comes back at twice its value instead of at zero. Nothing here detects it.
      !! The `A - B` matvec that frequency-dependent response needs is exactly that
      !! case, and it needs its own accumulation with the eight permutations written
      !! out rather than folded into `deg`.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h(:, :)          !! Core Hamiltonian, added to every set
      real(dp), intent(in) :: densities(:, :, :)  !! (n_ao, n_ao, n_set), each 2 C C^T
      real(dp), intent(in) :: bounds(:, :)     !! From `schwarz_bounds`
      real(dp), allocatable, intent(out) :: focks(:, :, :)
      type(direct_stats_t), intent(out) :: stats
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: screen_tol
      real(dp), intent(in), optional :: k_scale
         !! Fraction of exact exchange, one by default. A pure density
         !! functional's response operator has none; a hybrid's has its
         !! mixing fraction. The Coulomb term is the same either way.

      real(dp), allocatable :: buf(:), g(:, :, :), g_local(:, :, :), d_half(:, :, :)
      real(dp), allocatable :: bq(:, :)
      type(eri_shell_table_t) :: tab
      type(c_ptr) :: opt
      integer :: s1, s2, s3, s4
      integer :: d1, d2, d3, d4, o1, o2, o3, o4
      integer :: shls(4)
      integer :: f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, block_max, n, iset, n_set
      integer :: ij, kl, npair, ipair
      integer, allocatable :: pair_i(:), pair_j(:), dims(:), offs(:), order(:)
      integer :: itask
      integer(int64) :: n_total, n_computed, n_screened
      real(dp) :: tol, deg, value, scaled, kx

      n = mol%nao
      n_set = size(densities, 3)
      if (size(h, 1) /= n .or. size(densities, 1) /= n .or. size(densities, 2) /= n) then
         call error%set(ERROR_VALIDATION, "direct Fock: matrix dimensions do not match the basis")
         return
      end if
      if (n_set < 1) then
         call error%set(ERROR_VALIDATION, "direct Fock: no densities to contract against")
         return
      end if
      allocate (focks(n, n, n_set))

      tol = DEFAULT_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol

      ! The shells the quartet loop runs over: the fused-sp view when the
      ! molecule carries one, its split shells otherwise. libfint's `int2e` is
      ! the one driver that reads a fused L shell, and it is exactly what this
      ! loop calls, so the view is safe here and nowhere looser -- and since
      ! libfint learned to carry a derivative's tensor component through an L
      ! shell, the SCF gradient's int2e_ip1 loop takes the view too. The
      ! Schwarz bounds arrive per split shell either way, straight from
      ! `schwarz_bounds`, and are collapsed to match. Dimensions
      ! and offsets are copied out up front: both are needed inside the
      ! parallel region, and looking them up there would mean every thread
      ! reaching into `bas` for something that does not change.
      call eri_shell_table(mol, tab)
      call eri_schwarz_collapse(mol, bounds, bq)
      dims = tab%dims
      offs = tab%offs(1:tab%nbas)
      block_max = tab%block_max

      ! The quartet loop, flattened onto one index so it can be handed out.
      !
      ! The nested form -- s3 up to s1, and s4 up to s2 only when s3 equals s1 --
      ! is exactly "every shell pair (s3,s4) at or before (s1,s2)" in the
      ! canonical pair ordering, so enumerating pairs and taking kl <= ij covers
      ! the same quartets once each. Flattening is what makes the outer loop
      ! divisible; the triangular nest is not.
      npair = tab%nbas*(tab%nbas + 1)/2
      allocate (pair_i(npair), pair_j(npair))
      ipair = 0
      do s1 = 1, tab%nbas
         do s2 = 1, s1
            ipair = ipair + 1
            pair_i(ipair) = s1
            pair_j(ipair) = s2
         end do
      end do
      call pair_work_order(pair_i, pair_j, dims, order)

      allocate (g(n, n, n_set), d_half(n, n, n_set))
      g = 0.0_dp

      ! The six-update form below is written for D = C_occ C_occ^T, the density
      ! without its factor of two. `build_density` produces D = 2 C_occ C_occ^T,
      ! so halve it here. Skipping this makes both J and K exactly twice too
      ! large -- an error that still converges, to a badly wrong energy.
      kx = 1.0_dp
      if (present(k_scale)) kx = k_scale
      d_half = 0.5_dp*densities

      ! Created once and reused for every quartet in this build.
      opt = c_null_ptr
      call two_electron_optimizer(mol%cartesian, opt, mol%atm, mol%natm, tab%bas, &
                                  tab%nbas, tab%env)

      n_total = 0_int64
      n_computed = 0_int64
      n_screened = 0_int64

      ! Threaded over bra pairs. libcint carries no mutable state across calls
      ! -- the 2e path has no static globals, and `opt` is written once here and
      ! only read inside -- so the integrals themselves need nothing but a
      ! private buffer each.
      !
      ! The accumulator is the part that does need care. The six updates below
      ! scatter into `g` at positions that depend on all four shells, so two
      ! threads holding different quartets can land on the same element. Each
      ! thread therefore fills its own copy and adds it in once at the end.
      ! Atomics on the innermost statement would be correct too and far slower:
      ! six of them per integral, on the hottest line in the program.
      !
      ! The cost of that is one n*n array per thread. At forty threads it is 18
      ! MB for the 237-function case here, and it grows as n^2 -- worth watching
      ! if this ever meets a few thousand functions, where a blocked scheme
      ! along the lines of GTFock's would be the answer.
      !
      ! `schedule(dynamic)`: pair ij does ij quartets, so the last chunk is
      ! thousands of times the first, and a static split would leave most
      ! threads idle waiting for the tail.
      !$omp parallel default(none) &
      !$omp    shared(kx, mol, tab, bq, d_half, g, dims, offs, pair_i, pair_j, order, npair, tol, opt, n, &
      !$omp           block_max, n_set) &
      !$omp    private(itask, ij, kl, s1, s2, s3, s4, d1, d2, d3, d4, o1, o2, o3, o4, &
      !$omp            shls, f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, deg, value, scaled, &
      !$omp            buf, g_local, iset) &
      !$omp    reduction(+:n_total, n_computed, n_screened)
      allocate (buf(block_max**4))
      allocate (g_local(n, n, n_set))
      g_local = 0.0_dp

      !$omp do schedule(dynamic)
      do itask = 1, npair
         ij = order(itask)
         s1 = pair_i(ij)
         s2 = pair_j(ij)
         d1 = dims(s1)
         o1 = offs(s1)
         d2 = dims(s2)
         o2 = offs(s2)

         do kl = 1, ij
            s3 = pair_i(kl)
            s4 = pair_j(kl)
            d3 = dims(s3)
            o3 = offs(s3)
            d4 = dims(s4)
            o4 = offs(s4)

            n_total = n_total + 1_int64

            if (bq(s1, s2)*bq(s3, s4) < tol) then
               n_screened = n_screened + 1_int64
               cycle
            end if

            shls = [s1 - 1, s2 - 1, s3 - 1, s4 - 1]
            ret = two_electron_block(mol%cartesian, buf, shls, mol%atm, mol%natm, &
                                     tab%bas, tab%nbas, tab%env, opt)
            if (ret == 0) then
               n_screened = n_screened + 1_int64
               cycle
            end if
            n_computed = n_computed + 1_int64

            deg = pair_degeneracy(s1, s2, s3, s4)

            do f4 = 1, d4
               b4 = o4 + f4
               do f3 = 1, d3
                  b3 = o3 + f3
                  do f2 = 1, d2
                     b2 = o2 + f2
                     do f1 = 1, d1
                        b1 = o1 + f1

                        idx = f1 + (f2 - 1)*d1 + (f3 - 1)*d1*d2 + (f4 - 1)*d1*d2*d3
                        value = buf(idx)
                        scaled = value*deg

                        ! Two Coulomb and four exchange contributions per set. g
                        ! is not symmetric as it stands; symmetrising at the end is
                        ! what makes these six updates equivalent to the full sum
                        ! over all eight permutations.
                        !
                        ! This inner loop is the entire point of the routine: the
                        ! integral above was computed once and every set reuses it.
                        do iset = 1, n_set
                           g_local(b1, b2, iset) = g_local(b1, b2, iset) &
                                                   + d_half(b3, b4, iset)*scaled
                           g_local(b3, b4, iset) = g_local(b3, b4, iset) &
                                                   + d_half(b1, b2, iset)*scaled
                           g_local(b1, b3, iset) = g_local(b1, b3, iset) &
                                                   - kx*0.25_dp*d_half(b2, b4, iset)*scaled
                           g_local(b2, b4, iset) = g_local(b2, b4, iset) &
                                                   - kx*0.25_dp*d_half(b1, b3, iset)*scaled
                           g_local(b1, b4, iset) = g_local(b1, b4, iset) &
                                                   - kx*0.25_dp*d_half(b2, b3, iset)*scaled
                           g_local(b2, b3, iset) = g_local(b2, b3, iset) &
                                                   - kx*0.25_dp*d_half(b1, b4, iset)*scaled
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end do

      !$omp critical(mqc_direct_fock_accumulate)
      g = g + g_local
      !$omp end critical(mqc_direct_fock_accumulate)

      deallocate (buf, g_local)
      !$omp end parallel

      stats%quartets_total = n_total
      stats%quartets_computed = n_computed
      stats%quartets_screened = n_screened

      call libcint_del_optimizer(opt)

      do iset = 1, n_set
         focks(:, :, iset) = h + 0.5_dp*(g(:, :, iset) + transpose(g(:, :, iset)))
      end do

      deallocate (g, d_half, dims, offs, pair_i, pair_j)
   end subroutine build_fock_direct_many

   subroutine build_fock_direct_nosym(mol, h, densities, bounds, focks, stats, error, &
                                      screen_tol)
      !! J - K/2 for densities of **any** symmetry, over one pass of the integrals
      !!
      !! `build_fock_direct_many` is faster and cannot be used here. It folds three
      !! of the eightfold permutations into a `deg` factor -- doubling for
      !! `s1 /= s2` rather than ever touching `D(nu,mu)` -- which is only the same
      !! thing when `D` is symmetric. Hand it an antisymmetric density and those
      !! permutations reinforce where they must cancel, so the Coulomb term comes
      !! back at twice its size instead of at zero. Nothing detects that, and it is
      !! not repairable by changing the final symmetrisation, because the error is
      !! already in the accumulation.
      !!
      !! So this routine writes the permutations out. For each computed integral it
      !! generates the distinct index tuples that the block enumeration does not
      !! already cover -- exactly the ones `deg` was standing in for, under the same
      !! three conditions -- and applies to each
      !!
      !!     g(p,q) += V D(r,s)          the Coulomb contribution
      !!     g(p,r) -= V D(q,s) / 2      the exchange contribution
      !!
      !! with the real density elements rather than a multiplicity. **There is no
      !! final symmetrisation**, and there must not be: for an antisymmetric density
      !! `J` vanishes and `G = -K/2` is itself antisymmetric, so
      !! `0.5*(g + transpose(g))` would annihilate the entire result.
      !!
      !! For a symmetric density this reduces term by term to what
      !! `build_fock_direct_many` computes, which is asserted in the tests rather
      !! than asserted here.
      !!
      !! **Why this is the wanted operator.** The frequency-dependent
      !! coupled-perturbed equations need `A - B` as well as `A + B`, and both are
      !! this same contraction: symmetrise the response density and get `A + B`,
      !! antisymmetrise it and get `A - B`. The Coulomb term dropping out of the
      !! second is not a special case to code around -- it happens by itself,
      !! because the two-electron integral is symmetric in its ket pair while the
      !! density is not.
      !!
      !! **Cost.** Up to eight tuples times two updates against six updates, so
      !! roughly 2.7x the contraction work per integral. Since the batched build is
      !! contraction-bound rather than integral-bound past a few densities, prefer
      !! two passes -- symmetric densities through `build_fock_direct_many` and
      !! antisymmetric ones through here -- over routing everything through this one.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h(:, :)             !! Added to every set
      real(dp), intent(in) :: densities(:, :, :)  !! (n_ao, n_ao, n_set), any symmetry
      real(dp), intent(in) :: bounds(:, :)        !! From `schwarz_bounds`
      real(dp), allocatable, intent(out) :: focks(:, :, :)
      type(direct_stats_t), intent(out) :: stats
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: screen_tol

      real(dp), allocatable :: buf(:), g(:, :, :), g_local(:, :, :)
      real(dp), allocatable :: bq(:, :)
      type(eri_shell_table_t) :: tab
      type(c_ptr) :: opt
      integer :: s1, s2, s3, s4
      integer :: d1, d2, d3, d4, o1, o2, o3, o4
      integer :: shls(4)
      integer :: f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, block_max, n, iset, n_set
      integer :: ij, kl, npair, ipair, i1, i2, i3, pp, qq, rr, ss, tmp
      integer, allocatable :: pair_i(:), pair_j(:), dims(:), offs(:), order(:)
      integer :: itask
      integer(int64) :: n_total, n_computed, n_screened
      real(dp) :: tol, value
      logical :: swap_bra, swap_ket, swap_pairs

      n = mol%nao
      n_set = size(densities, 3)
      if (size(h, 1) /= n .or. size(densities, 1) /= n .or. size(densities, 2) /= n) then
         call error%set(ERROR_VALIDATION, "direct Fock: matrix dimensions do not match the basis")
         return
      end if
      if (n_set < 1) then
         call error%set(ERROR_VALIDATION, "direct Fock: no densities to contract against")
         return
      end if
      allocate (focks(n, n, n_set))

      tol = DEFAULT_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol

      ! The shells the quartet loop runs over: the fused-sp view when the
      ! molecule carries one, its split shells otherwise. libfint's `int2e` is
      ! the one driver that reads a fused L shell, and it is exactly what this
      ! loop calls, so the view is safe here and nowhere looser -- and since
      ! libfint learned to carry a derivative's tensor component through an L
      ! shell, the SCF gradient's int2e_ip1 loop takes the view too. The
      ! Schwarz bounds arrive per split shell either way, straight from
      ! `schwarz_bounds`, and are collapsed to match. Dimensions
      ! and offsets are copied out up front: both are needed inside the
      ! parallel region, and looking them up there would mean every thread
      ! reaching into `bas` for something that does not change.
      call eri_shell_table(mol, tab)
      call eri_schwarz_collapse(mol, bounds, bq)
      dims = tab%dims
      offs = tab%offs(1:tab%nbas)
      block_max = tab%block_max

      npair = tab%nbas*(tab%nbas + 1)/2
      allocate (pair_i(npair), pair_j(npair))
      ipair = 0
      do s1 = 1, tab%nbas
         do s2 = 1, s1
            ipair = ipair + 1
            pair_i(ipair) = s1
            pair_j(ipair) = s2
         end do
      end do
      call pair_work_order(pair_i, pair_j, dims, order)

      allocate (g(n, n, n_set))
      g = 0.0_dp

      opt = c_null_ptr
      call two_electron_optimizer(mol%cartesian, opt, mol%atm, mol%natm, tab%bas, &
                                  tab%nbas, tab%env)

      n_total = 0_int64
      n_computed = 0_int64
      n_screened = 0_int64

      !$omp parallel default(none) &
      !$omp    shared(mol, tab, bq, densities, g, dims, offs, pair_i, pair_j, order, npair, tol, &
      !$omp           opt, n, block_max, n_set) &
      !$omp    private(itask, ij, kl, s1, s2, s3, s4, d1, d2, d3, d4, o1, o2, o3, o4, &
      !$omp            shls, f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, value, &
      !$omp            buf, g_local, iset, i1, i2, i3, pp, qq, rr, ss, tmp, &
      !$omp            swap_bra, swap_ket, swap_pairs) &
      !$omp    reduction(+:n_total, n_computed, n_screened)
      allocate (buf(block_max**4))
      allocate (g_local(n, n, n_set))
      g_local = 0.0_dp

      !$omp do schedule(dynamic)
      do itask = 1, npair
         ij = order(itask)
         s1 = pair_i(ij)
         s2 = pair_j(ij)
         d1 = dims(s1)
         o1 = offs(s1)
         d2 = dims(s2)
         o2 = offs(s2)

         do kl = 1, ij
            s3 = pair_i(kl)
            s4 = pair_j(kl)
            d3 = dims(s3)
            o3 = offs(s3)
            d4 = dims(s4)
            o4 = offs(s4)

            n_total = n_total + 1_int64

            if (bq(s1, s2)*bq(s3, s4) < tol) then
               n_screened = n_screened + 1_int64
               cycle
            end if

            shls = [s1 - 1, s2 - 1, s3 - 1, s4 - 1]
            ret = two_electron_block(mol%cartesian, buf, shls, mol%atm, mol%natm, &
                                     tab%bas, tab%nbas, tab%env, opt)
            if (ret == 0) then
               n_screened = n_screened + 1_int64
               cycle
            end if
            n_computed = n_computed + 1_int64

            ! Which permutations the block enumeration leaves uncovered. Shell
            ! equality, not function equality: a block with s1 == s2 already runs
            ! over both orderings of its function pair as separate elements, so
            ! generating the swap as well would count it twice. These are the same
            ! three conditions that produce `deg` in `build_fock_direct_many`.
            swap_bra = s1 /= s2
            swap_ket = s3 /= s4
            swap_pairs = .not. (s1 == s3 .and. s2 == s4)

            do f4 = 1, d4
               b4 = o4 + f4
               do f3 = 1, d3
                  b3 = o3 + f3
                  do f2 = 1, d2
                     b2 = o2 + f2
                     do f1 = 1, d1
                        b1 = o1 + f1

                        idx = f1 + (f2 - 1)*d1 + (f3 - 1)*d1*d2 + (f4 - 1)*d1*d2*d3
                        value = buf(idx)

                        ! The orbit, generated by three independent swaps. Skipping
                        ! a swap when its shells coincide is what keeps each
                        ! distinct tuple appearing exactly once.
                        do i1 = 1, 2
                           if (i1 == 2 .and. .not. swap_bra) cycle
                           do i2 = 1, 2
                              if (i2 == 2 .and. .not. swap_ket) cycle
                              do i3 = 1, 2
                                 if (i3 == 2 .and. .not. swap_pairs) cycle

                                 pp = b1
                                 qq = b2
                                 rr = b3
                                 ss = b4
                                 if (i1 == 2) then
                                    tmp = pp
                                    pp = qq
                                    qq = tmp
                                 end if
                                 if (i2 == 2) then
                                    tmp = rr
                                    rr = ss
                                    ss = tmp
                                 end if
                                 if (i3 == 2) then
                                    tmp = pp
                                    pp = rr
                                    rr = tmp
                                    tmp = qq
                                    qq = ss
                                    ss = tmp
                                 end if

                                 do iset = 1, n_set
                                    g_local(pp, qq, iset) = g_local(pp, qq, iset) &
                                                            + densities(rr, ss, iset)*value
                                    g_local(pp, rr, iset) = g_local(pp, rr, iset) &
                                                            - 0.5_dp*densities(qq, ss, iset)*value
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end do

      !$omp critical(mqc_direct_nosym_accumulate)
      g = g + g_local
      !$omp end critical(mqc_direct_nosym_accumulate)

      deallocate (buf, g_local)
      !$omp end parallel

      stats%quartets_total = n_total
      stats%quartets_computed = n_computed
      stats%quartets_screened = n_screened

      call libcint_del_optimizer(opt)

      ! No symmetrisation. See the note at the top: for an antisymmetric density
      ! the result is antisymmetric, and symmetrising would return zero.
      do iset = 1, n_set
         focks(:, :, iset) = h + g(:, :, iset)
      end do

      deallocate (g, dims, offs, pair_i, pair_j)
   end subroutine build_fock_direct_nosym

   subroutine build_fock_direct_uhf(mol, h, d_alpha, d_beta, bounds, fock_a, fock_b, stats, error, &
                                    screen_tol, k_scale, j_scale, omega)
      !! F_sigma = H + J(D_alpha + D_beta) - K(D_sigma), without the tensor
      !!
      !! The same quartet loop as the closed-shell build, differing only in what
      !! the six updates read. Coulomb draws on half the total density and
      !! exchange on the same-spin density -- and in the closed-shell case those
      !! are the same matrix, D/2, which is why one form covers both.
      !!
      !! Kept as its own routine rather than folded into `build_fock_direct`
      !! with a flag: the closed-shell path would then carry two accumulators
      !! and two sets of updates to do one spin's work, and that path is the hot
      !! one. The cost is that the two loops have to stay in step; the
      !! elementwise check against the in-core UHF build is what enforces it.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h(:, :)          !! Core Hamiltonian
      real(dp), intent(in) :: d_alpha(:, :)    !! D_alpha = C_alpha C_alpha^T
      real(dp), intent(in) :: d_beta(:, :)     !! D_beta  = C_beta  C_beta^T
      real(dp), intent(in) :: bounds(:, :)     !! From `schwarz_bounds`
      real(dp), intent(out) :: fock_a(:, :), fock_b(:, :)
      type(direct_stats_t), intent(out) :: stats
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: screen_tol
      real(dp), intent(in), optional :: k_scale
         !! Fraction of exact exchange to keep, as in the closed-shell build. One
         !! is Hartree-Fock and the default; a hybrid Kohn-Sham build wants less.
      real(dp), intent(in), optional :: j_scale
         !! Fraction of the Coulomb term, default one. Zero is what the second
         !! pass of a range-separated build wants.
      real(dp), intent(in), optional :: omega
         !! Range-separation parameter, through `env(PTR_RANGE_OMEGA)`. Zero, the
         !! default, is the full Coulomb kernel. See the closed-shell build for
         !! why the Schwarz bounds need no rebuilding for an attenuated pass.

      real(dp), allocatable :: buf(:), ga(:, :), gb(:, :), ga_local(:, :), gb_local(:, :)
      real(dp), allocatable :: d_coul(:, :)
      real(dp), allocatable :: dsh(:, :)
      real(dp), allocatable :: env_local(:)
      real(dp), allocatable :: bq(:, :)
      type(eri_shell_table_t) :: tab
      real(dp) :: kq, jf
      type(c_ptr) :: opt
      integer :: s1, s2, s3, s4
      integer :: d1, d2, d3, d4, o1, o2, o3, o4
      integer :: shls(4)
      integer :: f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, block_max, n
      integer :: ij, kl, npair, ipair
      integer, allocatable :: pair_i(:), pair_j(:), dims(:), offs(:), order(:)
      integer :: itask
      integer(int64) :: n_total, n_computed, n_screened, n_schwarz, n_density
      real(dp) :: schwarz
      real(dp) :: tol, deg, value, scaled

      n = mol%nao
      if (size(h, 1) /= n .or. size(d_alpha, 1) /= n .or. size(d_beta, 1) /= n &
          .or. size(fock_a, 1) /= n .or. size(fock_b, 1) /= n) then
         call error%set(ERROR_VALIDATION, "direct UHF Fock: matrix dimensions do not match the basis")
         return
      end if

      tol = DEFAULT_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol

      ! The shells the quartet loop runs over: the fused-sp view when the
      ! molecule carries one, its split shells otherwise. libfint's `int2e` is
      ! the one driver that reads a fused L shell, and it is exactly what this
      ! loop calls, so the view is safe here and nowhere looser -- and since
      ! libfint learned to carry a derivative's tensor component through an L
      ! shell, the SCF gradient's int2e_ip1 loop takes the view too. The
      ! Schwarz bounds arrive per split shell either way, straight from
      ! `schwarz_bounds`, and are collapsed to match. Dimensions
      ! and offsets are copied out up front: both are needed inside the
      ! parallel region, and looking them up there would mean every thread
      ! reaching into `bas` for something that does not change.
      call eri_shell_table(mol, tab)
      call eri_schwarz_collapse(mol, bounds, bq)
      dims = tab%dims
      offs = tab%offs(1:tab%nbas)
      block_max = tab%block_max

      ! The quartet loop, flattened onto one index so it can be handed out.
      !
      ! The nested form -- s3 up to s1, and s4 up to s2 only when s3 equals s1 --
      ! is exactly "every shell pair (s3,s4) at or before (s1,s2)" in the
      ! canonical pair ordering, so enumerating pairs and taking kl <= ij covers
      ! the same quartets once each. Flattening is what makes the outer loop
      ! divisible; the triangular nest is not.
      npair = tab%nbas*(tab%nbas + 1)/2
      allocate (pair_i(npair), pair_j(npair))
      ipair = 0
      do s1 = 1, tab%nbas
         do s2 = 1, s1
            ipair = ipair + 1
            pair_i(ipair) = s1
            pair_j(ipair) = s2
         end do
      end do
      call pair_work_order(pair_i, pair_j, dims, order)

      allocate (ga(n, n), gb(n, n), d_coul(n, n))
      ga = 0.0_dp
      gb = 0.0_dp

      ! The eight exchange updates below carry a quarter each, so the exchange
      ! fraction folds in here and scales all eight at no per-quartet cost.
      kq = 0.25_dp
      if (present(k_scale)) kq = 0.25_dp*k_scale
      jf = 1.0_dp
      if (present(j_scale)) jf = j_scale

      ! Coulomb sees half the total density, exchange sees the same-spin one.
      !
      ! The closed-shell build passes D/2 to both, and that is not a coincidence
      ! worth losing: there D_alpha = D_beta = D/2, so half the total density is
      ! D/2 and the same-spin density is D/2. The two forms are one form.
      !
      ! `j_scale` rides on this matrix rather than on the two Coulomb updates,
      ! which is the same arithmetic one level out.
      d_coul = jf*0.5_dp*(d_alpha + d_beta)

      ! Elementwise largest of the three matrices the updates read. Taking the
      ! maximum rather than one of them keeps the bound a bound for every term:
      ! the Coulomb updates read `d_coul` and the exchange updates read the two
      ! spin densities, and any of the three can be the big one.
      call block_density_max(max(abs(d_coul), abs(d_alpha), abs(d_beta)), tab%nbas, offs, dims, dsh)

      ! A copy, because the range-separation parameter lives in `env` and the
      ! molecule is read-only here. See the closed-shell build for why the slot
      ! index carries a `+ 1`.
      env_local = tab%env
      if (present(omega)) env_local(LIBCINT_PTR_RANGE_OMEGA + 1) = omega

      ! Created once and reused for every quartet in this build.
      opt = c_null_ptr
      call two_electron_optimizer(mol%cartesian, opt, mol%atm, mol%natm, tab%bas, &
                                  tab%nbas, env_local)

      n_total = 0_int64
      n_computed = 0_int64
      n_screened = 0_int64
      n_schwarz = 0_int64
      n_density = 0_int64

      ! Threaded over bra pairs. libcint carries no mutable state across calls
      ! -- the 2e path has no static globals, and `opt` is written once here and
      ! only read inside -- so the integrals themselves need nothing but a
      ! private buffer each.
      !
      ! The accumulator is the part that does need care. The six updates below
      ! scatter into `g` at positions that depend on all four shells, so two
      ! threads holding different quartets can land on the same element. Each
      ! thread therefore fills its own copy and adds it in once at the end.
      ! Atomics on the innermost statement would be correct too and far slower:
      ! six of them per integral, on the hottest line in the program.
      !
      ! The cost of that is one n*n array per thread. At forty threads it is 18
      ! MB for the 237-function case here, and it grows as n^2 -- worth watching
      ! if this ever meets a few thousand functions, where a blocked scheme
      ! along the lines of GTFock's would be the answer.
      !
      ! `schedule(dynamic)`: pair ij does ij quartets, so the last chunk is
      ! thousands of times the first, and a static split would leave most
      ! threads idle waiting for the tail.
      !$omp parallel default(none) &
      !$omp    shared(mol, tab, bq, dsh, d_coul, d_alpha, d_beta, ga, gb, dims, offs, pair_i, pair_j, &
      !$omp            order, &
      !$omp            npair, tol, opt, n, block_max, kq, env_local) &
      !$omp    private(itask, ij, kl, s1, s2, s3, s4, d1, d2, d3, d4, o1, o2, o3, o4, &
      !$omp            shls, f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, deg, value, scaled, schwarz, &
      !$omp            buf, ga_local, gb_local) &
      !$omp    reduction(+:n_total, n_computed, n_screened, n_schwarz, n_density)
      allocate (buf(block_max**4))
      allocate (ga_local(n, n), gb_local(n, n))
      ga_local = 0.0_dp
      gb_local = 0.0_dp

      !$omp do schedule(dynamic)
      do itask = 1, npair
         ij = order(itask)
         s1 = pair_i(ij)
         s2 = pair_j(ij)
         d1 = dims(s1)
         o1 = offs(s1)
         d2 = dims(s2)
         o2 = offs(s2)

         do kl = 1, ij
            s3 = pair_i(kl)
            s4 = pair_j(kl)
            d3 = dims(s3)
            o3 = offs(s3)
            d4 = dims(s4)
            o4 = offs(s4)

            n_total = n_total + 1_int64

            ! One rather than `jf` for the Coulomb weight: this path folded `jf`
            ! into `d_coul` before the loop, so `dsh` already carries it and
            ! passing it again would square it.
            ! One criterion, not two. The quantity that must be small is the
            ! *contribution*, so that is what the test compares; screening on the
            ! bare Schwarz bound as well would discard quartets whose weight
            ! exceeds one -- `deg` reaches eight -- and those are exactly the ones
            ! that matter most. The two counters only attribute the decision:
            ! whether the Schwarz bound alone would have been enough to reach it.
            deg = pair_degeneracy(s1, s2, s3, s4)
            schwarz = bq(s1, s2)*bq(s3, s4)
            if (schwarz*density_weight(dsh, s1, s2, s3, s4, 1.0_dp, kq, deg) < tol) then
               n_screened = n_screened + 1_int64
               if (schwarz < tol) then
                  n_schwarz = n_schwarz + 1_int64
               else
                  n_density = n_density + 1_int64
               end if
               cycle
            end if

            shls = [s1 - 1, s2 - 1, s3 - 1, s4 - 1]
            ret = two_electron_block(mol%cartesian, buf, shls, mol%atm, mol%natm, &
                                     tab%bas, tab%nbas, env_local, opt)
            if (ret == 0) then
               n_screened = n_screened + 1_int64
               cycle
            end if
            n_computed = n_computed + 1_int64

            do f4 = 1, d4
               b4 = o4 + f4
               do f3 = 1, d3
                  b3 = o3 + f3
                  do f2 = 1, d2
                     b2 = o2 + f2
                     do f1 = 1, d1
                        b1 = o1 + f1

                        idx = f1 + (f2 - 1)*d1 + (f3 - 1)*d1*d2 + (f4 - 1)*d1*d2*d3
                        value = buf(idx)
                        scaled = value*deg

                        ! Two Coulomb and four exchange contributions. g is
                        ! not symmetric as it stands; symmetrising at the
                        ! end is what makes these six updates equivalent to
                        ! the full sum over all eight permutations.
                        ga_local(b1, b2) = ga_local(b1, b2) + d_coul(b3, b4)*scaled
                        ga_local(b3, b4) = ga_local(b3, b4) + d_coul(b1, b2)*scaled
                        ga_local(b1, b3) = ga_local(b1, b3) - kq*d_alpha(b2, b4)*scaled
                        ga_local(b2, b4) = ga_local(b2, b4) - kq*d_alpha(b1, b3)*scaled
                        ga_local(b1, b4) = ga_local(b1, b4) - kq*d_alpha(b2, b3)*scaled
                        ga_local(b2, b3) = ga_local(b2, b3) - kq*d_alpha(b1, b4)*scaled

                        gb_local(b1, b2) = gb_local(b1, b2) + d_coul(b3, b4)*scaled
                        gb_local(b3, b4) = gb_local(b3, b4) + d_coul(b1, b2)*scaled
                        gb_local(b1, b3) = gb_local(b1, b3) - kq*d_beta(b2, b4)*scaled
                        gb_local(b2, b4) = gb_local(b2, b4) - kq*d_beta(b1, b3)*scaled
                        gb_local(b1, b4) = gb_local(b1, b4) - kq*d_beta(b2, b3)*scaled
                        gb_local(b2, b3) = gb_local(b2, b3) - kq*d_beta(b1, b4)*scaled
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end do

      !$omp critical(mqc_direct_fock_uhf_accumulate)
      ga = ga + ga_local
      gb = gb + gb_local
      !$omp end critical(mqc_direct_fock_uhf_accumulate)

      deallocate (buf, ga_local, gb_local)
      !$omp end parallel

      stats%quartets_total = n_total
      stats%quartets_computed = n_computed
      stats%quartets_screened = n_screened
      stats%screened_schwarz = n_schwarz
      stats%screened_density = n_density

      call libcint_del_optimizer(opt)

      fock_a = h + 0.5_dp*(ga + transpose(ga))
      fock_b = h + 0.5_dp*(gb + transpose(gb))

      deallocate (ga, gb, d_coul, dims, offs, pair_i, pair_j, env_local)
   end subroutine build_fock_direct_uhf

end module mqc_libcint_direct
