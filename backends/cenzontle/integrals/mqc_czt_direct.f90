module mqc_czt_direct
   !! Direct Fock construction: integrals consumed as produced, never stored
   !!
   !! The in-core path in `mqc_czt_integrals` holds all n^4 integrals -- 65 GB
   !! at three hundred basis functions -- and is a reference implementation rather
   !! than a method. This builds the Fock matrix without ever forming the tensor.
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
   !! without being computed, which turns the formal n^4 into roughly n^2 for
   !! extended systems.
   !!
   !! **libcint's optimizer.** CINTOpt precomputes per-shell-pair data reused
   !! across quartet calls. Created once per Fock build.
   !!
   !! **The degeneracy factors test shell equality, not function equality.** A
   !! quartet where M == N already has both (mu nu|..) and (nu mu|..) inside its
   !! own block, so that permutation must not be counted again. Getting it wrong
   !! produces a Fock matrix wrong by a factor of two on its diagonal blocks
   !! only; `check_direct` compares against the in-core build elementwise.
   use pic_types, only: dp, int64, int_index
   use pic_sorting, only: sort_index
   use pic_timer, only: timer_type
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, shell_dim, &
                                two_electron_block, two_electron_optimizer, &
                                eri_shell_table_t, eri_shell_table, &
                                eri_schwarz_collapse
   use libcint_fortran, only: libcint_del_optimizer, LIBCINT_PTR_RANGE_OMEGA
   implicit none
   private

   public :: schwarz_bounds
   public :: shell_density_max
   public :: block_density_max
   public :: pair_degeneracy
   public :: pair_work_order
   public :: omp_threads
   public :: build_fock_direct
   public :: build_fock_direct_many
   public :: build_fock_direct_nosym
   public :: build_fock_direct_uhf
   public :: direct_stats_t
   public :: DEFAULT_SCREEN_TOL

   real(dp), parameter :: DEFAULT_SCREEN_TOL = 1.0e-11_dp
   integer, parameter :: FOCK_TILE_FUNCS = 16
      !! Functions per shell tile in the batched build's quartet loop. About
      !! one heavy atom in a double-zeta basis; the twelve blocks a tile
      !! quartet touches are `12 * FOCK_TILE_FUNCS**2 * n_set` doubles.
      !! Quartets whose Schwarz bound falls below this are skipped.
      !!
      !! The GTFock paper's value. The count of surviving quartets is
      !! insensitive to the threshold over a decade or so.

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
         !! Slowest thread's time on the quartet loop divided by the mean. One is
         !! perfect balance. A ratio within a single run, so contention on a
         !! shared node largely divides out.
      integer(int64) :: screened_density = 0
         !! Skipped because the *contribution* is negligible although the
         !! integral is not -- the extra reach the density weighting buys, and
         !! exactly reproducible where a wall-time measurement of it is not.
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
      !! Costs n_shells^2 quartet evaluations, done once per geometry: the bounds
      !! depend on the basis and the positions, not the density, so a whole SCF
      !! reuses one set.
      type(czt_molecule_t), intent(in) :: mol
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
      !! Threads available to the NEXT parallel region, or one without OpenMP
      !!
      !! `omp_get_max_threads()`, so this is an upper bound and not a count of
      !! what ran: dynamic adjustment or a `num_threads` clause can give fewer.
      !! Right for sizing per-thread storage before a region, wrong for
      !! normalising anything measured inside one -- count those in the region.
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
      !! thread idle behind it. Sorting on the estimated cost fixes both the
      !! trend and the variation within it.
      !!
      !! The cost is the unscreened function-quartet count, so screening makes
      !! this an estimate; the ordering only has to be roughly right.
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
      !! A quartet multiplies six density elements, and if all six are negligible
      !! it contributes nothing however large the integral is. The Schwarz bound
      !! alone cannot see that.
      !!
      !! Blocked to shells because the screening decision is made per shell
      !! quartet, and it has to be the largest in the block or the bound would
      !! not be one. Rebuilt every Fock build, unlike the Schwarz bounds: this
      !! one is a property of the density.
      type(czt_molecule_t), intent(in) :: mol
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
      !! split shell for any caller blocked that way.
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
      !! The six elements are the six the Fock updates read, with the weights
      !! they are multiplied by, so `bound * denmax` bounds the largest change any
      !! one Fock element can see from this quartet.
      !!
      !! `deg` belongs in here: the permutational factor multiplies the
      !! contribution, so leaving it out would make the bound smaller than the
      !! thing it bounds.
      real(dp), intent(in) :: dsh(:, :)
      integer, intent(in) :: s1, s2, s3, s4
      real(dp), intent(in) :: jf, kq, deg
      real(dp) :: denmax

      denmax = deg*max(jf*dsh(s1, s2), jf*dsh(s3, s4), &
                       kq*dsh(s1, s3), kq*dsh(s1, s4), &
                       kq*dsh(s2, s3), kq*dsh(s2, s4))
   end function density_weight

   pure subroutine shell_tiles(dims, budget, first, last, ntile)
      !! Consecutive shells grouped into tiles of at most `budget` functions
      !!
      !! A shell wider than the budget is a tile of its own. `first(t)` and
      !! `last(t)` are shell indices, inclusive.
      integer, intent(in) :: dims(:)
      integer, intent(in) :: budget
      integer, allocatable, intent(out) :: first(:), last(:)
      integer, intent(out) :: ntile

      integer :: s, nbas, used
      integer, allocatable :: f(:), l(:)

      nbas = size(dims)
      allocate (f(nbas), l(nbas))
      ntile = 0
      used = 0
      do s = 1, nbas
         if (ntile == 0 .or. used + dims(s) > budget) then
            ntile = ntile + 1
            f(ntile) = s
            used = 0
         end if
         l(ntile) = s
         used = used + dims(s)
      end do
      first = f(1:ntile)
      last = l(1:ntile)
   end subroutine shell_tiles

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
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h(:, :)          !! Core Hamiltonian
      real(dp), intent(in) :: density(:, :)    !! D = 2 C_occ C_occ^T
      real(dp), intent(in) :: bounds(:, :)     !! From `schwarz_bounds`
      real(dp), intent(out) :: fock(:, :)
      type(direct_stats_t), intent(out) :: stats
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: screen_tol
      real(dp), intent(in), optional :: k_scale
         !! Fraction of exact exchange to keep. One is Hartree-Fock and the
         !! default; zero is pure density-functional exchange; a hybrid is
         !! between.
      real(dp), intent(in), optional :: j_scale
         !! Fraction of the Coulomb term, default one. Zero is what a
         !! range-separated hybrid's *second* pass wants: a long-range exchange
         !! matrix and nothing else.
      real(dp), intent(in), optional :: omega
         !! Range-separation parameter. Zero, the default, is the full Coulomb
         !! kernel. Positive gives the long-range erf-attenuated one and negative
         !! the short-range complement, through `env(PTR_RANGE_OMEGA)`.
         !!
         !! The Schwarz bounds passed in are the full-kernel ones and stay valid
         !! here without being rebuilt: erf(omega r)/r <= 1/r pointwise, so an
         !! attenuated quartet is never larger than the bound screening it, and
         !! the screening can only become conservative rather than wrong.
      logical, intent(in), optional :: density_screen
         !! Weight the Schwarz bound by the density it multiplies before
         !! screening (default true). An SCF wants it.
         !!
         !! **A CPHF response build must pass false.** Its `density` is a Krylov
         !! trial vector the solver drives towards zero, so a screen keyed on the
         !! density's magnitude tightens as the solve proceeds: the operator
         !! stops being the same linear map from one matvec to the next and the
         !! iteration cannot converge. False recovers the plain Schwarz screen,
         !! which depends on the basis alone, and is bit-for-bit equal to
         !! `build_fock_direct_many` with one density.

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
      integer :: n_thread_ran   !! Threads that actually entered the region
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
      ! one of the two drivers that read a fused L shell, and is exactly what
      ! this loop calls. The Schwarz bounds arrive per split shell either way and
      ! are collapsed to match. Dimensions and offsets are copied out up front:
      ! both are read inside the parallel region.
      call eri_shell_table(mol, tab)
      call eri_schwarz_collapse(mol, bounds, bq)
      dims = tab%dims
      offs = tab%offs(1:tab%nbas)
      block_max = tab%block_max

      ! The quartet loop, flattened onto one index so it can be handed out.
      ! Enumerating shell pairs and taking `kl <= ij` covers the same quartets
      ! once each as the triangular nest, and is divisible where that is not.
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
      ! updates below multiply. Skipped when the caller has opted out: `dsh` is
      ! then never read, and a zero-size placeholder keeps the `shared` clause
      ! legal.
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
      ! molecule is read-only here.
      !
      ! **`+ 1` because the `env` pointer constants are libcint's own 0-based
      ! offsets and are not converted by the Fortran interface**, unlike the
      ! `atm`/`bas` column constants, which are. Getting this wrong is silent:
      ! slot 8 is `PTR_RINV_ZETA`, which a plain two-electron integral ignores,
      ! so the attenuated build returns full-range exchange and the SCF converges
      ! several Hartree out.
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
      n_thread_ran = 0

      ! Threaded over bra pairs. libcint carries no mutable state across calls --
      ! the 2e path has no static globals, and `opt` is written once here and
      ! only read inside -- so the integrals need nothing but a private buffer
      ! each.
      !
      ! The accumulator does need care. The updates below scatter at positions
      ! that depend on all four shells, so two threads holding
      ! different quartets can land on the same element. Each thread fills its
      ! own copy and adds it in once at the end; atomics on the innermost
      ! statement would be correct too and far slower. The cost is one n*n array
      ! per thread, growing as n^2 -- **worth watching at a few thousand
      ! functions**, where a blocked scheme like GTFock's would be the answer.
      !
      ! `schedule(dynamic)`: pair ij does ij quartets, so the last chunk is
      ! thousands of times the first.
      !$omp parallel default(none) &
      !$omp    shared(mol, tab, bq, dsh, d_half, g, dims, offs, pair_i, pair_j, order, npair, tol, opt, n, &
      !$omp           block_max, kq, jf, env_local, weight_density) &
      !$omp    private(itask, ij, kl, s1, s2, s3, s4, d1, d2, d3, d4, o1, o2, o3, o4, &
      !$omp            shls, f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, deg, value, scaled, schwarz, &
      !$omp            buf, g_local, thread_clock, t_thread) &
      !$omp    reduction(+:n_total, n_computed, n_screened, n_schwarz, n_density) &
      !$omp    reduction(max:t_thread_max) reduction(+:t_thread_sum, n_thread_ran)
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
            ! exceeds one -- `deg` reaches eight. The two counters only attribute
            ! the decision.
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
      ! `nowait`, so what is timed is this thread's share of the loop and not the
      ! wait for the slowest one, which is the thing being measured.
      call thread_clock%stop()
      t_thread = thread_clock%get_elapsed_time()
      t_thread_max = max(t_thread_max, t_thread)
      t_thread_sum = t_thread_sum + t_thread
      n_thread_ran = n_thread_ran + 1

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
      ! Counted inside the region rather than read from `omp_threads()`, which
      ! reports the maximum available for the NEXT region. `t_thread_sum` covers
      ! the threads that actually entered this one, and dynamic adjustment can
      ! make that fewer -- dividing by the larger number would understate the
      ! mean and overstate the imbalance.
      if (t_thread_sum > 0.0_dp .and. n_thread_ran > 0) then
         stats%thread_imbalance = t_thread_max/(t_thread_sum/real(n_thread_ran, dp))
      end if

      call libcint_del_optimizer(opt)

      fock = h + 0.5_dp*(g + transpose(g))

      deallocate (g, d_half, dims, offs, pair_i, pair_j, env_local)
   end subroutine build_fock_direct

   subroutine build_fock_direct_many(mol, h, densities, bounds, focks, stats, error, &
                                     screen_tol, k_scale, j_scale, omega, density_screen)
      !! F = H + J - K/2 for many densities, over one pass of the integrals
      !!
      !! In a direct scheme the integral evaluation dominates and the contractions
      !! against it are nearly free, so one quartet contracted against every
      !! density in hand is the difference between one integral pass and N of
      !! them.
      !!
      !! **Screening is shared, deliberately.** The Schwarz bound depends on the
      !! basis and not on any density, so every set sees exactly the same
      !! quartets. That keeps a batch bit-for-bit equal to the same densities
      !! passed one at a time, which is what makes the single-density wrapper
      !! above safe.
      !!
      !! **The batch can be wide.** The accumulator is one `n_set * n^2` array
      !! with the set index fastest, shared by every thread and written under
      !! a lock per tile pair, and the quartet loop is tiled so the blocks a
      !! tile quartet touches stay in cache whatever the width; a pass over
      !! 222 densities costs a few times a pass over 12. Memory is the
      !! densities, the accumulator and the result, three copies of
      !! `n_set * n^2`, whatever the thread count.
      !!
      !! **Symmetric densities only.** The six updates below assume `D` is
      !! symmetric in two places: the factor of two for `s1 /= s2` stands in for
      !! the `mu <-> nu` permutation without ever adding `D(nu,mu)`, and likewise
      !! for `s3 /= s4`. Hand this an antisymmetric density and those
      !! permutations double where they should cancel, so the Coulomb term comes
      !! back at twice its value instead of at zero, and nothing here detects it.
      !! Use `build_fock_direct_nosym` for that case.
!$    use omp_lib, only: omp_lock_kind, omp_init_lock, omp_set_lock, omp_unset_lock, omp_destroy_lock
      type(czt_molecule_t), intent(in) :: mol
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
      real(dp), intent(in), optional :: j_scale
         !! Fraction of Coulomb, one by default. A long-range exchange pass
         !! passes zero: the full-range pass has already supplied it.
      real(dp), intent(in), optional :: omega
         !! Range separation, through `env(PTR_RANGE_OMEGA)`.
      logical, intent(in), optional :: density_screen
         !! Also skip a quartet when the largest density element it multiplies,
         !! over every set, makes its contribution negligible. Off by default,
         !! which keeps a batch bit-for-bit what the same densities give one at
         !! a time; on, the bound is `Q_ij Q_kl max|D|` over the six elements
         !! the updates read, rigorous for any symmetric matrix, and what it
         !! buys is locality -- a response to one atom's displacement is small
         !! far from that atom.

      real(dp), allocatable :: buf(:), g(:, :, :), d_half(:, :, :)
      real(dp), allocatable :: bq(:, :), dsh(:, :), dsh_set(:, :), dsh_all(:, :, :)
      real(dp), allocatable :: l12(:, :, :), l34(:, :, :), l13(:, :, :), l24(:, :, :), l14(:, :, :), l23(:, :, :)
      integer, allocatable :: active(:), hits(:)
      logical, allocatable :: hit_ts(:), hit_task(:)
      logical :: all_ts, all_task
      integer :: nhit, fp, fq, fr, fs, wp, wq, wr, ws, tw_max
      integer, allocatable :: fo(:), tw(:)
!$    integer(omp_lock_kind), allocatable :: locks(:, :)
      integer :: na, ia, p12, p34, p13, p24, p14, p23
      real(dp), allocatable :: dd12(:, :), dd34(:, :), dd13(:, :), dd24(:, :), dd14(:, :), dd23(:, :)
      real(dp), allocatable :: gg12(:, :), gg34(:, :), gg13(:, :), gg24(:, :), gg14(:, :), gg23(:, :)
      real(dp) :: qq, kq
      logical :: weight_density
      type(eri_shell_table_t) :: tab
      type(c_ptr) :: opt
      integer :: s1, s2, s3, s4
      integer :: d1, d2, d3, d4, o1, o2, o3, o4
      integer :: shls(4)
      integer :: f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, block_max, n, iset, n_set
      integer :: ij, kl
      integer, allocatable :: dims(:), offs(:)
      integer :: itask
      integer :: ntile, ntp, tp, tq, tr, ts
      integer, allocatable :: tile_first(:), tile_last(:), tp_p(:), tp_q(:), tp_r(:)
      integer(int64) :: n_total, n_computed, n_screened
      real(dp) :: tol, deg, value, scaled, jscaled, kscaled, kx, jxm
      real(dp), allocatable :: env_many(:)

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
      ! one of the two drivers that read a fused L shell, and is exactly what
      ! this loop calls. The Schwarz bounds arrive per split shell either way and
      ! are collapsed to match. Dimensions and offsets are copied out up front:
      ! both are read inside the parallel region.
      call eri_shell_table(mol, tab)
      call eri_schwarz_collapse(mol, bounds, bq)
      dims = tab%dims
      offs = tab%offs(1:tab%nbas)
      block_max = tab%block_max

      ! **The quartet loop is tiled for the cache.** With the set index fastest
      ! a block of the accumulator or the density is `n_set` doubles per
      ! matrix element, so at a few hundred sets the six blocks one shell
      ! quartet touches are hundreds of kilobytes and a task that walks every
      ! ket pair below its bra pair streams them from memory once per quartet.
      ! Measured on cholesterol at 222 sets: the loop took the same 30 s on
      ! 16 threads and on 128, with 93 per cent of the samples in the six
      ! updates, which is a pass bound by memory bandwidth and not by the
      ! integrals. Shells are grouped into tiles of consecutive shells of
      ! about `FOCK_TILE_FUNCS` functions -- roughly one heavy atom -- and a
      ! task is one bra tile pair against every ket tile pair at or below it,
      ! shell quartets innermost. The twelve blocks a tile quartet touches fit
      ! the cache and every shell quartet inside it reuses them.
      !
      ! Coverage: a quartet `(ij, kl)` with `kl <= ij` has `s3 <= s1`, so its
      ! ket tile pair has `tile(s3) <= tile(s1)`; walking every `(r, t)` with
      ! `r <= p` and `t <= r` and keeping `kl <= ij` visits each unique
      ! quartet exactly once, the same set as the triangular nest.
      call shell_tiles(dims, FOCK_TILE_FUNCS, tile_first, tile_last, ntile)
      ! First function of each tile, one before it, and its width in functions.
      allocate (fo(ntile), tw(ntile))
      do tp = 1, ntile
         fo(tp) = offs(tile_first(tp))
         tw(tp) = offs(tile_last(tp)) + dims(tile_last(tp)) - fo(tp)
      end do
      tw_max = maxval(tw)
      ! A task is a bra tile pair and one ket bra tile: the blocks it
      ! holds on the bra side and on `tr` stay put while `ts` walks.
      ntp = 0
      do tp = 1, ntile
         ntp = ntp + tp*tp
      end do
      allocate (tp_p(ntp), tp_q(ntp), tp_r(ntp))
      itask = 0
      do tp = 1, ntile
         do tq = 1, tp
            do tr = 1, tp
               itask = itask + 1
               tp_p(itask) = tp
               tp_q(itask) = tq
               tp_r(itask) = tr
            end do
         end do
      end do

      ! **Set index fastest**, in the densities and in the accumulator both.
      ! The six updates per integral element land on six matrix elements, and
      ! with the set slowest each of those is `n_set` scattered read-modify-
      ! writes a stride of `n^2` apart -- a different page per set. With the
      ! set fastest they are six contiguous vectors of length `n_set`, which is
      ! what makes a wide batch cost little more than a narrow one. Measured on
      ! the Hessian's coupled-perturbed solve, where the cost of a pass grew
      ! almost linearly with the batch in the other layout.
      allocate (d_half(n_set, n, n))

      ! The six-update form below is written for D = C_occ C_occ^T, the density
      ! without its factor of two. `build_density` produces D = 2 C_occ C_occ^T,
      ! so halve it here. Skipping this makes both J and K exactly twice too
      ! large -- an error that still converges, to a badly wrong energy.
      kx = 1.0_dp
      if (present(k_scale)) kx = k_scale
      jxm = 1.0_dp
      if (present(j_scale)) jxm = j_scale
      do iset = 1, n_set
         d_half(iset, :, :) = 0.5_dp*densities(:, :, iset)
      end do

      ! The largest halved density element per shell pair, over every set,
      ! which is what the screen below multiplies the Schwarz bound by.
      weight_density = .false.
      if (present(density_screen)) weight_density = density_screen
      if (weight_density) then
         allocate (dsh(tab%nbas, tab%nbas), dsh_all(n_set, tab%nbas, tab%nbas))
         dsh = 0.0_dp
         do iset = 1, n_set
            call block_density_max(0.5_dp*densities(:, :, iset), tab%nbas, offs, dims, dsh_set)
            dsh = max(dsh, dsh_set)
            dsh_all(iset, :, :) = dsh_set
            deallocate (dsh_set)
         end do
      else
         allocate (dsh(0, 0), dsh_all(0, 0, 0))
      end if
      kq = kx*0.25_dp

      ! Created once and reused for every quartet in this build.
      opt = c_null_ptr
      ! As `build_fock_direct`: a local environment so the attenuated pass is
      ! the same quartets with the range-separation slot set.
      env_many = tab%env
      if (present(omega)) env_many(LIBCINT_PTR_RANGE_OMEGA + 1) = omega

      call two_electron_optimizer(mol%cartesian, opt, mol%atm, mol%natm, tab%bas, &
                                  tab%nbas, env_many)

      n_total = 0_int64
      n_computed = 0_int64
      n_screened = 0_int64

      ! Threaded over tile tasks. libcint carries no mutable state across calls --
      ! the 2e path has no static globals, and `opt` is written once here and
      ! only read inside -- so the integrals need nothing but a private buffer
      ! each.
      !
      ! **The accumulator is one shared array, written under a lock per tile
      ! pair.** The updates below scatter at positions that depend on all four
      ! shells, so two threads holding different quartets can land on the same
      ! element. A copy per thread was the first answer and its memory was the
      ! thread count times `n_set * n^2`: 60 GB for 128 threads on cholesterol
      ! at 222 sets, and it would not have fit a 256 GB node at 256 threads.
      ! Instead a thread accumulates the six blocks of the tile quartet in hand
      ! in six tile-sized local blocks -- the same blocks the tiling already
      ! keeps in cache -- and adds each into the shared array once the tile
      ! quartet is done, holding that tile pair's lock while it does. The three
      ! blocks fixed across a task, `(p,q)`, `(p,r)` and `(q,r)`, wait until the
      ! task ends, since every task with the same bra tile pair wants `(p,q)`
      ! and they would otherwise queue on it. Only the sets a block touched
      ! are added and zeroed, so a nearly screened tile quartet costs nearly
      ! nothing to flush. The memory traffic is what the per-thread copies
      ! already paid, since their blocks were evicted once per tile quartet
      ! too; what is gone is the reduction over the copies.
      !
      ! `schedule(dynamic)`: task `(p, q, r)` walks `r` ket tiles and the
      ! screen leaves the distant ones nearly empty, so the tasks are uneven.
      allocate (g(n_set, n, n))
!$    allocate (locks(ntile, ntile))
!$    do tq = 1, ntile
!$       do tp = 1, ntile
!$          call omp_init_lock(locks(tp, tq))
!$       end do
!$    end do

      !$omp parallel default(none) &
      !$omp    shared(kx, jxm, env_many, mol, tab, bq, dsh, weight_density, d_half, g, dims, offs, &
      !$omp           tile_first, tile_last, tp_p, tp_q, tp_r, ntp, tol, opt, n, block_max, n_set, &
      !$omp           focks, h, dsh_all, kq, locks, fo, tw, tw_max) &
      !$omp    private(itask, ij, kl, tp, tq, tr, ts, s1, s2, s3, s4, d1, d2, d3, d4, o1, o2, o3, o4, &
      !$omp            shls, f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, deg, value, scaled, &
      !$omp            jscaled, kscaled, buf, iset, active, na, ia, qq, &
      !$omp            p12, p34, p13, p24, p14, p23, dd12, dd34, dd13, dd24, dd14, dd23, &
      !$omp            gg12, gg34, gg13, gg24, gg14, gg23, &
      !$omp            l12, l34, l13, l24, l14, l23, hits, hit_ts, hit_task, all_ts, all_task, nhit, &
      !$omp            fp, fq, fr, fs, wp, wq, wr, ws) &
      !$omp    reduction(+:n_total, n_computed, n_screened)
      allocate (buf(block_max**4), active(n_set), hits(n_set), hit_ts(n_set), hit_task(n_set))
      active = [(iset, iset=1, n_set)]
      allocate (dd12(n_set, block_max**2), dd34(n_set, block_max**2), dd13(n_set, block_max**2), &
                dd24(n_set, block_max**2), dd14(n_set, block_max**2), dd23(n_set, block_max**2))
      allocate (gg12(n_set, block_max**2), gg34(n_set, block_max**2), gg13(n_set, block_max**2), &
                gg24(n_set, block_max**2), gg14(n_set, block_max**2), gg23(n_set, block_max**2))
      ! The six tile-local blocks, zero on entry and left zero by every flush.
      allocate (l12(n_set, tw_max, tw_max), l34(n_set, tw_max, tw_max), l13(n_set, tw_max, tw_max), &
                l24(n_set, tw_max, tw_max), l14(n_set, tw_max, tw_max), l23(n_set, tw_max, tw_max))
      l12 = 0.0_dp
      l34 = 0.0_dp
      l13 = 0.0_dp
      l24 = 0.0_dp
      l14 = 0.0_dp
      l23 = 0.0_dp

      ! The shared accumulator, zeroed by every thread so its pages are spread.
      !$omp do
      do b2 = 1, n
         g(:, :, b2) = 0.0_dp
      end do
      !$omp end do

      ! Descending: the later tile pairs carry the most ket tile pairs, and
      ! `schedule(dynamic)` wants the big tasks first.
      !$omp do schedule(dynamic)
      do itask = ntp, 1, -1
         tp = tp_p(itask)
         tq = tp_q(itask)
         tr = tp_r(itask)
         fp = fo(tp)
         fq = fo(tq)
         fr = fo(tr)
         wp = tw(tp)
         wq = tw(tq)
         wr = tw(tr)
         hit_task = .false.
         all_task = .false.
         do ts = 1, tr
            fs = fo(ts)
            ws = tw(ts)
            hit_ts = .false.
            all_ts = .false.
            do s1 = tile_first(tp), tile_last(tp)
               d1 = dims(s1)
               o1 = offs(s1)
               do s2 = tile_first(tq), min(tile_last(tq), s1)
                  d2 = dims(s2)
                  o2 = offs(s2)
                  ij = s1*(s1 - 1)/2 + s2
                  do s3 = tile_first(tr), min(tile_last(tr), s1)
                     d3 = dims(s3)
                     o3 = offs(s3)
                     do s4 = tile_first(ts), min(tile_last(ts), s3)
                        kl = s3*(s3 - 1)/2 + s4
                        if (kl > ij) cycle
                        d4 = dims(s4)
                        o4 = offs(s4)

                        n_total = n_total + 1_int64

                        if (bq(s1, s2)*bq(s3, s4) < tol) then
                           n_screened = n_screened + 1_int64
                           cycle
                        end if
                        deg = pair_degeneracy(s1, s2, s3, s4)
                        if (weight_density) then
                           if (bq(s1, s2)*bq(s3, s4)*density_weight(dsh, s1, s2, s3, s4, jxm, &
                                                                    kx*0.25_dp, deg) < tol) then
                              n_screened = n_screened + 1_int64
                              cycle
                           end if
                        end if

                        shls = [s1 - 1, s2 - 1, s3 - 1, s4 - 1]
                        ! `env_many`, not `tab%env`: the omega slot lives in the copy, and
                        ! handing the optimizer the attenuated environment while the quartet
                        ! call reads the unattenuated one gives full-range integrals scaled
                        ! by `k_scale` -- a long-range pass that is not long-range at all.
                        ret = two_electron_block(mol%cartesian, buf, shls, mol%atm, mol%natm, &
                                                 tab%bas, tab%nbas, env_many, opt)
                        if (ret == 0) then
                           n_screened = n_screened + 1_int64
                           cycle
                        end if
                        n_computed = n_computed + 1_int64

                        ! Which sets this quartet can touch above the tolerance: the same
                        ! bound as the quartet screen, per set. A response density is
                        ! local to its atom, so with a wide batch most sets fail it and
                        ! the six updates run over the survivors alone.
                        if (weight_density) then
                           na = 0
                           qq = deg*bq(s1, s2)*bq(s3, s4)
                           do iset = 1, n_set
                              if (qq*max(jxm*dsh_all(iset, s1, s2), jxm*dsh_all(iset, s3, s4), &
                                         kq*dsh_all(iset, s1, s3), kq*dsh_all(iset, s1, s4), &
                                         kq*dsh_all(iset, s2, s3), kq*dsh_all(iset, s2, s4)) >= tol) then
                                 na = na + 1
                                 active(na) = iset
                              end if
                           end do
                        else
                           na = n_set
                        end if
                        if (na == 0) cycle

                        ! Three shapes of the same six updates. At least half the sets
                        ! active, which is the usual case: contiguous vectors down the
                        ! whole set index straight into the accumulator, the inactive
                        ! sets included -- their contributions are below tolerance and
                        ! cost less as part of a vector than a scalar loop over the
                        ! survivors would, a measured ten times less per quartet. A
                        ! small sparse quartet: the updates indexed over the survivors,
                        ! since a handful of elements does not repay gathering blocks.
                        ! A large sparse quartet: the density blocks compacted to the
                        ! survivors once, the updates contiguous over them, and the six
                        ! output blocks scattered back once.
                        if (2*na >= n_set) na = n_set
                        if (na == n_set) then
                           all_ts = .true.
                           all_task = .true.
                        else
                           do ia = 1, na
                              hit_ts(active(ia)) = .true.
                              hit_task(active(ia)) = .true.
                           end do
                        end if
                        if (na == n_set .or. d1*d2*d3*d4 < 32) then
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
                                       jscaled = jxm*scaled
                                       kscaled = kx*0.25_dp*scaled

                                       ! Two Coulomb and four exchange contributions. g is
                                       ! not symmetric as it stands; symmetrising at the end
                                       ! is what makes these six updates equivalent to the
                                       ! full sum over all eight permutations.
                                       if (na == n_set) then
                                          l12(:, b1 - fp, b2 - fq) = l12(:, b1 - fp, b2 - fq) + jscaled*d_half(:, b3, b4)
                                          l34(:, b3 - fr, b4 - fs) = l34(:, b3 - fr, b4 - fs) + jscaled*d_half(:, b1, b2)
                                          l13(:, b1 - fp, b3 - fr) = l13(:, b1 - fp, b3 - fr) - kscaled*d_half(:, b2, b4)
                                          l24(:, b2 - fq, b4 - fs) = l24(:, b2 - fq, b4 - fs) - kscaled*d_half(:, b1, b3)
                                          l14(:, b1 - fp, b4 - fs) = l14(:, b1 - fp, b4 - fs) - kscaled*d_half(:, b2, b3)
                                          l23(:, b2 - fq, b3 - fr) = l23(:, b2 - fq, b3 - fr) - kscaled*d_half(:, b1, b4)
                                       else
                                          do ia = 1, na
                                             iset = active(ia)
                                  l12(iset, b1 - fp, b2 - fq) = l12(iset, b1 - fp, b2 - fq) + jscaled*d_half(iset, b3, b4)
                                  l34(iset, b3 - fr, b4 - fs) = l34(iset, b3 - fr, b4 - fs) + jscaled*d_half(iset, b1, b2)
                                  l13(iset, b1 - fp, b3 - fr) = l13(iset, b1 - fp, b3 - fr) - kscaled*d_half(iset, b2, b4)
                                  l24(iset, b2 - fq, b4 - fs) = l24(iset, b2 - fq, b4 - fs) - kscaled*d_half(iset, b1, b3)
                                  l14(iset, b1 - fp, b4 - fs) = l14(iset, b1 - fp, b4 - fs) - kscaled*d_half(iset, b2, b3)
                                  l23(iset, b2 - fq, b3 - fr) = l23(iset, b2 - fq, b3 - fr) - kscaled*d_half(iset, b1, b4)
                                          end do
                                       end if
                                    end do
                                 end do
                              end do
                           end do
                           cycle
                        end if

                        ! The six density blocks this quartet reads, compacted to the
                        ! active sets as `(na, d_a*d_b)` with the first function fastest.
                        ! The six output blocks accumulate in the same shape and go back
                        ! to the accumulator once per quartet, so the updates per element
                        ! are contiguous whatever the sets, and land in a few kilobytes.
                        do f2 = 1, d2
                           do f1 = 1, d1
                              p12 = f1 + (f2 - 1)*d1
                              do ia = 1, na
                                 dd12(ia, p12) = d_half(active(ia), o1 + f1, o2 + f2)
                              end do
                           end do
                        end do
                        do f4 = 1, d4
                           do f3 = 1, d3
                              p34 = f3 + (f4 - 1)*d3
                              do ia = 1, na
                                 dd34(ia, p34) = d_half(active(ia), o3 + f3, o4 + f4)
                              end do
                           end do
                        end do
                        do f3 = 1, d3
                           do f1 = 1, d1
                              p13 = f1 + (f3 - 1)*d1
                              do ia = 1, na
                                 dd13(ia, p13) = d_half(active(ia), o1 + f1, o3 + f3)
                              end do
                           end do
                        end do
                        do f4 = 1, d4
                           do f2 = 1, d2
                              p24 = f2 + (f4 - 1)*d2
                              do ia = 1, na
                                 dd24(ia, p24) = d_half(active(ia), o2 + f2, o4 + f4)
                              end do
                           end do
                        end do
                        do f4 = 1, d4
                           do f1 = 1, d1
                              p14 = f1 + (f4 - 1)*d1
                              do ia = 1, na
                                 dd14(ia, p14) = d_half(active(ia), o1 + f1, o4 + f4)
                              end do
                           end do
                        end do
                        do f3 = 1, d3
                           do f2 = 1, d2
                              p23 = f2 + (f3 - 1)*d2
                              do ia = 1, na
                                 dd23(ia, p23) = d_half(active(ia), o2 + f2, o3 + f3)
                              end do
                           end do
                        end do
                        gg12(1:na, 1:d1*d2) = 0.0_dp
                        gg34(1:na, 1:d3*d4) = 0.0_dp
                        gg13(1:na, 1:d1*d3) = 0.0_dp
                        gg24(1:na, 1:d2*d4) = 0.0_dp
                        gg14(1:na, 1:d1*d4) = 0.0_dp
                        gg23(1:na, 1:d2*d3) = 0.0_dp

                        do f4 = 1, d4
                           do f3 = 1, d3
                              p34 = f3 + (f4 - 1)*d3
                              do f2 = 1, d2
                                 p24 = f2 + (f4 - 1)*d2
                                 p23 = f2 + (f3 - 1)*d2
                                 do f1 = 1, d1
                                    p12 = f1 + (f2 - 1)*d1
                                    p13 = f1 + (f3 - 1)*d1
                                    p14 = f1 + (f4 - 1)*d1

                                    idx = f1 + (f2 - 1)*d1 + (f3 - 1)*d1*d2 + (f4 - 1)*d1*d2*d3
                                    value = buf(idx)
                                    scaled = value*deg
                                    jscaled = jxm*scaled
                                    kscaled = kx*0.25_dp*scaled

                                    ! Two Coulomb and four exchange contributions, every
                                    ! active set at once down the fast index. g is not
                                    ! symmetric as it stands; symmetrising at the end is
                                    ! what makes these six updates equivalent to the full
                                    ! sum over all eight permutations.
                                    !
                                    ! These six lines are the entire point of the routine:
                                    ! the integral above was computed once and every set
                                    ! reuses it.
                                    gg12(1:na, p12) = gg12(1:na, p12) + jscaled*dd34(1:na, p34)
                                    gg34(1:na, p34) = gg34(1:na, p34) + jscaled*dd12(1:na, p12)
                                    gg13(1:na, p13) = gg13(1:na, p13) - kscaled*dd24(1:na, p24)
                                    gg24(1:na, p24) = gg24(1:na, p24) - kscaled*dd13(1:na, p13)
                                    gg14(1:na, p14) = gg14(1:na, p14) - kscaled*dd23(1:na, p23)
                                    gg23(1:na, p23) = gg23(1:na, p23) - kscaled*dd14(1:na, p14)
                                 end do
                              end do
                           end do
                        end do

                        ! Added, not assigned: two of the six blocks coincide whenever the
                        ! quartet has a repeated shell pair.
                        do f2 = 1, d2
                           do f1 = 1, d1
                              p12 = f1 + (f2 - 1)*d1
                              do ia = 1, na
                 l12(active(ia), o1 + f1 - fp, o2 + f2 - fq) = l12(active(ia), o1 + f1 - fp, o2 + f2 - fq) + gg12(ia, p12)
                              end do
                           end do
                        end do
                        do f4 = 1, d4
                           do f3 = 1, d3
                              p34 = f3 + (f4 - 1)*d3
                              do ia = 1, na
                 l34(active(ia), o3 + f3 - fr, o4 + f4 - fs) = l34(active(ia), o3 + f3 - fr, o4 + f4 - fs) + gg34(ia, p34)
                              end do
                           end do
                        end do
                        do f3 = 1, d3
                           do f1 = 1, d1
                              p13 = f1 + (f3 - 1)*d1
                              do ia = 1, na
                 l13(active(ia), o1 + f1 - fp, o3 + f3 - fr) = l13(active(ia), o1 + f1 - fp, o3 + f3 - fr) + gg13(ia, p13)
                              end do
                           end do
                        end do
                        do f4 = 1, d4
                           do f2 = 1, d2
                              p24 = f2 + (f4 - 1)*d2
                              do ia = 1, na
                 l24(active(ia), o2 + f2 - fq, o4 + f4 - fs) = l24(active(ia), o2 + f2 - fq, o4 + f4 - fs) + gg24(ia, p24)
                              end do
                           end do
                        end do
                        do f4 = 1, d4
                           do f1 = 1, d1
                              p14 = f1 + (f4 - 1)*d1
                              do ia = 1, na
                 l14(active(ia), o1 + f1 - fp, o4 + f4 - fs) = l14(active(ia), o1 + f1 - fp, o4 + f4 - fs) + gg14(ia, p14)
                              end do
                           end do
                        end do
                        do f3 = 1, d3
                           do f2 = 1, d2
                              p23 = f2 + (f3 - 1)*d2
                              do ia = 1, na
                 l23(active(ia), o2 + f2 - fq, o3 + f3 - fr) = l23(active(ia), o2 + f2 - fq, o3 + f3 - fr) + gg23(ia, p23)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do

            ! The three blocks that change with the ket tile go in now.
            !
            ! **Why one lock per tile PAIR is enough, and why the pair is not
            ! sorted.** `add_tile_block` puts `blk` at `g(:, fa+1:, fb+1:)`, so
            ! the block a call writes is exactly the ordered pair it is given,
            ! and every call below takes the lock of that same ordered pair.
            ! The task enumeration constrains `tq <= tp` and `tr <= tp` but
            ! leaves `tq` and `tr` unordered against each other, so both
            ! `(tq, tr)` and `(tr, tq)` occur across tasks -- they address
            ! different halves of a full `n x n` accumulator and so are
            ! different memory under different locks. Canonicalising the index
            ! would be wrong here, not merely unnecessary: it would put two
            ! distinct blocks behind one lock.
            !
            ! Each lock is also taken and released before the next, so there is
            ! no hold-and-wait and no order to get wrong.
            call touched_sets(all_ts, hit_ts, hits, nhit)
            if (nhit > 0) then
!$             call omp_set_lock(locks(tr, ts))
               call add_tile_block(g, l34, wr, ws, fr, fs, nhit, hits)
!$             call omp_unset_lock(locks(tr, ts))
!$             call omp_set_lock(locks(tq, ts))
               call add_tile_block(g, l24, wq, ws, fq, fs, nhit, hits)
!$             call omp_unset_lock(locks(tq, ts))
!$             call omp_set_lock(locks(tp, ts))
               call add_tile_block(g, l14, wp, ws, fp, fs, nhit, hits)
!$             call omp_unset_lock(locks(tp, ts))
            end if
         end do

         ! The three the task held throughout.
         call touched_sets(all_task, hit_task, hits, nhit)
         if (nhit > 0) then
!$          call omp_set_lock(locks(tp, tq))
            call add_tile_block(g, l12, wp, wq, fp, fq, nhit, hits)
!$          call omp_unset_lock(locks(tp, tq))
!$          call omp_set_lock(locks(tp, tr))
            call add_tile_block(g, l13, wp, wr, fp, fr, nhit, hits)
!$          call omp_unset_lock(locks(tp, tr))
!$          call omp_set_lock(locks(tq, tr))
            call add_tile_block(g, l23, wq, wr, fq, fr, nhit, hits)
!$          call omp_unset_lock(locks(tq, tr))
         end if
      end do
      !$omp end do

      ! Back to one matrix per set, symmetrised on the way. The barrier the
      ! loop above ends with is what makes every block of `g` complete here.
      !$omp do
      do b2 = 1, n
         do iset = 1, n_set
            do b1 = 1, n
               focks(b1, b2, iset) = h(b1, b2) + 0.5_dp*(g(iset, b1, b2) + g(iset, b2, b1))
            end do
         end do
      end do
      !$omp end do

      deallocate (buf, active, hits, hit_ts, hit_task, dd12, dd34, dd13, dd24, dd14, dd23, &
                  gg12, gg34, gg13, gg24, gg14, gg23, l12, l34, l13, l24, l14, l23)
      !$omp end parallel

      stats%quartets_total = n_total
      stats%quartets_computed = n_computed
      stats%quartets_screened = n_screened

      call libcint_del_optimizer(opt)

!$    do tq = 1, ntile
!$       do tp = 1, ntile
!$          call omp_destroy_lock(locks(tp, tq))
!$       end do
!$    end do
!$    deallocate (locks)
      deallocate (g, d_half, dims, offs, tile_first, tile_last, tp_p, tp_q, tp_r, dsh, dsh_all, fo, tw)
   end subroutine build_fock_direct_many

   pure subroutine touched_sets(all, hit, hits, nhit)
      !! The list of sets a tile block was written in
      !!
      !! `all` short-circuits the mask: every set, and `hits` is not read by
      !! the flush in that case.
      logical, intent(in) :: all
      logical, intent(in) :: hit(:)
      integer, intent(inout) :: hits(:)
      integer, intent(out) :: nhit

      integer :: iset

      if (all) then
         nhit = size(hit)
         return
      end if
      nhit = 0
      do iset = 1, size(hit)
         if (hit(iset)) then
            nhit = nhit + 1
            hits(nhit) = iset
         end if
      end do
   end subroutine touched_sets

   subroutine add_tile_block(g, blk, wa, wb, fa, fb, nhit, hits)
      !! Adds a tile-local block into the shared accumulator and clears it
      !!
      !! `blk(:, 1:wa, 1:wb)` lands at `g(:, fa+1:fa+wa, fb+1:fb+wb)`. Only the
      !! `nhit` sets in `hits` are added and zeroed -- the rest were never
      !! written and are zero already -- and `nhit == size(g, 1)` means every
      !! set, contiguously, without reading `hits`. The caller holds the lock
      !! of the tile pair `(fa, fb)` belongs to.
      real(dp), intent(inout) :: g(:, :, :)
      real(dp), intent(inout) :: blk(:, :, :)
      integer, intent(in) :: wa, wb, fa, fb, nhit
      integer, intent(in) :: hits(:)

      integer :: a, b, ih, iset

      if (nhit == size(g, 1)) then
         do b = 1, wb
            do a = 1, wa
               g(:, fa + a, fb + b) = g(:, fa + a, fb + b) + blk(:, a, b)
               blk(:, a, b) = 0.0_dp
            end do
         end do
      else
         do b = 1, wb
            do a = 1, wa
               do ih = 1, nhit
                  iset = hits(ih)
                  g(iset, fa + a, fb + b) = g(iset, fa + a, fb + b) + blk(iset, a, b)
                  blk(iset, a, b) = 0.0_dp
               end do
            end do
         end do
      end if
   end subroutine add_tile_block

   subroutine build_fock_direct_nosym(mol, h, densities, bounds, focks, stats, error, &
                                      screen_tol)
      !! J - K/2 for densities of **any** symmetry, over one pass of the integrals
      !!
      !! `build_fock_direct_many` is faster and cannot be used here: it folds
      !! three of the eightfold permutations into a `deg` factor, which is only
      !! the same thing when `D` is symmetric. The error is in the accumulation,
      !! so no change to the final symmetrisation repairs it.
      !!
      !! This routine writes the permutations out instead. For each computed
      !! integral it generates the distinct index tuples the block enumeration
      !! does not already cover -- exactly the ones `deg` was standing in for --
      !! and applies to each
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
      !! The frequency-dependent coupled-perturbed equations need `A - B` as well
      !! as `A + B`, and both are this same contraction: symmetrise the response
      !! density for `A + B`, antisymmetrise it for `A - B`. The Coulomb term
      !! dropping out of the second happens by itself, because the two-electron
      !! integral is symmetric in its ket pair while the density is not.
      !!
      !! **Cost.** Up to eight tuples times two updates against six updates, so
      !! roughly 2.7x the contraction work per integral. Prefer two passes --
      !! symmetric densities through `build_fock_direct_many`, antisymmetric ones
      !! through here -- over routing everything through this one.
      type(czt_molecule_t), intent(in) :: mol
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
      ! one of the two drivers that read a fused L shell, and is exactly what
      ! this loop calls. The Schwarz bounds arrive per split shell either way and
      ! are collapsed to match. Dimensions and offsets are copied out up front:
      ! both are read inside the parallel region.
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
      !! the updates read: Coulomb draws on half the total density and exchange
      !! on the same-spin density, which in the closed-shell case are the same
      !! matrix, D/2.
      !!
      !! **The two loops have to stay in step** -- this one is a copy rather than
      !! a flag on `build_fock_direct`, so the hot closed-shell path carries one
      !! accumulator. The elementwise check against the in-core UHF build is what
      !! enforces the agreement.
      type(czt_molecule_t), intent(in) :: mol
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
      ! one of the two drivers that read a fused L shell, and is exactly what
      ! this loop calls. The Schwarz bounds arrive per split shell either way and
      ! are collapsed to match. Dimensions and offsets are copied out up front:
      ! both are read inside the parallel region.
      call eri_shell_table(mol, tab)
      call eri_schwarz_collapse(mol, bounds, bq)
      dims = tab%dims
      offs = tab%offs(1:tab%nbas)
      block_max = tab%block_max

      ! The quartet loop, flattened onto one index so it can be handed out.
      ! Enumerating shell pairs and taking `kl <= ij` covers the same quartets
      ! once each as the triangular nest, and is divisible where that is not.
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

      ! Coulomb sees half the total density, exchange sees the same-spin one. The
      ! closed-shell build passes D/2 to both because there D_alpha = D_beta =
      ! D/2, so the two forms are one form.
      !
      ! `j_scale` rides on this matrix rather than on the Coulomb updates, which
      ! is the same arithmetic one level out.
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

      ! Threaded over bra pairs. libcint carries no mutable state across calls --
      ! the 2e path has no static globals, and `opt` is written once here and
      ! only read inside -- so the integrals need nothing but a private buffer
      ! each.
      !
      ! The accumulator does need care. The updates below scatter at positions
      ! that depend on all four shells, so two threads holding
      ! different quartets can land on the same element. Each thread fills its
      ! own copy and adds it in once at the end; atomics on the innermost
      ! statement would be correct too and far slower. The cost is one n*n array
      ! per thread, growing as n^2 -- **worth watching at a few thousand
      ! functions**, where a blocked scheme like GTFock's would be the answer.
      !
      ! `schedule(dynamic)`: pair ij does ij quartets, so the last chunk is
      ! thousands of times the first.
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
            ! exceeds one -- `deg` reaches eight. The two counters only attribute
            ! the decision.
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

                        ! Two Coulomb and four exchange contributions per spin.
                        ! `ga` and `gb` are not symmetric as they stand;
                        ! symmetrising at the end is what makes these updates
                        ! equivalent to the full sum over all eight permutations.
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

end module mqc_czt_direct
