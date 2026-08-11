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
   use pic_types, only: dp, int64
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim, &
                                    two_electron_block, two_electron_optimizer
   use libcint_fortran, only: libcint_del_optimizer
   implicit none
   private

   public :: schwarz_bounds
   public :: build_fock_direct
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
      integer(int64) :: quartets_screened = 0  !! Skipped on the Schwarz bound
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

   subroutine build_fock_direct(mol, h, density, bounds, fock, stats, error, screen_tol)
      !! F = H + J - K/2, without forming the integral tensor
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h(:, :)          !! Core Hamiltonian
      real(dp), intent(in) :: density(:, :)    !! D = 2 C_occ C_occ^T
      real(dp), intent(in) :: bounds(:, :)     !! From `schwarz_bounds`
      real(dp), intent(out) :: fock(:, :)
      type(direct_stats_t), intent(out) :: stats
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: screen_tol

      real(dp), allocatable :: buf(:), g(:, :), g_local(:, :), d_half(:, :)
      type(c_ptr) :: opt
      integer :: s1, s2, s3, s4
      integer :: d1, d2, d3, d4, o1, o2, o3, o4
      integer :: shls(4)
      integer :: f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, block_max, n
      integer :: ij, kl, npair, ipair
      integer, allocatable :: pair_i(:), pair_j(:), dims(:), offs(:)
      integer(int64) :: n_total, n_computed, n_screened
      real(dp) :: tol, deg, value, scaled

      n = mol%nao
      if (size(h, 1) /= n .or. size(density, 1) /= n .or. size(fock, 1) /= n) then
         call error%set(ERROR_VALIDATION, "direct Fock: matrix dimensions do not match the basis")
         return
      end if

      tol = DEFAULT_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol

      ! Shell dimensions and offsets up front. Both are needed inside the
      ! parallel region, and looking them up there would mean every thread
      ! reaching into `bas` for something that does not change.
      allocate (dims(mol%nbas), offs(mol%nbas))
      block_max = 1
      do s1 = 1, mol%nbas
         dims(s1) = shell_dim(mol%cartesian, s1 - 1, mol%bas)
         offs(s1) = mol%shell_offset(s1)
         block_max = max(block_max, dims(s1))
      end do

      ! The quartet loop, flattened onto one index so it can be handed out.
      !
      ! The nested form -- s3 up to s1, and s4 up to s2 only when s3 equals s1 --
      ! is exactly "every shell pair (s3,s4) at or before (s1,s2)" in the
      ! canonical pair ordering, so enumerating pairs and taking kl <= ij covers
      ! the same quartets once each. Flattening is what makes the outer loop
      ! divisible; the triangular nest is not.
      npair = mol%nbas*(mol%nbas + 1)/2
      allocate (pair_i(npair), pair_j(npair))
      ipair = 0
      do s1 = 1, mol%nbas
         do s2 = 1, s1
            ipair = ipair + 1
            pair_i(ipair) = s1
            pair_j(ipair) = s2
         end do
      end do

      allocate (g(n, n), d_half(n, n))
      g = 0.0_dp

      ! The six-update form below is written for D = C_occ C_occ^T, the density
      ! without its factor of two. `build_density` produces D = 2 C_occ C_occ^T,
      ! so halve it here. Skipping this makes both J and K exactly twice too
      ! large -- an error that still converges, to a badly wrong energy.
      d_half = 0.5_dp*density

      ! Created once and reused for every quartet in this build.
      opt = c_null_ptr
      call two_electron_optimizer(mol%cartesian, opt, mol%atm, mol%natm, mol%bas, &
                                  mol%nbas, mol%env)

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
      !$omp    shared(mol, bounds, d_half, g, dims, offs, pair_i, pair_j, npair, tol, opt, n, block_max) &
      !$omp    private(ij, kl, s1, s2, s3, s4, d1, d2, d3, d4, o1, o2, o3, o4, &
      !$omp            shls, f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, deg, value, scaled, &
      !$omp            buf, g_local) &
      !$omp    reduction(+:n_total, n_computed, n_screened)
      allocate (buf(block_max**4))
      allocate (g_local(n, n))
      g_local = 0.0_dp

      !$omp do schedule(dynamic)
      do ij = 1, npair
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

            if (bounds(s1, s2)*bounds(s3, s4) < tol) then
               n_screened = n_screened + 1_int64
               cycle
            end if

            shls = [s1 - 1, s2 - 1, s3 - 1, s4 - 1]
            ret = two_electron_block(mol%cartesian, buf, shls, mol%atm, mol%natm, &
                                     mol%bas, mol%nbas, mol%env, opt)
            if (ret == 0) then
               n_screened = n_screened + 1_int64
               cycle
            end if
            n_computed = n_computed + 1_int64

            ! Shell equality, not function equality: a block with s1 == s2
            ! already contains both orderings of its function pair.
            deg = 1.0_dp
            if (s1 /= s2) deg = deg*2.0_dp
            if (s3 /= s4) deg = deg*2.0_dp
            if (.not. (s1 == s3 .and. s2 == s4)) deg = deg*2.0_dp

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
                        g_local(b1, b2) = g_local(b1, b2) + d_half(b3, b4)*scaled
                        g_local(b3, b4) = g_local(b3, b4) + d_half(b1, b2)*scaled
                        g_local(b1, b3) = g_local(b1, b3) - 0.25_dp*d_half(b2, b4)*scaled
                        g_local(b2, b4) = g_local(b2, b4) - 0.25_dp*d_half(b1, b3)*scaled
                        g_local(b1, b4) = g_local(b1, b4) - 0.25_dp*d_half(b2, b3)*scaled
                        g_local(b2, b3) = g_local(b2, b3) - 0.25_dp*d_half(b1, b4)*scaled
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

      fock = h + 0.5_dp*(g + transpose(g))

      deallocate (g, d_half, dims, offs, pair_i, pair_j)
   end subroutine build_fock_direct

   subroutine build_fock_direct_uhf(mol, h, d_alpha, d_beta, bounds, fock_a, fock_b, stats, error, screen_tol)
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

      real(dp), allocatable :: buf(:), ga(:, :), gb(:, :), ga_local(:, :), gb_local(:, :)
      real(dp), allocatable :: d_coul(:, :)
      type(c_ptr) :: opt
      integer :: s1, s2, s3, s4
      integer :: d1, d2, d3, d4, o1, o2, o3, o4
      integer :: shls(4)
      integer :: f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, block_max, n
      integer :: ij, kl, npair, ipair
      integer, allocatable :: pair_i(:), pair_j(:), dims(:), offs(:)
      integer(int64) :: n_total, n_computed, n_screened
      real(dp) :: tol, deg, value, scaled

      n = mol%nao
      if (size(h, 1) /= n .or. size(d_alpha, 1) /= n .or. size(d_beta, 1) /= n &
          .or. size(fock_a, 1) /= n .or. size(fock_b, 1) /= n) then
         call error%set(ERROR_VALIDATION, "direct UHF Fock: matrix dimensions do not match the basis")
         return
      end if

      tol = DEFAULT_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol

      ! Shell dimensions and offsets up front. Both are needed inside the
      ! parallel region, and looking them up there would mean every thread
      ! reaching into `bas` for something that does not change.
      allocate (dims(mol%nbas), offs(mol%nbas))
      block_max = 1
      do s1 = 1, mol%nbas
         dims(s1) = shell_dim(mol%cartesian, s1 - 1, mol%bas)
         offs(s1) = mol%shell_offset(s1)
         block_max = max(block_max, dims(s1))
      end do

      ! The quartet loop, flattened onto one index so it can be handed out.
      !
      ! The nested form -- s3 up to s1, and s4 up to s2 only when s3 equals s1 --
      ! is exactly "every shell pair (s3,s4) at or before (s1,s2)" in the
      ! canonical pair ordering, so enumerating pairs and taking kl <= ij covers
      ! the same quartets once each. Flattening is what makes the outer loop
      ! divisible; the triangular nest is not.
      npair = mol%nbas*(mol%nbas + 1)/2
      allocate (pair_i(npair), pair_j(npair))
      ipair = 0
      do s1 = 1, mol%nbas
         do s2 = 1, s1
            ipair = ipair + 1
            pair_i(ipair) = s1
            pair_j(ipair) = s2
         end do
      end do

      allocate (ga(n, n), gb(n, n), d_coul(n, n))
      ga = 0.0_dp
      gb = 0.0_dp

      ! Coulomb sees half the total density, exchange sees the same-spin one.
      !
      ! The closed-shell build passes D/2 to both, and that is not a coincidence
      ! worth losing: there D_alpha = D_beta = D/2, so half the total density is
      ! D/2 and the same-spin density is D/2. The two forms are one form.
      d_coul = 0.5_dp*(d_alpha + d_beta)

      ! Created once and reused for every quartet in this build.
      opt = c_null_ptr
      call two_electron_optimizer(mol%cartesian, opt, mol%atm, mol%natm, mol%bas, &
                                  mol%nbas, mol%env)

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
      !$omp    shared(mol, bounds, d_coul, d_alpha, d_beta, ga, gb, dims, offs, pair_i, pair_j, &
      !$omp            npair, tol, opt, n, block_max) &
      !$omp    private(ij, kl, s1, s2, s3, s4, d1, d2, d3, d4, o1, o2, o3, o4, &
      !$omp            shls, f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, deg, value, scaled, &
      !$omp            buf, ga_local, gb_local) &
      !$omp    reduction(+:n_total, n_computed, n_screened)
      allocate (buf(block_max**4))
      allocate (ga_local(n, n), gb_local(n, n))
      ga_local = 0.0_dp
      gb_local = 0.0_dp

      !$omp do schedule(dynamic)
      do ij = 1, npair
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

            if (bounds(s1, s2)*bounds(s3, s4) < tol) then
               n_screened = n_screened + 1_int64
               cycle
            end if

            shls = [s1 - 1, s2 - 1, s3 - 1, s4 - 1]
            ret = two_electron_block(mol%cartesian, buf, shls, mol%atm, mol%natm, &
                                     mol%bas, mol%nbas, mol%env, opt)
            if (ret == 0) then
               n_screened = n_screened + 1_int64
               cycle
            end if
            n_computed = n_computed + 1_int64

            ! Shell equality, not function equality: a block with s1 == s2
            ! already contains both orderings of its function pair.
            deg = 1.0_dp
            if (s1 /= s2) deg = deg*2.0_dp
            if (s3 /= s4) deg = deg*2.0_dp
            if (.not. (s1 == s3 .and. s2 == s4)) deg = deg*2.0_dp

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
                        ga_local(b1, b3) = ga_local(b1, b3) - 0.25_dp*d_alpha(b2, b4)*scaled
                        ga_local(b2, b4) = ga_local(b2, b4) - 0.25_dp*d_alpha(b1, b3)*scaled
                        ga_local(b1, b4) = ga_local(b1, b4) - 0.25_dp*d_alpha(b2, b3)*scaled
                        ga_local(b2, b3) = ga_local(b2, b3) - 0.25_dp*d_alpha(b1, b4)*scaled

                        gb_local(b1, b2) = gb_local(b1, b2) + d_coul(b3, b4)*scaled
                        gb_local(b3, b4) = gb_local(b3, b4) + d_coul(b1, b2)*scaled
                        gb_local(b1, b3) = gb_local(b1, b3) - 0.25_dp*d_beta(b2, b4)*scaled
                        gb_local(b2, b4) = gb_local(b2, b4) - 0.25_dp*d_beta(b1, b3)*scaled
                        gb_local(b1, b4) = gb_local(b1, b4) - 0.25_dp*d_beta(b2, b3)*scaled
                        gb_local(b2, b3) = gb_local(b2, b3) - 0.25_dp*d_beta(b1, b4)*scaled
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

      call libcint_del_optimizer(opt)

      fock_a = h + 0.5_dp*(ga + transpose(ga))
      fock_b = h + 0.5_dp*(gb + transpose(gb))

      deallocate (ga, gb, d_coul, dims, offs, pair_i, pair_j)
   end subroutine build_fock_direct_uhf

end module mqc_libcint_direct
