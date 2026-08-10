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

      real(dp), allocatable :: buf(:), g(:, :), d_half(:, :)
      type(c_ptr) :: opt
      integer :: s1, s2, s3, s4, s4_max
      integer :: d1, d2, d3, d4, o1, o2, o3, o4
      integer :: shls(4)
      integer :: f1, f2, f3, f4, b1, b2, b3, b4, idx, ret, block_max, n
      real(dp) :: tol, deg, value, scaled

      n = mol%nao
      if (size(h, 1) /= n .or. size(density, 1) /= n .or. size(fock, 1) /= n) then
         call error%set(ERROR_VALIDATION, "direct Fock: matrix dimensions do not match the basis")
         return
      end if

      tol = DEFAULT_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol

      block_max = 1
      do s1 = 1, mol%nbas
         block_max = max(block_max, mol%shell_offset(s1 + 1) - mol%shell_offset(s1))
      end do

      allocate (buf(block_max**4))
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

      do s1 = 1, mol%nbas
         d1 = shell_dim(mol%cartesian, s1 - 1, mol%bas)
         o1 = mol%shell_offset(s1)
         do s2 = 1, s1
            d2 = shell_dim(mol%cartesian, s2 - 1, mol%bas)
            o2 = mol%shell_offset(s2)
            do s3 = 1, s1
               d3 = shell_dim(mol%cartesian, s3 - 1, mol%bas)
               o3 = mol%shell_offset(s3)

               ! The pair (s3,s4) must not run past (s1,s2), or quartets would
               ! be counted twice under the bra-ket swap.
               s4_max = s3
               if (s3 == s1) s4_max = s2

               do s4 = 1, s4_max
                  d4 = shell_dim(mol%cartesian, s4 - 1, mol%bas)
                  o4 = mol%shell_offset(s4)

                  stats%quartets_total = stats%quartets_total + 1

                  if (bounds(s1, s2)*bounds(s3, s4) < tol) then
                     stats%quartets_screened = stats%quartets_screened + 1
                     cycle
                  end if

                  shls = [s1 - 1, s2 - 1, s3 - 1, s4 - 1]
                  ret = two_electron_block(mol%cartesian, buf, shls, mol%atm, mol%natm, &
                                           mol%bas, mol%nbas, mol%env, opt)
                  if (ret == 0) then
                     stats%quartets_screened = stats%quartets_screened + 1
                     cycle
                  end if
                  stats%quartets_computed = stats%quartets_computed + 1

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
                              g(b1, b2) = g(b1, b2) + d_half(b3, b4)*scaled
                              g(b3, b4) = g(b3, b4) + d_half(b1, b2)*scaled
                              g(b1, b3) = g(b1, b3) - 0.25_dp*d_half(b2, b4)*scaled
                              g(b2, b4) = g(b2, b4) - 0.25_dp*d_half(b1, b3)*scaled
                              g(b1, b4) = g(b1, b4) - 0.25_dp*d_half(b2, b3)*scaled
                              g(b2, b3) = g(b2, b3) - 0.25_dp*d_half(b1, b4)*scaled
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      call libcint_del_optimizer(opt)

      fock = h + 0.5_dp*(g + transpose(g))

      deallocate (buf, g, d_half)
   end subroutine build_fock_direct

end module mqc_libcint_direct
