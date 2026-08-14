!! Charge-penetration screening parameters, by fitting the electrostatic potential
module mqc_libcint_screening
   !! The two screening blocks a potential carries are *fitted*, not computed: the
   !! damped classical multipole potential is matched to the quantum one on a grid
   !! of fused scaled van der Waals spheres, and the damping exponents are what
   !! comes out. So this module is a grid, an objective and a minimizer, and every
   !! constant in it was read out of GAMESS rather than guessed -- the exponents in
   !! a reference potential are the check.
   !!
   !! **Two blocks, two damping functions.** GAMESS names them:
   !!
   !!     SCREEN2   1 - exp(-alpha r)     exponential, EFP2 fragment-fragment
   !!     SCREEN    1 - exp(-alpha r^2)   Gaussian, EFP1 ab initio-fragment
   !!
   !! and the linear coefficient in front is frozen at one in both, because
   !! `1 - A exp(...)` only vanishes at the origin when `A = 1`, which is why GAMESS freezes it there.
   !! That is why the first column of every screening record is `1.000000000`.
   !!
   !! **The grid, from the source.** Radii are Gavezzotti and Spackman's, selected
   !! because the point scheme is geodesic, over Gavezzotti and Spackman's radii;
   !! spheres sit on *every* expansion centre, since GAMESS passes the full centre count to
   !! `PDCPTS`; a bond midpoint takes the mean of its two atoms' radii
   !! ; layers run from `VDWSCL = 0.7` in steps of `0.1`, and the
   !! points per layer grow as the square of the scale, which is a constant surface
   !! density and is also the only weighting there is.
   !!
   !! Getting that last point wrong is instructive: a fixed count per layer with an
   !! explicit `1/(layer+1)` weight -- weighting the inner layers up, the opposite
   !! way round -- leaves GAMESS's own fitted exponents scoring *worse* than the
   !! initial guess they started from, which a converged search cannot do.
   !!
   !! **`alpha = 10` is an absorbing state and the reference shows it.** At the upper
   !! bound the damping is switched off and the objective is flat in that exponent,
   !! so a search that wanders there gets nothing back. Water's reference has one
   !! bond midpoint at `10.0` and the other at `2.48`, and nothing distinguishes them
   !! physically. So a fitted set will not match a reference digit for digit -- part
   !! of what is in a reference is search history. What matches is the objective.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_atomic_radii, only: vdw_radius_geodesic
   use mqc_calculation_defaults, only: VDW_SCALE => DEFAULT_VDW_SCALE
   use mqc_physical_constants, only: ANGSTROM_TO_BOHR, HARTREE_TO_KCALMOL, PI
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_libcint_esp, only: esp_contract
   use mqc_libcint_dma, only: dma_result_t, N_QUAD
   implicit none
   private

   public :: fit_screening
   public :: screening_grid
   public :: screening_target_t
   public :: SCREEN_EXPONENTIAL, SCREEN_GAUSSIAN

   !> The part of a screening fit that does not depend on the damping function
   !>
   !> Both blocks are fitted to the *same* quantum potential on the *same* grid;
   !> only the objective's damping term differs. Computing them once and handing
   !> them to both fits halves the stage, and the expensive half is the one being
   !> shared -- the grid runs to tens of thousands of points and every one of them
   !> costs an integral over every shell pair.
   type :: screening_target_t
      logical :: ready = .false.
      real(dp), allocatable :: grid(:, :)     !! (3, n_points), Bohr
      real(dp), allocatable :: quantum(:)     !! Electronic potential at each point
   contains
      procedure :: destroy => screening_target_destroy
   end type screening_target_t

   !> Which damping function to fit.
   integer, parameter :: SCREEN_EXPONENTIAL = 1   !! `1 - exp(-alpha r)`, the SCREEN2 block
   integer, parameter :: SCREEN_GAUSSIAN = 2      !! `1 - exp(-alpha r^2)`, the SCREEN block

   !> Layer schedule: `VDWSCL`, `VDWINC` and the layer count.
   real(dp), parameter :: VDW_STEP = 0.1_dp
   integer, parameter :: N_LAYER = 25

   !> Points on the innermost sphere; outer layers scale by area from this.
   integer, parameter :: N_ANGULAR = 110

   !> Bounds on an exponent, and the starting values: two on an atom, four on a
   !> bond midpoint, which is how GAMESS initialises the search.
   real(dp), parameter :: ALPHA_MIN = 0.5_dp
   real(dp), parameter :: ALPHA_MAX = 10.0_dp
   real(dp), parameter :: ALPHA_ATOM = 2.0_dp
   real(dp), parameter :: ALPHA_MIDPOINT = 4.0_dp

   !> Sweeps of the minimizer, and the bracket it stops refining at.
   integer, parameter :: MAX_SWEEPS = 12
   real(dp), parameter :: ALPHA_TOL = 1.0e-7_dp

contains

   subroutine screening_grid(dma, atomic_numbers, grid, error)
      !! Points on fused scaled spheres about every expansion centre
      type(dma_result_t), intent(in) :: dma
      integer, intent(in) :: atomic_numbers(:)
      real(dp), allocatable, intent(out) :: grid(:, :)   !! (3, n_points), Bohr
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: radii(:), unit(:, :), kept(:, :)
      real(dp) :: scale, distance
      integer :: n_centre, i, j, k, layer, n_point, total, n_unit
      logical :: inside

      n_centre = size(dma%points, 2)
      allocate (radii(n_centre))
      do i = 1, n_centre
         if (dma%labels(i) (1:2) == "BO") then
            ! A midpoint takes the mean of the two atoms its label names.
            read (dma%labels(i) (3:3), *) j
            read (dma%labels(i) (4:4), *) k
            radii(i) = 0.5_dp*(vdw_radius_geodesic(atomic_numbers(j)) &
                               + vdw_radius_geodesic(atomic_numbers(k)))
         else
            radii(i) = vdw_radius_geodesic(atomic_numbers(i))
         end if
      end do

      ! Two passes: count, then fill, so the grid is exactly sized.
      total = 0
      do layer = 1, N_LAYER
         scale = VDW_SCALE + VDW_STEP*real(layer - 1, dp)
         n_unit = int((scale/VDW_SCALE)**2*real(N_ANGULAR, dp))
         total = total + n_centre*n_unit
      end do
      allocate (kept(3, total))

      n_point = 0
      do layer = 1, N_LAYER
         scale = VDW_SCALE + VDW_STEP*real(layer - 1, dp)
         n_unit = int((scale/VDW_SCALE)**2*real(N_ANGULAR, dp))
         call fibonacci_sphere(n_unit, unit)
         do i = 1, n_centre
            do j = 1, n_unit
               inside = .false.
               do k = 1, n_centre
                  if (k == i) cycle
                  distance = norm2(dma%points(:, i) &
                                   + radii(i)*ANGSTROM_TO_BOHR*scale*unit(:, j) &
                                   - dma%points(:, k))
                  if (distance <= radii(k)*ANGSTROM_TO_BOHR*scale) then
                     inside = .true.
                     exit
                  end if
               end do
               if (inside) cycle
               n_point = n_point + 1
               kept(:, n_point) = dma%points(:, i) &
                                  + radii(i)*ANGSTROM_TO_BOHR*scale*unit(:, j)
            end do
         end do
         deallocate (unit)
      end do

      if (n_point < 1) then
         call error%set(ERROR_VALIDATION, "screening: every grid point fell inside "// &
                        "another sphere, which cannot happen for a real molecule")
         return
      end if
      allocate (grid(3, n_point))
      grid = kept(:, 1:n_point)
      deallocate (kept, radii)
   end subroutine screening_grid

   subroutine fibonacci_sphere(n, points)
      !! `n` roughly equal-area points on the unit sphere
      !!
      !! Not the geodesic tessellation GAMESS builds, and it does not need to be:
      !! the objective is an average over thousands of points, and what matters is
      !! that the density is uniform and scales with area. The fitted exponents move
      !! in the last digit or two between the two constructions.
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: points(:, :)
      real(dp) :: polar, azimuth, offset
      integer :: i
      allocate (points(3, n))
      do i = 1, n
         offset = real(i, dp) - 0.5_dp
         polar = acos(1.0_dp - 2.0_dp*offset/real(n, dp))
         azimuth = PI*(1.0_dp + sqrt(5.0_dp))*offset
         points(1, i) = cos(azimuth)*sin(polar)
         points(2, i) = sin(azimuth)*sin(polar)
         points(3, i) = cos(polar)
      end do
   end subroutine fibonacci_sphere

   subroutine precompute(dma, grid, quantum, base, monopole, argument, kind)
      !! Split the objective into the part that depends on alpha and the part that does not
      !!
      !! The exponents enter only through the damped monopole, so the dipole and
      !! quadrupole contributions and the quantum target can be folded into one
      !! array per grid point once, and each objective evaluation then costs one
      !! exponential per point per centre instead of a multipole sum. Without this
      !! the fit is minutes rather than seconds, since a line search calls it
      !! thousands of times.
      type(dma_result_t), intent(in) :: dma
      real(dp), intent(in) :: grid(:, :), quantum(:)
      real(dp), allocatable, intent(out) :: base(:), monopole(:, :), argument(:, :)
      integer, intent(in) :: kind

      integer, parameter :: PACK_A(N_QUAD) = [1, 2, 3, 1, 1, 2]
      integer, parameter :: PACK_B(N_QUAD) = [1, 2, 3, 2, 3, 3]
      real(dp) :: d(3), q(3, 3)
      real(dp) :: r
      integer :: g, k, c, a, b, n_grid, n_centre

      n_grid = size(grid, 2)
      n_centre = size(dma%points, 2)
      allocate (base(n_grid), monopole(n_grid, n_centre), argument(n_grid, n_centre))

      !$omp parallel do default(shared) private(g, k, c, a, b, d, r, q)
      do g = 1, n_grid
         base(g) = -quantum(g)
         do k = 1, n_centre
            d = grid(:, g) - dma%points(:, k)
            r = norm2(d)
            monopole(g, k) = dma%electronic(k)/r
            ! r for the exponential form, r^2 for the Gaussian one.
            argument(g, k) = r
            if (kind == SCREEN_GAUSSIAN) argument(g, k) = r*r
            base(g) = base(g) + dot_product(dma%dipole(:, k), d)/r**3
            q = 0.0_dp
            do c = 1, N_QUAD
               a = PACK_A(c)
               b = PACK_B(c)
               q(a, b) = dma%quadrupole(c, k)
               q(b, a) = dma%quadrupole(c, k)
            end do
            base(g) = base(g) + 1.5_dp*dot_product(d, matmul(q, d))/r**5 &
                      - 0.5_dp*(q(1, 1) + q(2, 2) + q(3, 3))/r**3
         end do
      end do
      !$omp end parallel do
   end subroutine precompute

   function objective(base, monopole, argument, alpha) result(rms)
      !! Root-mean-square miss of the damped classical potential, kcal/mol
      real(dp), intent(in) :: base(:), monopole(:, :), argument(:, :), alpha(:)
      real(dp) :: rms
      real(dp) :: total, residual
      integer :: g, k

      rms = 0.0_dp
      !$omp parallel do default(shared) private(g, k, total, residual) reduction(+:rms)
      do g = 1, size(base)
         total = base(g)
         do k = 1, size(alpha)
            total = total + monopole(g, k)*(1.0_dp - exp(-alpha(k)*argument(g, k)))
         end do
         residual = total*HARTREE_TO_KCALMOL
         rms = rms + residual*residual
      end do
      !$omp end parallel do
      rms = sqrt(rms/real(size(base), dp))
   end function objective

   subroutine hold_others(base, monopole, argument, alpha, skip, rest)
      !! The damped potential with one centre's contribution left out
      !!
      !! A line search moves one exponent and leaves the rest where they are, so
      !! everything but that centre is constant across the whole bracket. Summed
      !! once here, the search then costs one exponential per grid point instead of
      !! one per point per centre -- which on eighteen centres is the difference
      !! between the fit being the longest stage of a potential and being noise.
      !!
      !! Summed rather than subtracted from a running total. Subtracting would be
      !! one operation instead of eighteen, but the total is a residual that the
      !! fit is driving toward zero, so the cancellation gets worse exactly as the
      !! search converges.
      real(dp), intent(in) :: base(:), monopole(:, :), argument(:, :), alpha(:)
      integer, intent(in) :: skip
      real(dp), intent(out) :: rest(:)

      real(dp) :: total
      integer :: g, k

      !$omp parallel do default(shared) private(g, k, total)
      do g = 1, size(base)
         total = base(g)
         do k = 1, size(alpha)
            if (k == skip) cycle
            total = total + monopole(g, k)*(1.0_dp - exp(-alpha(k)*argument(g, k)))
         end do
         rest(g) = total
      end do
      !$omp end parallel do
   end subroutine hold_others

   function objective_one(rest, monopole, argument, alpha) result(rms)
      !! The objective with every other centre already summed into `rest`
      !!
      !! Identical to `objective` term for term; what differs is that seventeen of
      !! the eighteen exponentials were evaluated once by `hold_others` rather than
      !! again for every probe of the bracket.
      real(dp), intent(in) :: rest(:), monopole(:), argument(:)
      real(dp), intent(in) :: alpha
      real(dp) :: rms
      real(dp) :: residual
      integer :: g

      rms = 0.0_dp
      !$omp parallel do default(shared) private(g, residual) reduction(+:rms)
      do g = 1, size(rest)
         residual = (rest(g) + monopole(g)*(1.0_dp - exp(-alpha*argument(g)))) &
                    *KCAL_PER_HARTREE
         rms = rms + residual*residual
      end do
      !$omp end parallel do
      rms = sqrt(rms/real(size(rest), dp))
   end function objective_one

   subroutine screening_target_destroy(this)
      !! Release the shared grid and potential
      class(screening_target_t), intent(inout) :: this

      if (allocated(this%grid)) deallocate (this%grid)
      if (allocated(this%quantum)) deallocate (this%quantum)
      this%ready = .false.
   end subroutine screening_target_destroy

   subroutine fit_screening(mol, density, dma, atomic_numbers, kind, alpha, error, &
                            residual, grid_size, target)
      !! Fit one damping exponent per expansion centre
      !!
      !! Minimized by repeated bracketed line searches over one exponent at a time,
      !! swept to convergence. GAMESS uses Powell's direction-set method; this is the
      !! same objective and the same bounds, and on five smooth parameters both land
      !! on the same minimum -- which is what is being reproduced, since the search
      !! *path* is exactly the part of a reference that is not physics.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      type(dma_result_t), intent(in) :: dma
      integer, intent(in) :: atomic_numbers(:)
      integer, intent(in) :: kind
      real(dp), allocatable, intent(out) :: alpha(:)
      type(error_t), intent(inout) :: error
      real(dp), intent(out), optional :: residual
      integer, intent(out), optional :: grid_size
      type(screening_target_t), intent(inout), optional :: target
         !! The grid and quantum potential, built here if it is not ready and kept
         !! for the next fit. Both damping forms are fitted to the same target, so
         !! passing this across the pair halves the stage. Absent means build it,
         !! use it and drop it, which is what a lone fit wants.

      real(dp), allocatable :: grid(:, :), potential(:, :, :), quantum(:)
      real(dp), allocatable :: base(:), monopole(:, :), argument(:, :), rest(:)
      real(dp) :: best, trial, low, high, mid, value_low, value_high, keep
      integer :: n_centre, i, sweep, step, g
      logical :: shared

      if (kind /= SCREEN_EXPONENTIAL .and. kind /= SCREEN_GAUSSIAN) then
         call error%set(ERROR_VALIDATION, "screening: unknown damping function")
         return
      end if

      shared = .false.
      if (present(target)) shared = target%ready

      if (shared) then
         grid = target%grid
         quantum = target%quantum
      else
         call screening_grid(dma, atomic_numbers, grid, error)
         if (error%has_error()) return

         ! The quantum target: the electronic potential, nuclei excluded, because the
         ! classical side above carries the electronic monopoles alone.
         ! Contracted inside the integral loop. Holding the whole
         ! `(n_ao, n_ao, n_grid)` tensor to form this one vector is 786 MB at 58
         ! orbitals and 3.1 GB at 115, for a grid of about thirty thousand points.
         call esp_contract(mol, grid, density, quantum, error)
         if (error%has_error()) return

         if (present(target)) then
            target%grid = grid
            target%quantum = quantum
            target%ready = .true.
         end if
      end if

      n_centre = size(dma%points, 2)
      allocate (alpha(n_centre))
      do i = 1, n_centre
         if (dma%labels(i) (1:2) == "BO") then
            alpha(i) = ALPHA_MIDPOINT
         else
            alpha(i) = ALPHA_ATOM
         end if
      end do

      call precompute(dma, grid, quantum, base, monopole, argument, kind)
      allocate (rest(size(base)))
      best = objective(base, monopole, argument, alpha)
      do sweep = 1, MAX_SWEEPS
         keep = best
         do i = 1, n_centre
            ! Everything but centre `i` is fixed for the whole bracket below, so it
            ! is summed once here and the sixty-one probes that follow each cost a
            ! single exponential per grid point.
            call hold_others(base, monopole, argument, alpha, i, rest)

            ! Golden-section on one exponent, inside the bounds.
            low = ALPHA_MIN
            high = ALPHA_MAX
            do step = 1, 30
               if (high - low < ALPHA_TOL) exit
               mid = low + 0.381966_dp*(high - low)
               trial = high - 0.381966_dp*(high - low)
               value_low = objective_one(rest, monopole(:, i), argument(:, i), mid)
               value_high = objective_one(rest, monopole(:, i), argument(:, i), trial)
               if (value_low < value_high) then
                  high = trial
               else
                  low = mid
               end if
            end do
            mid = 0.5_dp*(low + high)
            trial = objective_one(rest, monopole(:, i), argument(:, i), mid)
            if (trial < best) then
               best = trial
               alpha(i) = mid
            end if
         end do
         if (keep - best < 1.0e-10_dp) exit
      end do

      if (present(residual)) residual = best
      if (present(grid_size)) grid_size = size(grid, 2)
      deallocate (grid, quantum, base, monopole, argument, rest)
   end subroutine fit_screening

end module mqc_libcint_screening
