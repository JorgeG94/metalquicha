module mqc_dft_grid
   !! Molecular integration grid for DFT
   !!
   !! Assembles the three pieces into one quadrature. Each atom gets a spherical
   !! product grid -- a Treutler-Ahlrichs radial mesh times a Lebedev sphere --
   !! and the Becke partition then decides how much of each point belongs to the
   !! atom that produced it, so the atomic grids together cover space once.
   !!
   !! The weight of a point is
   !!
   !!    w = 4*pi * r^2 * dr * w_lebedev * w_becke
   !!
   !! where w_lebedev sums to 1 over the sphere. The 4*pi therefore appears
   !! exactly once, here; `mqc_dft_radial` deliberately leaves it out of its
   !! weights so it cannot be applied twice.
   !!
   !! Grid size follows a `level` from 0 to 9, with radial and angular counts
   !! chosen per element period. Level 3 is the default.
   use pic_types, only: dp, int64, int_index
   use pic_sorting, only: sort_index
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_lebedev, only: lebedev_grid
   use mqc_dft_radial, only: treutler_ahlrichs_radial
   use mqc_dft_partition, only: becke_partition_weights, PARTITION_BECKE, ADJUST_TREUTLER
   use mqc_dft_prune, only: prune_angular_orders, PRUNE_NONE, PRUNE_NWCHEM
   use mqc_physical_constants, only: PI
   implicit none
   private

   public :: dft_grid_t
   public :: build_dft_grid
   public :: grid_level_radial, grid_level_angular
   public :: MIN_GRID_LEVEL, MAX_GRID_LEVEL, DEFAULT_GRID_LEVEL
   public :: DEFAULT_PRUNE

   integer, parameter :: N_DIM = 3
   integer, parameter :: N_PERIODS = 7
   integer, parameter :: MIN_GRID_LEVEL = 0
   integer, parameter :: MAX_GRID_LEVEL = 9
   integer, parameter :: DEFAULT_GRID_LEVEL = 3

   integer, parameter :: DEFAULT_PRUNE = PRUNE_NWCHEM
      !! Angular pruning is on unless a caller asks for `PRUNE_NONE`.

   integer, parameter :: PERIOD_LAST_Z(N_PERIODS) = [2, 10, 18, 36, 54, 86, 118]
      !! Highest Z in each period, used to place an element in the level tables

   integer, parameter :: RAD_GRIDS(MIN_GRID_LEVEL:MAX_GRID_LEVEL, N_PERIODS) = &
                         transpose(reshape([ &
                                           10, 15, 20, 30, 35, 40, 50, &
                                           30, 40, 50, 60, 65, 70, 75, &
                                           40, 60, 65, 75, 80, 85, 90, &
                                           50, 75, 80, 90, 95, 100, 105, &
                                           60, 90, 95, 105, 110, 115, 120, &
                                           70, 105, 110, 120, 125, 130, 135, &
                                           80, 120, 125, 135, 140, 145, 150, &
                                           90, 135, 140, 150, 155, 160, 165, &
                                           100, 150, 155, 165, 170, 175, 180, &
                                           200, 200, 200, 200, 200, 200, 200], [N_PERIODS, MAX_GRID_LEVEL + 1]))
      !! Radial shells per (level, period)

   integer, parameter :: ANG_POINTS(MIN_GRID_LEVEL:MAX_GRID_LEVEL, N_PERIODS) = &
                         transpose(reshape([ &
                                           50, 86, 110, 110, 110, 110, 110, &
                                           110, 194, 194, 194, 194, 194, 194, &
                                           194, 302, 302, 302, 302, 302, 302, &
                                           302, 302, 434, 434, 434, 434, 434, &
                                           434, 590, 590, 590, 590, 590, 590, &
                                           590, 770, 770, 770, 770, 770, 770, &
                                           770, 974, 974, 974, 974, 974, 974, &
                                           974, 1202, 1202, 1202, 1202, 1202, 1202, &
                                           1202, 1202, 1202, 1202, 1202, 1202, 1202, &
                                           1454, 1454, 1454, 1454, 1454, 1454, 1454], &
                                           [N_PERIODS, MAX_GRID_LEVEL + 1]))
      !! Lebedev points per (level, period), as point counts rather than the
      !! degrees the standard tables list, since that is what `lebedev_grid`
      !! takes.

   real(dp), parameter :: GRID_SORT_CELL = 3.0_dp
      !! Side of the cubes the points are binned into before blocking, Bohr.
      !! About the spacing of bonded heavy atoms: a block of a few hundred
      !! points then sits within one or two cubes.

   type :: dft_grid_t
      !! A molecular quadrature: points, weights, and which atom made each point
      integer :: n_points = 0
      real(dp), allocatable :: coords(:, :)  !! (3, n_points), Bohr
      real(dp), allocatable :: weights(:)    !! (n_points), full volume weight
      real(dp), allocatable :: quad_weights(:)
         !! (n_points), the radial times angular weight *before* the Becke
         !! partition multiplies it in. Only the partition carries a nuclear
         !! derivative, so a gradient needs the two factors apart, and
         !! recovering this by dividing `weights` by the partition would divide
         !! by zero exactly where the partition underflows.
      integer, allocatable :: atom(:)        !! (n_points), owning atom
      ! What the partition was built with: a gradient has to differentiate the
      ! same partition the energy integrated over.
      integer :: scheme = PARTITION_BECKE
      integer :: adjust = ADJUST_TREUTLER
      integer, allocatable :: numbers(:)     !! (n_atoms), Z, for the partition radii
   contains
      procedure :: destroy => dft_grid_destroy
   end type dft_grid_t

contains

   pure function element_period(atomic_number) result(period)
      !! Period of an element, 1 to 7, for indexing the level tables
      integer, intent(in) :: atomic_number
      integer :: period
      integer :: i

      period = N_PERIODS
      do i = 1, N_PERIODS
         if (atomic_number <= PERIOD_LAST_Z(i)) then
            period = i
            exit
         end if
      end do
   end function element_period

   pure function grid_level_radial(atomic_number, level) result(n_radial)
      !! Radial shells this element gets at this level
      integer, intent(in) :: atomic_number, level
      integer :: n_radial
      n_radial = RAD_GRIDS(clamp_level(level), element_period(atomic_number))
   end function grid_level_radial

   pure function grid_level_angular(atomic_number, level) result(n_angular)
      !! Lebedev points this element gets at this level
      integer, intent(in) :: atomic_number, level
      integer :: n_angular
      n_angular = ANG_POINTS(clamp_level(level), element_period(atomic_number))
   end function grid_level_angular

   pure function clamp_level(level) result(clamped)
      integer, intent(in) :: level
      integer :: clamped
      clamped = max(MIN_GRID_LEVEL, min(MAX_GRID_LEVEL, level))
   end function clamp_level

   subroutine build_dft_grid(atom_coords, atomic_numbers, grid, error, &
                             level, scheme, adjust, n_radial, n_angular, prune)
      !! Build the molecular grid
      !!
      !! `level` picks per-element sizes from the standard tables. `n_radial`
      !! and `n_angular` override it for every atom; supplying one without the
      !! other is refused rather than silently half-applied.
      ! TODO(mqc): the override is still pruned. With the default
      ! `PRUNE_NWCHEM` an explicit `n_angular` is the valence-zone order only,
      ! so a convergence study asking for 302 points gets 50/86/266/302/266 by
      ! zone unless it also passes `prune = PRUNE_NONE`.
      real(dp), intent(in) :: atom_coords(:, :)   !! (3, n_atoms), Bohr
      integer, intent(in) :: atomic_numbers(:)    !! Z per atom
      type(dft_grid_t), intent(out) :: grid
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: level      !! 0 to 9, default 3
      integer, intent(in), optional :: scheme     !! Partition cutoff
      integer, intent(in), optional :: adjust     !! Size adjustment
      integer, intent(in), optional :: n_radial   !! Override, all atoms
      integer, intent(in), optional :: n_angular  !! Override, all atoms
      integer, intent(in), optional :: prune      !! PRUNE_NONE or PRUNE_NWCHEM

      real(dp), allocatable :: r(:), dr(:), sphere(:, :), w_ang(:)
      integer, allocatable :: shell_order(:)
      integer :: n_atoms, used_level, used_scheme, used_adjust, used_prune
      integer :: ia, i, j, k, nr, na, total, offset, cached_order

      n_atoms = size(atom_coords, 2)
      if (size(atomic_numbers) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "grid: atomic_numbers does not match atom_coords")
         return
      end if
      if (n_atoms < 1) then
         call error%set(ERROR_VALIDATION, "grid: no atoms")
         return
      end if

      used_level = DEFAULT_GRID_LEVEL
      if (present(level)) used_level = level
      used_scheme = PARTITION_BECKE
      if (present(scheme)) used_scheme = scheme
      used_adjust = ADJUST_TREUTLER
      if (present(adjust)) used_adjust = adjust
      used_prune = DEFAULT_PRUNE
      if (present(prune)) used_prune = prune

      if (present(n_radial) .neqv. present(n_angular)) then
         call error%set(ERROR_VALIDATION, &
                        "grid: n_radial and n_angular must be given together")
         return
      end if

      ! Size the whole grid first so it can be allocated once. With pruning the
      ! count is no longer nr*na: each shell carries its own order, and those
      ! depend on the radii, so the radial mesh is built here as well.
      total = 0
      do ia = 1, n_atoms
         call atom_sizes(atomic_numbers(ia), used_level, n_radial, n_angular, nr, na)
         call treutler_ahlrichs_radial(nr, atomic_numbers(ia), r, dr, error)
         if (error%has_error()) return
         allocate (shell_order(nr))
         call prune_angular_orders(used_prune, atomic_numbers(ia), r, na, shell_order, error)
         if (error%has_error()) return
         total = total + sum(shell_order)
         deallocate (r, dr, shell_order)
      end do

      grid%n_points = total
      allocate (grid%coords(N_DIM, total))
      allocate (grid%weights(total))
      allocate (grid%atom(total))

      offset = 0
      do ia = 1, n_atoms
         call atom_sizes(atomic_numbers(ia), used_level, n_radial, n_angular, nr, na)

         call treutler_ahlrichs_radial(nr, atomic_numbers(ia), r, dr, error)
         if (error%has_error()) return

         allocate (shell_order(nr))
         call prune_angular_orders(used_prune, atomic_numbers(ia), r, na, shell_order, error)
         if (error%has_error()) return

         ! Shells sharing an order are not contiguous -- the outer zone repeats
         ! the order of the third -- so hold the last sphere built and reuse it
         ! when the order does not change.
         cached_order = 0
         do i = 1, nr
            if (shell_order(i) /= cached_order) then
               if (allocated(sphere)) deallocate (sphere, w_ang)
               call lebedev_grid(shell_order(i), sphere, w_ang, error)
               if (error%has_error()) return
               cached_order = shell_order(i)
            end if

            do j = 1, shell_order(i)
               k = offset + j
               grid%coords(:, k) = atom_coords(:, ia) + r(i)*sphere(:, j)
               ! 4*pi appears once, here: w_ang sums to 1 and dr excludes it.
               grid%weights(k) = 4.0_dp*PI*r(i)*r(i)*dr(i)*w_ang(j)
               grid%atom(k) = ia
            end do
            offset = offset + shell_order(i)
         end do

         deallocate (r, dr, shell_order)
         if (allocated(sphere)) deallocate (sphere, w_ang)
      end do

      grid%scheme = used_scheme
      grid%adjust = used_adjust
      allocate (grid%numbers(n_atoms), source=atomic_numbers)

      call sort_grid_spatially(grid)
      call apply_partition(grid, atom_coords, atomic_numbers, used_scheme, used_adjust, error)
   end subroutine build_dft_grid

   subroutine sort_grid_spatially(grid)
      !! Reorder the points so that consecutive ones are close in space
      !!
      !! The loop above lays the points down atom by atom and radial shell by
      !! radial shell, so a block of consecutive points is a whole sphere around
      !! one atom -- and for the outer shells a sphere of ten Bohr or more,
      !! whose bounding sphere spans the molecule. Every basis function is
      !! significant on such a block and its screening does nothing. Binned
      !! into cubes of `GRID_SORT_CELL` Bohr and ordered along a Morton curve,
      !! consecutive points share a cube or a neighbouring one, and a block
      !! reaches only the functions near it. The quadrature is a sum, so the
      !! order changes nothing but roundoff.
      type(dft_grid_t), intent(inout) :: grid

      integer(int64), allocatable :: key(:)
      integer(int_index), allocatable :: perm(:)
      real(dp) :: lo(N_DIM)
      integer :: k, ix(N_DIM), n
      ! TEMP(mqc): switch for the before/after measurement, remove before committing.
      character(len=16) :: envbuf
      integer :: envlen, envstat

      call get_environment_variable("MQC_GRID_NOSORT", envbuf, envlen, envstat)
      if (envstat == 0 .and. envlen > 0) return

      n = grid%n_points
      if (n < 2) return
      allocate (key(n), perm(n))
      lo = minval(grid%coords, dim=2)
      do k = 1, n
         ix = int((grid%coords(:, k) - lo)/GRID_SORT_CELL)
         ix = min(max(ix, 0), 1023)
         key(k) = morton_key(ix(1), ix(2), ix(3))
      end do
      call sort_index(key, perm)
      grid%coords = grid%coords(:, perm)
      grid%weights = grid%weights(perm)
      grid%atom = grid%atom(perm)
   end subroutine sort_grid_spatially

   pure function morton_key(ix, iy, iz) result(key)
      !! The three cell indices interleaved bit by bit, ten bits each
      integer, intent(in) :: ix, iy, iz
      integer(int64) :: key

      integer :: b

      key = 0_int64
      do b = 0, 9
         if (btest(ix, b)) key = ibset(key, 3*b)
         if (btest(iy, b)) key = ibset(key, 3*b + 1)
         if (btest(iz, b)) key = ibset(key, 3*b + 2)
      end do
   end function morton_key

   pure subroutine atom_sizes(atomic_number, level, n_radial, n_angular, nr, na)
      !! Radial and angular counts for one atom, honouring any override
      integer, intent(in) :: atomic_number, level
      integer, intent(in), optional :: n_radial, n_angular
      integer, intent(out) :: nr, na

      if (present(n_radial)) then
         nr = n_radial
         na = n_angular
      else
         nr = grid_level_radial(atomic_number, level)
         na = grid_level_angular(atomic_number, level)
      end if
   end subroutine atom_sizes

   subroutine apply_partition(grid, atom_coords, atomic_numbers, scheme, adjust, error)
      !! Multiply the product weights by the Becke cell weights
      type(dft_grid_t), intent(inout) :: grid
      real(dp), intent(in) :: atom_coords(:, :)
      integer, intent(in) :: atomic_numbers(:)
      integer, intent(in) :: scheme, adjust
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: cell(:)

      allocate (cell(grid%n_points))
      call becke_partition_weights(grid%coords, atom_coords, atomic_numbers, &
                                   grid%atom, scheme, adjust, cell, error)
      if (error%has_error()) return

      ! Before, not after: this is the quadrature weight on its own.
      if (allocated(grid%quad_weights)) deallocate (grid%quad_weights)
      allocate (grid%quad_weights(grid%n_points), source=grid%weights)

      grid%weights = grid%weights*cell
   end subroutine apply_partition

   pure subroutine dft_grid_destroy(this)
      !! Release the grid
      class(dft_grid_t), intent(inout) :: this

      this%n_points = 0
      if (allocated(this%coords)) deallocate (this%coords)
      if (allocated(this%weights)) deallocate (this%weights)
      if (allocated(this%quad_weights)) deallocate (this%quad_weights)
      if (allocated(this%numbers)) deallocate (this%numbers)
      if (allocated(this%atom)) deallocate (this%atom)
   end subroutine dft_grid_destroy

end module mqc_dft_grid
