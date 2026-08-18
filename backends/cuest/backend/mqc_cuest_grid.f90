!! Numerical integration grids for exchange-correlation
module mqc_cuest_grid
   !! Builds the per-atom quadrature grids a DFT calculation integrates the
   !! exchange-correlation energy over.
   !!
   !! The scheme is an unpruned direct product: a Treutler-Ahlrichs M4 radial
   !! quadrature scaled by an element-dependent radius, crossed with the same
   !! Lebedev angular order at every radial shell. cuEST combines the per-atom
   !! grids into a molecular grid, applying the atomic partition weights.
   !!
   !! Unpruned means every radial shell carries the full angular order, so the
   !! cost is `n_radial * n_angular` points per atom with no savings near the
   !! nucleus. It is the scheme cuEST's own samples use, and it keeps the grid
   !! reproducible; pruning would be the first place to look for speed.
   !!
   !! Reference: O. Treutler and R. Ahlrichs, J. Chem. Phys. 102, 346 (1995).
   !!
   !! Host/device note: radial nodes, radial weights and the angular-point
   !! counts are all CPU arrays, unlike most cuEST buffers.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_int64_t, &
                                                                             c_double, c_loc, c_associated
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use pic_io, only: to_char
   use mqc_cuest_runtime, only: cuest_status_check
   use cuest, only: cuestAtomGridCreate, cuestAtomGridDestroy, &
                    cuestParametersCreate, cuestParametersDestroy, &
                    CUEST_ATOMGRID_PARAMETERS
   implicit none
   private

   public :: atom_grid_set_t      !! Owned collection of cuEST atom grids
   public :: build_atom_grids     !! Build one atom grid per atom
   public :: ahlrichs_radius      !! Element radius used to scale the radial mesh
   public :: ahlrichs_radial_quadrature  !! Treutler-Ahlrichs M4 radial mesh

   real(dp), parameter :: PI = 3.14159265358979323846_dp

   real(dp), parameter :: M4_ALPHA = 0.6_dp
      !! Exponent of the (1+x) factor in the Treutler-Ahlrichs M4 mapping

   integer, parameter :: MAX_TABULATED_Z = 36  !! Radii are tabulated H..Kr

   real(dp), parameter :: AHLRICHS_RADII(0:MAX_TABULATED_Z) = [ &
                          1.00_dp, &
                          0.80_dp, 0.90_dp, &
                          1.80_dp, 1.40_dp, 1.30_dp, 1.10_dp, 0.90_dp, 0.90_dp, 0.90_dp, 0.90_dp, &
                          1.40_dp, 1.30_dp, 1.30_dp, 1.20_dp, 1.10_dp, 1.00_dp, 1.00_dp, 1.00_dp, &
                          1.50_dp, 1.40_dp, 1.30_dp, 1.20_dp, 1.20_dp, 1.20_dp, 1.20_dp, 1.20_dp, &
                          1.20_dp, 1.10_dp, 1.10_dp, 1.10_dp, 1.10_dp, 1.00_dp, 0.90_dp, 0.90_dp, &
                          0.90_dp, 0.90_dp]
      !! Ahlrichs radii in Bohr, indexed by atomic number (index 0 is a dummy)

   type :: atom_grid_set_t
      !! The per-atom cuEST grid handles for one molecule
      !!
      !! Handles are owned here and released by `destroy`. The set only needs
      !! to outlive the molecular grid built from it.
      integer :: n_atoms = 0
      integer :: n_radial = 0
      integer :: n_angular = 0
      type(c_ptr), allocatable :: grids(:)  !! cuestAtomGrid_t per atom
   contains
      procedure :: destroy => atom_grid_set_destroy
   end type atom_grid_set_t

contains

   pure function ahlrichs_radius(atomic_number) result(radius)
      !! Element radius used to scale the radial quadrature, in Bohr
      !!
      !! Falls back to 1.0 beyond krypton, which is what cuEST's own helper
      !! does; heavier elements therefore get an unoptimized radial mesh.
      integer, intent(in) :: atomic_number
      real(dp) :: radius

      if (atomic_number >= 1 .and. atomic_number <= MAX_TABULATED_Z) then
         radius = AHLRICHS_RADII(atomic_number)
      else
         radius = 1.0_dp
      end if
   end function ahlrichs_radius

   pure subroutine ahlrichs_radial_quadrature(n_points, radius, nodes, weights)
      !! Treutler-Ahlrichs M4 radial quadrature on `n_points` points
      !!
      !! Filled back to front so the nodes come out in ascending order.
      integer, intent(in) :: n_points
      real(dp), intent(in) :: radius       !! Element radius, Bohr
      real(dp), intent(out) :: nodes(:)    !! Radial nodes, ascending
      real(dp), intent(out) :: weights(:)  !! Quadrature weights, r^2 included

      real(dp) :: angle, cos_angle, sin_angle, log_term, map_term, r, w, n_plus_one
      integer :: i

      n_plus_one = real(n_points, dp) + 1.0_dp
      do i = 1, n_points
         angle = real(i, dp)*PI/n_plus_one
         cos_angle = cos(angle)
         sin_angle = sin(angle)
         log_term = log((1.0_dp - cos_angle)/2.0_dp)
         map_term = (1.0_dp + cos_angle)**M4_ALPHA/log(2.0_dp)
         r = -radius*map_term*log_term
         w = PI/n_plus_one*sin_angle*radius*map_term &
             *(-M4_ALPHA*log_term/(1.0_dp + cos_angle) + 1.0_dp/(1.0_dp - cos_angle))*r*r
         nodes(n_points - i + 1) = r
         weights(n_points - i + 1) = w
      end do
   end subroutine ahlrichs_radial_quadrature

   subroutine build_atom_grids(handle, atomic_numbers, n_radial, n_angular, grid_set, error, &
                               shard_index, shard_count)
      !! Build one unpruned direct-product grid per atom
      !!
      !! With `shard_count` > 1 each atom keeps only every `shard_count`-th
      !! radial shell, starting at `shard_index`. Every atom is still present
      !! and every atom still gets a grid, so the molecular grid built from the
      !! result covers the whole molecule at a lower radial density; the Becke
      !! partition weight at a surviving point is unchanged, because it is a
      !! function of that point and of all the nuclear positions and knows
      !! nothing about which other points exist. Summing the slices therefore
      !! reproduces the unsharded quadrature exactly, term for term -- this is
      !! a partition of the sum, not an approximation to it.
      !!
      !! Striding rather than blocking on purpose: consecutive radial shells
      !! carry wildly different weight (the mesh is dense near the nucleus and
      !! sparse in the tail), so a contiguous block would hand one rank the
      !! core and another the vacuum. Striding gives every rank the same mix.
      type(c_ptr), intent(in) :: handle          !! cuEST handle
      integer, intent(in) :: atomic_numbers(:)   !! Z per atom
      integer, intent(in) :: n_radial            !! Radial points per atom, before sharding
      integer, intent(in) :: n_angular           !! Lebedev order per radial shell
      type(atom_grid_set_t), intent(out) :: grid_set
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: shard_index  !! 0-based slice to keep
      integer, intent(in), optional :: shard_count  !! Number of slices; 1 keeps everything

      real(dp), allocatable :: nodes(:), weights(:)
      real(dp), allocatable :: my_nodes(:), my_weights(:)
      integer(c_int64_t), allocatable :: angular_points(:)
      type(c_ptr) :: grid_params
      integer(c_int) :: status
      integer :: iatom, n_atoms, n_shards, my_shard, n_kept, i, k

      n_shards = 1
      my_shard = 0
      if (present(shard_count)) n_shards = max(shard_count, 1)
      if (present(shard_index)) my_shard = shard_index

      ! A shard with no radial shells would leave this rank without an XC plan,
      ! and the XC plan is what the DF plan asks for its exact-exchange
      ! fraction -- so the rank would go on to build a *differently defined*
      ! Fock operator rather than merely contribute nothing. Refused here,
      ! where the cause is still visible.
      if (n_shards > n_radial) then
         call error%set(ERROR_VALIDATION, &
                        "multi-GPU XC grid: "//to_char(n_shards)//" ranks share "// &
                        to_char(n_radial)//" radial shells, so some rank would get none. "// &
                        "Raise keywords.dft.radial_points to at least the rank count, "// &
                        "or run on fewer ranks.")
         return
      end if

      n_kept = (n_radial - my_shard - 1)/n_shards + 1

      n_atoms = size(atomic_numbers)
      grid_set%n_atoms = n_atoms
      grid_set%n_radial = n_kept
      grid_set%n_angular = n_angular

      allocate (grid_set%grids(n_atoms))
      grid_set%grids = c_null_ptr

      grid_params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_ATOMGRID_PARAMETERS, grid_params), &
                              "cuestParametersCreate(atom grid)", error)
      if (error%has_error()) return

      ! One angular order for every radial shell -- this is what makes the
      ! grid "unpruned".
      allocate (nodes(n_radial), weights(n_radial))
      allocate (my_nodes(n_kept), my_weights(n_kept), angular_points(n_kept))
      angular_points = int(n_angular, c_int64_t)

      do iatom = 1, n_atoms
         call ahlrichs_radial_quadrature(n_radial, ahlrichs_radius(atomic_numbers(iatom)), &
                                         nodes, weights)
         k = 0
         do i = my_shard + 1, n_radial, n_shards
            k = k + 1
            my_nodes(k) = nodes(i)
            my_weights(k) = weights(i)
         end do
         call create_atom_grid(handle, int(n_kept, c_int64_t), my_nodes, my_weights, &
                               angular_points, grid_params, grid_set%grids(iatom), error)
         if (error%has_error()) exit
      end do

      status = cuestParametersDestroy(CUEST_ATOMGRID_PARAMETERS, grid_params)
      call cuest_status_check(status, "cuestParametersDestroy(atom grid)", error)

      deallocate (nodes, weights, my_nodes, my_weights, angular_points)
      if (error%has_error()) then
         call error%add_context("mqc_cuest_grid:build_atom_grids")
         call grid_set%destroy()
      end if
   end subroutine build_atom_grids

   subroutine create_atom_grid(handle, n_points, nodes, weights, angular_points, &
                               grid_params, grid, error)
      !! Create one cuEST atom grid from host quadrature arrays
      !!
      !! Split out so the node/weight arrays get the TARGET attribute C_LOC
      !! needs. All three arrays are HOST memory, not device buffers.
      type(c_ptr), intent(in) :: handle, grid_params
      integer(c_int64_t), intent(in) :: n_points
      real(dp), intent(in), target :: nodes(:), weights(:)
      integer(c_int64_t), intent(in) :: angular_points(:)
      type(c_ptr), intent(inout) :: grid
      type(error_t), intent(inout) :: error

      call cuest_status_check(cuestAtomGridCreate(handle, n_points, c_loc(nodes), &
                                                  c_loc(weights), angular_points, &
                                                  grid_params, grid), &
                              "cuestAtomGridCreate", error)
   end subroutine create_atom_grid

   subroutine atom_grid_set_destroy(this)
      !! Release every atom grid handle
      class(atom_grid_set_t), intent(inout) :: this
      integer(c_int) :: status
      integer :: i

      if (allocated(this%grids)) then
         do i = 1, size(this%grids)
            if (c_associated(this%grids(i))) status = cuestAtomGridDestroy(this%grids(i))
         end do
         deallocate (this%grids)
      end if
      this%n_atoms = 0
      this%n_radial = 0
      this%n_angular = 0
   end subroutine atom_grid_set_destroy

end module mqc_cuest_grid
