!! Electrostatic potential integrals at arbitrary points
module mqc_libcint_esp
   !! `<chi_u| 1/|r - R_g| |chi_v>` for a list of points `R_g`.
   !!
   !! Contracted against a density this gives the electronic electrostatic
   !! potential the molecule produces at each point, which is what a charge-
   !! penetration screening fit is fitted *to*: GAMESS optimizes its damping
   !! exponents by matching a damped classical multipole potential to the quantum
   !! one on a grid, so the quantum potential is the target and this is how it is
   !! computed. The same integral is what a CPU polarizable-continuum model would
   !! need; only the GPU path has one today.
   !!
   !! **libcint's grid integral is shaped differently from its other one-electron
   !! integrals, in three ways that all have to be right at once.**
   !!
   !!   * The points live in `env`, at the offset stored in `env(PTR_GRIDS)` --
   !!     an offset into `env` itself, not a coordinate. `PTR_GRIDS` is 12 and, like
   !!     `PTR_COMMON_ORIG` and `PTR_RANGE_OMEGA`, is one of libcint's own 0-based
   !!     slot constants that the Fortran interface does *not* convert, so it
   !!     carries a `+ 1` here. It is also not among the constants that interface
   !!     exports, hence the local parameter.
   !!   * The grid range is passed in `shls`, which is four long rather than two:
   !!     `[ish, jsh, first_grid, last_grid]`, and libcint takes
   !!     `ngrids = shls(4) - shls(3)`. So the shell pair and the grid block travel
   !!     in the same array.
   !!   * The signature has ten arguments, not seven -- `dims`, an optimizer and a
   !!     cache join the usual list. All three are passed null: null `dims` means
   !!     the natural block shape, there is no optimizer for this integral
   !!     (`int1e_grids_optimizer` sets it to null itself), and a null cache makes
   !!     libcint allocate and free its own.
   use pic_types, only: dp
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_loc, c_null_ptr
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim
   implicit none
   private

   public :: esp_matrices

   !> `PTR_GRIDS` from libcint's `cint.h`. Not exported by the Fortran interface,
   !> and 0-based like every other `PTR_*`, so it is used as `+ 1`.
   integer, parameter :: LIBCINT_PTR_GRIDS = 12

   interface
      function cint1e_grids_sph(buf, dims, shls, atm, natm, bas, nbas, env, opt, &
                                cache) result(ret) bind(C, name="int1e_grids_sph")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value, intent(in) :: buf, dims, shls, atm
         integer(c_int), value, intent(in) :: natm
         type(c_ptr), value, intent(in) :: bas
         integer(c_int), value, intent(in) :: nbas
         type(c_ptr), value, intent(in) :: env, opt, cache
         integer(c_int) :: ret
      end function cint1e_grids_sph

      function cint1e_grids_cart(buf, dims, shls, atm, natm, bas, nbas, env, opt, &
                                 cache) result(ret) bind(C, name="int1e_grids_cart")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value, intent(in) :: buf, dims, shls, atm
         integer(c_int), value, intent(in) :: natm
         type(c_ptr), value, intent(in) :: bas
         integer(c_int), value, intent(in) :: nbas
         type(c_ptr), value, intent(in) :: env, opt, cache
         integer(c_int) :: ret
      end function cint1e_grids_cart
   end interface

contains

   subroutine esp_matrices(mol, grids, matrices, error)
      !! `1/|r - R_g|` over every shell pair, for each point
      !!
      !! Returns `(n_ao, n_ao, n_points)`. Contract with a density and negate to get
      !! the electronic potential; the nuclear part is `sum_A Z_A/|R_g - R_A|` and is
      !! the caller's, since a screening fit wants the electronic part alone.
      type(libcint_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: grids(:, :)          !! (3, n_points), Bohr
      real(dp), allocatable, intent(out) :: matrices(:, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable, target :: buf(:), env_local(:)
      integer, target :: shls(4)
      integer :: ish, jsh, di, dj, i, j, g, ret, io, jo, block_max, n_grid
      integer :: grid_offset

      if (size(grids, 1) /= 3) then
         call error%set(ERROR_VALIDATION, "ESP integrals: grid points must be (3, n)")
         return
      end if
      n_grid = size(grids, 2)
      if (n_grid < 1) then
         call error%set(ERROR_VALIDATION, "ESP integrals: no grid points")
         return
      end if

      ! The points are appended to a copy of env, and the offset to them recorded in
      ! the PTR_GRIDS slot. libcint reads that slot as an index into env, so what
      ! goes there is the 0-based offset of the first coordinate.
      grid_offset = size(mol%env)
      allocate (env_local(grid_offset + 3*n_grid))
      env_local(1:grid_offset) = mol%env
      do g = 1, n_grid
         env_local(grid_offset + 3*(g - 1) + 1:grid_offset + 3*g) = grids(:, g)
      end do
      env_local(LIBCINT_PTR_GRIDS + 1) = real(grid_offset, dp)

      block_max = 1
      do ish = 1, mol%nbas
         block_max = max(block_max, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do

      allocate (matrices(mol%nao, mol%nao, n_grid))
      matrices = 0.0_dp
      allocate (buf(block_max*block_max*n_grid))

      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, mol%nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            ! All the points in one call: the shell pair and the grid range share
            ! this array, and 0 to n_grid is every point.
            shls = [ish - 1, jsh - 1, 0, n_grid]

            if (mol%cartesian) then
               ret = cint1e_grids_cart(c_loc(buf), c_null_ptr, c_loc(shls), &
                                       c_loc(mol%atm), mol%natm, c_loc(mol%bas), &
                                       mol%nbas, c_loc(env_local), c_null_ptr, &
                                       c_null_ptr)
            else
               ret = cint1e_grids_sph(c_loc(buf), c_null_ptr, c_loc(shls), &
                                      c_loc(mol%atm), mol%natm, c_loc(mol%bas), &
                                      mol%nbas, c_loc(env_local), c_null_ptr, &
                                      c_null_ptr)
            end if
            if (ret == 0) cycle      ! screened away, leave the block zero

            ! The grid index runs fastest, then the bra, then the ket -- the
            ! opposite nesting from the multipole blocks, where the component is
            ! slowest. Getting it backwards gives a matrix that is plausible at one
            ! point and wrong at every other, so the check against PySCF is
            ! elementwise over all points rather than on a contracted potential.
            do j = 1, dj
               do i = 1, di
                  do g = 1, n_grid
                     matrices(io + i, jo + j, g) = &
                        buf(g + (i - 1)*n_grid + (j - 1)*n_grid*di)
                  end do
               end do
            end do
         end do
      end do

      deallocate (buf, env_local)
   end subroutine esp_matrices

end module mqc_libcint_esp
