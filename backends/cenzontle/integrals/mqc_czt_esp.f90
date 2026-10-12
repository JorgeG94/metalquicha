!! Electrostatic potential integrals at arbitrary points
module mqc_czt_esp
   !! `<chi_u| 1/|r - R_g| |chi_v>` for a list of points `R_g`.
   !!
   !! Contracted against a density these give the electronic electrostatic
   !! potential the molecule produces at each point, which is what a
   !! charge-penetration screening fit is fitted to.
   !!
   !! **libcint's grid integral is shaped differently from its other one-electron
   !! integrals, in two ways that both have to be right at once.**
   !!
   !!   * The points live in `env`, at the offset stored in `env(PTR_GRIDS)` --
   !!     an offset into `env` itself, not a coordinate. `PTR_GRIDS` is 12 and,
   !!     like `PTR_COMMON_ORIG` and `PTR_RANGE_OMEGA`, is one of libcint's own
   !!     0-based slot constants that the Fortran interface does not convert, so it
   !!     carries a `+ 1`. `LIBCINT_PTR_GRIDS` comes from `mqc_czt_cint_compat`.
   !!   * The grid range is passed in `shls`, which is four long rather than two:
   !!     `[ish, jsh, first_grid, last_grid]`, and libcint takes
   !!     `ngrids = shls(4) - shls(3)`.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, shell_dim
   use libcint_fortran, only: LIBCINT_PTR_RINV_ORIG
   use mqc_czt_cint_compat, only: LIBCINT_PTR_GRIDS, &
                                  libcint_1e_grids_cart, libcint_1e_grids_sph, &
                                  libcint_1e_drinv_cart, libcint_1e_drinv_sph, &
                                  libcint_1e_ipiprinv_cart, libcint_1e_ipiprinv_sph, &
                                  libcint_1e_iprinvip_cart, libcint_1e_iprinvip_sph
   implicit none
   private

   public :: esp_matrices
   public :: esp_contract
   public :: drinv_matrices
   public :: ddrinv_matrices

contains

   subroutine esp_contract(mol, grids, density, values, error)
      !! The electronic potential a density produces at each point
      !!
      !! `values(g) = -sum_uv D_uv <chi_u| 1/|r - R_g| |chi_v>`, contracted inside the
      !! integral loop rather than after it.
      !!
      !! `esp_matrices` returns the whole `(n_ao, n_ao, n_grid)` tensor, which a
      !! charge-penetration grid of tens of thousands of points cannot hold;
      !! contracting as the blocks come out needs the grid vector and nothing
      !! else. Use the matrices when they are needed more than once.
      use omp_lib, only: omp_get_max_threads
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: grids(:, :)          !! (3, n_points), Bohr
      real(dp), intent(in) :: density(:, :)        !! AO density
      real(dp), allocatable, intent(out) :: values(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: env_local(:), buf(:)
      real(dp), allocatable :: local(:)
      integer, allocatable :: pair_i(:), pair_j(:)
      integer :: shls(4)
      integer :: ish, jsh, di, dj, i, j, g, ret, io, jo, block_max, n_grid
      integer :: grid_offset, n_pair, p
      real(dp) :: dij, wt

      if (size(grids, 1) /= 3) then
         call error%set(ERROR_VALIDATION, "ESP contraction: grid points must be (3, n)")
         return
      end if
      n_grid = size(grids, 2)
      if (n_grid < 1) then
         call error%set(ERROR_VALIDATION, "ESP contraction: no grid points")
         return
      end if

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

      allocate (values(n_grid))
      values = 0.0_dp

      ! Half the shell pairs: the integral is symmetric in the bra and ket for
      ! real functions, and so is the density contracted against it, so an
      ! off-diagonal block stands for itself and its transpose. Walked as a flat
      ! pair list rather than as a nested triangle so the schedule balances.
      n_pair = mol%nbas*(mol%nbas + 1)/2
      allocate (pair_i(n_pair), pair_j(n_pair))
      p = 0
      do ish = 1, mol%nbas
         do jsh = 1, ish
            p = p + 1
            pair_i(p) = ish
            pair_j(p) = jsh
         end do
      end do

      ! Threaded over shell pairs, each thread accumulating its own grid vector.
      !$omp parallel default(none) &
      !$omp    shared(mol, env_local, density, values, n_grid, block_max, &
      !$omp           n_pair, pair_i, pair_j) &
      !$omp    private(ish, jsh, di, dj, i, j, g, ret, io, jo, shls, buf, local, dij, p, wt)
      allocate (buf(block_max*block_max*n_grid), local(n_grid))
      local = 0.0_dp
      !$omp do schedule(dynamic)
      do p = 1, n_pair
         ish = pair_i(p)
         jsh = pair_j(p)
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
         jo = mol%shell_offset(jsh)
         ! A diagonal block already holds both `(i,j)` and `(j,i)`; an off-diagonal
         ! one holds only its own half of the pair.
         wt = 2.0_dp
         if (ish == jsh) wt = 1.0_dp
         shls = [ish - 1, jsh - 1, 0, n_grid]

         if (mol%cartesian) then
            ret = libcint_1e_grids_cart(buf, shls, mol%atm, mol%natm, &
                                        mol%bas, mol%nbas, env_local)
         else
            ret = libcint_1e_grids_sph(buf, shls, mol%atm, mol%natm, &
                                       mol%bas, mol%nbas, env_local)
         end if
         if (ret == 0) cycle

         do j = 1, dj
            do i = 1, di
               dij = wt*density(io + i, jo + j)
               if (dij == 0.0_dp) cycle
               do g = 1, n_grid
                  local(g) = local(g) &
                             + dij*buf(g + (i - 1)*n_grid + (j - 1)*n_grid*di)
               end do
            end do
         end do
      end do
      !$omp end do
      !$omp critical(mqc_esp_contract_accumulate)
      values = values + local
      !$omp end critical(mqc_esp_contract_accumulate)
      deallocate (buf, local)
      !$omp end parallel

      ! The electronic potential is negative where the density is positive.
      values = -values
      deallocate (env_local, pair_i, pair_j)
   end subroutine esp_contract

   subroutine drinv_matrices(mol, centres, matrices, error)
      !! The gradient of `1/|r - C|` with respect to `C`, over every shell pair
      !!
      !! Returns `(n_ao, n_ao, 3, n_centres)`. This is what a point *dipole*'s
      !! potential needs: the field of a dipole `d` at `C` is
      !! `d . grad_C (1/|r - C|)`.
      !!
      !! **The operator carries the derivative, not the basis.** libcint's
      !! `nabla-rinv` is `grad` applied to `1/r_C` itself, with the Gaussians left
      !! alone, so there is no integration by parts and no pairing with a
      !! transpose.
      !!
      !! **The sign convention is the test's.** The gradient could be with respect
      !! to `r` or to `C`, and those differ by a minus; `test_mqc_czt_esp`
      !! differences `esp_matrices` about each centre and requires this to
      !! reproduce it, which is what defines the convention for every caller.
      !!
      !! Unlike the grid integral this takes **one centre per call**, in
      !! `env(PTR_RINV_ORIG)`, so `env` is rewritten per centre.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: centres(:, :)        !! (3, n_centres), Bohr
      real(dp), allocatable, intent(out) :: matrices(:, :, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: buf(:), env_local(:)
      integer :: shls(2)
      integer :: ish, jsh, di, dj, i, j, c, x, ret, io, jo, block_max, n_centre

      if (size(centres, 1) /= 3) then
         call error%set(ERROR_VALIDATION, "drinv integrals: centres must be (3, n)")
         return
      end if
      n_centre = size(centres, 2)
      if (n_centre < 1) then
         call error%set(ERROR_VALIDATION, "drinv integrals: no centres")
         return
      end if

      block_max = 1
      do ish = 1, mol%nbas
         block_max = max(block_max, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do

      allocate (matrices(mol%nao, mol%nao, 3, n_centre))
      matrices = 0.0_dp
      allocate (buf(3*block_max*block_max))
      allocate (env_local(size(mol%env)))
      env_local = mol%env

      do c = 1, n_centre
         env_local(LIBCINT_PTR_RINV_ORIG + 1:LIBCINT_PTR_RINV_ORIG + 3) = centres(:, c)
         do ish = 1, mol%nbas
            di = shell_dim(mol%cartesian, ish - 1, mol%bas)
            io = mol%shell_offset(ish)
            do jsh = 1, mol%nbas
               dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
               jo = mol%shell_offset(jsh)
               shls = [ish - 1, jsh - 1]

               if (mol%cartesian) then
                  ret = libcint_1e_drinv_cart(buf, shls, mol%atm, mol%natm, &
                                              mol%bas, mol%nbas, env_local)
               else
                  ret = libcint_1e_drinv_sph(buf, shls, mol%atm, mol%natm, &
                                             mol%bas, mol%nbas, env_local)
               end if
               if (ret == 0) cycle      ! screened away, leave the block zero

               ! The Cartesian component is *slowest* here, the opposite of the
               ! grid integral where the point runs fastest.
               do x = 1, 3
                  do j = 1, dj
                     do i = 1, di
                        matrices(io + i, jo + j, x, c) = &
                           buf(i + (j - 1)*di + (x - 1)*di*dj)
                     end do
                  end do
               end do
            end do
         end do
      end do

      deallocate (buf, env_local)
   end subroutine drinv_matrices

   subroutine ddrinv_matrices(mol, centres, matrices, error)
      !! The second centre-derivative of `1/|r - C|`, over every shell pair
      !!
      !! Returns `(n_ao, n_ao, 3, 3, n_centres)`. This is what a point *quadrupole*'s
      !! potential needs, one gradient above `drinv_matrices`.
      !!
      !! **libcint has no `ddrinv`,** so unlike the first derivative this one is
      !! assembled from basis derivatives. Writing `f = 1/|r - C|` and using
      !! `grad_C f = -grad_r f` with an integration by parts whose boundary term
      !! vanishes,
      !!
      !!     d/dC_a I    = <d_a mu| f |nu> + <mu| f |d_a nu>
      !!     d2/dC_a dC_b I = <d_a d_b mu| f |nu> + <d_a mu| f |d_b nu>
      !!                    + <d_b mu| f |d_a nu> + <mu| f |d_a d_b nu>
      !!
      !! which is `int1e_ipiprinv` (`nabla nabla | rinv |`) plus its transpose, plus
      !! `int1e_iprinvip` (`nabla | rinv | nabla`) plus *its* index-swapped self.
      !!
      !! **Two of libcint's conventions cancel here rather than needing to be
      !! known.** Whether its `nabla` carries a minus does not matter, since every
      !! term above carries exactly two of them; and whether its nine components
      !! run `a` fastest or `b` fastest does not matter either, since the
      !! expression is symmetrized over `(a,b)`. Both questions are live for
      !! `drinv_matrices`, which has one derivative.
      !!
      !! `test_mqc_czt_esp` requires this to reproduce a second difference of
      !! `esp_matrices`.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: centres(:, :)        !! (3, n_centres), Bohr
      real(dp), allocatable, intent(out) :: matrices(:, :, :, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: buf_pp(:), buf_pip(:), env_local(:)
      integer :: shls(2)
      integer :: ish, jsh, di, dj, i, j, c, a, b, io, jo, block_max, n_centre
      integer :: slot, ret_pp, ret_pip

      if (size(centres, 1) /= 3) then
         call error%set(ERROR_VALIDATION, "ddrinv integrals: centres must be (3, n)")
         return
      end if
      n_centre = size(centres, 2)
      if (n_centre < 1) then
         call error%set(ERROR_VALIDATION, "ddrinv integrals: no centres")
         return
      end if

      block_max = 1
      do ish = 1, mol%nbas
         block_max = max(block_max, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do

      allocate (matrices(mol%nao, mol%nao, 3, 3, n_centre))
      matrices = 0.0_dp
      allocate (buf_pp(9*block_max*block_max))
      allocate (buf_pip(9*block_max*block_max))
      allocate (env_local(size(mol%env)))
      env_local = mol%env

      do c = 1, n_centre
         env_local(LIBCINT_PTR_RINV_ORIG + 1:LIBCINT_PTR_RINV_ORIG + 3) = centres(:, c)
         do ish = 1, mol%nbas
            di = shell_dim(mol%cartesian, ish - 1, mol%bas)
            io = mol%shell_offset(ish)
            do jsh = 1, mol%nbas
               dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
               jo = mol%shell_offset(jsh)
               shls = [ish - 1, jsh - 1]

               ! Both are asked for unconditionally and each zeroes its own buffer when
               ! screened away. Chaining them on one return code would drop a live
               ! contribution whenever the other one happened to be screened.
               if (mol%cartesian) then
                  ret_pp = libcint_1e_ipiprinv_cart(buf_pp, shls, mol%atm, mol%natm, &
                                                    mol%bas, mol%nbas, env_local)
                  ret_pip = libcint_1e_iprinvip_cart(buf_pip, shls, mol%atm, mol%natm, &
                                                     mol%bas, mol%nbas, env_local)
               else
                  ret_pp = libcint_1e_ipiprinv_sph(buf_pp, shls, mol%atm, mol%natm, &
                                                   mol%bas, mol%nbas, env_local)
                  ret_pip = libcint_1e_iprinvip_sph(buf_pip, shls, mol%atm, mol%natm, &
                                                    mol%bas, mol%nbas, env_local)
               end if
               if (ret_pp == 0 .and. ret_pip == 0) cycle
               if (ret_pp == 0) buf_pp = 0.0_dp
               if (ret_pip == 0) buf_pip = 0.0_dp

               do b = 1, 3
                  do a = 1, 3
                     slot = (a - 1) + (b - 1)*3
                     do j = 1, dj
                        do i = 1, di
                           ! <d_a d_b mu|f|nu> + <d_a mu|f|d_b nu>, the two terms whose
                           ! derivatives sit on this shell pair in this order; the other
                           ! two arrive when the (jsh, ish) block is visited, since this
                           ! loop covers every ordered pair.
                           matrices(io + i, jo + j, a, b, c) = &
                              matrices(io + i, jo + j, a, b, c) &
                              + buf_pp(i + (j - 1)*di + slot*di*dj) &
                              + buf_pip(i + (j - 1)*di + slot*di*dj)
                           ! ... and their partners, placed transposed: <mu|f|d_a d_b nu>
                           ! and <d_b mu|f|d_a nu> as seen from the other block.
                           matrices(jo + j, io + i, a, b, c) = &
                              matrices(jo + j, io + i, a, b, c) &
                              + buf_pp(i + (j - 1)*di + slot*di*dj) &
                              + buf_pip(i + (j - 1)*di + slot*di*dj)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      deallocate (buf_pp, buf_pip, env_local)
   end subroutine ddrinv_matrices

   subroutine esp_matrices(mol, grids, matrices, error)
      !! `1/|r - R_g|` over every shell pair, for each point
      !!
      !! Returns `(n_ao, n_ao, n_points)`. Contract with a density and negate to get
      !! the electronic potential; the nuclear part is `sum_A Z_A/|R_g - R_A|` and is
      !! the caller's, since a screening fit wants the electronic part alone.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: grids(:, :)          !! (3, n_points), Bohr
      real(dp), allocatable, intent(out) :: matrices(:, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: buf(:), env_local(:)
      integer :: shls(4)
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
               ret = libcint_1e_grids_cart(buf, shls, mol%atm, mol%natm, &
                                           mol%bas, mol%nbas, env_local)
            else
               ret = libcint_1e_grids_sph(buf, shls, mol%atm, mol%natm, &
                                          mol%bas, mol%nbas, env_local)
            end if
            if (ret == 0) cycle      ! screened away, leave the block zero

            ! The grid index runs fastest, then the bra, then the ket -- the
            ! opposite nesting from the multipole blocks, where the component is
            ! slowest. Backwards it gives a matrix that is plausible at one point
            ! and wrong at every other, so the check against PySCF is elementwise
            ! over all points rather than on a contracted potential.
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

end module mqc_czt_esp
