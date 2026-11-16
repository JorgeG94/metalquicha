!! Chain-of-spheres exchange: the exchange matrix by quadrature
module mqc_czt_cosx
   !! Seminumerical exchange. One electron of each two-electron integral is
   !! integrated on a grid and the other analytically:
   !!
   !!     (mu kappa | nu lambda) ~ sum_g w_g chi_mu(r_g) chi_kappa(r_g) A_g(nu, lambda)
   !!     A_g(nu, lambda) = <nu| 1/|r - r_g| |lambda>
   !!
   !! so that `K(mu, nu) = sum_(kappa lambda) P(kappa, lambda) (mu kappa | nu lambda)`
   !! is three steps per batch of points:
   !!
   !!     F(g, lambda) = w_g sum_kappa chi_kappa(r_g) P(kappa, lambda)
   !!     G(g, nu)     = sum_lambda A_g(nu, lambda) F(g, lambda)
   !!     K(mu, nu)   += sum_g chi_mu(r_g) G(g, nu)
   !!
   !! `K` is the exchange matrix itself, with no factor: a closed-shell Fock
   !! matrix takes `-K/2` of the total density, and the exchange energy is
   !! `-Tr(P K)/4`.
   !!
   !! Every array is indexed (point, function), points in the order given. The
   !! steps are separate public routines over plain arrays, and
   !! `cosx_esp_contract` is the only one that evaluates an integral.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, shell_dim
   use mqc_czt_ao, only: eval_ao_block
   use libcint_fortran, only: LIBCINT_PTR_GRIDS, &
                              libcint_1e_grids_cart, libcint_1e_grids_sph
   implicit none
   private

   public :: cosx_exchange
   public :: cosx_batch_exchange
   public :: cosx_density_on_grid
   public :: cosx_esp_contract
   public :: cosx_accumulate
   public :: COSX_DEFAULT_BATCH

   integer, parameter :: COSX_DEFAULT_BATCH = 128
      !! Points per batch when the caller does not say

contains

   subroutine cosx_exchange(mol, coords, weights, density, k, error, batch_size, symmetrize)
      !! The exchange matrix of a density, by quadrature over the points given
      !!
      !! Converges to the analytic exchange matrix as the grid grows. Threaded
      !! over batches, each calling the BLAS from its own thread.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)       !! (3, n_points), Bohr
      real(dp), intent(in) :: weights(:)         !! (n_points), quadrature weights
      real(dp), intent(in) :: density(:, :)      !! (n_ao, n_ao), AO basis
      real(dp), allocatable, intent(out) :: k(:, :)
         !! (n_ao, n_ao), `K(mu, nu) = sum P(kappa, lambda) (mu kappa | nu lambda)`
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: batch_size
         !! Points per batch, default `COSX_DEFAULT_BATCH`. Changes how the sum
         !! is grouped and not what it sums to.
      logical, intent(in), optional :: symmetrize
         !! Return `(K + K^T)/2`, default true. The quadrature treats the two
         !! electrons differently, so `K` is symmetric only in the limit of an
         !! exact grid even when the density is. Pass false for a density that
         !! is not symmetric, whose exchange matrix is not either.

      real(dp), allocatable :: k_local(:, :)
      type(error_t) :: local_error
      integer :: n_points, n_ao, step, g0, g1
      logical :: fold, failed

      n_ao = mol%nao
      if (size(coords, 1) /= 3) then
         call error%set(ERROR_VALIDATION, "COSX: grid points must be (3, n)")
         return
      end if
      n_points = size(coords, 2)
      if (n_points < 1) then
         call error%set(ERROR_VALIDATION, "COSX: no grid points")
         return
      end if
      if (size(weights) /= n_points) then
         call error%set(ERROR_VALIDATION, "COSX: one weight per grid point is required")
         return
      end if
      if (size(density, 1) /= n_ao .or. size(density, 2) /= n_ao) then
         call error%set(ERROR_VALIDATION, "COSX: the density does not match the basis")
         return
      end if

      step = COSX_DEFAULT_BATCH
      if (present(batch_size)) step = batch_size
      if (step < 1) then
         call error%set(ERROR_VALIDATION, "COSX: the batch size must be positive")
         return
      end if
      fold = .true.
      if (present(symmetrize)) fold = symmetrize

      allocate (k(n_ao, n_ao))
      k = 0.0_dp
      failed = .false.

      !$omp parallel default(none) &
      !$omp    shared(mol, coords, weights, density, k, error, failed, &
      !$omp           n_points, n_ao, step) &
      !$omp    private(g0, g1, k_local) &
      !$omp    firstprivate(local_error)
      allocate (k_local(n_ao, n_ao))
      k_local = 0.0_dp
      !$omp do schedule(dynamic)
      do g0 = 1, n_points, step
         if (failed) cycle
         g1 = min(g0 + step - 1, n_points)
         call cosx_batch_exchange(mol, coords(:, g0:g1), weights(g0:g1), density, &
                                  k_local, local_error)
         if (local_error%has_error()) then
            !$omp critical(mqc_cosx_exchange_failure)
            if (.not. failed) then
               failed = .true.
               error = local_error
            end if
            !$omp end critical(mqc_cosx_exchange_failure)
         end if
      end do
      !$omp end do
      !$omp critical(mqc_cosx_exchange_accumulate)
      k = k + k_local
      !$omp end critical(mqc_cosx_exchange_accumulate)
      deallocate (k_local)
      !$omp end parallel

      if (failed) then
         deallocate (k)
         return
      end if
      if (fold) k = 0.5_dp*(k + transpose(k))
   end subroutine cosx_exchange

   subroutine cosx_batch_exchange(mol, coords, weights, density, k, error)
      !! Add one batch of points to an exchange matrix
      !!
      !! The whole method for the points given, with no threading of its own.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)     !! (3, n_batch), Bohr
      real(dp), intent(in) :: weights(:)       !! (n_batch)
      real(dp), intent(in) :: density(:, :)    !! (n_ao, n_ao)
      real(dp), intent(inout) :: k(:, :)       !! (n_ao, n_ao), accumulated into
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: ao(:, :), f(:, :), g(:, :)

      call eval_ao_block(mol, coords, ao, error)
      if (error%has_error()) return
      call cosx_density_on_grid(ao, weights, density, f)
      call cosx_esp_contract(mol, coords, f, g, error)
      if (error%has_error()) return
      call cosx_accumulate(ao, g, k)
   end subroutine cosx_batch_exchange

   subroutine cosx_density_on_grid(ao, weights, density, f)
      !! `F(g, lambda) = w_g sum_kappa chi_kappa(r_g) P(kappa, lambda)`
      real(dp), intent(in) :: ao(:, :)         !! (n_batch, n_ao), function values
      real(dp), intent(in) :: weights(:)       !! (n_batch)
      real(dp), intent(in) :: density(:, :)    !! (n_ao, n_ao)
      real(dp), allocatable, intent(out) :: f(:, :)   !! (n_batch, n_ao)

      integer :: ig, mu

      allocate (f(size(ao, 1), size(ao, 2)))
      call pic_gemm(ao, density, f)
      do mu = 1, size(f, 2)
         do ig = 1, size(f, 1)
            f(ig, mu) = weights(ig)*f(ig, mu)
         end do
      end do
   end subroutine cosx_density_on_grid

   subroutine cosx_esp_contract(mol, coords, f, g, error)
      !! `G(g, nu) = sum_lambda <nu| 1/|r - r_g| |lambda> F(g, lambda)`
      !!
      !! Point `g` of `f` is contracted with the potential integral at point
      !! `g` of `coords` and no other. The integral over all points for all
      !! function pairs is never held: each shell pair's block is contracted as
      !! it is produced.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)     !! (3, n_batch), Bohr
      real(dp), intent(in) :: f(:, :)          !! (n_batch, n_ao)
      real(dp), allocatable, intent(out) :: g(:, :)   !! (n_batch, n_ao)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: env_local(:), buf(:)
      integer :: shls(4)
      integer :: ish, jsh, di, dj, i, j, ig, io, jo, ret, block_max
      integer :: n_batch, grid_offset, base

      n_batch = size(coords, 2)
      if (size(coords, 1) /= 3 .or. n_batch < 1) then
         call error%set(ERROR_VALIDATION, "COSX: grid points must be (3, n) with n > 0")
         return
      end if
      if (size(f, 1) /= n_batch .or. size(f, 2) /= mol%nao) then
         call error%set(ERROR_VALIDATION, "COSX: F must be (n_points, n_ao)")
         return
      end if

      ! The points ride at the end of a copy of env, and PTR_GRIDS holds the
      ! 0-based offset of the first; see `mqc_czt_esp`.
      grid_offset = size(mol%env)
      allocate (env_local(grid_offset + 3*n_batch))
      env_local(1:grid_offset) = mol%env
      do ig = 1, n_batch
         env_local(grid_offset + 3*(ig - 1) + 1:grid_offset + 3*ig) = coords(:, ig)
      end do
      env_local(LIBCINT_PTR_GRIDS + 1) = real(grid_offset, dp)

      block_max = 1
      do ish = 1, mol%nbas
         block_max = max(block_max, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do
      allocate (buf(block_max*block_max*n_batch))

      allocate (g(n_batch, mol%nao))
      g = 0.0_dp

      ! The integral is symmetric in its two functions, so the lower triangle
      ! of shell pairs serves both `G(:, i) += A F(:, j)` and its mirror.
      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, ish
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            shls = [ish - 1, jsh - 1, 0, n_batch]
            if (mol%cartesian) then
               ret = libcint_1e_grids_cart(buf, shls, mol%atm, mol%natm, &
                                           mol%bas, mol%nbas, env_local)
            else
               ret = libcint_1e_grids_sph(buf, shls, mol%atm, mol%natm, &
                                          mol%bas, mol%nbas, env_local)
            end if
            if (ret == 0) cycle

            ! The buffer runs point fastest, then the bra, then the ket.
            do j = 1, dj
               do i = 1, di
                  base = (i - 1)*n_batch + (j - 1)*n_batch*di
                  do ig = 1, n_batch
                     g(ig, io + i) = g(ig, io + i) + buf(base + ig)*f(ig, jo + j)
                  end do
                  if (ish == jsh) cycle
                  do ig = 1, n_batch
                     g(ig, jo + j) = g(ig, jo + j) + buf(base + ig)*f(ig, io + i)
                  end do
               end do
            end do
         end do
      end do
   end subroutine cosx_esp_contract

   subroutine cosx_accumulate(ao, g, k)
      !! `K(mu, nu) += sum_g chi_mu(r_g) G(g, nu)`
      real(dp), intent(in) :: ao(:, :)        !! (n_batch, n_ao)
      real(dp), intent(in) :: g(:, :)         !! (n_batch, n_ao)
      real(dp), intent(inout) :: k(:, :)      !! (n_ao, n_ao)

      call pic_gemm(ao, g, k, transa="T", beta=1.0_dp)
   end subroutine cosx_accumulate

end module mqc_czt_cosx
