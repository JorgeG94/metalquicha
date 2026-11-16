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
   !! A grid small enough to be worth using integrates products of basis
   !! functions poorly, so the grids here carry an overlap fit: a matrix,
   !! applied to the left index of `K`, that makes the same quadrature exact
   !! for the overlap. See `cosx_overlap_fit`.
   !!
   !! Every array is indexed (point, function), points in the order given. The
   !! steps are separate public routines over plain arrays, and
   !! `cosx_esp_contract` is the only one that evaluates an integral.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_calculation_defaults, only: DEFAULT_COSX_GRID
   use mqc_czt_integrals, only: czt_molecule_t, shell_dim
   use mqc_czt_ao, only: eval_ao_block
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_dft_grid, only: dft_grid_t, build_dft_grid, element_period
   use libcint_fortran, only: LIBCINT_PTR_GRIDS, &
                              libcint_1e_grids_cart, libcint_1e_grids_sph
   implicit none
   private

   public :: cosx_grid_t
   public :: build_cosx_grid
   public :: cosx_exchange
   public :: cosx_overlap_fit
   public :: cosx_numerical_overlap
   public :: cosx_batch_exchange
   public :: cosx_density_on_grid
   public :: cosx_esp_contract
   public :: cosx_accumulate
   public :: COSX_DEFAULT_BATCH
   public :: COSX_DEFAULT_SCREEN
   public :: COSX_DEFAULT_GRID_LEVEL, COSX_MIN_GRID_LEVEL, COSX_MAX_GRID_LEVEL
   public :: COSX_ROLE_SCF, COSX_ROLE_FINAL
   public :: cosx_radial_points, cosx_angular_points

   integer, parameter :: COSX_DEFAULT_BATCH = 128
      !! Points per batch when the caller does not say
   real(dp), parameter :: COSX_DEFAULT_SCREEN = 1.0e-10_dp
      !! A shell pair's potential integral is skipped for a batch when its
      !! Schwarz bound times the largest element of `F` it would multiply is
      !! below this
   real(dp), parameter :: COSX_FIT_NULL = 1.0e-10_dp
      !! Eigenvalues of the numerical overlap below this fraction of the
      !! largest are left out of its inverse: directions of the basis the grid
      !! cannot resolve, which a plain inverse would amplify
   integer, parameter :: COSX_MIN_GRID_LEVEL = 1
   integer, parameter :: COSX_MAX_GRID_LEVEL = 3
   integer, parameter :: COSX_DEFAULT_GRID_LEVEL = DEFAULT_COSX_GRID
      !! Level of the exchange grids, 1 (smallest) to 3. Not the DFT grid's
      !! scale: an exchange grid is far smaller than one for a functional.

   integer, parameter :: COSX_ROLE_SCF = 1     !! The grid the iterations run on
   integer, parameter :: COSX_ROLE_FINAL = 2   !! The larger one the final energy is taken on

   integer, parameter :: COSX_RADIAL_BASE(2, COSX_MAX_GRID_LEVEL) = &
                         reshape([5, 20, 10, 20, 10, 25], [2, COSX_MAX_GRID_LEVEL])
      !! (role, level). Radial shells are this plus `COSX_RADIAL_PER_PERIOD`
      !! times the period. The iteration grid is radially much coarser than
      !! the final one: its own exchange energy is poor, but the density it
      !! converges is not, and the energy is taken on the final grid.
   integer, parameter :: COSX_RADIAL_PER_PERIOD = 5
   integer, parameter :: COSX_ANGULAR(2, COSX_MAX_GRID_LEVEL) = &
                         reshape([50, 86, 50, 194, 86, 302], [2, COSX_MAX_GRID_LEVEL])
      !! (role, level). Lebedev points in the valence zone; the pruning steps
      !! down from it towards the nucleus and far from it.

   type :: cosx_grid_t
      !! The points the exchange is integrated over
      !!
      !! Unallocated means the exchange is not seminumerical, which is how a
      !! Fock build reads it.
      integer :: n_points = 0
      integer :: level = -1                   !! The level it was built at
      real(dp), allocatable :: coords(:, :)   !! (3, n_points), Bohr
      real(dp), allocatable :: weights(:)     !! (n_points)
      real(dp), allocatable :: fit(:, :)
         !! (n_ao, n_ao), the overlap-fitting matrix from `cosx_overlap_fit`.
         !! Unallocated means the quadrature is used as it stands.
      real(dp), allocatable :: bounds(:, :)
         !! (n_shell, n_shell) Schwarz bounds, `sqrt(max |(MN|MN)|)`, that the
         !! contraction screens shell pairs with. Unallocated means every pair
         !! is evaluated.
   contains
      procedure :: active => cosx_grid_active
      procedure :: destroy => cosx_grid_destroy
   end type cosx_grid_t

contains

   pure function cosx_grid_active(this) result(active)
      !! Whether a grid was built
      class(cosx_grid_t), intent(in) :: this
      logical :: active

      active = allocated(this%coords)
   end function cosx_grid_active

   pure subroutine cosx_grid_destroy(this)
      !! Release the points
      class(cosx_grid_t), intent(inout) :: this

      if (allocated(this%coords)) deallocate (this%coords)
      if (allocated(this%weights)) deallocate (this%weights)
      if (allocated(this%fit)) deallocate (this%fit)
      if (allocated(this%bounds)) deallocate (this%bounds)
      this%n_points = 0
      this%level = -1
   end subroutine cosx_grid_destroy

   pure function cosx_radial_points(atomic_number, level, role) result(n_radial)
      !! Radial shells an element gets on an exchange grid
      integer, intent(in) :: atomic_number
      integer, intent(in) :: level   !! Clamped to the levels that exist
      integer, intent(in) :: role    !! `COSX_ROLE_SCF` or `COSX_ROLE_FINAL`
      integer :: n_radial

      n_radial = COSX_RADIAL_BASE(role, clamp_level(level)) &
                 + COSX_RADIAL_PER_PERIOD*element_period(atomic_number)
   end function cosx_radial_points

   pure function cosx_angular_points(level, role) result(n_angular)
      !! Lebedev points in the valence zone of an exchange grid
      integer, intent(in) :: level
      integer, intent(in) :: role
      integer :: n_angular

      n_angular = COSX_ANGULAR(role, clamp_level(level))
   end function cosx_angular_points

   pure function clamp_level(level) result(clamped)
      integer, intent(in) :: level
      integer :: clamped

      clamped = max(COSX_MIN_GRID_LEVEL, min(COSX_MAX_GRID_LEVEL, level))
   end function clamp_level

   subroutine build_cosx_grid(mol, grid, error, level, role, n_radial, n_angular, overlap, &
                              screen)
      !! The molecular grid for the exchange
      !!
      !! Sized from the exchange grid tables by `level` and `role`, or from
      !! `n_radial` and `n_angular` when those are given. Ghost centres get
      !! points: they carry basis functions, and the quadrature integrates
      !! products of basis functions.
      type(czt_molecule_t), intent(in) :: mol
      type(cosx_grid_t), intent(out) :: grid
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: level   !! Default `COSX_DEFAULT_GRID_LEVEL`
      integer, intent(in), optional :: role    !! Default `COSX_ROLE_SCF`
      integer, intent(in), optional :: n_radial
         !! Radial shells on every atom, overriding the tables. Given with
         !! `n_angular` or not at all.
      integer, intent(in), optional :: n_angular
         !! Lebedev points in the valence zone of every atom
      real(dp), intent(in), optional :: overlap(:, :)
         !! The analytic overlap. Present, the grid carries the overlap fit
         !! and every exchange built on it is fitted.
      logical, intent(in), optional :: screen
         !! Build the Schwarz bounds the contraction screens with. Default
         !! false: one quartet per shell pair, worth paying only for a grid an
         !! SCF will reuse.

      type(dft_grid_t) :: full
      integer, allocatable :: numbers(:), radial(:), angular(:)
      integer :: used_role, ia

      grid%level = COSX_DEFAULT_GRID_LEVEL
      if (present(level)) grid%level = level
      used_role = COSX_ROLE_SCF
      if (present(role)) used_role = role
      if (used_role /= COSX_ROLE_SCF .and. used_role /= COSX_ROLE_FINAL) then
         call error%set(ERROR_VALIDATION, "COSX: unknown grid role")
         return
      end if
      if (grid%level < COSX_MIN_GRID_LEVEL .or. grid%level > COSX_MAX_GRID_LEVEL) then
         call error%set(ERROR_VALIDATION, "COSX: the grid level is "// &
                        to_char(grid%level)//" and must be "// &
                        to_char(COSX_MIN_GRID_LEVEL)//" to "//to_char(COSX_MAX_GRID_LEVEL))
         return
      end if

      ! The element and not the charge it presents, as the XC grid does: the
      ! charge is reduced by an ECP's core and is zero on a ghost.
      allocate (numbers(mol%natm))
      if (allocated(mol%atomic_numbers)) then
         numbers = mol%atomic_numbers
      else
         numbers = nint(mol%charges) + mol%core_electrons
      end if

      if (present(n_radial) .or. present(n_angular)) then
         call build_dft_grid(mol%coords, numbers, full, error, &
                             n_radial=n_radial, n_angular=n_angular)
         grid%level = -1
      else
         allocate (radial(mol%natm), angular(mol%natm))
         do ia = 1, mol%natm
            radial(ia) = cosx_radial_points(numbers(ia), grid%level, used_role)
            angular(ia) = cosx_angular_points(grid%level, used_role)
         end do
         call build_dft_grid(mol%coords, numbers, full, error, &
                             radial_by_atom=radial, angular_by_atom=angular)
      end if
      if (error%has_error()) return
      grid%n_points = full%n_points
      call move_alloc(full%coords, grid%coords)
      call move_alloc(full%weights, grid%weights)
      call full%destroy()
      if (present(overlap)) then
         call cosx_overlap_fit(mol, grid%coords, grid%weights, overlap, grid%fit, error)
         if (error%has_error()) return
      end if
      if (present(screen)) then
         if (screen) call schwarz_bounds(mol, grid%bounds, error)
      end if
   end subroutine build_cosx_grid

   subroutine cosx_exchange(mol, coords, weights, density, k, error, batch_size, symmetrize, fit, &
                            bounds, screen_tol)
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
      real(dp), intent(in), optional :: fit(:, :)
         !! (n_ao, n_ao) from `cosx_overlap_fit` on the same points. Applied
         !! to the left index of `K` before any symmetrisation, which makes the
         !! quadrature exact for the overlap and far more accurate for the
         !! exchange on a small grid. Absent is the plain quadrature.
      real(dp), intent(in), optional :: bounds(:, :)
         !! (n_shell, n_shell) Schwarz bounds from `schwarz_bounds`. Present,
         !! a shell pair is skipped for a batch whose `F` makes its
         !! contribution negligible. Absent, every pair is evaluated.
      real(dp), intent(in), optional :: screen_tol
         !! The threshold for that, default `COSX_DEFAULT_SCREEN`. Zero keeps
         !! every pair. Read only with `bounds`.

      real(dp), allocatable :: k_local(:, :), k_plain(:, :)
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
      if (present(fit)) then
         if (size(fit, 1) /= n_ao .or. size(fit, 2) /= n_ao) then
            call error%set(ERROR_VALIDATION, "COSX: the overlap fit does not match the basis")
            return
         end if
      end if

      allocate (k(n_ao, n_ao))
      k = 0.0_dp
      failed = .false.

      if (present(bounds)) then
         if (size(bounds, 1) /= mol%nbas .or. size(bounds, 2) /= mol%nbas) then
            call error%set(ERROR_VALIDATION, "COSX: the Schwarz bounds do not match the basis")
            deallocate (k)
            return
         end if
      end if

      ! `default(shared)`: `bounds` and `screen_tol` are optional dummies, and
      ! naming an absent one in a data-sharing clause is not portable.
      !$omp parallel default(shared) &
      !$omp    private(g0, g1, k_local) &
      !$omp    firstprivate(local_error)
      allocate (k_local(n_ao, n_ao))
      k_local = 0.0_dp
      !$omp do schedule(dynamic)
      do g0 = 1, n_points, step
         if (failed) cycle
         g1 = min(g0 + step - 1, n_points)
         call cosx_batch_exchange(mol, coords(:, g0:g1), weights(g0:g1), density, &
                                  k_local, local_error, bounds=bounds, &
                                  screen_tol=screen_tol)
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
      if (present(fit)) then
         call move_alloc(k, k_plain)
         allocate (k(n_ao, n_ao))
         call pic_gemm(fit, k_plain, k)
         deallocate (k_plain)
      end if
      if (fold) k = 0.5_dp*(k + transpose(k))
   end subroutine cosx_exchange

   subroutine cosx_numerical_overlap(mol, coords, weights, s_grid, error, batch_size)
      !! The overlap matrix as the grid integrates it
      !!
      !! `S(mu, nu) = sum_g w_g chi_mu(r_g) chi_nu(r_g)`. Equal to the analytic
      !! overlap in the limit of the grid, and how far it is from it on a given
      !! grid is how well that grid integrates a product of basis functions.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)       !! (3, n_points), Bohr
      real(dp), intent(in) :: weights(:)         !! (n_points)
      real(dp), allocatable, intent(out) :: s_grid(:, :)   !! (n_ao, n_ao)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: batch_size   !! Default `COSX_DEFAULT_BATCH`

      real(dp), allocatable :: s_local(:, :), ao(:, :), wao(:, :)
      type(error_t) :: local_error
      integer :: n_points, n_ao, step, g0, g1, ig, mu
      logical :: failed

      n_ao = mol%nao
      n_points = size(coords, 2)
      if (size(coords, 1) /= 3 .or. n_points < 1 .or. size(weights) /= n_points) then
         call error%set(ERROR_VALIDATION, "COSX: the grid must be (3, n) points with "// &
                        "n > 0 and one weight each")
         return
      end if
      step = COSX_DEFAULT_BATCH
      if (present(batch_size)) step = batch_size
      if (step < 1) then
         call error%set(ERROR_VALIDATION, "COSX: the batch size must be positive")
         return
      end if

      allocate (s_grid(n_ao, n_ao))
      s_grid = 0.0_dp
      failed = .false.

      !$omp parallel default(none) &
      !$omp    shared(mol, coords, weights, s_grid, error, failed, n_points, n_ao, step) &
      !$omp    private(g0, g1, ig, mu, s_local, ao, wao) &
      !$omp    firstprivate(local_error)
      allocate (s_local(n_ao, n_ao))
      s_local = 0.0_dp
      !$omp do schedule(dynamic)
      do g0 = 1, n_points, step
         if (failed) cycle
         g1 = min(g0 + step - 1, n_points)
         call eval_ao_block(mol, coords(:, g0:g1), ao, local_error)
         if (local_error%has_error()) then
            !$omp critical(mqc_cosx_overlap_failure)
            if (.not. failed) then
               failed = .true.
               error = local_error
            end if
            !$omp end critical(mqc_cosx_overlap_failure)
            cycle
         end if
         if (allocated(wao)) deallocate (wao)
         allocate (wao(size(ao, 1), n_ao))
         do mu = 1, n_ao
            do ig = 1, size(ao, 1)
               wao(ig, mu) = weights(g0 + ig - 1)*ao(ig, mu)
            end do
         end do
         call pic_gemm(ao, wao, s_local, transa="T", beta=1.0_dp)
      end do
      !$omp end do
      !$omp critical(mqc_cosx_overlap_accumulate)
      s_grid = s_grid + s_local
      !$omp end critical(mqc_cosx_overlap_accumulate)
      deallocate (s_local)
      !$omp end parallel

      if (failed) deallocate (s_grid)
   end subroutine cosx_numerical_overlap

   subroutine cosx_overlap_fit(mol, coords, weights, overlap, fit, error, n_dropped)
      !! The matrix that makes a grid integrate the overlap exactly
      !!
      !! `fit = S_analytic S_grid^-1`, with `S_grid` from
      !! `cosx_numerical_overlap`. Replacing the function values on one side of
      !! the quadrature by `chi fit^T` turns the numerical overlap into the
      !! analytic one, and applied to the exchange that is `K <- fit K`.
      !!
      !! It depends on the grid and the basis and not on the density, so it is
      !! built once per grid. It is the identity in the limit of the grid.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)       !! (3, n_points), Bohr
      real(dp), intent(in) :: weights(:)         !! (n_points)
      real(dp), intent(in) :: overlap(:, :)      !! (n_ao, n_ao), analytic
      real(dp), allocatable, intent(out) :: fit(:, :)   !! (n_ao, n_ao)
      type(error_t), intent(inout) :: error
      integer, intent(out), optional :: n_dropped
         !! Directions of the numerical overlap left out of its inverse as
         !! unresolved by the grid. Zero on any grid adequate for the basis.

      real(dp), allocatable :: s_grid(:, :), values(:), scaled(:, :), inverse(:, :)
      integer :: n_ao, info, i, dropped
      real(dp) :: floor_value

      n_ao = mol%nao
      if (size(overlap, 1) /= n_ao .or. size(overlap, 2) /= n_ao) then
         call error%set(ERROR_VALIDATION, "COSX: the overlap does not match the basis")
         return
      end if
      call cosx_numerical_overlap(mol, coords, weights, s_grid, error)
      if (error%has_error()) return

      ! A pseudo-inverse through the eigenvectors: on a coarse grid the
      ! numerical overlap of a large basis can be singular where the analytic
      ! one is not.
      allocate (values(n_ao))
      call pic_syev(s_grid, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "COSX: the numerical overlap could not be "// &
                        "diagonalized (info = "//to_char(info)//")")
         return
      end if
      floor_value = COSX_FIT_NULL*maxval(abs(values))
      allocate (scaled(n_ao, n_ao), inverse(n_ao, n_ao))
      dropped = 0
      do i = 1, n_ao
         if (values(i) > floor_value) then
            scaled(:, i) = s_grid(:, i)/values(i)
         else
            scaled(:, i) = 0.0_dp
            dropped = dropped + 1
         end if
      end do
      if (present(n_dropped)) n_dropped = dropped
      call pic_gemm(scaled, s_grid, inverse, transb="T")
      allocate (fit(n_ao, n_ao))
      call pic_gemm(overlap, inverse, fit)
   end subroutine cosx_overlap_fit

   subroutine cosx_batch_exchange(mol, coords, weights, density, k, error, bounds, screen_tol)
      !! Add one batch of points to an exchange matrix
      !!
      !! The whole method for the points given, with no threading of its own.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)     !! (3, n_batch), Bohr
      real(dp), intent(in) :: weights(:)       !! (n_batch)
      real(dp), intent(in) :: density(:, :)    !! (n_ao, n_ao)
      real(dp), intent(inout) :: k(:, :)       !! (n_ao, n_ao), accumulated into
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: bounds(:, :)   !! See `cosx_exchange`
      real(dp), intent(in), optional :: screen_tol     !! See `cosx_exchange`

      real(dp), allocatable :: ao(:, :), f(:, :), g(:, :)

      call eval_ao_block(mol, coords, ao, error)
      if (error%has_error()) return
      call cosx_density_on_grid(ao, weights, density, f)
      call cosx_esp_contract(mol, coords, f, g, error, bounds=bounds, screen_tol=screen_tol)
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

   subroutine cosx_esp_contract(mol, coords, f, g, error, bounds, screen_tol, n_skipped)
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
      real(dp), intent(in), optional :: bounds(:, :)
         !! (n_shell, n_shell) Schwarz bounds. Present, the half of a shell
         !! pair that would add `A F(:, shell)` is skipped when the bound times
         !! the largest `|F|` of that shell over the batch is below
         !! `screen_tol`, and the integral is not evaluated when both halves are.
      real(dp), intent(in), optional :: screen_tol
         !! Default `COSX_DEFAULT_SCREEN`; zero keeps everything
      integer, intent(out), optional :: n_skipped
         !! Shell pairs, of the lower triangle, whose integral was not
         !! evaluated for this batch

      real(dp), allocatable :: env_local(:), buf(:), f_max(:)
      real(dp) :: tol
      logical :: screening, into_i, into_j
      integer :: skipped
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

      ! The largest element of F each shell carries over this batch. F falls
      ! off with the distance of a shell from the batch as the density matrix
      ! does, which is where the linear scaling comes from.
      skipped = 0
      screening = present(bounds)
      tol = COSX_DEFAULT_SCREEN
      if (present(screen_tol)) tol = screen_tol
      if (screening) then
         allocate (f_max(mol%nbas))
         do ish = 1, mol%nbas
            di = shell_dim(mol%cartesian, ish - 1, mol%bas)
            io = mol%shell_offset(ish)
            f_max(ish) = maxval(abs(f(:, io + 1:io + di)))
         end do
      end if

      ! The integral is symmetric in its two functions, so the lower triangle
      ! of shell pairs serves both `G(:, i) += A F(:, j)` and its mirror.
      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, ish
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            ! Which halves are worth having: G(:, i) takes F(:, j) and its
            ! mirror the reverse.
            into_i = .true.
            into_j = ish /= jsh
            if (screening) then
               into_i = bounds(ish, jsh)*f_max(jsh) >= tol
               into_j = into_j .and. bounds(ish, jsh)*f_max(ish) >= tol
               if (.not. (into_i .or. into_j)) then
                  skipped = skipped + 1
                  cycle
               end if
            end if
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
                  if (into_i) then
                     do ig = 1, n_batch
                        g(ig, io + i) = g(ig, io + i) + buf(base + ig)*f(ig, jo + j)
                     end do
                  end if
                  if (into_j) then
                     do ig = 1, n_batch
                        g(ig, jo + j) = g(ig, jo + j) + buf(base + ig)*f(ig, io + i)
                     end do
                  end if
               end do
            end do
         end do
      end do
      if (present(n_skipped)) n_skipped = skipped
   end subroutine cosx_esp_contract

   subroutine cosx_accumulate(ao, g, k)
      !! `K(mu, nu) += sum_g chi_mu(r_g) G(g, nu)`
      real(dp), intent(in) :: ao(:, :)        !! (n_batch, n_ao)
      real(dp), intent(in) :: g(:, :)         !! (n_batch, n_ao)
      real(dp), intent(inout) :: k(:, :)      !! (n_ao, n_ao)

      call pic_gemm(ao, g, k, transa="T", beta=1.0_dp)
   end subroutine cosx_accumulate

end module mqc_czt_cosx
