!! MBIS atomic charges, by a stockholder partition fitted to the molecule
module mqc_czt_mbis
   !! Minimal Basis Iterative Stockholder (Verstraelen et al., J. Chem. Theory
   !! Comput. 12, 3894 (2016)). A stockholder partition, as Hirshfeld's is,
   !!
   !!     w_A(r) = rho_A^0(r) / sum_B rho_B^0(r)
   !!
   !! but the pro-atoms are not free atoms. Each is a sum of exponential shells,
   !!
   !!     rho_A^0(r) = sum_i N_Ai a_Ai^3 / (8 pi) exp(-a_Ai |r - R_A|)
   !!
   !! one shell per row of the periodic table the element reaches, and the
   !! populations and exponents are fitted to the molecule's own density by the
   !! fixed point
   !!
   !!     N_Ai = int rho w_Ai,      a_Ai = 3 N_Ai / int rho w_Ai |r - R_A|
   !!
   !! with w_Ai the shell's share, rho_Ai^0 / rho^0. At convergence each shell
   !! holds what it is given and its exponent matches the mean radius of what it
   !! holds, and q_A = Z_A - sum_i N_Ai.
   !!
   !! So there is no free-atom calculation and nothing to choose about one, and
   !! no element table beyond the shell count: the fit starts from the shells
   !! filled in aufbau order, with exponents from a Slater-style screened
   !! charge. The
   !! partition is the density's own, and a charged molecule's atoms come out
   !! charged rather than being measured against neutral ones -- the reason MBIS
   !! charges are larger than Hirshfeld's.
   !!
   !! **The weights depend on the density.** Two densities of one molecule are
   !! divided by two different partitions, unlike Hirshfeld's.
   !!
   !! The populations are renormalised to the density's exact electron count,
   !! Tr(D S), as Hirshfeld's are, so the charges sum to the molecular charge
   !! whatever the grid missed of a diffuse density.
   !!
   !! A lone atom needs no fit: every weight is one, and its charge is Z less
   !! the integrated density. Elsewhere a shell the fit gives under 1e-4
   !! electrons is retired and two shells of one atom whose exponents meet are
   !! merged, after which the fit goes on with the rest.
   !!
   !! Refused: an effective core potential, whose missing core electrons the
   !! shell count does not know about.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_string_utils, only: int_to_text
   use mqc_physical_constants, only: PI
   use mqc_dft_grid, only: dft_grid_t, build_dft_grid
   use mqc_diis, only: diis_state_t
   use mqc_population_analysis, only: renormalise_populations, GRID_DEFICIT_WARNING
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_ao, only: eval_ao_block, eval_rho, shell_extents, block_significant_aos, &
                         AO_POINT_BLOCK
   implicit none
   private

   public :: mbis_model_t        !! The fitted pro-model
   public :: mbis_populations    !! Charges, and spin populations for an open shell
   public :: mbis_shell_count    !! Shells an element is given
   public :: MBIS_GRID_LEVEL

   integer, parameter :: MBIS_GRID_LEVEL = 5
      !! Quadrature level when the caller names none.
   integer, parameter :: MAX_ITER = 1000
   integer, parameter :: DIIS_VECTORS = 8
   real(dp), parameter :: CONVERGENCE = 1.0e-10_dp
      !! Largest change of any shell population, electrons, and of any exponent
      !! relative to itself, between two iterations.
   real(dp), parameter :: MIN_SHELL_POPULATION = 1.0e-4_dp
      !! A shell given fewer electrons than this is retired from the model.
   real(dp), parameter :: MERGE_TOLERANCE = 1.0e-4_dp
      !! Two shells of one atom whose exponents differ by less than this
      !! fraction of their mean are merged into one.
   ! Both are denspart's reduction, at its thresholds: a shell holding nothing
   ! has an exponent that is a ratio of two vanishing integrals, and two shells
   ! at one exponent are one exponential whose split is undetermined.
   real(dp), parameter :: SHELL_FLOOR = 1.0e-14_dp
      !! A shell is zero beyond the radius where it falls below this, e/Bohr^3.
      !! Where another shell is larger its share there is below this relative
      !! to it, and where every shell is this small the density being divided
      !! is too.
   real(dp), parameter :: AO_SCREEN = 1.0e-12_dp
      !! AO value below which a shell is dropped from a grid block.
   integer, parameter :: N_ROWS = 7
   integer, parameter :: ROW_END(N_ROWS) = [2, 10, 18, 36, 54, 86, 118]
      !! Last element of each row of the periodic table.
   integer, parameter :: SHELL_CAPACITY(N_ROWS) = [2, 8, 8, 18, 18, 32, 32]
      !! Electrons the starting model puts in each shell before the next.
   real(dp), parameter :: SAME_SHELL_SCREENING = 0.35_dp
      !! Slater's screening by another electron of the same shell.

   type :: mbis_model_t
      !! The fitted pro-model, shell by shell
      integer :: n_shells = 0
      integer, allocatable :: atom(:)            !! (n_shells) owning atom, 1-based
      real(dp), allocatable :: population(:)     !! (n_shells) N_Ai, electrons
      real(dp), allocatable :: exponent(:)       !! (n_shells) a_Ai, 1/Bohr
      integer :: iterations = 0
      logical :: converged = .false.
      integer :: retired = 0     !! Shells dropped for holding nothing
      integer :: merged = 0      !! Shells folded into a neighbour with their exponent
   contains
      procedure :: destroy => mbis_model_destroy
   end type mbis_model_t

contains

   pure function mbis_shell_count(atomic_number) result(n)
      !! One shell per row of the periodic table the element reaches
      integer, intent(in) :: atomic_number
      integer :: n

      do n = 1, size(ROW_END)
         if (atomic_number <= ROW_END(n)) return
      end do
      n = size(ROW_END)
   end function mbis_shell_count

   subroutine mbis_populations(mol, density, charges, error, spin_density, &
                               spin_populations, grid_level, model)
      !! q_A = Z_A - sum_i N_Ai, from a pro-model fitted to `density`
      !!
      !! `density` is the total, both spins. With `spin_density` (alpha minus
      !! beta) the converged weights give `spin_populations`. A ghost centre has
      !! no shells and a charge of zero. An unconverged fit is an error. For a
      !! lone atom `model` is the unfitted start, with `iterations` zero.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)              !! (nao, nao)
      real(dp), allocatable, intent(out) :: charges(:)   !! (natm), e
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: spin_density(:, :)
      real(dp), allocatable, intent(out), optional :: spin_populations(:)
      integer, intent(in), optional :: grid_level        !! Default `MBIS_GRID_LEVEL`
      type(mbis_model_t), intent(out), optional :: model

      type(dft_grid_t) :: grid
      type(mbis_model_t) :: fit
      type(diis_state_t) :: diis
      real(dp), allocatable :: rho(:), spin_rho(:), extents(:), moment(:), spins(:)
      real(dp), allocatable :: share_population(:), electrons(:), overlap(:, :)
      real(dp) :: deficit
      integer :: level, k
      logical :: want_spin

      if (error%has_error()) return
      want_spin = present(spin_density) .and. present(spin_populations)
      if (allocated(mol%core_electrons)) then
         if (any(mol%core_electrons > 0)) then
            call error%set(ERROR_VALIDATION, "MBIS charges: not available with an "// &
                           "effective core potential")
            return
         end if
      end if

      level = MBIS_GRID_LEVEL
      if (present(grid_level)) level = grid_level
      call build_dft_grid(mol%coords, mol%atomic_numbers, grid, error, level=level)
      if (error%has_error()) return
      call shell_extents(mol, AO_SCREEN, extents)

      ! The density does not change while the pro-model is fitted to it, so it
      ! is evaluated once and every iteration reads it back.
      call density_on_grid(mol, grid, extents, density, rho, error)
      if (error%has_error()) return
      if (want_spin) then
         call density_on_grid(mol, grid, extents, spin_density, spin_rho, error)
         if (error%has_error()) return
      end if

      call initial_model(mol, fit)
      ! With one atom to share between, every weight is one whatever the
      ! pro-atom, so its charge needs no fit -- and the fit would be the
      ! flattest there is: a bare Li+, an ordinary monomer in an ion-water
      ! cluster, has two shells for one 1s pair and crawled along the
      ! split between them to the iteration cap.
      fit%converged = count(mol%charges > 0.0_dp) <= 1
      allocate (moment(fit%n_shells), share_population(fit%n_shells))
      if (.not. fit%converged) call diis%init(DIIS_VECTORS, 2*fit%n_shells, 2*fit%n_shells)
      do k = 1, MAX_ITER
         if (fit%converged) exit
         call shell_moments(mol, grid, rho, fit, share_population, moment)
         fit%iterations = k
         call advance_model(fit, share_population, moment, diis, fit%converged)
         ! Reduced only at a converged model, as denspart reduces: shells can
         ! pass close to one another on the way to a minimum that keeps them
         ! apart, and merging them in passing would decide that minimum away.
         if (.not. fit%converged) cycle
         if (.not. reduce_model(fit, fit%population)) cycle
         ! The parameter vector changes length, so the DIIS history goes with
         ! it, and the reduced model is fitted again from where it stands.
         deallocate (moment, share_population)
         allocate (moment(fit%n_shells), share_population(fit%n_shells))
         call diis%destroy()
         if (fit%n_shells > 0) call diis%init(DIIS_VECTORS, 2*fit%n_shells, 2*fit%n_shells)
         fit%converged = fit%n_shells == 0
      end do
      call diis%destroy()
      if (.not. fit%converged) then
         call error%set(ERROR_VALIDATION, "MBIS charges: the pro-model did not converge in "// &
                        int_to_text(MAX_ITER)//" iterations")
         return
      end if

      call logger%verbose("  mbis: "//to_char(fit%n_shells)//" shells fitted in "// &
                          to_char(fit%iterations)//" iterations on "// &
                          to_char(grid%n_points)//" grid points")
      if (fit%retired + fit%merged > 0) then
         call logger%verbose("  mbis: "//to_char(fit%retired)//" empty shells retired, "// &
                             to_char(fit%merged)//" merged into a neighbour with "// &
                             "the same exponent")
      end if

      ! Each atom's electrons are its weighted share of the density. At a
      ! converged fit that is the sum of its shell populations; it is taken
      ! from the weights so the unfitted lone atom gets it the same way.
      call atom_shares(mol, grid, rho, fit, electrons)
      if (want_spin) call atom_shares(mol, grid, spin_rho, fit, spins)
      call grid%destroy()

      ! Summed to Tr(D S) rather than to what the grid integrated, as
      ! Hirshfeld's are: a diffuse anion's tail is what the grid misses.
      call mol%overlap(overlap)
      if (want_spin) then
         call renormalise_populations(electrons, sum(density*overlap), deficit, &
                                      spins, sum(spin_density*overlap))
         call move_alloc(spins, spin_populations)
      else
         call renormalise_populations(electrons, sum(density*overlap), deficit)
      end if
      call logger%verbose("  mbis: the grid missed "//to_char(deficit)// &
                          " electrons, shared out by population")
      if (abs(deficit) > GRID_DEFICIT_WARNING) then
         call logger%warning("  mbis: the integration grid missed "//to_char(deficit)// &
                             " electrons; the charges are renormalised, but a grid this "// &
                             "coarse for the density is worth a finer one")
      end if

      allocate (charges(mol%natm))
      where (mol%charges > 0.0_dp)
         charges = mol%charges - electrons
      else where
         charges = 0.0_dp
      end where
      if (present(model)) model = fit
   end subroutine mbis_populations

   subroutine initial_model(mol, fit)
      ! Each real atom's shells, filled in aufbau order -- 2, 8, 8, 18, ...
      ! electrons, the last shell taking the rest -- with exponents from a
      ! Slater-style screened charge, 2 (Z - sigma) / n: sigma counts every
      ! inner electron in full and 0.35 for each other electron of the shell.
      !
      ! The start decides which stationary point the fixed point reaches. An
      ! even split with exponents from 2Z down to 2 left lithium's valence
      ! shell too full and too tight in LiF, and the fit collapsed it onto
      ! the core: Li at +1.000 where the minimum denspart finds is +0.940.
      ! From this start the fit reaches denspart's.
      type(czt_molecule_t), intent(in) :: mol
      type(mbis_model_t), intent(out) :: fit

      integer :: ia, i, n, k, z
      real(dp) :: inner, held, screened

      fit%n_shells = 0
      do ia = 1, mol%natm
         if (mol%charges(ia) <= 0.0_dp) cycle
         fit%n_shells = fit%n_shells + mbis_shell_count(mol%atomic_numbers(ia))
      end do
      allocate (fit%atom(fit%n_shells), fit%population(fit%n_shells), &
                fit%exponent(fit%n_shells))

      k = 0
      do ia = 1, mol%natm
         if (mol%charges(ia) <= 0.0_dp) cycle
         z = mol%atomic_numbers(ia)
         n = mbis_shell_count(z)
         inner = 0.0_dp
         do i = 1, n
            k = k + 1
            if (i < n) then
               held = min(real(SHELL_CAPACITY(i), dp), real(z, dp) - inner)
            else
               held = real(z, dp) - inner
            end if
            screened = inner + SAME_SHELL_SCREENING*max(held - 1.0_dp, 0.0_dp)
            fit%atom(k) = ia
            fit%population(k) = held
            fit%exponent(k) = 2.0_dp*(real(z, dp) - screened)/real(i, dp)
            inner = inner + held
         end do
         ! Electrons, not nuclear charge: a ghost has none and is skipped
         ! above, and the populations must sum to what the atom holds.
         fit%population(k - n + 1:k) = fit%population(k - n + 1:k)*mol%charges(ia)/real(z, dp)
      end do
   end subroutine initial_model

   function reduce_model(fit, share_population) result(reduced)
      ! Retire every shell given under `MIN_SHELL_POPULATION` electrons, and
      ! merge each shell into the previous one of the same atom when their
      ! exponents agree to `MERGE_TOLERANCE`: population summed, exponent the
      ! population-weighted mean. True when anything changed. An atom's
      ! shells are contiguous, so only neighbours need comparing.
      type(mbis_model_t), intent(inout) :: fit
      real(dp), intent(in) :: share_population(:)
      logical :: reduced

      logical, allocatable :: keep(:)
      integer :: k, prev
      real(dp) :: total

      allocate (keep(fit%n_shells))
      keep = share_population >= MIN_SHELL_POPULATION
      fit%retired = fit%retired + count(.not. keep)

      prev = 0
      do k = 1, fit%n_shells
         if (.not. keep(k)) cycle
         if (prev > 0) then
            if (fit%atom(prev) == fit%atom(k) .and. &
                abs(fit%exponent(k) - fit%exponent(prev)) < &
                MERGE_TOLERANCE*0.5_dp*(fit%exponent(k) + fit%exponent(prev))) then
               total = fit%population(prev) + fit%population(k)
               fit%exponent(prev) = (fit%population(prev)*fit%exponent(prev) + &
                                     fit%population(k)*fit%exponent(k))/total
               fit%population(prev) = total
               keep(k) = .false.
               fit%merged = fit%merged + 1
               cycle
            end if
         end if
         prev = k
      end do

      reduced = .not. all(keep)
      if (.not. reduced) return
      fit%atom = pack(fit%atom, keep)
      fit%population = pack(fit%population, keep)
      fit%exponent = pack(fit%exponent, keep)
      fit%n_shells = size(fit%atom)
   end function reduce_model

   subroutine mbis_model_destroy(this)
      !! Release the shells and reset to the unfitted state
      class(mbis_model_t), intent(inout) :: this

      if (allocated(this%atom)) deallocate (this%atom)
      if (allocated(this%population)) deallocate (this%population)
      if (allocated(this%exponent)) deallocate (this%exponent)
      this%n_shells = 0
      this%iterations = 0
      this%retired = 0
      this%merged = 0
      this%converged = .false.
   end subroutine mbis_model_destroy

   subroutine advance_model(fit, share_population, moment, diis, converged)
      ! One step of the fixed point -- each shell's population becomes its
      ! share and its exponent 3 N / <r> -- extrapolated by DIIS, which turns
      ! the map's linear convergence into a few dozen iterations: 154 plain
      ! against 37 for a 46-atom peptide, 72 against 19 for water. The
      ! parameters are the populations and the *logarithms* of the exponents: the residual is then dimensionless
      ! throughout, a change in ln a being a relative change in a, and an
      ! extrapolated exponent cannot come out negative. A population can, and
      ! a step that makes one is replaced by the plain update.
      type(mbis_model_t), intent(inout) :: fit
      real(dp), intent(in) :: share_population(:)
      real(dp), intent(in) :: moment(:)
      type(diis_state_t), intent(inout) :: diis
      logical, intent(out) :: converged

      real(dp), allocatable :: mapped(:), plain(:), residual(:)
      integer :: n
      logical :: ok

      n = fit%n_shells
      allocate (mapped(2*n), residual(2*n))
      mapped(1:n) = share_population
      mapped(n + 1:2*n) = log(3.0_dp*share_population/moment)
      residual(1:n) = mapped(1:n) - fit%population
      residual(n + 1:2*n) = mapped(n + 1:2*n) - log(fit%exponent)
      converged = maxval(abs(residual)) < CONVERGENCE

      if (.not. converged) then
         plain = mapped
         call diis%push(mapped, residual)
         call diis%extrapolate(mapped, ok)
         if (any(mapped(1:n) <= 0.0_dp)) then
            mapped = plain
            ok = .false.
         end if
         if (.not. ok .and. diis%count() >= 2) then
            ! A singular or unphysical subspace is restarted from the newest
            ! entry rather than kept, as the spherical-atom SCF does.
            call diis%destroy()
            call diis%init(DIIS_VECTORS, 2*n, 2*n)
            call diis%push(plain, residual)
         end if
      end if

      fit%population = mapped(1:n)
      fit%exponent = exp(mapped(n + 1:2*n))
   end subroutine advance_model

   subroutine shell_moments(mol, grid, rho, fit, pop, moment)
      ! One pass of the fixed point: each shell's share of the density and the
      ! first radial moment of that share, under the current model.
      type(czt_molecule_t), intent(in) :: mol
      type(dft_grid_t), intent(in) :: grid
      real(dp), intent(in) :: rho(:)
      type(mbis_model_t), intent(in) :: fit
      real(dp), intent(out) :: pop(:)
      real(dp), intent(out) :: moment(:)

      real(dp), allocatable :: share(:, :), dist(:, :), local_n(:), local_m(:)
      integer, allocatable :: near(:)
      integer :: g0, g1, n_near, i

      pop = 0.0_dp
      moment = 0.0_dp

      !$omp parallel default(shared) &
      !$omp    private(g0, g1, n_near, i, share, dist, near, local_n, local_m)
      allocate (near(fit%n_shells), local_n(fit%n_shells), local_m(fit%n_shells))
      local_n = 0.0_dp
      local_m = 0.0_dp
      !$omp do schedule(dynamic)
      do g0 = 1, grid%n_points, AO_POINT_BLOCK
         g1 = min(g0 + AO_POINT_BLOCK - 1, grid%n_points)
         call shell_shares(mol, grid%coords(:, g0:g1), fit, near, n_near, share, dist)
         do i = 1, n_near
            local_n(near(i)) = local_n(near(i)) + &
                               sum(grid%weights(g0:g1)*rho(g0:g1)*share(:, i))
            local_m(near(i)) = local_m(near(i)) + &
                               sum(grid%weights(g0:g1)*rho(g0:g1)*share(:, i)*dist(:, i))
         end do
      end do
      !$omp end do
      !$omp critical (mqc_mbis_sum)
      pop = pop + local_n
      moment = moment + local_m
      !$omp end critical (mqc_mbis_sum)
      !$omp end parallel
   end subroutine shell_moments

   subroutine atom_shares(mol, grid, values, fit, per_atom)
      ! The converged weights applied to another function on the grid, summed
      ! over each atom's shells.
      type(czt_molecule_t), intent(in) :: mol
      type(dft_grid_t), intent(in) :: grid
      real(dp), intent(in) :: values(:)
      type(mbis_model_t), intent(in) :: fit
      real(dp), allocatable, intent(out) :: per_atom(:)

      real(dp), allocatable :: share(:, :), dist(:, :), local_a(:)
      integer, allocatable :: near(:)
      integer :: g0, g1, n_near, i

      allocate (per_atom(mol%natm))
      per_atom = 0.0_dp

      !$omp parallel default(shared) private(g0, g1, n_near, i, share, dist, near, local_a)
      allocate (near(fit%n_shells), local_a(mol%natm))
      local_a = 0.0_dp
      !$omp do schedule(dynamic)
      do g0 = 1, grid%n_points, AO_POINT_BLOCK
         g1 = min(g0 + AO_POINT_BLOCK - 1, grid%n_points)
         call shell_shares(mol, grid%coords(:, g0:g1), fit, near, n_near, share, dist)
         do i = 1, n_near
            local_a(fit%atom(near(i))) = local_a(fit%atom(near(i))) + &
                                         sum(grid%weights(g0:g1)*values(g0:g1)*share(:, i))
         end do
      end do
      !$omp end do
      !$omp critical (mqc_mbis_atoms)
      per_atom = per_atom + local_a
      !$omp end critical (mqc_mbis_atoms)
      !$omp end parallel
   end subroutine atom_shares

   subroutine shell_shares(mol, coords, fit, near, n_near, share, dist)
      ! Each shell's share of every point in a block, for the shells that reach
      ! it, and the distance from the shell's atom. Zero at a point no shell
      ! reaches.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)
      type(mbis_model_t), intent(in) :: fit
      integer, intent(out) :: near(:)
      integer, intent(out) :: n_near
      real(dp), allocatable, intent(out) :: share(:, :)   ! (n_points, n_near)
      real(dp), allocatable, intent(out) :: dist(:, :)    ! (n_points, n_near)

      real(dp), allocatable :: total(:)
      real(dp) :: centre(3)
      real(dp) :: radius, cutoff, prefactor
      integer :: k, ig, i, np

      np = size(coords, 2)
      centre = sum(coords, dim=2)/real(np, dp)
      radius = 0.0_dp
      do ig = 1, np
         radius = max(radius, norm2(coords(:, ig) - centre))
      end do

      n_near = 0
      do k = 1, fit%n_shells
         cutoff = shell_cutoff(fit%population(k), fit%exponent(k))
         if (norm2(mol%coords(:, fit%atom(k)) - centre) > cutoff + radius) cycle
         n_near = n_near + 1
         near(n_near) = k
      end do

      allocate (share(np, n_near), dist(np, n_near), total(np))
      total = 0.0_dp
      do i = 1, n_near
         k = near(i)
         ! An atom's shells are contiguous, so its distances are worked out
         ! for the first of them and copied for the rest.
         if (i > 1) then
            if (fit%atom(near(i - 1)) == fit%atom(k)) then
               dist(:, i) = dist(:, i - 1)
            else
               call distances(coords, mol%coords(:, fit%atom(k)), dist(:, i))
            end if
         else
            call distances(coords, mol%coords(:, fit%atom(k)), dist(:, i))
         end if
         prefactor = fit%population(k)*fit%exponent(k)**3/(8.0_dp*PI)
         share(:, i) = prefactor*exp(-fit%exponent(k)*dist(:, i))
         total = total + share(:, i)
      end do
      where (total > 0.0_dp)
         total = 1.0_dp/total
      else where
         total = 0.0_dp
      end where
      do i = 1, n_near
         share(:, i) = share(:, i)*total
      end do
   end subroutine shell_shares

   pure subroutine distances(coords, centre, dist)
      ! |r_g - centre| for every point of a block.
      real(dp), intent(in) :: coords(:, :)
      real(dp), intent(in) :: centre(3)
      real(dp), intent(out) :: dist(:)

      integer :: ig

      do ig = 1, size(coords, 2)
         dist(ig) = sqrt((coords(1, ig) - centre(1))**2 + (coords(2, ig) - centre(2))**2 + &
                         (coords(3, ig) - centre(3))**2)
      end do
   end subroutine distances

   pure function shell_cutoff(population, exponent) result(r)
      ! Where N a^3/(8 pi) exp(-a r) falls to `SHELL_FLOOR`.
      real(dp), intent(in) :: population, exponent
      real(dp) :: r

      r = max(0.0_dp, log(population*exponent**3/(8.0_dp*PI*SHELL_FLOOR))/exponent)
   end function shell_cutoff

   subroutine density_on_grid(mol, grid, extents, density, rho, error)
      ! rho(r_g) at every grid point, block by block, screened as the
      ! exchange-correlation quadrature screens.
      type(czt_molecule_t), intent(in) :: mol
      type(dft_grid_t), intent(in) :: grid
      real(dp), intent(in) :: extents(:)
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: rho(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: ao(:, :), block_rho(:), d_sig(:, :)
      integer, allocatable :: ao_list(:), ao_offset(:)
      logical, allocatable :: shell_mask(:)
      type(error_t) :: local_error
      integer :: g0, g1, n_sig, ia, ja
      logical :: failed

      allocate (rho(grid%n_points))
      rho = 0.0_dp
      failed = .false.

      ! `local_error` firstprivate rather than private, for the reason
      ! `xc_add_potential` gives: a private copy need not be default-initialised.
      !$omp parallel default(shared) &
      !$omp    private(g0, g1, n_sig, ia, ja, ao, block_rho, d_sig, ao_list, ao_offset, &
      !$omp            shell_mask) &
      !$omp    firstprivate(local_error)
      allocate (d_sig(mol%nao, mol%nao), shell_mask(mol%nbas), ao_offset(mol%nbas), &
                ao_list(mol%nao))
      !$omp do schedule(dynamic)
      do g0 = 1, grid%n_points, AO_POINT_BLOCK
         if (failed) cycle
         g1 = min(g0 + AO_POINT_BLOCK - 1, grid%n_points)
         call block_significant_aos(mol, grid%coords(:, g0:g1), extents, shell_mask, &
                                    ao_list, ao_offset, n_sig)
         if (n_sig == 0) cycle
         do ja = 1, n_sig
            do ia = 1, n_sig
               d_sig(ia, ja) = density(ao_list(ia), ao_list(ja))
            end do
         end do
         call eval_ao_block(mol, grid%coords(:, g0:g1), ao, local_error, &
                            shell_mask=shell_mask, ao_offset=ao_offset, n_ao_out=n_sig)
         if (local_error%has_error()) then
            !$omp critical (mqc_mbis_error)
            if (.not. failed) error = local_error
            failed = .true.
            !$omp end critical (mqc_mbis_error)
            cycle
         end if
         call eval_rho(ao, d_sig(1:n_sig, 1:n_sig), block_rho)
         rho(g0:g1) = block_rho
      end do
      !$omp end do
      !$omp end parallel
   end subroutine density_on_grid

end module mqc_czt_mbis
