!! Chain-of-spheres exchange against the analytic exchange matrix
module test_mqc_czt_cosx
   !! COSX is a quadrature, so nothing it produces equals the analytic exchange
   !! matrix and a tolerance alone cannot tell a correct implementation on a
   !! modest grid from a subtly wrong one on any grid. Three things are checked
   !! that between them can:
   !!
   !!   * **The parts that are exact are held to rounding.** The contraction
   !!     with the potential integral is compared against the dense integral
   !!     tensor from `esp_matrices`, and the batching against itself at a
   !!     different batch size. Neither involves a grid error.
   !!   * **The error falls as the grid grows.** An approximation that works
   !!     converges; a wrong factor, a transposed index or a misplaced weight
   !!     leaves an error that no grid removes. The finest grid is then held to
   !!     a bound far below anything a bug would leave.
   !!   * **The reference is the direct four-index build**, on the same
   !!     density, with density screening off so it is the exchange matrix and
   !!     not an SCF's approximation to it.
   !!
   !! The density is a converged RHF one. Any symmetric matrix would do for the
   !! algebra, but the size of the quadrature error is only meaningful for a
   !! density a calculation would actually hand over.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_calculation_defaults, only: DEFAULT_COSX_GRID
   use mqc_czt_direct, only: schwarz_bounds, build_fock_direct, build_fock_direct_nosym, &
                             direct_stats_t
   use mqc_czt_esp, only: esp_matrices
   use mqc_czt_ao, only: eval_ao_block
   use mqc_czt_cosx, only: cosx_exchange, cosx_exchange_trace, cosx_batch_exchange, &
                           cosx_esp_contract, &
                           cosx_density_on_grid, cosx_accumulate, cosx_grid_t, &
                           build_cosx_grid, cosx_overlap_fit, cosx_numerical_overlap, &
                           cosx_radial_points, cosx_angular_points, &
                           COSX_ROLE_SCF, COSX_ROLE_FINAL, &
                           COSX_SCREEN_ITERATION, COSX_SCREEN_FINAL
   use mqc_dft_grid, only: dft_grid_t, build_dft_grid
   implicit none
   private

   public :: collect_mqc_czt_cosx_tests

   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [3, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]

contains

   subroutine collect_mqc_czt_cosx_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("cosx_contraction_matches_the_dense_integral", test_contract), &
                  new_unittest("cosx_steps_compose_to_the_batch", test_steps), &
                  new_unittest("cosx_does_not_depend_on_the_batch_size", test_batching), &
                  new_unittest("cosx_converges_to_the_analytic_exchange", test_convergence), &
                  new_unittest("cosx_handles_a_density_of_no_symmetry", test_nosym), &
                  new_unittest("cosx_in_a_spherical_basis", test_spherical), &
                  new_unittest("cosx_overlap_fit_makes_the_overlap_exact", test_fit_overlap), &
                  new_unittest("cosx_overlap_fit_improves_a_small_grid", test_fit_exchange), &
                  new_unittest("cosx_grid_levels_are_ordered", test_grid_levels), &
                  new_unittest("cosx_screening_drops_only_what_is_negligible", test_screening), &
                  new_unittest("cosx_batch_integral_gives_the_same_exchange", test_batch_kernel), &
                  new_unittest("cosx_trace_is_the_trace_of_the_exchange", test_trace), &
                  new_unittest("cosx_refuses_malformed_input", test_refusals) &
                  ]
   end subroutine collect_mqc_czt_cosx_tests

   subroutine water_density(basis, mol, density, error)
      !! Water and its converged RHF density in the basis named
      character(len=*), intent(in) :: basis
      type(czt_molecule_t), intent(out) :: mol
      real(dp), allocatable, intent(out) :: density(:, :)
      type(error_type), allocatable, intent(out) :: error

      type(rhf_result_t) :: scf
      type(error_t) :: err

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, basis, mol, err)
      call check(error,.not. err%has_error(), "the molecule should build: "//err%get_full_trace())
      if (allocated(error)) return
      call run_czt_rhf(mol, 10, 60, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
      call check(error,.not. err%has_error(), "the SCF should run: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, scf%converged, "the SCF should converge")
      if (allocated(error)) return
      density = scf%density
   end subroutine water_density

   subroutine analytic_exchange(mol, density, k, error)
      !! `K(mu, nu) = sum P(kappa, lambda) (mu kappa | nu lambda)` from the direct build
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: k(:, :)
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: bounds(:, :), zero(:, :)
      type(direct_stats_t) :: stats
      type(error_t) :: err

      call schwarz_bounds(mol, bounds, err)
      call check(error,.not. err%has_error(), "Schwarz bounds failed: "//err%get_full_trace())
      if (allocated(error)) return
      allocate (zero(mol%nao, mol%nao), k(mol%nao, mol%nao))
      zero = 0.0_dp
      ! With no core Hamiltonian and no Coulomb term the "Fock matrix" is -K/2.
      call build_fock_direct(mol, zero, density, bounds, k, stats, err, &
                             screen_tol=1.0e-14_dp, k_scale=1.0_dp, j_scale=0.0_dp, &
                             density_screen=.false.)
      call check(error,.not. err%has_error(), "the direct build failed: "//err%get_full_trace())
      if (allocated(error)) return
      k = -2.0_dp*k
   end subroutine analytic_exchange

   subroutine test_contract(error)
      !! `cosx_esp_contract` against the dense tensor, point by point
      !!
      !! `F` is arbitrary and different at every point, so a contraction that
      !! paired point `g` of `F` with the integral at another point, or that
      !! dropped the mirror half of a shell pair, cannot pass.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N_PTS = 5
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: dense(:, :, :), f(:, :), g(:, :)
      real(dp) :: points(3, N_PTS)
      real(dp) :: expect, worst, scale
      integer :: ig, mu, nu

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "6-31g*", mol, err)
      call check(error,.not. err%has_error(), "the molecule should build: "//err%get_full_trace())
      if (allocated(error)) return

      points = reshape([0.7_dp, 0.4_dp, 1.3_dp, &
                        -1.1_dp, 0.9_dp, -0.6_dp, &
                        0.0_dp, 0.1_dp, 0.2_dp, &
                        2.5_dp, -1.7_dp, 0.3_dp, &
                        0.0_dp, -1.4_dp, 1.0_dp], [3, N_PTS])

      allocate (f(N_PTS, mol%nao))
      do mu = 1, mol%nao
         do ig = 1, N_PTS
            f(ig, mu) = sin(1.7_dp*real(ig, dp) + 0.37_dp*real(mu, dp))
         end do
      end do

      call esp_matrices(mol, points, dense, err)
      call check(error,.not. err%has_error(), "esp_matrices failed: "//err%get_full_trace())
      if (allocated(error)) return
      call cosx_esp_contract(mol, points, f, g, err)
      call check(error,.not. err%has_error(), "the contraction failed: "//err%get_full_trace())
      if (allocated(error)) return

      worst = 0.0_dp
      scale = 0.0_dp
      do nu = 1, mol%nao
         do ig = 1, N_PTS
            expect = dot_product(dense(nu, :, ig), f(ig, :))
            worst = max(worst, abs(g(ig, nu) - expect))
            scale = max(scale, abs(expect))
         end do
      end do
      call check(error, worst < 1.0e-12_dp, &
                 "the contraction differs from the dense potential integral")
      if (allocated(error)) return
      call check(error, scale > 0.1_dp, "the contracted values are too small to test")
      call mol%destroy()
   end subroutine test_contract

   subroutine test_steps(error)
      !! The three public steps, called by hand, are what the batch routine does
      !!
      !! The steps are public so another backend can replace one of them; this
      !! pins that the batch routine is their composition and nothing more.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(dft_grid_t) :: grid
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), ao(:, :), f(:, :), g(:, :)
      real(dp), allocatable :: by_hand(:, :), batch(:, :)
      integer, parameter :: N_USE = 40

      call water_density("sto-3g", mol, density, error)
      if (allocated(error)) return
      call build_dft_grid(mol%coords, mol%atomic_numbers, grid, err, level=0)
      call check(error,.not. err%has_error(), "the grid should build: "//err%get_full_trace())
      if (allocated(error)) return

      allocate (by_hand(mol%nao, mol%nao), batch(mol%nao, mol%nao))
      by_hand = 0.0_dp
      batch = 0.0_dp

      call eval_ao_block(mol, grid%coords(:, 1:N_USE), ao, err)
      call cosx_density_on_grid(ao, grid%weights(1:N_USE), density, f)
      call cosx_esp_contract(mol, grid%coords(:, 1:N_USE), f, g, err)
      call cosx_accumulate(ao, g, by_hand)
      call check(error,.not. err%has_error(), "a step failed: "//err%get_full_trace())
      if (allocated(error)) return

      call cosx_batch_exchange(mol, grid%coords(:, 1:N_USE), grid%weights(1:N_USE), &
                               density, batch, err)
      call check(error,.not. err%has_error(), "the batch failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, maxval(abs(batch - by_hand)) < 1.0e-14_dp, &
                 "the batch routine is not the composition of its steps")
      if (allocated(error)) return
      call check(error, maxval(abs(batch)) > 1.0e-6_dp, "the batch contributed nothing")
      call mol%destroy()
   end subroutine test_steps

   subroutine test_batching(error)
      !! Batches of 7 and one batch of everything sum to the same matrix
      !!
      !! Seven divides nothing convenient, so the last batch is short and the
      !! boundaries fall mid-atom.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(dft_grid_t) :: grid
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), small(:, :), whole(:, :)

      call water_density("sto-3g", mol, density, error)
      if (allocated(error)) return
      call build_dft_grid(mol%coords, mol%atomic_numbers, grid, err, level=0)
      call check(error,.not. err%has_error(), "the grid should build: "//err%get_full_trace())
      if (allocated(error)) return

      call cosx_exchange(mol, grid%coords, grid%weights, density, small, err, batch_size=7)
      call cosx_exchange(mol, grid%coords, grid%weights, density, whole, err, &
                         batch_size=grid%n_points)
      call check(error,.not. err%has_error(), "COSX failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, maxval(abs(small - whole)) < 1.0e-11_dp, &
                 "the exchange matrix depends on the batch size")
      call mol%destroy()
   end subroutine test_batching

   subroutine test_convergence(error)
      !! The error against the analytic matrix falls with the grid, to a small bound
      !!
      !! Water in 6-31G*, grid levels 1, 4 and 9. Measured when this was
      !! written, on 10124, 59676 and 489828 points: the largest element of the
      !! difference went 8.1e-5, 6.5e-8, 3.1e-10 and the exchange energy error
      !! 1.8e-5, 4.5e-8, 1.5e-10 hartree. The bounds below sit well above the
      !! last of those and far below the first. A mistake in the
      !! algebra leaves an error of the order of the matrix itself, which is
      !! about 1, on every grid.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: LEVELS(3) = [1, 4, 9]
      type(czt_molecule_t) :: mol
      type(dft_grid_t) :: grid
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), exact(:, :), k(:, :)
      real(dp) :: worst(3), e_err(3)
      real(dp) :: e_exact, asym
      integer :: il

      call water_density("6-31g*", mol, density, error)
      if (allocated(error)) return
      call check(error, mol%cartesian, "6-31G* was expected to be Cartesian")
      if (allocated(error)) return
      call analytic_exchange(mol, density, exact, error)
      if (allocated(error)) return
      e_exact = -0.25_dp*sum(density*exact)

      do il = 1, 3
         call build_dft_grid(mol%coords, mol%atomic_numbers, grid, err, level=LEVELS(il))
         call check(error,.not. err%has_error(), "the grid should build: "//err%get_full_trace())
         if (allocated(error)) return
         call cosx_exchange(mol, grid%coords, grid%weights, density, k, err)
         call check(error,.not. err%has_error(), "COSX failed: "//err%get_full_trace())
         if (allocated(error)) return
         worst(il) = maxval(abs(k - exact))
         e_err(il) = abs(-0.25_dp*sum(density*k) - e_exact)
         asym = maxval(abs(k - transpose(k)))
         call check(error, asym < 1.0e-14_dp, "the symmetrised matrix is not symmetric")
         if (allocated(error)) return
         call grid%destroy()
      end do

      call check(error, worst(2) < worst(1) .and. worst(3) < worst(2), &
                 "the exchange matrix error does not fall as the grid grows")
      if (allocated(error)) return
      call check(error, e_err(2) < e_err(1) .and. e_err(3) < e_err(2), &
                 "the exchange energy error does not fall as the grid grows")
      if (allocated(error)) return
      call check(error, worst(3) < 1.0e-8_dp, &
                 "the exchange matrix on the finest grid is too far from the analytic one")
      if (allocated(error)) return
      call check(error, e_err(3) < 1.0e-8_dp, &
                 "the exchange energy on the finest grid is too far from the analytic one")
      if (allocated(error)) return
      ! Guards the guard: the exchange energy of water is about -9 hartree.
      call check(error, abs(e_exact) > 5.0_dp, "the reference exchange energy is implausible")
      call mol%destroy()
   end subroutine test_convergence

   subroutine test_nosym(error)
      !! A density with no symmetry, unsymmetrised, against the permutation-explicit build
      !!
      !! A response density is not symmetric and neither is its exchange
      !! matrix, so this is the case that fixes which index of the density
      !! meets the grid: with the density transposed the result is the
      !! transpose, which for a symmetric density is invisible.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(dft_grid_t) :: grid
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: density(:, :), skew(:, :, :), bounds(:, :), zero(:, :)
      real(dp), allocatable :: exact(:, :, :), k(:, :)
      real(dp) :: diff, diff_transposed
      integer :: mu, nu

      call water_density("6-31g*", mol, density, error)
      if (allocated(error)) return

      ! The SCF density plus an antisymmetric part of comparable size.
      allocate (skew(mol%nao, mol%nao, 1))
      do nu = 1, mol%nao
         do mu = 1, mol%nao
            skew(mu, nu, 1) = density(mu, nu) &
                              + 0.3_dp*sin(real(mu - nu, dp))*sqrt(abs(density(mu, mu)*density(nu, nu)))
         end do
      end do

      call schwarz_bounds(mol, bounds, err)
      allocate (zero(mol%nao, mol%nao))
      zero = 0.0_dp
      call build_fock_direct_nosym(mol, zero, skew, bounds, exact, stats, err, &
                                   screen_tol=1.0e-14_dp, k_scale=1.0_dp, j_scale=0.0_dp)
      call check(error,.not. err%has_error(), "the direct build failed: "//err%get_full_trace())
      if (allocated(error)) return
      exact = -2.0_dp*exact

      call build_dft_grid(mol%coords, mol%atomic_numbers, grid, err, level=9)
      call cosx_exchange(mol, grid%coords, grid%weights, skew(:, :, 1), k, err, &
                         symmetrize=.false.)
      call check(error,.not. err%has_error(), "COSX failed: "//err%get_full_trace())
      if (allocated(error)) return

      diff = maxval(abs(k - exact(:, :, 1)))
      diff_transposed = maxval(abs(transpose(k) - exact(:, :, 1)))
      ! Measured 3.0e-10 on this grid; the transpose is 0.94 away.
      call check(error, diff < 1.0e-8_dp, &
                 "the exchange matrix of an unsymmetric density is wrong")
      if (allocated(error)) return
      ! Guards the guard: the transpose must be visibly different, or the
      ! density was not unsymmetric enough to tell the two apart.
      call check(error, diff_transposed > 1.0e-2_dp, &
                 "the test density does not distinguish K from its transpose")
      call mol%destroy()
   end subroutine test_nosym

   subroutine test_spherical(error)
      !! A spherical d shell goes through the spherical integral
      !!
      !! 6-31G*, used above, is Cartesian in this repository's basis file and
      !! cc-pVDZ is spherical, so this is the spherical half. The form is
      !! asserted, so a basis file that changed would fail here and not silently
      !! test the same path twice. Measured 1.4e-8 on the level 6 grid.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(dft_grid_t) :: grid
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), exact(:, :), k(:, :)

      call water_density("cc-pvdz", mol, density, error)
      if (allocated(error)) return
      call analytic_exchange(mol, density, exact, error)
      if (allocated(error)) return
      call build_dft_grid(mol%coords, mol%atomic_numbers, grid, err, level=6)
      call cosx_exchange(mol, grid%coords, grid%weights, density, k, err)
      call check(error,.not. err%has_error(), "COSX failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error,.not. mol%cartesian, "cc-pVDZ was expected to be spherical")
      if (allocated(error)) return
      call check(error, maxval(abs(k - exact)) < 1.0e-6_dp, &
                 "the exchange matrix in the spherical form is wrong")
      call mol%destroy()
   end subroutine test_spherical

   subroutine test_fit_overlap(error)
      !! The fit turns the grid's overlap into the analytic one, and is the identity on a fine grid
      !!
      !! `fit S_grid = S` is the definition, so it holds to rounding on any
      !! grid, including one too coarse to be used. On a fine grid the two
      !! overlaps already agree and the fit has nothing to do.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(cosx_grid_t) :: coarse
      type(dft_grid_t) :: fine
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), s_grid(:, :), fit(:, :), fitted(:, :)
      real(dp) :: raw, off_identity
      integer :: n_dropped, mu

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "6-31g*", mol, err)
      call check(error,.not. err%has_error(), "the molecule should build: "//err%get_full_trace())
      if (allocated(error)) return
      call mol%overlap(s)

      call build_cosx_grid(mol, coarse, err, level=1, role=COSX_ROLE_SCF)
      call cosx_numerical_overlap(mol, coarse%coords, coarse%weights, s_grid, err)
      call cosx_overlap_fit(mol, coarse%coords, coarse%weights, s, fit, err, n_dropped=n_dropped)
      call check(error,.not. err%has_error(), "the fit failed: "//err%get_full_trace())
      if (allocated(error)) return

      raw = maxval(abs(s_grid - s))
      fitted = matmul(fit, s_grid)
      call check(error, maxval(abs(fitted - s)) < 1.0e-10_dp, &
                 "the fitted overlap is not the analytic one")
      if (allocated(error)) return
      call check(error, n_dropped == 0, "the grid left a direction of this basis unresolved")
      if (allocated(error)) return
      ! Guards the guard: on a grid that already integrated the overlap the
      ! check above would pass with the fit doing nothing.
      call check(error, raw > 1.0e-6_dp, "the coarse grid is too good to test the fit on")
      if (allocated(error)) return

      call build_dft_grid(mol%coords, mol%atomic_numbers, fine, err, level=9)
      call cosx_overlap_fit(mol, fine%coords, fine%weights, s, fit, err)
      call check(error,.not. err%has_error(), "the fit failed: "//err%get_full_trace())
      if (allocated(error)) return
      off_identity = 0.0_dp
      do mu = 1, mol%nao
         fit(mu, mu) = fit(mu, mu) - 1.0_dp
      end do
      off_identity = maxval(abs(fit))
      call check(error, off_identity < 1.0e-7_dp, &
                 "the fit is not the identity on a grid that needs none")
      call mol%destroy()
   end subroutine test_fit_overlap

   subroutine test_fit_exchange(error)
      !! On the grids the SCF uses, the fitted exchange is closer to the analytic one
      !!
      !! This is what the fit is for. Water in 6-31G* on 25 radial shells and
      !! 50 angular points, sized here and not from the tables so the test does
      !! not move when they are retuned: measured when this was written, the exchange energy is off by
      !! 5.0e-4 hartree as the quadrature stands and 7.6e-5 fitted. A fit
      !! applied to the wrong index, or transposed, makes it worse and not
      !! better.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(cosx_grid_t) :: grid
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), exact(:, :), s(:, :), k(:, :)
      real(dp) :: e_exact, plain, fitted

      call water_density("6-31g*", mol, density, error)
      if (allocated(error)) return
      call analytic_exchange(mol, density, exact, error)
      if (allocated(error)) return
      e_exact = -0.25_dp*sum(density*exact)
      call mol%overlap(s)

      call build_cosx_grid(mol, grid, err, n_radial=25, n_angular=50, overlap=s)
      call check(error,.not. err%has_error(), "the grid should build: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, allocated(grid%fit), "a grid built with the overlap should carry the fit")
      if (allocated(error)) return

      call cosx_exchange(mol, grid%coords, grid%weights, density, k, err)
      plain = abs(-0.25_dp*sum(density*k) - e_exact)
      call cosx_exchange(mol, grid%coords, grid%weights, density, k, err, fit=grid%fit)
      call check(error,.not. err%has_error(), "COSX failed: "//err%get_full_trace())
      if (allocated(error)) return
      fitted = abs(-0.25_dp*sum(density*k) - e_exact)

      call check(error, fitted < 0.5_dp*plain, &
                 "the overlap fit did not improve the exchange energy on a small grid")
      if (allocated(error)) return
      call check(error, fitted < 3.0e-4_dp, &
                 "the fitted exchange energy on the smallest grid is too far off")
      if (allocated(error)) return
      call check(error, maxval(abs(k - transpose(k))) < 1.0e-14_dp, &
                 "the fitted exchange matrix is not symmetric")
      call mol%destroy()
   end subroutine test_fit_exchange

   subroutine test_grid_levels(error)
      !! The grid tables: larger with the level, with the role and with the period
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(cosx_grid_t) :: grid
      type(error_t) :: err
      integer :: level, n_scf, n_final, last_final

      ! Hydrogen, oxygen, sulfur: one more period, more radial shells.
      call check(error, cosx_radial_points(1, 2, COSX_ROLE_SCF) < &
                 cosx_radial_points(8, 2, COSX_ROLE_SCF) .and. &
                 cosx_radial_points(8, 2, COSX_ROLE_SCF) < &
                 cosx_radial_points(16, 2, COSX_ROLE_SCF), &
                 "the radial count does not grow with the period")
      if (allocated(error)) return

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return

      last_final = 0
      do level = 1, 3
         call check(error, cosx_angular_points(level, COSX_ROLE_FINAL) > &
                    cosx_angular_points(level, COSX_ROLE_SCF), &
                    "the final grid has no more angular points than the SCF grid")
         if (allocated(error)) return
         call build_cosx_grid(mol, grid, err, level=level, role=COSX_ROLE_SCF)
         n_scf = grid%n_points
         call build_cosx_grid(mol, grid, err, level=level, role=COSX_ROLE_FINAL)
         n_final = grid%n_points
         call check(error,.not. err%has_error(), "the grid should build: "//err%get_full_trace())
         if (allocated(error)) return
         call check(error, n_final > n_scf, "the final grid is not larger than the SCF grid")
         if (allocated(error)) return
         call check(error, n_final > last_final, "the final grid does not grow with the level")
         if (allocated(error)) return
         call check(error,.not. allocated(grid%fit), &
                    "a grid built without the overlap should carry no fit")
         if (allocated(error)) return
         last_final = n_final
      end do

      call build_cosx_grid(mol, grid, err, level=4)
      call check(error, err%has_error(), "a grid level that does not exist was accepted")
      if (allocated(error)) return
      call err%clear()
      call build_cosx_grid(mol, grid, err, level=0)
      call check(error, err%has_error(), "grid level zero was accepted")
      call mol%destroy()
   end subroutine test_grid_levels

   subroutine test_screening(error)
      !! Screening skips shell pairs, and what it skips does not change the answer
      !!
      !! Two waters 40 bohr apart, so the pairs with one shell on each have
      !! Schwarz bounds near zero and there is something to skip. Three
      !! things, because a screen has three ways to be wrong:
      !!
      !!   * a threshold of zero must keep everything, to rounding;
      !!   * the default threshold must skip pairs and leave the matrix where
      !!     it was;
      !!   * a threshold loose enough to cut real contributions must visibly
      !!     change it, or the threshold was never read.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: Z2(6) = [8, 1, 1, 8, 1, 1]
      character(len=2), parameter :: SYM2(6) = ["O ", "H ", "H ", "O ", "H ", "H "]
      real(dp), parameter :: APART = 40.0_dp
      type(czt_molecule_t) :: mol
      type(cosx_grid_t) :: grid
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp) :: xyz(3, 6)
      real(dp), allocatable :: k_all(:, :), k_zero(:, :), k_default(:, :), k_loose(:, :)
      real(dp), allocatable :: ao(:, :), f(:, :), g(:, :)
      integer :: n_skipped, n_pairs

      xyz(:, 1:3) = WATER
      xyz(:, 4:6) = WATER
      xyz(1, 4:6) = xyz(1, 4:6) + APART
      call build_czt_molecule(Z2, SYM2, xyz, "6-31g*", mol, err)
      call check(error,.not. err%has_error(), "the molecule should build: "//err%get_full_trace())
      if (allocated(error)) return
      call run_czt_rhf(mol, 20, 60, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the SCF should converge")
      if (allocated(error)) return

      call build_cosx_grid(mol, grid, err, level=DEFAULT_COSX_GRID, role=COSX_ROLE_FINAL, &
                           screen=.true.)
      call check(error,.not. err%has_error(), "the grid should build: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, allocated(grid%bounds), "a grid built to screen should carry the bounds")
      if (allocated(error)) return

      call cosx_exchange(mol, grid%coords, grid%weights, scf%density, k_all, err)
      call cosx_exchange(mol, grid%coords, grid%weights, scf%density, k_zero, err, &
                         bounds=grid%bounds, screen_tol=0.0_dp)
      call cosx_exchange(mol, grid%coords, grid%weights, scf%density, k_default, err, &
                         bounds=grid%bounds)
      call cosx_exchange(mol, grid%coords, grid%weights, scf%density, k_loose, err, &
                         bounds=grid%bounds, screen_tol=1.0e-2_dp)
      call check(error,.not. err%has_error(), "COSX failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, maxval(abs(k_zero - k_all)) < 1.0e-13_dp, &
                 "a threshold of zero did not keep every shell pair")
      if (allocated(error)) return
      call check(error, maxval(abs(k_default - k_all)) < 1.0e-9_dp, &
                 "the default screening changed the exchange matrix")
      if (allocated(error)) return
      call check(error, maxval(abs(k_loose - k_all)) > 1.0e-6_dp, &
                 "a loose threshold changed nothing, so the threshold is not read")
      if (allocated(error)) return

      ! That the default skipped anything at all, on the first batch: without
      ! this the second check above passes with the screen doing nothing.
      call eval_ao_block(mol, grid%coords(:, 1:64), ao, err)
      call cosx_density_on_grid(ao, grid%weights(1:64), scf%density, f)
      call cosx_esp_contract(mol, grid%coords(:, 1:64), f, g, err, bounds=grid%bounds, &
                             n_skipped=n_skipped)
      call check(error,.not. err%has_error(), "the contraction failed: "//err%get_full_trace())
      if (allocated(error)) return
      n_pairs = mol%nbas*(mol%nbas + 1)/2
      call check(error, n_skipped > n_pairs/4, &
                 "the screening skipped too few pairs of two molecules far apart")
      if (allocated(error)) return

      call cosx_exchange(mol, grid%coords, grid%weights, scf%density, k_zero, err, &
                         bounds=grid%bounds(1:2, 1:2))
      call check(error, err%has_error(), "bounds of the wrong size were accepted")
      call mol%destroy()
   end subroutine test_screening

   subroutine test_batch_kernel(error)
      !! The exchange matrix through the batch potential integral is the library's, to rounding
      !!
      !! A grid carries the shells for the batch evaluator whenever the basis
      !! allows, and an SCF then uses them. Unscreened, the two routes are the
      !! same sum, so they agree to rounding; `test_mqc_czt_esp_batch` checks
      !! the integral itself shell form by shell form, and this checks that the
      !! exchange build reaches it and comes back with the same matrix.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(cosx_grid_t) :: iteration, final
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), k_library(:, :), k_batch(:, :)

      call water_density("cc-pvdz", mol, density, error)
      if (allocated(error)) return
      call build_cosx_grid(mol, iteration, err, level=2, role=COSX_ROLE_SCF, screen=.true.)
      call build_cosx_grid(mol, final, err, level=2, role=COSX_ROLE_FINAL, screen=.true.)
      call check(error,.not. err%has_error(), "the grids should build: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, final%shells%n_shell > 0, &
                 "a grid over a basis the batch evaluator takes should carry its shells")
      if (allocated(error)) return

      call cosx_exchange(mol, final%coords, final%weights, density, k_library, err)
      call cosx_exchange(mol, final%coords, final%weights, density, k_batch, err, &
                         shells=final%shells)
      call check(error,.not. err%has_error(), "COSX failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, maxval(abs(k_batch - k_library)) < 1.0e-11_dp, &
                 "the exchange matrix depends on which potential integral built it")
      if (allocated(error)) return

      ! The thresholds a grid asks for: looser where only the density has to
      ! come out right.
      call check(error, iteration%screen_tol == COSX_SCREEN_ITERATION .and. &
                 final%screen_tol == COSX_SCREEN_FINAL .and. &
                 COSX_SCREEN_ITERATION > COSX_SCREEN_FINAL, &
                 "the grids do not carry the screening thresholds of their roles")
      call mol%destroy()
   end subroutine test_batch_kernel

   subroutine test_trace(error)
      !! `cosx_exchange_trace` is `Tr[D K]` of the fitted, symmetrised matrix
      !!
      !! The final energy of an SCF comes from this and never forms `K`, so it
      !! has to be the same number the matrix would have given: the fit on the
      !! right side, the symmetrisation accounted for, the weights once.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(cosx_grid_t) :: grid
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), s(:, :), k(:, :)
      real(dp) :: from_matrix, plain, fitted, screened

      call water_density("cc-pvdz", mol, density, error)
      if (allocated(error)) return
      call mol%overlap(s)
      call build_cosx_grid(mol, grid, err, level=1, role=COSX_ROLE_SCF, overlap=s, screen=.true.)
      call check(error,.not. err%has_error(), "the grid should build: "//err%get_full_trace())
      if (allocated(error)) return

      ! Without the fit.
      call cosx_exchange(mol, grid%coords, grid%weights, density, k, err, shells=grid%shells)
      from_matrix = sum(density*k)
      call cosx_exchange_trace(mol, grid%coords, grid%weights, density, plain, err, grid%shells)
      call check(error,.not. err%has_error(), "COSX failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, abs(plain - from_matrix) < 1.0e-10_dp, &
                 "the unfitted trace is not the trace of the unfitted matrix")
      if (allocated(error)) return

      ! With it, which is what an SCF uses. This grid is coarse enough for the
      ! fit to move the trace by far more than the tolerance.
      call cosx_exchange(mol, grid%coords, grid%weights, density, k, err, fit=grid%fit, &
                         shells=grid%shells)
      call check(error, abs(sum(density*k) - from_matrix) > 1.0e-6_dp, &
                 "the fit did not change the trace, so this case cannot tell if it is applied")
      if (allocated(error)) return
      from_matrix = sum(density*k)
      call cosx_exchange_trace(mol, grid%coords, grid%weights, density, fitted, err, &
                               grid%shells, fit=grid%fit)
      call cosx_exchange_trace(mol, grid%coords, grid%weights, density, screened, err, &
                               grid%shells, fit=grid%fit, bounds=grid%bounds)
      call check(error,.not. err%has_error(), "COSX failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, abs(fitted - from_matrix) < 1.0e-10_dp, &
                 "the fitted trace is not the trace of the fitted matrix")
      if (allocated(error)) return
      call check(error, abs(screened - from_matrix) < 1.0e-8_dp, &
                 "screening moved the trace")
      if (allocated(error)) return
      ! Guards the guard: Tr[D K] of water is about four times its exchange
      ! energy of -9 hartree, in magnitude.
      call check(error, abs(from_matrix) > 20.0_dp, "the reference trace is implausible")
      call mol%destroy()
   end subroutine test_trace

   subroutine test_refusals(error)
      !! Shapes that cannot be right are refused, not read out of bounds
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: k(:, :), density(:, :)
      real(dp) :: points(3, 4), weights(4)

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return
      allocate (density(mol%nao, mol%nao))
      density = 0.0_dp
      points = 0.5_dp
      weights = 1.0_dp

      call cosx_exchange(mol, points, weights(1:3), density, k, err)
      call check(error, err%has_error(), "a weight count that differs from the points passed")
      if (allocated(error)) return
      call err%clear()

      call cosx_exchange(mol, points, weights, density(1:2, 1:2), k, err)
      call check(error, err%has_error(), "a density of the wrong size passed")
      if (allocated(error)) return
      call err%clear()

      call cosx_exchange(mol, points, weights, density, k, err, batch_size=0)
      call check(error, err%has_error(), "a batch size of zero passed")
      if (allocated(error)) return
      call err%clear()

      call cosx_exchange(mol, points(:, 1:0), weights(1:0), density, k, err)
      call check(error, err%has_error(), "an empty grid passed")
      call mol%destroy()
   end subroutine test_refusals

end module test_mqc_czt_cosx

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_cosx, only: collect_mqc_czt_cosx_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_cosx", collect_mqc_czt_cosx_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
