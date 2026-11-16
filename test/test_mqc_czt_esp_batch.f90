!! The batch potential integral against the integral library's
module test_mqc_czt_esp_batch
   !! `esp_batch_contract` computes what `cosx_esp_contract` gets from libfint's
   !! `int1e_grids`, by a different algorithm written here. libfint is checked
   !! against libcint value for value, so agreeing with it to rounding on every
   !! angular momentum pair is the whole test: there is no tolerance for an
   !! approximation to hide in, because there is no approximation.
   !!
   !! The comparison is on `G = A F` with `F` arbitrary and different at every
   !! point, which reaches every element of `A`; a second pass with `F` a unit
   !! vector isolates single columns, so a compensating pair of errors cannot
   !! pass the first by cancelling.
   !!
   !! One basis per thing that can go wrong separately: s and p only, a
   !! Cartesian d, a spherical d, f, and g, the last being the highest the
   !! evaluator takes.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_cosx, only: cosx_esp_contract
   use mqc_czt_esp_batch, only: esp_shells_t, build_esp_shells, esp_batch_contract, &
                                esp_batch_energy, boys_values, ESP_BATCH_LMAX
   implicit none
   private

   public :: collect_mqc_czt_esp_batch_tests

   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [3, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   integer, parameter :: N_PTS = 23

contains

   subroutine collect_mqc_czt_esp_batch_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("boys_function_matches_its_series", test_boys), &
                  new_unittest("esp_batch_s_and_p", test_sp), &
                  new_unittest("esp_batch_cartesian_d", test_cart_d), &
                  new_unittest("esp_batch_spherical_d", test_sph_d), &
                  new_unittest("esp_batch_f", test_f), &
                  new_unittest("esp_batch_g", test_g), &
                  new_unittest("esp_batch_far_and_near_batches", test_regimes), &
                  new_unittest("esp_batch_screens_as_the_library_path_does", test_screen), &
                  new_unittest("esp_batch_energy_is_the_contraction_traced", test_energy), &
                  new_unittest("esp_batch_refuses_what_it_cannot_take", test_refusals) &
                  ]
   end subroutine collect_mqc_czt_esp_batch_tests

   subroutine probe_points(points)
      !! Points near, between and far from the atoms, none on a symmetry element
      real(dp), intent(out) :: points(3, N_PTS)

      integer :: ig
      real(dp) :: s

      do ig = 1, N_PTS
         s = real(ig, dp)
         points(1, ig) = 0.31_dp*sin(1.3_dp*s) + 0.07_dp*s - 0.6_dp
         points(2, ig) = 1.9_dp*cos(0.7_dp*s + 0.2_dp)
         points(3, ig) = 0.9_dp*sin(2.1_dp*s + 1.0_dp) + 0.4_dp
      end do
      ! Far enough that the Boys argument leaves the tabulated range for every
      ! primitive pair, and close enough to a nucleus that it is near zero.
      points(:, N_PTS) = [14.0_dp, -9.0_dp, 11.0_dp]
      points(:, N_PTS - 1) = [1.0e-3_dp, -2.0e-3_dp, 1.5e-3_dp]
   end subroutine probe_points

   subroutine compare(basis, error, expect_l, cluster)
      !! Both evaluators on one basis, with an arbitrary F and then unit vectors
      character(len=*), intent(in) :: basis
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: expect_l   !! Highest angular momentum the basis must reach
      real(dp), intent(in), optional :: cluster(4)
         !! Centre and radius, Bohr: the points are drawn into this sphere
         !! instead of being spread from the nuclei outwards

      type(czt_molecule_t) :: mol
      type(esp_shells_t) :: shells
      type(error_t) :: err
      real(dp) :: points(3, N_PTS)
      real(dp), allocatable :: f(:, :), g_ref(:, :), g_new(:, :)
      real(dp) :: worst, scale
      integer :: ig, mu, column

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, basis, mol, err)
      call check(error,.not. err%has_error(), "the molecule should build: "//err%get_full_trace())
      if (allocated(error)) return
      call build_esp_shells(mol, shells, err)
      call check(error,.not. err%has_error(), "the shells should build: "//err%get_full_trace())
      if (allocated(error)) return
      ! Guards the guard: a basis that did not reach the angular momentum the
      ! case is named for would pass while testing something lower.
      call check(error, shells%l_max == expect_l, &
                 "the basis does not reach the angular momentum this case is for")
      if (allocated(error)) return

      call probe_points(points)
      if (present(cluster)) then
         ! The spread points are within about 25 Bohr of the origin; scaled by
         ! that they sit inside the unit sphere.
         do ig = 1, N_PTS
            points(:, ig) = cluster(1:3) + cluster(4)*points(:, ig)/25.0_dp
         end do
      end if
      allocate (f(N_PTS, mol%nao))
      do mu = 1, mol%nao
         do ig = 1, N_PTS
            f(ig, mu) = sin(1.7_dp*real(ig, dp) + 0.37_dp*real(mu, dp))
         end do
      end do

      call cosx_esp_contract(mol, points, f, g_ref, err)
      call esp_batch_contract(shells, points, f, g_new, err)
      call check(error,.not. err%has_error(), "a contraction failed: "//err%get_full_trace())
      if (allocated(error)) return
      worst = maxval(abs(g_new - g_ref))
      scale = maxval(abs(g_ref))
      call check(error, worst < 1.0e-11_dp*max(1.0_dp, scale), &
                 "the batch integral differs from the library's on "//basis)
      if (allocated(error)) return
      call check(error, scale > 1.0e-3_dp, "the contracted values are too small to test")
      if (allocated(error)) return

      ! Single columns of the integral, first, middle and last function.
      do column = 1, mol%nao, max(1, (mol%nao - 1)/2)
         f = 0.0_dp
         f(:, column) = 1.0_dp
         call cosx_esp_contract(mol, points, f, g_ref, err)
         call esp_batch_contract(shells, points, f, g_new, err)
         call check(error,.not. err%has_error(), "a contraction failed: "//err%get_full_trace())
         if (allocated(error)) return
         call check(error, maxval(abs(g_new - g_ref)) < 1.0e-11_dp*max(1.0_dp, maxval(abs(g_ref))), &
                    "a column of the batch integral differs from the library's on "//basis)
         if (allocated(error)) return
      end do
      call mol%destroy()
   end subroutine compare

   subroutine test_sp(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("6-31g", error, 1)
   end subroutine test_sp

   subroutine test_cart_d(error)
      !! 6-31G* is Cartesian here, so the block is used as it stands
      type(error_type), allocatable, intent(out) :: error
      call compare("6-31g*", error, 2)
   end subroutine test_cart_d

   subroutine test_sph_d(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("cc-pvdz", error, 2)
   end subroutine test_sph_d

   subroutine test_f(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("def2-tzvp", error, 3)
   end subroutine test_f

   subroutine test_g(error)
      type(error_type), allocatable, intent(out) :: error
      call compare("def2-qzvp", error, ESP_BATCH_LMAX)
   end subroutine test_g

   subroutine test_regimes(error)
      !! Batches that lie wholly beyond the Boys table, and wholly inside it
      !!
      !! The evaluator picks one way of seeding the Coulomb integrals for a
      !! whole batch when the batch's bounding sphere allows it. The spread
      !! points of the other cases straddle the table for most primitive pairs,
      !! so they exercise the point-by-point choice; these two exercise the
      !! shortcuts. def2-TZVP has both the tight and the diffuse exponents.
      type(error_type), allocatable, intent(out) :: error

      call compare("def2-tzvp", error, 3, cluster=[14.0_dp, -9.0_dp, 11.0_dp, 0.5_dp])
      if (allocated(error)) return
      call compare("def2-tzvp", error, 3, cluster=[0.05_dp, -0.1_dp, 0.2_dp, 0.2_dp])
   end subroutine test_regimes

   subroutine test_boys(error)
      !! `F_n(x)` against the series that defines it, across both branches
      !!
      !! `F_n(x) = exp(-x) sum_k (2x)^k / ((2n+1)(2n+3)...(2n+2k+1))`, summed
      !! here to convergence. The arguments straddle the mesh nodes, sit on
      !! them, and cross into the asymptotic range.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N_X = 12, N_MAX = 2*ESP_BATCH_LMAX
      real(dp), parameter :: X(N_X) = [0.0_dp, 1.0e-9_dp, 0.049_dp, 0.05_dp, 0.3_dp, &
                                       2.7183_dp, 9.95_dp, 17.0_dp, 29.999_dp, &
                                       35.97_dp, 36.5_dp, 250.0_dp]
      type(czt_molecule_t) :: mol
      type(esp_shells_t) :: shells
      type(error_t) :: err
      real(dp) :: values(N_X, 0:N_MAX), reference, term, total, worst
      integer :: ix, n, k

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call build_esp_shells(mol, shells, err)
      call check(error,.not. err%has_error(), "the shells should build: "//err%get_full_trace())
      if (allocated(error)) return
      call boys_values(shells, N_MAX, X, values)

      worst = 0.0_dp
      do n = 0, N_MAX
         do ix = 1, N_X
            term = 1.0_dp/real(2*n + 1, dp)
            total = term
            do k = 1, 2000
               term = term*2.0_dp*X(ix)/real(2*n + 2*k + 1, dp)
               total = total + term
               if (term < 1.0e-18_dp*total) exit
            end do
            reference = exp(-X(ix))*total
            ! Relative where the value is of order one, absolute where it has
            ! decayed: at x = 250 and n = 8 it is 1e-21.
            worst = max(worst, abs(values(ix, n) - reference)/max(reference, 1.0e-3_dp))
         end do
      end do
      call check(error, worst < 1.0e-12_dp, "the Boys function is off its defining series")
      call mol%destroy()
   end subroutine test_boys

   subroutine test_screen(error)
      !! With bounds, the same pairs are skipped and the same matrix comes back
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: Z2(6) = [8, 1, 1, 8, 1, 1]
      character(len=2), parameter :: SYM2(6) = ["O ", "H ", "H ", "O ", "H ", "H "]
      type(czt_molecule_t) :: mol
      type(esp_shells_t) :: shells
      type(error_t) :: err
      real(dp) :: xyz(3, 6), points(3, N_PTS)
      real(dp), allocatable :: f(:, :), g_ref(:, :), g_new(:, :), bounds(:, :)
      integer :: ig, mu, skipped_ref, skipped_new

      xyz(:, 1:3) = WATER
      xyz(:, 4:6) = WATER
      xyz(1, 4:6) = xyz(1, 4:6) + 40.0_dp
      ! A segmented basis, so that a shell is one contraction and the two
      ! evaluators count pairs the same way. On a generally contracted one the
      ! batch evaluator counts pairs of contractions.
      call build_czt_molecule(Z2, SYM2, xyz, "def2-svp", mol, err)
      call build_esp_shells(mol, shells, err)
      call schwarz_bounds(mol, bounds, err)
      call check(error,.not. err%has_error(), "setup failed: "//err%get_full_trace())
      if (allocated(error)) return

      call probe_points(points)
      allocate (f(N_PTS, mol%nao))
      ! F that lives on the first molecule only, as a density near these
      ! points would.
      f = 0.0_dp
      do mu = 1, mol%nao/2
         do ig = 1, N_PTS
            f(ig, mu) = sin(1.7_dp*real(ig, dp) + 0.37_dp*real(mu, dp))
         end do
      end do

      call cosx_esp_contract(mol, points, f, g_ref, err, bounds=bounds, &
                             screen_tol=1.0e-10_dp, n_skipped=skipped_ref)
      call esp_batch_contract(shells, points, f, g_new, err, bounds=bounds, &
                              screen_tol=1.0e-10_dp, n_skipped=skipped_new)
      call check(error,.not. err%has_error(), "a contraction failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, skipped_new == skipped_ref, &
                 "the two evaluators skipped different numbers of shell pairs")
      if (allocated(error)) return
      call check(error, skipped_new > 0, "nothing was screened, so nothing was tested")
      if (allocated(error)) return
      call check(error, maxval(abs(g_new - g_ref)) < 1.0e-11_dp, &
                 "the screened batch integral differs from the library's")
      call mol%destroy()
   end subroutine test_screen

   subroutine test_energy(error)
      !! `esp_batch_energy` is `sum Y G` over what `esp_batch_contract` returns
      !!
      !! Unscreened the two are the same sum in a different order, so they
      !! agree to rounding. Screened, the energy form may skip more -- it
      !! knows where `Y` is small as well as where `F` is -- and must still
      !! land on the same number.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(esp_shells_t) :: shells
      type(error_t) :: err
      real(dp) :: points(3, N_PTS)
      real(dp), allocatable :: f(:, :), y(:, :), g(:, :), bounds(:, :)
      real(dp) :: direct, traced, screened
      integer :: ig, mu, skipped_matrix, skipped_energy

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "def2-svp", mol, err)
      call build_esp_shells(mol, shells, err)
      call schwarz_bounds(mol, bounds, err)
      call check(error,.not. err%has_error(), "setup failed: "//err%get_full_trace())
      if (allocated(error)) return

      call probe_points(points)
      allocate (f(N_PTS, mol%nao), y(N_PTS, mol%nao))
      ! F on the first half of the functions and Y negligible on the second:
      ! a pair from the second half to the first carries F and so the matrix
      ! form evaluates it, while the energy form sees that nothing reads the
      ! result.
      f = 0.0_dp
      do mu = 1, mol%nao
         do ig = 1, N_PTS
            if (mu <= mol%nao/2) f(ig, mu) = sin(1.7_dp*real(ig, dp) + 0.37_dp*real(mu, dp))
            y(ig, mu) = cos(0.9_dp*real(ig, dp) - 0.53_dp*real(mu, dp))
         end do
      end do
      y(:, mol%nao/2 + 1:) = 1.0e-14_dp*y(:, mol%nao/2 + 1:)

      call esp_batch_contract(shells, points, f, g, err)
      direct = sum(y*g)
      call esp_batch_energy(shells, points, f, y, traced, err)
      call check(error,.not. err%has_error(), "a contraction failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, abs(traced - direct) < 1.0e-11_dp*max(1.0_dp, abs(direct)), &
                 "the traced contraction is not the trace of the contraction")
      if (allocated(error)) return
      call check(error, abs(direct) > 1.0e-3_dp, "the trace is too small to test")
      if (allocated(error)) return

      call esp_batch_contract(shells, points, f, g, err, bounds=bounds, &
                              screen_tol=1.0e-10_dp, n_skipped=skipped_matrix)
      call esp_batch_energy(shells, points, f, y, screened, err, bounds=bounds, &
                            screen_tol=1.0e-10_dp, n_skipped=skipped_energy)
      call check(error,.not. err%has_error(), "a contraction failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, abs(screened - direct) < 1.0e-9_dp, &
                 "screening moved the traced contraction")
      if (allocated(error)) return
      call check(error, skipped_energy > skipped_matrix, &
                 "the energy form skipped no more than the matrix form, so the "// &
                 "second factor is not screening anything")
      call mol%destroy()
   end subroutine test_energy

   subroutine test_refusals(error)
      !! Shapes that cannot be right are refused
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(esp_shells_t) :: shells
      type(error_t) :: err
      real(dp), allocatable :: f(:, :), g(:, :)
      real(dp) :: points(3, 4)

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call build_esp_shells(mol, shells, err)
      call check(error,.not. err%has_error(), "the shells should build")
      if (allocated(error)) return
      points = 0.5_dp
      allocate (f(3, mol%nao))
      f = 0.0_dp
      call esp_batch_contract(shells, points, f, g, err)
      call check(error, err%has_error(), "an F with the wrong number of points was accepted")
      if (allocated(error)) return
      call err%clear()
      deallocate (f)
      allocate (f(4, mol%nao + 1))
      f = 0.0_dp
      call esp_batch_contract(shells, points, f, g, err)
      call check(error, err%has_error(), "an F with the wrong number of functions was accepted")
      call mol%destroy()
   end subroutine test_refusals

end module test_mqc_czt_esp_batch

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_esp_batch, only: collect_mqc_czt_esp_batch_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_esp_batch", collect_mqc_czt_esp_batch_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
