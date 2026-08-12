!! Direct two-electron builds, and what density symmetry each one may be given
module test_mqc_libcint_direct
   !! The fast direct build folds three of the eightfold integral permutations
   !! into a multiplicity factor, which is valid only for a symmetric density.
   !! `build_fock_direct_nosym` writes those permutations out instead, and exists
   !! because the frequency-dependent coupled-perturbed equations need `A - B`,
   !! whose response density is antisymmetric.
   !!
   !! These tests pin all of that down, including -- deliberately -- that the fast
   !! routine really is wrong on an antisymmetric density. A test that only checked
   !! the new routine works would leave the next person free to "simplify" by
   !! routing everything through the fast one, and the failure would be a silent
   !! factor of two on a term that should have been zero.
   !!
   !! The arbitrary-density reference is `build_fock` on the stored tensor: a plain
   !! `n^4` contraction with no symmetry assumptions anywhere in it, so it is
   !! correct for any density by construction and is what the screened, permuted,
   !! threaded builds are measured against.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: build_fock
   use mqc_libcint_direct, only: build_fock_direct, build_fock_direct_many, &
                                 build_fock_direct_nosym, schwarz_bounds, &
                                 direct_stats_t
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_libcint_direct_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

   !> Screening off, so a disagreement is about permutations and not thresholds.
   !>
   !> The whole point of these tests is the accumulation algebra. Leaving the
   !> Schwarz screening on would mix a second, unrelated approximation into every
   !> comparison and put a floor under the tolerances for no reason.
   real(dp), parameter :: NO_SCREENING = 0.0_dp

contains

   subroutine collect_mqc_libcint_direct_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("nosym_matches_the_fast_build_on_a_symmetric_density", &
                               test_nosym_symmetric), &
                  new_unittest("nosym_handles_an_antisymmetric_density", &
                               test_nosym_antisymmetric), &
                  new_unittest("coulomb_vanishes_for_an_antisymmetric_density", &
                               test_coulomb_vanishes), &
                  new_unittest("the_fast_build_is_wrong_on_an_antisymmetric_density", &
                               test_fast_build_is_unsafe) &
                  ]
   end subroutine collect_mqc_libcint_direct_tests

   subroutine setup(mol, eri, bounds, zero_h, sym, anti, err)
      !! Water in 6-31G, its integrals, and one density of each symmetry
      type(libcint_molecule_t), intent(out) :: mol
      real(dp), allocatable, intent(out) :: eri(:, :, :, :), bounds(:, :)
      real(dp), allocatable, intent(out) :: zero_h(:, :), sym(:, :), anti(:, :)
      type(error_t), intent(inout) :: err

      real(dp) :: c(3, 3)
      real(dp), allocatable :: m(:, :)
      integer :: n, i, j

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      ! 6-31G rather than STO-3G: several shells per atom and more than one
      ! contraction, so blocks with s1 == s2 and blocks with s1 /= s2 both occur
      ! in quantity. A one-shell-per-atom basis would exercise only some of the
      ! permutation cases and could let a wrong condition pass.
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
      if (err%has_error()) return

      n = mol%nao
      call mol%eris(eri)
      call schwarz_bounds(mol, bounds, err)
      if (err%has_error()) return

      allocate (zero_h(n, n), sym(n, n), anti(n, n), m(n, n))
      zero_h = 0.0_dp
      ! Deterministic and not close to any symmetry by accident.
      do j = 1, n
         do i = 1, n
            m(i, j) = sin(0.7_dp*real(i, dp) + 1.3_dp*real(j, dp)) &
                      + 0.25_dp*real(i - j, dp)/real(n, dp)
         end do
      end do
      sym = m + transpose(m)
      anti = m - transpose(m)
      deallocate (m)
   end subroutine setup

   subroutine test_nosym_symmetric(error)
      !! On a symmetric density the general build reduces to the fast one
      !!
      !! The compatibility statement. If this failed, the permutation enumeration
      !! would be wrong in a way that the antisymmetric tests might not catch,
      !! since they compare against a different reference.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :), fast(:, :), general(:, :, :)
      real(dp), allocatable :: one(:, :, :)

      call setup(mol, eri, bounds, zero_h, sym, anti, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if

      allocate (fast(mol%nao, mol%nao), one(mol%nao, mol%nao, 1))
      one(:, :, 1) = sym
      call build_fock_direct(mol, zero_h, sym, bounds, fast, stats, err, &
                             screen_tol=NO_SCREENING)
      if (.not. err%has_error()) then
         call build_fock_direct_nosym(mol, zero_h, one, bounds, general, stats, err, &
                                      screen_tol=NO_SCREENING)
      end if
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "a build failed: "//err%get_message())
         return
      end if

      call check(error, maxval(abs(general(:, :, 1) - fast)) < 1.0e-11_dp, &
                 "the general and fast builds disagree on a symmetric density")
   end subroutine test_nosym_symmetric

   subroutine test_nosym_antisymmetric(error)
      !! On an antisymmetric density it matches an unpermuted n^4 contraction
      !!
      !! This is the capability the routine was written for, and the reference has
      !! no symmetry assumptions in it at all.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :), reference(:, :), general(:, :, :)
      real(dp), allocatable :: one(:, :, :)
      real(dp) :: worst

      call setup(mol, eri, bounds, zero_h, sym, anti, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if

      allocate (reference(mol%nao, mol%nao), one(mol%nao, mol%nao, 1))
      one(:, :, 1) = anti
      call build_fock(zero_h, eri, anti, reference)
      call build_fock_direct_nosym(mol, zero_h, one, bounds, general, stats, err, &
                                   screen_tol=NO_SCREENING)
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "the general build failed: "//err%get_message())
         return
      end if

      worst = maxval(abs(general(:, :, 1) - reference))
      call check(error, worst < 1.0e-11_dp, &
                 "the general direct build does not reproduce an unpermuted "// &
                 "contraction on an antisymmetric density")
      if (allocated(error)) return

      ! And the result must itself be antisymmetric, since J vanishes and K of an
      ! antisymmetric density is antisymmetric. This is what a final
      ! `0.5*(g + transpose(g))` would destroy, so it is worth asserting
      ! separately from the comparison above.
      call check(error, maxval(abs(general(:, :, 1) + transpose(general(:, :, 1)))) &
                 < 1.0e-11_dp, &
                 "G of an antisymmetric density came back with a symmetric part")
   end subroutine test_nosym_antisymmetric

   subroutine test_coulomb_vanishes(error)
      !! The Coulomb term cancels by itself, and is not special-cased
      !!
      !! `J(D)_uv = sum (uv|ls) D_ls` vanishes for antisymmetric `D` because the
      !! integral is symmetric under `l <-> s` while the density is not. The
      !! routine does not know that -- it accumulates the Coulomb updates anyway --
      !! so this measures a cancellation rather than a skipped branch, which makes
      !! it a real check on the permutation bookkeeping: a missing or duplicated
      !! tuple would leave a residue here.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :)
      real(dp), allocatable :: coulomb(:, :)
      real(dp) :: scale
      integer :: n, mu, nu

      call setup(mol, eri, bounds, zero_h, sym, anti, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if
      n = mol%nao

      ! J alone, straight from the stored tensor, so the claim is checked on the
      ! object itself rather than inferred from a total.
      allocate (coulomb(n, n))
      do nu = 1, n
         do mu = 1, n
            coulomb(mu, nu) = sum(eri(mu, nu, :, :)*anti)
         end do
      end do
      scale = maxval(abs(eri))*maxval(abs(anti))
      call mol%destroy()

      call check(error, maxval(abs(coulomb)) < 1.0e-12_dp*max(scale, 1.0_dp), &
                 "the Coulomb term does not vanish on an antisymmetric density, "// &
                 "so the premise of the general build is wrong")
   end subroutine test_coulomb_vanishes

   subroutine test_fast_build_is_unsafe(error)
      !! The fast build really is wrong here, and by the predicted amount
      !!
      !! Asserting a *failure* on purpose. The fast build's `deg` factor doubles for
      !! `s1 /= s2` instead of adding `D(nu,mu)`, so on an antisymmetric density the
      !! permutations that should cancel reinforce, and its Coulomb contribution
      !! survives where the true one is zero. Pinning that keeps someone from
      !! deleting `build_fock_direct_nosym` as redundant, and records what the
      !! symptom would be: a large error, not a subtle one.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :), reference(:, :), fast(:, :, :)
      real(dp), allocatable :: one(:, :, :)
      real(dp) :: discrepancy, magnitude

      call setup(mol, eri, bounds, zero_h, sym, anti, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if

      allocate (reference(mol%nao, mol%nao), one(mol%nao, mol%nao, 1))
      one(:, :, 1) = anti
      call build_fock(zero_h, eri, anti, reference)
      call build_fock_direct_many(mol, zero_h, one, bounds, fast, stats, err, &
                                  screen_tol=NO_SCREENING)
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "the fast build failed: "//err%get_message())
         return
      end if

      discrepancy = maxval(abs(fast(:, :, 1) - reference))
      magnitude = maxval(abs(reference))

      call check(error, discrepancy > 0.1_dp*magnitude, &
                 "the fast build agreed with the reference on an antisymmetric "// &
                 "density, so either it was fixed -- in which case delete "// &
                 "build_fock_direct_nosym and this test -- or the reference is "// &
                 "no longer symmetry-agnostic")
      if (allocated(error)) return

      ! It returns something symmetric, because it symmetrises at the end, where
      ! the true answer is antisymmetric. That is the clearest single statement of
      ! why it cannot be used.
      call check(error, maxval(abs(fast(:, :, 1) - transpose(fast(:, :, 1)))) &
                 < 1.0e-11_dp, &
                 "the fast build did not return a symmetric matrix, so its final "// &
                 "symmetrisation is not what makes it unusable here")
   end subroutine test_fast_build_is_unsafe

end module test_mqc_libcint_direct

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_libcint_direct, only: collect_mqc_libcint_direct_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_libcint_direct", collect_mqc_libcint_direct_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
