!! Direct two-electron builds, and what density symmetry each one may be given
module test_mqc_czt_direct
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
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, set_eri_path, &
                                eri_path_name, ROTAXIS_AVAILABLE, quartet_on_rotaxis, &
                                rotaxis_libfint_covers
   use mqc_czt_rhf, only: build_fock
   use mqc_czt_direct, only: build_fock_direct, build_fock_direct_many, &
                             build_fock_direct_nosym, schwarz_bounds, &
                             direct_stats_t
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_czt_direct_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

   !! Screening off, so a disagreement is about permutations and not thresholds.
   !!
   !! The whole point of these tests is the accumulation algebra. Leaving the
   !! Schwarz screening on would mix a second, unrelated approximation into every
   !! comparison and put a floor under the tolerances for no reason.
   real(dp), parameter :: NO_SCREENING = 0.0_dp

   !! Both limits of `erf(omega r)/r`, which is what an omega pass is.
   !!
   !! `omega` reaches libcint through a slot in `env`, not through a separate
   !! entry point, so handing the routine the shared environment instead of the
   !! local copy the omega was written into is completely silent: full-range
   !! integrals come back scaled by the long-range coefficient, nothing raises,
   !! and the Fock matrix is the right shape and the wrong operator. This build
   !! was one of the three places in the tree that did exactly that.
   !!
   !! Neither limit is sufficient alone. At `OMEGA_OFF` the kernel vanishes, so
   !! an ignored omega returns the full-range answer and fails; at `OMEGA_FULL`
   !! the kernel is `1/r`, which an ignored omega also satisfies -- but a
   !! short-range kernel, the sign convention inverted, fails it.
   real(dp), parameter :: OMEGA_OFF = 1.0e-6_dp
   real(dp), parameter :: OMEGA_FULL = 5.0e2_dp

contains

   subroutine collect_mqc_czt_direct_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("nosym_matches_the_fast_build_on_a_symmetric_density", &
                               test_nosym_symmetric), &
                  new_unittest("the_rotated_axis_path_builds_the_same_fock_matrix", &
                               test_rotaxis_fock), &
                  new_unittest("the_dispatch_agrees_with_libfint_on_what_the_path_covers", &
                               test_rotaxis_coverage), &
                  new_unittest("nosym_handles_an_antisymmetric_density", &
                               test_nosym_antisymmetric), &
                  new_unittest("the_fast_build_announced_antisymmetric_is_exact", &
                               test_fast_build_antisymmetric), &
                  new_unittest("a_wide_batch_contracts_through_the_blas_exactly", &
                               test_wide_batch), &
                  new_unittest("coulomb_vanishes_for_an_antisymmetric_density", &
                               test_coulomb_vanishes), &
                  new_unittest("the_fast_build_is_wrong_on_an_antisymmetric_density", &
                               test_fast_build_is_unsafe), &
                  new_unittest("an_omega_pass_is_actually_attenuated", &
                               test_attenuation_is_real) &
                  ]
   end subroutine collect_mqc_czt_direct_tests

   subroutine setup(mol, eri, bounds, zero_h, sym, anti, err)
      !! Water in 6-31G, its integrals, and one density of each symmetry
      type(czt_molecule_t), intent(out) :: mol
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
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
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

   subroutine test_rotaxis_coverage(error)
      !! `ROTAXIS_MAX_L` never sends libfint a quartet it does not cover
      !!
      !! Over every (i, j, i, j) quartet of water in cc-pVTZ, which has s, p,
      !! d and f shells. The dispatch may route *less* than the kernels cover
      !! -- it does, since the d classes measured slower than Rys -- but a
      !! quartet it routes that libfint refuses is an error stop in the middle
      !! of a Fock build, so that direction is held exactly.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp) :: c(3, 3)
      integer :: ish, jsh, n_on, n_off
      logical :: ours, theirs

      if (.not. ROTAXIS_AVAILABLE) then
         call check(error, .true.)
         return
      end if

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.4_dp, 1.1_dp, 0.0_dp, -1.4_dp, 1.1_dp], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "cc-pvtz", mol, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if

      n_on = 0
      n_off = 0
      do ish = 1, mol%nbas
         do jsh = 1, mol%nbas
            ours = quartet_on_rotaxis([ish - 1, jsh - 1, ish - 1, jsh - 1], mol%bas)
            theirs = rotaxis_libfint_covers([ish - 1, jsh - 1, ish - 1, jsh - 1], mol%bas, mol%nbas)
            if (ours .and. .not. theirs) then
               call check(error, .false., "the dispatch routes a quartet libfint does not cover")
               call mol%destroy()
               return
            end if
            if (ours) then
               n_on = n_on + 1
            else
               n_off = n_off + 1
            end if
         end do
      end do
      call check(error, n_on > 0 .and. n_off > 0, &
                 "cc-pVTZ water must have quartets on both sides of the limit")
      ! And the limit is what it says: an s/p pair is routed, a d one is not.
      if (.not. allocated(error)) then
         call check(error, quartet_on_rotaxis([0, 0, 0, 0], mol%bas) .and. &
                    .not. quartet_on_rotaxis([mol%nbas - 1, 0, mol%nbas - 1, 0], mol%bas), &
                    "the first shell of O must be routed and its last (f) must not")
      end if
      call mol%destroy()
   end subroutine test_rotaxis_coverage

   subroutine test_rotaxis_fock(error)
      !! The rotated-axis quartets agree with the Rys ones through a Fock build
      !!
      !! Water/6-31G is s, p and L shells, so every quartet takes the
      !! rotated-axis path when it is selected. A different algorithm, so not
      !! bit-identical: libfint holds the integrals to 1e-12 scaled, and the
      !! Fock matrix is a contraction over them. Skipped on a libcint build,
      !! which has the one path.
      !!
      !! `ROTAXIS_MAX_L` is checked against libfint's own answer for a d
      !! quartet and an f one, so the two ends of that agreement cannot drift
      !! apart silently: a wrong constant here would send a quartet the kernels
      !! lack into an error stop, or leave one they have on Rys.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :), rys(:, :), rotaxis(:, :)

      if (.not. ROTAXIS_AVAILABLE) then
         call set_eri_path("rotaxis", err)
         call check(error, err%has_error(), &
                    "a build without the path must refuse to be asked for it")
         return
      end if

      call setup(mol, eri, bounds, zero_h, sym, anti, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if

      allocate (rys(mol%nao, mol%nao), rotaxis(mol%nao, mol%nao))
      call set_eri_path("rys", err)
      call build_fock_direct(mol, zero_h, sym, bounds, rys, stats, err, &
                             screen_tol=NO_SCREENING)
      if (.not. err%has_error()) then
         call set_eri_path("rotaxis", err)
         call check(error, eri_path_name() == "rotaxis", "the path did not switch")
         if (.not. allocated(error)) then
            call build_fock_direct(mol, zero_h, sym, bounds, rotaxis, stats, err, &
                                   screen_tol=NO_SCREENING)
         end if
      end if
      ! Back to the default whatever happened, for the tests after this one.
      block
         type(error_t) :: reset
         call set_eri_path("rys", reset)
      end block
      call mol%destroy()
      if (allocated(error)) return
      if (err%has_error()) then
         call check(error, .false., "a build failed: "//err%get_message())
         return
      end if

      call check(error, maxval(abs(rotaxis - rys)) < 1.0e-10_dp, &
                 "the rotated-axis and Rys Fock matrices disagree")
   end subroutine test_rotaxis_fock

   subroutine test_nosym_symmetric(error)
      !! On a symmetric density the general build reduces to the fast one
      !!
      !! The compatibility statement. If this failed, the permutation enumeration
      !! would be wrong in a way that the antisymmetric tests might not catch,
      !! since they compare against a different reference.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
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

      type(czt_molecule_t) :: mol
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

   subroutine test_fast_build_antisymmetric(error)
      !! `build_fock_direct_many(antisymmetric=.true.)` matches the n^4 contraction
      !!
      !! The fast build's six updates fold the pair-swapped tuples into a final
      !! symmetrisation; announced antisymmetric it drops the Coulomb term and
      !! antisymmetrises instead. Held to the unpermuted reference and to
      !! `build_fock_direct_nosym`, on one density and on a batch that mixes
      !! magnitudes, since a batch is what the response solver hands it.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :), reference(:, :), fast(:, :, :)
      real(dp), allocatable :: general(:, :, :), batch(:, :, :)
      real(dp) :: worst

      call setup(mol, eri, bounds, zero_h, sym, anti, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if

      allocate (reference(mol%nao, mol%nao), batch(mol%nao, mol%nao, 3))
      batch(:, :, 1) = anti
      batch(:, :, 2) = 0.37_dp*anti
      batch(:, :, 3) = -2.5_dp*anti
      call build_fock(zero_h, eri, anti, reference)
      call build_fock_direct_many(mol, zero_h, batch, bounds, fast, stats, err, &
                                  screen_tol=NO_SCREENING, antisymmetric=.true.)
      if (err%has_error()) then
         call check(error, .false., "the fast build failed: "//err%get_message())
         call mol%destroy()
         return
      end if
      call build_fock_direct_nosym(mol, zero_h, batch, bounds, general, stats, err, &
                                   screen_tol=NO_SCREENING)
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "the general build failed: "//err%get_message())
         return
      end if

      worst = maxval(abs(fast(:, :, 1) - reference))
      call check(error, worst < 1.0e-11_dp, &
                 "the fast build announced antisymmetric does not reproduce the "// &
                 "unpermuted contraction")
      if (allocated(error)) return
      worst = max(maxval(abs(fast(:, :, 2) - 0.37_dp*reference)), &
                  maxval(abs(fast(:, :, 3) + 2.5_dp*reference)))
      call check(error, worst < 1.0e-10_dp, &
                 "the fast build announced antisymmetric is wrong on a batch")
      if (allocated(error)) return
      worst = maxval(abs(fast - general))
      call check(error, worst < 1.0e-10_dp, &
                 "the fast and the general build disagree on an antisymmetric batch")
      if (allocated(error)) return
      call check(error, maxval(abs(fast(:, :, 1) + transpose(fast(:, :, 1)))) < 1.0e-11_dp, &
                 "G of an antisymmetric density came back with a symmetric part")
   end subroutine test_fast_build_antisymmetric

   subroutine test_wide_batch(error)
      !! A batch wide enough for the matrix-product contraction, both symmetries
      !!
      !! Below `GEMM_SETS` densities the fast build updates element by element;
      !! at or above it each quartet goes through six matrix products. The two
      !! must agree to rounding with each other and with the unpermuted
      !! reference, on a batch of scaled copies so that every set is distinct.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: WIDE = 24
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :), ref_sym(:, :), ref_anti(:, :)
      real(dp), allocatable :: batch(:, :, :), fast(:, :, :)
      real(dp) :: worst, scale
      integer :: m

      call setup(mol, eri, bounds, zero_h, sym, anti, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if
      allocate (ref_sym(mol%nao, mol%nao), ref_anti(mol%nao, mol%nao), batch(mol%nao, mol%nao, WIDE))
      call build_fock(zero_h, eri, sym, ref_sym)
      call build_fock(zero_h, eri, anti, ref_anti)

      do m = 1, WIDE
         batch(:, :, m) = (0.1_dp*m - 1.3_dp)*sym
      end do
      call build_fock_direct_many(mol, zero_h, batch, bounds, fast, stats, err, &
                                  screen_tol=NO_SCREENING)
      if (err%has_error()) then
         call check(error, .false., "the wide symmetric build failed: "//err%get_message())
         call mol%destroy()
         return
      end if
      worst = 0.0_dp
      do m = 1, WIDE
         scale = 0.1_dp*m - 1.3_dp
         worst = max(worst, maxval(abs(fast(:, :, m) - scale*ref_sym)))
      end do
      call check(error, worst < 1.0e-10_dp, "a wide symmetric batch through the BLAS is wrong")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      do m = 1, WIDE
         batch(:, :, m) = (0.1_dp*m - 1.3_dp)*anti
      end do
      call build_fock_direct_many(mol, zero_h, batch, bounds, fast, stats, err, &
                                  screen_tol=NO_SCREENING, antisymmetric=.true.)
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "the wide antisymmetric build failed: "//err%get_message())
         return
      end if
      worst = 0.0_dp
      do m = 1, WIDE
         scale = 0.1_dp*m - 1.3_dp
         worst = max(worst, maxval(abs(fast(:, :, m) - scale*ref_anti)))
      end do
      call check(error, worst < 1.0e-10_dp, "a wide antisymmetric batch through the BLAS is wrong")
   end subroutine test_wide_batch

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

      type(czt_molecule_t) :: mol
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

      type(czt_molecule_t) :: mol
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

   subroutine test_attenuation_is_real(error)
      !! A long-range exchange pass must be long-range
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :)
      real(dp), allocatable :: zero_h(:, :), sym(:, :), anti(:, :)
      real(dp), allocatable :: dens(:, :, :)
      real(dp), allocatable :: full(:, :, :), faint(:, :, :), most(:, :, :)
      real(dp) :: scale

      call setup(mol, eri, bounds, zero_h, sym, anti, err)
      call check(error,.not. err%has_error(), "the molecule failed to build")
      if (allocated(error)) return

      ! Exchange only, which is the shape a long-range pass is always called
      ! in -- `j_scale = 0` because the full-range pass has already supplied
      ! the Coulomb term. It also keeps the norms below measuring the thing
      ! under test rather than a Coulomb matrix that dwarfs it.
      allocate (dens(size(sym, 1), size(sym, 2), 1))
      dens(:, :, 1) = sym

      call build_fock_direct_many(mol, zero_h, dens, bounds, full, stats, err, &
                                  screen_tol=NO_SCREENING, j_scale=0.0_dp)
      call build_fock_direct_many(mol, zero_h, dens, bounds, faint, stats, err, &
                                  screen_tol=NO_SCREENING, j_scale=0.0_dp, omega=OMEGA_OFF)
      call build_fock_direct_many(mol, zero_h, dens, bounds, most, stats, err, &
                                  screen_tol=NO_SCREENING, j_scale=0.0_dp, omega=OMEGA_FULL)
      call check(error,.not. err%has_error(), "the Fock build failed")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      scale = maxval(abs(full))
      call check(error, scale > 1.0e-2_dp, "the full-range exchange build is empty")
      if (.not. allocated(error)) &
         call check(error, maxval(abs(faint)) < 1.0e-4_dp*scale, &
                    "an omega pass through build_fock_direct_many is not attenuated")
      if (.not. allocated(error)) &
         call check(error, maxval(abs(most - full)) < 1.0e-2_dp*scale, &
                    "build_fock_direct_many at large omega does not recover full-range exchange")
      call mol%destroy()
   end subroutine test_attenuation_is_real

end module test_mqc_czt_direct

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_direct, only: collect_mqc_czt_direct_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_direct", collect_mqc_czt_direct_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
