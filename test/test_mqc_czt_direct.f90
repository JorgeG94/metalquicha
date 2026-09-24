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
                                rotaxis_libfint_covers, HGP_AVAILABLE, quartet_on_hgp, &
                                hgp_libfint_covers, eri_grad_dispatch_t, &
                                build_eri_grad_dispatch, rotaxis_grad_cached, &
                                hgp_grad_cached, rotaxis_grad_libfint, hgp_grad_libfint
   use mqc_czt_rhf, only: build_fock
   use mqc_czt_direct, only: build_fock_direct, build_fock_direct_many, &
                             build_fock_direct_nosym, build_fock_direct_uhf, &
                             build_fock_direct_uhf_many, schwarz_bounds, &
                             direct_stats_t, coulomb_from_blocks
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

   !! CAM-B3LYP's own numbers, because they are the ones a response operator
   !! will actually pass: `alpha + beta = 0.19` on the full-range exchange and
   !! `-beta = 0.46` on the attenuated pass at `omega = 0.33`. Nothing here
   !! depends on the functional; they are a triple of coefficients no two of
   !! which are equal, so a routine that confused `k_scale` with `j_scale` or
   !! dropped one of them cannot pass by coincidence.
   real(dp), parameter :: CAM_K_FULL = 0.19_dp
   real(dp), parameter :: CAM_K_LR = 0.46_dp
   real(dp), parameter :: CAM_OMEGA = 0.33_dp

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
                  new_unittest("the_head_gordon_pople_path_builds_the_same_fock_matrix", &
                               test_hgp_fock), &
                  new_unittest("the_dispatch_agrees_with_libfint_on_what_hgp_covers", &
                               test_hgp_coverage), &
                  new_unittest("the_gradient_dispatch_agrees_with_libfint_quartet_by_quartet", &
                               test_grad_dispatch), &
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
                               test_attenuation_is_real), &
                  new_unittest("nosym_takes_the_same_scales_as_the_fast_build", &
                               test_nosym_scales), &
                  new_unittest("the_stored_build_takes_the_same_scales_as_the_direct_one", &
                               test_stored_scales), &
                  new_unittest("nosym_scales_exchange_on_an_antisymmetric_density", &
                               test_nosym_antisymmetric_scaled), &
                  new_unittest("the_uhf_batch_is_the_closed_shell_build_on_equal_spins", &
                               test_uhf_many_closed_shell), &
                  new_unittest("the_uhf_batch_matches_an_explicit_contraction", &
                               test_uhf_many_explicit), &
                  new_unittest("the_uhf_batch_handles_antisymmetric_pairs", &
                               test_uhf_many_antisymmetric), &
                  new_unittest("the_block_coulomb_build_is_the_bra_block_of_the_full_one", &
                               test_block_coulomb) &
                  ]
   end subroutine collect_mqc_czt_direct_tests

   subroutine test_block_coulomb(error)
      !! `coulomb_from_blocks` is the bra block of the full Coulomb build
      !!
      !! The layout an FMO near field hands it: glycine tripeptide and water
      !! cut into four pieces across the two C-alpha bonds, a group first and
      !! three neighbours after it, each with the detached atom its SCF held as
      !! a ghost. A ghost sits on the same centre as the real atom another
      !! block owns, so the supersystem carries duplicated centres, exactly as
      !! the FMO one does. Each neighbour gets a density of its own; the group
      !! gets none.
      !!
      !! The reference is `build_fock_direct` over the whole supersystem with
      !! exchange off, which is what the near field was built from before this
      !! routine existed, sliced to the group. Held to 1e-12 with the default
      !! screen, which both apply by the same test, and with none. STO-3G has
      !! fused L shells and 6-31G* adds Cartesian d, so both shell tables the
      !! quartet loop can run on are covered.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: BRA(8) = [3, 4, 9, 10, 13, 14, 15, 2]
      integer, parameter :: KET1(6) = [1, 2, 5, 6, 7, 8]
      integer, parameter :: KET2(12) = [11, 12, 16, 17, 18, 19, 20, 21, 22, 23, 24, 10]
      integer, parameter :: KET3(3) = [25, 26, 27]
      character(len=6), parameter :: BASES(2) = ["sto-3g", "6-31g*"]
      type(czt_molecule_t) :: mol, piece
      type(error_t) :: err
      type(direct_stats_t) :: stats, block_stats
      integer :: z_all(27), order(29), z(29), ket_first(4), b, ib, k, nao_k, i, j, lo
      character(len=2) :: sym_all(27), sym(29)
      real(dp) :: xyz_all(3, 27), xyz(3, 29), diff
      real(dp), allocatable :: d(:, :), c(:, :), zero_h(:, :), bounds(:, :), full(:, :)
      real(dp), allocatable :: j_block(:, :)
      real(dp) :: tols(2)

      call gly3_water_geometry(z_all, sym_all, xyz_all)
      order = [BRA, KET1, KET2, KET3]
      z = z_all(order)
      sym = sym_all(order)
      xyz = xyz_all(:, order)
      tols = [1.0e-11_dp, NO_SCREENING]

      do ib = 1, size(BASES)
         call build_czt_molecule(z, sym, xyz, trim(BASES(ib)), mol, err)
         if (err%has_error()) then
            call check(error, .false., "supersystem: "//err%get_message())
            return
         end if

         ! Each block's size from the piece on its own, the way FMO counts it.
         allocate (d(mol%nao, mol%nao), source=0.0_dp)
         lo = 0
         do k = 0, 3
            select case (k)
            case (0)
               call build_czt_molecule(z_all(BRA), sym_all(BRA), xyz_all(:, BRA), &
                                       trim(BASES(ib)), piece, err)
            case (1)
               call build_czt_molecule(z_all(KET1), sym_all(KET1), xyz_all(:, KET1), &
                                       trim(BASES(ib)), piece, err)
            case (2)
               call build_czt_molecule(z_all(KET2), sym_all(KET2), xyz_all(:, KET2), &
                                       trim(BASES(ib)), piece, err)
            case (3)
               call build_czt_molecule(z_all(KET3), sym_all(KET3), xyz_all(:, KET3), &
                                       trim(BASES(ib)), piece, err)
            end select
            if (err%has_error()) then
               call check(error, .false., "piece: "//err%get_message())
               return
            end if
            nao_k = piece%nao
            call piece%destroy()
            if (k == 0) then
               b = nao_k
               lo = nao_k
               cycle
            end if
            ket_first(k) = lo
            ! D = 2 C C^T over half the block: symmetric, positive semidefinite
            ! and of a density's size, without an SCF.
            allocate (c(nao_k, max(1, nao_k/2)))
            do j = 1, size(c, 2)
               do i = 1, nao_k
                  c(i, j) = 0.3_dp*sin(1.7_dp*i + 0.9_dp*j*k)/sqrt(real(nao_k, dp)) + &
                            merge(0.8_dp, 0.0_dp, i == j)
               end do
            end do
            d(lo + 1:lo + nao_k, lo + 1:lo + nao_k) = 2.0_dp*matmul(c, transpose(c))
            deallocate (c)
            lo = lo + nao_k
         end do
         ket_first(4) = lo
         call check(error, lo, mol%nao, "the pieces do not add up to the supersystem")
         if (allocated(error)) return

         call schwarz_bounds(mol, bounds, err)
         allocate (zero_h(mol%nao, mol%nao), source=0.0_dp)
         allocate (full(mol%nao, mol%nao))
         do k = 1, size(tols)
            if (.not. err%has_error()) then
               call build_fock_direct(mol, zero_h, d, bounds, full, stats, err, &
                                      screen_tol=tols(k), k_scale=0.0_dp, j_scale=1.0_dp)
            end if
            if (.not. err%has_error()) then
               call coulomb_from_blocks(mol, b, ket_first, d, j_block, err, &
                                        screen_tol=tols(k), stats=block_stats)
            end if
            if (err%has_error()) then
               call check(error, .false., "a build failed: "//err%get_message())
               return
            end if
            diff = maxval(abs(j_block - full(1:b, 1:b)))
            write (*, '(3x,a,a,a,es9.2,a,es10.3,a,i0,a,i0,a,f8.3)') "block Coulomb, ", &
               trim(BASES(ib)), ", screen ", tols(k), ": max |diff| ", diff, &
               ", quartets ", block_stats%quartets_computed, " of ", stats%quartets_computed, &
               ", max |J| ", maxval(abs(j_block))
            call check(error, diff <= 1.0e-12_dp, "the block Coulomb build is not the "// &
                       "bra block of the full one in "//trim(BASES(ib)))
            if (allocated(error)) return
            call check(error, maxval(abs(j_block)) > 1.0_dp, &
                       "the block Coulomb build came back without a field")
            if (allocated(error)) return
         end do
         deallocate (d, zero_h, full, bounds)
         call mol%destroy()
      end do
   end subroutine test_block_coulomb

   subroutine gly3_water_geometry(z, sym, xyz)
      !! `validation/inputs/sample_inputs/gly3_water_pair.xyz`, in Bohr
      integer, intent(out) :: z(27)
      character(len=2), intent(out) :: sym(27)
      real(dp), intent(out) :: xyz(3, 27)
      integer :: i

      z = [7, 6, 6, 8, 1, 1, 1, 1, 7, 6, 6, 8, 1, 1, 1, 7, 6, 6, 8, 1, 1, 1, 8, 1, 8, 1, 1]
      do i = 1, 27
         select case (z(i))
         case (1)
            sym(i) = "H "
         case (6)
            sym(i) = "C "
         case (7)
            sym(i) = "N "
         case default
            sym(i) = "O "
         end select
      end do
      xyz = ANG*reshape([ &
                        0.0171625298_dp, -0.4776667709_dp, -0.0077801388_dp, &
                        1.3251492481_dp, 0.1638239831_dp, 0.0713249069_dp, &
                        1.8818395599_dp, 0.1764813685_dp, 1.4667973423_dp, &
                        1.1563644386_dp, 0.4758564459_dp, 2.4030731780_dp, &
                        2.0041403197_dp, -0.3893217244_dp, -0.6156078332_dp, &
                        1.2933738676_dp, 1.2140808724_dp, -0.2903017566_dp, &
                        -0.6557592247_dp, -0.0682256808_dp, 0.6785523482_dp, &
                        -0.3826962098_dp, -0.2691894812_dp, -0.9506317163_dp, &
                        3.2093591995_dp, -0.0780774266_dp, 1.6702200732_dp, &
                        3.8489825798_dp, -0.0589263473_dp, 2.9842578467_dp, &
                        5.3502343581_dp, -0.0788662970_dp, 2.9476716562_dp, &
                        5.9543074560_dp, -0.1656759551_dp, 1.8893430618_dp, &
                        3.5421254604_dp, 0.8561169960_dp, 3.5393994122_dp, &
                        3.4986665918_dp, -0.9402544817_dp, 3.5643998498_dp, &
                        3.7845901118_dp, -0.3119789206_dp, 0.8286081985_dp, &
                        6.0352251963_dp, 0.0003525130_dp, 4.1282386693_dp, &
                        7.4955375902_dp, -0.0138802141_dp, 4.2014382315_dp, &
                        8.0730347718_dp, 0.0277800836_dp, 5.5909529457_dp, &
                        7.3557278976_dp, 0.0641983810_dp, 6.5759347789_dp, &
                        7.8694940865_dp, -0.9353711779_dp, 3.7021749317_dp, &
                        7.8868335534_dp, 0.8596348618_dp, 3.6344677391_dp, &
                        5.4670886620_dp, 0.0786510231_dp, 5.0034540291_dp, &
                        9.3768940878_dp, 0.0221621974_dp, 5.7818296269_dp, &
                        9.9376629532_dp, -0.0106298905_dp, 4.9380771002_dp, &
                        7.3635223930_dp, -0.3681902986_dp, -0.5795840740_dp, &
                        6.8902239587_dp, -0.3001739023_dp, 0.2496289275_dp, &
                        6.6876213700_dp, -0.2710584500_dp, -1.2503709636_dp &
                        ], [3, 27])
   end subroutine gly3_water_geometry

   subroutine setup(mol, eri, bounds, zero_h, sym, anti, err, basis)
      !! Water in 6-31G, its integrals, and one density of each symmetry
      type(czt_molecule_t), intent(out) :: mol
      real(dp), allocatable, intent(out) :: eri(:, :, :, :), bounds(:, :)
      real(dp), allocatable, intent(out) :: zero_h(:, :), sym(:, :), anti(:, :)
      type(error_t), intent(inout) :: err
      character(len=*), intent(in), optional :: basis
         !! Another basis instead, for a test that needs shells 6-31G lacks.

      real(dp) :: c(3, 3)
      real(dp), allocatable :: m(:, :)
      character(len=:), allocatable :: basis_name
      integer :: n, i, j

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      ! 6-31G rather than STO-3G: several shells per atom and more than one
      ! contraction, so blocks with s1 == s2 and blocks with s1 /= s2 both occur
      ! in quantity. A one-shell-per-atom basis would exercise only some of the
      ! permutation cases and could let a wrong condition pass.
      basis_name = "6-31g"
      if (present(basis)) basis_name = basis
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis_name, mol, err)
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

   subroutine test_hgp_coverage(error)
      !! `HGP_MAX_L` never sends libfint a quartet it does not cover
      !!
      !! The same contract as `test_rotaxis_coverage`, one shell higher: over
      !! every (i, j, i, j) quartet of water in cc-pVTZ, a quartet the dispatch
      !! routes that libfint refuses is an error stop inside a Fock build.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp) :: c(3, 3)
      integer :: ish, jsh, n_on, n_off
      logical :: ours, theirs

      if (.not. HGP_AVAILABLE) then
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
            ours = quartet_on_hgp([ish - 1, jsh - 1, ish - 1, jsh - 1], mol%bas)
            theirs = hgp_libfint_covers([ish - 1, jsh - 1, ish - 1, jsh - 1], mol%bas, mol%nbas)
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
      ! And the limit is one shell above the rotated-axis one: d is on the
      ! path, f is not, so the two constants cannot silently become the same.
      ! The f shell is found rather than assumed to be last -- shells run atom
      ! by atom, so the last one belongs to a hydrogen.
      if (.not. allocated(error)) then
         block
            integer :: n_only_hgp, n_neither
            ! Counted rather than indexed: shells run atom by atom, so which
            ! index carries the f shell is a property of the basis file.
            n_only_hgp = 0
            n_neither = 0
            do ish = 1, mol%nbas
               if (quartet_on_hgp([ish - 1, 0, ish - 1, 0], mol%bas) .and. &
                   .not. quartet_on_rotaxis([ish - 1, 0, ish - 1, 0], mol%bas)) then
                  n_only_hgp = n_only_hgp + 1
               end if
               if (.not. quartet_on_hgp([ish - 1, 0, ish - 1, 0], mol%bas)) then
                  n_neither = n_neither + 1
               end if
            end do
            call check(error, quartet_on_hgp([0, 0, 0, 0], mol%bas), &
                       "the first shell of O must be routed")
            if (.not. allocated(error)) then
               call check(error, n_only_hgp > 0, &
                          "the d shells must be on hgp and off rotaxis")
            end if
            if (.not. allocated(error)) then
               call check(error, n_neither > 0, &
                          "the f shell must be off both paths")
            end if
         end block
      end if
      call mol%destroy()
   end subroutine test_hgp_coverage

   subroutine test_grad_dispatch(error)
      !! The cached gradient dispatch answers exactly as libfint would, per quartet
      !!
      !! `eri_grad_dispatch_t` replaces a per-quartet call to
      !! `libcint_*_grad_supported` with a table looked up by the quartet's four
      !! shell kinds, because that call copies the whole shell table and the
      !! quartet loop is four deep. The replacement is only sound if the two
      !! agree everywhere, so this asks both about every (i, j, i, j) quartet of
      !! a molecule carrying s, p, L, d and f shells and requires them to match
      !! -- not merely that the cache is conservative. A cache that said yes
      !! where libfint says no is an error stop inside a gradient; one that said
      !! no where libfint says yes is a silent loss of the path.
      !!
      !! It also requires that *something* is routed and something is not, since
      !! a table of all-false would agree with nothing being dispatched and pass
      !! a weaker test while quietly disabling the feature.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(eri_grad_dispatch_t) :: disp
      type(error_t) :: err
      real(dp) :: c(3, 3)
      integer :: ish, jsh, shls(4), n_rot, n_hgp, n_off
      logical :: theirs

      if (.not. (ROTAXIS_AVAILABLE .and. HGP_AVAILABLE)) then
         call check(error, .true.)
         return
      end if

      ! cc-pVTZ for the f shells the gradient sets do not cover, and a geometry
      ! with no symmetry so nothing is zero by accident.
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.4_dp, 1.1_dp, 0.0_dp, -1.4_dp, 1.1_dp], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "cc-pvtz", mol, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if

      call build_eri_grad_dispatch(mol%bas, mol%nbas, disp)

      n_rot = 0
      n_hgp = 0
      n_off = 0
      do ish = 1, mol%nbas
         do jsh = 1, mol%nbas
            shls = [ish - 1, jsh - 1, ish - 1, jsh - 1]
            theirs = rotaxis_grad_cached(disp, shls)
            if (theirs .neqv. rotaxis_grad_libfint(shls, mol%bas, mol%nbas)) then
               call check(error, .false., "the cached rotated-axis gradient dispatch "// &
                          "disagrees with libfint")
               call mol%destroy()
               return
            end if
            if (theirs) n_rot = n_rot + 1

            theirs = hgp_grad_cached(disp, shls)
            if (theirs .neqv. hgp_grad_libfint(shls, mol%bas, mol%nbas)) then
               call check(error, .false., "the cached Head-Gordon-Pople gradient dispatch "// &
                          "disagrees with libfint")
               call mol%destroy()
               return
            end if
            if (theirs) n_hgp = n_hgp + 1
            if (.not. (rotaxis_grad_cached(disp, shls) .or. hgp_grad_cached(disp, shls))) then
               n_off = n_off + 1
            end if
         end do
      end do
      call mol%destroy()
      call disp%destroy()

      if (.not. allocated(error)) then
         call check(error, n_rot > 0, "no quartet reached the rotated-axis gradient")
      end if
      if (.not. allocated(error)) then
         call check(error, n_hgp > 0, "no quartet reached the Head-Gordon-Pople gradient")
      end if
      if (.not. allocated(error)) then
         call check(error, n_off > 0, "cc-pVTZ has f shells, which neither gradient covers")
      end if
   end subroutine test_grad_dispatch

   subroutine test_hgp_fock(error)
      !! The Head-Gordon-Pople and hybrid paths agree with Rys through a Fock build
      !!
      !! Water/6-31G* rather than 6-31G, because d is what this path adds and a
      !! basis without one would leave the new kernels untouched. Under
      !! `hybrid` the same molecule splits: the s/p quartets go rotated-axis
      !! and the d-touching ones Head-Gordon-Pople, so one build covers the
      !! branch that chooses between them. Different algorithms, so agreement
      !! is to 1e-10 rather than bitwise.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :), rys(:, :), hgp(:, :), hybrid(:, :)

      if (.not. HGP_AVAILABLE) then
         call set_eri_path("hgp", err)
         call check(error, err%has_error(), &
                    "a build without the path must refuse to be asked for it")
         return
      end if

      call setup(mol, eri, bounds, zero_h, sym, anti, err, basis="6-31g_st_")
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if

      allocate (rys(mol%nao, mol%nao), hgp(mol%nao, mol%nao), hybrid(mol%nao, mol%nao))
      call set_eri_path("rys", err)
      call build_fock_direct(mol, zero_h, sym, bounds, rys, stats, err, &
                             screen_tol=NO_SCREENING)
      if (.not. err%has_error()) then
         call set_eri_path("hgp", err)
         call check(error, eri_path_name() == "hgp", "the path did not switch")
         if (.not. allocated(error)) then
            call build_fock_direct(mol, zero_h, sym, bounds, hgp, stats, err, &
                                   screen_tol=NO_SCREENING)
         end if
      end if
      if (.not. err%has_error() .and. .not. allocated(error)) then
         call set_eri_path("hybrid", err)
         call check(error, eri_path_name() == "hybrid", "the path did not switch")
         if (.not. allocated(error)) then
            call build_fock_direct(mol, zero_h, sym, bounds, hybrid, stats, err, &
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

      call check(error, maxval(abs(hgp - rys)) < 1.0e-10_dp, &
                 "the Head-Gordon-Pople and Rys Fock matrices disagree")
      if (.not. allocated(error)) then
         call check(error, maxval(abs(hybrid - rys)) < 1.0e-10_dp, &
                    "the hybrid and Rys Fock matrices disagree")
      end if
   end subroutine test_hgp_fock

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

   subroutine test_nosym_scales(error)
      !! `k_scale`, `j_scale` and `omega` mean on the general build what they
      !! mean on the fast one
      !!
      !! A symmetric density is the only place the two builds can be compared,
      !! and it is enough: the coefficients multiply the same six contributions
      !! either way, so agreeing here pins them for any density. Both the
      !! full-range hybrid pass and the attenuated exchange-only one, since
      !! those are the two a range-separated response operator makes.
      !!
      !! The attenuated result is also required to *differ* from the
      !! unattenuated one. `omega` reaches libcint through a slot in `env`, so a
      !! build that forgot to copy the environment returns full-range integrals
      !! scaled by the long-range coefficient, silently, and would match the
      !! fast build's coefficients while answering the wrong operator.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :)
      real(dp), allocatable :: dens(:, :, :)
      real(dp), allocatable :: fast(:, :, :), general(:, :, :)
      real(dp), allocatable :: fast_lr(:, :, :), general_lr(:, :, :)

      call setup(mol, eri, bounds, zero_h, sym, anti, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if

      allocate (dens(mol%nao, mol%nao, 1))
      dens(:, :, 1) = sym

      call build_fock_direct_many(mol, zero_h, dens, bounds, fast, stats, err, &
                                  screen_tol=NO_SCREENING, k_scale=CAM_K_FULL)
      if (.not. err%has_error()) &
         call build_fock_direct_nosym(mol, zero_h, dens, bounds, general, stats, err, &
                                      screen_tol=NO_SCREENING, k_scale=CAM_K_FULL)
      if (.not. err%has_error()) &
         call build_fock_direct_many(mol, zero_h, dens, bounds, fast_lr, stats, err, &
                                     screen_tol=NO_SCREENING, k_scale=CAM_K_LR, &
                                     j_scale=0.0_dp, omega=CAM_OMEGA)
      if (.not. err%has_error()) &
         call build_fock_direct_nosym(mol, zero_h, dens, bounds, general_lr, stats, err, &
                                      screen_tol=NO_SCREENING, k_scale=CAM_K_LR, &
                                      j_scale=0.0_dp, omega=CAM_OMEGA)
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "a build failed: "//err%get_message())
         return
      end if

      call check(error, maxval(abs(general(:, :, 1) - fast(:, :, 1))) < 1.0e-11_dp, &
                 "the general build does not scale exact exchange the way the fast "// &
                 "build does")
      if (allocated(error)) return
      call check(error, maxval(abs(general_lr(:, :, 1) - fast_lr(:, :, 1))) < 1.0e-11_dp, &
                 "the general build disagrees with the fast one on an attenuated "// &
                 "exchange-only pass")
      if (allocated(error)) return
      ! The long-range pass is a fraction of the full-range one at this omega,
      ! so the two are nowhere near each other even after the coefficients.
      call check(error, maxval(abs(general_lr(:, :, 1) - (CAM_K_LR/CAM_K_FULL) &
                                   *general(:, :, 1))) > 1.0e-3_dp, &
                 "an omega pass through build_fock_direct_nosym is not attenuated")
   end subroutine test_nosym_scales

   subroutine test_stored_scales(error)
      !! `k_scale` and `j_scale` mean the same on the stored tensor as direct
      !!
      !! `build_fock` gained `j_scale` so the in-core route is not quietly
      !! Coulomb-carrying where its direct twin is not, and nothing on the
      !! response path reaches it: the core always asks for `direct = .true.`.
      !! An argument no caller exercises is an argument nobody finds out is
      !! wrong, so it is exercised here instead -- against
      !! `build_fock_direct`, which is the same six contributions computed a
      !! completely different way.
      !!
      !! Both coefficients at once and neither equal to the other or to one,
      !! and then the `j_scale = 0` case on its own: that is the shape the
      !! triplet response and a long-range exchange pass both use, and it is
      !! the branch where the Coulomb accumulation is now skipped rather than
      !! computed and multiplied by zero.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: J_SCALE = 0.63_dp
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :)
      real(dp), allocatable :: stored(:, :), direct(:, :)
      real(dp) :: scale

      call setup(mol, eri, bounds, zero_h, sym, anti, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if

      allocate (stored(mol%nao, mol%nao), direct(mol%nao, mol%nao))

      call build_fock(zero_h, eri, sym, stored, k_scale=CAM_K_FULL, j_scale=J_SCALE)
      call build_fock_direct(mol, zero_h, sym, bounds, direct, stats, err, &
                             screen_tol=NO_SCREENING, k_scale=CAM_K_FULL, &
                             j_scale=J_SCALE)
      if (err%has_error()) then
         call mol%destroy()
         call check(error, .false., "the direct build failed: "//err%get_message())
         return
      end if
      scale = maxval(abs(direct))
      call check(error, scale > 1.0e-2_dp, "the scaled build is empty, so the "// &
                 "comparison below would hold for any pair of coefficients")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      call check(error, maxval(abs(stored - direct)) < 1.0e-11_dp, &
                 "the stored-tensor build does not scale the two-electron terms "// &
                 "the way the direct build does")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      ! Exchange alone. Checked against the direct build rather than against
      ! the scaled result above, so a Coulomb term that survived the skip has
      ! nowhere to hide.
      call build_fock(zero_h, eri, sym, stored, k_scale=CAM_K_FULL, j_scale=0.0_dp)
      call build_fock_direct(mol, zero_h, sym, bounds, direct, stats, err, &
                             screen_tol=NO_SCREENING, k_scale=CAM_K_FULL, &
                             j_scale=0.0_dp)
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "the exchange-only direct build failed: "// &
                    err%get_message())
         return
      end if
      call check(error, maxval(abs(stored - direct)) < 1.0e-11_dp, &
                 "the stored-tensor build with j_scale = 0 does not match the "// &
                 "direct build's exchange-only pass")
      if (allocated(error)) return
      call check(error, maxval(abs(stored)) > 1.0e-2_dp, &
                 "the exchange-only build is empty, so matching it says nothing")
   end subroutine test_stored_scales

   subroutine test_nosym_antisymmetric_scaled(error)
      !! A scaled `A - B` pass: no Coulomb, and exchange at the coefficient asked
      !!
      !! The reference is written out here rather than taken from `build_fock`,
      !! because what is under test is the *coefficients* and `build_fock`
      !! applies its own. `J` is formed too, and required to vanish, which is
      !! what says the Coulomb term dropped out of its own accord rather than
      !! because `j_scale` was zero -- it is left at its default of one here on
      !! purpose.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: K_SCALE = 0.37_dp
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :)
      real(dp), allocatable :: dens(:, :, :), general(:, :, :)
      real(dp), allocatable :: j_ref(:, :), k_ref(:, :)
      integer :: n, a, b, c, d

      ! STO-3G, so the explicit n^4 reference below is seven functions rather
      ! than thirteen and the loop is a test and not a benchmark.
      call setup(mol, eri, bounds, zero_h, sym, anti, err, basis="sto-3g")
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if

      n = mol%nao
      allocate (dens(n, n, 1), j_ref(n, n), k_ref(n, n))
      dens(:, :, 1) = anti
      j_ref = 0.0_dp
      k_ref = 0.0_dp
      do d = 1, n
         do c = 1, n
            do b = 1, n
               do a = 1, n
                  j_ref(a, b) = j_ref(a, b) + eri(a, b, c, d)*anti(c, d)
                  k_ref(a, c) = k_ref(a, c) + eri(a, b, c, d)*anti(b, d)
               end do
            end do
         end do
      end do

      call build_fock_direct_nosym(mol, zero_h, dens, bounds, general, stats, err, &
                                   screen_tol=NO_SCREENING, k_scale=K_SCALE)
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "the general build failed: "//err%get_message())
         return
      end if

      call check(error, maxval(abs(j_ref)) < 1.0e-12_dp, &
                 "the Coulomb contraction of an antisymmetric density is not zero, "// &
                 "so this test cannot say what the general build left out")
      if (allocated(error)) return
      call check(error, maxval(abs(general(:, :, 1) + 0.5_dp*K_SCALE*k_ref)) < 1.0e-11_dp, &
                 "the general build does not reproduce a scaled exchange contraction "// &
                 "on an antisymmetric density")
   end subroutine test_nosym_antisymmetric_scaled

   subroutine uhf_reference(eri, zero_h, d_alpha, d_beta, k, g_a, g_b)
      !! `G_s = J[D_a + D_b] - k K[D_s]` from the stored tensor, by definition
      !!
      !! Two calls to `build_fock` per spin, which is the plain `n^4`
      !! contraction with no symmetry assumption in it: one with the Coulomb
      !! term alone on the total density, one with the exchange term alone on
      !! that spin's. `build_fock` carries `K/2`, so `k_scale = 2 k` is what
      !! makes the second call the **full** same-spin exchange an unrestricted
      !! build wants.
      real(dp), intent(in) :: eri(:, :, :, :), zero_h(:, :)
      real(dp), intent(in) :: d_alpha(:, :), d_beta(:, :)
      real(dp), intent(in) :: k
      real(dp), allocatable, intent(out) :: g_a(:, :), g_b(:, :)

      real(dp), allocatable :: coul(:, :), exch(:, :)
      integer :: n

      n = size(zero_h, 1)
      allocate (g_a(n, n), g_b(n, n), coul(n, n), exch(n, n))
      call build_fock(zero_h, eri, d_alpha + d_beta, coul, k_scale=0.0_dp)
      call build_fock(zero_h, eri, d_alpha, exch, k_scale=2.0_dp*k, j_scale=0.0_dp)
      g_a = coul + exch
      call build_fock(zero_h, eri, d_beta, exch, k_scale=2.0_dp*k, j_scale=0.0_dp)
      g_b = coul + exch
      deallocate (coul, exch)
   end subroutine uhf_reference

   subroutine test_uhf_many_closed_shell(error)
      !! Half the closed-shell density in each spin gives the closed-shell build
      !!
      !! This is the convention, stated as a test: the unrestricted build takes
      !! the **true** spin densities and returns full same-spin exchange, so
      !! `D_a = D_b = D/2` has to reproduce `J[D] - K[D]/2` -- the matrix the
      !! restricted batch returns for `D` -- element for element. Get the
      !! factor wrong in either direction and this fails by exactly two.
      !!
      !! A scaled batch of three, so the set index is exercised and a routine
      !! that quietly contracted every set against the first would show up.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: SETS = 3
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :)
      real(dp), allocatable :: batch(:, :, :), half(:, :, :)
      real(dp), allocatable :: closed(:, :, :), ga(:, :, :), gb(:, :, :)
      real(dp) :: worst
      integer :: m

      call setup(mol, eri, bounds, zero_h, sym, anti, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if
      allocate (batch(mol%nao, mol%nao, SETS), half(mol%nao, mol%nao, SETS))
      do m = 1, SETS
         batch(:, :, m) = (0.7_dp*m - 1.1_dp)*sym
         half(:, :, m) = 0.5_dp*batch(:, :, m)
      end do

      call build_fock_direct_many(mol, zero_h, batch, bounds, closed, stats, err, &
                                  screen_tol=NO_SCREENING, k_scale=CAM_K_FULL)
      if (.not. err%has_error()) then
         call build_fock_direct_uhf_many(mol, zero_h, half, half, bounds, ga, gb, &
                                         stats, err, screen_tol=NO_SCREENING, &
                                         k_scale=CAM_K_FULL)
      end if
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "a build failed: "//err%get_message())
         return
      end if

      worst = 0.0_dp
      do m = 1, SETS
         worst = max(worst, maxval(abs(ga(:, :, m) - closed(:, :, m))))
         worst = max(worst, maxval(abs(gb(:, :, m) - closed(:, :, m))))
      end do
      call check(error, worst < 1.0e-12_dp, "the unrestricted batch on two equal "// &
                 "half densities is not the closed-shell build")
   end subroutine test_uhf_many_closed_shell

   subroutine test_uhf_many_explicit(error)
      !! A genuinely spin-polarised pair, against the four-index contraction
      !!
      !! Equal spin densities cannot separate `J[D_a + D_b]` from `2 J[D_a]`,
      !! nor same-spin exchange from any other combination, so the case above
      !! passes for several wrong builds. Here the two spins are different
      !! matrices and every term is pinned on its own, against `build_fock` --
      !! which assumes no permutational symmetry and is correct by
      !! construction.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: SETS = 2
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :)
      real(dp), allocatable :: da(:, :, :), db(:, :, :), ga(:, :, :), gb(:, :, :)
      real(dp), allocatable :: ref_a(:, :), ref_b(:, :), single_a(:, :), single_b(:, :)
      real(dp) :: worst
      integer :: m, n, i, j

      call setup(mol, eri, bounds, zero_h, sym, anti, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if
      n = mol%nao
      allocate (da(n, n, SETS), db(n, n, SETS), single_a(n, n), single_b(n, n))
      ! Two symmetric matrices that are not multiples of each other, so the
      ! alpha and beta halves of every term are distinguishable.
      do j = 1, n
         do i = 1, n
            da(i, j, 1) = 0.3_dp*sym(i, j) + 0.05_dp*cos(real(2*i + j, dp))
            db(i, j, 1) = -0.4_dp*sym(i, j) + 0.02_dp*sin(real(i - 3*j, dp))
         end do
      end do
      da(:, :, 1) = 0.5_dp*(da(:, :, 1) + transpose(da(:, :, 1)))
      db(:, :, 1) = 0.5_dp*(db(:, :, 1) + transpose(db(:, :, 1)))
      da(:, :, 2) = db(:, :, 1)
      db(:, :, 2) = da(:, :, 1)

      call build_fock_direct_uhf_many(mol, zero_h, da, db, bounds, ga, gb, stats, err, &
                                      screen_tol=NO_SCREENING, k_scale=CAM_K_FULL)
      ! And the single-density build it is the batch of: the two loops are
      ! separate copies and this is what keeps them in step.
      if (.not. err%has_error()) then
         call build_fock_direct_uhf(mol, zero_h, da(:, :, 1), db(:, :, 1), bounds, &
                                    single_a, single_b, stats, err, &
                                    screen_tol=NO_SCREENING, k_scale=CAM_K_FULL)
      end if
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "a build failed: "//err%get_message())
         return
      end if

      worst = 0.0_dp
      do m = 1, SETS
         call uhf_reference(eri, zero_h, da(:, :, m), db(:, :, m), CAM_K_FULL, &
                            ref_a, ref_b)
         worst = max(worst, maxval(abs(ga(:, :, m) - ref_a)))
         worst = max(worst, maxval(abs(gb(:, :, m) - ref_b)))
         deallocate (ref_a, ref_b)
      end do
      call check(error, worst < 1.0e-10_dp, "the unrestricted batch disagrees with "// &
                 "an explicit four-index contraction")
      if (allocated(error)) return

      worst = max(maxval(abs(ga(:, :, 1) - single_a)), &
                  maxval(abs(gb(:, :, 1) - single_b)))
      call check(error, worst < 1.0e-10_dp, "the unrestricted batch disagrees with "// &
                 "the single-density unrestricted build it is a batch of")
   end subroutine test_uhf_many_explicit

   subroutine test_uhf_many_antisymmetric(error)
      !! Announced antisymmetric: no Coulomb, and the output antisymmetrised
      !!
      !! `(A - B)` hands this build an antisymmetric pair. What has to come
      !! back is exchange alone -- the Coulomb term vanishes for an
      !! antisymmetric density and the fold that stands in for the remaining
      !! permutations has to change sign with it. Checked against the
      !! four-index contraction, which needs no announcement, and the result is
      !! checked to be antisymmetric rather than assumed so.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: sym(:, :), anti(:, :)
      real(dp), allocatable :: da(:, :, :), db(:, :, :), ga(:, :, :), gb(:, :, :)
      real(dp), allocatable :: ref_a(:, :), ref_b(:, :)
      real(dp) :: worst
      integer :: n

      call setup(mol, eri, bounds, zero_h, sym, anti, err)
      if (err%has_error()) then
         call check(error, .false., "setup failed: "//err%get_message())
         return
      end if
      n = mol%nao
      allocate (da(n, n, 1), db(n, n, 1))
      da(:, :, 1) = anti
      db(:, :, 1) = -0.6_dp*anti

      call build_fock_direct_uhf_many(mol, zero_h, da, db, bounds, ga, gb, stats, err, &
                                      screen_tol=NO_SCREENING, k_scale=CAM_K_FULL, &
                                      antisymmetric=.true.)
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "the antisymmetric unrestricted build failed: "// &
                    err%get_message())
         return
      end if

      call uhf_reference(eri, zero_h, da(:, :, 1), db(:, :, 1), CAM_K_FULL, ref_a, ref_b)
      worst = max(maxval(abs(ga(:, :, 1) - ref_a)), maxval(abs(gb(:, :, 1) - ref_b)))
      call check(error, worst < 1.0e-10_dp, "the announced-antisymmetric unrestricted "// &
                 "batch disagrees with an explicit contraction")
      if (allocated(error)) return
      worst = max(maxval(abs(ga(:, :, 1) + transpose(ga(:, :, 1)))), &
                  maxval(abs(gb(:, :, 1) + transpose(gb(:, :, 1)))))
      call check(error, worst < 1.0e-12_dp, "an antisymmetric pair came back with a "// &
                 "symmetric part")
   end subroutine test_uhf_many_antisymmetric

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
