!! EFMO across covalent bonds, at the limits where its total is known
module test_mqc_czt_efmo_covalent
   !! A partition that cuts a covalent bond is run with `bond_breaking = "afo"`:
   !! every monomer, its potential and every near group is solved with the cut
   !! bonds detached by adjusted frozen orbitals, through FMO's own group
   !! assembly, and the detached atom's nucleus split so every fragment is
   !! neutral.
   !!
   !! What is asserted is what the expression guarantees whatever the
   !! partition did:
   !!
   !!   1. **Level equal to the fragment count, cutoff huge, is the
   !!      unfragmented RHF energy.** The in-vacuo series telescopes to the
   !!      supersystem's SCF, which holds every bond whole, and the induction
   !!      series telescopes to `E_pol^total` and cancels it. Checked on
   !!      propane cut once and cut twice, on propane numbered so one carbon is
   !!      the detached end of both bonds, and on methylcyclopropane cut at the
   !!      exocyclic bond.
   !!   2. **A covalently joined pair is quantum at any cutoff.** With a cutoff
   !!      far inside every separation, the two bonded pairs of propane are
   !!      still dimers, and a two-fragment cut is still the whole molecule.
   !!   3. **Every fragment's potential is neutral**, which is what splitting
   !!      the nucleus is for.
   !!   4. Correlation across a cut is refused by name.
   !!
   !! 6-31G throughout: s and p only, so the Cartesian molecules EFMO builds
   !! and the reference here are the same basis.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_czt_efmo, only: efmo_options_t, efmo_result_t, run_efmo, EFMO_CORR_MP2, &
                           efmo_pair_contribution
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_atomic_guess, only: build_restricted_guess
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_physical_constants, only: ANGSTROM_TO_BOHR
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_czt_efmo_covalent_tests

   real(dp), parameter :: ANG = ANGSTROM_TO_BOHR

   character(len=*), parameter :: BASIS = "6-31g"

   real(dp), parameter :: TOL = 1.0e-8_dp
      !! The identities are exact. What this clears is SCF convergence, 1e-10
      !! per SCF over the dozen or so a three-fragment run takes, and the
      !! frozen-virtual shift's back-transform noise, about 1e-13.

contains

   subroutine collect_mqc_czt_efmo_covalent_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("efmo_afo_propane_one_cut_is_the_molecule", &
                               test_propane_one_cut), &
                  new_unittest("efmo_afo_propane_two_cuts_full_level_is_the_molecule", &
                               test_propane_two_cuts), &
                  new_unittest("efmo_afo_doubly_detached_atom_full_level_is_the_molecule", &
                               test_propane_middle_first), &
                  new_unittest("efmo_afo_butane_four_pieces_full_level_is_the_molecule", &
                               test_butane_four), &
                  new_unittest("efmo_afo_methylcyclopropane_exocyclic_is_the_molecule", &
                               test_methylcyclopropane), &
                  new_unittest("efmo_afo_joined_pairs_are_quantum_at_any_cutoff", &
                               test_joined_pairs_quantum), &
                  new_unittest("efmo_afo_cut_potentials_are_neutral", test_neutral), &
                  new_unittest("efmo_afo_far_water_pair_against_the_quantum_pair", &
                               test_far_water), &
                  new_unittest("efmo_afo_correlation_is_refused", test_refuse_mp2) &
                  ]
   end subroutine collect_mqc_czt_efmo_covalent_tests

   subroutine test_propane_one_cut(error)
      !! Propane split at one C-C bond: methyl and ethyl
      !!
      !! Two fragments at level two is the full expansion, so the answer is
      !! the molecule's own energy.
      type(error_type), allocatable, intent(out) :: error
      integer :: z(11)
      character(len=2) :: sym(11)
      real(dp) :: xyz(3, 11)

      call propane(z, sym, xyz)
      call full_level_is_whole(z, sym, xyz, [1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2], 2, 1, error)
   end subroutine test_propane_one_cut

   subroutine test_propane_two_cuts(error)
      !! Propane in three pieces, CH3 | CH2 | CH3, at level three
      !!
      !! The dimer of the two ends is the interesting group: they are not
      !! bonded to each other, and it carries the boundary of each to the
      !! middle.
      type(error_type), allocatable, intent(out) :: error
      integer :: z(11)
      character(len=2) :: sym(11)
      real(dp) :: xyz(3, 11)

      call propane(z, sym, xyz)
      call full_level_is_whole(z, sym, xyz, [1, 2, 3, 1, 1, 1, 2, 2, 3, 3, 3], 3, 2, error)
   end subroutine test_propane_two_cuts

   subroutine test_propane_middle_first(error)
      !! The same three pieces, numbered so the middle carbon is detached twice
      !!
      !! The bond-detached atom of a cut is its lower-numbered end, so with the
      !! middle carbon first both end fragments carry a ghost of the *same*
      !! atom, and their dimer carries that ghost once with `+2` on it. Two
      !! fragments sharing a centre are joined, so that dimer must be quantum.
      type(error_type), allocatable, intent(out) :: error
      integer :: z(11)
      character(len=2) :: sym(11)
      real(dp) :: xyz(3, 11)

      call propane_middle_first(z, sym, xyz)
      call full_level_is_whole(z, sym, xyz, [1, 2, 3, 1, 1, 2, 2, 2, 3, 3, 3], 3, 2, error)
   end subroutine test_propane_middle_first

   subroutine test_butane_four(error)
      !! Anti butane cut at all three C-C bonds, one carbon per fragment, level 4
      !!
      !! Fifteen groups, and the two middle fragments are each the detached end
      !! of one bond and the attached end of another.
      type(error_type), allocatable, intent(out) :: error
      integer :: z(14)
      character(len=2) :: sym(14)
      real(dp) :: xyz(3, 14)

      call butane(z, sym, xyz)
      call full_level_is_whole(z, sym, xyz, [1, 2, 3, 4, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4], 4, 3, &
                               error)
   end subroutine test_butane_four

   subroutine test_methylcyclopropane(error)
      !! Methylcyclopropane, cut at the one bond not in the ring
      !!
      !! Cyclopropane itself cannot be used: any partition of a three-ring cuts
      !! two bonds between the same pair of fragments, which is refused.
      type(error_type), allocatable, intent(out) :: error
      integer :: z(12)
      character(len=2) :: sym(12)
      real(dp) :: xyz(3, 12)

      call methylcyclopropane(z, sym, xyz)
      call full_level_is_whole(z, sym, xyz, [1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 2], 2, 1, &
                               error)
   end subroutine test_methylcyclopropane

   subroutine test_joined_pairs_quantum(error)
      !! A cutoff far inside every separation still leaves bonded pairs quantum
      !!
      !! At `R_cut = 0.1` no pair of propane's fragments is near by distance --
      !! a C-C bond is about 0.45 in these units -- so without the rule every
      !! pair would be effective and the bond would be described by two
      !! multipole expansions sharing a centre.
      type(error_type), allocatable, intent(out) :: error
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      type(error_t) :: err
      integer :: z(11)
      character(len=2) :: sym(11)
      real(dp) :: xyz(3, 11), whole

      call propane(z, sym, xyz)
      call settings(opts)
      opts%rcut = 0.1_dp

      ! Two fragments: the one pair is bonded, so the answer is the molecule.
      call run_efmo(z, sym, xyz, [1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2], [0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, res%n_qm_pairs, 1, message="the bonded pair was not quantum")
      if (allocated(error)) return
      call whole_rhf(z, sym, xyz, whole, err)
      call check(error,.not. err%has_error(), "the whole-molecule RHF failed")
      if (allocated(error)) return
      write (*, *) "   two fragments at R_cut 0.1: EFMO - RHF =", res%energy - whole
      call check(error, res%energy, whole, thr=TOL, &
                 message="a bonded pair made quantum by rule is not the molecule")
      if (allocated(error)) return

      ! Three fragments. The ends are not bonded to each other, but the third
      ! carries a ghost of the middle carbon, a bond length from the first
      ! carbon: joined, so quantum too. Effective, that pair came out at -0.34
      ! Hartree against +0.27 quantum.
      call run_efmo(z, sym, xyz, [1, 2, 3, 1, 1, 1, 2, 2, 3, 3, 3], [0, 0, 0], opts, res, &
                    err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, res%n_qm_pairs, 3, message="a joined pair of propane is not quantum")
      if (allocated(error)) return

      ! Butane in four: the first and last share no centre and no bond, so
      ! that one pair alone is effective; its two numbers are reported.
      pair_one_four: block
         integer :: zb(14)
         character(len=2) :: sb(14)
         real(dp) :: xb(3, 14), efp_value, qm_value
         call butane(zb, sb, xb)
         call run_efmo(zb, sb, xb, [1, 2, 3, 4, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4], [0, 0, 0, 0], &
                       opts, res, err)
         call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
         if (allocated(error)) return
         call check(error, res%n_qm_pairs, 5, message="butane's joined pairs are not quantum")
         if (allocated(error)) return
         call check(error, res%n_efp_pairs, 1, message="butane's end pair is not effective")
         if (allocated(error)) return
         call far_and_near(zb, sb, xb, [1, 2, 3, 4, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4], [1, 4], &
                           "afo", efp_value, qm_value, err)
         call check(error,.not. err%has_error(), "the end pair failed: "// &
                    err%get_full_trace())
         if (allocated(error)) return
         write (*, *) "   butane end pair (1,4): efp", efp_value, " qm", qm_value
      end block pair_one_four
   end subroutine test_joined_pairs_quantum

   subroutine test_neutral(error)
      !! Every fragment's potential sums to the charge it was given
      !!
      !! Propane in three pieces, so the middle fragment both owns a detached
      !! atom (`Z-1`, hybrid empty) and carries a ghost of another (`+1`,
      !! hybrid full). With whole nuclei the two ends would be at +1 and -1
      !! and the middle at zero only by cancellation.
      type(error_type), allocatable, intent(out) :: error
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      type(error_t) :: err
      integer :: z(11), k
      character(len=2) :: sym(11)
      real(dp) :: xyz(3, 11)

      call propane(z, sym, xyz)
      call settings(opts)
      opts%level = 1
      call run_efmo(z, sym, xyz, [1, 2, 3, 1, 1, 1, 2, 2, 3, 3, 3], [0, 0, 0], opts, res, &
                    err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return
      write (*, *) "   potential charges =", res%potential_charge
      do k = 1, 3
         call check(error, abs(res%potential_charge(k)) < 1.0e-6_dp, &
                    "a cut fragment's potential is not neutral")
         if (allocated(error)) return
      end do
   end subroutine test_neutral

   subroutine test_far_water(error)
      !! A water beyond the cutoff, against a cut fragment and against the same
      !! molecule uncut
      !!
      !! Each far pair's four effective-fragment terms are compared with the
      !! same pair solved as a quantum dimer, `E_IJ^0 - E_I^0 - E_J^0 -
      !! E_IJ^pol`, which is what the pair contributes when it is near. The
      !! uncut molecule is the baseline: the effective-fragment approximation
      !! has an error of its own there, and what a cut must not do is make it
      !! materially worse.
      type(error_type), allocatable, intent(out) :: error
      integer :: z(14), k
      character(len=2) :: sym(14)
      real(dp) :: xyz(3, 14), gap(3), efp_w, qm_w, efp_c, qm_c, worst_whole, worst_cut
      real(dp), parameter :: SEPARATION(3) = [3.5_dp, 4.5_dp, 6.0_dp]
      type(error_t) :: err

      worst_whole = 0.0_dp
      worst_cut = 0.0_dp
      do k = 1, size(SEPARATION)
         call propane_water(SEPARATION(k), z, sym, xyz)
         ! Uncut: propane is fragment 1, the water fragment 2.
         call far_and_near(z, sym, xyz, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2], &
                           [1, 2], "none", efp_w, qm_w, err)
         call check(error,.not. err%has_error(), "uncut run failed: "// &
                    err%get_full_trace())
         if (allocated(error)) return
         ! Cut: the methyl facing the water is fragment 1, the ethyl 2.
         call far_and_near(z, sym, xyz, [1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3], &
                           [1, 3], "afo", efp_c, qm_c, err)
         call check(error,.not. err%has_error(), "cut run failed: "// &
                    err%get_full_trace())
         if (allocated(error)) return
         gap = [efp_w - qm_w, efp_c - qm_c, 0.0_dp]
         write (*, "(a,f5.2,a,4es13.4,a,2es11.3)") "   water at ", SEPARATION(k), &
            " A: whole efp/qm, cut efp/qm =", efp_w, qm_w, efp_c, qm_c, &
            "  errors whole, cut:", gap(1), gap(2)
         worst_whole = max(worst_whole, abs(gap(1)))
         worst_cut = max(worst_cut, abs(gap(2)))
      end do
      write (*, *) "   worst |efp - qm|: whole", worst_whole, " cut", worst_cut
      call check(error, worst_cut < 3.0_dp*worst_whole + 2.0e-4_dp, &
                 "a far pair with a cut fragment is much worse than the same pair uncut")
   end subroutine test_far_water

   subroutine far_and_near(z, sym, xyz, owner, pair, bond_breaking, efp_value, qm_value, err)
      !! One pair's contribution computed effective and computed quantum
      integer, intent(in) :: z(:), owner(:), pair(2)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: xyz(:, :)
      character(len=*), intent(in) :: bond_breaking
      real(dp), intent(out) :: efp_value, qm_value
      type(error_t), intent(inout) :: err

      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      integer, allocatable :: charges(:)
      integer :: k

      efp_value = 0.0_dp
      qm_value = 0.0_dp
      allocate (charges(maxval(owner)), source=0)
      call settings(opts)
      opts%bond_breaking = bond_breaking
      opts%induction_damping = 0.1_dp
      opts%rcut = 0.1_dp
      call run_efmo(z, sym, xyz, owner, charges, opts, res, err)
      if (err%has_error()) return
      do k = 1, size(res%pairs)
         if (res%pairs(k)%i == pair(1) .and. res%pairs(k)%j == pair(2)) then
            efp_value = efmo_pair_contribution(res, k)
         end if
      end do
      opts%rcut = 1.0e6_dp
      call run_efmo(z, sym, xyz, owner, charges, opts, res, err)
      if (err%has_error()) return
      do k = 1, size(res%pairs)
         if (res%pairs(k)%i == pair(1) .and. res%pairs(k)%j == pair(2)) then
            qm_value = efmo_pair_contribution(res, k)
         end if
      end do
   end subroutine far_and_near

   subroutine propane_water(distance, z, sym, xyz)
      !! Propane, and a water on the far side of its first methyl
      !!
      !! The oxygen sits `distance` Angstrom beyond the first carbon along the
      !! C-C axis, its hydrogens pointing away.
      real(dp), intent(in) :: distance
      integer, intent(out) :: z(14)
      character(len=2), intent(out) :: sym(14)
      real(dp), intent(out) :: xyz(3, 14)

      integer :: zp(11)
      character(len=2) :: sp(11)
      real(dp) :: xp(3, 11), x0

      call propane(zp, sp, xp)
      z(1:11) = zp
      sym(1:11) = sp
      xyz(:, 1:11) = xp
      z(12:14) = [8, 1, 1]
      sym(12:14) = ["O ", "H ", "H "]
      x0 = 1.5260_dp + distance
      xyz(:, 12) = [x0, 0.0_dp, 0.0_dp]*ANG
      xyz(:, 13) = [x0 + 0.5686_dp, 0.7725_dp, 0.0_dp]*ANG
      xyz(:, 14) = [x0 + 0.5686_dp, -0.7725_dp, 0.0_dp]*ANG
   end subroutine propane_water

   subroutine test_refuse_mp2(error)
      !! MP2 across a cut would correlate into the frozen virtual; refused
      type(error_type), allocatable, intent(out) :: error
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      type(error_t) :: err
      integer :: z(11)
      character(len=2) :: sym(11)
      real(dp) :: xyz(3, 11)

      call propane(z, sym, xyz)
      call settings(opts)
      opts%correlation = EFMO_CORR_MP2
      call run_efmo(z, sym, xyz, [1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2], [0, 0], opts, res, err)
      call check(error, err%has_error(), "MP2 across a detached bond should be refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "Hartree-Fock") > 0, &
                 "the refusal should say only Hartree-Fock runs: "//err%get_message())
   end subroutine test_refuse_mp2

   subroutine full_level_is_whole(z, sym, xyz, owner, level, n_cuts, error)
      !! Level `level` on a partition with that many fragments, cutoff huge
      integer, intent(in) :: z(:), owner(:), level, n_cuts
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: xyz(:, :)
      type(error_type), allocatable, intent(out) :: error

      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      type(error_t) :: err
      integer, allocatable :: charges(:)
      real(dp) :: whole

      call settings(opts)
      opts%rcut = 1.0e6_dp
      opts%level = level
      allocate (charges(maxval(owner)), source=0)
      call run_efmo(z, sym, xyz, owner, charges, opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, res%n_cuts, n_cuts, message="wrong number of detached bonds")
      if (allocated(error)) return

      call whole_rhf(z, sym, xyz, whole, err)
      call check(error,.not. err%has_error(), "the whole-molecule RHF failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      write (*, *) "   level", level, ": EFMO =", res%energy, " RHF =", whole, &
         " diff =", res%energy - whole
      write (*, *) "   induction correction - total =", &
         res%induction_correction - res%polarization_total
      call check(error, res%energy, whole, thr=TOL, &
                 message="EFMO at full level across detached bonds is not the "// &
                 "unfragmented RHF energy")
   end subroutine full_level_is_whole

   subroutine settings(opts)
      type(efmo_options_t), intent(out) :: opts

      opts%basis = BASIS
      opts%bond_breaking = "afo"
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-10_dp
      opts%scf_density_tol = 1.0e-8_dp
      opts%scf_grad_tol = 1.0e-8_dp
      opts%scf%grad_tol = 1.0e-8_dp
   end subroutine settings

   subroutine whole_rhf(z, sym, xyz, energy, err)
      !! The molecule's own RHF energy, Cartesian as EFMO's molecules are
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: xyz(:, :)
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: guess_density(:, :)
      integer :: guess_kind

      energy = 0.0_dp
      call build_czt_molecule(z, sym, xyz, BASIS, mol, err, force_cartesian=.true.)
      if (err%has_error()) return
      call build_restricted_guess(mol, "auto", guess_kind, guess_density, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, sum(z), 200, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err, &
                       guess=guess_kind, guess_density=guess_density, grad_tol=1.0e-8_dp)
      call mol%destroy()
      if (err%has_error()) return
      energy = scf%energy
   end subroutine whole_rhf

   subroutine propane(z, sym, xyz)
      !! Idealised propane, carbons first: C1 end, C2 middle, C3 end
      integer, intent(out) :: z(11)
      character(len=2), intent(out) :: sym(11)
      real(dp), intent(out) :: xyz(3, 11)

      z = [6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1]
      sym = ["C ", "C ", "C ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", "H "]
      xyz = reshape([1.5260_dp, 0.0000_dp, 0.0000_dp, &
                     0.0000_dp, 0.0000_dp, 0.0000_dp, &
                     -0.5716_dp, 1.4149_dp, 0.0000_dp, &
                     2.1553_dp, -0.8900_dp, 0.0000_dp, &
                     2.1553_dp, 0.4450_dp, -0.7707_dp, &
                     2.1553_dp, 0.4450_dp, 0.7707_dp, &
                     -0.3519_dp, -0.5217_dp, 0.8900_dp, &
                     -0.3519_dp, -0.5217_dp, -0.8900_dp, &
                     0.0178_dp, 2.3318_dp, 0.0000_dp, &
                     -1.2200_dp, 1.8317_dp, -0.7707_dp, &
                     -1.2200_dp, 1.8317_dp, 0.7707_dp], [3, 11])*ANG
   end subroutine propane

   subroutine propane_middle_first(z, sym, xyz)
      !! The same propane with the middle carbon numbered first
      integer, intent(out) :: z(11)
      character(len=2), intent(out) :: sym(11)
      real(dp), intent(out) :: xyz(3, 11)

      z = [6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1]
      sym = ["C ", "C ", "C ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", "H "]
      xyz = reshape([0.0000_dp, 0.0000_dp, 0.0000_dp, &     ! C middle
                     1.5260_dp, 0.0000_dp, 0.0000_dp, &     ! C end
                     -0.5716_dp, 1.4149_dp, 0.0000_dp, &    ! C end
                     -0.3519_dp, -0.5217_dp, 0.8900_dp, &   ! H on the middle
                     -0.3519_dp, -0.5217_dp, -0.8900_dp, &
                     2.1553_dp, -0.8900_dp, 0.0000_dp, &    ! H on the first end
                     2.1553_dp, 0.4450_dp, -0.7707_dp, &
                     2.1553_dp, 0.4450_dp, 0.7707_dp, &
                     0.0178_dp, 2.3318_dp, 0.0000_dp, &     ! H on the second end
                     -1.2200_dp, 1.8317_dp, -0.7707_dp, &
                     -1.2200_dp, 1.8317_dp, 0.7707_dp], [3, 11])*ANG
   end subroutine propane_middle_first

   subroutine butane(z, sym, xyz)
      !! Anti butane, C-C 1.53, CCC 112 degrees, C-H 1.09; carbons first
      integer, intent(out) :: z(14)
      character(len=2), intent(out) :: sym(14)
      real(dp), intent(out) :: xyz(3, 14)

      z = [6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
      sym = ["C ", "C ", "C ", "C ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", &
             "H ", "H "]
      xyz = reshape([0.0000_dp, 0.0000_dp, 0.0000_dp, &
                     1.2684_dp, 0.8556_dp, 0.0000_dp, &
                     2.5369_dp, 0.0000_dp, 0.0000_dp, &
                     3.8053_dp, 0.8556_dp, 0.0000_dp, &
                     0.2729_dp, -1.0553_dp, 0.0000_dp, &
                     -0.5889_dp, 0.2224_dp, -0.8898_dp, &
                     -0.5889_dp, 0.2224_dp, 0.8898_dp, &
                     1.2684_dp, 1.4847_dp, 0.8901_dp, &
                     1.2684_dp, 1.4847_dp, -0.8901_dp, &
                     2.5369_dp, -0.6291_dp, -0.8901_dp, &
                     2.5369_dp, -0.6291_dp, 0.8901_dp, &
                     3.5324_dp, 1.9108_dp, 0.0000_dp, &
                     4.3942_dp, 0.6331_dp, -0.8898_dp, &
                     4.3942_dp, 0.6331_dp, 0.8898_dp], [3, 14])*ANG
   end subroutine butane

   subroutine methylcyclopropane(z, sym, xyz)
      !! Methylcyclopropane, ring carbons C1-C3, the methyl carbon C4
      integer, intent(out) :: z(12)
      character(len=2), intent(out) :: sym(12)
      real(dp), intent(out) :: xyz(3, 12)

      z = [6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1]
      sym = ["C ", "C ", "C ", "C ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", "H "]
      xyz = reshape([0.8718_dp, 0.0000_dp, 0.0000_dp, &    ! C1, carries the methyl
                     -0.4359_dp, 0.7550_dp, 0.0000_dp, &   ! C2
                     -0.4359_dp, -0.7550_dp, 0.0000_dp, &  ! C3
                     1.7379_dp, 0.0000_dp, 1.2370_dp, &    ! C4, the methyl carbon
                     1.4970_dp, 0.0000_dp, -0.8929_dp, &   ! H on C1
                     -0.7485_dp, 1.2964_dp, 0.8929_dp, &   ! H on C2
                     -0.7485_dp, 1.2964_dp, -0.8929_dp, &
                     -0.7485_dp, -1.2964_dp, 0.8929_dp, &  ! H on C3
                     -0.7485_dp, -1.2964_dp, -0.8929_dp, &
                     1.9466_dp, 1.0274_dp, 1.5351_dp, &    ! H on C4
                     1.2178_dp, -0.5137_dp, 2.0454_dp, &
                     2.6755_dp, -0.5137_dp, 1.0248_dp], [3, 12])*ANG
   end subroutine methylcyclopropane

end module test_mqc_czt_efmo_covalent

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_efmo_covalent, only: collect_mqc_czt_efmo_covalent_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_efmo_covalent", collect_mqc_czt_efmo_covalent_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
