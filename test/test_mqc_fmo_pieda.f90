!! PIEDA: an FMO2 pair's interaction energy split into Ees, Eex and a residual
module test_mqc_fmo_pieda
   !! The gates PIEDA is held to:
   !!
   !! 1. the sum identity, asserted for every decomposed pair;
   !! 2. the total is bit-identical with the analysis on and off;
   !! 3. the cyclic water trimer against an independent PySCF reference
   !!    (`tools/cpu_validation/gen_pieda_refs.py`, ported from the scratch
   !!    probe) and GAMESS's printed kcal/mol;
   !! 4. glycine tripeptide with a water: `test_mqc_fmo_pieda_long`, labelled
   !!    `LONG` for its ten-minute runtime;
   !! 5. butane cut into two ethyls, water 4.5 A beyond C4, against GAMESS's
   !!    own `IPIEDA=1` table, `pieda_hl = "gamess"`;
   !! 6. `pieda_hl = "projected"` on the same cut system -- checked rather
   !!    than assumed for the frozen-virtual occupation, inside
   !!    `project_out_frozen_virtuals` itself -- and the two modes agreeing
   !!    wherever there is no frozen virtual to project;
   !! 7. a separated pair's Ees is its whole energy, Eex and Ect+mix exactly
   !!    zero;
   !! 9. the pivoted-Cholesky recovery: rank equals the occupied count, and a
   !!    mismatch is an error rather than a warning;
   !! 10. a connected pair comes back flagged and undecomposed;
   !! 11. `pieda_dispersion`: Edi on the water trimer against dftd4's and
   !!     s-dftd3's own programs, HF numbers bit-identical with it on and off,
   !!     nonzero on a separated pair, absent on a connected one. Skipped
   !!     where the build lacks the library, as `test_mqc_dispersion` is.
   !!
   !! Gate 8 (rank independence) needs MPI, which a test-drive unit test
   !! cannot reach; 1 and 2 ranks were compared by hand on the water trimer.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: to_bohr
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   use mqc_czt_pieda, only: cholesky_occupied_orbitals
   use mqc_physical_constants, only: HARTREE_TO_KCALMOL
   use mqc_dispersion_apply, only: dispersion_kind_available
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   implicit none
   private

   public :: collect_mqc_fmo_pieda

   real(dp), parameter :: REF_TOL = 2.0e-9_dp
      !! Against the PySCF reference: measured agreement is a few times
      !! 1e-10 (design gate 3), so this clears it with room for a different
      !! compiler's SCF to land on a slightly different point on the same
      !! tight-tolerance plateau.
   real(dp), parameter :: SUM_TOL = 1.0e-9_dp
      !! The sum identity is exact algebra (`combine_pieda_terms` builds
      !! `ect_mix` as the residual), so the gap is round-off; this only
      !! guards against a genuine bookkeeping slip.

contains

   subroutine collect_mqc_fmo_pieda(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("sum_identity_holds_for_every_decomposed_pair", &
                               test_sum_identity), &
                  new_unittest("total_energy_is_bit_identical_with_pieda_on_and_off", &
                               test_bit_identical), &
                  new_unittest("water_trimer_pieda_matches_an_independent_reference", &
                               test_water_trimer), &
                  new_unittest("a_separated_pairs_ees_is_its_whole_energy", &
                               test_separated_pair), &
                  new_unittest("a_connected_pair_is_reported_undecomposed", &
                               test_connected_pair), &
                  new_unittest("a_pair_next_to_a_cut_matches_gamess_ipieda", &
                               test_cut_pair_matches_gamess), &
                  new_unittest("two_cut_glycine_water_pairs_match_gamess_ipieda", &
                               test_two_cut_glycine_water_matches_gamess), &
                  new_unittest("pairs_across_a_doubly_cut_fragment_are_flagged", &
                               test_doubly_cut_neighbor_flag), &
                  new_unittest("pieda_hl_projected_is_safe_and_variational", &
                               test_projected_mode_is_safe_and_variational), &
                  new_unittest("pieda_hl_modes_coincide_with_no_frozen_virtual", &
                               test_pieda_hl_modes_coincide_with_no_frozen_virtual), &
                  new_unittest("cholesky_recovers_the_occupied_orbitals", &
                               test_cholesky_recovery), &
                  new_unittest("cholesky_refuses_a_rank_mismatch", &
                               test_cholesky_rank_mismatch), &
                  new_unittest("water_trimer_edi_matches_an_independent_reference", &
                               test_water_trimer_edi), &
                  new_unittest("pieda_dispersion_leaves_the_hf_numbers_bit_identical", &
                               test_pieda_dispersion_bit_identical), &
                  new_unittest("a_separated_pairs_edi_is_not_zero", &
                               test_separated_pair_edi), &
                  new_unittest("a_connected_pairs_edi_is_zero", &
                               test_connected_pair_edi) &
                  ]
   end subroutine collect_mqc_fmo_pieda

   subroutine test_sum_identity(error)
      !! Gate 1: Ees + Eex + Ect+mix = the pair's own energy, every pair
      type(error_type), allocatable, intent(out) :: error

      type(fmo_result_t) :: res
      integer :: p

      call water_trimer_pieda_run(res, error)
      if (allocated(error)) return

      do p = 1, size(res%pairs)
         call check(error, res%pairs(p)%pieda, "every unconnected water-trimer pair "// &
                    "should be decomposed")
         if (allocated(error)) return
         call check(error, abs(res%pairs(p)%energy - (res%pairs(p)%ees + res%pairs(p)%eex + &
                                                      res%pairs(p)%ect_mix)) < SUM_TOL, &
                    "Ees + Eex + Ect+mix does not close the pair's own energy")
         if (allocated(error)) return
      end do
   end subroutine test_sum_identity

   subroutine test_bit_identical(error)
      !! Gate 2: the analysis is read-only on the SCF path
      !!
      !! Forced to one OpenMP thread: `build_fock_direct`'s quartet loop
      !! merges over an unordered critical section, so two otherwise-identical
      !! runs already scatter by ~1e-10 under threading (see AGENTS.md,
      !! "OpenMP merge order is not a race") and bit-identity is only a
      !! meaningful question at one thread.
      type(error_type), allocatable, intent(out) :: error
      type(fmo_result_t) :: res_on, res_off
      type(error_t) :: err
      type(fmo_options_t) :: opts
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      integer :: saved

      saved = 1
!$    saved = omp_get_max_threads()
!$    call omp_set_num_threads(1)

      call cyclic_water_trimer(z, sym, xyz)
      call pieda_options(opts, pieda=.false.)
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res_off, err)
      call check(error,.not. err%has_error(), "the pieda-off run failed: "// &
                 err%get_message())
      if (allocated(error)) then
!$       call omp_set_num_threads(saved)
         return
      end if

      call pieda_options(opts, pieda=.true.)
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res_on, err)
      call check(error,.not. err%has_error(), "the pieda-on run failed: "// &
                 err%get_message())
!$    call omp_set_num_threads(saved)
      if (allocated(error)) return

      call check(error, res_on%energy == res_off%energy, &
                 "PIEDA must not change the total by so much as an ulp")
      if (allocated(error)) then
         write (*, *) "   pieda on  =", res_on%energy
         write (*, *) "   pieda off =", res_off%energy
      end if
   end subroutine test_bit_identical

   subroutine test_water_trimer(error)
      !! Gate 3: the cyclic water trimer, HF/6-31G, exact ESP
      !!
      !! Reference from `tools/cpu_validation/gen_pieda_refs.py`:
      !!
      !!   python3 tools/cpu_validation/gen_pieda_refs.py \
      !!     water3_cyclic.xyz 6-31g w3 '[[0,1,2],[3,4,5],[6,7,8]]'
      !!
      !! which agrees with GAMESS's PIEDA table to its printed kcal/mol: Eex 5.280/5.291/4.807,
      !! Ect+mix -2.624/-2.549/-2.381 for pairs 1-2/1-3/2-3.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: PAIR_I(3) = [1, 1, 2], PAIR_J(3) = [2, 3, 3]
      real(dp), parameter :: REF_EES(3) = &
                             [-0.018035578699_dp, -0.018240724203_dp, -0.016095874680_dp]
      real(dp), parameter :: REF_EEX(3) = &
                             [0.008414961460_dp, 0.008432488467_dp, 0.007661073300_dp]
      real(dp), parameter :: REF_ECT(3) = &
                             [-0.004181798363_dp, -0.004061890874_dp, -0.003794447605_dp]
      type(fmo_result_t) :: res
      integer :: p, k

      call water_trimer_pieda_run(res, error)
      if (allocated(error)) return

      call check(error, size(res%pairs), 3, "the water trimer should have three pairs")
      if (allocated(error)) return

      do p = 1, 3
         k = findloc(res%pairs%i == PAIR_I(p) .and. res%pairs%j == PAIR_J(p), .true., dim=1)
         call check(error, k > 0, "a pair is missing")
         if (allocated(error)) return
         call check(error, abs(res%pairs(k)%ees - REF_EES(p)) < REF_TOL, &
                    "Ees does not match the reference")
         if (allocated(error)) then
            write (*, *) "   pair", PAIR_I(p), PAIR_J(p), "Ees - ref =", &
               res%pairs(k)%ees - REF_EES(p)
            return
         end if
         call check(error, abs(res%pairs(k)%eex - REF_EEX(p)) < REF_TOL, &
                    "Eex does not match the reference")
         if (allocated(error)) then
            write (*, *) "   pair", PAIR_I(p), PAIR_J(p), "Eex - ref =", &
               res%pairs(k)%eex - REF_EEX(p)
            return
         end if
         call check(error, abs(res%pairs(k)%ect_mix - REF_ECT(p)) < REF_TOL, &
                    "Ect+mix does not match the reference")
         if (allocated(error)) then
            write (*, *) "   pair", PAIR_I(p), PAIR_J(p), "Ect+mix - ref =", &
               res%pairs(k)%ect_mix - REF_ECT(p)
            return
         end if
      end do
   end subroutine test_water_trimer

   subroutine test_separated_pair(error)
      !! Gate 7: beyond `resdim`, Ees is the pair's whole energy and Eex,
      !! Ect+mix are exactly zero -- no pair SCF, no union state to build
      type(error_type), allocatable, intent(out) :: error
      type(fmo_result_t) :: res
      type(error_t) :: err
      type(fmo_options_t) :: opts
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: ang(3, 6)

      z = [8, 1, 1, 8, 1, 1]
      sym = ["O ", "H ", "H ", "O ", "H ", "H "]
      ! Two waters 15 A apart -- far past `resdim`'s default vdW-scaled 2.0 --
      ! so the only pair is separated.
      ang = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                     0.0_dp, -0.7572_dp, 0.5865_dp, &
                     0.0_dp, 0.7572_dp, 0.5865_dp, &
                     0.0_dp, 0.0_dp, 15.0_dp, &
                     0.0_dp, -0.7572_dp, 15.5865_dp, &
                     0.0_dp, 0.7572_dp, 15.5865_dp], [3, 6])

      opts%basis = "sto-3g"
      opts%pieda = .true.
      opts%resdim = 2.0_dp
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      call run_fmo2(z, sym, to_bohr(ang), [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "two distant waters failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, size(res%pairs), 1, "there should be exactly one pair")
      if (allocated(error)) return
      call check(error, res%pairs(1)%separated, "the pair should be separated")
      if (allocated(error)) return
      call check(error, res%pairs(1)%pieda, "a separated pair is still decomposed")
      if (allocated(error)) return

      ! Exact, not to round-off: a separated pair's term is taken whole as
      ! its Ees, and nothing is subtracted from it.
      call check(error, res%pairs(1)%ees == res%pairs(1)%energy, &
                 "a separated pair's Ees must equal its energy exactly")
      if (allocated(error)) return
      call check(error, res%pairs(1)%eex == 0.0_dp, &
                 "a separated pair's Eex must be exactly zero")
      if (allocated(error)) return
      call check(error, res%pairs(1)%ect_mix == 0.0_dp, &
                 "a separated pair's Ect+mix must be exactly zero")
   end subroutine test_separated_pair

   subroutine test_connected_pair(error)
      !! Gate 10: a pair joined by a detached bond is reported undecomposed,
      !! with no refusal -- `assemble_group` makes the bond whole again
      !! inside the dimer that holds both ends, so this pair holds no frozen
      !! virtual of its own
      type(error_type), allocatable, intent(out) :: error
      type(fmo_result_t) :: res
      type(error_t) :: err
      type(fmo_options_t) :: opts
      integer :: z(8)
      character(len=2) :: sym(8)
      real(dp) :: xyz(3, 8)

      call ethane(z, sym, xyz)
      opts%basis = "sto-3g"
      opts%esp = "exact"
      opts%expansion = "fmo"
      opts%bond_breaking = "afo"
      opts%pieda = .true.
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp

      call run_fmo2(z, sym, xyz, [1, 1, 1, 1, 2, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "ethane across a detached bond failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, size(res%pairs), 1, "ethane has exactly one pair")
      if (allocated(error)) return
      call check(error, res%pairs(1)%connected, "the pair should be connected")
      if (allocated(error)) return
      call check(error,.not. res%pairs(1)%pieda, &
                 "a connected pair must be reported undecomposed")
      if (allocated(error)) return
      call check(error, res%pairs(1)%ees == 0.0_dp .and. res%pairs(1)%eex == 0.0_dp .and. &
                 res%pairs(1)%ect_mix == 0.0_dp, &
                 "an undecomposed pair's term columns must stay at their default")
   end subroutine test_connected_pair

   subroutine test_cut_pair_matches_gamess(error)
      !! Gate 5: a pair next to a cut -- not itself connected, but one of its
      !! monomers still carries a frozen virtual from a cut elsewhere --
      !! against GAMESS's own PIEDA table, `pieda_hl = "gamess"` (the default)
      !!
      !! GAMESS 2026, `erloc_xbut2w_er.inp`/`.log`
      !! (`~/dev/mqc_worktrees/er_localize/er_gamess/`): butane cut into two
      !! ethyls, water 4.5 A beyond C4, RHF/STO-3G, `LOCAL=RUEDNBRG`,
      !! `$FMO RESPPC=2.0 RESDIM=0`, `$FMOPRP IPIEDA=1`. Printed table,
      !! kcal/mol (`I J`, `Ees`, `Eex`, `Ect+mix`):
      !!
      !!   3 1  -0.158  -0.000  -0.000
      !!   3 2   0.303   0.000  -0.001
      !!
      !! The connected pair 1-2 is not decomposed here; its own total,
      !! -9122.246 kcal/mol, is what `test_gamess_butane_water_exact_er` in
      !! `test_mqc_afo_fmo.f90` already checks `energy` against.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: PAIR_I(2) = [1, 2], PAIR_J(2) = [3, 3]
      real(dp), parameter :: GAMESS_EES(2) = [-0.158_dp, 0.303_dp]
      real(dp), parameter :: GAMESS_EEX(2) = [-0.000_dp, 0.000_dp]
      real(dp), parameter :: GAMESS_ECT(2) = [-0.000_dp, -0.001_dp]
      real(dp), parameter :: GAMESS_TOL = 2.0e-3_dp
         !! GAMESS prints three decimals in kcal/mol; this clears the
         !! rounding with room for its own ~1e-6 Ha model-SCF floor (design
         !! risk 6) besides.
      type(fmo_result_t) :: res
      integer :: p, k

      call cut_pair_run("gamess", res, error)
      if (allocated(error)) return

      do p = 1, 2
         k = findloc(res%pairs%i == PAIR_I(p) .and. res%pairs%j == PAIR_J(p), .true., dim=1)
         call check(error, k > 0, "a pair is missing")
         if (allocated(error)) return
         call check(error, res%pairs(k)%pieda, "a pair next to a cut should be decomposed")
         if (allocated(error)) return
         call check(error, &
                    abs(res%pairs(k)%ees*HARTREE_TO_KCALMOL - GAMESS_EES(p)) < GAMESS_TOL, &
                    "Ees does not match GAMESS's PIEDA")
         if (allocated(error)) then
            write (*, *) "   pair", PAIR_I(p), PAIR_J(p), "Ees =", &
               res%pairs(k)%ees*HARTREE_TO_KCALMOL
            return
         end if
         call check(error, &
                    abs(res%pairs(k)%eex*HARTREE_TO_KCALMOL - GAMESS_EEX(p)) < GAMESS_TOL, &
                    "Eex does not match GAMESS's PIEDA")
         if (allocated(error)) then
            write (*, *) "   pair", PAIR_I(p), PAIR_J(p), "Eex =", &
               res%pairs(k)%eex*HARTREE_TO_KCALMOL
            return
         end if
         call check(error, &
                    abs(res%pairs(k)%ect_mix*HARTREE_TO_KCALMOL - GAMESS_ECT(p)) < GAMESS_TOL, &
                    "Ect+mix does not match GAMESS's PIEDA")
         if (allocated(error)) then
            write (*, *) "   pair", PAIR_I(p), PAIR_J(p), "Ect+mix =", &
               res%pairs(k)%ect_mix*HARTREE_TO_KCALMOL
            return
         end if
      end do
   end subroutine test_cut_pair_matches_gamess

   subroutine test_two_cut_glycine_water_matches_gamess(error)
      !! Gate 5 where it bites: glycine tripeptide and water, two C-alpha--C
      !! cuts, so fragments 1, 2 and 3 all carry frozen orbitals and pair 3-4
      !! has a large exchange term
      !!
      !! GAMESS 2026, `erloc_xg3w_er.inp`/`.log`
      !! (`~/dev/mqc_worktrees/er_localize/er_gamess/`), RHF/STO-3G,
      !! `LOCAL=RUEDNBRG`, `$FMO NBODY=2 RAFO(1)=1,1,1 RESPPC=2.0 RESDIM=0
      !! RESPAP=0`, `$FMOPRP IPIEDA=1`. Its unconnected PIEDA rows, kcal/mol
      !! (`I J`, `Ees`, `Eex`, `Ect+mix`):
      !!
      !!   3 1   1.872  -0.000   0.027
      !!   4 1   0.075  -0.000   0.000
      !!   4 2   1.265   0.012  -0.012
      !!   4 3  -6.795   6.048  -4.728
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: PAIR_I(4) = [1, 1, 2, 3], PAIR_J(4) = [3, 4, 4, 4]
      real(dp), parameter :: GAMESS_EES(4) = [1.872_dp, 0.075_dp, 1.265_dp, -6.795_dp]
      real(dp), parameter :: GAMESS_EEX(4) = [-0.000_dp, -0.000_dp, 0.012_dp, 6.048_dp]
      real(dp), parameter :: GAMESS_ECT(4) = [0.027_dp, 0.000_dp, -0.012_dp, -4.728_dp]
      real(dp), parameter :: GAMESS_TOL = 2.0e-3_dp
         !! Three printed decimals, plus GAMESS's ~1e-6 Hartree model-SCF floor
      type(fmo_result_t) :: res
      type(error_t) :: err
      type(fmo_options_t) :: opts
      integer :: z(27)
      character(len=2) :: sym(27)
      real(dp) :: xyz(3, 27)
      real(dp) :: ours(3), ref(3)
      integer :: p, k, c

      call gly3_water_cut(z, sym, xyz)
      opts%basis = "sto-3g"
      opts%esp = "exact"
      opts%expansion = "fmo"
      opts%bond_breaking = "afo"
      opts%afo_localization = "er"
      opts%resppc = 2.0_dp
      opts%resdim = 0.0_dp
      opts%scf_energy_tol = 1.0e-10_dp
      opts%scf_density_tol = 1.0e-8_dp
      opts%outer_tol = 1.0e-9_dp
      opts%pieda = .true.
      call run_fmo2(z, sym, xyz, [1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 3, 3, 2, 2, 2, &
                                  3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4], opts, res, err)
      call check(error,.not. err%has_error(), "two-cut glycine and water failed: "// &
                 err%get_message())
      if (allocated(error)) return

      do p = 1, 4
         k = findloc(res%pairs%i == PAIR_I(p) .and. res%pairs%j == PAIR_J(p), .true., dim=1)
         call check(error, k > 0, "a pair is missing")
         if (allocated(error)) return
         call check(error, res%pairs(k)%pieda, "an unconnected pair should be decomposed")
         if (allocated(error)) return
         ours = [res%pairs(k)%ees, res%pairs(k)%eex, res%pairs(k)%ect_mix]*HARTREE_TO_KCALMOL
         ref = [GAMESS_EES(p), GAMESS_EEX(p), GAMESS_ECT(p)]
         write (*, "(a,2i3,3f12.4)") "   pair, Ees Eex Ect+mix (kcal/mol):", &
            PAIR_I(p), PAIR_J(p), ours
         do c = 1, 3
            call check(error, abs(ours(c) - ref(c)) < GAMESS_TOL, &
                       "a PIEDA term does not match GAMESS's IPIEDA=1 table")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_two_cut_glycine_water_matches_gamess

   subroutine test_doubly_cut_neighbor_flag(error)
      !! Butane as CH3 | CH2CH2 | CH3, with the middle fragment holding both
      !! detached atoms (C2 and C3, bonded to each other): the pair of methyls
      !! either side of it is flagged, and neither pair touching it is
      type(error_type), allocatable, intent(out) :: error
      type(fmo_result_t) :: res
      type(error_t) :: err
      type(fmo_options_t) :: opts
      integer :: z(17)
      character(len=2) :: sym(17)
      real(dp) :: xyz(3, 17)
      integer :: k

      call butane_water(z, sym, xyz)
      opts%basis = "sto-3g"
      opts%esp = "exact"
      opts%expansion = "fmo"
      opts%bond_breaking = "afo"
      opts%detached = [2, 3]
      opts%resdim = 0.0_dp
      opts%pieda = .true.
      opts%scf_energy_tol = 1.0e-10_dp
      opts%scf_density_tol = 1.0e-8_dp
      call run_fmo2(z(1:14), sym(1:14), xyz(:, 1:14), &
                    [1, 2, 2, 3, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "butane in three failed: "//err%get_message())
      if (allocated(error)) return

      do k = 1, size(res%pairs)
         write (*, "(a,2i3,2l3)") "   pair, connected, doubly_cut_neighbor:", &
            res%pairs(k)%i, res%pairs(k)%j, res%pairs(k)%connected, &
            res%pairs(k)%doubly_cut_neighbor
         call check(error, res%pairs(k)%doubly_cut_neighbor .eqv. &
                    (res%pairs(k)%i == 1 .and. res%pairs(k)%j == 3), &
                    "only the methyl-methyl pair should be flagged")
         if (allocated(error)) return
      end do
   end subroutine test_doubly_cut_neighbor_flag

   subroutine test_projected_mode_is_safe_and_variational(error)
      !! Gate 6: `pieda_hl = "projected"` on the same cut system
      !!
      !! `project_out_frozen_virtuals` checks its own residual occupation
      !! internally (design gate 6's "0 to 1e-12", measured at 2-4e-17 on
      !! this system and on Gly3+water's two cuts) and would refuse rather
      !! than silently carry on if it ever failed, so a clean run already
      !! covers that half. What is asserted here is the other half: on a
      !! pair whose union state variationally cannot do better than the
      !! constrained pair SCF, `Ect+mix` is not meaningfully positive -- the
      !! measured values are -0.0002 and -0.0002 Hartree here, comfortably
      !! inside the tolerance.
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: NOT_POSITIVE_TOL = 1.0e-6_dp
      type(fmo_result_t) :: res
      integer :: k

      call cut_pair_run("projected", res, error)
      if (allocated(error)) return

      do k = 1, size(res%pairs)
         if (.not. res%pairs(k)%pieda) cycle
         call check(error, res%pairs(k)%ect_mix < NOT_POSITIVE_TOL, &
                    "a projected pair's Ect+mix should not be meaningfully positive")
         if (allocated(error)) then
            write (*, *) "   pair", res%pairs(k)%i, res%pairs(k)%j, "Ect+mix =", &
               res%pairs(k)%ect_mix
            return
         end if
      end do
   end subroutine test_projected_mode_is_safe_and_variational

   subroutine test_pieda_hl_modes_coincide_with_no_frozen_virtual(error)
      !! Gate 6: the two `pieda_hl` modes are identical wherever a pair holds
      !! no frozen virtual -- the water trimer has no cut at all, so every
      !! pair is such a pair, and `"projected"`'s extra step never runs
      !!
      !! Forced to one OpenMP thread, as `test_bit_identical` is: the two
      !! modes run the identical code path here, so any difference between
      !! them would only be the Fock build's own unordered-critical-section
      !! scatter under threading (AGENTS.md, "OpenMP merge order is not a
      !! race"), which is not what this gate is asking about.
      type(error_type), allocatable, intent(out) :: error
      type(fmo_result_t) :: gamess_res, projected_res
      type(error_t) :: err
      type(fmo_options_t) :: opts
      integer :: z(9), p
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      integer :: saved

      saved = 1
!$    saved = omp_get_max_threads()
!$    call omp_set_num_threads(1)

      call cyclic_water_trimer(z, sym, xyz)

      call pieda_options(opts, pieda=.true.)
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, gamess_res, err)
      call check(error,.not. err%has_error(), "the gamess-mode run failed: "// &
                 err%get_message())
      if (allocated(error)) then
!$       call omp_set_num_threads(saved)
         return
      end if

      opts%pieda_hl = "projected"
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, projected_res, err)
      call check(error,.not. err%has_error(), "the projected-mode run failed: "// &
                 err%get_message())
!$    call omp_set_num_threads(saved)
      if (allocated(error)) return

      call check(error, projected_res%energy == gamess_res%energy, &
                 "the two modes must give the same total with no frozen virtual anywhere")
      if (allocated(error)) return
      do p = 1, size(gamess_res%pairs)
         call check(error, projected_res%pairs(p)%eex == gamess_res%pairs(p)%eex .and. &
                    projected_res%pairs(p)%ect_mix == gamess_res%pairs(p)%ect_mix, &
                    "the two modes must agree pair by pair with no frozen virtual anywhere")
         if (allocated(error)) return
      end do
   end subroutine test_pieda_hl_modes_coincide_with_no_frozen_virtual

   subroutine cut_pair_run(pieda_hl, res, error)
      !! Butane cut into two ethyls, water 4.5 A beyond C4 -- the system
      !! `cut_pair_matches_gamess` and `projected_mode_is_safe_and_variational`
      !! both run, at the `pieda_hl` named
      character(len=*), intent(in) :: pieda_hl
      type(fmo_result_t), intent(out) :: res
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(fmo_options_t) :: opts
      integer :: z(17)
      character(len=2) :: sym(17)
      real(dp) :: xyz(3, 17)

      call butane_water(z, sym, xyz)
      opts%basis = "sto-3g"
      opts%esp = "exact"
      opts%expansion = "fmo"
      opts%bond_breaking = "afo"
      opts%afo_localization = "er"
      opts%resppc = 2.0_dp
      opts%resdim = 0.0_dp
      opts%pieda = .true.
      opts%pieda_hl = pieda_hl
      opts%scf_energy_tol = 1.0e-10_dp
      opts%scf_density_tol = 1.0e-8_dp
      opts%outer_tol = 1.0e-9_dp

      call run_fmo2(z, sym, xyz, [1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3], &
                    opts, res, err)
      call check(error,.not. err%has_error(), "butane cut into two ethyls, with water, "// &
                 "failed: "//err%get_message())
   end subroutine cut_pair_run

   subroutine test_cholesky_recovery(error)
      !! Gate 9: `C_occ` recovered from `D = 2 C C^T` reproduces `D`, and the
      !! factorisation's rank is the occupied count
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: c_true(4, 2), density(4, 4)
      real(dp), allocatable :: c_occ(:, :)
      real(dp) :: rebuilt(4, 4)

      ! Two orthonormal columns of a 4-dimensional space -- an arbitrary but
      ! idempotent closed-shell density, `S = I` so orthonormal is enough.
      c_true = reshape([0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, &
                        0.5_dp, 0.5_dp, -0.5_dp, -0.5_dp], [4, 2])
      density = 2.0_dp*matmul(c_true, transpose(c_true))

      call cholesky_occupied_orbitals(density, 2, "test density", c_occ, err)
      call check(error,.not. err%has_error(), "recovery failed: "//err%get_message())
      if (allocated(error)) return
      call check(error, size(c_occ, 2), 2, "the recovered rank should be the occupied count")
      if (allocated(error)) return

      rebuilt = 2.0_dp*matmul(c_occ, transpose(c_occ))
      call check(error, maxval(abs(rebuilt - density)) < 1.0e-12_dp, &
                 "2 C_occ C_occ^T should reproduce the original density")
   end subroutine test_cholesky_recovery

   subroutine test_cholesky_rank_mismatch(error)
      !! Gate 9: a rank that disagrees with the declared occupied count is an
      !! error, not a warning
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: c_true(4, 2), density(4, 4)
      real(dp), allocatable :: c_occ(:, :)

      c_true = reshape([0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, &
                        0.5_dp, 0.5_dp, -0.5_dp, -0.5_dp], [4, 2])
      density = 2.0_dp*matmul(c_true, transpose(c_true))

      call cholesky_occupied_orbitals(density, 1, "test density", c_occ, err)
      call check(error, err%has_error(), "a rank-2 density declared as one occupied "// &
                 "orbital must be refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "test density") > 0, &
                 "the refusal should name the density it was given")
   end subroutine test_cholesky_rank_mismatch

   subroutine pieda_options(opts, pieda)
      !! The tight-tolerance settings gate 3/4 need, exact field, every pair
      !! solved -- matching `gen_pieda_refs.py`'s own treatment of the field
      logical, intent(in) :: pieda
      type(fmo_options_t), intent(out) :: opts

      opts%basis = "6-31g"
      opts%esp = "exact"
      opts%expansion = "fmo"
      opts%resppc = -1.0_dp
      opts%resdim = 0.0_dp
      opts%outer_tol = 1.0e-12_dp
      opts%max_outer = 100
      opts%scf_energy_tol = 1.0e-12_dp
      opts%scf_density_tol = 1.0e-10_dp
      opts%pieda = pieda
   end subroutine pieda_options

   subroutine edi_options(opts, kind)
      !! PIEDA on, STO-3G, exact field, every pair solved, default SCF
      !! tolerances, and `pieda_dispersion = kind` -- the cheap settings the
      !! Edi tests use, where the HF part only has to finish
      character(len=*), intent(in) :: kind
      type(fmo_options_t), intent(out) :: opts

      opts%basis = "sto-3g"
      opts%esp = "exact"
      opts%expansion = "fmo"
      opts%resdim = 0.0_dp
      opts%pieda = .true.
      opts%pieda_dispersion = kind
   end subroutine edi_options

   subroutine water_trimer_pieda_run(res, error)
      !! The cyclic water trimer, PIEDA on, at gate 3's tolerances
      type(fmo_result_t), intent(out) :: res
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(fmo_options_t) :: opts
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)

      call cyclic_water_trimer(z, sym, xyz)
      call pieda_options(opts, pieda=.true.)
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "the water trimer PIEDA run failed: "// &
                 err%get_message())
   end subroutine water_trimer_pieda_run

   subroutine test_water_trimer_edi(error)
      !! Edi on the cyclic water trimer, D4 and D3(BJ), against the two
      !! libraries' own command-line programs
      !!
      !! Independent reference: xtb's standalone builds, dftd4 3.7.0 and
      !! s-dftd3 1.2.1 (not the libraries this program links), run on the
      !! three monomers and three dimers of `water3_cyclic.xyz` split as
      !! `cyclic_water_trimer` splits it, each written as its own .xyz in
      !! Angstrom:
      !!
      !!   dftd4 -c 0 -f hf  monoK.xyz / dimerIJ.xyz
      !!   s-dftd3 --bj hf   monoK.xyz / dimerIJ.xyz
      !!
      !! and `Edi = E(dimerIJ) - E(monoI) - E(monoJ)` from the printed
      !! "Dispersion energy" lines. The pinned dftd4 4.2.0 and s-dftd3 1.4.0
      !! programs print the same energies to all fourteen digits.
      !!
      !! STO-3G and loose SCF tolerances: Edi depends on the geometry and the
      !! charges alone, so the HF part only has to finish.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: PAIR_I(3) = [1, 1, 2], PAIR_J(3) = [2, 3, 3]
      real(dp), parameter :: REF_D4(3) = &
                             [-1.9162492975e-3_dp, -1.9169965590e-3_dp, -1.8870430470e-3_dp]
      real(dp), parameter :: REF_D3BJ(3) = &
                             [-3.0888832909e-3_dp, -3.0630954919e-3_dp, -2.9867146243e-3_dp]
      real(dp), parameter :: EDI_TOL = 1.0e-12_dp
         !! The reference is quoted to 1e-13 Hartree; the programs read
         !! Angstrom and convert with their own constant, which moves Edi by
         !! well under this. A wrong atom list or charge moves it by 1e-5 and
         !! more.
      character(len=4), parameter :: KINDS(2) = ["d4  ", "d3bj"]
      real(dp) :: ref(3)
      type(fmo_result_t) :: res
      type(error_t) :: err
      type(fmo_options_t) :: opts
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      integer :: p, k, m

      call cyclic_water_trimer(z, sym, xyz)
      do m = 1, size(KINDS)
         if (.not. dispersion_kind_available(trim(KINDS(m)))) cycle
         if (m == 1) then
            ref = REF_D4
         else
            ref = REF_D3BJ
         end if

         call edi_options(opts, trim(KINDS(m)))
         call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
         call check(error,.not. err%has_error(), "the water trimer Edi run ("// &
                    trim(KINDS(m))//") failed: "//err%get_message())
         if (allocated(error)) return
         call check(error, size(res%pairs), 3, "the water trimer should have three pairs")
         if (allocated(error)) return

         do p = 1, 3
            k = findloc(res%pairs%i == PAIR_I(p) .and. res%pairs%j == PAIR_J(p), .true., &
                        dim=1)
            call check(error, k > 0, "a pair is missing")
            if (allocated(error)) return
            call check(error, res%pairs(k)%pieda, "every unconnected pair should be decomposed")
            if (allocated(error)) return
            call check(error, abs(res%pairs(k)%edi - ref(p)) < EDI_TOL, &
                       "Edi ("//trim(KINDS(m))//") does not match the library's own program")
            if (allocated(error)) then
               write (*, *) "   ", trim(KINDS(m)), " pair", PAIR_I(p), PAIR_J(p), &
                  "Edi =", res%pairs(k)%edi, " Edi - ref =", res%pairs(k)%edi - ref(p)
               return
            end if
         end do
      end do
   end subroutine test_water_trimer_edi

   subroutine test_pieda_dispersion_bit_identical(error)
      !! `pieda_dispersion` changes no HF number by so much as an ulp:
      !! the total, every pair's own `energy`, and Ees/Eex/Ect+mix are
      !! bit-identical with it on and off. Edi is read beside them, never
      !! folded into any of them.
      !!
      !! Forced to one OpenMP thread, as `test_bit_identical` is. STO-3G,
      !! since bit-identity does not care which basis.
      type(error_type), allocatable, intent(out) :: error
      type(fmo_result_t) :: res_on, res_off
      type(error_t) :: err
      type(fmo_options_t) :: opts
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      integer :: p, saved

      if (.not. dispersion_kind_available("d4")) return

      saved = 1
!$    saved = omp_get_max_threads()
!$    call omp_set_num_threads(1)

      call cyclic_water_trimer(z, sym, xyz)
      call edi_options(opts, "none")
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res_off, err)
      call check(error,.not. err%has_error(), "the dispersion-off run failed: "// &
                 err%get_message())
      if (allocated(error)) then
!$       call omp_set_num_threads(saved)
         return
      end if

      call edi_options(opts, "d4")
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res_on, err)
!$    call omp_set_num_threads(saved)
      call check(error,.not. err%has_error(), "the dispersion-on run failed: "// &
                 err%get_message())
      if (allocated(error)) return

      call check(error, res_on%energy == res_off%energy, &
                 "pieda_dispersion must not change the total by so much as an ulp")
      if (allocated(error)) return
      call check(error, size(res_on%pairs), size(res_off%pairs), "the pair lists differ")
      if (allocated(error)) return
      do p = 1, size(res_off%pairs)
         call check(error, res_on%pairs(p)%energy == res_off%pairs(p)%energy .and. &
                    res_on%pairs(p)%ees == res_off%pairs(p)%ees .and. &
                    res_on%pairs(p)%eex == res_off%pairs(p)%eex .and. &
                    res_on%pairs(p)%ect_mix == res_off%pairs(p)%ect_mix, &
                    "pieda_dispersion must not change a pair's HF terms")
         if (allocated(error)) return
         call check(error, res_off%pairs(p)%edi == 0.0_dp, &
                    "Edi must stay at its default with pieda_dispersion off")
         if (allocated(error)) return
         call check(error, res_on%pairs(p)%edi /= 0.0_dp, &
                    "Edi should be nonzero on the water trimer with pieda_dispersion on")
         if (allocated(error)) return
      end do
   end subroutine test_pieda_dispersion_bit_identical

   subroutine test_separated_pair_edi(error)
      !! A separated pair still gets a nonzero Edi: dispersion does not
      !! vanish at the range `resdim` cuts the pair SCF off at
      !!
      !! Only nonzero is asserted, not negative. With the three-body term
      !! "-D4" includes, dftd4's own program gives this pair a tiny positive
      !! Edi: `dftd4 -c 0 -f hf` on the two waters and their union prints
      !! -9.8088920013418E-04 twice and -1.9616201546051E-03, so
      !! Edi = +1.58e-7 Hartree. s-dftd3 (`--bj hf`, two-body only) gives
      !! -8.7e-8.
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: REF_D4 = 1.5824566326e-7_dp
      type(fmo_result_t) :: res
      type(error_t) :: err
      type(fmo_options_t) :: opts
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: ang(3, 6)

      if (.not. dispersion_kind_available("d4")) return

      z = [8, 1, 1, 8, 1, 1]
      sym = ["O ", "H ", "H ", "O ", "H ", "H "]
      ang = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                     0.0_dp, -0.7572_dp, 0.5865_dp, &
                     0.0_dp, 0.7572_dp, 0.5865_dp, &
                     0.0_dp, 0.0_dp, 15.0_dp, &
                     0.0_dp, -0.7572_dp, 15.5865_dp, &
                     0.0_dp, 0.7572_dp, 15.5865_dp], [3, 6])

      opts%basis = "sto-3g"
      opts%pieda = .true.
      opts%pieda_dispersion = "d4"
      opts%resdim = 2.0_dp
      call run_fmo2(z, sym, to_bohr(ang), [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "two distant waters failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, size(res%pairs), 1, "there should be exactly one pair")
      if (allocated(error)) return
      call check(error, res%pairs(1)%separated, "the pair should be separated")
      if (allocated(error)) return
      call check(error, res%pairs(1)%pieda, "a separated pair is still decomposed")
      if (allocated(error)) return
      call check(error, res%pairs(1)%edi /= 0.0_dp, &
                 "a separated pair's Edi should be nonzero")
      if (allocated(error)) return
      call check(error, abs(res%pairs(1)%edi - REF_D4) < 1.0e-12_dp, &
                 "a separated pair's Edi does not match dftd4's own program")
      if (allocated(error)) write (*, *) "   Edi =", res%pairs(1)%edi
   end subroutine test_separated_pair_edi

   subroutine test_connected_pair_edi(error)
      !! A connected pair's Edi stays at its default: its term already
      !! carries the bond itself, so nothing here decomposes it further
      !!
      !! Ethane's only pair is that connected one, so `pieda_dispersion_term`
      !! is never called and this run needs no dispersion library at all --
      !! unlike the other Edi tests, this one is not skipped on a build
      !! without dftd4, and asking for `d4` here must not be refused either.
      type(error_type), allocatable, intent(out) :: error
      type(fmo_result_t) :: res
      type(error_t) :: err
      type(fmo_options_t) :: opts
      integer :: z(8)
      character(len=2) :: sym(8)
      real(dp) :: xyz(3, 8)

      call ethane(z, sym, xyz)
      opts%basis = "sto-3g"
      opts%esp = "exact"
      opts%expansion = "fmo"
      opts%bond_breaking = "afo"
      opts%pieda = .true.
      opts%pieda_dispersion = "d4"
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp

      call run_fmo2(z, sym, xyz, [1, 1, 1, 1, 2, 2, 2, 2], opts, res, err)

      call check(error,.not. err%has_error(), "ethane across a detached bond failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, size(res%pairs), 1, "ethane has exactly one pair")
      if (allocated(error)) return
      call check(error, res%pairs(1)%connected, "the pair should be connected")
      if (allocated(error)) return
      call check(error, res%pairs(1)%edi == 0.0_dp, &
                 "a connected pair's Edi must stay at its default")
   end subroutine test_connected_pair_edi

   subroutine cyclic_water_trimer(z, sym, xyz)
      !! GAMESS's `3h2o.pieda.inp`, RHF/6-31G* optimised, Angstrom; fragment
      !! `k` is atoms `3k-2..3k`; `sample_inputs/water3_cyclic.xyz`.
      integer, intent(out) :: z(9)
      character(len=2), intent(out) :: sym(9)
      real(dp), intent(out) :: xyz(3, 9)
      real(dp) :: ang(3, 9)

      z = [8, 1, 1, 8, 1, 1, 8, 1, 1]
      sym = ["O ", "H ", "H ", "O ", "H ", "H ", "O ", "H ", "H "]
      ang = reshape([-0.8309834169_dp, 1.3972845867_dp, -0.2245777899_dp, &
                     -1.7252647239_dp, 1.0653378323_dp, -0.1532773281_dp, &
                     -0.7392089607_dp, 2.0451809087_dp, 0.4613985156_dp, &
                     -0.3320345076_dp, -1.3821619786_dp, 0.2567741911_dp, &
                     -0.1833783620_dp, -0.4480121188_dp, 0.1145770742_dp, &
                     0.1095285385_dp, -1.8294391791_dp, -0.4524662378_dp, &
                     -3.0234372871_dp, -0.3756747342_dp, 0.2555351867_dp, &
                     -2.2955882864_dp, -0.9880884300_dp, 0.3537585219_dp, &
                     -3.5920719939_dp, -0.7444718872_dp, -0.4067791338_dp], [3, 9])
      xyz = to_bohr(ang)
   end subroutine cyclic_water_trimer

   subroutine ethane(z, sym, xyz)
      !! Ethane, split at the C-C bond into two methyls
      integer, intent(out) :: z(8)
      character(len=2), intent(out) :: sym(8)
      real(dp), intent(out) :: xyz(3, 8)
      real(dp) :: ang(3, 8)
      integer :: i

      z = [6, 1, 1, 1, 6, 1, 1, 1]
      do i = 1, 8
         if (z(i) == 6) then
            sym(i) = "C "
         else
            sym(i) = "H "
         end if
      end do
      ang = reshape([0.000_dp, 0.000_dp, 0.768_dp, &
                     -1.019_dp, 0.000_dp, 1.157_dp, &
                     0.510_dp, 0.883_dp, 1.157_dp, &
                     0.510_dp, -0.883_dp, 1.157_dp, &
                     0.000_dp, 0.000_dp, -0.768_dp, &
                     1.019_dp, 0.000_dp, -1.157_dp, &
                     -0.510_dp, -0.883_dp, -1.157_dp, &
                     -0.510_dp, 0.883_dp, -1.157_dp], [3, 8])
      xyz = to_bohr(ang)
   end subroutine ethane

   subroutine butane_water(z, sym, xyz)
      !! Anti butane, cut into two ethyls, and a water 4.5 A beyond C4, as
      !! the GAMESS deck `test_mqc_afo_fmo.f90` matches has them
      integer, intent(out) :: z(17)
      character(len=2), intent(out) :: sym(17)
      real(dp), intent(out) :: xyz(3, 17)
      real(dp) :: ang(3, 17)

      z = [6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1]
      sym = ["C ", "C ", "C ", "C ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", &
             "H ", "H ", "O ", "H ", "H "]
      ang = reshape([0.000000_dp, 0.000000_dp, 0.000000_dp, &
                     1.268400_dp, 0.855600_dp, 0.000000_dp, &
                     2.536900_dp, 0.000000_dp, 0.000000_dp, &
                     3.805300_dp, 0.855600_dp, 0.000000_dp, &
                     0.272900_dp, -1.055300_dp, 0.000000_dp, &
                     -0.588900_dp, 0.222400_dp, -0.889800_dp, &
                     -0.588900_dp, 0.222400_dp, 0.889800_dp, &
                     1.268400_dp, 1.484700_dp, 0.890100_dp, &
                     1.268400_dp, 1.484700_dp, -0.890100_dp, &
                     2.536900_dp, -0.629100_dp, -0.890100_dp, &
                     2.536900_dp, -0.629100_dp, 0.890100_dp, &
                     3.532400_dp, 1.910800_dp, 0.000000_dp, &
                     4.394200_dp, 0.633100_dp, -0.889800_dp, &
                     4.394200_dp, 0.633100_dp, 0.889800_dp, &
                     7.535896_dp, 3.372076_dp, 0.000000_dp, &
                     7.575283_dp, 4.330466_dp, 0.000000_dp, &
                     8.439273_dp, 3.049628_dp, 0.000000_dp], [3, 17])
      xyz = to_bohr(ang)
   end subroutine butane_water

   subroutine gly3_water_cut(z, sym, xyz)
      !! `validation/inputs/sample_inputs/gly3_water_pair.xyz`, as `test_mqc_afo_fmo` builds it
      integer, intent(out) :: z(27)
      character(len=2), intent(out) :: sym(27)
      real(dp), intent(out) :: xyz(3, 27)
      real(dp) :: ang(3, 27)
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
      ang = reshape([ &
                    0.0171625298_dp, -0.4776667709_dp, -0.0077801388_dp, &   ! N
                    1.3251492481_dp, 0.1638239831_dp, 0.0713249069_dp, &   ! C
                    1.8818395599_dp, 0.1764813685_dp, 1.4667973423_dp, &   ! C
                    1.1563644386_dp, 0.4758564459_dp, 2.4030731780_dp, &   ! O
                    2.0041403197_dp, -0.3893217244_dp, -0.6156078332_dp, &   ! H
                    1.2933738676_dp, 1.2140808724_dp, -0.2903017566_dp, &   ! H
                    -0.6557592247_dp, -0.0682256808_dp, 0.6785523482_dp, &   ! H
                    -0.3826962098_dp, -0.2691894812_dp, -0.9506317163_dp, &   ! H
                    3.2093591995_dp, -0.0780774266_dp, 1.6702200732_dp, &   ! N
                    3.8489825798_dp, -0.0589263473_dp, 2.9842578467_dp, &   ! C
                    5.3502343581_dp, -0.0788662970_dp, 2.9476716562_dp, &   ! C
                    5.9543074560_dp, -0.1656759551_dp, 1.8893430618_dp, &   ! O
                    3.5421254604_dp, 0.8561169960_dp, 3.5393994122_dp, &   ! H
                    3.4986665918_dp, -0.9402544817_dp, 3.5643998498_dp, &   ! H
                    3.7845901118_dp, -0.3119789206_dp, 0.8286081985_dp, &   ! H
                    6.0352251963_dp, 0.0003525130_dp, 4.1282386693_dp, &   ! N
                    7.4955375902_dp, -0.0138802141_dp, 4.2014382315_dp, &   ! C
                    8.0730347718_dp, 0.0277800836_dp, 5.5909529457_dp, &   ! C
                    7.3557278976_dp, 0.0641983810_dp, 6.5759347789_dp, &   ! O
                    7.8694940865_dp, -0.9353711779_dp, 3.7021749317_dp, &   ! H
                    7.8868335534_dp, 0.8596348618_dp, 3.6344677391_dp, &   ! H
                    5.4670886620_dp, 0.0786510231_dp, 5.0034540291_dp, &   ! H
                    9.3768940878_dp, 0.0221621974_dp, 5.7818296269_dp, &   ! O
                    9.9376629532_dp, -0.0106298905_dp, 4.9380771002_dp, &   ! H
                    7.3635223930_dp, -0.3681902986_dp, -0.5795840740_dp, &   ! O
                    6.8902239587_dp, -0.3001739023_dp, 0.2496289275_dp, &   ! H
                    6.6876213700_dp, -0.2710584500_dp, -1.2503709636_dp &   ! H
                    ], [3, 27])
      xyz = to_bohr(ang)
   end subroutine gly3_water_cut

end module test_mqc_fmo_pieda

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fmo_pieda, only: collect_mqc_fmo_pieda
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fmo_pieda", collect_mqc_fmo_pieda)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
