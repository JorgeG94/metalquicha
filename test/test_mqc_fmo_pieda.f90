!! PIEDA: an FMO2 pair's interaction energy split into Ees, Eex and a residual
module test_mqc_fmo_pieda
   !! Phase A of PIEDA Layer 3 -- FMO only, no cut bonds -- against
   !! `PIEDA_LAYER3_DESIGN.md`'s gates:
   !!
   !! 1. the sum identity, asserted for every decomposed pair;
   !! 2. the total is bit-identical with the analysis on and off;
   !! 3. the cyclic water trimer against an independent PySCF reference
   !!    (`tools/cpu_validation/gen_pieda_refs.py`, ported from the scratch
   !!    probe) and GAMESS's printed kcal/mol;
   !! 4. glycine tripeptide with a water: `test_mqc_fmo_pieda_long`, labelled
   !!    `LONG` for its ten-minute runtime;
   !! 7. a separated pair's Ees is its whole energy, Eex and Ect+mix exactly
   !!    zero;
   !! 9. the pivoted-Cholesky recovery: rank equals the occupied count, and a
   !!    mismatch is an error rather than a warning;
   !! 10. a connected pair comes back flagged and undecomposed, and a pair
   !!     that is not itself connected but still holds a frozen virtual --
   !!     one of its monomers cut elsewhere -- is refused by name, since
   !!     Phase A does not implement PIEDA next to a cut.
   !!
   !! Gate 5/6 (cut systems, `pieda_hl = "projected"`) and gate 8 (rank
   !! independence, which needs MPI and so a test-drive unit test cannot
   !! reach) are outside Phase A; the refusal that stands in for 5 is covered
   !! here.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: to_bohr
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   use mqc_czt_pieda, only: cholesky_occupied_orbitals
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
                  new_unittest("a_pair_next_to_a_cut_is_refused_by_name", &
                               test_cut_pair_refused), &
                  new_unittest("cholesky_recovers_the_occupied_orbitals", &
                               test_cholesky_recovery), &
                  new_unittest("cholesky_refuses_a_rank_mismatch", &
                               test_cholesky_rank_mismatch) &
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
      !! which agrees with GAMESS's PIEDA table (`PIEDA_LAYER3_DESIGN.md`
      !! section 2) to its printed kcal/mol: Eex 5.280/5.291/4.807,
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

   subroutine test_cut_pair_refused(error)
      !! Gate 5/10: a pair next to a cut -- not itself connected, but one of
      !! its monomers still carries a frozen virtual from a cut elsewhere --
      !! is refused by name, Phase A's stand-in for the `pieda_hl` support
      !! that would otherwise be needed to give it a number
      type(error_type), allocatable, intent(out) :: error
      type(fmo_result_t) :: res
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
      opts%pieda = .true.
      opts%scf_energy_tol = 1.0e-10_dp
      opts%scf_density_tol = 1.0e-8_dp

      call run_fmo2(z, sym, xyz, [1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3], &
                    opts, res, err)
      call check(error, err%has_error(), &
                 "PIEDA next to a cut must be refused, not silently wrong")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "not implemented yet") > 0, &
                 "the refusal should say this is not implemented yet, not that it "// &
                 "is impossible")
   end subroutine test_cut_pair_refused

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

   subroutine cyclic_water_trimer(z, sym, xyz)
      !! GAMESS's `3h2o.pieda.inp`, RHF/6-31G* optimised, Angstrom; fragment
      !! `k` is atoms `3k-2..3k`. The geometry `gen_pieda_refs.py` and the
      !! design's own reference table (`PIEDA_LAYER3_DESIGN.md` section 2)
      !! were both taken against.
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
