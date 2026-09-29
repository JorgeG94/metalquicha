!! MP2, RI-MP2 and SCS-MP2 as correlation on embedded Hartree-Fock, against the supermolecule
module test_mqc_fmo_mp2
   !! The MP2-family counterpart of `mqc_fmo_dft`: with the fragment count
   !! equal to the expansion level, FMO and EE-MBE are exact by inclusion and
   !! exclusion, so a correlated total has to reproduce an ordinary
   !! Hartree-Fock-plus-MP2 calculation on the whole molecule -- the same
   !! identity `test_mqc_fmo_dft` checks for Kohn-Sham, carried over to MP2,
   !! RI-MP2 and SCS-MP2 (see `mqc_docs/source/developer_fragment_solver.rst`,
   !! phase 3: correlation is added per fragment and per n-mer on top of the
   !! embedded Hartree-Fock, as GAMESS FMO-MP2 does).
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: to_bohr
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: run_czt_rhf, rhf_result_t
   use mqc_czt_mp2, only: mp2_result_t, run_czt_mp2, run_czt_ri_mp2
   use mqc_elements, only: core_orbital_count
   use mqc_method_config, only: correlation_config_t
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   implicit none
   private

   public :: collect_mqc_fmo_mp2

   real(dp), parameter :: TOL = 1.0e-8_dp
      !! Both sides are the same Hartree-Fock-plus-MP2 calculation on the same
      !! atoms, in the same basis, so what separates them is rounding in the
      !! fragment assembly rather than any physics.

   character(len=*), parameter :: BASIS = "sto-3g"
   character(len=*), parameter :: AUX = "cc-pvdz-rifit"
      !! RI-MP2's fitting basis. Not matched to `BASIS`, exactly as
      !! `test_mqc_czt_efmo` uses it: the identity checked here is FMO against
      !! the supermolecule at the *same* fitting basis on both sides, not
      !! RI-MP2 against conventional MP2.

contains

   subroutine collect_mqc_fmo_mp2(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fmo2_water_dimer_mp2_is_the_supermolecule", test_dimer_mp2), &
                  new_unittest("fmo2_water_dimer_ri_mp2_is_the_supermolecule", &
                               test_dimer_ri_mp2), &
                  new_unittest("fmo2_water_dimer_scs_mp2_is_the_supermolecule", &
                               test_dimer_scs_mp2), &
                  new_unittest("fmo3_water_trimer_mp2_is_the_supermolecule", test_trimer_mp2), &
                  new_unittest("eembe2_water_dimer_mp2_is_the_supermolecule", &
                               test_eembe_dimer_mp2), &
                  new_unittest("fmo2_water_dimer_mp2_without_a_frozen_core", &
                               test_dimer_mp2_no_frozen_core), &
                  new_unittest("eembe2_water_trimer_mp2_family_matches_pyscf", &
                               test_eembe2_trimer_pyscf), &
                  new_unittest("fmo2_and_eembe2_water_trimer_mp2_against_the_supermolecule", &
                               test_trimer_pair_truncated_mp2), &
                  new_unittest("plain_mbe2_water_dimer_mp2_is_the_supermolecule", &
                               test_plain_mbe_dimer_mp2), &
                  new_unittest("fmo2_water_trimer_cyclic_mp2_matches_gamess", &
                               test_trimer_cyclic_mp2_gamess), &
                  new_unittest("fmo2_water_trimer_cyclic_separated_pair_carries_no_correlation", &
                               test_trimer_cyclic_separated_pair) &
                  ]
   end subroutine collect_mqc_fmo_mp2

   subroutine test_dimer_mp2(error)
      !! Two waters, FMO2 (the default field), MP2, frozen core
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)
      real(dp) :: whole

      call water_dimer(z, sym, xyz)
      call mp2_options(opts, use_ri=.false., use_scs=.false., freeze_core=.true.)
      opts%level = 2

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO2 MP2 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_mp2(z, sym, xyz, use_ri=.false., use_scs=.false., freeze_core=.true., &
                             nelec=20, energy=whole, error=err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call check(error, abs(res%energy - whole) < TOL, &
                 "FMO2 MP2 on the water dimer does not reproduce the supermolecule")
      if (allocated(error)) then
         write (*, *) "   fmo2        =", res%energy
         write (*, *) "   supermolecule=", whole
         write (*, *) "   difference  =", res%energy - whole
      end if
   end subroutine test_dimer_mp2

   subroutine test_dimer_ri_mp2(error)
      !! The same dimer, RI-MP2
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)
      real(dp) :: whole

      call water_dimer(z, sym, xyz)
      call mp2_options(opts, use_ri=.true., use_scs=.false., freeze_core=.true.)
      opts%level = 2

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO2 RI-MP2 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_mp2(z, sym, xyz, use_ri=.true., use_scs=.false., freeze_core=.true., &
                             nelec=20, energy=whole, error=err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      call check(error, abs(res%energy - whole) < TOL, &
                 "FMO2 RI-MP2 on the water dimer does not reproduce the supermolecule")
      if (allocated(error)) then
         write (*, *) "   fmo2        =", res%energy
         write (*, *) "   supermolecule=", whole
         write (*, *) "   difference  =", res%energy - whole
      end if
   end subroutine test_dimer_ri_mp2

   subroutine test_dimer_scs_mp2(error)
      !! The same dimer, SCS-MP2 -- the spin-component scales this build reads
      !! the deck with (`correlation_config_t`'s own defaults)
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)
      real(dp) :: whole

      call water_dimer(z, sym, xyz)
      call mp2_options(opts, use_ri=.false., use_scs=.true., freeze_core=.true.)
      opts%level = 2

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO2 SCS-MP2 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_mp2(z, sym, xyz, use_ri=.false., use_scs=.true., freeze_core=.true., &
                             nelec=20, energy=whole, error=err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      call check(error, abs(res%energy - whole) < TOL, &
                 "FMO2 SCS-MP2 on the water dimer does not reproduce the supermolecule")
      if (allocated(error)) then
         write (*, *) "   fmo2        =", res%energy
         write (*, *) "   supermolecule=", whole
         write (*, *) "   difference  =", res%energy - whole
      end if
   end subroutine test_dimer_scs_mp2

   subroutine test_trimer_mp2(error)
      !! Three waters, FMO3 (level equal to the fragment count), MP2
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      real(dp) :: whole

      call water_trimer(z, sym, xyz)
      call mp2_options(opts, use_ri=.false., use_scs=.false., freeze_core=.true.)
      opts%level = 3

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO3 MP2 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_mp2(z, sym, xyz, use_ri=.false., use_scs=.false., freeze_core=.true., &
                             nelec=30, energy=whole, error=err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      call check(error, abs(res%energy - whole) < TOL, &
                 "FMO3 MP2 on the water trimer does not reproduce the supermolecule")
      if (allocated(error)) then
         write (*, *) "   fmo3        =", res%energy
         write (*, *) "   supermolecule=", whole
         write (*, *) "   difference  =", res%energy - whole
      end if
   end subroutine test_trimer_mp2

   subroutine test_eembe_dimer_mp2(error)
      !! Two waters, EE-MBE (esp = "ptc", expansion = "mbe"), MP2
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)
      real(dp) :: whole

      call water_dimer(z, sym, xyz)
      call mp2_options(opts, use_ri=.false., use_scs=.false., freeze_core=.true.)
      opts%esp = "ptc"
      opts%expansion = "mbe"
      opts%level = 2

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the EE-MBE MP2 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_mp2(z, sym, xyz, use_ri=.false., use_scs=.false., freeze_core=.true., &
                             nelec=20, energy=whole, error=err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      call check(error, abs(res%energy - whole) < TOL, &
                 "EE-MBE MP2 on the water dimer does not reproduce the supermolecule")
      if (allocated(error)) then
         write (*, *) "   ee-mbe      =", res%energy
         write (*, *) "   supermolecule=", whole
         write (*, *) "   difference  =", res%energy - whole
      end if
   end subroutine test_eembe_dimer_mp2

   subroutine test_dimer_mp2_no_frozen_core(error)
      !! The dimer again, every orbital correlated: `freeze_core = .false.`
      !! has to match its own (all-electron) supermolecule just as the frozen
      !! default matches its own in `test_dimer_mp2`.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)
      real(dp) :: whole

      call water_dimer(z, sym, xyz)
      call mp2_options(opts, use_ri=.false., use_scs=.false., freeze_core=.false.)
      opts%level = 2

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO2 MP2 run (no frozen core) failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_mp2(z, sym, xyz, use_ri=.false., use_scs=.false., &
                             freeze_core=.false., nelec=20, energy=whole, error=err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      call check(error, abs(res%energy - whole) < TOL, &
                 "FMO2 MP2 on the water dimer with every orbital correlated does not "// &
                 "reproduce the supermolecule")
      if (allocated(error)) then
         write (*, *) "   fmo2        =", res%energy
         write (*, *) "   supermolecule=", whole
         write (*, *) "   difference  =", res%energy - whole
      end if
   end subroutine test_dimer_mp2_no_frozen_core

   subroutine test_plain_mbe_dimer_mp2(error)
      !! `esp = "none"`, `expansion = "mbe"`: plain MBE, no field, MP2
      !!
      !! The monomer's one solve is also its only correlated one (there is no
      !! outer loop to converge first), and it has to still be the identity's
      !! only run: `add_monomer_correlation`'s `isolated` branch.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)
      real(dp) :: whole

      call water_dimer(z, sym, xyz)
      call mp2_options(opts, use_ri=.false., use_scs=.false., freeze_core=.true.)
      opts%esp = "none"
      opts%expansion = "mbe"
      opts%level = 2

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the plain-MBE MP2 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_mp2(z, sym, xyz, use_ri=.false., use_scs=.false., freeze_core=.true., &
                             nelec=20, energy=whole, error=err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call check(error, abs(res%energy - whole) < TOL, &
                 "plain-MBE MP2 on the water dimer does not reproduce the supermolecule")
      if (allocated(error)) then
         write (*, *) "   mbe2        =", res%energy
         write (*, *) "   supermolecule=", whole
         write (*, *) "   difference  =", res%energy - whole
      end if
   end subroutine test_plain_mbe_dimer_mp2

   subroutine test_trimer_cyclic_mp2_gamess(error)
      !! FMO2-MP2, cyclic water trimer, 6-31G, exact field, frozen core,
      !! against GAMESS
      !!
      !! GAMESS deck `tools/fmo_validation/gamess/w3_mp2.inp` (`mplevl=2`,
      !! `$fmo ... respap=0 resppc=0 resdim=1000`, GAMESS's default frozen
      !! core: 1 core orbital per water). The total, "The best FMO energy", is
      !! -228.375055246 Hartree. The pair terms are read off the log's
      !! per-fragment/per-dimer `EFMOc` (correlated internal energy) and `Tr`
      !! lines at full precision, `EFMOc(IJ) - EFMOc(I) - EFMOc(J) + Tr`:
      !! -0.015671582 (1-2), -0.015689099 (1-3), -0.014026446 (2-3) Hartree.
      !! This build reproduces the total to 2.1e-8 and every pair to 9e-10 --
      !! roughly the same agreement HF gets against GAMESS elsewhere in this
      !! suite, since the MP2 family adds correlation on the already-embedded
      !! Hartree-Fock reference (Decision 1, `developer_fragment_solver.rst`)
      !! rather than introducing a new source of cross-code noise.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      real(dp), parameter :: REF_TOTAL = -228.375055246_dp
      real(dp), parameter :: REF_PAIR(3) = [-0.015671582_dp, -0.015689099_dp, -0.014026446_dp]
      real(dp), parameter :: TOL_TOTAL = 1.0e-7_dp
      real(dp), parameter :: TOL_PAIR = 5.0e-9_dp
      integer :: k

      call water_trimer_cyclic(z, sym, xyz)
      call mp2_options(opts, use_ri=.false., use_scs=.false., freeze_core=.true.)
      opts%basis = "6-31g"
      opts%level = 2
      opts%resppc = -1.0_dp  ! exact field everywhere, matching GAMESS's resppc=0/respap=0
      opts%resdim = 0.0_dp
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      opts%outer_tol = 1.0e-10_dp

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO2 MP2 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call check(error, abs(res%energy - REF_TOTAL) < TOL_TOTAL, &
                 "FMO2-MP2 water trimer total does not match GAMESS")
      if (allocated(error)) then
         write (*, *) "   mqc   =", res%energy
         write (*, *) "   GAMESS=", REF_TOTAL
         write (*, *) "   diff  =", res%energy - REF_TOTAL
         return
      end if

      do k = 1, size(res%pairs)
         call check(error, abs(res%pairs(k)%energy - REF_PAIR(k)) < TOL_PAIR, &
                    "FMO2-MP2 water trimer pair IFIE does not match GAMESS")
         if (allocated(error)) then
            write (*, *) "   pair  =", res%pairs(k)%i, res%pairs(k)%j
            write (*, *) "   mqc   =", res%pairs(k)%energy
            write (*, *) "   GAMESS=", REF_PAIR(k)
            return
         end if
      end do
   end subroutine test_trimer_cyclic_mp2_gamess

   subroutine test_trimer_cyclic_separated_pair(error)
      !! A pair beyond `resdim` carries no MP2 correlation, HF or MP2 alike
      !!
      !! `water_trimer` (stacked 2.9 A apart, not the hydrogen-bonded ring) at
      !! GAMESS's FMO2 default `RESDIM` (2.0 -- `fmo_options_t%resdim` itself
      !! defaults to 0, solve every pair, and only the JSON deck layer applies
      !! GAMESS's FMO2 default when a deck says nothing) separates its 1-3
      !! pair (5.27 A apart, `check_fmo`'s own kind of measurement). Per
      !! Decision 1
      !! (`developer_fragment_solver.rst`) and GAMESS's own ES-dimer formula
      !! (`mqc_docs/source/fmo.rst`), a separated pair's term is the
      !! electrostatic interaction of its two converged Hartree-Fock monomers
      !! with no response term and no correlation -- the same number whether
      !! `model.method` is `hf` or `mp2`. This checks that identity directly,
      !! since it is what the design claims and does not itself need a fresh
      !! GAMESS run: `mqc_docs/source/fmo.rst` already pins mqc's separated-pair
      !! formula against GAMESS's `RESDIM=2.0` to 5e-9 on a twenty-water
      !! system, and that agreement does not depend on which reference method
      !! produced the monomer densities.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts_hf, opts_mp2
      type(fmo_result_t) :: res_hf, res_mp2
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      real(dp), parameter :: TOL = 1.0e-8_dp
      integer :: k, k13

      call water_trimer(z, sym, xyz)

      opts_hf%basis = "6-31g"
      opts_hf%level = 2
      opts_hf%resdim = 2.0_dp  ! GAMESS's FMO2 default RESDIM, which
      ! fmo_options_t itself does not default to
      opts_hf%scf_max_iter = 200
      opts_hf%scf_energy_tol = 1.0e-11_dp
      opts_hf%scf_density_tol = 1.0e-9_dp
      opts_hf%outer_tol = 1.0e-10_dp

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts_hf, res_hf, err)
      call check(error,.not. err%has_error(), "the FMO2 HF run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call mp2_options(opts_mp2, use_ri=.false., use_scs=.false., freeze_core=.true.)
      opts_mp2%basis = "6-31g"
      opts_mp2%level = 2
      opts_mp2%resdim = 2.0_dp
      opts_mp2%outer_tol = 1.0e-10_dp

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts_mp2, res_mp2, err)
      call check(error,.not. err%has_error(), "the FMO2 MP2 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      k13 = 0
      do k = 1, size(res_hf%pairs)
         if (res_hf%pairs(k)%i == 1 .and. res_hf%pairs(k)%j == 3) k13 = k
      end do
      call check(error, k13 > 0, "pair (1,3) was not found")
      if (allocated(error)) return

      call check(error, res_hf%pairs(k13)%separated, &
                 "pair (1,3) was expected to be separated beyond resdim")
      if (allocated(error)) return
      call check(error, res_mp2%pairs(k13)%separated, &
                 "pair (1,3) was expected to be separated beyond resdim under MP2 too")
      if (allocated(error)) return

      call check(error, abs(res_hf%pairs(k13)%energy - res_mp2%pairs(k13)%energy) < TOL, &
                 "a separated pair's term should not move when correlation is added")
      if (allocated(error)) then
         write (*, *) "   HF term to  =", res_hf%pairs(k13)%energy
         write (*, *) "   MP2 term to =", res_mp2%pairs(k13)%energy
         write (*, *) "   difference  =", res_hf%pairs(k13)%energy - res_mp2%pairs(k13)%energy
      end if
   end subroutine test_trimer_cyclic_separated_pair

   subroutine test_trimer_pair_truncated_mp2(error)
      !! The water trimer truncated at pairs, MP2, 6-31G -- printed, not
      !! asserted
      !!
      !! Not a full-level identity: with three fragments, level two omits the
      !! three-body term, so FMO2-MP2 and EE-MBE2-MP2 both differ from the
      !! supermolecule by that term (plus, for EE-MBE, the point-charge
      !! approximation). Printed to twelve decimals so an independent
      !! reference -- PySCF fed this repository's own 6-31G basis JSON, at the
      !! geometry `water_trimer` below, frozen core on (`core_orbital_count`:
      !! one 1s orbital per oxygen) -- can be checked against these numbers
      !! without rerunning this build.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res_fmo, res_eembe
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      real(dp) :: whole

      call water_trimer(z, sym, xyz)

      call mp2_options(opts, use_ri=.false., use_scs=.false., freeze_core=.true.)
      opts%basis = "6-31g"
      opts%level = 2
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res_fmo, err)
      call check(error,.not. err%has_error(), "the FMO2 MP2 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      opts%esp = "ptc"
      opts%expansion = "mbe"
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res_eembe, err)
      call check(error,.not. err%has_error(), "the EE-MBE2 MP2 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_mp2(z, sym, xyz, use_ri=.false., use_scs=.false., freeze_core=.true., &
                             nelec=30, energy=whole, error=err, basis_override="6-31g")
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      write (*, "(a)") "   water trimer, MP2/6-31G, frozen core on, geometry from "// &
         "water_trimer() below:"
      write (*, "(a, f20.12)") "   supermolecule       =", whole
      write (*, "(a, f20.12)") "   FMO2-MP2            =", res_fmo%energy
      write (*, "(a, f20.12)") "   FMO2-MP2 - super     =", res_fmo%energy - whole
      write (*, "(a, f20.12)") "   EE-MBE2-MP2          =", res_eembe%energy
      write (*, "(a, f20.12)") "   EE-MBE2-MP2 - super  =", res_eembe%energy - whole
   end subroutine test_trimer_pair_truncated_mp2

   subroutine test_eembe2_trimer_pyscf(error)
      !! EE-MBE2 on the water trimer, MP2, RI-MP2 and SCS-MP2, against PySCF
      !!
      !! The references come from `tools/fmo_validation/eembe_pyscf.py`, an
      !! independent reimplementation of EE-MBE fed this repository's 6-31G
      !! and def2-universal-jkfit JSON. It reproduces this code's Hartree-Fock
      !! total to 3e-11 Eh. Unlike a full-level identity, which the embedding
      !! cancels out of, this checks MP2 on embedded monomers and pairs. The
      !! tolerance is set by the outer charge loop's convergence (`outer_tol`,
      !! 1e-7 on the monomer sum), which the two codes stop at independently.
      !! Frozen core on.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: TOL = 5.0e-8_dp
      real(dp), parameter :: REF_MP2 = -228.354417062659_dp
      real(dp), parameter :: REF_RI_MP2 = -228.354367521908_dp
      real(dp), parameter :: REF_SCS_MP2 = -228.353055036908_dp

      call eembe2_trimer_against(.false., .false., REF_MP2, TOL, "MP2", error)
      if (allocated(error)) return
      call eembe2_trimer_against(.true., .false., REF_RI_MP2, TOL, "RI-MP2", error)
      if (allocated(error)) return
      call eembe2_trimer_against(.false., .true., REF_SCS_MP2, TOL, "SCS-MP2", error)
   end subroutine test_eembe2_trimer_pyscf

   subroutine eembe2_trimer_against(use_ri, use_scs, reference, tol, label, error)
      !! One EE-MBE2 water-trimer run in 6-31G, checked against `reference`
      logical, intent(in) :: use_ri, use_scs
      real(dp), intent(in) :: reference, tol
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)

      call water_trimer(z, sym, xyz)
      call mp2_options(opts, use_ri=use_ri, use_scs=use_scs, freeze_core=.true.)
      opts%basis = "6-31g"
      opts%method%aux_basis_set = "def2-universal-jkfit"
      opts%level = 2
      opts%esp = "ptc"
      opts%expansion = "mbe"
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "the EE-MBE2 "//label//" run failed: "// &
                 err%get_message())
      if (allocated(error)) return
      write (*, "(a, a, a, es12.4)") "   EE-MBE2-", label, " - PySCF = ", res%energy - reference
      call check(error, abs(res%energy - reference) < tol, &
                 "EE-MBE2 "//label//" is not the PySCF reference")
   end subroutine eembe2_trimer_against

   subroutine mp2_options(opts, use_ri, use_scs, freeze_core)
      !! What every FMO/EE-MBE run in this file shares: the SCF settings
      !! `test_mqc_fmo_dft` uses, plus MP2 (or one of its variants) as the
      !! fragment method
      type(fmo_options_t), intent(out) :: opts
      logical, intent(in) :: use_ri, use_scs, freeze_core

      type(correlation_config_t) :: default_corr

      opts%basis = BASIS
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      opts%method%run_mp2 = .true.
      opts%method%corr_density_fitting = use_ri
      opts%method%aux_basis_set = AUX
      opts%method%freeze_core = freeze_core
      if (use_scs) then
         opts%method%scs_ss = default_corr%scs_ss
         opts%method%scs_os = default_corr%scs_os
      end if
   end subroutine mp2_options

   subroutine supermolecule_mp2(z, sym, xyz, use_ri, use_scs, freeze_core, nelec, energy, &
                                error, basis_override)
      !! An ordinary Hartree-Fock-plus-MP2 energy on the whole system, for the
      !! full-level identity to be checked against
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: xyz(:, :)
      logical, intent(in) :: use_ri, use_scs, freeze_core
      integer, intent(in) :: nelec
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error
      character(len=*), intent(in), optional :: basis_override   !! Defaults to `BASIS`

      type(czt_molecule_t) :: mol, aux_mol
      type(rhf_result_t) :: scf
      type(mp2_result_t) :: mp2
      type(correlation_config_t) :: default_corr
      character(len=:), allocatable :: this_basis
      real(dp) :: ss, os
      integer :: frozen

      energy = 0.0_dp
      this_basis = BASIS
      if (present(basis_override)) this_basis = basis_override
      ss = 1.0_dp
      os = 1.0_dp
      if (use_scs) then
         ss = default_corr%scs_ss
         os = default_corr%scs_os
      end if
      frozen = 0
      if (freeze_core) frozen = core_orbital_count(z)

      call build_czt_molecule(z, sym, xyz, this_basis, mol, error)
      if (error%has_error()) return
      call run_czt_rhf(mol, nelec, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, error)
      if (error%has_error()) return

      if (use_ri) then
         call build_czt_molecule(z, sym, xyz, AUX, aux_mol, error)
         if (error%has_error()) return
         call run_czt_ri_mp2(mol, aux_mol, scf%orbitals, scf%orbital_energies, nelec/2, &
                             scf%energy, mp2, error, n_frozen=frozen)
         call aux_mol%destroy()
      else
         call run_czt_mp2(mol, scf%orbitals, scf%orbital_energies, nelec/2, scf%energy, &
                          mp2, error, n_frozen=frozen)
      end if
      if (error%has_error()) return
      energy = scf%energy + ss*mp2%same_spin + os*mp2%opposite_spin
   end subroutine supermolecule_mp2

   subroutine water_dimer(z, sym, xyz)
      !! Two waters, stacked 2.9 A apart -- the first pair of `water_trimer`,
      !! as in `test_mqc_fmo_dft`
      integer, intent(out) :: z(6)
      character(len=2), intent(out) :: sym(6)
      real(dp), intent(out) :: xyz(3, 6)
      real(dp) :: ang(3, 6)

      z = [8, 1, 1, 8, 1, 1]
      sym = ["O ", "H ", "H ", "O ", "H ", "H "]
      ang = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                     0.0_dp, -0.7572_dp, 0.5865_dp, &
                     0.0_dp, 0.7572_dp, 0.5865_dp, &
                     0.0_dp, 0.0_dp, 2.9_dp, &
                     0.0_dp, -0.7572_dp, 3.4865_dp, &
                     0.0_dp, 0.7572_dp, 3.4865_dp], [3, 6])
      xyz = to_bohr(ang)
   end subroutine water_dimer

   subroutine water_trimer(z, sym, xyz)
      !! `sample_inputs/w3.xyz`, three waters stacked 2.9 A apart, as in
      !! `test_mqc_fmo_dft`
      integer, intent(out) :: z(9)
      character(len=2), intent(out) :: sym(9)
      real(dp), intent(out) :: xyz(3, 9)
      real(dp) :: ang(3, 9)

      z = [8, 1, 1, 8, 1, 1, 8, 1, 1]
      sym = ["O ", "H ", "H ", "O ", "H ", "H ", "O ", "H ", "H "]
      ang = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                     0.0_dp, -0.7572_dp, 0.5865_dp, &
                     0.0_dp, 0.7572_dp, 0.5865_dp, &
                     0.0_dp, 0.0_dp, 2.9_dp, &
                     0.0_dp, -0.7572_dp, 3.4865_dp, &
                     0.0_dp, 0.7572_dp, 3.4865_dp, &
                     0.0_dp, 0.0_dp, 5.8_dp, &
                     0.0_dp, -0.7572_dp, 6.3865_dp, &
                     0.0_dp, 0.7572_dp, 6.3865_dp], [3, 9])
      xyz = to_bohr(ang)
   end subroutine water_trimer

   subroutine water_trimer_cyclic(z, sym, xyz)
      !! `sample_inputs/water3_cyclic.xyz`, the hydrogen-bonded ring GAMESS's
      !! own `3h2o.pieda.inp` optimised, RHF/6-31G*, as in `test_mqc_fmo_dft`
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
   end subroutine water_trimer_cyclic

end module test_mqc_fmo_mp2

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fmo_mp2, only: collect_mqc_fmo_mp2
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fmo_mp2", collect_mqc_fmo_mp2)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
