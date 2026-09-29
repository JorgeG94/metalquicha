!! Kohn-Sham fragments through the fragment solver, against the supermolecule
module test_mqc_fmo_dft
   !! The DFT counterpart of `mqc_afo_fmo`: with the fragment count equal to
   !! the expansion level, FMO and EE-MBE are exact by inclusion and
   !! exclusion whatever the field or the partition did (see
   !! `mqc_docs/source/developer_fragment_solver.rst`), so a Kohn-Sham total
   !! has to reproduce an ordinary Kohn-Sham calculation on the whole
   !! molecule -- the same identity `test_mqc_afo_fmo` checks for
   !! Hartree-Fock, carried over to a GGA, a global hybrid, a range-separated
   !! hybrid, EE-MBE and a detached covalent bond.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: to_bohr
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: run_czt_rhf, rhf_result_t
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   implicit none
   private

   public :: collect_mqc_fmo_dft

   real(dp), parameter :: TOL = 1.0e-8_dp
      !! Both sides are the same Kohn-Sham SCF on the same atoms, in the same
      !! basis and on the same grid, so what separates them is rounding in
      !! the fragment assembly rather than any physics.

   integer, parameter :: GRID_LEVEL = 3
      !! Passed identically to every fragment solve (through `opts%method`)
      !! and to every supermolecule reference built here, so the two sides
      !! integrate on the same grid. A full-level identity failing at 1e-6
      !! rather than holding to 1e-8 usually means the two grids differed,
      !! not that the physics did.

   integer, parameter :: GAMESS_GRID_LEVEL = 5
      !! For the tests checked against GAMESS below, where the other side is
      !! GAMESS's own grid (`nrad=200 nleb=1202`) rather than a grid this
      !! build also controls -- so a larger level than `GRID_LEVEL` narrows
      !! the quadrature disagreement instead of matching it exactly. 770
      !! angular points on oxygen against GAMESS's 1202.

contains

   subroutine collect_mqc_fmo_dft(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fmo2_water_dimer_pbe_is_the_supermolecule", test_dimer_pbe), &
                  new_unittest("fmo2_water_dimer_cam_b3lyp_is_the_supermolecule", &
                               test_dimer_cam_b3lyp), &
                  new_unittest("fmo3_water_trimer_b3lyp_is_the_supermolecule", &
                               test_trimer_b3lyp), &
                  new_unittest("eembe_water_dimer_pbe_is_the_supermolecule", test_eembe_dimer), &
                  new_unittest("propane_cut_across_one_bond_pbe_is_the_supermolecule", &
                               test_propane_afo), &
                  new_unittest("propane_in_three_fragments_pbe_freezes_a_kohn_sham_orbital", &
                               test_propane_afo_three), &
                  new_unittest("fmo2_water_trimer_pbe_difference_from_the_supermolecule", &
                               test_trimer_fmo2_difference), &
                  new_unittest("fmo2_water_trimer_cyclic_pbe_matches_gamess", &
                               test_trimer_cyclic_pbe_gamess), &
                  new_unittest("fmo2_water_trimer_cyclic_b3lyp_matches_gamess", &
                               test_trimer_cyclic_b3lyp_gamess) &
                  ]
   end subroutine collect_mqc_fmo_dft

   subroutine test_dimer_pbe(error)
      !! Two waters, FMO2 (the default field), PBE
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)
      real(dp) :: whole

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call water_dimer(z, sym, xyz)

      opts%basis = "sto-3g"
      opts%level = 2
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      opts%method%functional = "pbe"
      opts%method%grid_level = GRID_LEVEL

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO2 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_energy(z, sym, xyz, "sto-3g", "pbe", GRID_LEVEL, 20, whole, err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      call check(error, abs(res%energy - whole) < TOL, &
                 "FMO2 PBE on the water dimer does not reproduce the supermolecule")
      if (allocated(error)) then
         write (*, *) "   fmo2        =", res%energy
         write (*, *) "   supermolecule=", whole
         write (*, *) "   difference  =", res%energy - whole
      end if
   end subroutine test_dimer_pbe

   subroutine test_dimer_cam_b3lyp(error)
      !! The same dimer, CAM-B3LYP: the long-range exchange has to reach
      !! through the embedding operator the same way the short-range part does
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)
      real(dp) :: whole

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call water_dimer(z, sym, xyz)

      opts%basis = "sto-3g"
      opts%level = 2
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      opts%method%functional = "cam-b3lyp"
      opts%method%grid_level = GRID_LEVEL

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO2 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_energy(z, sym, xyz, "sto-3g", "cam-b3lyp", GRID_LEVEL, 20, whole, err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      call check(error, abs(res%energy - whole) < TOL, &
                 "FMO2 CAM-B3LYP on the water dimer does not reproduce the supermolecule")
      if (allocated(error)) then
         write (*, *) "   fmo2        =", res%energy
         write (*, *) "   supermolecule=", whole
         write (*, *) "   difference  =", res%energy - whole
      end if
   end subroutine test_dimer_cam_b3lyp

   subroutine test_trimer_b3lyp(error)
      !! Three waters, FMO3 (level equal to the fragment count), B3LYP
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      real(dp) :: whole

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call water_trimer(z, sym, xyz)

      opts%basis = "sto-3g"
      opts%level = 3
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      opts%method%functional = "b3lyp"
      opts%method%grid_level = GRID_LEVEL

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO3 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_energy(z, sym, xyz, "sto-3g", "b3lyp", GRID_LEVEL, 30, whole, err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      call check(error, abs(res%energy - whole) < TOL, &
                 "FMO3 B3LYP on the water trimer does not reproduce the supermolecule")
      if (allocated(error)) then
         write (*, *) "   fmo3        =", res%energy
         write (*, *) "   supermolecule=", whole
         write (*, *) "   difference  =", res%energy - whole
      end if
   end subroutine test_trimer_b3lyp

   subroutine test_eembe_dimer(error)
      !! Two waters, EE-MBE (esp = "ptc", expansion = "mbe"), PBE
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)
      real(dp) :: whole

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call water_dimer(z, sym, xyz)

      opts%basis = "sto-3g"
      opts%esp = "ptc"
      opts%expansion = "mbe"
      opts%level = 2
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      opts%method%functional = "pbe"
      opts%method%grid_level = GRID_LEVEL

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the EE-MBE run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_energy(z, sym, xyz, "sto-3g", "pbe", GRID_LEVEL, 20, whole, err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      call check(error, abs(res%energy - whole) < TOL, &
                 "EE-MBE PBE on the water dimer does not reproduce the supermolecule")
      if (allocated(error)) then
         write (*, *) "   ee-mbe      =", res%energy
         write (*, *) "   supermolecule=", whole
         write (*, *) "   difference  =", res%energy - whole
      end if
   end subroutine test_eembe_dimer

   subroutine test_propane_afo(error)
      !! Propane cut into a methyl and an ethyl across one C-C bond, PBE
      !!
      !! Two fragments across a detached bond is exact for the same reason
      !! `two_fragments_across_a_cut_bond_are_exact` is in `test_mqc_afo_fmo`:
      !! the dimer holds both ends of the cut, so it carries no ghost, no
      !! frozen orbital and no electron shift, and the monomer terms cancel
      !! against the pair correction exactly.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(11)
      character(len=2) :: sym(11)
      real(dp) :: xyz(3, 11)
      real(dp) :: whole

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call propane(z, sym, xyz)

      opts%basis = "sto-3g"
      opts%bond_breaking = "afo"
      opts%level = 2
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      opts%method%functional = "pbe"
      opts%method%grid_level = GRID_LEVEL

      ! Methyl (atom 1 and its three hydrogens) against ethyl (atoms 2-3 and
      ! their six hydrogens), cutting the C1-C2 bond.
      call run_fmo2(z, sym, xyz, [1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the frozen-orbital expansion failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_energy(z, sym, xyz, "sto-3g", "pbe", GRID_LEVEL, 26, whole, err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      call check(error, abs(res%energy - whole) < TOL, &
                 "the detached-bond expansion in PBE does not reproduce the whole "// &
                 "molecule")
      if (allocated(error)) then
         write (*, *) "   fmo2        =", res%energy
         write (*, *) "   supermolecule=", whole
         write (*, *) "   difference  =", res%energy - whole
      end if
   end subroutine test_propane_afo

   subroutine test_propane_afo_three(error)
      !! Propane in three fragments, two detached bonds, FMO2-PBE: the frozen
      !! orbitals come from a model solved at PBE
      !!
      !! Below full level, so the embedding and the frozen orbitals matter and
      !! the result is not the supermolecule. The model system under each cut
      !! is solved at the deck's functional, as GAMESS does, where it used to be
      !! Hartree-Fock by construction. `E_HF_MODEL` is what this same run gave
      !! before that change, obtained by leaving `afo_opts%method` unallocated
      !! in `build_afo_context`; `E_KS_MODEL` is the value after it, recorded
      !! here. The two differ because the frozen orbital does. Neither is close
      !! to the supermolecule: the middle fragment is a CH2 held by two cuts,
      !! which is a poor partition at level two whatever the model. The
      !! comparison is between the models and not against the whole molecule.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(11)
      character(len=2) :: sym(11)
      real(dp) :: xyz(3, 11)
      real(dp) :: whole
      real(dp), parameter :: E_HF_MODEL = -117.31265747731965_dp
      real(dp), parameter :: E_KS_MODEL = -117.30068114551975_dp
      real(dp), parameter :: TOL_RECORDED = 1.0e-7_dp
      real(dp), parameter :: MIN_CHANGE = 1.0e-3_dp

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call propane(z, sym, xyz)

      opts%basis = "sto-3g"
      opts%bond_breaking = "afo"
      opts%level = 2
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      opts%method%functional = "pbe"
      opts%method%grid_level = GRID_LEVEL

      ! One carbon and its hydrogens per fragment, cutting C1-C2 and C2-C3.
      call run_fmo2(z, sym, xyz, [1, 2, 3, 1, 1, 1, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "the two-cut frozen-orbital expansion failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if
      call check(error, res%converged, "the two-cut PBE expansion did not converge")
      if (allocated(error)) return

      call supermolecule_energy(z, sym, xyz, "sto-3g", "pbe", GRID_LEVEL, 26, whole, err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      call check(error, abs(res%energy - E_KS_MODEL) < TOL_RECORDED, &
                 "the two-cut PBE expansion does not reproduce its recorded value")
      if (allocated(error)) then
         write (*, *) "   fmo2     =", res%energy
         write (*, *) "   recorded =", E_KS_MODEL
         return
      end if
      call check(error, abs(res%energy - E_HF_MODEL) > MIN_CHANGE, &
                 "the two-cut PBE expansion is what a Hartree-Fock model gave: the "// &
                 "model's method was ignored")
      if (allocated(error)) return

      write (*, *) "   FMO2 PBE, three fragments      =", res%energy
      write (*, *) "   supermolecule                  =", whole
      write (*, *) "   error against the supermolecule =", res%energy - whole
      write (*, *) "   error with an HF model          =", E_HF_MODEL - whole
      write (*, *) "   change from the HF model        =", res%energy - E_HF_MODEL
   end subroutine test_propane_afo_three

   subroutine test_trimer_fmo2_difference(error)
      !! The water trimer truncated at pairs, PBE -- printed, not asserted
      !!
      !! Not a full-level identity: with three fragments, level two omits the
      !! three-body term, so this is the approximation `three_fragment_error`
      !! in `test_mqc_afo_fmo` measures for Hartree-Fock, carried over to a
      !! functional.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      real(dp) :: whole

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call water_trimer(z, sym, xyz)

      opts%basis = "sto-3g"
      opts%level = 2
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      opts%method%functional = "pbe"
      opts%method%grid_level = GRID_LEVEL

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO2 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call supermolecule_energy(z, sym, xyz, "sto-3g", "pbe", GRID_LEVEL, 30, whole, err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      write (*, *) "   FMO2 - supermolecule (water trimer, PBE, sto-3g) =", res%energy - whole
   end subroutine test_trimer_fmo2_difference

   subroutine test_trimer_cyclic_pbe_gamess(error)
      !! FMO2-PBE, cyclic water trimer, 6-31G, exact field, against GAMESS
      !!
      !! GAMESS deck `tools/fmo_validation/gamess/w3_pbe.inp`
      !! (`$fmo ... respap=0 resppc=0 resdim=0`, `$dft nrad=200 nleb=1202`,
      !! `dfttyp=pbe`), same geometry as `water_trimer_cyclic` below. GAMESS's
      !! total is -228.942512413 Hartree ("The best FMO energy"); the three
      !! pair terms, `EFMOu(IJ) - EFMOu(I) - EFMOu(J) + Tr`, read off the
      !! log's per-fragment/per-dimer `EFMOu`/`Tr` lines at full precision, are
      !! -0.018353617 (1-2), -0.018277205 (1-3) and -0.016561905 (2-3) Hartree.
      !! At `GAMESS_GRID_LEVEL` this build reproduces the total to 1.2e-8 and
      !! every pair to 1.4e-7 -- grid-quadrature noise, not a disagreement in
      !! the physics (see
      !! `mqc_docs/source/developer_fragment_solver.rst`, phase 2, and the
      !! module docstring's remark that a full-level identity is blind to the
      !! embedding).
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      real(dp), parameter :: REF_TOTAL = -228.942512413_dp
      real(dp), parameter :: REF_PAIR(3) = [-0.018353617_dp, -0.018277205_dp, -0.016561905_dp]
         !! Pairs (1,2), (1,3), (2,3), the order `run_fmo2` enumerates them
      real(dp), parameter :: TOL_TOTAL = 2.0e-7_dp
      real(dp), parameter :: TOL_PAIR = 5.0e-7_dp
      integer :: k

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call water_trimer_cyclic(z, sym, xyz)

      opts%basis = "6-31g"
      opts%level = 2
      opts%resppc = -1.0_dp  ! exact field everywhere, matching GAMESS's resppc=0/respap=0
      opts%resdim = 0.0_dp
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      opts%outer_tol = 1.0e-10_dp
      opts%method%functional = "pbe"
      opts%method%grid_level = GAMESS_GRID_LEVEL

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO2 PBE run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call check(error, abs(res%energy - REF_TOTAL) < TOL_TOTAL, &
                 "FMO2-PBE water trimer total does not match GAMESS")
      if (allocated(error)) then
         write (*, *) "   mqc   =", res%energy
         write (*, *) "   GAMESS=", REF_TOTAL
         write (*, *) "   diff  =", res%energy - REF_TOTAL
         return
      end if

      do k = 1, size(res%pairs)
         call check(error, abs(res%pairs(k)%energy - REF_PAIR(k)) < TOL_PAIR, &
                    "FMO2-PBE water trimer pair IFIE does not match GAMESS")
         if (allocated(error)) then
            write (*, *) "   pair  =", res%pairs(k)%i, res%pairs(k)%j
            write (*, *) "   mqc   =", res%pairs(k)%energy
            write (*, *) "   GAMESS=", REF_PAIR(k)
            return
         end if
      end do
   end subroutine test_trimer_cyclic_pbe_gamess

   subroutine test_trimer_cyclic_b3lyp_gamess(error)
      !! FMO2, cyclic water trimer, 6-31G, exact field, against GAMESS B3LYP
      !!
      !! GAMESS's plain `dfttyp=b3lyp` uses VWN formula V for the local
      !! correlation (GAMESS's functional table, `dftxca.src`), not
      !! VWN-RPA; libxc's `hyb_gga_xc_b3lyp` is the
      !! VWN-RPA variant (`XC_LDA_C_VWN_RPA` in `hyb_gga_xc_b3lyp.c`) and
      !! `hyb_gga_xc_b3lyp5` is the VWN5 one GAMESS's name means
      !! (`XC_LDA_C_VWN`, "B3LYP with VWN functional 5 instead of RPA"). This
      !! test asks for `hyb_gga_xc_b3lyp5` to match GAMESS's `B3LYP`, not
      !! mqc's own `"b3lyp"` alias.
      !!
      !! GAMESS deck `tools/fmo_validation/gamess/w3_b3lyp.inp`, same
      !! settings and geometry as `test_trimer_cyclic_pbe_gamess`. GAMESS's
      !! total is -229.088616246 Hartree; the pair terms are -0.017106306
      !! (1-2), -0.017073134 (1-3) and -0.015386181 (2-3) Hartree. This build
      !! reproduces the total to 5.7e-8 and every pair to 1.4e-7.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      real(dp), parameter :: REF_TOTAL = -229.088616246_dp
      real(dp), parameter :: REF_PAIR(3) = [-0.017106306_dp, -0.017073134_dp, -0.015386181_dp]
      real(dp), parameter :: TOL_TOTAL = 2.0e-7_dp
      real(dp), parameter :: TOL_PAIR = 5.0e-7_dp
      integer :: k

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call water_trimer_cyclic(z, sym, xyz)

      opts%basis = "6-31g"
      opts%level = 2
      opts%resppc = -1.0_dp
      opts%resdim = 0.0_dp
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      opts%outer_tol = 1.0e-10_dp
      opts%method%functional = "hyb_gga_xc_b3lyp5"
      opts%method%grid_level = GAMESS_GRID_LEVEL

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO2 B3LYP5 run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call check(error, abs(res%energy - REF_TOTAL) < TOL_TOTAL, &
                 "FMO2-B3LYP5 water trimer total does not match GAMESS's B3LYP")
      if (allocated(error)) then
         write (*, *) "   mqc   =", res%energy
         write (*, *) "   GAMESS=", REF_TOTAL
         write (*, *) "   diff  =", res%energy - REF_TOTAL
         return
      end if

      do k = 1, size(res%pairs)
         call check(error, abs(res%pairs(k)%energy - REF_PAIR(k)) < TOL_PAIR, &
                    "FMO2-B3LYP5 water trimer pair IFIE does not match GAMESS's B3LYP")
         if (allocated(error)) then
            write (*, *) "   pair  =", res%pairs(k)%i, res%pairs(k)%j
            write (*, *) "   mqc   =", res%pairs(k)%energy
            write (*, *) "   GAMESS=", REF_PAIR(k)
            return
         end if
      end do
   end subroutine test_trimer_cyclic_b3lyp_gamess

   subroutine supermolecule_energy(z, sym, xyz, basis, functional, level, nelec, energy, error)
      !! An ordinary restricted Kohn-Sham energy on the whole system, for the
      !! full-level identity to be checked against
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: xyz(:, :)
      character(len=*), intent(in) :: basis
      character(len=*), intent(in) :: functional
      integer, intent(in) :: level
      integer, intent(in) :: nelec
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: mol
      type(xc_context_t) :: xc
      type(rhf_result_t) :: scf

      energy = 0.0_dp
      call build_czt_molecule(z, sym, xyz, basis, mol, error)
      if (error%has_error()) return
      call xc_context_create(mol, functional, xc, error, level=level, polarized=.false.)
      if (error%has_error()) return
      call run_czt_rhf(mol, nelec, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, error, xc=xc)
      call xc%destroy()
      if (error%has_error()) return
      energy = scf%energy
   end subroutine supermolecule_energy

   subroutine water_dimer(z, sym, xyz)
      !! Two waters, stacked 2.9 A apart -- the first pair of `water_trimer`
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
      !! `test_mqc_fmo_pairs`
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
      !! own `3h2o.pieda.inp` optimised, RHF/6-31G*: the geometry every
      !! GAMESS-referenced FMO2 test in this file and in `test_mqc_fmo_mp2`
      !! is run on.
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

   subroutine propane(z, sym, xyz)
      !! Idealised propane, carbons first -- as in `test_mqc_afo_fmo`
      integer, intent(out) :: z(11)
      character(len=2), intent(out) :: sym(11)
      real(dp), intent(out) :: xyz(3, 11)
      real(dp) :: ang(3, 11)
      integer :: i

      z = [6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1]
      do i = 1, 11
         if (z(i) == 6) then
            sym(i) = "C "
         else
            sym(i) = "H "
         end if
      end do
      ang = reshape([1.5260_dp, 0.0000_dp, 0.0000_dp, &
                     0.0000_dp, 0.0000_dp, 0.0000_dp, &
                     -0.5716_dp, 1.4149_dp, 0.0000_dp, &
                     2.1553_dp, -0.8900_dp, 0.0000_dp, &
                     2.1553_dp, 0.4450_dp, -0.7707_dp, &
                     2.1553_dp, 0.4450_dp, 0.7707_dp, &
                     -0.3519_dp, -0.5217_dp, 0.8900_dp, &
                     -0.3519_dp, -0.5217_dp, -0.8900_dp, &
                     0.0178_dp, 2.3318_dp, 0.0000_dp, &
                     -1.2200_dp, 1.8317_dp, -0.7707_dp, &
                     -1.2200_dp, 1.8317_dp, 0.7707_dp], [3, 11])
      xyz = to_bohr(ang)
   end subroutine propane

end module test_mqc_fmo_dft

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fmo_dft, only: collect_mqc_fmo_dft
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fmo_dft", collect_mqc_fmo_dft)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
