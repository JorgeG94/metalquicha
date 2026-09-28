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
                  new_unittest("fmo2_water_trimer_pbe_difference_from_the_supermolecule", &
                               test_trimer_fmo2_difference) &
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
