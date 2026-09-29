!! Kohn-Sham FMO across two detached peptide bonds, glycine tripeptide and water
module test_mqc_fmo_dft_long
   !! Split from `test_mqc_fmo_dft` and labelled `LONG`: fifteen PBE n-mer
   !! SCFs over the tripeptide take a few minutes, against seconds for every
   !! other Kohn-Sham FMO case.
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

   public :: collect_mqc_fmo_dft_long

   real(dp), parameter :: TOL = 1.0e-8_dp
      !! As in `test_mqc_fmo_dft`: the same Kohn-Sham SCF on the same atoms and
      !! grid on both sides

   integer, parameter :: GRID_LEVEL = 3
      !! Passed identically to every fragment solve and to the supermolecule

contains

   subroutine collect_mqc_fmo_dft_long(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("gly3_water_two_cuts_pbe_full_level_is_the_supermolecule", &
                               test_gly3_water_afo) &
                  ]
   end subroutine collect_mqc_fmo_dft_long

   subroutine test_gly3_water_afo(error)
      !! Glycine tripeptide cut at both C-alpha--C(=O) bonds, with a water, PBE
      !!
      !! Four fragments in STO-3G, the partition of `gly3_water` in
      !! `test_mqc_afo_fmo` and of GAMESS's `gly3w_afo_pbe.inp`. At pairs, each
      !! cut bond's model system is solved at PBE and must converge -- the case
      !! on which GAMESS's own PBE model SCF diverges, which is why cut
      !! FMO-DFT has no GAMESS reference -- and the pair error against the
      !! molecule is printed, not asserted. At full order the tetramer holds
      !! both cuts whole again, so the expansion is the molecule to `TOL`.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(27), owner(27), run
      character(len=2) :: sym(27)
      real(dp) :: xyz(3, 27)
      real(dp) :: whole
      integer, parameter :: LEVEL(2) = [2, 4]

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check

      call gly3_water(z, sym, xyz)
      owner = [1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 3, 3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4]

      call supermolecule_energy(z, sym, xyz, "sto-3g", "pbe", GRID_LEVEL, sum(z), whole, err)
      call check(error,.not. err%has_error(), "the supermolecule reference failed")
      if (allocated(error)) return

      do run = 1, size(LEVEL)
         opts%basis = "sto-3g"
         opts%bond_breaking = "afo"
         opts%level = LEVEL(run)
         opts%scf_max_iter = 200
         opts%scf_energy_tol = 1.0e-11_dp
         opts%scf_density_tol = 1.0e-9_dp
         opts%method%functional = "pbe"
         opts%method%grid_level = GRID_LEVEL
         call run_fmo2(z, sym, xyz, owner, opts, res, err)
         call check(error,.not. err%has_error(), "the two-cut PBE expansion failed")
         if (allocated(error)) then
            write (*, *) "   level ", LEVEL(run), ": ", trim(err%get_message())
            return
         end if
         call check(error, res%converged, "the two-cut PBE expansion did not converge")
         if (allocated(error)) return
         write (*, "(a,i2,a,es12.4)") "    gly3 + water, PBE, level", LEVEL(run), &
            ": error against the molecule (hartree) =", res%energy - whole
      end do

      call check(error, abs(res%energy - whole) < TOL, &
                 "full order over two detached peptide bonds in PBE is not the molecule")
      if (allocated(error)) then
         write (*, *) "   fmo4         =", res%energy
         write (*, *) "   supermolecule=", whole
      end if
   end subroutine test_gly3_water_afo

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

   subroutine gly3_water(z, sym, xyz)
      !! `validation/inputs/sample_inputs/gly3_water_pair.xyz`
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
   end subroutine gly3_water

end module test_mqc_fmo_dft_long

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fmo_dft_long, only: collect_mqc_fmo_dft_long
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fmo_dft_long", collect_mqc_fmo_dft_long)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
