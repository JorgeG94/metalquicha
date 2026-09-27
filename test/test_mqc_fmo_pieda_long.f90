!! PIEDA on glycine tripeptide with a water, against PySCF
module test_mqc_fmo_pieda_long
   !! Gate 4 of `PIEDA_LAYER3_DESIGN.md`, split from `test_mqc_fmo_pieda` and
   !! labelled `LONG`: the tripeptide's monomer SCFs at tight tolerances take
   !! several minutes, against seconds for every other PIEDA case.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: to_bohr
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   implicit none
   private

   public :: collect_mqc_fmo_pieda_long

   real(dp), parameter :: REF_TOL = 2.0e-9_dp
      !! Against the PySCF reference, as in `test_mqc_fmo_pieda`

contains

   subroutine collect_mqc_fmo_pieda_long(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("glycine_tripeptide_and_water_pieda_matches_the_reference", &
                               test_glycine_water) &
                  ]
   end subroutine collect_mqc_fmo_pieda_long

   subroutine test_glycine_water(error)
      !! Gate 4: glycine tripeptide with a water, no cut, HF/6-31G
      !!
      !! Reference from `gen_pieda_refs.py` on the exact geometry `gly3_water`
      !! below builds -- **not** `sample_inputs/gly3_water.xyz`, whose water is
      !! placed differently; this one is `sample_inputs/gly3_water_pair.xyz`,
      !! the file `test_mqc_fmo_pairs.f90`'s own `gly3_water` docstring names:
      !!
      !!   python3 tools/cpu_validation/gen_pieda_refs.py gly3_water_pair.xyz 6-31g glyw \
      !!     '[[0,...,23],[24,25,26]]'
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: REF_EES = -0.017416031579_dp
      real(dp), parameter :: REF_EEX = 0.008601981502_dp
      real(dp), parameter :: REF_ECT = -0.002649651634_dp
      type(fmo_result_t) :: res
      type(error_t) :: err
      type(fmo_options_t) :: opts
      integer :: z(27)
      character(len=2) :: sym(27)
      real(dp) :: xyz(3, 27)
      integer :: owner(27)

      call gly3_water(z, sym, xyz)
      owner(1:24) = 1
      owner(25:27) = 2
      opts%basis = "6-31g"
      opts%esp = "exact"
      opts%expansion = "fmo"
      opts%resppc = -1.0_dp
      opts%resdim = 0.0_dp
      ! Above the outer loop's OpenMP noise floor, which is a few 1e-10 here:
      ! at 1e-12 it wanders there for another ten iterations and gains nothing.
      opts%outer_tol = 1.0e-10_dp
      opts%max_outer = 100
      opts%scf_energy_tol = 1.0e-12_dp
      opts%scf_density_tol = 1.0e-10_dp
      opts%pieda = .true.
      call run_fmo2(z, sym, xyz, owner, opts, res, err)
      call check(error,.not. err%has_error(), "glycine tripeptide with a water failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, size(res%pairs), 1, "there should be exactly one pair")
      if (allocated(error)) return
      call check(error, res%pairs(1)%pieda, "the pair should be decomposed")
      if (allocated(error)) return

      call check(error, abs(res%pairs(1)%ees - REF_EES) < REF_TOL, &
                 "Ees does not match the reference")
      if (allocated(error)) then
         write (*, *) "   Ees - ref =", res%pairs(1)%ees - REF_EES
         return
      end if
      call check(error, abs(res%pairs(1)%eex - REF_EEX) < REF_TOL, &
                 "Eex does not match the reference")
      if (allocated(error)) then
         write (*, *) "   Eex - ref =", res%pairs(1)%eex - REF_EEX
         return
      end if
      call check(error, abs(res%pairs(1)%ect_mix - REF_ECT) < REF_TOL, &
                 "Ect+mix does not match the reference")
      if (allocated(error)) write (*, *) "   Ect+mix - ref =", res%pairs(1)%ect_mix - REF_ECT
   end subroutine test_glycine_water

   subroutine gly3_water(z, sym, xyz)
      !! `sample_inputs/gly3_water.xyz`, in Bohr -- glycine tripeptide with a
      !! water donating an H-bond to the middle carbonyl O; the same geometry
      !! `test_mqc_fmo_pairs.f90`'s own `gly3_water` builds
      integer, intent(out) :: z(27)
      character(len=2), intent(out) :: sym(27)
      real(dp), intent(out) :: xyz(3, 27)
      real(dp) :: ang(3, 27)
      integer :: i

      sym = ["N ", "C ", "C ", "O ", "H ", "H ", "H ", "H ", "N ", "C ", "C ", "O ", &
             "H ", "H ", "H ", "N ", "C ", "C ", "O ", "H ", "H ", "H ", "O ", "H ", &
             "O ", "H ", "H "]
      do i = 1, 27
         select case (sym(i))
         case ("N ")
            z(i) = 7
         case ("C ")
            z(i) = 6
         case ("O ")
            z(i) = 8
         case default
            z(i) = 1
         end select
      end do
      ang = reshape([ &
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
                    6.6876213700_dp, -0.2710584500_dp, -1.2503709636_dp], [3, 27])
      xyz = to_bohr(ang)
   end subroutine gly3_water

end module test_mqc_fmo_pieda_long

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fmo_pieda_long, only: collect_mqc_fmo_pieda_long
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fmo_pieda_long", collect_mqc_fmo_pieda_long)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
