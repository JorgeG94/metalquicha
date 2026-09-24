!! The per-pair terms an FMO run exports, against the sums they came out of
module test_mqc_fmo_pairs
   !! `fmo_result_t%pairs` keeps each two-member term of the expansion instead
   !! of summing it away. Nothing about the energy changes, so the tests are
   !! bookkeeping identities rather than physics:
   !!
   !! * at level two the pairs are every n-mer term there is, so they add up
   !!   to `pair_sum` -- summed in the order the expansion summed them, which
   !!   makes this exact rather than a tolerance;
   !! * at level three they do not, and the shortfall is the three-body sum
   !!   exactly. Asserted, so a pair map read off a level-three run is known to
   !!   be missing what it is missing;
   !! * the total is still `monomer_sum + pair_sum`, bit for bit, and still
   !!   the number the validation suite pins;
   !! * on a peptide with a water, the two pairs a detached bond joins are
   !!   flagged and every ligand row is present and unflagged.
   !!
   !! Distances are checked by hand on the water trimer, where the closest
   !! approach can be written down.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: to_bohr
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   implicit none
   private

   public :: collect_mqc_fmo_pairs

   real(dp), parameter :: SUM_TOL = 1.0e-12_dp
      !! How far a re-summed pair list may sit from the sum the expansion
      !! formed. Same terms, same order, so the measured gap is zero; the
      !! bound is there so a reordering that costs an ulp does not fail CI.

contains

   subroutine collect_mqc_fmo_pairs(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("level_two_pairs_sum_to_the_pair_sum", test_level_two), &
                  new_unittest("level_three_pairs_fall_short_by_the_three_body_sum", &
                               test_level_three), &
                  new_unittest("a_peptide_flags_its_joined_pairs_and_keeps_the_ligand", &
                               test_peptide) &
                  ]
   end subroutine collect_mqc_fmo_pairs

   subroutine test_level_two(error)
      !! Water trimer, FMO2 with the exact field
      type(error_type), allocatable, intent(out) :: error

      type(fmo_result_t) :: res
      real(dp) :: total, expected_r

      call water_trimer_run(2, res, error)
      if (allocated(error)) return

      call check(error, allocated(res%pairs), "no pairs were exported")
      if (allocated(error)) return
      call check(error, size(res%pairs), 3, "a trimer has three pairs")
      if (allocated(error)) return

      write (*, "(a,es12.4)") "   pairs - pair_sum: ", sum(res%pairs%energy) - res%pair_sum
      call check(error, sum(res%pairs%energy), res%pair_sum, thr=SUM_TOL, &
                 message="the pairs do not add up to the pair sum at level two")
      if (allocated(error)) return

      ! Enumeration order, lower fragment first, numbered from one.
      call check(error, all(res%pairs%i == [1, 1, 2]) .and. &
                 all(res%pairs%j == [2, 3, 3]), "pairs are not (1,2), (1,3), (2,3)")
      if (allocated(error)) return
      call check(error,.not. any(res%pairs%connected), "a water pair was flagged "// &
                 "as covalently connected")
      if (allocated(error)) return

      ! Closest approach of waters 1 and 2 is a hydrogen of the first to the
      ! oxygen of the second: 0.7572 across, 2.9 - 0.5865 up.
      expected_r = sqrt(0.7572_dp**2 + (2.9_dp - 0.5865_dp)**2)
      call check(error, res%pairs(1)%distance, expected_r, thr=1.0e-10_dp, &
                 message="pair (1,2) distance is not the closest H...O approach")
      if (allocated(error)) return

      ! The total is formed exactly as before the pairs were kept.
      total = res%monomer_sum + res%pair_sum
      call check(error, res%energy == total, "the total is no longer "// &
                 "monomer_sum + pair_sum exactly")
      if (allocated(error)) return
      write (*, "(a,f22.15)") "   FMO2 total: ", res%energy

      call check(error, sum(res%pairs%response), res%response_sum, thr=SUM_TOL, &
                 message="the pair responses do not add up to response_sum at level two")
   end subroutine test_level_two

   subroutine test_level_three(error)
      !! The same trimer at FMO3, where the pair map is not the whole expansion
      type(error_type), allocatable, intent(out) :: error

      type(fmo_result_t) :: res
      real(dp) :: shortfall

      call water_trimer_run(3, res, error)
      if (allocated(error)) return

      call check(error, size(res%pairs), 3, "level three must still export "// &
                 "exactly the three pairs, and nothing larger")
      if (allocated(error)) return
      call check(error, size(res%level_sum), 3, "level_sum is not one slot per level")
      if (allocated(error)) return

      shortfall = res%pair_sum - sum(res%pairs%energy)
      write (*, "(a,es14.6,a,es14.6,a,es10.2)") "   shortfall: ", shortfall, &
         "  three-body sum: ", res%level_sum(3), "  difference: ", &
         shortfall - res%level_sum(3)

      ! Not hidden: there is a three-body term and the pairs do not hold it.
      call check(error, abs(res%level_sum(3)) > 1.0e-7_dp, "the three-body term "// &
                 "vanished, so this case no longer shows the shortfall")
      if (allocated(error)) return
      call check(error, shortfall, res%level_sum(3), thr=SUM_TOL, &
                 message="the pairs fall short of pair_sum by something other than "// &
                 "the three-body sum")
      if (allocated(error)) return
      call check(error, sum(res%pairs%energy), res%level_sum(2), thr=SUM_TOL, &
                 message="the pairs do not add up to the two-body level sum")
      if (allocated(error)) return
      call check(error, res%level_sum(1), res%monomer_sum, thr=0.0_dp, &
                 message="level_sum(1) is not the monomer sum")
   end subroutine test_level_three

   subroutine test_peptide(error)
      !! Glycine tripeptide cut at both C-alpha--C(=O) bonds, plus one water
      !!
      !! The geometry is `sample_inputs/gly3_water_pair.xyz`. Fragments 1-2 and
      !! 2-3 share a detached bond and must be flagged; 1-3 and every row with
      !! fragment 4, the water, must not. STO-3G, point charges, FMO2.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(27), owner(27), k
      character(len=2) :: sym(27)
      real(dp) :: xyz(3, 27)
      logical :: joined, seen(3)

      call gly3_water(z, sym, xyz)
      ! The deck's fragments, 0-based there, as owners here.
      owner = 0
      owner([0, 1, 4, 5, 6, 7] + 1) = 1
      owner([2, 3, 8, 9, 12, 13, 14] + 1) = 2
      owner([10, 11, 15, 16, 17, 18, 19, 20, 21, 22, 23] + 1) = 3
      owner([24, 25, 26] + 1) = 4

      opts%basis = "sto-3g"
      opts%esp = "ptc"
      opts%expansion = "fmo"
      opts%bond_breaking = "afo"
      opts%level = 2

      call run_fmo2(z, sym, xyz, owner, opts, res, err)
      call check(error,.not. err%has_error(), "the peptide run failed")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call check(error, size(res%pairs), 6, "four fragments have six pairs")
      if (allocated(error)) return

      seen = .false.
      do k = 1, size(res%pairs)
         joined = (res%pairs(k)%i == 1 .and. res%pairs(k)%j == 2) .or. &
                  (res%pairs(k)%i == 2 .and. res%pairs(k)%j == 3)
         write (*, "(3x,i0,'-',i0,f9.3,es20.10,l3)") res%pairs(k)%i, res%pairs(k)%j, &
            res%pairs(k)%distance, res%pairs(k)%energy, res%pairs(k)%connected
         call check(error, res%pairs(k)%connected .eqv. joined, "pair "// &
                    char(48 + res%pairs(k)%i)//"-"//char(48 + res%pairs(k)%j)// &
                    " has the wrong connected flag")
         if (allocated(error)) return
         if (res%pairs(k)%j == 4) seen(res%pairs(k)%i) = .true.
      end do
      call check(error, all(seen), "a ligand-residue row is missing")
      if (allocated(error)) return

      ! The water is hydrogen-bonded to the carbonyl fragment 3 holds, at 1.893
      ! Angstrom; its row has to say so and read as an interaction, not a bond.
      do k = 1, size(res%pairs)
         if (res%pairs(k)%i /= 3 .or. res%pairs(k)%j /= 4) cycle
         call check(error, res%pairs(k)%distance, 1.893_dp, thr=1.0e-3_dp, &
                    message="the 3-4 row is not the 1.893 A hydrogen bond")
         if (allocated(error)) return
         call check(error, abs(res%pairs(k)%energy) < 0.05_dp, "the 3-4 term is "// &
                    "not the size of a hydrogen bond")
         if (allocated(error)) return
      end do

      write (*, "(a,es12.4)") "   pairs - pair_sum: ", sum(res%pairs%energy) - res%pair_sum
      call check(error, sum(res%pairs%energy), res%pair_sum, thr=SUM_TOL, &
                 message="the pairs do not add up to the pair sum at level two")
   end subroutine test_peptide

   subroutine water_trimer_run(level, res, error)
      !! `sample_inputs/w3.xyz`, three waters stacked 2.9 A apart, STO-3G
      integer, intent(in) :: level
      type(fmo_result_t), intent(out) :: res
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(fmo_options_t) :: opts
      real(dp) :: ang(3, 9)
      integer :: z(9)
      character(len=2) :: sym(9)

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

      opts%basis = "sto-3g"
      opts%level = level
      call run_fmo2(z, sym, to_bohr(ang), [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "the water trimer run failed")
      if (allocated(error)) write (*, *) "   message: ", trim(err%get_message())
   end subroutine water_trimer_run

   subroutine gly3_water(z, sym, xyz)
      !! `sample_inputs/gly3_water_pair.xyz`, in Bohr
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

end module test_mqc_fmo_pairs

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fmo_pairs, only: collect_mqc_fmo_pairs
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fmo_pairs", collect_mqc_fmo_pairs)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
