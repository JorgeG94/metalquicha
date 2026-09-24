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
                               test_peptide), &
                  new_unittest("separated_water_pairs_match_gamess_es_dimers", &
                               test_separated_waters) &
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

   subroutine test_separated_waters(error)
      !! Twenty waters, FMO2 with separated pairs, against GAMESS
      !!
      !! `sample_inputs/w20_isomer1.xyz`, a water per fragment, RHF/STO-3G in
      !! the exact field. GAMESS 2026 `$FMO NBODY=2 RESDIM=2.0` with its other
      !! defaults, which are `RESPPC=2.0` and `RESPAP=0`, `$SCF CONV=1D-8`: 64
      !! of the 190 pairs separated, and -1499.603170330. Each separated pair's
      !! `E"IJ-E"I-E"J`, printed to eight decimals, is listed below; it has no
      !! response term, so it is the whole of the pair's energy.
      !!
      !! With `RESDIM=0` GAMESS gives -1499.603176697, which is what this code
      !! gave before separated pairs existed, and still gives at `resdim` zero.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N_ES = 64
      integer, parameter :: ES_I(N_ES) = [2, 4, 4, 5, 6, 2, 5, 6, 2, 5, 5, 6, 8, 9, 7, 10, &
                                          5, 6, 4, 7, 8, 9, 10, 12, 4, 6, 8, 9, 11, 13, 4, 6, &
                                          8, 9, 11, 4, 7, 10, 12, 14, 15, 2, 5, 6, 13, 14, 15, 16, &
                                          4, 8, 9, 17, 5, 6, 8, 9, 11, 13, 16, 17, 10, 14, 15, 19]
      integer, parameter :: ES_J(N_ES) = [4, 5, 6, 7, 7, 8, 8, 8, 9, 9, 10, 10, 10, 10, 11, 11, &
                                          12, 12, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 15, 15, &
                                          15, 15, 15, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, &
                                          18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20]
      real(dp), parameter :: ES_E(N_ES) = [-0.00026283_dp, 0.00010071_dp, -0.00012010_dp, 0.00004412_dp, &
                                           0.00014818_dp, 0.00019615_dp, -0.00014373_dp, -0.00000383_dp, &
                                           0.00009161_dp, -0.00004211_dp, -0.00008450_dp, -0.00014765_dp, &
                                           -0.00013899_dp, 0.00018447_dp, -0.00029092_dp, -0.00014250_dp, &
                                           -0.00008697_dp, -0.00027936_dp, -0.00001436_dp, -0.00010456_dp, &
                                           -0.00002696_dp, -0.00011125_dp, 0.00006651_dp, 0.00005598_dp, &
                                           -0.00012878_dp, -0.00015869_dp, 0.00005897_dp, -0.00001511_dp, &
                                           -0.00013229_dp, 0.00051170_dp, -0.00003657_dp, 0.00012355_dp, &
                                           0.00001580_dp, 0.00009799_dp, 0.00007224_dp, 0.00018470_dp, &
                                           -0.00014432_dp, 0.00008857_dp, 0.00028586_dp, 0.00007309_dp, &
                                           -0.00009295_dp, -0.00012413_dp, 0.00001655_dp, 0.00014023_dp, &
                                           0.00005259_dp, -0.00002717_dp, -0.00014763_dp, -0.00013955_dp, &
                                           -0.00000044_dp, 0.00000881_dp, 0.00030791_dp, -0.00031591_dp, &
                                           0.00011565_dp, 0.00009779_dp, 0.00013096_dp, -0.00011885_dp, &
                                           0.00014039_dp, -0.00015564_dp, -0.00014492_dp, -0.00001086_dp, &
                                           0.00026863_dp, 0.00003409_dp, -0.00007594_dp, -0.00024506_dp]
      real(dp), parameter :: GAMESS_RESDIM2 = -1499.603170330_dp
      real(dp), parameter :: GAMESS_RESDIM0 = -1499.603176697_dp
      type(fmo_result_t) :: res
      integer :: p, k
      real(dp) :: worst

      call water20_run(2.0_dp, res, error)
      if (allocated(error)) return
      write (*, *) "   FMO2 w20, resdim 2 - GAMESS =", res%energy - GAMESS_RESDIM2
      call check(error, abs(res%energy - GAMESS_RESDIM2) < 1.0e-7_dp, &
                 "FMO2 with separated pairs does not match GAMESS")
      if (allocated(error)) return
      call check(error, count(res%pairs%separated), N_ES, &
                 "a different set of pairs is separated than GAMESS's")
      if (allocated(error)) return

      worst = 0.0_dp
      do p = 1, N_ES
         k = findloc(res%pairs%i == ES_I(p) .and. res%pairs%j == ES_J(p), .true., dim=1)
         call check(error, k > 0, "a pair GAMESS separates is missing")
         if (allocated(error)) return
         call check(error, res%pairs(k)%separated, "a pair GAMESS separates was solved")
         if (allocated(error)) return
         call check(error, res%pairs(k)%response == 0.0_dp, &
                    "a separated pair carries a response term")
         if (allocated(error)) return
         worst = max(worst, abs(res%pairs(k)%energy - ES_E(p)))
      end do
      write (*, *) "   worst separated pair - GAMESS =", worst
      ! Eight printed decimals, so half a unit of the last is rounding.
      call check(error, worst < 1.0e-8_dp, &
                 "a separated pair's energy does not match GAMESS's ES dimer")
      if (allocated(error)) return

      call water20_run(0.0_dp, res, error)
      if (allocated(error)) return
      write (*, *) "   FMO2 w20, resdim 0 - GAMESS =", res%energy - GAMESS_RESDIM0
      call check(error,.not. any(res%pairs%separated), "resdim 0 separated a pair")
      if (allocated(error)) return
      call check(error, abs(res%energy - GAMESS_RESDIM0) < 1.0e-7_dp, &
                 "FMO2 solving every pair does not match GAMESS at RESDIM=0")
   end subroutine test_separated_waters

   subroutine water20_run(resdim, res, error)
      !! `sample_inputs/w20_isomer1.xyz`, a water per fragment, STO-3G, exact field
      real(dp), intent(in) :: resdim
      type(fmo_result_t), intent(out) :: res
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(fmo_options_t) :: opts
      real(dp) :: ang(3, 60)
      integer :: z(60), owner(60), k
      character(len=2) :: sym(60)

      do k = 1, 60
         owner(k) = (k - 1)/3 + 1
         if (mod(k - 1, 3) == 0) then
            z(k) = 8
            sym(k) = "O "
         else
            z(k) = 1
            sym(k) = "H "
         end if
      end do
      ang = reshape([ &
                    2.19756413_dp, 2.30645657_dp, 2.12398648_dp, &
                    1.30374396_dp, 2.67852998_dp, 2.14887738_dp, &
                    2.77602649_dp, 2.98882961_dp, 2.48710012_dp, &
                    -0.446642250_dp, 3.18157983_dp, -0.553412318_dp, &
                    -0.294790566_dp, 4.09050608_dp, -0.841292858_dp, &
                    0.405002892_dp, 2.72651124_dp, -0.635274589_dp, &
                    -0.371120274_dp, 3.40200686_dp, 2.21183610_dp, &
                    -0.776223004_dp, 2.65618086_dp, 2.67454791_dp, &
                    -0.552393317_dp, 3.32637596_dp, 1.26568747_dp, &
                    -0.344809651_dp, 2.19895148_dp, 6.16216040_dp, &
                    -0.799917400_dp, 1.92929876_dp, 5.35409737_dp, &
                    -0.216623634_dp, 3.15416622_dp, 6.11619949_dp, &
                    2.97632003_dp, 4.44283104_dp, -1.21666169_dp, &
                    2.13466644_dp, 4.88871813_dp, -1.37117922_dp, &
                    2.76599479_dp, 3.50445747_dp, -1.12249064_dp, &
                    4.16477871_dp, 5.80204582_dp, 0.813909113_dp, &
                    3.81304598_dp, 5.27708006_dp, 7.27024004E-02_dp, &
                    4.95197010_dp, 6.23447180_dp, 0.493922502_dp, &
                    -1.12275231_dp, 1.15836358_dp, 3.63779879_dp, &
                    -1.67120850_dp, 0.587195158_dp, 3.08564949_dp, &
                    -0.245801628_dp, 0.747786403_dp, 3.66580105_dp, &
                    0.566819966_dp, 4.92840958_dp, 6.15334225_dp, &
                    1.46272981_dp, 4.54431629_dp, 6.13453388_dp, &
                    0.518503189_dp, 5.47743940_dp, 6.93166018_dp, &
                    2.99093223_dp, 3.72968745_dp, 5.88431597_dp, &
                    3.32716727_dp, 3.93724680_dp, 5.00251436_dp, &
                    2.86334229_dp, 2.77451253_dp, 5.92706394_dp, &
                    0.381184578_dp, -1.53543997_dp, 1.63772929_dp, &
                    -0.556514800_dp, -1.34652114_dp, 1.76359665_dp, &
                    0.843591511_dp, -1.10289204_dp, 2.36837602_dp, &
                    3.40574169_dp, 4.54636431_dp, 3.20606375_dp, &
                    2.53101349_dp, 4.95308876_dp, 3.30422902_dp, &
                    3.81790257_dp, 4.97632933_dp, 2.44633842_dp, &
                    1.49123096_dp, 0.106784083_dp, 3.62084031_dp, &
                    1.85748398_dp, 0.819141865_dp, 3.07536340_dp, &
                    1.79046202_dp, 0.273016363_dp, 4.52292490_dp, &
                    0.409543365_dp, 5.76687193_dp, -1.21077549_dp, &
                    9.09162611E-02_dp, 6.45845699_dp, -1.78449512_dp, &
                    0.723045111_dp, 6.20045519_dp, -0.396752506_dp, &
                    0.965811729_dp, -0.716764033_dp, -0.878285229_dp, &
                    0.826097012_dp, -1.06112885_dp, 2.22627297E-02_dp, &
                    1.20147121_dp, -1.46186376_dp, -1.42444563_dp, &
                    -1.57260728_dp, 0.603498220_dp, -0.564730942_dp, &
                    -0.767514944_dp, 0.135007307_dp, -0.818493128_dp, &
                    -1.36242235_dp, 1.54356337_dp, -0.635614216_dp, &
                    1.48347354_dp, 6.85378027_dp, 1.03128028_dp, &
                    2.42572093_dp, 6.64860678_dp, 1.00056386_dp, &
                    1.17456794_dp, 6.55655336_dp, 1.89760959_dp, &
                    2.09537220_dp, 1.00540781_dp, 6.21120834_dp, &
                    1.20502126_dp, 1.39716053_dp, 6.27762365_dp, &
                    2.25407863_dp, 0.541421950_dp, 7.02924109_dp, &
                    2.02927136_dp, 1.86724663_dp, -0.597557843_dp, &
                    1.77290356_dp, 0.952067673_dp, -0.764892220_dp, &
                    2.20549679_dp, 1.94485366_dp, 0.352207899_dp, &
                    -2.32855034_dp, -0.538468540_dp, 1.77467108_dp, &
                    -3.16010594_dp, -0.993369758_dp, 1.67372715_dp, &
                    -2.13997388_dp, -0.103332303_dp, 0.923533261_dp, &
                    0.792866886_dp, 5.56255627_dp, 3.42352605_dp, &
                    0.296446860_dp, 4.84693289_dp, 2.99636626_dp, &
                    0.622788429_dp, 5.48711681_dp, 4.36984348_dp &
                    ], [3, 60])

      opts%basis = "sto-3g"
      opts%esp = "exact"
      opts%expansion = "fmo"
      opts%resppc = 2.0_dp
      opts%resdim = resdim
      opts%level = 2
      call run_fmo2(z, sym, to_bohr(ang), owner, opts, res, err)
      call check(error,.not. err%has_error(), "the twenty-water run failed")
      if (allocated(error)) write (*, *) "   message: ", trim(err%get_message())
   end subroutine water20_run

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
