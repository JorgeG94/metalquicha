!! FMO2's interaction energy of one fragment: every monomer, then only its pairs
module test_mqc_fmo_interaction
   !! `fmo_options_t%reference_fragment` keeps the monomer loop whole -- each
   !! monomer is solved in the field of all the others, so none can be left
   !! out -- and then solves only the pairs that hold the reference. A pair's
   !! term needs only its own SCF and its two monomers, so every row the
   !! reduced run reports must be the row the full FMO2 run reports for the
   !! same pair. That is what is checked here, for every choice of reference:
   !!
   !! 1. the water trimer `sample_inputs/w3.xyz`, HF/6-31G, PIEDA on: pair
   !!    energy, response and Ees/Eex/Ect+mix against the full run's rows;
   !! 2. the glycine tripeptide with a water, `sample_inputs/gly3_water.xyz`,
   !!    cut at both C-alpha--C(=O) bonds with adjusted frozen orbitals under
   !!    the exact field, HF/STO-3G, PIEDA on, `resdim` 2.0: two of the pairs
   !!    are connected (undecomposed, carrying the bond) and one is separated
   !!    (electrostatics, no pair SCF), and both kinds must match too;
   !! 3. the work: the same outer iterations, so the same monomer SCFs, and
   !!    `n - 1` pairs holding the reference with at most `n - 1` pair SCFs;
   !!    no total energy;
   !! 4. the refusals: above level two, under the embedded MBE expansion, and
   !!    a reference out of range.
   !!
   !! The deck path -- the driver, the JSON `interaction_energy` block and the
   !! refusals of `check_interaction_energy_support` -- is checked in
   !! `test_mqc_config_roundtrip`, `test_mqc_config_adapter` and the
   !! validation deck `interaction_energy/ie_fmo2_w3_pieda.json`.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: to_bohr
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   implicit none
   private

   public :: collect_mqc_fmo_interaction

   integer, parameter :: N_W3 = 9
      !! Atoms in the water trimer
   integer, parameter :: N_GLY = 27
      !! Atoms in the glycine tripeptide with a water
   real(dp), parameter :: ROW_TOL = 1.0e-9_dp
      !! A reduced run's row against the full run's, Hartree. The two solve
      !! the same SCFs from the same start, so what separates them is the
      !! order a threaded Fock build merges in -- around 1e-11 measured here.
   real(dp), parameter :: CONNECTED_RESPONSE_TOL = 1.0e-7_dp
      !! The response of a pair joined by a detached bond, which moves by
      !! 1.7e-8 Hartree between two identical threaded full runs of the
      !! glycine case while the pair's energy agrees to 1e-11

contains

   subroutine collect_mqc_fmo_interaction(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("water_trimer_reference_pairs_match_the_full_run", &
                               test_water_trimer), &
                  new_unittest("cut_glycine_water_reference_pairs_match_the_full_run", &
                               test_gly3_water), &
                  new_unittest("a_reference_is_refused_outside_fmo2", test_refusals) &
                  ]
   end subroutine collect_mqc_fmo_interaction

   subroutine test_water_trimer(error)
      !! Three stacked waters, one fragment each, HF/6-31G, PIEDA on
      type(error_type), allocatable, intent(out) :: error

      integer :: z(N_W3)
      character(len=2) :: sym(N_W3)
      real(dp) :: xyz(3, N_W3)
      type(fmo_options_t) :: opts

      call water_trimer(z, sym, xyz)
      call tight_options(opts, "6-31g")
      call compare_with_full_run(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], 3, opts, &
                                 "w3", error)
   end subroutine test_water_trimer

   subroutine test_gly3_water(error)
      !! The tripeptide cut at both C-alpha--C(=O) bonds, and the water
      !!
      !! Fragments as `mqc_docs/source/fmo.rst` cuts a backbone, the water
      !! fourth: atoms 0-1, 4-7 | 2-3, 8-9, 12-14 | 10-11, 15-23 | 24-26,
      !! 0-based. The peptide bonds stay whole inside a fragment, which is
      !! what makes the cut a single bond the frozen orbital can detach.
      type(error_type), allocatable, intent(out) :: error

      integer :: z(N_GLY)
      character(len=2) :: sym(N_GLY)
      real(dp) :: xyz(3, N_GLY)
      type(fmo_options_t) :: opts

      call gly3_water(z, sym, xyz)
      call tight_options(opts, "sto-3g")
      opts%bond_breaking = "afo"
      opts%afo_localization = "er"
      opts%resppc = 2.0_dp
      opts%resdim = 2.0_dp
      call compare_with_full_run(z, sym, xyz, &
                                 [1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 3, 3, 2, 2, 2, &
                                  3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4], 4, opts, &
                                 "gly3_water", error)
   end subroutine test_gly3_water

   subroutine compare_with_full_run(z, sym, xyz, owner, n_frag, base, label, error)
      !! The full FMO2 run, then one reduced run per reference, row by row
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: owner(:)
      integer, intent(in) :: n_frag
      type(fmo_options_t), intent(in) :: base
      character(len=*), intent(in) :: label

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: full, red
      type(error_t) :: err
      real(dp) :: worst_row, worst_monomer, expected_sum
      integer :: ref, p, k, n_pair_scf

      opts = base
      call run_fmo2(z, sym, xyz, owner, opts, full, err)
      call check(error,.not. err%has_error(), label//": the full FMO2 run failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, size(full%pairs), n_frag*(n_frag - 1)/2, &
                 label//": the full run should report every pair")
      if (allocated(error)) return
      write (*, "(a,a,a,i0,a,i0,a,i0,a)") "   ", label, ": full run, ", n_frag, &
         " monomers, ", full%outer_iterations, " outer iterations, ", &
         count(.not. full%pairs%separated), " pair SCFs"

      worst_row = 0.0_dp
      worst_monomer = 0.0_dp
      do ref = 1, n_frag
         opts = base
         opts%reference_fragment = ref
         call run_fmo2(z, sym, xyz, owner, opts, red, err)
         call check(error,.not. err%has_error(), label//": the reduced run failed: "// &
                    err%get_message())
         if (allocated(error)) return

         ! The monomer phase is the full run's, untouched.
         call check(error, red%converged .and. red%outer_iterations == full%outer_iterations, &
                    label//": the monomer loop should run exactly as in the full run")
         if (allocated(error)) return
         worst_monomer = max(worst_monomer, &
                             maxval(abs(red%monomer_energy - full%monomer_energy)))

         ! Only the reference's pairs, and no total.
         call check(error, size(red%pairs), n_frag - 1, &
                    label//": a reduced run should report the n - 1 pairs of its reference")
         if (allocated(error)) return
         call check(error, all(red%pairs%i == ref .or. red%pairs%j == ref), &
                    label//": every pair of a reduced run should hold the reference")
         if (allocated(error)) return
         n_pair_scf = count(.not. red%pairs%separated)
         call check(error, n_pair_scf <= n_frag - 1, &
                    label//": no more than n - 1 pair SCFs")
         if (allocated(error)) return
         call check(error, red%energy == 0.0_dp, &
                    label//": a reduced run must not report a total energy")
         if (allocated(error)) return
         write (*, "(a,a,a,i0,a,i0,a,i0,a,i0,a)") "   ", label, ": reference ", ref, ", ", &
            n_frag, " monomers x ", red%outer_iterations, " passes, ", n_pair_scf, &
            " pair SCFs"

         expected_sum = 0.0_dp
         do p = 1, size(red%pairs)
            k = findloc(full%pairs%i == red%pairs(p)%i .and. &
                        full%pairs%j == red%pairs(p)%j, .true., dim=1)
            call check(error, k > 0, label//": a reduced pair is missing from the full run")
            if (allocated(error)) return
            associate (a => red%pairs(p), b => full%pairs(k))
               call check(error, (a%connected .eqv. b%connected) .and. &
                          (a%separated .eqv. b%separated) .and. (a%pieda .eqv. b%pieda), &
                          label//": a pair's flags differ from the full run's")
               if (allocated(error)) return
               call check(error, a%distance == b%distance, &
                          label//": a pair's distance differs from the full run's")
               if (allocated(error)) return
               write (*, "(a,a,a,i0,a,i0,a,l1,a,l1,a,5es10.2)") "   ", label, ": pair ", &
                  a%i, "-", a%j, " connected ", a%connected, " separated ", a%separated, &
                  " d(E, resp, ees, eex, ect) =", a%energy - b%energy, &
                  a%response - b%response, a%ees - b%ees, a%eex - b%eex, &
                  a%ect_mix - b%ect_mix
               worst_row = max(worst_row, abs(a%energy - b%energy), abs(a%ees - b%ees), &
                               abs(a%eex - b%eex), abs(a%ect_mix - b%ect_mix))
               ! A connected pair's response is not reproducible run to run --
               ! see the TODO beside `response(t)` in `calculate_polymers` --
               ! though its energy is; held to the size of that wobble.
               if (a%connected) then
                  call check(error, abs(a%response - b%response) < CONNECTED_RESPONSE_TOL, &
                             label//": a connected pair's response differs from the "// &
                             "full run's beyond its run-to-run wobble")
                  if (allocated(error)) return
               else
                  worst_row = max(worst_row, abs(a%response - b%response))
               end if
               expected_sum = expected_sum + b%energy
            end associate
         end do
         call check(error, abs(red%pair_sum - expected_sum) < ROW_TOL*n_frag, &
                    label//": the reduced pair sum is not the full run's sum over "// &
                    "the reference's pairs")
         if (allocated(error)) return
         call check(error, abs(red%level_sum(2) - red%pair_sum) < epsilon(1.0_dp)* &
                    max(1.0_dp, abs(red%pair_sum)), &
                    label//": level_sum(2) should be the pair sum")
         if (allocated(error)) return
      end do

      write (*, "(a,a,a,es10.3,a,es10.3)") "   ", label, &
         ": worst row difference ", worst_row, " Eh, worst monomer ", worst_monomer
      call check(error, worst_row < ROW_TOL, &
                 label//": a reduced row differs from the full run's beyond SCF noise")
      if (allocated(error)) return
      call check(error, worst_monomer < ROW_TOL, &
                 label//": a monomer energy differs from the full run's beyond SCF noise")
      if (allocated(error)) return

      ! PIEDA ran on the reduced pairs it would have run on in the full one.
      call check(error, any(full%pairs%pieda), label//": PIEDA decomposed nothing")
   end subroutine compare_with_full_run

   subroutine test_refusals(error)
      !! Above level two, under the embedded MBE expansion, and out of range
      type(error_type), allocatable, intent(out) :: error

      integer :: z(N_W3)
      character(len=2) :: sym(N_W3)
      real(dp) :: xyz(3, N_W3)
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(error_t) :: err

      call water_trimer(z, sym, xyz)

      call tight_options(opts, "sto-3g")
      opts%level = 3
      opts%reference_fragment = 1
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error, err%has_error(), "a reference at level 3 must be refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "pairs of FMO2") > 0, &
                 "level 3 was refused for another reason: "//err%get_message())
      if (allocated(error)) return

      call err%clear()
      call tight_options(opts, "sto-3g")
      opts%expansion = "mbe"
      opts%esp = "ptc"
      opts%reference_fragment = 1
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error, err%has_error(), "a reference under EE-MBE must be refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "pairs of FMO2") > 0, &
                 "EE-MBE was refused for another reason: "//err%get_message())
      if (allocated(error)) return

      call err%clear()
      call tight_options(opts, "sto-3g")
      opts%reference_fragment = 4
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error, err%has_error(), "a reference past the last fragment must be refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "out of range") > 0, &
                 "an out-of-range reference was refused for another reason: "// &
                 err%get_message())
   end subroutine test_refusals

   subroutine tight_options(opts, basis)
      !! FMO2, exact field, PIEDA on, tight enough that rows agree to 1e-10
      character(len=*), intent(in) :: basis
      type(fmo_options_t), intent(out) :: opts

      opts%basis = basis
      opts%esp = "exact"
      opts%expansion = "fmo"
      opts%level = 2
      opts%resdim = 0.0_dp
      opts%outer_tol = 1.0e-10_dp
      opts%max_outer = 100
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      opts%pieda = .true.
   end subroutine tight_options

   subroutine water_trimer(z, sym, xyz)
      !! `validation/inputs/sample_inputs/w3.xyz`: three waters stacked 2.9 A
      !! apart along z; fragment `k` is atoms `3k-2..3k`
      integer, intent(out) :: z(N_W3)
      character(len=2), intent(out) :: sym(N_W3)
      real(dp), intent(out) :: xyz(3, N_W3)
      real(dp) :: ang(3, N_W3)

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

   subroutine gly3_water(z, sym, xyz)
      !! `validation/inputs/sample_inputs/gly3_water.xyz`: the glycine
      !! tripeptide with a water donating a hydrogen bond to the middle
      !! carbonyl oxygen, Angstrom
      integer, intent(out) :: z(N_GLY)
      character(len=2), intent(out) :: sym(N_GLY)
      real(dp), intent(out) :: xyz(3, N_GLY)
      real(dp) :: ang(3, N_GLY)
      integer :: i

      z = [7, 6, 6, 8, 1, 1, 1, 1, 7, 6, 6, 8, 1, 1, 1, 7, 6, 6, 8, 1, 1, 1, 8, 1, 8, 1, 1]
      do i = 1, N_GLY
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
                    5.9540000000_dp, -3.0660000000_dp, 1.8890000000_dp, &
                    5.9540000000_dp, -2.1060000000_dp, 1.8890000000_dp, &
                    6.8834217348_dp, -3.3063648039_dp, 1.8890000000_dp], [3, 27])
      xyz = to_bohr(ang)
   end subroutine gly3_water

end module test_mqc_fmo_interaction

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fmo_interaction, only: collect_mqc_fmo_interaction
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fmo_interaction", collect_mqc_fmo_interaction)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
