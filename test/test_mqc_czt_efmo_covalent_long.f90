!! EFMO across covalent bonds: the far pairs, measured at every separation
module test_mqc_czt_efmo_covalent_long
   !! Split from `test_mqc_czt_efmo_covalent` and labelled `LONG`. Each value
   !! here is a pair computed twice, effective and as a quantum dimer, and each
   !! of those is a full EFMO run with a MAKEFP per monomer. Together they were
   !! most of that file's time, while the routine file keeps every assertion:
   !! its far-water case runs the one separation this file's assertion turns
   !! on.
   !!
   !! What is here:
   !!
   !!   1. A water beyond the cutoff, against a cut fragment and against the
   !!      same molecule uncut, at 3.5, 4.5 and 6.0 Angstrom. These are the
   !!      far-pair numbers `efmo.rst` quotes.
   !!   2. Butane in four, whose end pair is the one pair that is effective;
   !!      its effective and quantum values are reported.
   !!
   !! 6-31G throughout, the basis `efmo.rst` quotes; the routine file runs
   !! STO-3G.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_czt_efmo, only: efmo_options_t, efmo_result_t, run_efmo, efmo_pair_contribution
   use mqc_physical_constants, only: ANGSTROM_TO_BOHR
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_czt_efmo_covalent_long_tests

   real(dp), parameter :: ANG = ANGSTROM_TO_BOHR

   character(len=*), parameter :: BASIS = "6-31g"

contains

   subroutine collect_mqc_czt_efmo_covalent_long_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("efmo_afo_far_water_pairs_at_every_separation", &
                               test_far_water_sweep), &
                  new_unittest("efmo_afo_butane_end_pair_effective_and_quantum", &
                               test_butane_end_pair) &
                  ]
   end subroutine collect_mqc_czt_efmo_covalent_long_tests

   subroutine test_far_water_sweep(error)
      !! A water beyond the cutoff, against a cut fragment and against the same
      !! molecule uncut, at three separations
      !!
      !! Each far pair's four effective-fragment terms are compared with the
      !! same pair solved as a quantum dimer, `E_IJ^0 - E_I^0 - E_J^0 -
      !! E_IJ^pol`, which is what the pair contributes when it is near. The
      !! uncut molecule is the baseline: the effective-fragment approximation
      !! has an error of its own there, and what a cut must not do is make it
      !! materially worse. The routine file asserts the same at 3.5 Angstrom,
      !! where both errors are largest.
      type(error_type), allocatable, intent(out) :: error
      integer :: z(14), k
      character(len=2) :: sym(14)
      real(dp) :: xyz(3, 14), gap(3), efp_w, qm_w, efp_c, qm_c, worst_whole, worst_cut
      real(dp), parameter :: SEPARATION(3) = [3.5_dp, 4.5_dp, 6.0_dp]
      type(error_t) :: err

      worst_whole = 0.0_dp
      worst_cut = 0.0_dp
      do k = 1, size(SEPARATION)
         call propane_water(SEPARATION(k), z, sym, xyz)
         ! Uncut: propane is fragment 1, the water fragment 2.
         call far_and_near(z, sym, xyz, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2], &
                           [1, 2], "none", efp_w, qm_w, err)
         call check(error,.not. err%has_error(), "uncut run failed: "// &
                    err%get_full_trace())
         if (allocated(error)) return
         ! Cut: the methyl facing the water is fragment 1, the ethyl 2.
         call far_and_near(z, sym, xyz, [1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3], &
                           [1, 3], "afo", efp_c, qm_c, err)
         call check(error,.not. err%has_error(), "cut run failed: "// &
                    err%get_full_trace())
         if (allocated(error)) return
         gap = [efp_w - qm_w, efp_c - qm_c, 0.0_dp]
         write (*, "(a,f5.2,a,4es13.4,a,2es11.3)") "   water at ", SEPARATION(k), &
            " A: whole efp/qm, cut efp/qm =", efp_w, qm_w, efp_c, qm_c, &
            "  errors whole, cut:", gap(1), gap(2)
         worst_whole = max(worst_whole, abs(gap(1)))
         worst_cut = max(worst_cut, abs(gap(2)))
      end do
      write (*, *) "   worst |efp - qm|: whole", worst_whole, " cut", worst_cut
      call check(error, worst_cut < 3.0_dp*worst_whole + 2.0e-4_dp, &
                 "a far pair with a cut fragment is much worse than the same pair uncut")
   end subroutine test_far_water_sweep

   subroutine test_butane_end_pair(error)
      !! Butane in four pieces, its end pair effective and quantum
      !!
      !! The first and last fragments share no centre and no bond, so at a
      !! cutoff far inside every separation that pair alone is effective. Two
      !! cuts apart and at `R_IJ` about 0.85 it is a poor effective pair --
      !! almost all electrostatics -- and the two values are reported so the
      !! number `efmo.rst` quotes can be reproduced.
      type(error_type), allocatable, intent(out) :: error
      integer :: z(14)
      character(len=2) :: sym(14)
      real(dp) :: xyz(3, 14), efp_value, qm_value
      type(error_t) :: err

      call butane(z, sym, xyz)
      call far_and_near(z, sym, xyz, [1, 2, 3, 4, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4], [1, 4], &
                        "afo", efp_value, qm_value, err)
      call check(error,.not. err%has_error(), "the end pair failed: "//err%get_full_trace())
      if (allocated(error)) return
      write (*, *) "   butane end pair (1,4): efp", efp_value, " qm", qm_value
      call check(error, efp_value /= 0.0_dp .and. qm_value /= 0.0_dp, &
                 "the end pair was not found in both runs")
   end subroutine test_butane_end_pair

   subroutine far_and_near(z, sym, xyz, owner, pair, bond_breaking, efp_value, qm_value, err)
      !! One pair's contribution computed effective and computed quantum
      integer, intent(in) :: z(:), owner(:), pair(2)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: xyz(:, :)
      character(len=*), intent(in) :: bond_breaking
      real(dp), intent(out) :: efp_value, qm_value
      type(error_t), intent(inout) :: err

      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      integer, allocatable :: charges(:)
      integer :: k

      efp_value = 0.0_dp
      qm_value = 0.0_dp
      allocate (charges(maxval(owner)), source=0)
      call settings(opts)
      opts%bond_breaking = bond_breaking
      opts%induction_damping = 0.1_dp
      opts%rcut = 0.1_dp
      call run_efmo(z, sym, xyz, owner, charges, opts, res, err)
      if (err%has_error()) return
      do k = 1, size(res%pairs)
         if (res%pairs(k)%i == pair(1) .and. res%pairs(k)%j == pair(2)) then
            efp_value = efmo_pair_contribution(res, k)
         end if
      end do
      opts%rcut = 1.0e6_dp
      call run_efmo(z, sym, xyz, owner, charges, opts, res, err)
      if (err%has_error()) return
      do k = 1, size(res%pairs)
         if (res%pairs(k)%i == pair(1) .and. res%pairs(k)%j == pair(2)) then
            qm_value = efmo_pair_contribution(res, k)
         end if
      end do
   end subroutine far_and_near

   subroutine settings(opts)
      !! The routine file's settings
      type(efmo_options_t), intent(out) :: opts

      opts%basis = BASIS
      opts%bond_breaking = "afo"
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-10_dp
      opts%scf_density_tol = 1.0e-8_dp
      opts%scf_grad_tol = 1.0e-8_dp
      opts%scf%grad_tol = 1.0e-8_dp
   end subroutine settings

   subroutine propane_water(distance, z, sym, xyz)
      !! Propane, and a water on the far side of its first methyl
      !!
      !! The oxygen sits `distance` Angstrom beyond the first carbon along the
      !! C-C axis, its hydrogens pointing away.
      real(dp), intent(in) :: distance
      integer, intent(out) :: z(14)
      character(len=2), intent(out) :: sym(14)
      real(dp), intent(out) :: xyz(3, 14)

      integer :: zp(11)
      character(len=2) :: sp(11)
      real(dp) :: xp(3, 11), x0

      call propane(zp, sp, xp)
      z(1:11) = zp
      sym(1:11) = sp
      xyz(:, 1:11) = xp
      z(12:14) = [8, 1, 1]
      sym(12:14) = ["O ", "H ", "H "]
      x0 = 1.5260_dp + distance
      xyz(:, 12) = [x0, 0.0_dp, 0.0_dp]*ANG
      xyz(:, 13) = [x0 + 0.5686_dp, 0.7725_dp, 0.0_dp]*ANG
      xyz(:, 14) = [x0 + 0.5686_dp, -0.7725_dp, 0.0_dp]*ANG
   end subroutine propane_water

   subroutine propane(z, sym, xyz)
      !! Idealised propane, carbons first: C1 end, C2 middle, C3 end
      integer, intent(out) :: z(11)
      character(len=2), intent(out) :: sym(11)
      real(dp), intent(out) :: xyz(3, 11)

      z = [6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1]
      sym = ["C ", "C ", "C ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", "H "]
      xyz = reshape([1.5260_dp, 0.0000_dp, 0.0000_dp, &
                     0.0000_dp, 0.0000_dp, 0.0000_dp, &
                     -0.5716_dp, 1.4149_dp, 0.0000_dp, &
                     2.1553_dp, -0.8900_dp, 0.0000_dp, &
                     2.1553_dp, 0.4450_dp, -0.7707_dp, &
                     2.1553_dp, 0.4450_dp, 0.7707_dp, &
                     -0.3519_dp, -0.5217_dp, 0.8900_dp, &
                     -0.3519_dp, -0.5217_dp, -0.8900_dp, &
                     0.0178_dp, 2.3318_dp, 0.0000_dp, &
                     -1.2200_dp, 1.8317_dp, -0.7707_dp, &
                     -1.2200_dp, 1.8317_dp, 0.7707_dp], [3, 11])*ANG
   end subroutine propane

   subroutine butane(z, sym, xyz)
      !! Anti butane, C-C 1.53, CCC 112 degrees, C-H 1.09; carbons first
      integer, intent(out) :: z(14)
      character(len=2), intent(out) :: sym(14)
      real(dp), intent(out) :: xyz(3, 14)

      z = [6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
      sym = ["C ", "C ", "C ", "C ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", &
             "H ", "H "]
      xyz = reshape([0.0000_dp, 0.0000_dp, 0.0000_dp, &
                     1.2684_dp, 0.8556_dp, 0.0000_dp, &
                     2.5369_dp, 0.0000_dp, 0.0000_dp, &
                     3.8053_dp, 0.8556_dp, 0.0000_dp, &
                     0.2729_dp, -1.0553_dp, 0.0000_dp, &
                     -0.5889_dp, 0.2224_dp, -0.8898_dp, &
                     -0.5889_dp, 0.2224_dp, 0.8898_dp, &
                     1.2684_dp, 1.4847_dp, 0.8901_dp, &
                     1.2684_dp, 1.4847_dp, -0.8901_dp, &
                     2.5369_dp, -0.6291_dp, -0.8901_dp, &
                     2.5369_dp, -0.6291_dp, 0.8901_dp, &
                     3.5324_dp, 1.9108_dp, 0.0000_dp, &
                     4.3942_dp, 0.6331_dp, -0.8898_dp, &
                     4.3942_dp, 0.6331_dp, 0.8898_dp], [3, 14])*ANG
   end subroutine butane

end module test_mqc_czt_efmo_covalent_long

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_efmo_covalent_long, only: collect_mqc_czt_efmo_covalent_long_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_efmo_covalent_long", &
                               collect_mqc_czt_efmo_covalent_long_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
