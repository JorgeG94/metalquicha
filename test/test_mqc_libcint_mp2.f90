module test_mqc_libcint_mp2
   !! Pins conventional MP2 against what a reference energy cannot check.
   !!
   !! The validation manifest already compares total energies to PySCF across
   !! bases, molecules and frozen-core counts, so this suite is for the
   !! structural properties:
   !!
   !!   * the four-index transform has to preserve the integrals' permutational
   !!     symmetry, (ia|jb) == (jb|ia). An index transposed in one of the four
   !!     quarter-steps breaks that while still producing a plausible energy;
   !!   * the correlation energy has to be negative for a stable closed-shell
   !!     reference, and the two spin components separately so;
   !!   * freezing orbitals has to raise it, monotonically, and by roughly what
   !!     a core contributes rather than by an arbitrary amount;
   !!   * the spin scaling has to be applied by `total` and not folded into the
   !!     components, so a scaled result still reports what it scaled.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_result_types, only: mp2_energy_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_mp2, only: mp2_result_t, run_libcint_mp2
   implicit none
   private
   public :: collect_mqc_libcint_mp2_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: E_TOL = 1.0e-11_dp
   real(dp), parameter :: D_TOL = 1.0e-9_dp

contains

   subroutine collect_mqc_libcint_mp2_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("mp2_correlation_is_negative", test_negative), &
                  new_unittest("mp2_components_sum_to_the_correlation", test_sum), &
                  new_unittest("mp2_freezing_core_raises_the_energy", test_frozen), &
                  new_unittest("mp2_rejects_freezing_everything", test_bad_frozen), &
                  new_unittest("mp2_scaling_is_applied_by_total_only", test_scaling), &
                  new_unittest("mp2_scs_matches_its_published_factors", test_scs) &
                  ]
   end subroutine collect_mqc_libcint_mp2_tests

   subroutine water_mp2(basis, frozen, mp2, err)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: frozen
      type(mp2_result_t), intent(out) :: mp2
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.757_dp*ANG, 0.587_dp*ANG, &
                   0.0_dp, -0.757_dp*ANG, 0.587_dp*ANG], [3, 3])
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 200, E_TOL, D_TOL, .false., scf, err)
      if (err%has_error()) return
      call run_libcint_mp2(mol, scf%orbitals, scf%orbital_energies, 5, scf%energy, &
                           mp2, err, n_frozen=frozen)
      call mol%destroy()
   end subroutine water_mp2

   subroutine test_negative(error)
      !! Both spin components are negative for a stable reference
      !!
      !! Every denominator is e_occ - e_virt, which is negative below the gap,
      !! and both numerators are squares or near-squares. A positive component
      !! means either the reference is not a minimum or a sign is wrong, and
      !! the two are worth being able to tell apart.
      type(error_type), allocatable, intent(out) :: error
      type(mp2_result_t) :: mp2
      type(error_t) :: err

      call water_mp2("cc-pvdz", 0, mp2, err)
      call check(error,.not. err%has_error(), "water MP2 must run")
      if (allocated(error)) return
      call check(error, mp2%opposite_spin < 0.0_dp, "E_OS must be negative")
      if (allocated(error)) return
      call check(error, mp2%same_spin < 0.0_dp, "E_SS must be negative")
      if (allocated(error)) return
      call check(error, mp2%total < mp2%scf_energy, &
                 "correlation must lower the reference energy")
   end subroutine test_negative

   subroutine test_sum(error)
      !! The components add up, and opposite-spin dominates
      !!
      !! The ratio is the useful half. E_OS is around three times E_SS for a
      !! closed-shell molecule, and a transform that mixed the two -- swapping
      !! (ib|ja) for (ia|jb) in the same-spin term, say -- would still sum to
      !! something negative and still look like a correlation energy.
      type(error_type), allocatable, intent(out) :: error
      type(mp2_result_t) :: mp2
      type(error_t) :: err

      call water_mp2("cc-pvdz", 0, mp2, err)
      call check(error, abs(mp2%correlation - (mp2%same_spin + mp2%opposite_spin)) &
                 < 1.0e-14_dp, "the components must sum to the correlation energy")
      if (allocated(error)) return
      call check(error, abs(mp2%total - (mp2%scf_energy + mp2%correlation)) < 1.0e-14_dp, &
                 "the total must be the reference plus the correlation")
      if (allocated(error)) return
      call check(error, mp2%opposite_spin < mp2%same_spin, &
                 "opposite-spin correlation must be the larger of the two")
      if (allocated(error)) return
      call check(error, mp2%opposite_spin/mp2%same_spin > 2.0_dp .and. &
                 mp2%opposite_spin/mp2%same_spin < 4.0_dp, &
                 "the OS/SS ratio should be around three for a closed shell")
   end subroutine test_sum

   subroutine test_frozen(error)
      !! Correlating fewer orbitals recovers less correlation
      type(error_type), allocatable, intent(out) :: error
      type(mp2_result_t) :: all_e, frozen
      type(error_t) :: err

      call water_mp2("cc-pvdz", 0, all_e, err)
      call water_mp2("cc-pvdz", 1, frozen, err)
      call check(error,.not. err%has_error(), "both must run")
      if (allocated(error)) return

      call check(error, frozen%n_occupied, all_e%n_occupied - 1)
      if (allocated(error)) return
      call check(error, frozen%n_virtual, all_e%n_virtual, &
                 "freezing the core must not change the virtual space")
      if (allocated(error)) return
      call check(error, frozen%correlation > all_e%correlation, &
                 "a frozen core must recover less correlation")
      if (allocated(error)) return
      ! Oxygen's 1s is worth a couple of millihartree here, not tens.
      call check(error, abs(frozen%correlation - all_e%correlation) < 0.01_dp, &
                 "freezing one core orbital should cost only a few millihartree")
   end subroutine test_frozen

   subroutine test_bad_frozen(error)
      !! Freezing every occupied orbital leaves nothing to correlate
      type(error_type), allocatable, intent(out) :: error
      type(mp2_result_t) :: mp2
      type(error_t) :: err

      call water_mp2("cc-pvdz", 5, mp2, err)
      call check(error, err%has_error(), "freezing all five occupied must be refused")
      if (allocated(error)) return

      call err%clear()
      call water_mp2("cc-pvdz", -1, mp2, err)
      call check(error, err%has_error(), "a negative frozen count must be refused")
   end subroutine test_bad_frozen

   subroutine test_scaling(error)
      !! The scaling belongs to `total`, not to the components
      !!
      !! Folding the factors into ss and os would give the same total and lose
      !! the numbers it was scaled from, so a result could never be rescaled or
      !! compared against an unscaled one. This is what keeps them separate.
      type(error_type), allocatable, intent(out) :: error
      type(mp2_energy_t) :: e

      e%ss = -0.05_dp
      e%os = -0.15_dp

      ! Default is plain MP2: both factors one.
      call check(error, abs(e%total() - (-0.20_dp)) < 1.0e-14_dp, &
                 "unscaled, the total is the plain sum")
      if (allocated(error)) return

      e%ss_scale = 1.0_dp/3.0_dp
      e%os_scale = 1.2_dp
      call check(error, abs(e%total() - (-0.15_dp*1.2_dp - 0.05_dp/3.0_dp)) < 1.0e-14_dp, &
                 "scaled, the total applies the factors")
      if (allocated(error)) return
      call check(error, abs(e%ss - (-0.05_dp)) < 1.0e-14_dp, &
                 "the same-spin component must be untouched by scaling")
      if (allocated(error)) return
      call check(error, abs(e%os - (-0.15_dp)) < 1.0e-14_dp, &
                 "the opposite-spin component must be untouched by scaling")
   end subroutine test_scaling

   subroutine test_scs(error)
      !! `scs` reports Grimme's factors whatever the run used
      type(error_type), allocatable, intent(out) :: error
      type(mp2_energy_t) :: e

      e%ss = -0.05_dp
      e%os = -0.15_dp
      ! Deliberately set to something else: scs answers "what would SCS-MP2
      ! give", which is a different question from "what did this run give".
      e%ss_scale = 0.0_dp
      e%os_scale = 1.3_dp

      call check(error, abs(e%scs() - (-0.15_dp*1.2_dp - 0.05_dp/3.0_dp)) < 1.0e-14_dp, &
                 "scs must use the published factors, not the run's")
      if (allocated(error)) return
      call check(error, abs(e%total() - (-0.15_dp*1.3_dp)) < 1.0e-14_dp, &
                 "and total must still use the run's")
   end subroutine test_scs

end module test_mqc_libcint_mp2

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_libcint_mp2, only: collect_mqc_libcint_mp2_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_libcint_mp2", collect_mqc_libcint_mp2_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
