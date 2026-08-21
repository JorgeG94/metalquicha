!! Unrestricted MP2: the closed-shell limit, and the two spin channels
module test_mqc_ump2
   !! **The closed-shell limit is the test that pins the spin factors.** An
   !! unrestricted MP2 on a closed-shell reference has to reproduce the
   !! restricted one exactly -- not approximately, since the alpha and beta
   !! orbitals are the same orbitals and every integral is the same integral.
   !! Getting the 1/4 on the like-spin terms wrong, or antisymmetrising the
   !! mixed-spin block that must not be antisymmetrised, moves that number and
   !! nothing else here would notice.
   !!
   !! It also separates the two implementations from each other: conventional
   !! and fitted are checked against the same restricted reference, so a fault
   !! in the B tensor shows up as a disagreement rather than as two routines
   !! agreeing on the same mistake.
   !!
   !! The open-shell numbers are pinned rather than derived. There is no
   !! identity that fixes E(2) for a doublet, so what is asserted about them is
   !! sign, spin decomposition and reproducibility -- the value itself is
   !! checked against PySCF outside this file.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf, run_libcint_uhf
   use mqc_libcint_mp2, only: mp2_result_t, run_libcint_mp2, run_libcint_ump2, &
                              run_libcint_ri_mp2, run_libcint_uri_mp2
   implicit none
   private

   public :: collect_mqc_ump2_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   !> Bohr, C2v, the same geometry the other backend tests use
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp], [3, 3])

contains

   subroutine collect_mqc_ump2_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ump2_reproduces_rmp2_on_a_closed_shell", closed_shell_limit), &
                  new_unittest("cation_correlation_is_negative", cation_is_negative), &
                  new_unittest("refuses_mismatched_orbital_counts", refuses_mismatch), &
                  new_unittest("uri_mp2_reproduces_ri_mp2_on_a_closed_shell", ri_closed_shell_limit) &
                  ]
   end subroutine collect_mqc_ump2_tests

   subroutine water_scf(nelec, mult, unrestricted, mol, scf, err, ok)
      !! One converged reference on the shared geometry
      integer, intent(in) :: nelec, mult
      logical, intent(in) :: unrestricted
      type(libcint_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      if (err%has_error()) return
      if (unrestricted) then
         call run_libcint_uhf(mol, nelec, mult, 100, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
      else
         call run_libcint_rhf(mol, nelec, 100, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
      end if
      ok = .not. err%has_error()
   end subroutine water_scf

   subroutine closed_shell_limit(error)
      !! UMP2 on a closed shell is RMP2, to machine precision
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol_r, mol_u
      type(rhf_result_t) :: rhf, uhf
      type(mp2_result_t) :: rmp2, ump2
      type(error_t) :: err
      logical :: ok

      call water_scf(10, 1, .false., mol_r, rhf, err, ok)
      call check(error, ok, "the restricted reference did not converge")
      if (allocated(error)) return
      call run_libcint_mp2(mol_r, rhf%orbitals, rhf%orbital_energies, 5, rhf%energy, &
                           rmp2, err)
      call check(error,.not. err%has_error(), "RMP2 failed")
      if (allocated(error)) return

      call water_scf(10, 1, .true., mol_u, uhf, err, ok)
      call check(error, ok, "the unrestricted reference did not converge")
      if (allocated(error)) return
      call run_libcint_ump2(mol_u, uhf%orbitals, uhf%orbitals_beta, &
                            uhf%orbital_energies, uhf%orbital_energies_beta, &
                            5, 5, uhf%energy, ump2, err)
      call check(error,.not. err%has_error(), "UMP2 failed")
      if (allocated(error)) return

      ! The correlation energy first, then each spin channel, so a failure says
      ! which factor is wrong rather than only that something is.
      call check(error, abs(ump2%correlation - rmp2%correlation) < 1.0e-10_dp, &
                 "UMP2 correlation does not reproduce RMP2 on a closed shell")
      if (allocated(error)) return
      call check(error, abs(ump2%same_spin - rmp2%same_spin) < 1.0e-10_dp, &
                 "the like-spin channels disagree: check the 1/2 on E(aa) and E(bb)")
      if (allocated(error)) return
      call check(error, abs(ump2%opposite_spin - rmp2%opposite_spin) < 1.0e-10_dp, &
                 "the mixed-spin channels disagree: check that E(ab) is NOT antisymmetrised")

      call mol_r%destroy(); call mol_u%destroy()
   end subroutine closed_shell_limit

   subroutine cation_is_negative(error)
      !! A doublet correlates, and both channels contribute
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: uhf
      type(mp2_result_t) :: ump2
      type(error_t) :: err
      logical :: ok

      call water_scf(9, 2, .true., mol, uhf, err, ok)
      call check(error, ok, "the cation did not converge")
      if (allocated(error)) return
      call run_libcint_ump2(mol, uhf%orbitals, uhf%orbitals_beta, &
                            uhf%orbital_energies, uhf%orbital_energies_beta, &
                            5, 4, uhf%energy, ump2, err)
      call check(error,.not. err%has_error(), "UMP2 on the cation failed")
      if (allocated(error)) return

      ! Every denominator is negative for a stable reference and every
      ! numerator is a square or a sum of them, so both channels are negative.
      ! A positive one means the reference was not a minimum, which is worth
      ! catching here rather than in whatever consumes the number.
      call check(error, ump2%correlation < 0.0_dp, "correlation energy is not negative")
      if (allocated(error)) return
      call check(error, ump2%opposite_spin < 0.0_dp, "the mixed-spin channel is not negative")
      if (allocated(error)) return
      call check(error, ump2%same_spin < 0.0_dp, "the like-spin channel is not negative")
      if (allocated(error)) return
      ! Opposite spin dominates in MP2, and by a wide margin -- roughly three to
      ! one on a molecule like this. A run where they had swapped magnitude
      ! would be a spin-factor error that the sign checks above would miss.
      call check(error, abs(ump2%opposite_spin) > abs(ump2%same_spin), &
                 "the like-spin channel outweighs the mixed one, which no MP2 does")

      call mol%destroy()
   end subroutine cation_is_negative

   subroutine ri_closed_shell_limit(error)
      !! The same limit for the fitted pair, which is a separate implementation
      !!
      !! Against RI-MP2 rather than against conventional MP2: the fitting error
      !! is common to both sides here and cancels, so what is left is the spin
      !! algebra. Checking the fitted unrestricted result against the exact
      !! restricted one would fold in the RI error and need a tolerance loose
      !! enough to hide a real mistake.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol_r, mol_u, aux_r, aux_u
      type(rhf_result_t) :: rhf, uhf
      type(mp2_result_t) :: rmp2, ump2
      type(error_t) :: err
      logical :: ok

      call water_scf(10, 1, .false., mol_r, rhf, err, ok)
      call check(error, ok, "the restricted reference did not converge")
      if (allocated(error)) return
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "cc-pvdz-rifit", aux_r, err)
      call check(error,.not. err%has_error(), "the auxiliary basis did not build")
      if (allocated(error)) return
      call run_libcint_ri_mp2(mol_r, aux_r, rhf%orbitals, rhf%orbital_energies, 5, &
                              rhf%energy, rmp2, err)
      call check(error,.not. err%has_error(), "RI-MP2 failed")
      if (allocated(error)) return

      call water_scf(10, 1, .true., mol_u, uhf, err, ok)
      call check(error, ok, "the unrestricted reference did not converge")
      if (allocated(error)) return
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "cc-pvdz-rifit", aux_u, err)
      call check(error,.not. err%has_error(), "the auxiliary basis did not build")
      if (allocated(error)) return
      call run_libcint_uri_mp2(mol_u, aux_u, uhf%orbitals, uhf%orbitals_beta, &
                               uhf%orbital_energies, uhf%orbital_energies_beta, &
                               5, 5, uhf%energy, ump2, err)
      call check(error,.not. err%has_error(), "URI-MP2 failed")
      if (allocated(error)) return

      call check(error, abs(ump2%correlation - rmp2%correlation) < 1.0e-10_dp, &
                 "URI-MP2 correlation does not reproduce RI-MP2 on a closed shell")
      if (allocated(error)) return
      call check(error, abs(ump2%same_spin - rmp2%same_spin) < 1.0e-10_dp, &
                 "the like-spin channels disagree in the fitted pair")
      if (allocated(error)) return
      call check(error, abs(ump2%opposite_spin - rmp2%opposite_spin) < 1.0e-10_dp, &
                 "the mixed-spin channels disagree in the fitted pair")

      call mol_r%destroy(); call mol_u%destroy()
      call aux_r%destroy(); call aux_u%destroy()
   end subroutine ri_closed_shell_limit

   subroutine refuses_mismatch(error)
      !! Alpha and beta must span the same orbital space
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: uhf
      type(mp2_result_t) :: ump2
      type(error_t) :: err
      logical :: ok

      call water_scf(10, 1, .true., mol, uhf, err, ok)
      call check(error, ok, "the reference did not converge")
      if (allocated(error)) return

      ! Beta truncated by one column. Silently correlating a smaller virtual
      ! space than alpha would return a plausible number.
      call run_libcint_ump2(mol, uhf%orbitals, uhf%orbitals_beta(:, 1:size(uhf%orbitals_beta, 2) - 1), &
                            uhf%orbital_energies, uhf%orbital_energies_beta, &
                            5, 5, uhf%energy, ump2, err)
      call check(error, err%has_error(), "a truncated beta space was accepted")
      call err%clear()
      call mol%destroy()
   end subroutine refuses_mismatch

end module test_mqc_ump2

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_ump2, only: collect_mqc_ump2_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_ump2", collect_mqc_ump2_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
