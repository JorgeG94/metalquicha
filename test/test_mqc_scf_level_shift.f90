module test_mqc_scf_level_shift
   !! Level shifting changes how the SCF gets there, not where it arrives
   !!
   !! **These are invariance tests, and that is deliberate.** A convergence aid
   !! is normally shown off on a case that fails without it, and there is no
   !! such case here yet -- so what is pinned instead is the property that
   !! actually risks being broken. A shift raises the virtual block before each
   !! diagonalisation, so every virtual orbital energy is high by that amount
   !! while it is in force. Those energies are not decoration once the SCF
   !! returns: they are the denominators of MP2 and coupled cluster, the weights
   !! of the gradient's energy-weighted density, and `eps_occ` and the response
   !! poles of a fragment potential. A shift left in at exit moves all of them
   !! by an amount nothing downstream could recognise as a shift.
   !!
   !! `run_czt_rhf` prevents that by making "the shift is off" part of the
   !! convergence test rather than something checked afterwards, so the last
   !! diagonalisation is always unshifted. These tests are what says that
   !! mechanism works, and they are written against the consequences -- the
   !! energy, the spectrum, and an MP2 built on it -- rather than against the
   !! mechanism, so that they still hold if it is ever rewritten.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, run_czt_uhf
   use mqc_czt_mp2, only: run_czt_mp2, mp2_result_t
   implicit none
   private

   public :: collect_level_shift_tests

   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [3, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]

   real(dp), parameter :: E_TOL = 1.0e-12_dp, D_TOL = 1.0e-10_dp
   real(dp), parameter :: SHIFT = 0.5_dp
      !! Large enough to change the iteration path outright -- half a Hartree on
      !! every virtual -- so that "the answer did not move" is a statement about
      !! the fixed point and not about a shift too small to do anything.

contains

   subroutine collect_level_shift_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("a shift does not move the RHF energy", &
                               shift_keeps_the_energy), &
                  new_unittest("a shift does not move the orbital energies", &
                               shift_keeps_the_spectrum), &
                  new_unittest("a shift does not move an MP2 built on it", &
                               shift_keeps_mp2), &
                  new_unittest("a shift does not move the UHF energy", &
                               shift_keeps_the_uhf_energy), &
                  new_unittest("a negative shift is refused", &
                               negative_shift_is_refused) &
                  ]
   end subroutine collect_level_shift_tests

   subroutine converge(error, shift, scf, mol, basis)
      !! One closed-shell SCF, optionally shifted
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: shift
      type(rhf_result_t), intent(out) :: scf
      type(czt_molecule_t), intent(out) :: mol
      character(len=*), intent(in), optional :: basis

      type(error_t) :: err
      character(len=16) :: set

      set = "sto-3g"
      if (present(basis)) set = basis

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, trim(set), mol, err)
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return

      call run_czt_rhf(mol, 10, 300, E_TOL, D_TOL, .false., scf, err, &
                       level_shift=shift)
      call check(error,.not. err%has_error(), "the SCF should run")
      if (allocated(error)) return
      call check(error, scf%converged, "the SCF should converge")
   end subroutine converge

   subroutine shift_keeps_the_energy(error)
      !! Same fixed point, reached along a different path
      type(error_type), allocatable, intent(out) :: error
      type(rhf_result_t) :: plain, shifted
      type(czt_molecule_t) :: mol

      call converge(error, 0.0_dp, plain, mol)
      if (allocated(error)) return
      call converge(error, SHIFT, shifted, mol)
      if (allocated(error)) return

      ! To the convergence threshold, not to machine precision: the two took
      ! different routes and stopped when the density stopped moving, so they
      ! agree to about as well as that test is tight.
      call check(error, abs(plain%energy - shifted%energy) < 1.0e-9_dp, &
                 "a level shift must not move the converged energy")
   end subroutine shift_keeps_the_energy

   subroutine shift_keeps_the_spectrum(error)
      !! The test the whole design exists to pass
      !!
      !! Were the shift still in force on the iteration that converged, every
      !! virtual eigenvalue here would be high by `SHIFT` -- half a Hartree,
      !! which is nine orders above the bound below. Nothing else in the suite
      !! would notice.
      type(error_type), allocatable, intent(out) :: error
      type(rhf_result_t) :: plain, shifted
      type(czt_molecule_t) :: mol
      real(dp) :: worst

      call converge(error, 0.0_dp, plain, mol)
      if (allocated(error)) return
      call converge(error, SHIFT, shifted, mol)
      if (allocated(error)) return

      call check(error, size(plain%orbital_energies) == size(shifted%orbital_energies), &
                 "the two runs should describe the same number of orbitals")
      if (allocated(error)) return

      worst = maxval(abs(plain%orbital_energies - shifted%orbital_energies))
      call check(error, worst < 1.0e-8_dp, &
                 "a level shift must not survive into the reported orbital energies")
   end subroutine shift_keeps_the_spectrum

   subroutine shift_keeps_mp2(error)
      !! The same guard through a consumer that would not complain on its own
      !!
      !! `run_czt_mp2` takes the orbital energies as its denominators, so a
      !! residual shift arrives here as a correlation energy that is wrong and
      !! entirely plausible.
      type(error_type), allocatable, intent(out) :: error
      type(rhf_result_t) :: plain, shifted
      type(czt_molecule_t) :: mol
      type(mp2_result_t) :: mp2_plain, mp2_shifted
      type(error_t) :: err

      call converge(error, 0.0_dp, plain, mol, basis="6-31g")
      if (allocated(error)) return
      call run_czt_mp2(mol, plain%orbitals, plain%orbital_energies, &
                       plain%n_occupied, plain%energy, mp2_plain, err)
      call check(error,.not. err%has_error(), "the reference MP2 should run")
      if (allocated(error)) return

      call converge(error, SHIFT, shifted, mol, basis="6-31g")
      if (allocated(error)) return
      call run_czt_mp2(mol, shifted%orbitals, shifted%orbital_energies, &
                       shifted%n_occupied, shifted%energy, mp2_shifted, err)
      call check(error,.not. err%has_error(), "the shifted-run MP2 should run")
      if (allocated(error)) return

      call check(error, abs(mp2_plain%correlation - mp2_shifted%correlation) < 1.0e-9_dp, &
                 "a level shift must not move the correlation energy built on it")
   end subroutine shift_keeps_mp2

   subroutine shift_keeps_the_uhf_energy(error)
      !! The unrestricted path shifts each spin against its own virtuals
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: plain, shifted
      type(error_t) :: err

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return

      ! The cation, so there is an unpaired electron for the two spins to
      ! disagree about and the two Fock matrices are genuinely different.
      call run_czt_uhf(mol, 9, 2, 300, E_TOL, D_TOL, .false., plain, err)
      call check(error,.not. err%has_error() .and. plain%converged, &
                 "the unshifted UHF should converge")
      if (allocated(error)) return

      call run_czt_uhf(mol, 9, 2, 300, E_TOL, D_TOL, .false., shifted, err, &
                       level_shift=SHIFT)
      call check(error,.not. err%has_error() .and. shifted%converged, &
                 "the shifted UHF should converge")
      if (allocated(error)) return

      call check(error, abs(plain%energy - shifted%energy) < 1.0e-9_dp, &
                 "a level shift must not move the converged UHF energy")
   end subroutine shift_keeps_the_uhf_energy

   subroutine negative_shift_is_refused(error)
      !! Refused rather than clamped to zero, because it is a typo and not a request
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return

      call run_czt_rhf(mol, 10, 300, E_TOL, D_TOL, .false., scf, err, &
                       level_shift=-0.5_dp)
      call check(error, err%has_error(), "a negative level shift must be refused")
      if (allocated(error)) return
      call check(error, len_trim(err%get_message()) > 0, &
                 "and it must say why")
   end subroutine negative_shift_is_refused

end module test_mqc_scf_level_shift

program test_scf_level_shift_driver
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_scf_level_shift, only: collect_level_shift_tests
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   testsuites = [new_testsuite("scf_level_shift", collect_level_shift_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat)
   if (stat > 0) then
      write (error_unit, '(i0, a)') stat, " test(s) failed!"
      error stop
   end if
end program test_scf_level_shift_driver
