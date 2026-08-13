!! The CPU analytic SCF gradient, against finite differences and translation
program check_scf_gradient
   !! Three checks, and they fail in different ways on purpose.
   !!
   !! **Finite difference of this program's own energy** is the one that
   !! catches a wrong gradient of a right energy. It compares the analytic
   !! derivative against the function it claims to differentiate, so it cannot
   !! be fooled by an error shared between the two.
   !!
   !! **Translational invariance**, `sum_A dE/dR_A = 0`, is free and catches
   !! exactly one thing: a missing or double-counted Hellmann-Feynman term.
   !! Move every nucleus and every basis function together and the energy
   !! cannot change, so the gradient must sum to zero whatever the geometry.
   !!
   !! **The printed values** are for comparison against PySCF fed this
   !! repository's own basis JSON. Not asserted here -- a reference that has to
   !! be regenerated when a basis file changes belongs in the validation
   !! manifest, not compiled in.
   !!
   !! Not a ctest case: an SCF per displaced geometry is seconds, not
   !! milliseconds, and the unit suite is meant to stay instant.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf, run_libcint_uhf
   use mqc_libcint_gradient, only: libcint_scf_gradient
   implicit none

   integer, parameter :: N_DIM = 3
   integer :: n_bad

   n_bad = 0

   write (*, "(a)") "== CPU SCF gradients"

   call check_case("H2 / sto-3g", [1, 1], ["H", "H"], &
                   reshape([0.0_dp, 0.0_dp, -0.7_dp, &
                            0.0_dp, 0.0_dp, 0.7_dp], [N_DIM, 2]), &
                   "sto-3g", 2, 1, n_bad)

   call check_case("H2O / sto-3g", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", 10, 1, n_bad)

   call check_case("H2O / 6-31g", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "6-31g", 10, 1, n_bad)

   ! Asymmetric on purpose: a symmetric molecule can hide a term that is wrong
   ! in a way symmetry cancels.
   call check_case("HCN / sto-3g", [1, 6, 7], ["H", "C", "N"], &
                   reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                   "sto-3g", 14, 1, n_bad)

   call check_case("CH3 doublet / sto-3g (UHF)", [6, 1, 1, 1], ["C", "H", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            2.05_dp, 0.0_dp, 0.0_dp, &
                            -1.02_dp, 1.78_dp, 0.0_dp, &
                            -1.02_dp, -1.78_dp, 0.3_dp], [N_DIM, 4]), &
                   "sto-3g", 9, 2, n_bad)

   write (*, "(a)") ""
   if (n_bad == 0) then
      write (*, "(a)") "all gradient checks passed"
   else
      write (*, "(a,i0,a)") "FAILED: ", n_bad, " check(s)"
      error stop 1
   end if

contains

   subroutine check_case(label, numbers, symbols, coords, basis, nelec, multiplicity, n_bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec, multiplicity
      integer, intent(inout) :: n_bad

      real(dp), allocatable :: analytic(:, :), numeric(:, :)
      real(dp) :: worst, translation(3)
      integer :: natm, ia, ic
      type(error_t) :: error

      natm = size(numbers)

      write (*, "(a)") ""
      write (*, "(a,a)") "== ", label

      call gradient_at(numbers, symbols, coords, basis, nelec, multiplicity, &
                       analytic, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      call numeric_gradient(numbers, symbols, coords, basis, nelec, multiplicity, &
                            numeric, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL (finite difference): ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      write (*, "(a)") "  atom        analytic (x, y, z)                      finite difference"
      do ia = 1, natm
         write (*, "(i6,3f14.9,a,3f14.9)") ia, analytic(:, ia), "   ", numeric(:, ia)
      end do

      worst = maxval(abs(analytic - numeric))
      write (*, "(a,es12.4)") "  largest deviation from finite difference: ", worst
      if (worst > 2.0e-5_dp) then
         write (*, "(a)") "  FAIL: the analytic gradient does not differentiate this energy"
         n_bad = n_bad + 1
      end if

      do ic = 1, 3
         translation(ic) = sum(analytic(ic, :))
      end do
      write (*, "(a,es12.4)") "  |sum over atoms| (should be zero): ", maxval(abs(translation))
      if (maxval(abs(translation)) > 1.0e-8_dp) then
         write (*, "(a)") "  FAIL: not translationally invariant"
         n_bad = n_bad + 1
      end if
   end subroutine check_case

   subroutine gradient_at(numbers, symbols, coords, basis, nelec, multiplicity, &
                          gradient, error)
      !! Converge an SCF at this geometry and differentiate it
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec, multiplicity
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf

      call build_libcint_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) return

      if (multiplicity == 1) then
         call run_libcint_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, error)
         if (error%has_error()) return
         call libcint_scf_gradient(mol, scf%density, &
                                   orbitals=scf%orbitals, &
                                   orbital_energies=scf%orbital_energies, &
                                   n_occupied=scf%n_occupied, &
                                   gradient=gradient, error=error)
      else
         call run_libcint_uhf(mol, nelec, multiplicity, 200, 1.0e-12_dp, 1.0e-10_dp, &
                              .false., scf, error)
         if (error%has_error()) return
         call libcint_scf_gradient(mol, scf%density, &
                                   density_beta=scf%density_beta, &
                                   orbitals=scf%orbitals, &
                                   orbitals_beta=scf%orbitals_beta, &
                                   orbital_energies=scf%orbital_energies, &
                                   orbital_energies_beta=scf%orbital_energies_beta, &
                                   n_occupied=scf%n_occupied, &
                                   n_occupied_beta=scf%n_occupied_beta, &
                                   gradient=gradient, error=error)
      end if
   end subroutine gradient_at

   subroutine energy_at(numbers, symbols, coords, basis, nelec, multiplicity, energy, error)
      !! One converged SCF energy, for the finite difference
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec, multiplicity
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf

      energy = 0.0_dp
      call build_libcint_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) return

      if (multiplicity == 1) then
         call run_libcint_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, error)
      else
         call run_libcint_uhf(mol, nelec, multiplicity, 200, 1.0e-12_dp, 1.0e-10_dp, &
                              .false., scf, error)
      end if
      if (error%has_error()) return
      energy = scf%energy
   end subroutine energy_at

   subroutine numeric_gradient(numbers, symbols, coords, basis, nelec, multiplicity, &
                               gradient, error)
      !! Central differences on the SCF energy
      !!
      !! The step is a compromise the tolerance above is set from: too large
      !! and the neglected third derivative shows, too small and the SCF's own
      !! convergence noise does. 0.002 Bohr against a 1e-12 energy threshold
      !! puts both below 1e-5.
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec, multiplicity
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      real(dp), parameter :: step = 0.002_dp
      real(dp), allocatable :: shifted(:, :)
      real(dp) :: e_plus, e_minus
      integer :: natm, ia, ic

      natm = size(numbers)
      allocate (gradient(3, natm))
      allocate (shifted(3, natm))
      gradient = 0.0_dp

      do ia = 1, natm
         do ic = 1, 3
            shifted = coords
            shifted(ic, ia) = coords(ic, ia) + step
            call energy_at(numbers, symbols, shifted, basis, nelec, multiplicity, e_plus, error)
            if (error%has_error()) return

            shifted = coords
            shifted(ic, ia) = coords(ic, ia) - step
            call energy_at(numbers, symbols, shifted, basis, nelec, multiplicity, e_minus, error)
            if (error%has_error()) return

            gradient(ic, ia) = (e_plus - e_minus)/(2.0_dp*step)
         end do
      end do
   end subroutine numeric_gradient

end program check_scf_gradient
