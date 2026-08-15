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
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_available
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

   ! Density fitting. Only the two-electron term changes -- a fitted SCF is
   ! still variational -- so these exercise the three-centre and metric
   ! derivatives against the fitted energy they belong to.
   call check_case("H2O / sto-3g + cc-pvdz-rifit (DF)", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", 10, 1, n_bad, aux_basis="cc-pvdz-rifit")

   call check_case("H2O / 6-31g + cc-pvdz-rifit (DF)", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "6-31g", 10, 1, n_bad, aux_basis="cc-pvdz-rifit")

   call check_case("HCN / sto-3g + cc-pvdz-rifit (DF)", [1, 6, 7], ["H", "C", "N"], &
                   reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                   "sto-3g", 14, 1, n_bad, aux_basis="cc-pvdz-rifit")

   call check_unrestricted_df_channels(n_bad)

   ! Kohn-Sham. SVWN is Slater exchange with VWN correlation -- a pure local
   ! functional, so the exchange-correlation gradient is exercised with no Fock
   ! exchange derivative mixed into it. Skipped rather than failed on a build
   ! without libxc, which is the default.
   if (xc_available()) then
      call check_case("H2O / sto-3g, SVWN (RKS)", [8, 1, 1], ["O", "H", "H"], &
                      reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                               0.0_dp, -1.4308_dp, 1.1078_dp, &
                               0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                      "sto-3g", 10, 1, n_bad, functional="svwn")

      call check_case("HCN / sto-3g, SVWN (RKS)", [1, 6, 7], ["H", "C", "N"], &
                      reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                               0.0_dp, 0.0_dp, 0.0_dp, &
                               0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                      "sto-3g", 14, 1, n_bad, functional="svwn")

      call check_case("CH3 doublet / sto-3g, SVWN (UKS)", [6, 1, 1, 1], &
                      ["C", "H", "H", "H"], &
                      reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                               2.05_dp, 0.0_dp, 0.0_dp, &
                               -1.02_dp, 1.78_dp, 0.0_dp, &
                               -1.02_dp, -1.78_dp, 0.3_dp], [N_DIM, 4]), &
                      "sto-3g", 9, 2, n_bad, functional="svwn")
   else
      write (*, "(a)") ""
      write (*, "(a)") "== Kohn-Sham cases skipped: this build has no libxc"
   end if

   write (*, "(a)") ""
   if (n_bad == 0) then
      write (*, "(a)") "all gradient checks passed"
   else
      write (*, "(a,i0,a)") "FAILED: ", n_bad, " check(s)"
      error stop 1
   end if

contains

   subroutine check_case(label, numbers, symbols, coords, basis, nelec, multiplicity, &
                         n_bad, aux_basis, functional)
      character(len=*), intent(in) :: label
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec, multiplicity
      integer, intent(inout) :: n_bad
      character(len=*), intent(in), optional :: aux_basis
         !! Present fits J and K, and differentiates the fitted energy
      character(len=*), intent(in), optional :: functional
         !! Present makes it Kohn-Sham, and adds the exchange-correlation term

      real(dp), allocatable :: analytic(:, :), numeric(:, :)
      real(dp) :: translation(3)
      real(dp) :: worst
      integer :: natm, ia, ic
      type(error_t) :: error

      natm = size(numbers)

      write (*, "(a)") ""
      write (*, "(a,a)") "== ", label

      call gradient_at(numbers, symbols, coords, basis, nelec, multiplicity, &
                       analytic, error, aux_basis, functional)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      call numeric_gradient(numbers, symbols, coords, basis, nelec, multiplicity, &
                            numeric, error, aux_basis, functional)
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

   subroutine check_unrestricted_df_channels(n_bad)
      !! The unrestricted density-fitted branch, on a closed shell
      !!
      !! mqc has no unrestricted density-fitted SCF -- `run_libcint_uhf` says
      !! so itself -- so there is no open-shell fitted energy to differentiate
      !! and no finite difference to check the unrestricted branch against.
      !!
      !! What can still be checked is the part that is easy to get wrong: the
      !! channel weights. Feed a closed-shell system through the unrestricted
      !! path as two identical half-filled spin channels, and the result must
      !! equal the restricted one exactly. The restricted expression carries a
      !! factor of two on the exchange density and one on the metric term; the
      !! unrestricted one carries one and a half per channel. If either is
      !! wrong the two disagree, and this is the only thing that would notice.
      integer, intent(inout) :: n_bad

      integer, parameter :: numbers(3) = [8, 1, 1]
      character(len=1), parameter :: symbols(3) = ["O", "H", "H"]
      real(dp) :: coords(N_DIM, 3)
      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: error
      real(dp), allocatable :: restricted(:, :), unrestricted(:, :)
      real(dp), allocatable :: half_density(:, :)
      real(dp) :: worst

      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, -1.4308_dp, 1.1078_dp, &
                        0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3])

      write (*, "(a)") ""
      write (*, "(a)") "== unrestricted DF branch against the restricted one (H2O / sto-3g)"

      call build_libcint_molecule(numbers, symbols, coords, "sto-3g", mol, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
         return
      end if
      call build_libcint_molecule(numbers, symbols, coords, "cc-pvdz-rifit", aux, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      call run_libcint_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, error, aux=aux)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      call libcint_scf_gradient(mol, scf%density, &
                                orbitals=scf%orbitals, &
                                orbital_energies=scf%orbital_energies, &
                                n_occupied=scf%n_occupied, &
                                gradient=restricted, error=error, aux=aux)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      ! One electron per spin orbital rather than two, the same orbitals in
      ! both channels: physically the same closed shell, expressed the way an
      ! unrestricted reference would.
      allocate (half_density(mol%nao, mol%nao))
      half_density = 0.5_dp*scf%density

      call libcint_scf_gradient(mol, half_density, &
                                density_beta=half_density, &
                                orbitals=scf%orbitals, &
                                orbitals_beta=scf%orbitals, &
                                orbital_energies=scf%orbital_energies, &
                                orbital_energies_beta=scf%orbital_energies, &
                                n_occupied=scf%n_occupied, &
                                n_occupied_beta=scf%n_occupied, &
                                gradient=unrestricted, error=error, aux=aux)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      worst = maxval(abs(restricted - unrestricted))
      write (*, "(a,es12.4)") "  largest difference between the two paths: ", worst
      if (worst > 1.0e-11_dp) then
         write (*, "(a)") "  FAIL: the unrestricted channel weights do not reproduce the restricted result"
         n_bad = n_bad + 1
      end if
   end subroutine check_unrestricted_df_channels

   subroutine gradient_at(numbers, symbols, coords, basis, nelec, multiplicity, &
                          gradient, error, aux_basis, functional)
      !! Converge an SCF at this geometry and differentiate it
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec, multiplicity
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error
      character(len=*), intent(in), optional :: aux_basis
      character(len=*), intent(in), optional :: functional

      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(xc_context_t) :: xc

      call build_libcint_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) return
      if (present(aux_basis)) then
         call build_libcint_molecule(numbers, symbols, coords, aux_basis, aux, error)
         if (error%has_error()) return
      end if

      ! Kohn-Sham. Kept apart from the Hartree-Fock branches below rather than
      ! threaded through them: the SCF and the gradient both need the same
      ! context, and building it twice would integrate the energy on one grid
      ! and differentiate it on another.
      if (present(functional)) then
         call xc_context_create(mol, functional, xc, error, &
                                polarized=(multiplicity /= 1))
         if (error%has_error()) return
         if (multiplicity == 1) then
            call run_libcint_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, &
                                 error, xc=xc)
            if (error%has_error()) return
            call libcint_scf_gradient(mol, scf%density, &
                                      orbitals=scf%orbitals, &
                                      orbital_energies=scf%orbital_energies, &
                                      n_occupied=scf%n_occupied, &
                                      gradient=gradient, error=error, xc=xc)
         else
            call run_libcint_uhf(mol, nelec, multiplicity, 200, 1.0e-12_dp, 1.0e-10_dp, &
                                 .false., scf, error, xc=xc)
            if (error%has_error()) return
            call libcint_scf_gradient(mol, scf%density, &
                                      density_beta=scf%density_beta, &
                                      orbitals=scf%orbitals, &
                                      orbitals_beta=scf%orbitals_beta, &
                                      orbital_energies=scf%orbital_energies, &
                                      orbital_energies_beta=scf%orbital_energies_beta, &
                                      n_occupied=scf%n_occupied, &
                                      n_occupied_beta=scf%n_occupied_beta, &
                                      gradient=gradient, error=error, xc=xc)
         end if
         return
      end if

      if (multiplicity == 1) then
         if (present(aux_basis)) then
            call run_libcint_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, &
                                 error, aux=aux)
            if (error%has_error()) return
            call libcint_scf_gradient(mol, scf%density, &
                                      orbitals=scf%orbitals, &
                                      orbital_energies=scf%orbital_energies, &
                                      n_occupied=scf%n_occupied, &
                                      gradient=gradient, error=error, aux=aux)
            return
         end if
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

   subroutine energy_at(numbers, symbols, coords, basis, nelec, multiplicity, energy, &
                        error, aux_basis, functional)
      !! One converged SCF energy, for the finite difference
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec, multiplicity
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error
      character(len=*), intent(in), optional :: aux_basis
      character(len=*), intent(in), optional :: functional

      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(xc_context_t) :: xc

      energy = 0.0_dp
      call build_libcint_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) return
      if (present(aux_basis)) then
         call build_libcint_molecule(numbers, symbols, coords, aux_basis, aux, error)
         if (error%has_error()) return
      end if

      if (present(functional)) then
         ! The same grid the analytic gradient differentiates. `xc_context_create`
         ! builds the grid from the geometry it is handed, so a displaced point
         ! gets a displaced grid -- which is what makes the finite difference a
         ! test of the grid response rather than of nothing.
         call xc_context_create(mol, functional, xc, error, &
                                polarized=(multiplicity /= 1))
         if (error%has_error()) return
         if (multiplicity == 1) then
            call run_libcint_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, &
                                 error, xc=xc)
         else
            call run_libcint_uhf(mol, nelec, multiplicity, 200, 1.0e-12_dp, 1.0e-10_dp, &
                                 .false., scf, error, xc=xc)
         end if
      else if (multiplicity == 1) then
         if (present(aux_basis)) then
            call run_libcint_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, &
                                 error, aux=aux)
         else
            call run_libcint_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, error)
         end if
      else
         call run_libcint_uhf(mol, nelec, multiplicity, 200, 1.0e-12_dp, 1.0e-10_dp, &
                              .false., scf, error)
      end if
      if (error%has_error()) return
      energy = scf%energy
   end subroutine energy_at

   subroutine numeric_gradient(numbers, symbols, coords, basis, nelec, multiplicity, &
                               gradient, error, aux_basis, functional)
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
      character(len=*), intent(in), optional :: aux_basis
      character(len=*), intent(in), optional :: functional

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
            call energy_at(numbers, symbols, shifted, basis, nelec, multiplicity, e_plus, &
                           error, aux_basis, functional)
            if (error%has_error()) return

            shifted = coords
            shifted(ic, ia) = coords(ic, ia) - step
            call energy_at(numbers, symbols, shifted, basis, nelec, multiplicity, e_minus, &
                           error, aux_basis, functional)
            if (error%has_error()) return

            gradient(ic, ia) = (e_plus - e_minus)/(2.0_dp*step)
         end do
      end do
   end subroutine numeric_gradient

end program check_scf_gradient
