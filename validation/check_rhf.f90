!! Manual check that the CPU RHF gives the right total energy
!!
!!     cmake -B build -DMQC_ENABLE_CZT=ON && ./build/check_rhf
!!
!! Two systems, both with answers that exist independently of this program.
!! The references below are PySCF 2.14 on the same geometries and the same
!! basis data -- not the textbook roundings the tolerances were first written
!! against. Fed the same BSE numbers, the two codes agree to about 1e-11,
!! which is the SCF convergence floor rather than a difference of method.
!!
!!   * H2 / STO-3G at R = 1.4 bohr, -1.1167 Ha. s functions only, so it tests
!!     the SCF itself: orthogonalisation, the Fock build, the density, the
!!     energy expression.
!!   * H2O / STO-3G at the standard geometry, -74.9659 Ha. This one has p
!!     shells, so it also tests angular momentum handling and the ordering
!!     libcint returns spherical functions in -- which the first system
!!     cannot, having none.
!!
!! Both matter. A code that gets H2 right and water wrong has an offset bug
!! in the shell blocking, and that is exactly the failure that looks like a
!! converged answer.
!!
!! The basis is written out here rather than read through the BSE reader, as
!! in check_libcint: a wrong number should be this backend's fault and not the
!! reader's. Wiring the reader in is the next step.
program check_rhf
   use pic_types, only: dp
   use mqc_cgto, only: molecular_basis_type
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_error, only: error_t
   implicit none

   integer :: failures

   failures = 0
   call run_h2()
   call run_water()

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[rhf] all ok -- CPU Hartree-Fock reproduces both references"
   else
      write (*, "(A,I0,A)") "[rhf] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine run_h2()
      !! H2, STO-3G, R = 1.4 bohr. Reference -1.1167 Ha.
      type(molecular_basis_type) :: basis
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: error
      real(dp) :: coords(3, 2)
      integer :: z(2)

      z = [1, 1]
      coords = 0.0_dp
      coords(3, 2) = 1.4_dp

      call basis%allocate_elements(2)
      call hydrogen_sto3g(basis, 1)
      call hydrogen_sto3g(basis, 2)

      call mol%build(z, coords, basis, error)
      call expect(.not. error%has_error(), "H2 packed: "//error%get_message())
      call expect(mol%nao == 2, "H2/STO-3G has two basis functions")

      call run_czt_rhf(mol, 2, 100, 1.0e-10_dp, 1.0e-8_dp, .false., scf, error)
      call expect(.not. error%has_error(), "H2 SCF: "//error%get_message())
      call expect(scf%converged, "H2 SCF converged")

      write (*, "(A,F16.10,A,I0,A)") "H2   /STO-3G  E = ", scf%energy, &
         "  (", scf%iterations, " iterations)"
      write (*, "(A,F16.10,A,F16.10)") "              electronic ", scf%electronic, &
         "   nuclear ", scf%nuclear_repulsion
      call expect(abs(scf%nuclear_repulsion - 1.0_dp/1.4_dp) < 1.0e-12_dp, &
                  "H2 nuclear repulsion is 1/R")
      ! PySCF 2.14, same geometry and basis: -1.1167143251
      call expect(abs(scf%energy - (-1.1167143251_dp)) < 1.0e-9_dp, &
                  "H2 total energy matches PySCF to 1e-9")

      call mol%destroy()
      call basis%destroy()
   end subroutine run_h2

   subroutine run_water()
      !! H2O, STO-3G, the standard geometry. Reference -74.9659 Ha.
      !!
      !! r(OH) = 0.9689 angstrom, angle 103.9 degrees, placed in the yz plane.
      type(molecular_basis_type) :: basis
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf, plain
      type(error_t) :: error
      real(dp) :: coords(3, 3)
      integer :: z(3)

      z = [8, 1, 1]
      coords = 0.0_dp
      coords(:, 1) = [0.0_dp, 0.0_dp, -0.1364652_dp]
      coords(:, 2) = [0.0_dp, 1.4304924_dp, 1.0826636_dp]
      coords(:, 3) = [0.0_dp, -1.4304924_dp, 1.0826636_dp]

      call basis%allocate_elements(3)
      call oxygen_sto3g(basis, 1)
      call hydrogen_sto3g(basis, 2)
      call hydrogen_sto3g(basis, 3)

      call mol%build(z, coords, basis, error)
      call expect(.not. error%has_error(), "H2O packed: "//error%get_message())
      ! O contributes 1s + 2s + 2p = 5, each H one. The count is worth
      ! asserting: a p shell mis-sized would give 3 or 7 and every later
      ! number would still look plausible.
      call expect(mol%nao == 7, "H2O/STO-3G has seven basis functions")

      call run_czt_rhf(mol, 10, 100, 1.0e-10_dp, 1.0e-8_dp, .false., scf, error)
      call expect(.not. error%has_error(), "H2O SCF: "//error%get_message())
      call expect(scf%converged, "H2O SCF converged")

      ! The same SCF with the extrapolation switched off. DIIS is allowed to
      ! change how many iterations this takes and nothing else; if the two
      ! answers differ, the extrapolation is altering the fixed point rather
      ! than finding it faster, which no iteration count would reveal.
      call run_czt_rhf(mol, 10, 200, 1.0e-10_dp, 1.0e-8_dp, .false., plain, error, &
                       diis_vectors=0)
      call expect(.not. error%has_error(), "H2O SCF without DIIS: "//error%get_message())
      call expect(plain%converged, "H2O SCF converged without DIIS")
      call expect(abs(scf%energy - plain%energy) < 1.0e-9_dp, &
                  "DIIS and plain iteration reach the same energy")
      write (*, "(A,I0,A,I0,A)") "              DIIS ", scf%iterations, &
         " iterations against ", plain%iterations, " without"
      call expect(scf%iterations < plain%iterations, &
                  "DIIS actually converges in fewer iterations")

      write (*, "(A,F16.10,A,I0,A)") "H2O  /STO-3G  E = ", scf%energy, &
         "  (", scf%iterations, " iterations)"
      write (*, "(A,F16.10,A,F16.10)") "              electronic ", scf%electronic, &
         "   nuclear ", scf%nuclear_repulsion
      ! PySCF 2.14, same geometry and basis: -74.9658162796
      call expect(abs(scf%energy - (-74.9658162796_dp)) < 1.0e-9_dp, &
                  "H2O total energy matches PySCF to 1e-9")

      call mol%destroy()
      call basis%destroy()
   end subroutine run_water

   subroutine hydrogen_sto3g(basis, iatom)
      !! One s shell
      type(molecular_basis_type), intent(inout) :: basis
      integer, intent(in) :: iatom

      basis%elements(iatom)%element = "H"
      call basis%elements(iatom)%allocate_shells(1)
      call set_shell(basis, iatom, 1, 0, &
                     [3.42525091_dp, 0.62391373_dp, 0.16885540_dp], &
                     [0.15432897_dp, 0.53532814_dp, 0.44463454_dp])
   end subroutine hydrogen_sto3g

   subroutine oxygen_sto3g(basis, iatom)
      !! 1s, 2s and 2p
      type(molecular_basis_type), intent(inout) :: basis
      integer, intent(in) :: iatom

      basis%elements(iatom)%element = "O"
      call basis%elements(iatom)%allocate_shells(3)
      call set_shell(basis, iatom, 1, 0, &
                     [130.7093200_dp, 23.8088610_dp, 6.4436083_dp], &
                     [0.15432897_dp, 0.53532814_dp, 0.44463454_dp])
      call set_shell(basis, iatom, 2, 0, &
                     [5.0331513_dp, 1.1695961_dp, 0.3803890_dp], &
                     [-0.09996723_dp, 0.39951283_dp, 0.70011547_dp])
      call set_shell(basis, iatom, 3, 1, &
                     [5.0331513_dp, 1.1695961_dp, 0.3803890_dp], &
                     [0.15591627_dp, 0.60768372_dp, 0.39195739_dp])
   end subroutine oxygen_sto3g

   subroutine set_shell(basis, iatom, ishell, ang, exponents, coefficients)
      type(molecular_basis_type), intent(inout) :: basis
      integer, intent(in) :: iatom, ishell, ang
      real(dp), intent(in) :: exponents(:), coefficients(:)

      basis%elements(iatom)%shells(ishell)%ang_mom = ang
      call basis%elements(iatom)%shells(ishell)%allocate_arrays(size(exponents))
      basis%elements(iatom)%shells(ishell)%exponents = exponents
      basis%elements(iatom)%shells(ishell)%coefficients = coefficients
   end subroutine set_shell

   subroutine expect(condition, what)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: what

      if (.not. condition) then
         write (*, "(A)") "  FAILED: "//trim(what)
         failures = failures + 1
      end if
   end subroutine expect

end program check_rhf
