!! Manual check that the coupled-perturbed equations and the polarizability agree
!! with PySCF
!!
!!     cmake -B build -DMQC_ENABLE_LIBCINT=ON
!!     ./build/check_cphf && python3 validation/check_cphf.py
!!
!! Milestone 3 of the EFP plan. `test/test_mqc_libcint_cphf.f90` already pins
!! this against a finite field through our own SCF, which fixes the sign, the
!! factor and the coupling without any outside help -- so what is added here is
!! an independent *implementation*: PySCF's integrals, PySCF's SCF, PySCF's
!! response solver.
!!
!! **aug-cc-pVDZ is the case that matters.** A polarizability is a measure of how
!! far the density will move, so it needs functions out where the density is
!! going. STO-3G gives water an out-of-plane polarizability of 0.04 Bohr^3
!! against a true value near 9 -- the minimal basis has no out-of-plane
!! flexibility at all. That makes the small bases a fine test of the *equations*
!! and a useless test of the physics, and the diffuse basis is what shows both
!! codes agree where the number is large and the virtual space is big enough for
!! the coupling to do real work.
program check_cphf
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_cphf, only: static_polarizability
   use mqc_error, only: error_t
   implicit none

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   integer :: failures

   failures = 0
   call one_case(1, "sto-3g", "water")
   call one_case(2, "cc-pvdz", "water")
   call one_case(3, "aug-cc-pvdz", "water")
   ! HCN in 6-31G*, which is where GAMESS's own AO-basis CPHF gives up:
   ! `TOO MANY ITERATIONS IN AOCPCG`, residual stalling around 3e-3 through its
   ! 299-iteration grace limit, and MAKEFP aborts before it reaches the
   ! polarizabilities. Its diagnostic blames the wavefunction -- "CHANGE SCFTYP,
   ! OR CHECK HOMO/LUMO FILLING ORDER" -- but the reference is a closed shell with
   ! a 0.685 Hartree gap, and the orbital Hessian is positive definite with a
   ! preconditioned condition number of 6.3. So this case is here permanently: it
   ! is the cheapest evidence that the solver, not the physics, is what fails
   ! there, and it is the molecule the EFP plan wants for the dummy-point
   ! convention.
   call one_case(4, "6-31g*", "hcn")

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[cphf] wrote every case; run validation/check_cphf.py"
   else
      write (*, "(A,I0,A)") "[cphf] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine one_case(case_index, basis, molecule)
      integer, intent(in) :: case_index
      character(len=*), intent(in) :: basis
      character(len=*), intent(in) :: molecule   !! "water" or "hcn"

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp) :: c(3, 3), alpha(3, 3)
      real(dp) :: worst_asym, isotropic
      integer :: unit, k, l, iters, nelec
      integer :: z(3)
      character(len=2) :: symbols(3)
      character(len=48) :: name

      select case (molecule)
      case ("water")
         z = [8, 1, 1]
         symbols = ["O ", "H ", "H "]
         nelec = 10
         c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                      0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                      0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      case default
         ! Linear, so every off-diagonal polarizability is zero by symmetry and
         ! the two perpendicular components must come out equal -- a free check
         ! that water, being merely planar, cannot provide.
         z = [7, 6, 1]
         symbols = ["N ", "C ", "H "]
         nelec = 14
         c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                      0.0_dp, 0.0_dp, 1.1560_dp*ANG, &
                      0.0_dp, 0.0_dp, 2.2230_dp*ANG], [3, 3])
      end select

      call build_libcint_molecule(z, symbols, c, basis, mol, err)
      if (err%has_error()) then
         write (*, "(A,A)") "[cphf] basis failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call run_libcint_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, &
                           in_core=.true.)
      if (err%has_error() .or. .not. scf%converged) then
         write (*, "(A)") "[cphf] SCF failed"
         failures = failures + 1
         call mol%destroy()
         return
      end if

      call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                 scf%n_occupied, alpha, err, iterations=iters)
      if (err%has_error()) then
         write (*, "(A,A)") "[cphf] CPHF failed: ", err%get_message()
         failures = failures + 1
         call mol%destroy()
         return
      end if

      ! Symmetry is not imposed anywhere, so its residual is a free measure of
      ! how well the three independent CG solves converged.
      worst_asym = 0.0_dp
      do l = 1, 3
         do k = 1, 3
            worst_asym = max(worst_asym, abs(alpha(k, l) - alpha(l, k)))
         end do
      end do
      isotropic = (alpha(1, 1) + alpha(2, 2) + alpha(3, 3))/3.0_dp

      write (name, "(A,I0,A)") "/tmp/mqc_cphf_", case_index, ".txt"
      open (newunit=unit, file=trim(name), status="replace", action="write")
      write (unit, "(A)") trim(basis)
      write (unit, "(A)") trim(molecule)
      write (unit, "(I0,1X,I0)") mol%nao, scf%n_occupied
      write (unit, "(es25.16e3)") scf%energy
      write (unit, "(3es25.16e3)") alpha
      close (unit)

      write (*, "(A,A6,1X,A12,A,I0,A,I0,A,F12.6,A,ES9.2)") &
         "  ", trim(molecule), trim(basis), "  nao ", mol%nao, "  CG iters ", iters, &
         "   alpha_iso ", isotropic, "   asymmetry ", worst_asym

      if (worst_asym > 1.0e-8_dp) then
         write (*, "(A)") "    FAIL alpha is not symmetric"
         failures = failures + 1
      end if

      call mol%destroy()
   end subroutine one_case

end program check_cphf
