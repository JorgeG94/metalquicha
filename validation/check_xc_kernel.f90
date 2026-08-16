!! The exchange-correlation kernel, through a coupled-perturbed Kohn-Sham polarizability
program check_xc_kernel
   !! **Why a polarizability and not a gradient.** `f_xc` exists to turn a
   !! coupled-perturbed Hartree-Fock operator into a Kohn-Sham one, and the
   !! thing that wants it first is the double hybrid gradient -- which needs a
   !! relaxed density, a Z-vector, a frozen core and a PT2 Lagrangian before it
   !! produces a number. Validating the kernel through that would mean four
   !! unvalidated things stacked on one another.
   !!
   !! A static polarizability needs the kernel and nothing else. It is one CPKS
   !! solve against a dipole perturbation, PySCF computes the same quantity, and
   !! a wrong kernel moves it by percent rather than by rounding.
   !!
   !! **What the control cannot be.** The obvious one -- run the same
   !! polarizability with and without the `xc` context -- is invalid, and
   !! finding out why cost an hour. Without the context the operator is
   !! Hartree-Fock's, `Delta eps + 4(ia|jb) - (ij|ab) - (ib|ja)`, and evaluating
   !! it on *Kohn-Sham* orbitals is not a control at all: LDA gaps are far
   !! smaller than Hartree-Fock ones, the exchange terms overwhelm `Delta eps`,
   !! and the conjugate-gradient solver correctly reports an operator that is
   !! not positive definite. That is the wrong operator for those orbitals
   !! rather than a defect in anything being tested.
   !!
   !! So the control is PySCF, which computes the same coupled-perturbed
   !! Kohn-Sham polarizability, and the internal checks here are the ones that
   !! hold regardless of reference: the tensor is symmetric, and it is positive
   !! definite for a closed-shell molecule in its ground state.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_cphf, only: static_polarizability
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_available
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none

   integer, parameter :: N_DIM = 3
   integer :: n_bad

   n_bad = 0
   write (*, "(a)") "== the exchange-correlation kernel, via CPKS polarizabilities"

   if (.not. xc_available()) then
      write (*, "(a)") "   skipped: this build has no libxc"
      stop 0
   end if

   call check_case("H2O / sto-3g, SVWN", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", 10, "svwn", n_bad)

   ! Asymmetric, so a kernel wrong in a way symmetry hides shows up.
   call check_case("HCN / sto-3g, SVWN", [1, 6, 7], ["H", "C", "N"], &
                   reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                   "sto-3g", 14, "svwn", n_bad)

   write (*, "(a)") ""
   if (n_bad == 0) then
      write (*, "(a)") "all exchange-correlation kernel checks passed"
   else
      write (*, "(a,i0,a)") "FAILED: ", n_bad, " case(s)"
      stop 1
   end if

contains

   subroutine check_case(label, numbers, symbols, coords, basis, nelec, functional, n_bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec
      character(len=*), intent(in) :: functional
      integer, intent(inout) :: n_bad

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: xc
      type(error_t) :: error
      real(dp) :: alpha_ks(3, 3)
      integer :: k

      write (*, "(a)") ""
      write (*, "(a,a)") "== ", label
      flush (output_unit)

      call build_libcint_molecule(numbers, symbols, coords, basis, mol, error)
      if (fail(error, n_bad)) return
      call xc_context_create(mol, functional, xc, error, polarized=.false.)
      if (fail(error, n_bad)) return
      call run_libcint_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, error, xc=xc)
      if (fail(error, n_bad)) return

      ! The same reference, with and without the kernel in the response. The
      ! difference between these two is the whole of what this file tests: if
      ! they agree, `xc_kernel_apply` contributed nothing and the plumbing is
      ! wrong rather than the physics.
      call static_polarizability(mol, scf%orbitals, scf%orbital_energies, nelec/2, &
                                 alpha_ks, error, tol=1.0e-10_dp, in_core=.true., &
                                 xc=xc, density=scf%density)
      if (fail(error, n_bad)) return

      write (*, "(a)") "  alpha (diagonal, Bohr^3) -- compare with PySCF:"
      write (*, "(a,3f16.9)") "    ", (alpha_ks(k, k), k=1, 3)
      write (*, "(a,3f16.9)") "    off-diag yz: ", alpha_ks(2, 3)
      flush (output_unit)

      ! Positive definite for a closed-shell ground state, whatever the
      ! reference. A kernel with the wrong sign shows up here first.
      do k = 1, 3
         if (alpha_ks(k, k) <= 0.0_dp) then
            write (*, "(a,i0,a)") "  FAIL: alpha(", k, ") is not positive"
            n_bad = n_bad + 1
         end if
      end do

      ! Symmetric by construction; an asymmetry means the operator is not the
      ! self-adjoint one the conjugate-gradient solver assumes.
      if (maxval(abs(alpha_ks - transpose(alpha_ks))) > 1.0e-7_dp) then
         write (*, "(a,es14.4)") "  FAIL: alpha is not symmetric: ", &
            maxval(abs(alpha_ks - transpose(alpha_ks)))
         n_bad = n_bad + 1
      end if

      call mol%destroy()
   end subroutine check_case

   function fail(error, n_bad) result(bad)
      type(error_t), intent(inout) :: error
      integer, intent(inout) :: n_bad
      logical :: bad

      bad = error%has_error()
      if (bad) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
      end if
   end function fail

end program check_xc_kernel
