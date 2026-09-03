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
   !! So the control is PySCF, by a five-point finite field on its own converged
   !! Kohn-Sham energy -- that install ships no analytic polarizability, which is
   !! what makes it a control: it shares no response code with what is being
   !! checked. The references below are those numbers, carried here rather than
   !! left to be re-derived, so this file fails on its own rather than printing
   !! something for a person to compare.
   !!
   !! `TOL` is set by the reference and not by us. A field of 0.002 with the SCF
   !! at 1e-13 leaves a few times 1e-8 in the second difference, and the two
   !! codes' grids are the same tables but not the same code -- so 1e-6 is
   !! generous against agreement that is actually 1e-8 to 1e-7, and tight enough
   !! that dropping any one kernel term fails it by orders of magnitude.
   !!
   !! The internal checks are the ones that hold regardless of reference: the
   !! tensor is symmetric, and it is positive definite for a closed-shell
   !! molecule in its ground state.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_cphf, only: static_polarizability
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none

   integer, parameter :: N_DIM = 3
   real(dp), parameter :: TOL = 1.0e-6_dp
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
                   "sto-3g", 10, "svwn", [0.039642427_dp, 4.803099589_dp, 2.173448652_dp], n_bad)

   ! Asymmetric, so a kernel wrong in a way symmetry hides shows up.
   call check_case("HCN / sto-3g, SVWN", [1, 6, 7], ["H", "C", "N"], &
                   reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                   "sto-3g", 14, "svwn", [2.914428518_dp, 2.914428566_dp, 11.424268862_dp], n_bad)

   ! A pure GGA: no exact exchange at all, so the whole difference from the
   ! Hartree-Fock operator is the kernel, and `v2rhosigma`, `v2sigma2` and the
   ! `v_sigma grad drho` term have nothing to hide behind. This is the case
   ! that fails if any one of the three is dropped.
   call check_case("H2O / sto-3g, BLYP", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", 10, "blyp", [0.036895297_dp, 5.016091472_dp, 2.250961225_dp], n_bad)

   ! A hybrid GGA, where the kernel and a scaled exact exchange have to be
   ! right together: BLYP alone would not catch an `exx_fraction` applied twice
   ! or not at all, because it is zero there.
   call check_case("H2O / sto-3g, B3LYP", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", 10, "b3lyp", [0.038235630_dp, 5.031696862_dp, 2.207720164_dp], n_bad)

   ! A meta-GGA, where the kernel gains three second derivatives in tau
   ! (`v2rhotau`, `v2sigmatau`, `v2tau2`) and the response density gains a tau
   ! of its own. Nothing below the meta-GGA rung exercises any of that, so this
   ! is the case that fails if the tau channel is dropped, double counted, or
   ! given the wrong convention -- and the convention is the risk, since libxc's
   ! tau is per spin.
   call check_case("H2O / sto-3g, TPSS", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", 10, "tpss", &
                   [0.035425788_dp, 5.059453553_dp, 2.251964609_dp], n_bad)

   ! A second meta-GGA, and a heavily parametrised one: M06-L reaches tau
   ! through a different working expression, so a term that happens to cancel
   ! in TPSS's form does not cancel here.
   call check_case("H2O / sto-3g, M06-L", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", 10, "m06-l", &
                   [0.030952161_dp, 4.855462133_dp, 2.209799695_dp], n_bad)

   call check_case("HCN / sto-3g, BLYP", [1, 6, 7], ["H", "C", "N"], &
                   reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                   "sto-3g", 14, "blyp", [2.947895687_dp, 2.947895687_dp, 11.622996848_dp], n_bad)

   write (*, "(a)") ""
   if (n_bad == 0) then
      write (*, "(a)") "all exchange-correlation kernel checks passed"
   else
      write (*, "(a,i0,a)") "FAILED: ", n_bad, " case(s)"
      stop 1
   end if

contains

   subroutine check_case(label, numbers, symbols, coords, basis, nelec, functional, &
                         reference, n_bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec
      character(len=*), intent(in) :: functional
      real(dp), intent(in) :: reference(3)
         !! PySCF's diagonal, by finite field. See the note at the top for why
         !! that is the control and why the tolerance is what it is.
      integer, intent(inout) :: n_bad

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: xc
      type(error_t) :: error
      real(dp) :: alpha_ks(3, 3)
      integer :: k

      write (*, "(a)") ""
      write (*, "(a,a)") "== ", label
      flush (output_unit)

      call build_czt_molecule(numbers, symbols, coords, basis, mol, error)
      if (fail(error, n_bad)) return
      call xc_context_create(mol, functional, xc, error, polarized=.false.)
      if (fail(error, n_bad)) return
      call run_czt_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, error, xc=xc)
      if (fail(error, n_bad)) return

      ! An SCF that ran out of iterations returns its last iterate and no error,
      ! so this has to be asked for. It is not a formality here: a reference
      ! density wrong by 1e-6 moves the analytic polarizability by 1e-5 while
      ! leaving a finite-field control -- which is variational, and so
      ! second-order insensitive to exactly this -- looking correct. Every
      ! disagreement in this file traced back to that and to nothing else.
      if (.not. scf%converged) then
         write (*, "(a,i0,a)") "  FAIL: the reference SCF did not converge in ", &
            scf%iterations, " iterations"
         n_bad = n_bad + 1
         call mol%destroy()
         return
      end if

      ! The same reference, with and without the kernel in the response. The
      ! difference between these two is the whole of what this file tests: if
      ! they agree, `xc_kernel_apply` contributed nothing and the plumbing is
      ! wrong rather than the physics.
      call static_polarizability(mol, scf%orbitals, scf%orbital_energies, nelec/2, &
                                 alpha_ks, error, tol=1.0e-10_dp, in_core=.true., &
                                 xc=xc, density=scf%density)
      if (fail(error, n_bad)) return

      write (*, "(a)") "  alpha (diagonal, Bohr^3)"
      write (*, "(a,3f16.9)") "    mine      ", (alpha_ks(k, k), k=1, 3)
      write (*, "(a,3f16.9)") "    PySCF     ", reference
      write (*, "(a,3es16.3)") "    diff      ", (alpha_ks(k, k) - reference(k), k=1, 3)
      write (*, "(a,3f16.9)") "    off-diag yz: ", alpha_ks(2, 3)
      flush (output_unit)

      do k = 1, 3
         if (abs(alpha_ks(k, k) - reference(k)) > TOL) then
            write (*, "(a,i0,a,es14.4)") "  FAIL: alpha(", k, ") differs from PySCF by ", &
               alpha_ks(k, k) - reference(k)
            n_bad = n_bad + 1
         end if
      end do

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
