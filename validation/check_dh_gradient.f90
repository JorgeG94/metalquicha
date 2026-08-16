!! The double hybrid gradient, against finite difference of its own energy
program check_dh_gradient
   !! **The check that catches the error this assembly is most likely to make.**
   !! A double hybrid's energy is a Kohn-Sham part plus a scaled perturbative
   !! one, and only the first is variational. Its gradient is therefore the
   !! Kohn-Sham gradient plus a correlation term that carries a Z-vector -- and
   !! the tempting mistake is to reuse the MP2 gradient whole, which returns the
   !! *total* MP2 gradient with a Hartree-Fock reference inside it. That gives a
   !! number of entirely plausible magnitude, with the wrong reference gradient
   !! buried in it.
   !!
   !! Finite difference of the energy this program itself computes is what
   !! distinguishes those. It answers "does the gradient belong to the energy" and
   !! nothing else -- no external code, no convention to align, and a doubled or
   !! missing reference shows up immediately rather than as a small discrepancy.
   !!
   !! **The difference is extrapolated, which is what makes this sharp.** A plain
   !! central difference at a usable step leaves its own truncation error --
   !! measured here at 2.4e-6, 5.9e-7 and 1.5e-7 for steps of 4e-3, 2e-3 and
   !! 1e-3, which is 4.00x per halving to three figures and therefore textbook
   !! `h^2`. A tolerance loose enough to admit that would also admit a real error
   !! of the same size. Taking two steps and combining them as
   !! `(4 D(h/2) - D(h)) / 3` removes the `h^2` term and leaves `h^4`, which
   !! brings the comparison down to where the SCF's own convergence is the floor
   !! -- and lets this check to 1e-8 rather than to 1e-6.
   !!
   !! It costs four converged SCFs per component instead of two. On systems this
   !! size that is nothing, and the alternative is a test that cannot distinguish
   !! a correct gradient from one wrong in its fifth digit.
   !!
   !! All electron. The perturbative term of a double hybrid is computed with no
   !! frozen core here (`mqc_method_dft` never sets one), so an all-electron
   !! gradient differentiates the energy that is actually produced.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf, SCF_GUESS_SAD
   use mqc_libcint_atomic_guess, only: build_atomic_guess
   use mqc_libcint_mp2, only: run_libcint_mp2, mp2_result_t
   use mqc_libcint_gradient, only: libcint_scf_gradient
   use mqc_libcint_mp2_gradient, only: libcint_mp2_gradient
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_available
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none

   integer, parameter :: N_DIM = 3
   real(dp), parameter :: STEP = 4.0e-3_dp
      !! The coarse step. The fine one is half of it, and the two are combined
      !! by Richardson. Deliberately not small: the pair has to straddle a range
      !! where `h^2` dominates, and a coarse step already down at the SCF's noise
      !! floor would extrapolate noise instead of removing truncation.
   real(dp), parameter :: TOL = 1.0e-8_dp
   real(dp), parameter :: DENSITY_TOL = 1.0e-9_dp
      !! The density residual to stop at, with the energy held to 1e-12.
      !!
      !! **Not a slackening, a floor.** Asking for 1e-10 here cost most of an
      !! afternoon: HCN/STO-3G plateaus at 8.4e-10 and then grinds to the
      !! iteration limit, which this harness reported as "did not converge" --
      !! and which looked in turn like a hard molecule, a bad guess and a DIIS
      !! that could not handle a degenerate pi pair. It was none of those. The
      !! exchange-correlation quadrature puts noise in the Fock matrix, so the
      !! commutator cannot go below it, and the same geometry through the driver
      !! converges in nine iterations to the same energy PySCF reaches. The
      !! energy is converged to 1e-13 either way, which is what the finite
      !! difference actually consumes.
   integer :: n_bad

   n_bad = 0
   write (*, "(a)") "== the double hybrid gradient"

   if (.not. xc_available()) then
      write (*, "(a)") "   skipped: this build has no libxc"
      stop 0
   end if

   call check_case("H2O / sto-3g, B2PLYP", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", 10, "b2plyp", n_bad)

   ! No symmetry left, and a different PT2 coefficient: B2GP-PLYP takes 0.36
   ! where B2PLYP takes 0.27, so a coefficient applied in the wrong place --
   ! twice, or before the Z-vector rather than after -- moves this case and that
   ! one by different amounts.
   call check_case("H2O distorted / sto-3g, B2GP-PLYP", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.31_dp, -1.52_dp, 1.03_dp, &
                            -0.09_dp, 1.35_dp, 1.21_dp], [N_DIM, 3]), &
                   "sto-3g", 10, "b2gp-plyp", n_bad)

   ! Linear, with a doubly degenerate pi HOMO, and the case that exposed the
   ! density tolerance below.
   call check_case("HCN / sto-3g, B2GP-PLYP", [1, 6, 7], ["H", "C", "N"], &
                   reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                   "sto-3g", 14, "b2gp-plyp", n_bad)

   ! A third coefficient, and a molecule whose heavy atom is not oxygen.
   call check_case("NH3 pyramidal / sto-3g, mPW2-PLYP", [7, 1, 1, 1], &
                   ["N", "H", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            1.77_dp, 0.11_dp, -0.51_dp, &
                            -0.83_dp, 1.58_dp, -0.47_dp, &
                            -0.91_dp, -1.61_dp, -0.44_dp], [N_DIM, 4]), &
                   "sto-3g", 10, "mpw2plyp", n_bad)

   write (*, "(a)") ""
   if (n_bad == 0) then
      write (*, "(a)") "all double hybrid gradient checks passed"
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

      type(error_t) :: error
      real(dp), allocatable :: analytic(:, :), numeric(:, :)
      real(dp), allocatable :: shifted(:, :)
      real(dp) :: coarse, fine, e0, worst
      integer :: natm, ia, comp

      natm = size(numbers)
      write (*, "(a)") ""
      write (*, "(a,a)") "== ", label
      flush (output_unit)

      call dh_gradient(numbers, symbols, coords, basis, nelec, functional, &
                       e0, analytic, error)
      if (fail(error, n_bad)) return

      allocate (numeric(3, natm), shifted(N_DIM, natm))
      do ia = 1, natm
         do comp = 1, 3
            call central(numbers, symbols, coords, basis, nelec, functional, &
                         ia, comp, STEP, coarse, error)
            if (fail(error, n_bad)) return
            call central(numbers, symbols, coords, basis, nelec, functional, &
                         ia, comp, 0.5_dp*STEP, fine, error)
            if (fail(error, n_bad)) return
            ! Richardson: the h^2 error cancels between the two, leaving h^4.
            numeric(comp, ia) = (4.0_dp*fine - coarse)/3.0_dp
         end do
      end do

      write (*, "(a,f20.12)") "   E(double hybrid) = ", e0
      write (*, "(a)") "   atom      analytic (x, y, z)                        "// &
         "finite difference"
      do ia = 1, natm
         write (*, "(i6,3f14.9,a,3f14.9)") ia, analytic(:, ia), "   ", numeric(:, ia)
      end do
      worst = maxval(abs(analytic - numeric))
      write (*, "(a,es14.4)") "   largest deviation: ", worst
      write (*, "(a,es14.4)") "   |sum over atoms|:  ", maxval(abs(sum(analytic, dim=2)))
      flush (output_unit)

      if (worst > TOL) then
         write (*, "(a)") "  FAIL: analytic and finite difference disagree"
         n_bad = n_bad + 1
      end if
      if (maxval(abs(sum(analytic, dim=2))) > 1.0e-7_dp) then
         write (*, "(a)") "  FAIL: does not sum to zero over atoms"
         n_bad = n_bad + 1
      end if
   end subroutine check_case

   subroutine converge(mol, nelec, xc, scf, error)
      !! One SCF, from superposed atomic densities
      !!
      !! **The guess is load-bearing here, which is not obvious.** A converged
      !! SCF is a converged SCF whatever it started from, so a harness has no
      !! reason to care -- until a displaced geometry fails to converge at all
      !! and the finite difference silently measures the iteration instead of
      !! the energy. HCN/STO-3G is that case: nine iterations from SAD, between
      !! 68 and 130 from Wolfsberg-Helmholz, and at some displacements it never
      !! arrives. This is the guess the driver uses, so it is the one a check of
      !! the driver's arithmetic should use too.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: nelec
      type(xc_context_t), intent(inout) :: xc
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: guess_a(:, :), guess_b(:, :)

      call build_atomic_guess(mol, SCF_GUESS_SAD, guess_a, guess_b, error)
      if (error%has_error()) return
      call run_libcint_rhf(mol, nelec, 200, 1.0e-12_dp, DENSITY_TOL, .false., scf, &
                           error, xc=xc, guess=SCF_GUESS_SAD, &
                           guess_density=guess_a + guess_b)
   end subroutine converge

   subroutine central(numbers, symbols, coords, basis, nelec, functional, &
                      ia, comp, h, deriv, error)
      !! One central difference of the double hybrid energy
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec
      character(len=*), intent(in) :: functional
      integer, intent(in) :: ia, comp
      real(dp), intent(in) :: h
      real(dp), intent(out) :: deriv
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: shifted(:, :)
      real(dp) :: plus, minus

      deriv = 0.0_dp
      allocate (shifted(size(coords, 1), size(coords, 2)))
      shifted = coords
      shifted(comp, ia) = coords(comp, ia) + h
      call dh_energy(numbers, symbols, shifted, basis, nelec, functional, plus, error)
      if (error%has_error()) return
      shifted(comp, ia) = coords(comp, ia) - h
      call dh_energy(numbers, symbols, shifted, basis, nelec, functional, minus, error)
      if (error%has_error()) return
      deriv = (plus - minus)/(2.0_dp*h)
   end subroutine central

   subroutine dh_gradient(numbers, symbols, coords, basis, nelec, functional, &
                          energy, gradient, error)
      !! The total double hybrid gradient: Kohn-Sham plus scaled correlation
      !!
      !! The two halves are separate calls on purpose. The Kohn-Sham part is
      !! variational and its gradient is the ordinary one; the perturbative part
      !! is not, and its gradient carries the Z-vector. Adding a *total* MP2
      !! gradient here instead would add a Hartree-Fock reference gradient to a
      !! Kohn-Sham one.
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec
      character(len=*), intent(in) :: functional
      real(dp), intent(out) :: energy
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: xc
      type(mp2_result_t) :: mp2
      real(dp), allocatable :: corr(:, :)

      call build_libcint_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) return
      call xc_context_create(mol, functional, xc, error, polarized=.false.)
      if (error%has_error()) return
      call converge(mol, nelec, xc, scf, error)
      if (error%has_error()) return

      call run_libcint_mp2(mol, scf%orbitals, scf%orbital_energies, nelec/2, &
                           scf%energy, mp2, error)
      if (error%has_error()) return
      energy = scf%energy + xc%pt2_fraction*mp2%correlation

      call libcint_scf_gradient(mol, scf%density, orbitals=scf%orbitals, &
                                orbital_energies=scf%orbital_energies, &
                                n_occupied=scf%n_occupied, gradient=gradient, &
                                error=error, xc=xc)
      if (error%has_error()) return

      call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, nelec/2, &
                                corr, error, xc=xc, scf_density=scf%density, &
                                pt2_scale=xc%pt2_fraction)
      if (error%has_error()) return
      gradient = gradient + corr

      call xc%destroy()
      call mol%destroy()
   end subroutine dh_gradient

   subroutine dh_energy(numbers, symbols, coords, basis, nelec, functional, &
                        energy, error)
      !! `E_KS + c E_PT2` at one geometry, converged from scratch
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec
      character(len=*), intent(in) :: functional
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: xc
      type(mp2_result_t) :: mp2

      energy = 0.0_dp
      call build_libcint_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) return
      call xc_context_create(mol, functional, xc, error, polarized=.false.)
      if (error%has_error()) return
      call converge(mol, nelec, xc, scf, error)
      if (error%has_error()) return
      ! A displaced point that did not converge would enter the finite
      ! difference as whatever the iteration happened to reach, and the
      ! comparison would then be against a quantity that is not the energy.
      if (.not. scf%converged) then
         call error%set(ERROR_VALIDATION, "a displaced SCF did not converge")
         return
      end if
      call run_libcint_mp2(mol, scf%orbitals, scf%orbital_energies, nelec/2, &
                           scf%energy, mp2, error)
      if (error%has_error()) return
      energy = scf%energy + xc%pt2_fraction*mp2%correlation

      call xc%destroy()
      call mol%destroy()
   end subroutine dh_energy

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

end program check_dh_gradient
