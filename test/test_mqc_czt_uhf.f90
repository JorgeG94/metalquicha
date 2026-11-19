module test_mqc_czt_uhf
   !! Pins the unrestricted SCF against the things only it can be checked by.
   !!
   !! The validation manifest already compares open-shell energies to PySCF, so
   !! this suite is for the properties a reference energy cannot see:
   !!
   !!   * a closed-shell system run unrestricted must reproduce the restricted
   !!     answer exactly -- the two Fock matrices collapse onto one, and any
   !!     factor of two in the spin densities shows up here and nowhere else,
   !!     because an open-shell reference would absorb it into "some other
   !!     stationary point";
   !!   * the direct build and the in-core build must agree, which is what
   !!     separates a wrong six-update form from a wrong integral;
   !!   * <S^2> must be exact for a one-electron doublet, where there is no
   !!     contamination to hide behind;
   !!   * an electron count and multiplicity that cannot be paired must be
   !!     refused rather than silently rounded into one that can.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_cgto, only: molecular_basis_type
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, build_df_tensor
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, run_czt_uhf, &
                          assemble_fock_uhf, SCF_GUESS_CORE
   implicit none
   private
   public :: collect_mqc_czt_uhf_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: E_TOL = 1.0e-10_dp
   real(dp), parameter :: D_TOL = 1.0e-8_dp

contains

   subroutine collect_mqc_czt_uhf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("uhf_closed_shell_reproduces_rhf", test_closed_shell), &
                  new_unittest("uhf_direct_matches_in_core", test_direct_vs_in_core), &
                  new_unittest("uhf_spin_squared_is_exact_for_one_electron", test_s2), &
                  new_unittest("uhf_rejects_impossible_multiplicity", test_bad_multiplicity), &
                  new_unittest("uhf_odd_electron_count_is_allowed", test_odd_allowed), &
                  new_unittest("uhf_diis_delay_finds_the_ground_state", test_diis_delay), &
                  new_unittest("fitted_uhf_closed_shell_reproduces_fitted_rhf", &
                               test_fitted_closed_shell), &
                  new_unittest("fitted_uhf_build_matches_the_fitted_integrals", &
                               test_fitted_build), &
                  new_unittest("fitted_uhf_refuses_in_core", test_fitted_refuses_in_core) &
                  ]
   end subroutine collect_mqc_czt_uhf_tests

   subroutine hydrogen(mol, err)
      !! A single hydrogen atom: one electron, a doublet, nothing to cancel
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      real(dp) :: c(3, 1)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1])
      call build_czt_molecule([1], ["H "], c, "cc-pvdz", mol, err)
   end subroutine hydrogen

   subroutine dihydrogen(mol, err)
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      real(dp) :: c(3, 2)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.7414_dp*ANG], [3, 2])
      call build_czt_molecule([1, 1], ["H ", "H "], c, "cc-pvdz", mol, err)
   end subroutine dihydrogen

   subroutine dihydrogen_aux(aux, err)
      !! The fitting basis for `dihydrogen`, on the same two centres
      type(czt_molecule_t), intent(out) :: aux
      type(error_t), intent(inout) :: err
      real(dp) :: c(3, 2)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.7414_dp*ANG], [3, 2])
      call build_czt_molecule([1, 1], ["H ", "H "], c, "cc-pvdz-rifit", aux, err)
   end subroutine dihydrogen_aux

   subroutine test_fitted_closed_shell(error)
      !! H2 as a singlet, both fitted: unrestricted must land on restricted
      !!
      !! The fitted unrestricted build is the restricted one called three times
      !! with its factors rearranged, so the factor of two between a spin
      !! density and the total has one more place to go wrong than on the
      !! exact path. Equal spin densities are where that shows.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol, aux
      type(error_t) :: err
      type(rhf_result_t) :: r, u, exact

      call dihydrogen(mol, err)
      call dihydrogen_aux(aux, err)
      call check(error,.not. err%has_error(), "H2 must build: "//err%get_full_trace())
      if (allocated(error)) return

      call run_czt_rhf(mol, 2, 100, E_TOL, D_TOL, .false., r, err, aux=aux)
      call check(error, r%converged, "fitted RHF must converge: "//err%get_full_trace())
      if (allocated(error)) return
      call run_czt_uhf(mol, 2, 1, 100, E_TOL, D_TOL, .false., u, err, aux=aux)
      call check(error, u%converged, "fitted UHF must converge: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, abs(u%energy - r%energy) < 1.0e-11_dp, &
                 "a closed-shell fitted UHF must reproduce the fitted RHF energy")
      if (allocated(error)) return

      ! And the fit must have been used at all: the exact energy is a
      ! different number, by far more than either SCF's tolerance.
      call run_czt_uhf(mol, 2, 1, 100, E_TOL, D_TOL, .false., exact, err)
      call check(error, abs(u%energy - exact%energy) > 1.0e-8_dp, &
                 "the fitted energy must differ from the exact one")
   end subroutine test_fitted_closed_shell

   subroutine test_fitted_build(error)
      !! The fitted build against the stored-integral build fed `B B^T`
      !!
      !! With `(mu nu|la si)` replaced by `sum_P B(mu nu,P) B(la si,P)` the
      !! in-core build computes exactly what the fitted one claims to, by a
      !! four-index loop that knows nothing of pseudo-orbitals or of which
      !! spin's exchange went where. The densities are a triplet's, so the two
      !! spins differ and beta is not alpha scaled.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol, aux
      type(error_t) :: err
      type(rhf_result_t) :: u
      real(dp), allocatable :: b(:, :), eri(:, :, :, :), no_eri(:, :, :, :), no_b(:, :)
      real(dp), allocatable :: bounds(:, :), h(:, :)
      real(dp), allocatable :: fa(:, :), fb(:, :), ga(:, :), gb(:, :)
      real(dp) :: e_fit, e_ref
      integer :: n, mu, nu, la, si

      call dihydrogen(mol, err)
      call dihydrogen_aux(aux, err)
      call check(error,.not. err%has_error(), "H2 must build: "//err%get_full_trace())
      if (allocated(error)) return

      ! The triplet supplies alpha's density, two electrons in two orbitals.
      ! Its own beta density is empty, so beta's is made by hand below.
      call run_czt_uhf(mol, 2, 3, 100, E_TOL, D_TOL, .false., u, err)
      call check(error, u%converged, "triplet H2 must converge: "//err%get_full_trace())
      if (allocated(error)) return

      n = mol%nao
      call build_df_tensor(mol, aux, b, err)
      call check(error,.not. err%has_error(), "the fit must build: "//err%get_full_trace())
      if (allocated(error)) return
      call mol%core_hamiltonian(h)

      allocate (eri(n, n, n, n))
      do si = 1, n
         do la = 1, n
            do nu = 1, n
               do mu = 1, n
                  eri(mu, nu, la, si) = dot_product(b(mu + (nu - 1)*n, :), b(la + (si - 1)*n, :))
               end do
            end do
         end do
      end do

      allocate (fa(n, n), fb(n, n), ga(n, n), gb(n, n))
      ! Beta is the bonding orbital alone -- a different density from alpha's
      ! two, and not a multiple of it.
      block
         real(dp), allocatable :: d_a(:, :), d_b(:, :)
         allocate (d_a(n, n), d_b(n, n))
         d_a = u%density
         d_b = matmul(u%orbitals(:, 1:1), transpose(u%orbitals(:, 1:1)))

         call assemble_fock_uhf(mol, h, d_a, d_b, no_eri, bounds, fock_a=fa, fock_b=fb, &
                                e_elec=e_fit, error=err, bmat=b)
         call check(error,.not. err%has_error(), "fitted build: "//err%get_full_trace())
         if (allocated(error)) return
         call assemble_fock_uhf(mol, h, d_a, d_b, eri, bounds, fock_a=ga, fock_b=gb, &
                                e_elec=e_ref, error=err, bmat=no_b)
         call check(error,.not. err%has_error(), "stored build: "//err%get_full_trace())
         if (allocated(error)) return
      end block

      call check(error, maxval(abs(fa - ga)) < 1.0e-11_dp, &
                 "the fitted alpha Fock matrix must match the stored-integral one")
      if (allocated(error)) return
      call check(error, maxval(abs(fb - gb)) < 1.0e-11_dp, &
                 "the fitted beta Fock matrix must match the stored-integral one")
      if (allocated(error)) return
      call check(error, abs(e_fit - e_ref) < 1.0e-11_dp, &
                 "the two builds must give one energy")
      if (allocated(error)) return
      ! The spins must actually differ, or the check above could not tell
      ! alpha's exchange from beta's.
      call check(error, maxval(abs(fa - fb)) > 1.0e-3_dp, &
                 "the two spin Fock matrices must differ for this test to mean anything")
   end subroutine test_fitted_build

   subroutine test_fitted_refuses_in_core(error)
      !! Stored integrals and a fit are two answers to one question
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol, aux
      type(error_t) :: err, build_err
      type(rhf_result_t) :: u

      call dihydrogen(mol, build_err)
      call dihydrogen_aux(aux, build_err)
      call check(error,.not. build_err%has_error(), "H2 must build")
      if (allocated(error)) return

      call run_czt_uhf(mol, 2, 3, 100, E_TOL, D_TOL, .false., u, err, in_core=.true., aux=aux)
      call check(error, err%has_error(), "in_core beside a fit must be refused")
   end subroutine test_fitted_refuses_in_core

   subroutine test_closed_shell(error)
      !! H2 as a singlet: unrestricted must land exactly on restricted
      !!
      !! The strongest check here, because it needs no external number. With
      !! n_alpha == n_beta the two spin densities are equal and the unrestricted
      !! equations reduce to the restricted ones, so any disagreement is a
      !! defect in the reduction rather than a different physical answer.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(rhf_result_t) :: r, u

      call dihydrogen(mol, err)
      call check(error,.not. err%has_error(), "H2 must build: "//err%get_full_trace())
      if (allocated(error)) return

      call run_czt_rhf(mol, 2, 100, E_TOL, D_TOL, .false., r, err)
      call check(error, r%converged, "RHF must converge")
      if (allocated(error)) return

      call run_czt_uhf(mol, 2, 1, 100, E_TOL, D_TOL, .false., u, err)
      call check(error, u%converged, "UHF must converge")
      if (allocated(error)) return

      call check(error, abs(u%energy - r%energy) < 1.0e-12_dp, &
                 "a closed-shell UHF must reproduce the RHF energy")
      if (allocated(error)) return
      call check(error, abs(u%spin_squared) < 1.0e-12_dp, &
                 "a closed-shell UHF must have <S^2> = 0")
      if (allocated(error)) return
      call check(error, u%n_occupied, u%n_occupied_beta)
   end subroutine test_closed_shell

   subroutine test_direct_vs_in_core(error)
      !! The two Fock builds are separate code and must agree
      !!
      !! The direct build folds the eightfold symmetry into six scattered
      !! updates; the in-core one contracts a stored tensor with no symmetry at
      !! all. A degeneracy factor or a spin swapped between the Coulomb and
      !! exchange sources fails here while still converging.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(rhf_result_t) :: direct, in_core

      call hydrogen(mol, err)
      call check(error,.not. err%has_error(), "H must build: "//err%get_full_trace())
      if (allocated(error)) return

      call run_czt_uhf(mol, 1, 2, 100, E_TOL, D_TOL, .false., direct, err)
      call run_czt_uhf(mol, 1, 2, 100, E_TOL, D_TOL, .false., in_core, err, in_core=.true.)
      call check(error, direct%converged .and. in_core%converged, "both must converge")
      if (allocated(error)) return
      call check(error, abs(direct%energy - in_core%energy) < 1.0e-11_dp, &
                 "the direct and in-core unrestricted builds must agree")
   end subroutine test_direct_vs_in_core

   subroutine test_s2(error)
      !! One electron is a pure doublet, so <S^2> is 3/4 with nothing to spare
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(rhf_result_t) :: u

      call hydrogen(mol, err)
      call run_czt_uhf(mol, 1, 2, 100, E_TOL, D_TOL, .false., u, err)
      call check(error, u%converged, "H must converge")
      if (allocated(error)) return
      call check(error, abs(u%spin_squared - 0.75_dp) < 1.0e-12_dp, &
                 "a one-electron doublet has <S^2> = 3/4 exactly")
      if (allocated(error)) return
      call check(error, u%n_occupied, 1)
      if (allocated(error)) return
      call check(error, u%n_occupied_beta, 0)
   end subroutine test_s2

   subroutine test_bad_multiplicity(error)
      !! Parity has to work out, and asking for more unpaired electrons than
      !! there are electrons has to be refused
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(rhf_result_t) :: u

      call dihydrogen(mol, err)

      ! Two electrons cannot make a doublet: n_alpha - n_beta would be 1 and
      ! the two would have to sum to 2, which needs a half electron.
      call run_czt_uhf(mol, 2, 2, 100, E_TOL, D_TOL, .false., u, err)
      call check(error, err%has_error(), "two electrons cannot form a doublet")
      if (allocated(error)) return

      call err%clear()
      call run_czt_uhf(mol, 2, 5, 100, E_TOL, D_TOL, .false., u, err)
      call check(error, err%has_error(), &
                 "a multiplicity needing more unpaired electrons than exist must be refused")
      if (allocated(error)) return

      call err%clear()
      call run_czt_uhf(mol, 2, 0, 100, E_TOL, D_TOL, .false., u, err)
      call check(error, err%has_error(), "multiplicity below one must be refused")
   end subroutine test_bad_multiplicity

   subroutine test_odd_allowed(error)
      !! The restricted path refuses an odd electron count; this one must not
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(rhf_result_t) :: r, u

      call hydrogen(mol, err)

      call run_czt_rhf(mol, 1, 100, E_TOL, D_TOL, .false., r, err)
      call check(error, err%has_error(), "RHF must refuse one electron")
      if (allocated(error)) return

      call err%clear()
      call run_czt_uhf(mol, 1, 2, 100, E_TOL, D_TOL, .false., u, err)
      call check(error,.not. err%has_error(), "UHF must accept one electron: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, u%converged, "and converge")
   end subroutine test_odd_allowed

   subroutine test_diis_delay(error)
      !! OH in def2-SVP, the case that made the DIIS start delay necessary
      !!
      !! Extrapolating from the first iteration converges it to a sigma-hole
      !! state 158 mHartree above the ground state, and converges it tidily --
      !! dE below 1e-9, no warning, an <S^2> of 0.7584 that looks like an
      !! ordinary doublet. The tell is in the orbital energies: the pi pair
      !! stays degenerate to six digits, because a symmetric guess extrapolated
      !! linearly cannot leave the symmetric subspace.
      !!
      !! Both halves are asserted. The default has to find the ground state,
      !! and starting at iteration one has to still find the wrong one --
      !! otherwise this passes for some unrelated reason and stops guarding
      !! the delay it exists for.
      !!
      !! **The second half names its guess, and that is load-bearing.** The
      !! pathology needs a symmetric starting point, and the core guess is the
      !! symmetric one -- `F = H` knows nothing about the electrons, so the pi
      !! pair enters exactly degenerate and linear extrapolation cannot split
      !! it. Wolfsberg-Helmholz enters already slightly split, and from there
      !! even `diis_start=1` finds the ground state. That was discovered by
      !! changing the default guess to GWH and watching this assertion fail:
      !! both runs came back at -75.325108417965, the ground state, and the
      !! test stopped guarding anything. Pinned to `core` it keeps demonstrating
      !! the failure the delay exists to prevent.
      !!
      !! Worth stating plainly, because it is the more useful half of that
      !! discovery: a starting guess does not only change how *fast* an SCF
      !! converges, it changes *which solution* it finds.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(rhf_result_t) :: good, bad
      real(dp) :: c(3, 2)
      real(dp), parameter :: GROUND = -75.325108417978_dp

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9697_dp*ANG], [3, 2])
      call build_czt_molecule([8, 1], ["O ", "H "], c, "def2-svp", mol, err)
      call check(error,.not. err%has_error(), "OH must build: "//err%get_full_trace())
      if (allocated(error)) return

      call run_czt_uhf(mol, 9, 2, 400, E_TOL, D_TOL, .false., good, err)
      call check(error, good%converged, "the default must converge")
      if (allocated(error)) return
      call check(error, abs(good%energy - GROUND) < 1.0e-8_dp, &
                 "the default DIIS delay must reach the pi-hole ground state")
      if (allocated(error)) return

      call run_czt_uhf(mol, 9, 2, 400, E_TOL, D_TOL, .false., bad, err, &
                       diis_start=1, guess=SCF_GUESS_CORE)
      call check(error, bad%converged, "extrapolating from the start still converges")
      if (allocated(error)) return
      call check(error, bad%energy > GROUND + 0.1_dp, &
                 "extrapolating from the start must still find the higher state; if it "// &
                 "does not, this test no longer guards the delay")
   end subroutine test_diis_delay

end module test_mqc_czt_uhf

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_uhf, only: collect_mqc_czt_uhf_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_uhf", collect_mqc_czt_uhf_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
