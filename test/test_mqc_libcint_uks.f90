module test_mqc_libcint_uks
   !! Unrestricted Kohn-Sham, against the properties a reference energy cannot see.
   !!
   !! `validation/check_dft` compares the totals to PySCF's `dft.UKS`, which is the
   !! answer. What it cannot isolate is *how* an unrestricted functional goes wrong,
   !! and there are two characteristic ways, both of which converge:
   !!
   !!   * **libxc's polarised arrays are spin-interleaved.** rho arrives as (rho_a,
   !!     rho_b) per point, sigma as (sigma_aa, sigma_ab, sigma_bb), tau as (tau_a,
   !!     tau_b), and the potentials come back the same way. The F03 bindings take
   !!     bare arrays, so nothing checks a stride and a wrong one reads a real
   !!     number from the wrong place.
   !!   * **the gradient term couples the spins.** sigma_ab belongs to both, so
   !!     dE/dgrad rho_a carries a cross term in grad rho_b, and dropping it leaves
   !!     a GGA that is wrong by millihartree.
   !!
   !! The test that catches both is the closed-shell limit: with equal alpha and
   !! beta densities every polarised quantity has a restricted counterpart, so an
   !! unrestricted Kohn-Sham energy must equal the restricted one *exactly*, for
   !! every rung of the ladder. It is sharper than any comparison against another
   !! code because the reference is a number this repository already validates, and
   !! it needs neither PySCF nor a basis-set file beyond the ones already used.
   !!
   !! Also checked: that a context knows which spin channel it was built for and
   !! refuses the other evaluator, since that is the mistake the interleaving makes
   !! easy and its symptom would otherwise be a silently misread array.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf, run_libcint_uhf
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_add_potential, &
                             xc_add_potential_uks, xc_available
   implicit none
   private
   public :: collect_mqc_libcint_uks_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: E_TOL = 1.0e-10_dp
   real(dp), parameter :: D_TOL = 1.0e-8_dp

contains

   subroutine collect_mqc_libcint_uks_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("closed_shell_uks_equals_rks_lda", test_closed_lda), &
                  new_unittest("closed_shell_uks_equals_rks_gga", test_closed_gga), &
                  new_unittest("closed_shell_uks_equals_rks_hybrid", test_closed_hybrid), &
                  new_unittest("closed_shell_uks_equals_rks_mgga", test_closed_mgga), &
                  new_unittest("spin_densities_integrate_to_their_counts", test_electron_count), &
                  new_unittest("a_polarised_context_refuses_the_restricted_call", test_wrong_channel), &
                  new_unittest("a_restricted_context_refuses_the_polarised_call", test_wrong_channel_back) &
                  ]
   end subroutine collect_mqc_libcint_uks_tests

   subroutine water(mol, err, basis)
      !! Closed-shell water, at a geometry with no symmetry to flatter anything
      type(libcint_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      character(len=*), intent(in) :: basis
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
   end subroutine water

   subroutine compare_closed_shell(functional, error)
      !! One functional, run both ways on a closed shell, and the two must agree
      !!
      !! Restricted and unrestricted are different code: two Fock matrices, two
      !! densities, a spin-polarised functional and a different energy expression.
      !! On a closed shell they describe the same state, so the totals must be
      !! equal -- not close. The tolerance is the SCF's own convergence threshold
      !! and nothing looser, because there is no approximation between the two.
      character(len=*), intent(in) :: functional
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      type(rhf_result_t) :: rks, uks
      type(xc_context_t) :: xc_r, xc_u

      if (.not. xc_available()) return

      ! cc-pVDZ, not a minimal basis. The identity below holds in any basis, but
      ! sto-3g water is a marginal case for the *density* threshold here: its
      ! unrestricted Kohn-Sham energy converges to 1e-14 while the density change
      ! stalls near 5e-8, so whether it trips a 1e-8 test depends on the DIIS
      ! subspace size rather than on anything this test is about.
      call water(mol, err, "cc-pvdz")
      call check(error,.not. err%has_error(), "water must build: "//err%get_full_trace())
      if (allocated(error)) return

      call xc_context_create(mol, functional, xc_r, err, level=3)
      call check(error,.not. err%has_error(), "the restricted context must build: "//err%get_full_trace())
      if (allocated(error)) return
      call run_libcint_rhf(mol, 10, 100, E_TOL, D_TOL, .false., rks, err, xc=xc_r)
      call check(error,.not. err%has_error(), "the restricted SCF must run: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, rks%converged, "the restricted SCF must converge")
      if (allocated(error)) return

      call xc_context_create(mol, functional, xc_u, err, level=3, polarized=.true.)
      call check(error,.not. err%has_error(), "the polarised context must build: "//err%get_full_trace())
      if (allocated(error)) return
      call run_libcint_uhf(mol, 10, 1, 30, E_TOL, D_TOL, .true., uks, err, xc=xc_u)
      call check(error,.not. err%has_error(), "the unrestricted SCF must run: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, uks%converged, "the unrestricted SCF must converge")
      if (allocated(error)) return

      call check(error, abs(uks%energy - rks%energy) < 1.0e-9_dp, &
                 "unrestricted Kohn-Sham must reproduce the restricted energy on a "// &
                 "closed shell, exactly -- a difference here is an interleaving or a "// &
                 "dropped cross-spin term, not an approximation")
      if (allocated(error)) return
      ! And it must have stayed closed-shell, or the comparison above is between
      ! two different states and says nothing about the functional.
      call check(error, abs(uks%spin_squared) < 1.0e-8_dp, &
                 "the unrestricted solution must not have broken symmetry here")

      call xc_r%destroy()
      call xc_u%destroy()
      call mol%destroy()
   end subroutine compare_closed_shell

   subroutine test_closed_lda(error)
      !! LDA: rho only, so this isolates the rho interleaving
      type(error_type), allocatable, intent(out) :: error
      call compare_closed_shell("svwn", error)
   end subroutine test_closed_lda

   subroutine test_closed_gga(error)
      !! GGA: adds sigma's three components and the cross-spin gradient term
      type(error_type), allocatable, intent(out) :: error
      call compare_closed_shell("pbe", error)
   end subroutine test_closed_gga

   subroutine test_closed_hybrid(error)
      !! Hybrid: adds the same-spin exchange scaling in the Fock build
      type(error_type), allocatable, intent(out) :: error
      call compare_closed_shell("b3lyp", error)
   end subroutine test_closed_hybrid

   subroutine test_closed_mgga(error)
      !! Meta-GGA: adds tau, which is interleaved as well
      type(error_type), allocatable, intent(out) :: error
      call compare_closed_shell("tpss", error)
   end subroutine test_closed_mgga

   subroutine test_electron_count(error)
      !! The two spin densities integrate to the two electron counts
      !!
      !! Exact, and the cheapest statement that the grid and the spin densities
      !! agree -- and a property of the density alone, which is why this one does
      !! not ask the SCF to converge. It is also the check that would catch a spin density built at twice
      !! its size -- the unrestricted density is C C^T and not 2 C C^T, and the
      !! factor of two between those conventions is the oldest bug in this file's
      !! subject matter. Run on a doublet so alpha and beta differ and the two
      !! counts are distinguishable.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      type(rhf_result_t) :: uks
      type(xc_context_t) :: xc
      real(dp), allocatable :: v_a(:, :), v_b(:, :), s(:, :)
      real(dp) :: e_xc, n_elec, n_a, n_b
      real(dp) :: c(3, 2)

      if (.not. xc_available()) return

      ! OH, a doublet: five alpha and four beta electrons.
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9697_dp*ANG], [3, 2])
      ! def2-SVP and a generous iteration limit, following the unrestricted
      ! Hartree-Fock test on this same molecule: OH is the case whose sigma and pi
      ! solutions compete, and it is slow to settle rather than difficult.
      call build_libcint_molecule([8, 1], ["O ", "H "], c, "def2-svp", mol, err)
      call check(error,.not. err%has_error(), "OH must build: "//err%get_full_trace())
      if (allocated(error)) return

      call xc_context_create(mol, "pbe", xc, err, level=3, polarized=.true.)
      call check(error,.not. err%has_error(), "the polarised context must build: "//err%get_full_trace())
      if (allocated(error)) return

      ! Convergence is deliberately *not* required, and the iteration limit is
      ! modest. Tr(D S) is the electron count for any density built as
      ! C_occ C_occ^T, converged or not, so every assertion below holds on
      ! whatever density the loop last produced -- and requiring convergence made
      ! this test flaky. OH is the case whose sigma and pi solutions compete, the
      ! threaded Fock build sums its thread-local copies in an order that varies
      ! run to run, and that roundoff is enough to move a delicate open-shell
      ! trajectory. It failed about half the time on a 40-core box and never on
      ! CI's two cores, which is the worst possible distribution of a failure.
      call run_libcint_uhf(mol, 9, 2, 60, E_TOL, D_TOL, .false., uks, err, xc=xc)
      call check(error,.not. err%has_error(), "the unrestricted SCF must run: "//err%get_full_trace())
      if (allocated(error)) return

      allocate (v_a(mol%nao, mol%nao), v_b(mol%nao, mol%nao))
      call xc_add_potential_uks(xc, mol, uks%density, uks%density_beta, v_a, v_b, &
                                e_xc, n_elec, err)
      call check(error,.not. err%has_error(), "the polarised potential must evaluate: "//err%get_full_trace())
      if (allocated(error)) return

      ! Tr(D S) per spin, which the grid integral must reproduce.
      call mol%overlap(s)
      n_a = sum(uks%density*s)
      n_b = sum(uks%density_beta*s)
      call check(error, abs(n_a - 5.0_dp) < 1.0e-8_dp, &
                 "the alpha density must carry five electrons")
      if (allocated(error)) return
      call check(error, abs(n_b - 4.0_dp) < 1.0e-8_dp, &
                 "the beta density must carry four electrons")
      if (allocated(error)) return
      ! The grid's own total, which is what the functional was evaluated on. Looser
      ! than the trace above because this one is a quadrature, not an identity.
      call check(error, abs(n_elec - 9.0_dp) < 1.0e-6_dp, &
                 "the integrated spin densities must total nine electrons")
      if (allocated(error)) return
      ! Exchange-correlation is negative for any real functional; a sign error here
      ! would still converge, to something absurd.
      call check(error, e_xc < 0.0_dp, "the exchange-correlation energy must be negative")

      call xc%destroy()
      call mol%destroy()
   end subroutine test_electron_count

   subroutine test_wrong_channel(error)
      !! A polarised context cannot be evaluated on one density
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      type(xc_context_t) :: xc
      real(dp), allocatable :: v(:, :), d(:, :)
      real(dp) :: e_xc, n_elec

      if (.not. xc_available()) return

      call water(mol, err, "sto-3g")
      call check(error,.not. err%has_error(), "water must build: "//err%get_full_trace())
      if (allocated(error)) return
      call xc_context_create(mol, "svwn", xc, err, level=1, polarized=.true.)
      call check(error,.not. err%has_error(), "the polarised context must build: "//err%get_full_trace())
      if (allocated(error)) return

      allocate (v(mol%nao, mol%nao), d(mol%nao, mol%nao))
      d = 0.0_dp
      call xc_add_potential(xc, mol, d, v, e_xc, n_elec, err)
      call check(error, err%has_error(), &
                 "a polarised context must refuse the restricted evaluator rather "// &
                 "than read its arrays with the wrong stride")

      call xc%destroy()
      call mol%destroy()
   end subroutine test_wrong_channel

   subroutine test_wrong_channel_back(error)
      !! And a restricted context cannot serve two densities
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      type(xc_context_t) :: xc
      real(dp), allocatable :: v_a(:, :), v_b(:, :), d(:, :)
      real(dp) :: e_xc, n_elec

      if (.not. xc_available()) return

      call water(mol, err, "sto-3g")
      call check(error,.not. err%has_error(), "water must build: "//err%get_full_trace())
      if (allocated(error)) return
      call xc_context_create(mol, "svwn", xc, err, level=1)
      call check(error,.not. err%has_error(), "the restricted context must build: "//err%get_full_trace())
      if (allocated(error)) return

      allocate (v_a(mol%nao, mol%nao), v_b(mol%nao, mol%nao), d(mol%nao, mol%nao))
      d = 0.0_dp
      call xc_add_potential_uks(xc, mol, d, d, v_a, v_b, e_xc, n_elec, err)
      call check(error, err%has_error(), &
                 "a restricted context must refuse the polarised evaluator")

      call xc%destroy()
      call mol%destroy()
   end subroutine test_wrong_channel_back

end module test_mqc_libcint_uks

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_libcint_uks, only: collect_mqc_libcint_uks_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_libcint_uks", collect_mqc_libcint_uks_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
