!! The unrestricted RI-MP2 gradient, against differences of its own energy
module test_mqc_czt_uri_mp2_gradient
   !! The reference is the energy this gradient claims to differentiate:
   !! `run_czt_uhf` followed by `run_czt_uri_mp2`, at geometries displaced
   !! along each Cartesian coordinate, by central differences. That energy is
   !! held to PySCF's elsewhere, so agreement here is agreement with it.
   !!
   !! Three systems, chosen for what each one reaches:
   !!
   !!   * the water cation, a doublet with no symmetry to hide a sign behind,
   !!     all electrons correlated;
   !!   * the same with a frozen core, which adds the frozen-correlated
   !!     rotations and their part of the right-hand side;
   !!   * triplet water, where the two spins differ by two orbitals and the
   !!     opposite-spin block is furthest from square.
   !!
   !! Beside the differences: the correlation energy recomputed from the
   !! three-index density must be the energy routine's, the gradient must sum
   !! to zero over the atoms, and a closed shell must give the restricted
   !! gradient.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, run_czt_uhf
   use mqc_czt_mp2, only: mp2_result_t, run_czt_uri_mp2
   use mqc_czt_ri_mp2_gradient, only: czt_ri_mp2_gradient
   use mqc_czt_uri_mp2_gradient, only: czt_uri_mp2_gradient
   implicit none
   private

   public :: collect_mqc_czt_uri_mp2_gradient_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: STEP = 1.0e-3_dp
      !! Bohr. Central differences at this step are 3e-7 off here and at half
      !! of it 8e-8, a clean factor of four: the error is the differences',
      !! and the two are extrapolated to remove it.
   real(dp), parameter :: FD_TOL = 2.0e-8_dp
   character(len=*), parameter :: BASIS = "6-31g"
   character(len=*), parameter :: AUX_BASIS = "cc-pvdz-rifit"

contains

   subroutine collect_mqc_czt_uri_mp2_gradient_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("doublet_matches_finite_differences", test_doublet), &
                  new_unittest("frozen_core_matches_finite_differences", test_frozen), &
                  new_unittest("triplet_matches_finite_differences", test_triplet), &
                  new_unittest("closed_shell_is_the_restricted_gradient", test_closed_shell) &
                  ]
   end subroutine collect_mqc_czt_uri_mp2_gradient_tests

   function water_coords() result(c)
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
   end function water_coords

   subroutine energy_at(coords, nelec, mult, frozen, energy, err, scf_out, mol_out, aux_out)
      !! UHF + RI-MP2 at one geometry, converged tightly enough to difference
      real(dp), intent(in) :: coords(3, 3)
      integer, intent(in) :: nelec, mult, frozen
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: err
      type(rhf_result_t), intent(out), optional :: scf_out
      type(czt_molecule_t), intent(out), optional :: mol_out, aux_out

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(mp2_result_t) :: mp2

      energy = 0.0_dp
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], coords, BASIS, mol, err)
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], coords, AUX_BASIS, aux, err)
      if (err%has_error()) return
      call run_czt_uhf(mol, nelec, mult, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err, &
                       grad_tol=1.0e-10_dp)
      if (err%has_error()) return
      call run_czt_uri_mp2(mol, aux, scf%orbitals, scf%orbitals_beta, &
                           scf%orbital_energies, scf%orbital_energies_beta, &
                           scf%n_occupied, scf%n_occupied_beta, scf%energy, mp2, err, &
                           n_frozen=frozen)
      energy = mp2%total
      if (present(scf_out)) scf_out = scf
      if (present(mol_out)) mol_out = mol
      if (present(aux_out)) aux_out = aux
   end subroutine energy_at

   subroutine compare(error, nelec, mult, frozen, label)
      !! The analytic gradient of one state against central differences
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: nelec, mult, frozen
      character(len=*), intent(in) :: label

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(mp2_result_t) :: mp2
      type(error_t) :: err
      real(dp), allocatable :: gradient(:, :)
      real(dp) :: coords(3, 3), moved(3, 3), numeric(3, 3), e0, ep, em, e_corr
      real(dp) :: coarse, fine
      integer :: iatom, comp

      coords = water_coords()
      call energy_at(coords, nelec, mult, frozen, e0, err, scf, mol, aux)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 label//": the reference must converge: "//err%get_full_trace())
      if (allocated(error)) return
      call run_czt_uri_mp2(mol, aux, scf%orbitals, scf%orbitals_beta, &
                           scf%orbital_energies, scf%orbital_energies_beta, &
                           scf%n_occupied, scf%n_occupied_beta, scf%energy, mp2, err, &
                           n_frozen=frozen)
      call czt_uri_mp2_gradient(mol, aux, scf%orbitals, scf%orbitals_beta, &
                                scf%orbital_energies, scf%orbital_energies_beta, &
                                scf%n_occupied, scf%n_occupied_beta, gradient, err, &
                                n_frozen=frozen, e_corr=e_corr)
      call check(error,.not. err%has_error(), label//": "//err%get_full_trace())
      if (allocated(error)) return

      ! The three-index density is built from the same amplitudes the energy
      ! is, so contracting it back with the fitted integrals has to return
      ! the energy routine's number.
      call check(error, abs(e_corr - mp2%correlation) < 1.0e-10_dp, &
                 label//": the three-index density does not give back the "// &
                 "correlation energy")
      if (allocated(error)) return

      do iatom = 1, 3
         do comp = 1, 3
            moved = coords
            moved(comp, iatom) = coords(comp, iatom) + STEP
            call energy_at(moved, nelec, mult, frozen, ep, err)
            moved(comp, iatom) = coords(comp, iatom) - STEP
            call energy_at(moved, nelec, mult, frozen, em, err)
            coarse = (ep - em)/(2.0_dp*STEP)
            moved(comp, iatom) = coords(comp, iatom) + 0.5_dp*STEP
            call energy_at(moved, nelec, mult, frozen, ep, err)
            moved(comp, iatom) = coords(comp, iatom) - 0.5_dp*STEP
            call energy_at(moved, nelec, mult, frozen, em, err)
            fine = (ep - em)/STEP
            ! Richardson: the leading error is quadratic in the step.
            numeric(comp, iatom) = (4.0_dp*fine - coarse)/3.0_dp
         end do
      end do
      call check(error,.not. err%has_error(), label//": a displaced energy failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      write (*, "(4x,a,a)") label, ": analytic, then analytic minus numeric"
      do iatom = 1, 3
         write (*, "(6x,3f14.9,3x,3es10.2)") gradient(:, iatom), &
            gradient(:, iatom) - numeric(:, iatom)
      end do
      call check(error, maxval(abs(gradient - numeric)) < FD_TOL, &
                 label//": the gradient does not match differences of the energy")
      if (allocated(error)) return
      call check(error, maxval(abs(sum(gradient, 2))) < 1.0e-8_dp, &
                 label//": the gradient does not sum to zero over the atoms")
      call mol%destroy()
      call aux%destroy()
   end subroutine compare

   subroutine test_doublet(error)
      type(error_type), allocatable, intent(out) :: error

      call compare(error, 9, 2, 0, "water cation")
   end subroutine test_doublet

   subroutine test_frozen(error)
      type(error_type), allocatable, intent(out) :: error

      call compare(error, 9, 2, 1, "water cation, frozen core")
   end subroutine test_frozen

   subroutine test_triplet(error)
      type(error_type), allocatable, intent(out) :: error

      call compare(error, 10, 3, 1, "triplet water, frozen core")
   end subroutine test_triplet

   subroutine test_closed_shell(error)
      !! With equal spins every term here is the restricted one counted per
      !! spin, so the two gradients must agree to the solvers' tolerances.
      !! The restricted routine is asked for its exact response, which is
      !! what this one always uses.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: r
      type(error_t) :: err
      real(dp), allocatable :: restricted(:, :), unrestricted(:, :)
      real(dp) :: coords(3, 3)

      coords = water_coords()
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], coords, BASIS, mol, err)
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], coords, AUX_BASIS, aux, err)
      call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., r, err, &
                       grad_tol=1.0e-10_dp)
      call check(error,.not. err%has_error() .and. r%converged, &
                 "water must converge: "//err%get_full_trace())
      if (allocated(error)) return

      call czt_ri_mp2_gradient(mol, aux, r%orbitals, r%orbital_energies, r%n_occupied, &
                               restricted, err, n_frozen=1, fit_response=.false.)
      call czt_uri_mp2_gradient(mol, aux, r%orbitals, r%orbitals, r%orbital_energies, &
                                r%orbital_energies, r%n_occupied, r%n_occupied, &
                                unrestricted, err, n_frozen=1)
      call check(error,.not. err%has_error(), "gradients: "//err%get_full_trace())
      if (allocated(error)) return
      write (*, "(4x,a,es10.2)") "closed shell: largest difference from restricted ", &
         maxval(abs(unrestricted - restricted))
      call check(error, maxval(abs(unrestricted - restricted)) < 1.0e-8_dp, &
                 "a closed shell's unrestricted gradient is not the restricted one")
      call mol%destroy()
      call aux%destroy()
   end subroutine test_closed_shell

end module test_mqc_czt_uri_mp2_gradient

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_uri_mp2_gradient, only: collect_mqc_czt_uri_mp2_gradient_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_uri_mp2_gradient", &
                               collect_mqc_czt_uri_mp2_gradient_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
