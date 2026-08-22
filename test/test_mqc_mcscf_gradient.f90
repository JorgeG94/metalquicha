module test_mqc_mcscf_gradient
   !! The analytic MCSCF nuclear gradient
   !!
   !! Two of these are structural rather than numerical, and they are the ones
   !! that catch the mistakes this gradient is actually prone to. Emptying the
   !! active space has to reproduce the closed-shell gradient *exactly*, because
   !! with no active orbitals every term this module adds on top of the
   !! separable one is identically zero -- so a disagreement there is a bug in
   !! the part shared with Hartree-Fock, isolated from the active-space part
   !! that would otherwise mask it. And the cumulant has to vanish for a
   !! determinant, because a single determinant's two-particle density *is* the
   !! mean-field product that the separable term already counted; if it does not
   !! vanish the active energy is being counted twice, which is invisible at
   !! `n_active = 0` and wrong everywhere else.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_gradient, only: libcint_scf_gradient
   use mqc_libcint_mcscf, only: run_libcint_casscf, casscf_result_t
   use mqc_libcint_mcscf_gradient, only: libcint_mcscf_gradient, &
                                         cumulant_two_particle_density
   implicit none
   private

   public :: collect_mcscf_gradient_tests

   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.190459_dp, &
                           0.0_dp, 1.459853_dp, -0.884001_dp, &
                           0.0_dp, -1.459853_dp, -0.884001_dp], [3, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]

contains

   subroutine collect_mcscf_gradient_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("an empty active space is the SCF gradient", &
                               empty_active_space_is_scf), &
                  new_unittest("the cumulant vanishes for a determinant", &
                               cumulant_vanishes_for_a_determinant), &
                  new_unittest("the gradient differences the energy", &
                               gradient_differences_the_energy), &
                  new_unittest("water 6-31g CAS(4,4) agrees with PySCF", &
                               water_631g_matches_pyscf), &
                  new_unittest("water 6-31g CAS(6,5) agrees with PySCF", &
                               water_631g_65_matches_pyscf), &
                  new_unittest("an ORMAS gradient differences its energy", &
                               ormas_gradient_differences_the_energy) &
                  ]
   end subroutine collect_mcscf_gradient_tests

   subroutine empty_active_space_is_scf(error)
      !! With nothing active, this is the closed-shell gradient
      !!
      !! Not to machine precision, and the reason is worth stating rather than
      !! hiding in a loose tolerance. The two build the Pulay term from
      !! different objects: the SCF weights the orbitals by their converged
      !! energies, while this weights them by a generalised Fock matrix rebuilt
      !! from the density through a Schwarz-screened direct build. Those agree
      !! to the screening, which is around 1e-11 here on a gradient of 7e-2.
      !!
      !! The bound is still tight enough to be worth something: a factor of two
      !! or a wrong sign on any term -- which is what this test exists to catch
      !! -- lands at 1e-2, nine orders above where it passes.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: g_scf(:, :), g_mc(:, :)
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      real(dp) :: worst
      integer :: i, j

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return

      call run_libcint_rhf(mol, 10, 300, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error, scf%converged, "the SCF should converge")
      if (allocated(error)) return

      call libcint_scf_gradient(mol, scf%density, orbitals=scf%orbitals, &
                                orbital_energies=scf%orbital_energies, &
                                n_occupied=scf%n_occupied, gradient=g_scf, error=err)
      call check(error,.not. err%has_error(), "the SCF gradient should build")
      if (allocated(error)) return

      ! No active orbitals: the densities are zero-sized, and every active term
      ! is skipped rather than being multiplied by zero.
      allocate (dm1(0, 0), dm2(0, 0, 0, 0))
      call libcint_mcscf_gradient(mol, scf%orbitals, scf%n_occupied, 0, dm1, dm2, &
                                  g_mc, err)
      call check(error,.not. err%has_error(), "the MCSCF gradient should build")
      if (allocated(error)) return

      worst = 0.0_dp
      do j = 1, size(g_scf, 2)
         do i = 1, 3
            worst = max(worst, abs(g_scf(i, j) - g_mc(i, j)))
         end do
      end do
      call check(error, worst < 1.0e-9_dp, &
                 "an empty active space should reproduce the SCF gradient")
   end subroutine empty_active_space_is_scf

   subroutine cumulant_vanishes_for_a_determinant(error)
      !! A closed-shell determinant's two-particle density is exactly the
      !! mean-field product the separable term already carries.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 4
      real(dp) :: dm1(N, N), dm2(N, N, N, N)
      real(dp), allocatable :: ddm2(:, :, :, :)
      integer :: t, u, v, w
      real(dp) :: worst

      ! Every active orbital doubly occupied: `dm1 = 2 delta`, and the
      ! two-particle density of that determinant follows from Wick's theorem.
      dm1 = 0.0_dp
      do t = 1, N
         dm1(t, t) = 2.0_dp
      end do
      do w = 1, N
         do v = 1, N
            do u = 1, N
               do t = 1, N
                  dm2(t, u, v, w) = dm1(t, u)*dm1(v, w) - 0.5_dp*dm1(t, v)*dm1(u, w)
               end do
            end do
         end do
      end do

      call cumulant_two_particle_density(dm1, dm2, ddm2)
      worst = maxval(abs(ddm2))
      call check(error, worst < 1.0e-14_dp, &
                 "the cumulant should vanish for a single determinant")
   end subroutine cumulant_vanishes_for_a_determinant

   subroutine gradient_differences_the_energy(error)
      !! Against this program's own energy, with no reference code involved
      !!
      !! The PySCF comparisons check the whole contraction against another
      !! implementation, which is the stronger test of correctness but shares
      !! a blind spot with it: both codes have to agree on *which* stationary
      !! point they converged to, and on a surface with more than one they need
      !! not. Water in STO-3G is exactly that case -- this optimiser finds a
      !! point about 1.7 mHartree below the one PySCF reaches from Hartree-Fock
      !! orbitals -- and a gradient compared across two different solutions
      !! disagrees for a reason that has nothing to do with the gradient.
      !!
      !! Differencing the energy here sidesteps that entirely: both numbers come
      !! from the same optimiser landing on the same point. It is the check that
      !! the analytic expression really is the derivative of the energy this
      !! program computes, rather than of the one another program computes.
      !!
      !! **The tolerance measures the difference formula, not the gradient.**
      !! Central differences carry an `O(h^2)` truncation term, and a step study
      !! on the largest component showed it behaving exactly as it should:
      !!
      !! ```
      !!   h        error      ratio
      !!   8e-3     9.75e-6      -
      !!   4e-3     2.44e-6    4.00
      !!   2e-3     6.09e-7    4.00
      !!   1e-3     1.52e-7    4.00
      !!   5e-4     3.82e-8    3.98
      !! ```
      !!
      !! Four consecutive factors of four: the residual is the formula's, and
      !! the analytic gradient is the limit it is converging to. So the bound
      !! below is set by `h`, and tightening it without also shrinking `h` would
      !! be testing arithmetic that is not in question. `h` is not shrunk
      !! further because the noise floor grows as `1/h` once the energy's own
      !! convergence starts to show.
      !!
      !! Every component is differenced rather than one, which costs eighteen
      !! optimisations and catches a term that is wrong only off the symmetry
      !! axis.
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: STEP = 2.0e-3_dp
      real(dp) :: plus, minus, numerical, analytic
      real(dp), allocatable :: gradient(:, :)
      real(dp) :: moved(3, 3)

      integer :: iatom, comp
      real(dp) :: worst

      call energy_at(error, WATER, gradient=gradient)
      if (allocated(error)) return

      worst = 0.0_dp
      do iatom = 1, 3
         do comp = 1, 3
            analytic = gradient(comp, iatom)
            moved = WATER
            moved(comp, iatom) = WATER(comp, iatom) + STEP
            call energy_at(error, moved, energy=plus)
            if (allocated(error)) return
            moved(comp, iatom) = WATER(comp, iatom) - STEP
            call energy_at(error, moved, energy=minus)
            if (allocated(error)) return
            numerical = (plus - minus)/(2.0_dp*STEP)
            worst = max(worst, abs(numerical - analytic))
         end do
      end do
      call check(error, worst < 1.0e-6_dp, &
                 "the analytic gradient should difference the energy")
   end subroutine gradient_differences_the_energy

   subroutine energy_at(error, coordinates, energy, gradient)
      !! A converged CAS(4,4) at one geometry, and optionally its gradient
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: coordinates(3, 3)
      real(dp), intent(out), optional :: energy
      real(dp), allocatable, intent(out), optional :: gradient(:, :)

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casscf_result_t) :: cas
      type(error_t) :: err

      call build_libcint_molecule(WATER_Z, WATER_SYM, coordinates, "6-31g", mol, err)
      call run_libcint_rhf(mol, 10, 300, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error, scf%converged, "the SCF should converge")
      if (allocated(error)) return

      call run_libcint_casscf(mol, scf%orbitals, 3, 4, 2, 2, cas, err, &
                              max_iterations=400, gradient_tol=1.0e-9_dp, &
                              verbose=.false.)
      call check(error, cas%converged, "the CASSCF should converge")
      if (allocated(error)) return

      if (present(energy)) energy = cas%energy
      if (present(gradient)) then
         call libcint_mcscf_gradient(mol, cas%orbitals, 3, 4, cas%dm1, cas%dm2, &
                                     gradient, err)
         call check(error,.not. err%has_error(), "the gradient should build")
      end if
   end subroutine energy_at

   subroutine water_631g_matches_pyscf(error)
      !! The same, with six virtual orbitals the STO-3G case does not have
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: REFERENCE(3, 3) = reshape( &
                             [0.0_dp, 0.0_dp, -1.714654279167e-02_dp, &
                              0.0_dp, -1.454057677107e-02_dp, 8.573295541402e-03_dp, &
                              0.0_dp, 1.454064057661e-02_dp, 8.573247250258e-03_dp], [3, 3])

      call compare_against(error, "6-31g", 3, 4, 2, 2, -76.037525832320_dp, REFERENCE)
   end subroutine water_631g_matches_pyscf

   subroutine water_631g_65_matches_pyscf(error)
      !! A different split of the same basis: two inactive rather than three
      !!
      !! Worth having beside CAS(4,4) because the inactive-active cross term is
      !! the one place a wrong factor can hide behind a correct total, and the
      !! two cases weight it differently.
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: REFERENCE(3, 3) = reshape( &
                             [0.0_dp, 0.0_dp, -1.924985750758e-02_dp, &
                              0.0_dp, -1.493070902144e-02_dp, 9.624928753936e-03_dp, &
                              0.0_dp, 1.493070902109e-02_dp, 9.624928753640e-03_dp], [3, 3])

      call compare_against(error, "6-31g", 2, 5, 3, 3, -76.039192258019_dp, REFERENCE)
   end subroutine water_631g_65_matches_pyscf

   subroutine ormas_gradient_differences_the_energy(error)
      !! The same expression, on a wave function with a restricted space
      !!
      !! Nothing here is specific to a complete active space: the gradient is
      !! built from `dm1` and `dm2` and does not ask where they came from, and
      !! `run_libcint_casscf` optimises a restricted space too -- its redundancy
      !! rule already knows that rotating one active orbital into another stops
      !! being redundant once the two fall in different subspaces. So this
      !! ought to work, and "ought to" is why it is tested rather than assumed.
      !!
      !! Differenced against the energy rather than against a reference code,
      !! because PySCF has no ORMAS to compare with.
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: STEP = 2.0e-3_dp
      real(dp) :: plus, minus, numerical, analytic
      real(dp), allocatable :: gradient(:, :)
      real(dp) :: moved(3, 3)

      call ormas_at(error, WATER, gradient=gradient)
      if (allocated(error)) return
      analytic = gradient(3, 2)

      moved = WATER
      moved(3, 2) = WATER(3, 2) + STEP
      call ormas_at(error, moved, energy=plus)
      if (allocated(error)) return

      moved(3, 2) = WATER(3, 2) - STEP
      call ormas_at(error, moved, energy=minus)
      if (allocated(error)) return

      numerical = (plus - minus)/(2.0_dp*STEP)
      call check(error, abs(numerical - analytic) < 1.0e-6_dp, &
                 "the ORMAS gradient should difference the ORMAS energy")
   end subroutine ormas_gradient_differences_the_energy

   subroutine ormas_at(error, coordinates, energy, gradient)
      !! A converged restricted-space optimisation at one geometry
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: coordinates(3, 3)
      real(dp), intent(out), optional :: energy
      real(dp), allocatable, intent(out), optional :: gradient(:, :)

      ! Active orbitals 4-7 in two subspaces, with at most one electron above
      ! the first: a hard enough restriction that the inter-subspace rotations
      ! carry a gradient worth differencing.
      integer, parameter :: SUBSPACES(2) = [1, 3]
      integer, parameter :: MIN_E(2) = [3, 0], MAX_E(2) = [4, 1]
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casscf_result_t) :: cas
      type(error_t) :: err

      call build_libcint_molecule(WATER_Z, WATER_SYM, coordinates, "sto-3g", mol, err)
      call run_libcint_rhf(mol, 10, 300, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error, scf%converged, "the SCF should converge")
      if (allocated(error)) return

      call run_libcint_casscf(mol, scf%orbitals, 3, 4, 2, 2, cas, err, &
                              max_iterations=400, gradient_tol=1.0e-9_dp, &
                              verbose=.false., subspaces=SUBSPACES, &
                              min_electrons=MIN_E, max_electrons=MAX_E)
      call check(error,.not. err%has_error(), "the ORMAS-SCF should run")
      if (allocated(error)) return
      call check(error, cas%converged, "the ORMAS-SCF should converge")
      if (allocated(error)) return

      if (present(energy)) energy = cas%energy
      if (present(gradient)) then
         call libcint_mcscf_gradient(mol, cas%orbitals, 3, 4, cas%dm1, cas%dm2, &
                                     gradient, err)
         call check(error,.not. err%has_error(), "the gradient should build")
      end if
   end subroutine ormas_at

   subroutine compare_against(error, basis, n_inactive, n_active, n_alpha, n_beta, &
                              reference_energy, reference)
      !! Converge a CASSCF here and check its gradient against PySCF's
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: basis
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(in) :: reference_energy
      real(dp), intent(in) :: reference(3, 3)

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casscf_result_t) :: cas
      type(error_t) :: err
      real(dp), allocatable :: gradient(:, :)
      real(dp) :: worst, translation

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, basis, mol, err)
      call run_libcint_rhf(mol, 10, 300, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error, scf%converged, "the SCF should converge")
      if (allocated(error)) return

      call run_libcint_casscf(mol, scf%orbitals, n_inactive, n_active, n_alpha, n_beta, &
                              cas, err, max_iterations=300, gradient_tol=1.0e-8_dp, &
                              verbose=.false.)
      call check(error,.not. err%has_error(), "the CASSCF should run")
      if (allocated(error)) return
      call check(error, cas%converged, "the CASSCF should converge")
      if (allocated(error)) return
      call check(error, abs(cas%energy - reference_energy) < 1.0e-8_dp, &
                 "the CASSCF energy should match PySCF's")
      if (allocated(error)) return

      call libcint_mcscf_gradient(mol, cas%orbitals, n_inactive, n_active, &
                                  cas%dm1, cas%dm2, gradient, err)
      call check(error,.not. err%has_error(), "the gradient should build")
      if (allocated(error)) return

      ! PySCF's own macro-iteration stops around 1e-7 on these, so that is the
      ! floor for the comparison rather than anything about this contraction.
      worst = maxval(abs(gradient - reference))
      call check(error, worst < 1.0e-6_dp, "the gradient should match PySCF's")
      if (allocated(error)) return

      ! Independent of the reference: a gradient that does not sum to zero has
      ! lost a Pulay term, which a comparison against one number per atom can
      ! miss if both are wrong the same way.
      translation = maxval(abs(sum(gradient, dim=2)))
      call check(error, translation < 1.0e-8_dp, &
                 "the gradient should be translationally invariant")
   end subroutine compare_against

end module test_mqc_mcscf_gradient

program test_mcscf_gradient_driver
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mcscf_gradient, only: collect_mcscf_gradient_tests
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   testsuites = [new_testsuite("mcscf_gradient", collect_mcscf_gradient_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat)
   if (stat > 0) then
      write (error_unit, '(i0, a)') stat, " test(s) failed!"
      error stop
   end if
end program test_mcscf_gradient_driver
