!! Atomic charges, against what a charge is obliged to do
module test_mqc_czt_charges
   !! Neither scheme here has a reference this test could use. Mulliken charges
   !! are basis-set dependent by construction, CHELPG charges depend on the grid
   !! they were fitted on, and both are compared against other programs in
   !! `validation/`. What can be checked without any of that is the arithmetic
   !! every population analysis is obliged to satisfy, and those checks catch
   !! the failures that a number comparison would not localise:
   !!
   !!   * **The charges sum to the molecular charge.** For Mulliken because the
   !!     trace of `D S` is the electron count; for CHELPG because the fit is
   !!     solved under that constraint as a Lagrange multiplier. A dropped
   !!     off-diagonal term or a mis-assigned basis function breaks the first;
   !!     a constraint row assembled wrongly breaks the second, and both still
   !!     produce charges of a plausible size.
   !!   * **Every basis function belongs to exactly one atom, in range.** That
   !!     mapping is the whole of `ao_to_atom`, and everything downstream --
   !!     populations, distributed multipoles, the atomic blocks of the energy
   !!     decomposition -- inherits it.
   !!   * **The fitting grid lies outside the molecule.** CHELPG excludes points
   !!     inside the van der Waals envelope, because the potential there is not
   !!     what point charges are meant to reproduce. A grid that leaked inside
   !!     would still fit, to a potential nobody wants fitted.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: run_czt_rhf, run_czt_uhf, rhf_result_t
   use mqc_czt_charges, only: ao_to_atom, mulliken_charges, chelpg_charges, &
                              mulliken_spin_populations, &
                              chelpg_grid
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_czt_charges_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_czt_charges_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ao_to_atom_is_a_partition", test_ao_map), &
                  new_unittest("mulliken_sums_to_the_molecular_charge", test_mulliken_sum), &
                  new_unittest("chelpg_sums_to_the_molecular_charge", test_chelpg_sum), &
                  new_unittest("chelpg_grid_avoids_the_molecule", test_grid_excludes), &
                  new_unittest("spin_populations_sum_to_the_spin", test_spin_populations), &
                  new_unittest("closed_shell_spin_populations_vanish", test_spin_closed_shell) &
                  ]
   end subroutine collect_mqc_czt_charges_tests

   subroutine water(mol, err)
      !! One water molecule in a small basis, which is all any of these need.
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err

      real(dp) :: c(3, 3)
      integer :: z(3)
      character(len=2) :: symbols(3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule(z, symbols, c, "sto-3g", mol, err)
   end subroutine water

   subroutine converged_water(mol, scf, err)
      !! The same molecule with a converged restricted density on it.
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      call water(mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, 10, 60, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
   end subroutine converged_water

   subroutine test_ao_map(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      integer, allocatable :: owner(:)
      integer :: iatom, mu
      logical :: seen

      call water(mol, err)
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return

      call ao_to_atom(mol, owner)
      call check(error, size(owner) == mol%nao, "one owner per basis function")
      if (allocated(error)) return

      ! In range, and every atom carries something: sto-3g puts functions on
      ! hydrogen too, so an empty atom would mean the map has lost a shell.
      do mu = 1, size(owner)
         call check(error, owner(mu) >= 1 .and. owner(mu) <= mol%natm, &
                    "a basis function is owned by an atom that does not exist")
         if (allocated(error)) return
      end do
      do iatom = 1, mol%natm
         seen = any(owner == iatom)
         call check(error, seen, "an atom owns no basis functions")
         if (allocated(error)) return
      end do
   end subroutine test_ao_map

   subroutine test_mulliken_sum(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :), charges(:)

      call converged_water(mol, scf, err)
      call check(error,.not. err%has_error(), "the reference SCF")
      if (allocated(error)) return
      call check(error, scf%converged, "the reference SCF did not converge")
      if (allocated(error)) return

      call mol%overlap(overlap)
      call mulliken_charges(mol, scf%density, overlap, charges, err)
      call check(error,.not. err%has_error(), "computing Mulliken charges")
      if (allocated(error)) return

      call check(error, size(charges) == mol%natm, "one charge per atom")
      if (allocated(error)) return
      ! Neutral water: the populations must account for every electron.
      call check(error, abs(sum(charges)) < 1.0e-10_dp, &
                 "Mulliken charges do not sum to the molecular charge")
      if (allocated(error)) return
      ! And oxygen must be the negative end of it. Not a tight number -- the
      ! point is that the sign convention is charge and not population.
      call check(error, charges(1) < 0.0_dp, &
                 "oxygen came out positive; the sign convention is inverted")
   end subroutine test_mulliken_sum

   subroutine test_chelpg_sum(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: charges(:)

      call converged_water(mol, scf, err)
      call check(error,.not. err%has_error(), "the reference SCF")
      if (allocated(error)) return

      call chelpg_charges(mol, scf%density, charges, err, total_charge=0.0_dp)
      call check(error,.not. err%has_error(), "fitting CHELPG charges")
      if (allocated(error)) return

      call check(error, size(charges) == mol%natm, "one charge per atom")
      if (allocated(error)) return
      ! The constraint is solved exactly rather than approached, so this is
      ! tight on purpose: a loose sum means the multiplier row is wrong.
      call check(error, abs(sum(charges)) < 1.0e-9_dp, &
                 "CHELPG charges do not sum to the constrained total")
      if (allocated(error)) return
      call check(error, charges(1) < 0.0_dp, &
                 "oxygen came out positive; the sign convention is inverted")
      if (allocated(error)) return
      ! Symmetry the geometry does not have exactly, so this is deliberately
      ! loose: the two hydrogens are inequivalent here only through the grid.
      call check(error, abs(charges(2) - charges(3)) < 5.0e-2_dp, &
                 "the two hydrogens of water disagree by more than the grid explains")
   end subroutine test_chelpg_sum

   subroutine test_spin_populations(error)
      !! The spin partition sums to n_alpha - n_beta, not to the charge
      !!
      !! The two sum rules are what tell the routines apart, and getting them
      !! crossed is the plausible bug: a spin population built by subtracting
      !! from the nuclear charge would look reasonable per atom and sum to
      !! something with no meaning.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :), populations(:), spin_density(:, :)

      call water(mol, err)
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return

      ! The cation as a doublet: nine electrons, five alpha and four beta.
      call run_czt_uhf(mol, 9, 2, 60, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
      call check(error,.not. err%has_error(), "the reference UHF")
      if (allocated(error)) return
      call check(error, scf%converged, "the reference UHF did not converge")
      if (allocated(error)) return

      call mol%overlap(overlap)
      spin_density = scf%density - scf%density_beta
      call mulliken_spin_populations(mol, spin_density, overlap, populations, err)
      call check(error,.not. err%has_error(), "computing spin populations")
      if (allocated(error)) return

      call check(error, size(populations) == mol%natm, "one population per atom")
      if (allocated(error)) return
      call check(error, abs(sum(populations) - 1.0_dp) < 1.0e-10_dp, &
                 "spin populations must sum to n_alpha - n_beta, here one")
   end subroutine test_spin_populations

   subroutine test_spin_closed_shell(error)
      !! A closed shell has no spin density, so every population is zero
      !!
      !! Worth its own test because it is the case where a sign error or a
      !! stray nuclear charge cannot hide behind a plausible-looking number:
      !! the answer is exactly zero on every atom.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :), populations(:), spin_density(:, :)

      call converged_water(mol, scf, err)
      call check(error,.not. err%has_error(), "the reference SCF")
      if (allocated(error)) return

      call mol%overlap(overlap)
      allocate (spin_density(size(scf%density, 1), size(scf%density, 2)), source=0.0_dp)
      call mulliken_spin_populations(mol, spin_density, overlap, populations, err)
      call check(error,.not. err%has_error(), "computing spin populations")
      if (allocated(error)) return
      call check(error, maxval(abs(populations)) < 1.0e-12_dp, &
                 "a closed shell must carry no spin population anywhere")
   end subroutine test_spin_closed_shell

   subroutine test_grid_excludes(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: points(:, :)
      real(dp) :: distance, closest
      integer :: ipoint, iatom

      call water(mol, err)
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return

      call chelpg_grid(mol, points, err)
      call check(error,.not. err%has_error(), "building the CHELPG grid")
      if (allocated(error)) return
      call check(error, size(points, 2) > 100, "the grid is suspiciously small")
      if (allocated(error)) return

      ! Every kept point is outside every atom's excluded radius. Half a Bohr
      ! is well inside any van der Waals radius, so a point closer than that to
      ! a nucleus means the exclusion did not run at all.
      closest = huge(1.0_dp)
      do ipoint = 1, size(points, 2)
         do iatom = 1, mol%natm
            distance = norm2(points(:, ipoint) - mol%coords(:, iatom))
            closest = min(closest, distance)
         end do
      end do
      call check(error, closest > 0.5_dp, &
                 "a fitting point sits on top of a nucleus")
   end subroutine test_grid_excludes

end module test_mqc_czt_charges

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_charges, only: collect_mqc_czt_charges_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_charges", collect_mqc_czt_charges_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
