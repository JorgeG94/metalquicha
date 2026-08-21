!! Distributed multipoles, against the molecular moments they must reproduce
module test_mqc_libcint_dma
   !! A distributed multipole analysis spreads one molecule's charge
   !! distribution over many expansion points -- every atom and every bond
   !! midpoint here -- and the individual site moments have no meaning on their
   !! own to compare against. What is fixed is what they must add back up to:
   !!
   !!   * **The site charges sum to the molecular charge.** The electronic part
   !!     sums to minus the electron count and the nuclear part to the total
   !!     atomic number, so for a neutral molecule the two cancel exactly. A
   !!     site whose share was computed with the wrong weight breaks this while
   !!     leaving every individual number looking reasonable.
   !!   * **Charges and dipoles together reproduce the molecular dipole.** The
   !!     total dipole about the origin is `sum_s (q_s r_s + mu_s)`, and it must
   !!     equal the one computed directly as `-Tr(D r) + sum_A Z_A R_A`. This is
   !!     the check that the site dipoles carry the right sign and the right
   !!     origin, which is exactly what a `.efp` file depends on and what a
   !!     GAMESS comparison one rung further up would only tell you about
   !!     indirectly.
   !!   * **The expansion points are the molecule's own geometry.** Every atom
   !!     appears, every bond midpoint lies between two atoms, and the nuclear
   !!     charge sits on atoms and nowhere else.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: run_libcint_rhf, rhf_result_t
   use mqc_libcint_multipole, only: multipole_matrices
   use mqc_libcint_dma, only: dma_result_t, distributed_multipoles, expansion_points
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_libcint_dma_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   integer, parameter :: Z(3) = [8, 1, 1]

contains

   subroutine collect_mqc_libcint_dma_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("expansion_points_are_atoms_and_midpoints", test_points), &
                  new_unittest("site_charges_sum_to_the_molecular_charge", test_charge_sum), &
                  new_unittest("sites_reproduce_the_molecular_dipole", test_dipole_sum) &
                  ]
   end subroutine collect_mqc_libcint_dma_tests

   subroutine geometry(c)
      real(dp), intent(out) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
   end subroutine geometry

   subroutine converged_water(mol, scf, err)
      type(libcint_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      real(dp) :: c(3, 3)
      character(len=2) :: symbols(3)

      symbols = ["O ", "H ", "H "]
      call geometry(c)
      call build_libcint_molecule(Z, symbols, c, "6-31g", mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 60, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
   end subroutine converged_water

   subroutine test_points(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp) :: c(3, 3)
      real(dp), allocatable :: points(:, :), nuclear(:)
      character(len=8), allocatable :: labels(:)
      real(dp) :: nearest, distance
      integer :: ipoint, iatom

      call geometry(c)
      call expansion_points(Z, c, points, labels, nuclear, err)
      call check(error,.not. err%has_error(), "building the expansion points")
      if (allocated(error)) return

      ! Three atoms, and the bonds between them: more points than atoms, and a
      ! label for every one.
      call check(error, size(points, 2) > 3, "no bond midpoints were added")
      if (allocated(error)) return
      call check(error, size(labels) == size(points, 2), "a point without a label")
      if (allocated(error)) return
      call check(error, size(nuclear) == size(points, 2), "a point without a nuclear charge")
      if (allocated(error)) return

      ! The first three points are the atoms themselves, carrying their charge.
      do iatom = 1, 3
         call check(error, norm2(points(:, iatom) - c(:, iatom)) < 1.0e-12_dp, &
                    "an atomic expansion point is not on its atom")
         if (allocated(error)) return
         call check(error, abs(nuclear(iatom) - real(Z(iatom), dp)) < 1.0e-12_dp, &
                    "an atomic site does not carry its nuclear charge")
         if (allocated(error)) return
      end do

      ! Every further point is a midpoint, so it lies inside the molecule and
      ! carries no nucleus.
      do ipoint = 4, size(points, 2)
         call check(error, abs(nuclear(ipoint)) < 1.0e-12_dp, &
                    "a bond midpoint carries nuclear charge")
         if (allocated(error)) return
         nearest = huge(1.0_dp)
         do iatom = 1, 3
            distance = norm2(points(:, ipoint) - c(:, iatom))
            nearest = min(nearest, distance)
         end do
         call check(error, nearest < 2.0_dp, "a bond midpoint is nowhere near the molecule")
         if (allocated(error)) return
      end do
   end subroutine test_points

   subroutine test_charge_sum(error)
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(dma_result_t) :: dma

      call converged_water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF")
      if (allocated(error)) return

      call distributed_multipoles(mol, scf%density, Z, dma, err)
      call check(error,.not. err%has_error(), "the distributed multipole analysis")
      if (allocated(error)) return

      ! Ten electrons, spread over however many sites.
      call check(error, abs(sum(dma%electronic) + 10.0_dp) < 1.0e-9_dp, &
                 "the electronic site charges do not account for every electron")
      if (allocated(error)) return
      call check(error, abs(sum(dma%nuclear) - 10.0_dp) < 1.0e-12_dp, &
                 "the nuclear site charges do not sum to the nuclear charge")
      if (allocated(error)) return
      call check(error, abs(sum(dma%electronic) + sum(dma%nuclear)) < 1.0e-9_dp, &
                 "the molecule did not come out neutral")
   end subroutine test_charge_sum

   subroutine test_dipole_sum(error)
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(dma_result_t) :: dma
      real(dp), allocatable :: dip(:, :, :)
      real(dp) :: c(3, 3), from_sites(3), direct(3)
      integer :: comp, ipoint, iatom, i, j

      call converged_water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF")
      if (allocated(error)) return
      call geometry(c)

      call distributed_multipoles(mol, scf%density, Z, dma, err)
      call check(error,.not. err%has_error(), "the distributed multipole analysis")
      if (allocated(error)) return

      ! What the distribution says the molecular dipole is: each site's charge
      ! at its position, plus each site's own dipole.
      from_sites = 0.0_dp
      do ipoint = 1, size(dma%points, 2)
         from_sites = from_sites &
                      + (dma%electronic(ipoint) + dma%nuclear(ipoint))*dma%points(:, ipoint) &
                      + dma%dipole(:, ipoint)
      end do

      ! What it is, from the density and the nuclei directly.
      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, err)
      call check(error,.not. err%has_error(), "the dipole matrices")
      if (allocated(error)) return
      direct = 0.0_dp
      do comp = 1, 3
         do j = 1, mol%nao
            do i = 1, mol%nao
               direct(comp) = direct(comp) - scf%density(i, j)*dip(i, j, comp)
            end do
         end do
         do iatom = 1, 3
            direct(comp) = direct(comp) + real(Z(iatom), dp)*c(comp, iatom)
         end do
      end do

      call check(error, maxval(abs(from_sites - direct)) < 1.0e-8_dp, &
                 "the distributed multipoles do not reproduce the molecular dipole")
   end subroutine test_dipole_sum

end module test_mqc_libcint_dma

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_libcint_dma, only: collect_mqc_libcint_dma_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_libcint_dma", collect_mqc_libcint_dma_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
