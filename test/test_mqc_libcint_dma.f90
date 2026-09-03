!! Distributed multipoles: the moments they must reproduce, and thread safety
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
   use mqc_libcint_screening, only: screening_grid
   use mqc_error, only: error_t
   use mqc_calculation_defaults, only: DEFAULT_VDW_SCALE
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   implicit none
   private

   public :: collect_mqc_libcint_dma_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   integer, parameter :: Z(3) = [8, 1, 1]

   !! A water dimer for the threading test, where the water above will not do.
   !!
   !! Two molecules give bond midpoints that several primitive pairs can claim,
   !! and a partition with something to partition -- which is where a missing
   !! reduction shows. One water is small enough that a broken thread split can
   !! still look right.
   integer, parameter :: Z_DIMER(6) = [8, 1, 1, 8, 1, 1]
   character(len=2), parameter :: SYM_DIMER(6) = ["O ", "H ", "H ", "O ", "H ", "H "]
   real(dp), parameter :: GEO_DIMER(3, 6) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 0.0_dp, 5.6_dp, &
                           0.0_dp, 1.4308_dp, 6.7078_dp, &
                           0.0_dp, -1.4308_dp, 6.7078_dp], [3, 6])

contains

   subroutine collect_mqc_libcint_dma_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("expansion_points_are_atoms_and_midpoints", test_points), &
                  new_unittest("bond_labels_past_nine_atoms_read_back", test_wide_labels), &
                  new_unittest("site_charges_sum_to_the_molecular_charge", test_charge_sum), &
                  new_unittest("sites_reproduce_the_molecular_dipole", test_dipole_sum), &
                  new_unittest("thread_count_does_not_change_the_multipoles", &
                               test_thread_invariance) &
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

   subroutine test_wide_labels(error)
      !! A bond whose atoms need two digits still names the right pair
      !!
      !! Midpoint labels were written `"BO"//i//j` at minimal width and read one
      !! character per index, so "BO123" was atoms 12 and 3 or atoms 1 and 23
      !! with nothing to tell them apart -- and the reader took 1 and 2. Where
      !! the higher atom was a multiple of ten the second index came back as 0
      !! and indexed `atomic_numbers(0)`.
      !!
      !! Three atoms cannot see any of that, which is why the fixture here is a
      !! ten-carbon chain: it puts a bond on atoms 10 and 9, the case the old
      !! encoding could not express, and runs it through the reader that broke.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NC = 10
      real(dp), parameter :: SPACING = 2.8_dp
         !! Bohr, so 1.48 Angstrom -- inside a C-C bond, while the next-nearest
         !! pair at 2.96 is well outside one. The chain bonds only adjacently.
      type(error_t) :: err
      type(dma_result_t) :: dma
      integer :: znum(NC)
      real(dp) :: c(3, NC)
      real(dp), allocatable :: points(:, :), nuclear(:), grid(:, :)
      character(len=8), allocatable :: labels(:)
      integer :: i, hi, lo, found

      znum = 6
      c = 0.0_dp
      do i = 1, NC
         c(3, i) = real(i - 1, dp)*SPACING
      end do

      call expansion_points(znum, c, points, labels, nuclear, err)
      call check(error,.not. err%has_error(), "building the expansion points")
      if (allocated(error)) return

      call check(error, size(labels) == NC + (NC - 1), &
                 "a ten-atom chain should give ten atoms and nine bond midpoints")
      if (allocated(error)) return

      found = 0
      do i = NC + 1, size(labels)
         if (labels(i) == "BO010009") found = i
      end do
      call check(error, found > 0, &
                 "no BO010009 label: the bond between atoms 10 and 9 is missing or "// &
                 "written at minimal width, where it collides with other pairs")
      if (allocated(error)) return

      ! The same substrings `screening_grid` reads.
      read (labels(found) (3:5), *) hi
      read (labels(found) (6:8), *) lo
      call check(error, hi == 10, "label characters 3:5 did not read back as atom 10")
      if (allocated(error)) return
      call check(error, lo == 9, "label characters 6:8 did not read back as atom 9")
      if (allocated(error)) return

      ! And through the reader itself, which is where the wrong pair -- or index
      ! zero -- used to be picked up.
      allocate (dma%points(3, size(points, 2)))
      dma%points = points
      allocate (dma%labels(size(labels)))
      dma%labels = labels
      call screening_grid(dma, znum, grid, err, vdw_scale=DEFAULT_VDW_SCALE)
      call check(error,.not. err%has_error(), "the screening grid could not be built")
      if (allocated(error)) return
      call check(error, size(grid, 2) > 0, "the screening grid came out empty")
   end subroutine test_wide_labels

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

   subroutine dimer_multipoles(dma, nelec, err, ok)
      type(dma_result_t), intent(out) :: dma
      integer, intent(out) :: nelec
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf

      ok = .false.
      nelec = sum(Z_DIMER)
      call build_libcint_molecule(Z_DIMER, SYM_DIMER, GEO_DIMER, "sto-3g", mol, err)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if
      call run_libcint_rhf(mol, nelec, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         call mol%destroy()
         return
      end if
      call distributed_multipoles(mol, scf%density, Z_DIMER, dma, err)
      ok = .not. err%has_error()
      call mol%destroy()
   end subroutine dimer_multipoles

   subroutine test_thread_invariance(error)
      !! One thread and many must give the same multipoles
      !!
      !! The reduction in the monopole loop is the part that can go wrong: every
      !! primitive pair may contribute to any expansion point, so a missing
      !! private copy shows up as a thread-count-dependent answer rather than as
      !! a crash. The dipole and above accumulate into one column per point and
      !! cannot race at all, which is worth stating because it is why only one
      !! of the three loops needed a reduction.
      type(error_type), allocatable, intent(out) :: error

      type(dma_result_t) :: one, many
      type(error_t) :: err
      integer :: nelec, saved
      logical :: ok
      real(dp) :: worst

      saved = omp_get_max_threads()
      if (saved < 2) then
         ! Nothing to compare against; not a failure, just no information.
         return
      end if

      call omp_set_num_threads(1)
      call dimer_multipoles(one, nelec, err, ok)
      call check(error, ok, "the single-threaded reference failed")
      if (allocated(error)) then
         call omp_set_num_threads(saved)
         return
      end if

      call omp_set_num_threads(saved)
      call dimer_multipoles(many, nelec, err, ok)
      call check(error, ok, "the threaded run failed: "//err%get_message())
      if (allocated(error)) return

      worst = max(maxval(abs(one%electronic - many%electronic)), &
                  maxval(abs(one%dipole - many%dipole)))
      worst = max(worst, maxval(abs(one%quadrupole - many%quadrupole)))
      worst = max(worst, maxval(abs(one%octopole - many%octopole)))

      call check(error, maxval(abs(one%electronic)) > 0.1_dp, &
                 "the multipoles are empty, so this compares nothing")
      if (allocated(error)) return
      ! Summation order differs between one thread and many, so this is a
      ! rounding-level bound rather than bit equality.
      call check(error, worst < 1.0e-12_dp, &
                 "the multipoles depend on the thread count")
   end subroutine test_thread_invariance

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
