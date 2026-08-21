!! Foster-Boys localization, against what a rotation cannot change
module test_mqc_libcint_localize
   !! Localized orbitals are not unique numbers to compare against -- any
   !! program's Boys orbitals differ by phase, by ordering, and by however far
   !! its sweeps ran. What is not free is everything the transformation is
   !! obliged to preserve, and those are exactly the properties that a broken
   !! Jacobi rotation quietly violates:
   !!
   !!   * **The occupied space is unchanged.** Localization is a rotation among
   !!     occupied orbitals, so the density matrix built from them is identical
   !!     to the canonical one. A rotation applied to the wrong pair of columns,
   !!     or one that leaks into a virtual, changes the density and therefore
   !!     the energy -- while still producing orbitals that look localized.
   !!   * **Orthonormality survives.** `C^T S C = I` over the occupied block.
   !!   * **The functional does not go down.** Boys maximizes the sum of squared
   !!     centroids, and each sweep is a sequence of exact two-by-two maxima, so
   !!     the localized value must be at least the canonical one. On any molecule
   !!     with a lone pair it is strictly greater, which is what says the sweep
   !!     did anything at all.
   !!   * **The centroids sit inside the molecule.** A bond orbital's centroid
   !!     lies between its atoms; a centroid several Bohr outside the nuclear
   !!     frame means the dipole matrices reached the sweep with the wrong sign
   !!     or the wrong origin.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: run_libcint_rhf, rhf_result_t
   use mqc_libcint_localize, only: boys_localize
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_libcint_localize_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_libcint_localize_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("boys_preserves_the_occupied_space", test_density_invariant), &
                  new_unittest("boys_orbitals_stay_orthonormal", test_orthonormal), &
                  new_unittest("boys_increases_its_own_functional", test_functional), &
                  new_unittest("boys_centroids_lie_in_the_molecule", test_centroids) &
                  ]
   end subroutine collect_mqc_libcint_localize_tests

   subroutine converged_water(mol, scf, err)
      type(libcint_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      real(dp) :: c(3, 3)
      integer :: z(3)
      character(len=2) :: symbols(3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_libcint_molecule(z, symbols, c, "6-31g", mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 60, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
   end subroutine converged_water

   pure function occupied_density(c, n_occ) result(d)
      !! The closed-shell density from the first `n_occ` columns.
      real(dp), intent(in) :: c(:, :)
      integer, intent(in) :: n_occ
      real(dp) :: d(size(c, 1), size(c, 1))

      d = 2.0_dp*matmul(c(:, 1:n_occ), transpose(c(:, 1:n_occ)))
   end function occupied_density

   pure function boys_functional(centroids) result(value)
      !! The quantity Boys maximizes: the sum of squared orbital centroids.
      real(dp), intent(in) :: centroids(:, :)
      real(dp) :: value

      value = sum(centroids**2)
   end function boys_functional

   subroutine test_density_invariant(error)
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: localized(:, :), centroids(:, :)
      real(dp), allocatable :: before(:, :), after(:, :)

      call converged_water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF")
      if (allocated(error)) return

      before = occupied_density(scf%orbitals, scf%n_occupied)
      call boys_localize(mol, scf%orbitals, scf%n_occupied, localized, centroids, err)
      call check(error,.not. err%has_error(), "localizing")
      if (allocated(error)) return

      after = occupied_density(localized, scf%n_occupied)
      call check(error, maxval(abs(after - before)) < 1.0e-10_dp, &
                 "localization changed the density, so it was not a rotation "// &
                 "within the occupied space")
   end subroutine test_density_invariant

   subroutine test_orthonormal(error)
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: localized(:, :), centroids(:, :), overlap(:, :)
      real(dp), allocatable :: metric(:, :)
      real(dp) :: worst
      integer :: i, j, n

      call converged_water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF")
      if (allocated(error)) return

      call mol%overlap(overlap)
      call boys_localize(mol, scf%orbitals, scf%n_occupied, localized, centroids, err)
      call check(error,.not. err%has_error(), "localizing")
      if (allocated(error)) return

      n = scf%n_occupied
      metric = matmul(transpose(localized(:, 1:n)), matmul(overlap, localized(:, 1:n)))
      worst = 0.0_dp
      do j = 1, n
         do i = 1, n
            if (i == j) then
               worst = max(worst, abs(metric(i, j) - 1.0_dp))
            else
               worst = max(worst, abs(metric(i, j)))
            end if
         end do
      end do
      call check(error, worst < 1.0e-10_dp, "the localized orbitals are not orthonormal")
   end subroutine test_orthonormal

   subroutine test_functional(error)
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: localized(:, :), centroids(:, :)
      real(dp), allocatable :: canonical(:, :), canonical_centroids(:, :)
      real(dp) :: before, after

      call converged_water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF")
      if (allocated(error)) return

      ! Zero sweeps returns the canonical orbitals and their centroids, which is
      ! the starting value of the functional -- asked for through the same
      ! routine so the two values are computed identically.
      call boys_localize(mol, scf%orbitals, scf%n_occupied, canonical, &
                         canonical_centroids, err, max_sweeps=0)
      call check(error,.not. err%has_error(), "the unlocalized reference")
      if (allocated(error)) return
      call boys_localize(mol, scf%orbitals, scf%n_occupied, localized, centroids, err)
      call check(error,.not. err%has_error(), "localizing")
      if (allocated(error)) return

      before = boys_functional(canonical_centroids)
      after = boys_functional(centroids)
      call check(error, after >= before - 1.0e-10_dp, &
                 "the Boys functional went down, which no sequence of exact "// &
                 "two-by-two maxima can do")
      if (allocated(error)) return
      ! Water's canonical orbitals are delocalized enough that localizing them
      ! must move the functional appreciably; equality would mean no rotation
      ! was applied at all.
      call check(error, after > before + 1.0e-3_dp, &
                 "localization left the functional where it started")
   end subroutine test_functional

   subroutine test_centroids(error)
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: localized(:, :), centroids(:, :)
      real(dp) :: nearest, distance
      integer :: iorb, iatom

      call converged_water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF")
      if (allocated(error)) return

      call boys_localize(mol, scf%orbitals, scf%n_occupied, localized, centroids, err)
      call check(error,.not. err%has_error(), "localizing")
      if (allocated(error)) return
      call check(error, size(centroids, 2) == scf%n_occupied, "one centroid per orbital")
      if (allocated(error)) return

      ! Every centroid within two Bohr of some nucleus: a core orbital sits on
      ! its atom, a bond between two, a lone pair just outside one. Anything
      ! further out is not a localized orbital of this molecule.
      do iorb = 1, size(centroids, 2)
         nearest = huge(1.0_dp)
         do iatom = 1, mol%natm
            distance = norm2(centroids(:, iorb) - mol%coords(:, iatom))
            nearest = min(nearest, distance)
         end do
         call check(error, nearest < 2.0_dp, &
                    "a Boys centroid lies outside the molecule")
         if (allocated(error)) return
      end do
   end subroutine test_centroids

end module test_mqc_libcint_localize

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_libcint_localize, only: collect_mqc_libcint_localize_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_libcint_localize", collect_mqc_libcint_localize_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
