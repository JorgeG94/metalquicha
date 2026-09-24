!! Foster-Boys and Edmiston-Ruedenberg localization, against what a rotation cannot change
module test_mqc_czt_localize
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
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: run_czt_rhf, rhf_result_t
   use mqc_czt_localize, only: boys_localize, er_localize, occupied_eri
   use mqc_czt_mp2, only: transform_block
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_czt_localize_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_czt_localize_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("boys_preserves_the_occupied_space", test_density_invariant), &
                  new_unittest("boys_orbitals_stay_orthonormal", test_orthonormal), &
                  new_unittest("boys_increases_its_own_functional", test_functional), &
                  new_unittest("boys_centroids_lie_in_the_molecule", test_centroids), &
                  new_unittest("occupied_eri_matches_the_ao_tensor", test_occupied_eri), &
                  new_unittest("er_preserves_the_occupied_space", test_er_projector), &
                  new_unittest("er_rises_with_every_sweep", test_er_monotone), &
                  new_unittest("er_does_not_depend_on_the_starting_rotation", test_er_start), &
                  new_unittest("er_matches_gamess_and_pyscf", test_er_references) &
                  ]
   end subroutine collect_mqc_czt_localize_tests

   subroutine converged_water(mol, scf, err)
      type(czt_molecule_t), intent(out) :: mol
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
      call build_czt_molecule(z, symbols, c, "6-31g", mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, 10, 60, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
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

      type(czt_molecule_t) :: mol
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

      type(czt_molecule_t) :: mol
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

      type(czt_molecule_t) :: mol
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

      type(czt_molecule_t) :: mol
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

   subroutine small_molecule(name, basis, mol, scf, err)
      !! One of the molecules the ER references were computed on
      !!
      !! Geometries are Angstrom times this file's `ANG`, and the GAMESS and
      !! PySCF runs were given the same Bohr coordinates, so no conversion
      !! constant separates them.
      character(len=*), intent(in) :: name, basis
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      real(dp), allocatable :: c(:, :)
      integer, allocatable :: z(:)
      character(len=2), allocatable :: symbols(:)

      select case (name)
      case ("water")
         z = [8, 1, 1]
         symbols = ["O ", "H ", "H "]
         c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                      0.0_dp, 0.0_dp, 0.9584_dp, &
                      0.9268_dp, 0.0_dp, -0.2400_dp], [3, 3])
      case ("ethane")
         ! Staggered: C-C 1.530, C-H 1.094, H-C-C 111.2 degrees.
         z = [6, 6, 1, 1, 1, 1, 1, 1]
         symbols = ["C ", "C ", "H ", "H ", "H ", "H ", "H ", "H "]
         c = reshape([0.0_dp, 0.0_dp, 0.765_dp, &
                      0.0_dp, 0.0_dp, -0.765_dp, &
                      1.019962238529770_dp, 0.0_dp, 1.160617279669809_dp, &
                      -0.509981119264885_dp, 0.883313209467624_dp, 1.160617279669809_dp, &
                      -0.509981119264886_dp, -0.883313209467624_dp, 1.160617279669809_dp, &
                      0.509981119264885_dp, 0.883313209467624_dp, -1.160617279669809_dp, &
                      -1.019962238529770_dp, 0.0_dp, -1.160617279669809_dp, &
                      0.509981119264885_dp, -0.883313209467624_dp, -1.160617279669809_dp], &
                     [3, 8])
      case default
         z = [6, 8, 1, 1]
         symbols = ["C ", "O ", "H ", "H "]
         c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                      0.0_dp, 0.0_dp, 1.205_dp, &
                      0.0_dp, 0.9429_dp, -0.5876_dp, &
                      0.0_dp, -0.9429_dp, -0.5876_dp], [3, 4])
      end select
      c = c*ANG
      call build_czt_molecule(z, symbols, c, basis, mol, err)
      if (err%has_error()) return
      ! The commutator bound stated, not derived: sum (ii|ii) is first order in
      ! the orbitals, and sqrt(energy_tol) left it 2.4e-8 off on water/6-31G.
      call run_czt_rhf(mol, sum(z), 200, 1.0e-11_dp, 1.0e-10_dp, .false., scf, err, &
                       grad_tol=1.0e-9_dp)
   end subroutine small_molecule

   subroutine test_occupied_eri(error)
      !! The direct occupied transform against the packed AO tensor
      !!
      !! Formaldehyde twice: STO-3G, whose L shells the integral loop sees
      !! fused, and 6-31G*, which puts a d shell in every heavy-atom quartet.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: v(:, :), ao(:, :), full(:, :, :, :), c(:, :)
      real(dp) :: worst
      integer :: i, j, k, l, n, run
      character(len=8), parameter :: BASES(2) = ["sto-3g  ", "6-31g*  "]

      do run = 1, 2
         call small_molecule("formaldehyde", trim(BASES(run)), mol, scf, err)
         call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF")
         if (allocated(error)) return
         n = scf%n_occupied
         c = scf%orbitals(:, 1:n)
         call occupied_eri(mol, c, v, err)
         call check(error,.not. err%has_error(), "the direct transform")
         if (allocated(error)) return
         call mol%eris_packed(ao)
         call transform_block(ao, c, c, c, c, full)
         worst = 0.0_dp
         do l = 1, n
            do k = 1, n
               do j = 1, n
                  do i = 1, n
                     worst = max(worst, abs(full(i, j, k, l) - &
                                            v(pair(i, j), pair(k, l))))
                  end do
               end do
            end do
         end do
         write (*, "(a,a,a,es10.2)") "    ", BASES(run), " max |direct - in-core| =", worst
         call check(error, worst < 1.0e-12_dp, &
                    "the direct occupied transform disagrees with the AO tensor")
         if (allocated(error)) return
         call mol%destroy()
      end do
   end subroutine test_occupied_eri

   pure integer function pair(i, j)
      integer, intent(in) :: i, j
      pair = max(i, j)*(max(i, j) - 1)/2 + min(i, j)
   end function pair

   subroutine test_er_projector(error)
      !! ER is a rotation among the occupied orbitals: `C C^T` does not move
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: localized(:, :), centroids(:, :), overlap(:, :), metric(:, :)
      real(dp) :: worst
      integer :: n, i

      call converged_water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF")
      if (allocated(error)) return
      n = scf%n_occupied
      call er_localize(mol, scf%orbitals, n, localized, centroids, err)
      call check(error,.not. err%has_error(), "localizing")
      if (allocated(error)) return

      worst = maxval(abs(occupied_density(localized, n) - occupied_density(scf%orbitals, n)))
      write (*, "(a,es10.2)") "    max |P_er - P_canonical| =", 0.5_dp*worst
      call check(error, 0.5_dp*worst <= 1.0e-12_dp, &
                 "ER changed the occupied projector, so it was not a rotation within it")
      if (allocated(error)) return

      call mol%overlap(overlap)
      metric = matmul(transpose(localized), matmul(overlap, localized))
      do i = 1, n
         metric(i, i) = metric(i, i) - 1.0_dp
      end do
      call check(error, maxval(abs(metric)) < 1.0e-12_dp, &
                 "the ER orbitals are not orthonormal")
   end subroutine test_er_projector

   subroutine test_er_monotone(error)
      !! `sum_i (ii|ii)` never falls, sweep by sweep
      !!
      !! Every rotation is the exact maximum of its pair, so the functional
      !! can only rise; a wrong sign in A or B still converges, but to a point
      !! that is not the maximum, and on the way it goes down.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: localized(:, :), centroids(:, :)
      real(dp) :: previous, now, first
      integer :: sweeps, taken, final_sweeps

      call converged_water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF")
      if (allocated(error)) return

      call er_localize(mol, scf%orbitals, scf%n_occupied, localized, centroids, err, &
                       sweeps_taken=final_sweeps, functional=now)
      call check(error,.not. err%has_error(), "localizing")
      if (allocated(error)) return

      previous = -huge(1.0_dp)
      do sweeps = 0, final_sweeps
         call er_localize(mol, scf%orbitals, scf%n_occupied, localized, centroids, err, &
                          max_sweeps=sweeps, sweeps_taken=taken, functional=now)
         call check(error,.not. err%has_error(), "localizing")
         if (allocated(error)) return
         if (sweeps == 0) first = now
         call check(error, now >= previous - 1.0e-12_dp, &
                    "the ER functional went down during a sweep")
         if (allocated(error)) return
         previous = now
      end do
      write (*, "(a,f16.10,a,f16.10,a,i0,a)") "    sum (ii|ii):", first, " ->", now, &
         " in ", final_sweeps, " sweeps"
      call check(error, now > first + 1.0e-2_dp, &
                 "localization left the ER functional where it started")
   end subroutine test_er_monotone

   subroutine test_er_start(error)
      !! The same localized set from the canonical orbitals and from a rotation of them
      !!
      !! Compared through the overlap: each localized orbital from the rotated
      !! start must be, up to sign, one from the canonical start.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: a(:, :), b(:, :), ca(:, :), cb(:, :), overlap(:, :)
      real(dp), allocatable :: mixed(:, :), cross(:, :), u(:, :)
      real(dp) :: da, db, worst, angle, cs, sn, x, y
      integer :: n, i, j, k, m

      call converged_water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF")
      if (allocated(error)) return
      n = scf%n_occupied

      ! An arbitrary orthogonal mix: a Givens rotation on every pair.
      allocate (u(n, n))
      u = 0.0_dp
      do i = 1, n
         u(i, i) = 1.0_dp
      end do
      k = 0
      do j = 2, n
         do i = 1, j - 1
            k = k + 1
            angle = 0.37_dp*k
            cs = cos(angle)
            sn = sin(angle)
            do m = 1, n
               x = u(m, i)
               y = u(m, j)
               u(m, i) = cs*x - sn*y
               u(m, j) = sn*x + cs*y
            end do
         end do
      end do
      mixed = matmul(scf%orbitals(:, 1:n), u)

      call er_localize(mol, scf%orbitals, n, a, ca, err, functional=da)
      call check(error,.not. err%has_error(), "localizing the canonical orbitals")
      if (allocated(error)) return
      call er_localize(mol, mixed, n, b, cb, err, functional=db)
      call check(error,.not. err%has_error(), "localizing the rotated orbitals")
      if (allocated(error)) return

      call mol%overlap(overlap)
      cross = abs(matmul(transpose(a), matmul(overlap, b)))
      worst = 0.0_dp
      do i = 1, n
         worst = max(worst, abs(1.0_dp - maxval(cross(i, :))))
      end do
      write (*, "(a,es10.2,a,es10.2)") "    functional difference", db - da, &
         ",  worst 1 - max|<a_i|b_j>|", worst
      call check(error, abs(db - da) < 1.0e-10_dp, &
                 "the rotated start reached a different ER functional")
      if (allocated(error)) return
      call check(error, worst < 1.0e-8_dp, &
                 "the rotated start reached a different set of localized orbitals")
   end subroutine test_er_start

   subroutine test_er_references(error)
      !! Water, ethane and formaldehyde in STO-3G and 6-31G, against two codes
      !!
      !! GAMESS 2026 (`../mgga/gamess`), `LOCAL=RUEDNBRG`, `$LOCAL FCORE=.F.
      !! CVGLOC=1D-11`, RHF converged to 1e-10, `ICUT=12 ITOL=30`,
      !! `$TRANS CUTTRF=1D-14`, each basis spelt out from this repository's
      !! JSON: its printed `DIAGONAL SUM D`. PySCF 2.14's second-order
      !! `lo.ER`, started from GAMESS's localized orbitals and converged to
      !! 1e-12: its sum, and its centroids, which carry more digits than the
      !! nine GAMESS punches. The two agree to 5e-10.
      !!
      !! PySCF started from the canonical orbitals instead stops at lower
      !! stationary points on five of the six -- 8.1146 on water/STO-3G against
      !! 8.4550 -- which is why it is not the reference for which maximum.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: localized(:, :), centroids(:, :), ref(:, :)
      real(dp) :: d, nearest, worst
      integer :: run, i, j, n
      character(len=12), parameter :: NAMES(6) = [character(len=12) :: "water", "water", &
                                                                                     "ethane", "ethane", "formaldehyde", &
                                                                                                           "formaldehyde"]
      character(len=6), parameter :: BASES(6) = [character(len=6) :: "sto-3g", "6-31g", &
                                                                                     "sto-3g", "6-31g", "sto-3g", "6-31g"]
      real(dp), parameter :: D_GAMESS(6) = [8.4550244411_dp, 8.2686027165_dp, &
                                            11.9933311078_dp, 11.8830104139_dp, &
                                            13.2585043555_dp, 13.1035586880_dp]
      real(dp), parameter :: D_PYSCF(6) = [8.4550244412_dp, 8.2686027170_dp, &
                                           11.9933311079_dp, 11.8830104137_dp, &
                                           13.2585043558_dp, 13.1035586884_dp]

      do run = 1, 6
         call small_molecule(trim(NAMES(run)), trim(BASES(run)), mol, scf, err)
         call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF")
         if (allocated(error)) return
         n = scf%n_occupied
         call er_localize(mol, scf%orbitals, n, localized, centroids, err, functional=d)
         call check(error,.not. err%has_error(), "localizing")
         if (allocated(error)) return
         call reference_centroids(run, ref)

         worst = 0.0_dp
         do i = 1, n
            nearest = huge(1.0_dp)
            do j = 1, n
               nearest = min(nearest, norm2(centroids(:, i) - ref(:, j)))
            end do
            worst = max(worst, nearest)
         end do
         write (*, "(4x,a,1x,a,a,f15.10,a,es9.1,a,es9.1,a,es9.1)") NAMES(run), BASES(run), &
            " D =", d, "  - GAMESS", d - D_GAMESS(run), "  - PySCF", d - D_PYSCF(run), &
            "  centroids", worst
         call check(error, abs(d - D_GAMESS(run)) < 1.0e-8_dp, &
                    "sum (ii|ii) does not match GAMESS's ER maximum")
         if (allocated(error)) return
         call check(error, abs(d - D_PYSCF(run)) < 1.0e-8_dp, &
                    "sum (ii|ii) does not match PySCF's")
         if (allocated(error)) return
         call check(error, worst < 1.0e-7_dp, &
                    "the ER orbital centroids are not PySCF's")
         if (allocated(error)) return
         call mol%destroy()
      end do
   end subroutine test_er_references

   subroutine reference_centroids(run, ref)
      !! PySCF's ER centroids, Bohr, (3, n_occ), in `test_er_references`' order
      integer, intent(in) :: run
      real(dp), allocatable, intent(out) :: ref(:, :)

      select case (run)
      case (1)
         ref = reshape([ &
                       -0.00049268_dp, -0.00000000_dp, -0.00038274_dp, 1.05357179_dp, -0.00000000_dp, &
                       -0.25908345_dp, 0.01324508_dp, -0.00000000_dp, 1.08570400_dp, -0.22950825_dp, &
                       0.50608004_dp, -0.17752882_dp, -0.22950825_dp, -0.50608004_dp, -0.17752882_dp], &
                       [3, 5])
      case (2)
         ref = reshape([ &
                       -0.00042643_dp, -0.00000000_dp, -0.00033139_dp, 0.95764440_dp, -0.00000000_dp, &
                       -0.23413580_dp, 0.01340348_dp, 0.00000000_dp, 0.98635501_dp, -0.25202528_dp, &
                       -0.51289431_dp, -0.19494494_dp, -0.25202528_dp, 0.51289431_dp, -0.19494494_dp], &
                       [3, 5])
      case (3)
         ref = reshape([ &
                       -0.00000000_dp, 0.00000000_dp, 1.44548005_dp, -0.00000000_dp, -0.00000000_dp, &
                       -1.44548005_dp, 0.66234246_dp, 1.14721079_dp, -1.94599522_dp, -0.66234246_dp, &
                       -1.14721080_dp, 1.94599522_dp, 0.66234246_dp, -1.14721079_dp, -1.94599522_dp, &
                       -1.32468493_dp, -0.00000000_dp, -1.94599522_dp, 0.00000000_dp, 0.00000000_dp, &
                       0.00000000_dp, -0.66234246_dp, 1.14721079_dp, 1.94599522_dp, 1.32468492_dp, &
                       0.00000000_dp, 1.94599522_dp], [3, 9])
      case (4)
         ref = reshape([ &
                       0.00000000_dp, -0.00000000_dp, -1.44536500_dp, 0.00000000_dp, 0.00000000_dp, &
                       1.44536500_dp, -0.65194591_dp, -1.12920343_dp, 1.94225470_dp, -1.30389181_dp, &
                       0.00000000_dp, -1.94225470_dp, -0.65194589_dp, 1.12920343_dp, 1.94225470_dp, &
                       1.30389181_dp, -0.00000001_dp, 1.94225469_dp, -0.00000000_dp, 0.00000001_dp, &
                       0.00000001_dp, 0.65194590_dp, -1.12920342_dp, -1.94225470_dp, 0.65194590_dp, &
                       1.12920342_dp, -1.94225470_dp], [3, 9])
      case (5)
         ref = reshape([ &
                       -0.00000000_dp, -0.00000000_dp, 2.27750815_dp, -0.00000000_dp, 0.00000000_dp, &
                       0.00005042_dp, 0.49526324_dp, -0.00000001_dp, 1.28571494_dp, -0.00000000_dp, &
                       -1.21782149_dp, -0.76738797_dp, 0.00000000_dp, 1.21782149_dp, -0.76738797_dp, &
                       0.00000000_dp, 0.51542139_dp, 2.49143853_dp, -0.49526324_dp, 0.00000001_dp, &
                       1.28571494_dp, -0.00000001_dp, -0.51542139_dp, 2.49143854_dp], [3, 8])
      case default
         ref = reshape([ &
                       0.00000000_dp, 0.00000000_dp, 2.27734076_dp, 0.00000000_dp, -0.00000000_dp, &
                       -0.00009950_dp, 0.43233630_dp, -0.00000000_dp, 1.39638573_dp, -0.00000000_dp, &
                       1.19203529_dp, -0.75433519_dp, 0.00000000_dp, -1.19203529_dp, -0.75433519_dp, &
                       -0.00000001_dp, -0.55760912_dp, 2.51355017_dp, -0.43233630_dp, 0.00000000_dp, &
                       1.39638578_dp, -0.00000000_dp, 0.55760912_dp, 2.51355016_dp], [3, 8])
      end select
   end subroutine reference_centroids

end module test_mqc_czt_localize

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_localize, only: collect_mqc_czt_localize_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_localize", collect_mqc_czt_localize_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
