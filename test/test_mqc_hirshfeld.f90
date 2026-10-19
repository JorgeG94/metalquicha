module test_mqc_hirshfeld
   !! Pins Hirshfeld charges, spin populations and condensed Fukui indices.
   !!
   !! The reference numbers come from `tools/hirshfeld/hirshfeld_reference.py`,
   !! which follows the same recipe in PySCF -- polarized spherical wB97X free
   !! atoms in def2-TZVP -- but evaluates every free atom exactly at every point
   !! rather than from a radial table, and integrates on PySCF's own 150x974
   !! grid. Agreement therefore checks the table and the quadrature as well as
   !! the partition. What is checked:
   !!
   !!   * a free atom's tabulated density integrates to its electron count;
   !!   * water's charges, RHF/cc-pVDZ, and its cation's charges and spin
   !!     populations, UHF;
   !!   * the Hirshfeld-condensed Fukui indices of water, which are positive
   !!     where Mulliken's need not be, because one partition divides all
   !!     three states;
   !!   * the charges sum to the molecular charge to the accuracy of the grid,
   !!     with nothing renormalised;
   !!   * an element without a free atom is refused.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_constants, only: PI
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, run_czt_uhf
   use mqc_czt_xc, only: xc_available
   use mqc_czt_fukui, only: fukui_result_t, fukui_indices
   use mqc_czt_hirshfeld, only: hirshfeld_populations, free_atom_density, &
                                clear_free_atom_cache
   implicit none
   private
   public :: collect_mqc_hirshfeld_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp], [3, 3])
      !! Bohr, the water of test_mqc_fukui and of the reference script

   real(dp), parameter :: Q_NEUTRAL(3) = [-0.333925335674_dp, 0.166962667733_dp, &
                                          0.166962667733_dp]
   real(dp), parameter :: Q_CATION(3) = [0.335263944926_dp, 0.332368027475_dp, &
                                         0.332368027475_dp]
   real(dp), parameter :: S_CATION(3) = [0.956254823341_dp, 0.021872588392_dp, &
                                         0.021872588392_dp]
   real(dp), parameter :: F_PLUS(3) = [0.304597030610_dp, 0.347701484691_dp, &
                                       0.347701484691_dp]
   real(dp), parameter :: F_MINUS(3) = [0.669189280600_dp, 0.165405359742_dp, &
                                        0.165405359742_dp]

   real(dp), parameter :: CHARGE_TOL = 1.0e-7_dp
      !! Two quadratures, mqc's at its default level against PySCF's 150x974.
      !! Measured at 4e-8 or better, and flat at 1.5e-8 for water from grid
      !! level 4 to 9: what is left is not the molecular quadrature.
   real(dp), parameter :: SUM_TOL = 1.0e-7_dp
      !! How far the charges may sum from the molecular charge, with nothing
      !! renormalised. Measured at 1.4e-8 at the default level and 1.7e-7 at
      !! level 3.

contains

   subroutine collect_mqc_hirshfeld_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("free_atoms_hold_their_electrons", test_free_atoms), &
                  new_unittest("water_matches_reference", test_water), &
                  new_unittest("cation_spins_match_reference", test_cation), &
                  new_unittest("fukui_matches_reference", test_fukui), &
                  new_unittest("heavy_element_is_refused", test_refused), &
                  new_unittest("effective_core_is_refused", test_ecp_refused), &
                  new_unittest("cleared_cache_resolves_identically", test_cache_clear) &
                  ]
   end subroutine collect_mqc_hirshfeld_tests

   subroutine test_free_atoms(error)
      !! 4 pi int r^2 rho(r) dr = Z, on a radial grid uniform in ln r
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 4000
      type(error_t) :: err
      real(dp), allocatable :: r(:), rho(:)
      real(dp) :: h, total
      integer :: k, iz
      integer, parameter :: ELEMENTS(2) = [1, 8]

      if (.not. xc_available()) return
      allocate (r(N))
      h = (log(80.0_dp) - log(1.0e-7_dp))/real(N - 1, dp)
      do k = 1, N
         r(k) = 1.0e-7_dp*exp(real(k - 1, dp)*h)
      end do
      do iz = 1, size(ELEMENTS)
         call free_atom_density(ELEMENTS(iz), r, rho, err)
         call check(error,.not. err%has_error(), "free atom: "//err%get_full_trace())
         if (allocated(error)) return
         ! dr = r d(ln r): trapezoid in ln r, so the integrand carries r^3.
         total = 4.0_dp*PI*h*(sum(r**3*rho) - 0.5_dp*(r(1)**3*rho(1) + r(N)**3*rho(N)))
         call check(error, total, real(ELEMENTS(iz), dp), "free atom electron count", &
                    thr=1.0e-8_dp)
         if (allocated(error)) return
      end do
   end subroutine test_free_atoms

   subroutine water_scf(charge, multiplicity, scf, mol, err)
      integer, intent(in) :: charge, multiplicity
      type(rhf_result_t), intent(out) :: scf
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "cc-pvdz", mol, err)
      if (err%has_error()) return
      if (multiplicity == 1) then
         call run_czt_rhf(mol, 10 - charge, 100, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err)
      else
         call run_czt_uhf(mol, 10 - charge, multiplicity, 100, 1.0e-11_dp, 1.0e-9_dp, &
                          .false., scf, err)
      end if
   end subroutine water_scf

   subroutine test_water(error)
      !! Closed-shell charges, and their sum
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: q(:)
      integer :: i

      if (.not. xc_available()) return
      call water_scf(0, 1, scf, mol, err)
      call hirshfeld_populations(mol, scf%density, q, err)
      call check(error,.not. err%has_error(), "water: "//err%get_full_trace())
      if (allocated(error)) return

      do i = 1, 3
         call check(error, q(i), Q_NEUTRAL(i), "water charge", thr=CHARGE_TOL)
         if (allocated(error)) return
      end do
      call check(error, abs(sum(q)) < SUM_TOL, "water's charges must sum to zero")
      call mol%destroy()
   end subroutine test_water

   subroutine test_cation(error)
      !! Open-shell charges, and spin populations from the same weights
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: q(:), s(:)
      integer :: i

      if (.not. xc_available()) return
      call water_scf(1, 2, scf, mol, err)
      call hirshfeld_populations(mol, scf%density + scf%density_beta, q, err, &
                                 spin_density=scf%density - scf%density_beta, &
                                 spin_populations=s)
      call check(error,.not. err%has_error(), "cation: "//err%get_full_trace())
      if (allocated(error)) return

      do i = 1, 3
         call check(error, q(i), Q_CATION(i), "cation charge", thr=CHARGE_TOL)
         if (allocated(error)) return
         call check(error, s(i), S_CATION(i), "cation spin population", thr=CHARGE_TOL)
         if (allocated(error)) return
      end do
      call check(error, abs(sum(q) - 1.0_dp) < SUM_TOL, "the cation's charges must sum to one")
      if (allocated(error)) return
      call check(error, abs(sum(s) - 1.0_dp) < SUM_TOL, "one unpaired electron")
      call mol%destroy()
   end subroutine test_cation

   subroutine test_fukui(error)
      !! Condensed by one partition, the indices are positive and match
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(fukui_result_t) :: res
      integer :: i

      if (.not. xc_available()) return
      call water_scf(0, 1, scf, mol, err)
      call fukui_indices(mol, 10, 1, scf%density, scf%energy, "hirshfeld", 100, &
                         1.0e-11_dp, 1.0e-9_dp, res, err)
      call check(error,.not. err%has_error(), "fukui: "//err%get_full_trace())
      if (allocated(error)) return

      do i = 1, 3
         call check(error, res%f_plus(i), F_PLUS(i), "f+", thr=CHARGE_TOL)
         if (allocated(error)) return
         call check(error, res%f_minus(i), F_MINUS(i), "f-", thr=CHARGE_TOL)
         if (allocated(error)) return
      end do
      call check(error, all(res%f_plus > 0.0_dp) .and. all(res%f_minus > 0.0_dp), &
                 "Hirshfeld-condensed indices of water are all positive")
      call mol%destroy()
   end subroutine test_fukui

   subroutine test_refused(error)
      !! Xenon has no free-atom configuration, so it has no Hirshfeld share
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp), allocatable :: rho(:)

      call free_atom_density(54, [1.0_dp], rho, err)
      call check(error, err%has_error(), "xenon must be refused")
   end subroutine test_refused

   subroutine test_ecp_refused(error)
      !! An effective core removes electrons an all-electron free atom counts
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: q(:), d(:, :)
      real(dp) :: c(3, 2)

      if (.not. xc_available()) return
      c = 0.0_dp
      c(3, 2) = 3.0_dp
      call build_czt_molecule([1, 53], ["H ", "I "], c, "def2-svp", mol, err, &
                              ecp_name="def2-ecp")
      call check(error,.not. err%has_error(), "HI must build: "//err%get_full_trace())
      if (allocated(error)) return
      allocate (d(mol%nao, mol%nao))
      d = 0.0_dp
      call hirshfeld_populations(mol, d, q, err)
      call check(error, err%has_error(), "an effective core must be refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "effective core") > 0, &
                 "refused for the core, not for iodine being past krypton")
      call mol%destroy()
   end subroutine test_ecp_refused

   subroutine test_cache_clear(error)
      !! Clearing the cache forces a re-solve that reproduces the same table
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp), allocatable :: before(:), after(:)
      real(dp), parameter :: RADII(4) = [0.01_dp, 0.5_dp, 2.0_dp, 6.0_dp]

      if (.not. xc_available()) return
      call free_atom_density(1, RADII, before, err)
      call clear_free_atom_cache()
      call free_atom_density(1, RADII, after, err)
      call check(error,.not. err%has_error(), "re-solve: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, maxval(abs(after - before)) < 1.0e-14_dp*maxval(before), &
                 "a cleared cache must re-solve to the same free atom")
   end subroutine test_cache_clear

end module test_mqc_hirshfeld

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_hirshfeld, only: collect_mqc_hirshfeld_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_hirshfeld", collect_mqc_hirshfeld_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
