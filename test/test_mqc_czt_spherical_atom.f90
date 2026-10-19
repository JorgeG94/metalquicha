module test_mqc_czt_spherical_atom
   !! Pins the spherically averaged free atom against PySCF.
   !!
   !! The reference numbers come from
   !! `tools/spherical_atom/spherical_atom_reference.py`, which reads the
   !! orbital basis from this repository's own JSON so the two codes share both
   !! the basis data and the AO order. Its restricted cases are PySCF's own
   !! `atom_hf.AtomSphericAverageRHF` and `atom_ks.AtomSphericAverageRKS`; its
   !! polarized cases are PySCF's UHF and UKS with the same averaging and a
   !! Hund's-rule occupation, since PySCF has no spin-polarised spherical atom.
   !! What is checked:
   !!
   !!   * Hartree-Fock energies and density fingerprints, restricted and
   !!     polarized, H, C, O, Ne and Fe;
   !!   * Kohn-Sham energies for PBE, B3LYP and the range-separated wB97X,
   !!     restricted and polarized, on a grid fine enough that the quadrature
   !!     no longer limits the agreement;
   !!   * the density is spherical: equal across the components of each column
   !!     pair, nothing between different l, and the configuration's electrons
   !!     in each spin;
   !!   * two limits that need no PySCF: closed-shell neon is the ordinary RHF
   !!     atom, and polarized hydrogen is the ordinary UHF atom;
   !!   * the configuration table, and the refusals.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, subshell_layout
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, run_czt_uhf, SCF_GUESS_CORE
   use mqc_czt_xc, only: xc_available
   use mqc_czt_spherical_atom, only: spherical_atom_settings_t, spherical_atom_t, &
                                     solve_spherical_atom, spherical_atom_electrons
   implicit none
   private
   public :: collect_mqc_czt_spherical_atom_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   real(dp), parameter :: HF_TOL = 1.0e-10_dp
      !! Measured at 7e-13 or better. Both codes converge the commutator to
      !! 1e-10, which matters more than the energy threshold: the reference's
      !! density was only good to 1e-6 while it converged on the energy alone.
   real(dp), parameter :: DENSITY_TOL = 1.0e-9_dp
      !! Measured at 2e-10 or better, on |D| and the first diagonal element.
   integer, parameter :: KS_GRID_LEVEL = 7
      !! Fine enough that what is left against PySCF's 200x1202 grid is no
      !! longer the quadrature. wB97X is the slowest to get there -- its gap
      !! falls 3.7e-6, 2.5e-7, 1.1e-9, 1.7e-11 at levels 2, 3, 5, 7 -- for
      !! the reason `test_mqc_rsh_gradient` records: its semilocal part is
      !! stiff, and PySCF's own wB97X energy still moves with the grid. PBE
      !! and B3LYP are at 1e-13 from level 5.
   real(dp), parameter :: KS_TOL = 1.0e-9_dp
      !! Measured at 2e-11 or better at `KS_GRID_LEVEL`; the densities then
      !! agree as closely as Hartree-Fock's, so `DENSITY_TOL` serves both.

   type :: reference_t
      character(len=2) :: symbol = ""
      integer :: z = 0
      character(len=16) :: basis = ""
      character(len=8) :: functional = ""
      logical :: polarized = .false.
      real(dp) :: energy = 0.0_dp
      real(dp) :: n_alpha = 0.0_dp
      real(dp) :: n_beta = 0.0_dp
      real(dp) :: da11 = 0.0_dp
      real(dp) :: db11 = 0.0_dp
      real(dp) :: norm_total = 0.0_dp
   end type reference_t

contains

   subroutine collect_mqc_czt_spherical_atom_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("hf_restricted_matches_pyscf", test_hf_restricted), &
                  new_unittest("hf_polarized_matches_pyscf", test_hf_polarized), &
                  new_unittest("ks_matches_pyscf", test_ks), &
                  new_unittest("density_is_spherical", test_spherical), &
                  new_unittest("neon_is_the_rhf_atom", test_neon_rhf), &
                  new_unittest("polarized_hydrogen_is_the_uhf_atom", test_hydrogen_uhf), &
                  new_unittest("configuration_table", test_configuration), &
                  new_unittest("refusals", test_refusals) &
                  ]
   end subroutine collect_mqc_czt_spherical_atom_tests

   subroutine one_atom(z, symbol, basis, mol, err)
      integer, intent(in) :: z
      character(len=*), intent(in) :: symbol, basis
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      real(dp) :: c(3, 1)

      c = 0.0_dp
      call build_czt_molecule([z], [symbol], c, basis, mol, err)
   end subroutine one_atom

   function hf_restricted_references() result(refs)
      type(reference_t), allocatable :: refs(:)
      refs = [ &
             reference_t("H ", 1, "cc-pvdz", "", .false., -0.353687564421_dp, &
                         0.500000000000_dp, 0.500000000000_dp, 0.278176484601_dp, &
                         0.278176484601_dp, 0.631195203521_dp), &
             reference_t("C ", 6, "cc-pvdz", "", .false., -37.327421827178_dp, &
                         3.000000000000_dp, 3.000000000000_dp, 0.999463175488_dp, &
                         0.999463175488_dp, 3.144002909999_dp), &
             reference_t("O ", 8, "cc-pvdz", "", .false., -74.269187125634_dp, &
                         4.000000000000_dp, 4.000000000000_dp, 0.999637167996_dp, &
                         0.999637167996_dp, 3.592607801879_dp), &
             reference_t("Ne", 10, "cc-pvdz", "", .false., -128.488775551741_dp, &
                         5.000000000000_dp, 5.000000000000_dp, 0.999999959122_dp, &
                         0.999999959122_dp, 4.472134414124_dp), &
             reference_t("Fe", 26, "def2-svp", "", .false., -1261.250648997709_dp, &
                         13.000000000000_dp, 13.000000000000_dp, 1.118607382283_dp, &
                         1.118607382283_dp, 7.899502045392_dp) &
             ]
   end function hf_restricted_references

   function hf_polarized_references() result(refs)
      type(reference_t), allocatable :: refs(:)
      refs = [ &
             reference_t("H ", 1, "cc-pvdz", "", .true., -0.499278403420_dp, &
                         1.000000000000_dp, 0.000000000000_dp, 1.000001523866_dp, &
                         0.000000000000_dp, 1.000001523867_dp), &
             reference_t("C ", 6, "cc-pvdz", "", .true., -37.523403444387_dp, &
                         4.000000000000_dp, 2.000000000000_dp, 1.001109329638_dp, &
                         0.998434315986_dp, 3.073100139013_dp), &
             reference_t("O ", 8, "cc-pvdz", "", .true., -74.580300440122_dp, &
                         5.000000000000_dp, 3.000000000000_dp, 1.001589015667_dp, &
                         0.998561560864_dp, 3.635412390408_dp), &
             reference_t("Fe", 26, "def2-svp", "", .true., -1261.454485903318_dp, &
                         14.000000000000_dp, 12.000000000000_dp, 1.118609066878_dp, &
                         1.118583222392_dp, 7.894319230741_dp) &
             ]
   end function hf_polarized_references

   function ks_references() result(refs)
      type(reference_t), allocatable :: refs(:)
      refs = [ &
             reference_t("O ", 8, "cc-pvdz", "pbe", .false., -74.912917804104_dp, &
                         4.000000000000_dp, 4.000000000000_dp, 1.004545955854_dp, &
                         1.004545955854_dp, 3.642826168924_dp), &
             reference_t("O ", 8, "cc-pvdz", "pbe", .true., -74.968878340957_dp, &
                         5.000000000000_dp, 3.000000000000_dp, 1.005119331732_dp, &
                         1.004118489991_dp, 3.648361478017_dp), &
             reference_t("O ", 8, "cc-pvdz", "b3lyp", .true., -75.017035557546_dp, &
                         5.000000000000_dp, 3.000000000000_dp, 1.004491674529_dp, &
                         1.003180168352_dp, 3.653910996778_dp), &
             reference_t("O ", 8, "cc-pvdz", "wb97x", .true., -74.924372399935_dp, &
                         5.000000000000_dp, 3.000000000000_dp, 1.004141646589_dp, &
                         1.002517141059_dp, 3.635138674263_dp) &
             ]
   end function ks_references

   subroutine check_reference(error, ref, energy_tol, density_tol)
      !! Solve one reference atom and compare everything the script printed
      type(error_type), allocatable, intent(out) :: error
      type(reference_t), intent(in) :: ref
      real(dp), intent(in) :: energy_tol, density_tol

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(spherical_atom_settings_t) :: settings
      type(spherical_atom_t) :: atom
      real(dp), allocatable :: s(:, :)
      character(len=:), allocatable :: label

      label = trim(ref%symbol)//"/"//trim(ref%basis)
      if (len_trim(ref%functional) > 0) label = label//"/"//trim(ref%functional)
      if (ref%polarized) then
         label = label//" polarized"
      else
         label = label//" restricted"
      end if

      call one_atom(ref%z, ref%symbol, ref%basis, mol, err)
      call check(error,.not. err%has_error(), label//" must build: "//err%get_full_trace())
      if (allocated(error)) return

      settings%polarized = ref%polarized
      settings%energy_tol = 1.0e-12_dp
      settings%commutator_tol = 1.0e-10_dp
      if (len_trim(ref%functional) > 0) settings%functional = trim(ref%functional)
      if (len_trim(ref%functional) > 0) settings%grid_level = KS_GRID_LEVEL
      call solve_spherical_atom(mol, settings, atom, err)
      call check(error,.not. err%has_error(), label//" must solve: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, atom%converged, label//" must converge")
      if (allocated(error)) return

      call check(error, atom%energy, ref%energy, label//" energy", thr=energy_tol)
      if (allocated(error)) return

      call mol%overlap(s)
      call check(error, sum(atom%d_alpha*s), ref%n_alpha, label//" alpha electrons", &
                 thr=1.0e-10_dp)
      if (allocated(error)) return
      call check(error, sum(atom%d_beta*s), ref%n_beta, label//" beta electrons", &
                 thr=1.0e-10_dp)
      if (allocated(error)) return
      call check(error, atom%d_alpha(1, 1), ref%da11, label//" D_alpha(1,1)", thr=density_tol)
      if (allocated(error)) return
      call check(error, atom%d_beta(1, 1), ref%db11, label//" D_beta(1,1)", thr=density_tol)
      if (allocated(error)) return
      call check(error, norm2(atom%d_alpha + atom%d_beta), ref%norm_total, &
                 label//" |D|", thr=density_tol)
      call mol%destroy()
   end subroutine check_reference

   subroutine test_hf_restricted(error)
      !! PySCF's own spherically averaged restricted atom
      type(error_type), allocatable, intent(out) :: error
      type(reference_t), allocatable :: refs(:)
      integer :: i

      refs = hf_restricted_references()
      do i = 1, size(refs)
         call check_reference(error, refs(i), HF_TOL, DENSITY_TOL)
         if (allocated(error)) return
      end do
   end subroutine test_hf_restricted

   subroutine test_hf_polarized(error)
      !! Hund's rule within each l, against PySCF's UHF with the same averaging
      type(error_type), allocatable, intent(out) :: error
      type(reference_t), allocatable :: refs(:)
      integer :: i

      refs = hf_polarized_references()
      do i = 1, size(refs)
         call check_reference(error, refs(i), HF_TOL, DENSITY_TOL)
         if (allocated(error)) return
      end do
   end subroutine test_hf_polarized

   subroutine test_ks(error)
      !! Pure, hybrid and range-separated functionals, through the SCF's own Fock build
      type(error_type), allocatable, intent(out) :: error
      type(reference_t), allocatable :: refs(:)
      integer :: i

      if (.not. xc_available()) return  ! no libxc in this build: nothing to check
      refs = ks_references()
      do i = 1, size(refs)
         call check_reference(error, refs(i), KS_TOL, DENSITY_TOL)
         if (allocated(error)) return
      end do
   end subroutine test_ks

   subroutine test_spherical(error)
      !! Same in every component, nothing between l, the configuration's electrons
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(spherical_atom_settings_t) :: settings
      type(spherical_atom_t) :: atom
      integer, allocatable :: ang(:), first(:), ncomp(:)
      real(dp), allocatable :: d(:, :)
      integer :: n_sub, a, b, k, m
      real(dp) :: worst

      call one_atom(26, "Fe", "def2-svp", mol, err)
      call check(error,.not. err%has_error(), "Fe must build: "//err%get_full_trace())
      if (allocated(error)) return
      call solve_spherical_atom(mol, settings, atom, err)
      call check(error,.not. err%has_error(), "Fe must solve: "//err%get_full_trace())
      if (allocated(error)) return

      call subshell_layout(mol, ang, first, ncomp, n_sub)
      d = atom%d_alpha + atom%d_beta
      worst = 0.0_dp
      do b = 1, n_sub
         do a = 1, n_sub
            associate (blk => d(first(a) + 1:first(a) + ncomp(a), &
                                first(b) + 1:first(b) + ncomp(b)))
               if (ang(a) /= ang(b)) then
                  worst = max(worst, maxval(abs(blk)))
                  cycle
               end if
               do m = 1, ncomp(a)
                  do k = 1, ncomp(a)
                     if (k == m) then
                        worst = max(worst, abs(blk(k, k) - blk(1, 1)))
                     else
                        worst = max(worst, abs(blk(k, m)))
                     end if
                  end do
               end do
            end associate
         end do
      end do
      call check(error, worst < 1.0e-14_dp, "the density must be spherical by construction")
      call mol%destroy()
   end subroutine test_spherical

   subroutine test_neon_rhf(error)
      !! A closed shell has nothing fractional, so this is the ordinary RHF atom
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(spherical_atom_settings_t) :: settings
      type(spherical_atom_t) :: restricted, polarized
      type(rhf_result_t) :: scf

      call one_atom(10, "Ne", "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "Ne must build: "//err%get_full_trace())
      if (allocated(error)) return

      settings%polarized = .false.
      call solve_spherical_atom(mol, settings, restricted, err)
      settings%polarized = .true.
      call solve_spherical_atom(mol, settings, polarized, err)
      call run_czt_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-9_dp, .false., scf, err, &
                       guess=SCF_GUESS_CORE)
      call check(error,.not. err%has_error(), "Ne must solve: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, restricted%energy, scf%energy, "restricted Ne is RHF Ne", &
                 thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, polarized%energy, scf%energy, "polarized Ne is RHF Ne", &
                 thr=1.0e-9_dp)
      call mol%destroy()
   end subroutine test_neon_rhf

   subroutine test_hydrogen_uhf(error)
      !! One electron in one orbital: polarized is exactly the doublet UHF atom
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(spherical_atom_settings_t) :: settings
      type(spherical_atom_t) :: atom
      type(rhf_result_t) :: scf

      call one_atom(1, "H", "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "H must build: "//err%get_full_trace())
      if (allocated(error)) return

      call solve_spherical_atom(mol, settings, atom, err)
      call run_czt_uhf(mol, 1, 2, 100, 1.0e-12_dp, 1.0e-9_dp, .false., scf, err, &
                       guess=SCF_GUESS_CORE)
      call check(error,.not. err%has_error(), "H must solve: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, atom%energy, scf%energy, "polarized H is UHF H", thr=1.0e-9_dp)
      call mol%destroy()
   end subroutine test_hydrogen_uhf

   subroutine test_configuration(error)
      !! Hund's rule within each l, and the restricted halves
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: e(0:3, 2)

      call spherical_atom_electrons(8, .true., e, err)
      call check(error, all(e(:, 1) == [2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp]) .and. &
                 all(e(:, 2) == [2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]), "O polarized")
      if (allocated(error)) return
      call spherical_atom_electrons(8, .false., e, err)
      call check(error, all(e(:, 1) == [2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]) .and. &
                 all(e(:, 2) == e(:, 1)), "O restricted")
      if (allocated(error)) return
      ! NRSRHF's iron: 3d8 and an empty 4s.
      call spherical_atom_electrons(26, .true., e, err)
      call check(error, all(e(:, 1) == [3.0_dp, 6.0_dp, 5.0_dp, 0.0_dp]) .and. &
                 all(e(:, 2) == [3.0_dp, 6.0_dp, 3.0_dp, 0.0_dp]), "Fe polarized")
      if (allocated(error)) return
      call spherical_atom_electrons(1, .false., e, err)
      call check(error, e(0, 1) == 0.5_dp .and. e(0, 2) == 0.5_dp, "H restricted")
      if (allocated(error)) return
      call check(error,.not. err%has_error(), "H-Kr must be in the table")
      if (allocated(error)) return

      call spherical_atom_electrons(54, .true., e, err)
      call check(error, err%has_error(), "xenon is outside the table")
   end subroutine test_configuration

   subroutine test_refusals(error)
      !! Two centres, and a Cartesian d shell, are refused rather than averaged
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(spherical_atom_settings_t) :: settings
      type(spherical_atom_t) :: atom
      real(dp) :: c(3, 2)

      c = 0.0_dp
      c(3, 2) = 0.74_dp*ANG
      call build_czt_molecule([1, 1], ["H ", "H "], c, "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "H2 must build: "//err%get_full_trace())
      if (allocated(error)) return
      call solve_spherical_atom(mol, settings, atom, err)
      call check(error, err%has_error(), "two centres must be refused")
      if (allocated(error)) return
      call mol%destroy()
      call err%clear()

      call build_czt_molecule([8], ["O "], c(:, 1:1), "cc-pvdz", mol, err, &
                              force_cartesian=.true.)
      call check(error,.not. err%has_error(), "Cartesian O must build: "//err%get_full_trace())
      if (allocated(error)) return
      call solve_spherical_atom(mol, settings, atom, err)
      call check(error, err%has_error(), "a Cartesian d shell must be refused")
      call mol%destroy()
   end subroutine test_refusals

end module test_mqc_czt_spherical_atom

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_spherical_atom, only: collect_mqc_czt_spherical_atom_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_spherical_atom", collect_mqc_czt_spherical_atom_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
