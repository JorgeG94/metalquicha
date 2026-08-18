!! CASCI on real molecules, against PySCF
module test_mqc_casci
   !! Where the CI machinery meets a molecule. Everything upstream of this --
   !! strings, sigma, Davidson -- is validated on model Hamiltonians, which
   !! isolates it but cannot catch a mistake in the step that produces the
   !! integrals those models stand in for. That step is here.
   !!
   !! **Three numbers are checked, not one.** A CASCI total energy is the
   !! inactive constant plus the CI eigenvalue, and the constant is the larger
   !! of the two by an order of magnitude -- for water in cc-pVDZ it is -69.85
   !! against -6.18. So a total energy that agrees says almost nothing about the
   !! CI and a great deal about the inactive Fock, and the two have to be
   !! separated to mean anything. They are compared independently below.
   !!
   !! **But not to the same tolerance, and the reason is worth stating.** The
   !! total agrees with PySCF to 1e-12 while the split agrees only to 1e-9, and
   !! that is not sloppiness in the split -- it is what the two quantities are.
   !! Our orbitals differ from PySCF's by about 6e-10, measured as the
   !! off-diagonal of the overlap between the two sets. The constant is
   !! first-order sensitive to that: moving an occupied orbital across the
   !! inactive-active boundary changes it by of order its orbital energy. The
   !! total sees the same rotation only through the CASCI orbital gradient,
   !! which at a Hartree-Fock solution is small -- and making it zero is exactly
   !! what CASSCF will do next. So the tolerances differ by three orders of
   !! magnitude on purpose, and the loose one still catches an inactive Fock
   !! that is actually wrong by a wide margin.
   !!
   !! The references are taken *before* `mc.kernel()`, from `get_h1eff()` on the
   !! uncanonicalised orbitals. PySCF canonicalises inside `kernel`, which
   !! rotates orbitals within the inactive and active blocks and moves the split
   !! by 1e-8 -- ten times the disagreement being measured. Taken afterwards, as
   !! they were at first, they measure PySCF's canonicalisation.
   !!
   !! The SCF energy is checked first for the same reason: the active-space
   !! Hamiltonian is built from the SCF's orbitals, so orbitals that disagree
   !! with PySCF's would make every later comparison a statement about the SCF.
   !!
   !! References from `pyscf.mcscf.CASCI` on canonical RHF orbitals, driven with
   !! this repository's own basis JSON rather than PySCF's internal tables --
   !! those differ in the eighth decimal on some sets, which looks exactly like
   !! a bug in whichever code is being checked.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_casci, only: run_libcint_casci, casci_result_t, active_space_integrals
   implicit none
   private

   public :: collect_mqc_casci_tests

   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [3, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]

   real(dp), parameter :: NITROGEN(3, 2) = reshape( &
                          [0.0_dp, 0.0_dp, -1.0371_dp, &
                           0.0_dp, 0.0_dp, 1.0371_dp], [3, 2])
   integer, parameter :: NITROGEN_Z(2) = [7, 7]
   character(len=2), parameter :: NITROGEN_SYM(2) = ["N ", "N "]

contains

   subroutine collect_mqc_casci_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("water_cas44", test_water), &
                  new_unittest("nitrogen_cas66", test_nitrogen), &
                  new_unittest("water_full_valence", test_full_valence), &
                  new_unittest("empty_active_space_is_the_scf", test_no_active), &
                  new_unittest("refusals", test_refusals) &
                  ]
   end subroutine collect_mqc_casci_tests

   subroutine one_case(atomic_numbers, symbols, coordinates, basis, n_electrons, &
                       n_inactive, n_active, n_alpha, n_beta, scf_energy, result, err, ok)
      !! An SCF followed by a CASCI on its orbitals
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: n_electrons, n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(out) :: scf_energy
      type(casci_result_t), intent(out) :: result
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf

      ok = .false.
      scf_energy = 0.0_dp
      call build_libcint_molecule(atomic_numbers, symbols, coordinates, basis, mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, n_electrons, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      if (err%has_error() .or. .not. scf%converged) return
      scf_energy = scf%energy

      call run_libcint_casci(mol, scf%orbitals, n_inactive, n_active, n_alpha, n_beta, &
                             result, err, tolerance=1.0e-11_dp)
      ok = .not. err%has_error()
      call mol%destroy()
   end subroutine one_case

   subroutine test_water(error)
      !! Water, cc-pVDZ, CAS(4,4)
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(casci_result_t) :: result
      real(dp) :: scf_energy
      logical :: ok

      real(dp), parameter :: RHF = -76.026781904325_dp
      real(dp), parameter :: CORE = -69.849146671212_dp
      real(dp), parameter :: CI = -6.178181145247_dp
      real(dp), parameter :: TOTAL = -76.027327816459_dp

      call one_case(WATER_Z, WATER_SYM, WATER, "cc-pvdz", 10, 3, 4, 2, 2, &
                    scf_energy, result, err, ok)
      call check(error, ok, "the calculation should succeed")
      if (allocated(error)) return
      call check(error, result%converged, "the CI should converge")
      if (allocated(error)) return

      ! The orbitals first. Everything below is built from them.
      call check(error, scf_energy, RHF, "the RHF energy", thr=1.0e-9_dp)
      if (allocated(error)) return

      call check(error, result%n_determinants, 36, "six alpha strings by six beta")
      if (allocated(error)) return
      call check(error, result%core_energy, CORE, &
                 "the inactive plus nuclear constant", thr=1.0e-7_dp)
      if (allocated(error)) return
      call check(error, result%active_energy, CI, "the CI eigenvalue", thr=1.0e-7_dp)
      if (allocated(error)) return
      call check(error, result%energy, TOTAL, "the total CASCI energy", thr=1.0e-10_dp)
      if (allocated(error)) return

      ! CASCI is variational and contains the reference determinant, so it
      ! cannot be above the SCF. A trivially true statement that is trivially
      ! false if the constant and the eigenvalue have been added wrongly.
      call check(error, result%energy < scf_energy, &
                 "a CASCI containing the Hartree-Fock determinant cannot be "// &
                 "above the Hartree-Fock energy")
   end subroutine test_water

   subroutine test_nitrogen(error)
      !! N2, cc-pVDZ, CAS(6,6) -- the textbook active space
      !!
      !! Six electrons in the three bonding and three antibonding orbitals of
      !! the triple bond. Four hundred determinants, and the correlation
      !! recovered is 68 millihartree against water's 0.5, because a triple bond
      !! is where a single determinant is worst.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(casci_result_t) :: result
      real(dp) :: scf_energy
      logical :: ok

      real(dp), parameter :: RHF = -108.954139046666_dp
      real(dp), parameter :: CORE = -97.547027145109_dp
      real(dp), parameter :: CI = -11.474754117514_dp
      real(dp), parameter :: TOTAL = -109.021781262623_dp

      call one_case(NITROGEN_Z, NITROGEN_SYM, NITROGEN, "cc-pvdz", 14, 4, 6, 3, 3, &
                    scf_energy, result, err, ok)
      call check(error, ok, "the calculation should succeed")
      if (allocated(error)) return
      call check(error, result%converged, "the CI should converge")
      if (allocated(error)) return

      call check(error, scf_energy, RHF, "the RHF energy", thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, result%n_determinants, 400, "twenty alpha strings by twenty")
      if (allocated(error)) return
      call check(error, result%core_energy, CORE, "the constant", thr=1.0e-7_dp)
      if (allocated(error)) return
      call check(error, result%active_energy, CI, "the CI eigenvalue", thr=1.0e-7_dp)
      if (allocated(error)) return
      call check(error, result%energy, TOTAL, "the total", thr=1.0e-10_dp)
      if (allocated(error)) return

      call check(error, scf_energy - result%energy > 0.06_dp, &
                 "CAS(6,6) on a triple bond should recover tens of millihartree, "// &
                 "which is the reason to do it at all")
   end subroutine test_nitrogen

   subroutine test_full_valence(error)
      !! Water, STO-3G, all six valence orbitals active
      !!
      !! Eight electrons in six orbitals, 225 determinants, and unequal to the
      !! other two cases in a way that matters: the active space is nearly the
      !! whole basis, so the inactive block is a single orbital and the constant
      !! is doing much less of the work.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(casci_result_t) :: result
      real(dp) :: scf_energy
      logical :: ok

      real(dp), parameter :: RHF = -74.962986068617_dp
      real(dp), parameter :: CORE = -51.469677153406_dp
      real(dp), parameter :: CI = -23.542755319298_dp
      real(dp), parameter :: TOTAL = -75.012432472704_dp

      call one_case(WATER_Z, WATER_SYM, WATER, "sto-3g", 10, 1, 6, 4, 4, &
                    scf_energy, result, err, ok)
      call check(error, ok, "the calculation should succeed")
      if (allocated(error)) return
      call check(error, result%converged, "the CI should converge")
      if (allocated(error)) return

      call check(error, scf_energy, RHF, "the RHF energy", thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, result%n_determinants, 225, "fifteen strings each way")
      if (allocated(error)) return
      call check(error, result%core_energy, CORE, "the constant", thr=1.0e-7_dp)
      if (allocated(error)) return
      call check(error, result%active_energy, CI, "the CI eigenvalue", thr=1.0e-7_dp)
      if (allocated(error)) return
      call check(error, result%energy, TOTAL, "the total", thr=1.0e-10_dp)
   end subroutine test_full_valence

   subroutine test_no_active(error)
      !! With every occupied orbital inactive, the CASCI energy is the SCF energy
      !!
      !! A limit worth asserting because it ties the two ends together with no
      !! reference at all: the inactive constant is built by a completely
      !! different route from the SCF energy -- a Fock matrix traced against a
      !! density, rather than an eigenvalue sum -- and with no active electrons
      !! it must land on the same number.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(casci_result_t) :: result
      real(dp) :: scf_energy
      logical :: ok

      call one_case(WATER_Z, WATER_SYM, WATER, "cc-pvdz", 10, 5, 0, 0, 0, &
                    scf_energy, result, err, ok)
      call check(error, ok, "the calculation should succeed")
      if (allocated(error)) return
      call check(error, result%n_determinants, 1, "one determinant, and it is empty")
      if (allocated(error)) return
      call check(error, result%core_energy, scf_energy, &
                 "with no active electrons the constant is the whole energy, and "// &
                 "it must equal what the SCF reported", thr=1.0e-10_dp)
      if (allocated(error)) return
      call check(error, result%energy, scf_energy, "and so is the total", &
                 thr=1.0e-10_dp)
   end subroutine test_no_active

   subroutine test_refusals(error)
      !! Active spaces that do not exist
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(casci_result_t) :: result
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: h_eff(:, :), eri_act(:, :, :, :)
      real(dp) :: core

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call run_libcint_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error, scf%converged, "the SCF should converge")
      if (allocated(error)) return

      ! More orbitals than the basis has.
      call active_space_integrals(mol, scf%orbitals, 3, 20, h_eff, eri_act, core, err)
      call check(error, err%has_error(), &
                 "an active space larger than the basis should be refused")
      if (allocated(error)) return
      call err%clear()

      ! More electrons than the active orbitals can hold.
      call run_libcint_casci(mol, scf%orbitals, 1, 2, 3, 3, result, err)
      call check(error, err%has_error(), &
                 "three electrons of one spin in two orbitals should be refused")
      call err%clear()
      call mol%destroy()
   end subroutine test_refusals

end module test_mqc_casci

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_casci, only: collect_mqc_casci_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_casci", collect_mqc_casci_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
