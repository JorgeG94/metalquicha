!! Where a molecule reacts, by difference in the electron count
module test_mqc_fukui
   !! The condensed Fukui indices are differences of atomic charges between the
   !! neutral molecule and its two ions, so what can be asserted without a
   !! reference number is quite a lot: they sum to one because exactly one
   !! electron moved, the ionisation potential is positive because removing an
   !! electron costs energy, and water's oxygen is the nucleophilic site because
   !! that is where the lone pairs are.
   !!
   !! Reference values are deliberately not used. Condensed Fukui indices depend
   !! on the population scheme and on the basis, so a stored number would pin
   !! this implementation to one combination and say nothing about whether the
   !! quantity is right.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_fukui, only: fukui_result_t, fukui_indices
   implicit none
   private

   public :: collect_mqc_fukui_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   !> Bohr, C2v, with the two hydrogens exactly equivalent
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp], [3, 3])

contains

   subroutine collect_mqc_fukui_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("indices_sum_to_one_electron", indices_sum_to_one_electron), &
                  new_unittest("water_is_nucleophilic_at_oxygen", water_is_nucleophilic_at_oxygen), &
                  new_unittest("refuses_an_open_shell_neutral", refuses_an_open_shell_neutral) &
                  ]
   end subroutine collect_mqc_fukui_tests

   subroutine water_fukui(scheme, res, err, ok)
      !! Converge water and run the analysis on it
      character(len=*), intent(in) :: scheme
      type(fukui_result_t), intent(out) :: res
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 100, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if
      call fukui_indices(mol, 10, 1, scf%density, scf%energy, scheme, 100, &
                         1.0e-10_dp, 1.0e-8_dp, res, err)
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine water_fukui

   subroutine indices_sum_to_one_electron(error)
      !! Exactly one electron was added, and exactly one removed
      !!
      !! The population scheme divides the whole density, so the difference
      !! between two states differing by one electron has to account for one
      !! electron. That holds for either scheme and for any molecule, and it is
      !! what catches an ion run at the wrong charge, a density that was only
      !! one spin, or a fit that lost part of the molecule.
      type(error_type), allocatable, intent(out) :: error

      type(fukui_result_t) :: res
      type(error_t) :: err
      logical :: ok

      call water_fukui("chelpg", res, err, ok)
      call check(error, ok, "the analysis did not run: "//err%get_message())
      if (allocated(error)) return

      call check(error, sum(res%f_plus), 1.0_dp, thr=1.0e-6_dp, &
                 more="f+ does not account for the electron that was added")
      if (allocated(error)) return
      call check(error, sum(res%f_minus), 1.0_dp, thr=1.0e-6_dp, &
                 more="f- does not account for the electron that was removed")
      if (allocated(error)) return
      ! The dual descriptor is a difference of two things that each sum to one.
      call check(error, sum(res%dual), 0.0_dp, thr=1.0e-6_dp, &
                 more="the dual descriptor should sum to zero")
   end subroutine indices_sum_to_one_electron

   subroutine water_is_nucleophilic_at_oxygen(error)
      !! The chemistry, which is the only thing worth checking without a reference
      !!
      !! Water gives up charge from the oxygen lone pairs and accepts it at the
      !! hydrogens, so the dual descriptor is negative on the oxygen and
      !! positive on both hydrogens. Any implementation that had `f+` and `f-`
      !! the wrong way round would reproduce every sum rule above and get this
      !! exactly backwards.
      !!
      !! The two hydrogens are symmetry-equivalent in this geometry, so they
      !! must agree -- loosely, since the electrostatic fit samples a grid that
      !! does not share the molecule's symmetry.
      type(error_type), allocatable, intent(out) :: error

      type(fukui_result_t) :: res
      type(error_t) :: err
      logical :: ok

      call water_fukui("chelpg", res, err, ok)
      call check(error, ok, "the analysis did not run: "//err%get_message())
      if (allocated(error)) return

      call check(error, res%dual(1) < 0.0_dp, &
                 "water's oxygen should be the nucleophilic site")
      if (allocated(error)) return
      call check(error, res%dual(2) > 0.0_dp .and. res%dual(3) > 0.0_dp, &
                 "water's hydrogens should be the electrophilic sites")
      if (allocated(error)) return
      call check(error, res%f_minus(1) > res%f_minus(2), &
                 "the oxygen should give up charge more readily than a hydrogen")
      if (allocated(error)) return
      call check(error, abs(res%dual(2) - res%dual(3)) < 0.05_dp, &
                 "the two equivalent hydrogens disagree")
      if (allocated(error)) return

      ! Removing an electron costs energy. Nothing else here would notice the
      ! two ion energies being swapped.
      call check(error, res%ionisation_potential > 0.0_dp, &
                 "the ionisation potential came out negative")
   end subroutine water_is_nucleophilic_at_oxygen

   subroutine refuses_an_open_shell_neutral(error)
      !! An open-shell neutral is refused rather than guessed at
      !!
      !! Both ions of a closed shell are doublets. For anything else the two
      !! multiplicities are a chemical choice, and picking one silently would
      !! produce numbers that look exactly like the ones above.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(fukui_result_t) :: res
      type(error_t) :: err

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call run_libcint_rhf(mol, 10, 100, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
      call check(error,.not. err%has_error(), "could not converge water")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call fukui_indices(mol, 10, 3, scf%density, scf%energy, "chelpg", 100, &
                         1.0e-10_dp, 1.0e-8_dp, res, err)
      call check(error, err%has_error(), "an open-shell neutral was accepted")
      call mol%destroy()
   end subroutine refuses_an_open_shell_neutral

end module test_mqc_fukui

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fukui, only: collect_mqc_fukui_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fukui", collect_mqc_fukui_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
