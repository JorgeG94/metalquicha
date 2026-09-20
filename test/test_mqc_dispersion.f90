!! The empirical dispersion wiring: names, refusals, and the gradient's sign
module test_mqc_dispersion
   !! Two halves, and only the first runs on every build.
   !!
   !! The name mapping and the refusals are facts about two vocabularies and are
   !! compiled whether or not s-dftd3 was linked, so they are checked
   !! unconditionally. The energy and the gradient need the library; where it is
   !! absent those cases assert instead that the build says so, naming the
   !! option, rather than returning a zero that would read as "no dispersion
   !! here".
   !!
   !! The gradient case is the one that earns its keep. A dispersion energy is a
   !! few millihartree and a transposed index or a flipped sign is invisible in a
   !! total; central differences of the dispersion energy *alone*, on a geometry
   !! with nothing stationary about it, is what catches either. The pattern is
   !! the one in test_mqc_czt_soscf.f90.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp, default_int
   use mqc_error, only: error_t
   use mqc_dispersion_names, only: d3_functional_alias, dispersion_kind_is_known
   use mqc_dispersion, only: dispersion_available, dispersion_correction
   implicit none
   private

   public :: collect_dispersion_tests

   real(dp), parameter :: REFERENCE_ENERGY = -0.0019027135280896_dp
      !! s-dftd3's own answer for `methane_dimer` below with B3LYP-D3(BJ), read
      !! off its Python bindings rather than from this program.

contains

   subroutine collect_dispersion_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("functional_aliases", test_aliases), &
                  new_unittest("refusals", test_refusals), &
                  new_unittest("kinds", test_kinds), &
                  new_unittest("energy_is_attractive", test_energy), &
                  new_unittest("gradient_matches_finite_difference", test_gradient) &
                  ]
   end subroutine collect_dispersion_tests

   subroutine test_aliases(error)
      !! Our spelling in, s-dftd3's out, including the ones that differ
      type(error_type), allocatable, intent(out) :: error

      character(len=32) :: alias
      type(error_t) :: err

      call d3_functional_alias("b3lyp", alias, err)
      call check(error,.not. err%has_error(), "b3lyp must have D3 parameters")
      if (allocated(error)) return
      call check(error, trim(alias), "b3lyp")
      if (allocated(error)) return

      ! Case is ours to fold: a deck may shout.
      call d3_functional_alias("B3LYP", alias, err)
      call check(error, trim(alias), "b3lyp")
      if (allocated(error)) return

      ! The two vocabularies disagree about hyphens, and this is the whole
      ! reason the mapping exists rather than a pass-through.
      call d3_functional_alias("cam-b3lyp", alias, err)
      call check(error,.not. err%has_error(), "cam-b3lyp must have D3 parameters")
      if (allocated(error)) return
      call check(error, trim(alias), "camb3lyp")
      if (allocated(error)) return

      call d3_functional_alias("m06-l", alias, err)
      call check(error, trim(alias), "m06l")
      if (allocated(error)) return

      call d3_functional_alias("b2gp-plyp", alias, err)
      call check(error, trim(alias), "b2gpplyp")
   end subroutine test_aliases

   subroutine test_refusals(error)
      !! Three ways to have no parameters, each said differently
      type(error_type), allocatable, intent(out) :: error

      character(len=32) :: alias
      type(error_t) :: err

      ! A functional that already carries VV10. Adding D3 on top double counts,
      ! and the message has to say so rather than listing alternatives.
      call d3_functional_alias("wb97x-v", alias, err)
      call check(error, err%has_error(), "wb97x-v plus D3 must be refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "VV10") > 0, &
                 "the -V refusal must say why: it already has non-local correlation")
      if (allocated(error)) return

      ! wB97X without the -V is a different functional with its own D3(BJ)
      ! reparametrisation, and must not be swept up by the refusal above.
      call d3_functional_alias("wb97x", alias, err)
      call check(error,.not. err%has_error(), "wb97x has D3(BJ) parameters of its own")
      if (allocated(error)) return

      ! A functional this program supports that D3 was never fitted for.
      call d3_functional_alias("scan0", alias, err)
      call check(error, err%has_error(), "scan0 has no D3 parameters and must be refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "b3lyp") > 0, &
                 "a refusal must list what is available")
      if (allocated(error)) return

      ! A typo.
      call d3_functional_alias("b3lpy", alias, err)
      call check(error, err%has_error(), "an unknown functional must be refused")
   end subroutine test_refusals

   subroutine test_kinds(error)
      !! Which corrections exist, and what an unknown one does
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: energy
      type(error_t) :: err

      call check(error, dispersion_kind_is_known("d3bj"), "d3bj must be known")
      if (allocated(error)) return
      call check(error,.not. dispersion_kind_is_known("d4"), &
                 "d4 is not wired and must not claim to be")
      if (allocated(error)) return

      ! Refused before anything else is looked at, so this holds on a build with
      ! no library as much as on one with it.
      call dispersion_correction("d4", "b3lyp", [1_default_int, 1_default_int], &
                                 reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.4_dp], [3, 2]), &
                                 energy, error=err)
      call check(error, err%has_error(), "an unwired correction must be refused")
      if (allocated(error)) return
      call check(error, energy == 0.0_dp, "a refused correction must not invent a number")
   end subroutine test_kinds

   subroutine test_energy(error)
      !! A dispersion energy is negative, or the build says it has no library
      type(error_type), allocatable, intent(out) :: error

      integer(default_int) :: numbers(4)
      real(dp) :: coordinates(3, 4), energy
      type(error_t) :: err

      call methane_dimer(numbers, coordinates)

      call dispersion_correction("d3bj", "b3lyp", numbers, coordinates, energy, error=err)

      if (.not. dispersion_available()) then
         call check(error, err%has_error(), &
                    "a build without s-dftd3 must refuse, not return zero")
         if (allocated(error)) return
         call check(error, index(err%get_message(), "MQC_ENABLE_DFTD3") > 0, &
                    "the refusal must name the option that would fix it")
         return
      end if

      call check(error,.not. err%has_error(), "b3lyp-D3(BJ) should evaluate: "// &
                 err%get_message())
      if (allocated(error)) return
      ! London dispersion is attractive between any two atoms, so any geometry
      ! with more than one of them has a negative D3 energy. Weak as a number
      ! and strong as a wiring check: a units or ordering mistake that collapsed
      ! every distance to zero would not survive it.
      call check(error, energy < 0.0_dp, "a dispersion energy must be attractive")
      if (allocated(error)) return
      ! Millihartree, not hartree and not microhartree. A Bohr/Angstrom mix-up
      ! moves this by orders of magnitude.
      call check(error, abs(energy) > 1.0e-6_dp .and. abs(energy) < 1.0e-1_dp, &
                 "a four-atom D3 energy of this size is not physically plausible")
      if (allocated(error)) return

      ! And the number itself, from s-dftd3's own Python bindings on this
      ! geometry -- `RationalDampingParam(method="b3lyp")`, atm left at its
      ! default of False, coordinates in Bohr. Every digit agreed when this was
      ! recorded, which is what says the wiring adds nothing of its own.
      ! A mismatch here is this program's, or a changed pin; it is not a
      ! tolerance to widen.
      call check(error, abs(energy - REFERENCE_ENERGY) < 1.0e-14_dp, &
                 "the D3(BJ) energy does not match s-dftd3's own bindings")
   end subroutine test_energy

   subroutine test_gradient(error)
      !! dE/dR against central differences of E, on a geometry that is not a minimum
      type(error_type), allocatable, intent(out) :: error

      integer(default_int) :: numbers(4)
      real(dp) :: coordinates(3, 4), displaced(3, 4)
      real(dp) :: analytic(3, 4), energy, e_plus, e_minus, numeric, biggest
      real(dp), parameter :: STEP = 1.0e-4_dp
      real(dp), parameter :: TOL = 1.0e-8_dp
      type(error_t) :: err
      integer :: atom, xyz

      if (.not. dispersion_available()) then
         ! Nothing to difference. The refusal is checked in test_energy.
         return
      end if

      call methane_dimer(numbers, coordinates)

      call dispersion_correction("d3bj", "b3lyp", numbers, coordinates, energy, analytic, err)
      call check(error,.not. err%has_error(), "the D3 gradient should evaluate: "// &
                 err%get_message())
      if (allocated(error)) return

      biggest = 0.0_dp
      do atom = 1, 4
         do xyz = 1, 3
            displaced = coordinates
            displaced(xyz, atom) = coordinates(xyz, atom) + STEP
            call dispersion_correction("d3bj", "b3lyp", numbers, displaced, e_plus, error=err)
            displaced(xyz, atom) = coordinates(xyz, atom) - STEP
            call dispersion_correction("d3bj", "b3lyp", numbers, displaced, e_minus, error=err)
            numeric = (e_plus - e_minus)/(2.0_dp*STEP)
            biggest = max(biggest, abs(numeric - analytic(xyz, atom)))
         end do
      end do

      ! A sign error doubles the residual, a transposed index scatters it across
      ! components, and either is thousands of times this tolerance.
      call check(error, biggest < TOL, "the D3 gradient disagrees with central differences")
   end subroutine test_gradient

   pure subroutine methane_dimer(numbers, coordinates)
      !! Two carbons and two hydrogens, deliberately not a stationary geometry
      !!
      !! Four atoms rather than a real dimer: the point is a geometry where
      !! every Cartesian component of the gradient is nonzero and no symmetry
      !! makes a mistake cancel. Bohr, which is what both this program and
      !! s-dftd3 use.
      integer(default_int), intent(out) :: numbers(4)
      real(dp), intent(out) :: coordinates(3, 4)

      numbers = [6_default_int, 1_default_int, 6_default_int, 1_default_int]
      coordinates = reshape([ &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.13_dp, 1.97_dp, 0.41_dp, &
                            6.21_dp, 0.77_dp, -0.53_dp, &
                            7.05_dp, 2.31_dp, 1.09_dp], [3, 4])
   end subroutine methane_dimer

end module test_mqc_dispersion

program tester_dispersion
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_dispersion, only: collect_dispersion_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("dispersion", collect_dispersion_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_dispersion
