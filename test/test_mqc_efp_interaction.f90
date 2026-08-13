!! Electrostatic interaction between two fragments, against GAMESS's own answer
module test_mqc_efp_interaction
   !! **The references here came out of GAMESS.** Two copies of a water potential
   !! 3.0 Angstrom apart along x, and GAMESS's reported electrostatic energy for
   !! that dimer, read off `tools/efp_validation/dimer_energy.py`. What makes them
   !! usable without GAMESS installed is `zero_sections.py`: with the higher
   !! multipoles zeroed out of the potential, GAMESS's electrostatic energy *is* the
   !! truncated sum, so each rank has a reference of its own.
   !!
   !! And what makes those references reachable from a test is that `max_rank`
   !! truncates the same way zeroing the file did. A potential with its quadrupoles
   !! and octupoles zeroed, summed over every rank, equals the full potential summed
   !! through rank one -- the monopoles and dipoles are untouched either way. So the
   !! numbers below are GAMESS's, obtained on modified potentials, checked here
   !! against an unmodified one.
   !!
   !!   rank 0, charges only            0.005641619
   !!   rank 1, through the dipole      0.003913482
   !!   rank 2, through the quadrupole  0.005680360
   !!   rank 3, through the octupole    0.005736532
   !!
   !! Both are quoted to GAMESS's printed nine decimals, which is why the tolerance
   !! is 1e-9 and not tighter: that is the precision of the reference, not of the
   !! arithmetic.
   !!
   !! **The sign of the charge-dipole term is why this test exists.** The two
   !! possible signs differ by 4e-03 Hartree here -- not subtle, but nothing inside
   !! the program would notice, since both give a plausible interaction energy for a
   !! pair of water molecules. The first implementation had it backwards.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_efp_potential, only: efp_potential_t, make_efp_potential, &
                                write_efp_potential
   use mqc_efp_read, only: efp_fragment_t, read_efp_potential
   use mqc_efp_interaction, only: efp_system_t, build_efp_system, &
                                  electrostatic_energy
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_efp_interaction_tests

   !> GAMESS's Bohr, as everywhere else that compares against it.
   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp

   !> GAMESS's electrostatic energy for the dimer, per truncation rank.
   real(dp), parameter :: E_RANK0 = 0.005641619_dp
   real(dp), parameter :: E_RANK1 = 0.003913482_dp
   real(dp), parameter :: E_RANK2 = 0.005680360_dp
   real(dp), parameter :: E_RANK3 = 0.005736532_dp

   !> The references carry nine decimals.
   real(dp), parameter :: REF_TOL = 1.0e-9_dp

   !> The separation `dimer_energy.py` used, along x.
   real(dp), parameter :: SEPARATION = 3.0_dp

contains

   subroutine collect_mqc_efp_interaction_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("efp_charge_charge_vs_gamess", test_rank0), &
                  new_unittest("efp_through_dipole_vs_gamess", test_rank1), &
                  new_unittest("efp_through_quadrupole_vs_gamess", test_rank2), &
                  new_unittest("efp_through_octupole_vs_gamess", test_rank3), &
                  new_unittest("efp_translation_invariance", test_translation), &
                  new_unittest("efp_no_self_interaction", test_no_self) &
                  ]
   end subroutine collect_mqc_efp_interaction_tests

   subroutine dimer(system, err, shift)
      !! Two copies of one water potential, the second displaced along x
      type(efp_system_t), intent(out) :: system
      type(error_t), intent(inout) :: err
      real(dp), intent(in), optional :: shift
         !! Move *both* fragments by this much, which must not change anything.

      type(efp_potential_t) :: pot
      type(efp_fragment_t) :: frags(2)
      real(dp) :: c(3, 3), translations(3, 2)
      real(dp) :: offset
      integer :: z(3)
      character(len=2) :: symbols(3)
      character(len=*), parameter :: path = "test_efp_interaction.efp"

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      ! validation/inputs/sample_inputs/w1.xyz, which is the geometry the GAMESS
      ! reference energies above were generated from. The multipoles depend on it,
      ! so a different water -- even the one the other checks here use -- gives a
      ! different energy and the references stop meaning anything.
      c = reshape([0.00000000000000_dp, 0.00000000009155_dp, 0.10077199490609_dp, &
                   0.00000000000000_dp, 0.77250895271063_dp, -0.46780199741728_dp, &
                   0.00000000000000_dp, -0.77250895280218_dp, -0.46780199748881_dp], &
                  [3, 3])*ANG

      call make_efp_potential(z, symbols, c, "6-31g*", "WATER", pot, err)
      if (err%has_error()) return
      call write_efp_potential(pot, path, err)
      if (err%has_error()) return
      call read_efp_potential(path, frags(1), err)
      if (err%has_error()) return
      call read_efp_potential(path, frags(2), err)
      if (err%has_error()) return

      offset = 0.0_dp
      if (present(shift)) offset = shift
      translations = 0.0_dp
      translations(1, 1) = offset
      translations(1, 2) = SEPARATION*ANG + offset

      call build_efp_system(frags, translations, system, err)
      call frags(1)%destroy()
      call frags(2)%destroy()
      call pot%destroy()
      call delete(path)
   end subroutine dimer

   subroutine test_rank0(error)
      type(error_type), allocatable, intent(out) :: error

      type(efp_system_t) :: system
      type(error_t) :: err

      call dimer(system, err)
      call check(error,.not. err%has_error(), "building the dimer failed")
      if (allocated(error)) return
      call check(error, electrostatic_energy(system, 0), E_RANK0, thr=REF_TOL, &
                 message="the charge-charge energy disagrees with GAMESS")
      call system%destroy()
   end subroutine test_rank0

   subroutine test_rank1(error)
      type(error_type), allocatable, intent(out) :: error

      type(efp_system_t) :: system
      type(error_t) :: err

      call dimer(system, err)
      call check(error,.not. err%has_error(), "building the dimer failed")
      if (allocated(error)) return
      call check(error, electrostatic_energy(system, 1), E_RANK1, thr=REF_TOL, &
                 message="the energy through the dipole disagrees with GAMESS")
      call system%destroy()
   end subroutine test_rank1

   subroutine test_rank2(error)
      !! Three new terms at once -- charge-quadrupole, dipole-quadrupole and
      !! quadrupole-quadrupole -- so this one number is a weaker test than the
      !! lower rungs. Each factor was pinned separately against a potential with
      !! the other ranks zeroed, and this checks the three of them together.
      type(error_type), allocatable, intent(out) :: error

      type(efp_system_t) :: system
      type(error_t) :: err

      call dimer(system, err)
      call check(error,.not. err%has_error(), "building the dimer failed")
      if (allocated(error)) return
      call check(error, electrostatic_energy(system, 2), E_RANK2, thr=REF_TOL, &
                 message="the energy through the quadrupole disagrees with GAMESS")
      call system%destroy()
   end subroutine test_rank2

   subroutine test_rank3(error)
      !! The whole undamped multipole expansion, against GAMESS with its screening
      !! sections removed. Rank three adds charge-octupole alone: GAMESS's total
      !! carries no dipole-octupole, quadrupole-octupole or octupole-octupole term,
      !! which came from reading `FFELEC` rather than from fitting -- an energy
      !! cannot tell a missing term from a small one.
      type(error_type), allocatable, intent(out) :: error

      type(efp_system_t) :: system
      type(error_t) :: err

      call dimer(system, err)
      call check(error,.not. err%has_error(), "building the dimer failed")
      if (allocated(error)) return
      call check(error, electrostatic_energy(system, 3), E_RANK3, thr=REF_TOL, &
                 message="the full multipole energy disagrees with GAMESS")
      call system%destroy()
   end subroutine test_rank3

   subroutine test_translation(error)
      !! Moving the whole system must not change its internal energy
      !!
      !! Cheap, and it catches a class of mistake the reference numbers cannot: an
      !! energy that depends on absolute position rather than on separation. That
      !! would still match GAMESS at this one geometry and be wrong everywhere else.
      type(error_type), allocatable, intent(out) :: error

      type(efp_system_t) :: here, there
      type(error_t) :: err

      call dimer(here, err)
      call check(error,.not. err%has_error(), "building the dimer failed")
      if (allocated(error)) return
      call dimer(there, err, shift=17.0_dp)
      call check(error,.not. err%has_error(), "building the shifted dimer failed")
      if (allocated(error)) return

      call check(error, electrostatic_energy(there, 1), &
                 electrostatic_energy(here, 1), thr=1.0e-12_dp, &
                 message="the energy changed when the whole system moved")
      call here%destroy()
      call there%destroy()
   end subroutine test_translation

   subroutine test_no_self(error)
      !! A lone fragment interacts with nothing
      !!
      !! The pair loop skips same-fragment pairs, and if it did not, the answer
      !! would be enormous and obviously wrong -- but it would be enormous in the
      !! dimer too, where it might be mistaken for a units problem. One fragment
      !! localizes it.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(efp_fragment_t) :: frag(1)
      type(efp_system_t) :: system
      type(error_t) :: err
      real(dp) :: c(3, 3), translations(3, 1)
      integer :: z(3)
      character(len=2) :: symbols(3)
      character(len=*), parameter :: path = "test_efp_single.efp"

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      ! validation/inputs/sample_inputs/w1.xyz, which is the geometry the GAMESS
      ! reference energies above were generated from. The multipoles depend on it,
      ! so a different water -- even the one the other checks here use -- gives a
      ! different energy and the references stop meaning anything.
      c = reshape([0.00000000000000_dp, 0.00000000009155_dp, 0.10077199490609_dp, &
                   0.00000000000000_dp, 0.77250895271063_dp, -0.46780199741728_dp, &
                   0.00000000000000_dp, -0.77250895280218_dp, -0.46780199748881_dp], &
                  [3, 3])*ANG
      call make_efp_potential(z, symbols, c, "6-31g*", "WATER", pot, err)
      call check(error,.not. err%has_error(), "building the potential failed")
      if (allocated(error)) return
      call write_efp_potential(pot, path, err)
      call read_efp_potential(path, frag(1), err)
      call check(error,.not. err%has_error(), "reading the potential failed")
      if (allocated(error)) return

      translations = 0.0_dp
      call build_efp_system(frag, translations, system, err)
      call check(error,.not. err%has_error(), "building the one-fragment system failed")
      if (allocated(error)) return
      call check(error, electrostatic_energy(system, 3), 0.0_dp, thr=1.0e-14_dp, &
                 message="a lone fragment interacted with itself")

      call frag(1)%destroy()
      call system%destroy()
      call pot%destroy()
      call delete(path)
   end subroutine test_no_self

   subroutine delete(path)
      character(len=*), intent(in) :: path

      integer :: unit, stat

      open (newunit=unit, file=path, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete")
   end subroutine delete

end module test_mqc_efp_interaction

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_efp_interaction, only: collect_mqc_efp_interaction_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_efp_interaction", collect_mqc_efp_interaction_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
