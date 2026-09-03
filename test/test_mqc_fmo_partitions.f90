!! Which partitions an FMO run refuses, and which it must not
!!
!! Fragmenting across a covalent bond is not an approximation, it is a
!! different molecule. The guard that stops it is worth a test of its own
!! because the obvious safety net does not catch the interesting case:
!!
!!   * cut **one** bond and both fragments come back with an odd electron
!!     count, which the closed-shell check refuses on its own
!!   * cut an **even** number per fragment -- a ring, a double bond -- and
!!     every count stays even. Nothing but a connectivity check stands between
!!     that and a confident wrong answer. Cyclopropane split into three CH2
!!     used to come back 0.28 Hartree low: 176 kcal/mol wearing the shape of an
!!     answer
!!
!! And the other direction matters as much. A hydrogen bond is not a covalent
!! one, so a guard tightened until it refuses clusters has been tightened too
!! far -- every water dimer in this suite would stop working.
!!
!! A deck cannot express any of this: a refusal is not a number, so there is
!! nothing for a validation case to compare against. It belongs here, and it
!! costs milliseconds -- the refusals happen before any SCF starts.
module test_mqc_fmo_partitions
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   use mqc_error, only: error_t
   use mqc_physical_constants, only: ANGSTROM_TO_BOHR
   implicit none
   private

   public :: collect_fmo_partitions

   integer, parameter :: N_DIM = 3
   integer, parameter :: N_ETHANE = 8
   integer, parameter :: N_RING = 9
   integer, parameter :: N_DIMER = 6

contains

   subroutine collect_fmo_partitions(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("a_single_cut_bond_is_refused", test_ethane), &
                  new_unittest("an_even_cut_is_refused_too", test_cyclopropane), &
                  new_unittest("a_hydrogen_bond_is_not_a_cut", test_close_dimer) &
                  ]
   end subroutine collect_fmo_partitions

   subroutine test_ethane(error)
      !! Ethane split at the C-C bond
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: ethane(N_DIM, N_ETHANE)

      ethane = reshape([0.000_dp, 0.000_dp, 0.768_dp, &
                        -1.019_dp, 0.000_dp, 1.157_dp, &
                        0.510_dp, 0.883_dp, 1.157_dp, &
                        0.510_dp, -0.883_dp, 1.157_dp, &
                        0.000_dp, 0.000_dp, -0.768_dp, &
                        1.019_dp, 0.000_dp, -1.157_dp, &
                        -0.510_dp, -0.883_dp, -1.157_dp, &
                        -0.510_dp, 0.883_dp, -1.157_dp], [N_DIM, N_ETHANE])

      call must_refuse(error, "ethane cut at C-C", [6, 1, 1, 1, 6, 1, 1, 1], &
                       ethane, [1, 1, 1, 1, 2, 2, 2, 2])
   end subroutine test_ethane

   subroutine test_cyclopropane(error)
      !! Three CH2, eight electrons each -- every count even
      !!
      !! This is the one the electron-count check cannot see. Two bonds are cut
      !! per fragment, so the parity that catches ethane is restored and only
      !! connectivity is left to object.
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: ring(N_DIM, N_RING)

      ring = reshape([0.8718_dp, 0.0000_dp, 0.0000_dp, &
                      1.3818_dp, 0.0000_dp, 0.9500_dp, &
                      1.3818_dp, 0.0000_dp, -0.9500_dp, &
                      -0.4359_dp, 0.7550_dp, 0.0000_dp, &
                      -0.6909_dp, 1.1967_dp, 0.9500_dp, &
                      -0.6909_dp, 1.1967_dp, -0.9500_dp, &
                      -0.4359_dp, -0.7550_dp, 0.0000_dp, &
                      -0.6909_dp, -1.1967_dp, 0.9500_dp, &
                      -0.6909_dp, -1.1967_dp, -0.9500_dp], [N_DIM, N_RING])

      call must_refuse(error, "cyclopropane cut into three CH2", &
                       [6, 1, 1, 6, 1, 1, 6, 1, 1], ring, [1, 1, 1, 2, 2, 2, 3, 3, 3])
   end subroutine test_cyclopropane

   subroutine test_close_dimer(error)
      !! Two waters 2.5 Angstrom apart must pass
      !!
      !! Tighter than any hydrogen bond the rest of the suite uses. A guard that
      !! refused this would refuse every water cluster here, so this is what
      !! keeps the check from being tightened into uselessness.
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: dimer(N_DIM, N_DIMER)

      dimer = reshape([0.0000_dp, 0.0000_dp, 0.0000_dp, &
                       0.0000_dp, -0.7572_dp, 0.5865_dp, &
                       0.0000_dp, 0.7572_dp, 0.5865_dp, &
                       0.0000_dp, 0.0000_dp, 2.5000_dp, &
                       0.0000_dp, -0.7572_dp, 3.0865_dp, &
                       0.0000_dp, 0.7572_dp, 3.0865_dp], [N_DIM, N_DIMER])

      call must_allow(error, "two waters 2.5 A apart", [8, 1, 1, 8, 1, 1], &
                      dimer, [1, 1, 1, 2, 2, 2])
   end subroutine test_close_dimer

   subroutine must_refuse(error, label, z, coords_ang, owner)
      !! This partition must be rejected before anything is computed
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      integer, intent(in) :: z(:), owner(:)
      real(dp), intent(in) :: coords_ang(:, :)

      type(error_t) :: err

      call attempt(z, coords_ang, owner, err)
      call check(error, err%has_error(), &
                 label//" was answered rather than refused")
      call err%clear()
   end subroutine must_refuse

   subroutine must_allow(error, label, z, coords_ang, owner)
      !! This partition must be accepted
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      integer, intent(in) :: z(:), owner(:)
      real(dp), intent(in) :: coords_ang(:, :)

      type(error_t) :: err

      call attempt(z, coords_ang, owner, err)
      call check(error,.not. err%has_error(), &
                 label//" was refused: "//err%get_message())
      call err%clear()
   end subroutine must_allow

   subroutine attempt(z, coords_ang, owner, error)
      !! Run FMO2 on this partition and report whatever came back
      integer, intent(in) :: z(:), owner(:)
      real(dp), intent(in) :: coords_ang(:, :)
      type(error_t), intent(inout) :: error

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      character(len=2), allocatable :: symbols(:)
      real(dp), allocatable :: coords(:, :)
      integer :: i

      allocate (symbols(size(z)))
      do i = 1, size(z)
         select case (z(i))
         case (6)
            symbols(i) = "C "
         case (7)
            symbols(i) = "N "
         case (8)
            symbols(i) = "O "
         case default
            symbols(i) = "H "
         end select
      end do
      coords = coords_ang*ANGSTROM_TO_BOHR

      opts%basis = "sto-3g"
      call run_fmo2(z, symbols, coords, owner, opts, res, error)
   end subroutine attempt

end module test_mqc_fmo_partitions

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fmo_partitions, only: collect_fmo_partitions
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fmo_partitions", collect_fmo_partitions)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
