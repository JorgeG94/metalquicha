!! Manual check that libcint is wired up and gives the right integrals
!!
!!     cmake -B build -DMQC_ENABLE_CZT=ON && ./build/check_libcint
!!
!! H2 in STO-3G, at the textbook geometry, because it has answers that can be
!! checked rather than merely admired: the overlap diagonal is 1 if and only
!! if the contraction coefficients were normalised on the way in, and the
!! off-diagonal is 0.6593 at R = 1.4 bohr, a number that predates all of us.
!!
!! Getting the normalisation wrong is the standard way to misuse libcint: it
!! wants coefficients already multiplied by each primitive's norm, and if they
!! are not, everything still runs and every number is quietly wrong by a
!! constant per shell. The diagonal is what catches it, so it is checked
!! first and separately.
!!
!! The basis is written out by hand here rather than read through the BSE
!! reader. That is deliberate for a first check: it tests libcint and the
!! packing alone, so a failure cannot be the basis reader's. Wiring the real
!! reader in is the next step and has its own way of being wrong.
program check_libcint
   use pic_types, only: dp
   use libcint_fortran, only: libcint_1e_ovlp_sph, libcint_cgto_sph, &
                              libcint_tot_cgto_sph, libcint_gto_norm, &
                              LIBCINT_ATM_SLOTS, LIBCINT_BAS_SLOTS, &
                              LIBCINT_CHARGE_OF, LIBCINT_PTR_COORD, &
                              LIBCINT_ATOM_OF, LIBCINT_ANG_OF, LIBCINT_NPRIM_OF, &
                              LIBCINT_NCTR_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_COEFF, &
                              LIBCINT_PTR_ENV_START
   implicit none

   integer, parameter :: NATM = 2, NBAS = 2, NPRIM = 3
   real(dp), parameter :: R = 1.4_dp          !! bohr, the classic H2 separation
   real(dp), parameter :: S12_REFERENCE = 0.6593_dp
   integer, parameter :: ENV_SIZE = 64
      !! Coordinates, exponents and coefficients for this molecule, past
      !! libcint's reserved header. Ample for two hydrogens.

   !! The slot constants from `libcint_fortran` are already 1-based; the C
   !  header's are 0-based and the Fortran interface converts them. Adding one
   !  here as well writes into the neighbouring slot and leaves the real one
   !  uninitialised, which shows up as a crash inside libcint rather than as a
   !  wrong number.

   !! STO-3G hydrogen, as tabulated
   real(dp), parameter :: EXPONENTS(NPRIM) = &
                          [3.42525091_dp, 0.62391373_dp, 0.16885540_dp]
   real(dp), parameter :: COEFFS(NPRIM) = &
                          [0.15432897_dp, 0.53532814_dp, 0.44463454_dp]

   integer :: atm(LIBCINT_ATM_SLOTS, NATM)
   integer :: bas(LIBCINT_BAS_SLOTS, NBAS)
   real(dp) :: env(ENV_SIZE)
   real(dp), allocatable :: buf(:)
   real(dp) :: s(2, 2)
   integer :: shls(2)
   integer :: off, i, j, di, dj, ret, nao, failures

   failures = 0
   ! Zeroed, not merely filled: libcint reads slots this program never sets --
   ! NUC_MOD_OF, PTR_ZETA, KAPPA_OF -- and stack garbage in a pointer slot is
   ! a segfault inside the library with nothing to say where it came from.
   atm = 0
   bas = 0
   env = 0.0_dp
   off = LIBCINT_PTR_ENV_START

   ! ---- two hydrogens on the z axis ----------------------------------------
   do i = 1, NATM
      atm(LIBCINT_CHARGE_OF, i) = 1
      atm(LIBCINT_PTR_COORD, i) = off
      env(off + 1) = 0.0_dp
      env(off + 2) = 0.0_dp
      env(off + 3) = real(i - 1, dp)*R
      off = off + 3
   end do

   ! ---- one s shell on each ------------------------------------------------
   do i = 1, NBAS
      bas(LIBCINT_ATOM_OF, i) = i - 1        ! libcint counts atoms from 0
      bas(LIBCINT_ANG_OF, i) = 0             ! s
      bas(LIBCINT_NPRIM_OF, i) = NPRIM
      bas(LIBCINT_NCTR_OF, i) = 1
      bas(LIBCINT_PTR_EXP, i) = off
      env(off + 1:off + NPRIM) = EXPONENTS
      off = off + NPRIM
      bas(LIBCINT_PTR_COEFF, i) = off
      ! The normalisation libcint expects to have been applied already.
      do j = 1, NPRIM
         env(off + j) = COEFFS(j)*libcint_gto_norm(0, EXPONENTS(j))
      end do
      off = off + NPRIM
   end do

   nao = libcint_tot_cgto_sph(bas, NBAS)
   call expect(nao == 2, "two basis functions for two hydrogens")

   ! ---- the overlap, shell pair by shell pair ------------------------------
   allocate (buf(16))
   do i = 1, NBAS
      di = libcint_cgto_sph(i - 1, bas)
      do j = 1, NBAS
         dj = libcint_cgto_sph(j - 1, bas)
         shls = [i - 1, j - 1]
         ret = libcint_1e_ovlp_sph(buf, shls, atm, NATM, bas, NBAS, env)
         call expect(ret /= 0, "libcint returned a shell pair")
         s(i, j) = buf(1)
      end do
   end do

   write (*, "(A)") "H2, STO-3G, R = 1.4 bohr"
   write (*, "(A,2F12.8)") "   S row 1: ", s(1, 1), s(1, 2)
   write (*, "(A,2F12.8)") "   S row 2: ", s(2, 1), s(2, 2)

   ! Diagonal first: it is the normalisation, and every other number is wrong
   ! by the same factor if it is off.
   !
   ! 1e-6 rather than machine precision, and the reason is the table not the
   ! code: STO-3G contraction coefficients are published to eight decimals, so
   ! the contraction they define is normalised only to about 1e-8. Demanding
   ! more asks the constants to be exact. It still discriminates by orders of
   ! magnitude -- omitting libcint_gto_norm entirely does not move this
   ! diagonal by 1e-6, it moves it off 1 completely.
   call expect(abs(s(1, 1) - 1.0_dp) < 1.0e-6_dp, &
               "overlap diagonal is 1 -- contraction coefficients normalised")
   call expect(abs(s(2, 2) - 1.0_dp) < 1.0e-6_dp, "second diagonal is 1")
   call expect(abs(s(1, 2) - s(2, 1)) < 1.0e-14_dp, "overlap is symmetric")
   call expect(abs(s(1, 2) - S12_REFERENCE) < 5.0e-5_dp, &
               "off-diagonal matches the tabulated 0.6593")

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A,F10.6,A)") "[libcint] all ok -- S12 = ", s(1, 2), &
         " against a reference of 0.6593"
   else
      write (*, "(A,I0,A)") "[libcint] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine expect(condition, what)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: what

      if (.not. condition) then
         write (*, "(A)") "  FAILED: "//what
         failures = failures + 1
      end if
   end subroutine expect

end program check_libcint
