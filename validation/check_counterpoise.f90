!! What the counterpoise correction is actually worth
!!
!!     cmake --build build --target check_counterpoise
!!     ./build/check_counterpoise
!!
!! A measurement rather than a test: the numbers move with basis and separation
!! by design, so there is nothing here to assert. What it shows is the size of
!! the thing, which is the argument for doing the correction at all.
!!
!! Basis-set superposition error is monomer A borrowing B's basis functions and
!! so describing itself better than its own basis allows. In a plain expansion
!! that borrowing happens only in the pair term, so the pair looks more bound
!! than it is and the error survives into the total. Computing each monomer in
!! the *pair's* basis -- ghost centres on the partner's atoms -- puts the same
!! borrowing on both sides of the difference, where it cancels.
!!
!! The 3 Angstrom row is checkable against something independent: SAPT's
!! counterpoise-corrected supermolecular Hartree-Fock on the same dimer and
!! basis, which reaches the number through the dimer-centred basis rather than
!! through any of this, reports 0.009315543671.
program check_counterpoise
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_error, only: error_t
   use pic_timer, only: timer_type
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_convergence_report, only: convergence_header, convergence_footer
   implicit none

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   integer, parameter :: N_ATOMS = 6            !! Two waters
   integer, parameter :: N_MONOMER = 3          !! Atoms in one of them
   integer, parameter :: N_SEP = 3              !! Separations scanned
   integer, parameter :: NELEC_MONOMER = 10
   integer, parameter :: NELEC_DIMER = 20
   integer, parameter :: MAX_ITER = 200
   real(dp), parameter :: E_TOL = 1.0e-11_dp
   real(dp), parameter :: D_TOL = 1.0e-9_dp
   character(len=*), parameter :: BASIS = "6-31g"

   real(dp), parameter :: SEPS(N_SEP) = [3.0_dp, 4.0_dp, 6.0_dp]
   logical, parameter :: GHOST_B(N_ATOMS) = &
                         [.false., .false., .false., .true., .true., .true.]

   integer, parameter :: TABLE_WIDTH = 78

   type(timer_type) :: total_clock, row_clock
   integer :: z(N_ATOMS)
   character(len=2) :: sym(N_ATOMS)
   real(dp) :: c(3, N_ATOMS)
   integer :: k

   z = [8, 1, 1, 8, 1, 1]
   sym = ["O ", "H ", "H ", "O ", "H ", "H "]

   call convergence_header(.true., "Counterpoise on a water dimer, "//BASIS//", RHF", &
                           "    R/Ang     BSSE/monomer         raw dE"// &
                           "     CP-corrected dE      wall/s", TABLE_WIDTH)
   call total_clock%start()
   do k = 1, N_SEP
      call one_separation(SEPS(k))
   end do
   call total_clock%stop()
   call logger%info("  "//repeat("-", TABLE_WIDTH))
   call logger%info("  three SCFs per row, nine in "// &
                    trim(seconds(total_clock%get_elapsed_time())))
   call logger%info("")
   call logger%info("  BSSE falls away with overlap, so the two dE columns must")
   call logger%info("  converge as R grows -- and the correction must not invent")
   call logger%info("  anything where there is no overlap left to correct.")

contains

   subroutine one_separation(sep)
      !! One row: the BSSE on a monomer, and the interaction energy both ways
      real(dp), intent(in) :: sep

      type(libcint_molecule_t) :: alone, in_pair, dimer
      type(rhf_result_t) :: scf_alone, scf_pair, scf_dimer
      type(error_t) :: err
      real(dp) :: bsse, raw, corrected
      character(len=128) :: line
      integer :: i

      c(:, 1) = [0.0_dp, 0.0_dp, 0.10077199_dp]
      c(:, 2) = [0.0_dp, 0.77250895_dp, -0.46780200_dp]
      c(:, 3) = [0.0_dp, -0.77250895_dp, -0.46780200_dp]
      do i = 1, N_MONOMER
         c(:, i + N_MONOMER) = c(:, i)
         c(1, i + N_MONOMER) = c(1, i) + sep
      end do
      c = c*ANG
      call row_clock%start()

      ! Monomer A in its own basis
      call build_libcint_molecule(z(1:N_MONOMER), sym(1:N_MONOMER), &
                                  c(:, 1:N_MONOMER), BASIS, alone, err)
      call run_libcint_rhf(alone, NELEC_MONOMER, MAX_ITER, E_TOL, D_TOL, &
                           .false., scf_alone, err)
      if (bailed(err, sep)) return

      ! The same monomer, same electrons, in the pair's basis
      call build_libcint_molecule(z, sym, c, BASIS, in_pair, err, ghost=GHOST_B)
      call run_libcint_rhf(in_pair, NELEC_MONOMER, MAX_ITER, E_TOL, D_TOL, &
                           .false., scf_pair, err)
      if (bailed(err, sep)) return

      call build_libcint_molecule(z, sym, c, BASIS, dimer, err)
      call run_libcint_rhf(dimer, NELEC_DIMER, MAX_ITER, E_TOL, D_TOL, &
                           .false., scf_dimer, err)
      if (bailed(err, sep)) return

      ! Symmetric dimer, so both monomers carry the same BSSE and one suffices.
      bsse = scf_alone%energy - scf_pair%energy
      raw = scf_dimer%energy - 2.0_dp*scf_alone%energy
      corrected = scf_dimer%energy - 2.0_dp*scf_pair%energy
      call row_clock%stop()
      write (line, "(F9.2,3ES19.8,A12)") sep, bsse, raw, corrected, &
         trim(seconds(row_clock%get_elapsed_time()))
      call logger%info("  "//trim(line))

      call alone%destroy()
      call in_pair%destroy()
      call dimer%destroy()
   end subroutine one_separation

   function bailed(err, sep) result(did_fail)
      !! Report and skip a row rather than stopping the scan
      type(error_t), intent(inout) :: err
      real(dp), intent(in) :: sep
      logical :: did_fail

      did_fail = err%has_error()
      if (did_fail) then
         call logger%info("  "//trim(to_char(sep))//"   failed: "//err%get_message())
         call err%clear()
      end if
   end function bailed

   function seconds(t) result(text)
      !! A wall time at a fixed width, so the column stays a column
      real(dp), intent(in) :: t
      character(len=16) :: text

      write (text, "(F10.2,A)") t, " s"
      text = adjustl(text)
   end function seconds

end program check_counterpoise
