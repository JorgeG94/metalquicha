!! Manual check that the multipole integrals are right
!!
!!     cmake -B build -DMQC_ENABLE_LIBCINT=ON
!!     ./build/check_multipole && python3 validation/check_multipole.py
!!
!! Milestone 1 of the EFP plan. Two things are being established, and they fail
!! differently:
!!
!!   * the integrals themselves, elementwise against `mol.intor('int1e_r')` and
!!     friends. A wrong component ordering or a wrong block offset shows up here
!!     and essentially nowhere else -- a transposed quadrupole still contracts to
!!     a plausible number against a density.
!!   * the origin actually reaching libcint. It travels through `env`, whose slot
!!     constants are 0-based and unconverted, and the failure mode is silent: the
!!     moments come back about the wrong point. The check for it is a *shifted*
!!     origin, because about the default origin a wrong slot gives the right
!!     answer by accident.
!!
!! The dipole from `Tr(D r)` is also printed. That number does not currently
!! exist on the CPU path at all -- only the GPU backend computes a dipole -- so
!! this is the first time it can be compared with anything.
program check_multipole
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_multipole, only: multipole_matrices
   use mqc_error, only: error_t
   implicit none

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   integer :: failures

   failures = 0
   call one_case(1, "sto-3g", [0.0_dp, 0.0_dp, 0.0_dp])
   call one_case(2, "cc-pvdz", [0.0_dp, 0.0_dp, 0.0_dp])
   ! A shifted origin, which is the only arrangement that can catch the origin
   ! never arriving: at (0,0,0) a dropped origin is indistinguishable from a
   ! delivered one.
   call one_case(3, "cc-pvdz", [1.3_dp, -0.7_dp, 2.1_dp])

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[multipole] wrote every case; run validation/check_multipole.py"
   else
      write (*, "(A,I0,A)") "[multipole] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine one_case(case_index, basis, origin)
      !! Files are named by index, not by the origin. Writing a float into a
      !! filename invites exactly the mismatch it caused: `F0.1` drops the
      !! leading zero, so 0.0 became ".0" and the Python side looked for a file
      !! that was never written. The basis and origin go *inside* the file.
      integer, intent(in) :: case_index
      character(len=*), intent(in) :: basis
      real(dp), intent(in) :: origin(3)

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: dip(:, :, :), quad(:, :, :), oct(:, :, :)
      real(dp) :: c(3, 3), mu(3), nuc(3)
      integer :: unit, k, iatom
      character(len=64) :: name

      ! Water, off-centre and unsymmetric, so no component is zero by accident.
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])

      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
      if (err%has_error()) then
         write (*, "(A,A)") "[multipole] basis failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call multipole_matrices(mol, origin, 1, dip, err)
      if (.not. err%has_error()) call multipole_matrices(mol, origin, 2, quad, err)
      if (.not. err%has_error()) call multipole_matrices(mol, origin, 3, oct, err)
      if (err%has_error()) then
         write (*, "(A,A)") "[multipole] integrals failed: ", err%get_message()
         failures = failures + 1
         call mol%destroy()
         return
      end if

      call run_libcint_rhf(mol, 10, 100, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err, &
                           in_core=.true.)
      if (err%has_error() .or. .not. scf%converged) then
         write (*, "(A)") "[multipole] SCF failed"
         failures = failures + 1
         call mol%destroy()
         return
      end if

      ! The electronic dipole is -Tr(D r); the nuclei contribute +sum Z (R - O).
      do k = 1, 3
         mu(k) = -sum(scf%density*dip(:, :, k))
      end do
      nuc = 0.0_dp
      do iatom = 1, mol%natm
         do k = 1, 3
            nuc(k) = nuc(k) + mol%charges(iatom)*(mol%coords(k, iatom) - origin(k))
         end do
      end do

      write (name, "(A,I0,A)") "/tmp/mqc_multipole_", case_index, ".txt"
      open (newunit=unit, file=trim(name), status="replace", action="write")
      write (unit, "(A)") trim(basis)
      write (unit, "(I0)") mol%nao
      write (unit, "(3es25.16e3)") origin
      write (unit, "(3es25.16e3)") mu + nuc
      write (unit, "(es25.16e3)") scf%energy
      ! The matrices themselves, for the elementwise comparison.
      write (unit, "(es25.16e3)") dip
      write (unit, "(es25.16e3)") quad
      write (unit, "(es25.16e3)") oct
      ! The overlap too. A basis function that is not exactly normalised is
      ! invisible to the SCF energy -- scaling a function does not change the
      ! space it spans, so the coefficients absorb it -- but it scales every
      ! per-AO quantity, including these moments. The overlap diagonal is where
      ! it shows.
      block
         real(dp), allocatable :: ovl(:, :)
         call mol%overlap(ovl)
         write (unit, "(es25.16e3)") ovl
         deallocate (ovl)
      end block
      close (unit)

      write (*, "(A,A,A,3F7.2,A,I0,A,3F14.8)") "  ", trim(basis), "  origin", origin, &
         "  nao ", mol%nao, "   dipole", mu + nuc
      call mol%destroy()
      deallocate (dip, quad, oct)
   end subroutine one_case

end program check_multipole
