!! Manual check that Boys localization is right
!!
!!     cmake -B build -DMQC_ENABLE_LIBCINT=ON
!!     ./build/check_localize && python3 validation/check_localize.py
!!
!! Milestone 2 of the EFP plan. Three internal checks that need no reference at
!! all, then a dump for comparison against `pyscf.lo.Boys`:
!!
!!   * **the energy is invariant.** A rotation among the occupied orbitals cannot
!!     change the density, so `Tr(D H) + ...` is untouched. This is the check that
!!     catches a rotation that is not orthogonal, which is the way a Jacobi sweep
!!     usually breaks, and it is exact rather than approximate.
!!   * **the functional rises monotonically.** The derivation says each pair's
!!     gain is `sqrt(P^2+Q^2) - P >= 0`, so no sweep may lower L. A wrong sign in
!!     the angle still converges -- to a stationary point that is a minimum or a
!!     saddle -- and the energy check would not notice, because any orthogonal
!!     rotation preserves it.
!!   * **the orbitals stay orthonormal**, `C^T S C = 1`.
!!
!! Then water's centroids, which should be recognisably two bonds and two lone
!! pairs, and are what the polarizabilities will later be placed on.
program check_localize
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_localize, only: boys_localize
   use mqc_error, only: error_t
   implicit none

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   integer :: failures

   failures = 0
   call one_case(1, "sto-3g")
   call one_case(2, "cc-pvdz")

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[localize] wrote every case; run validation/check_localize.py"
   else
      write (*, "(A,I0,A)") "[localize] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine one_case(case_index, basis)
      integer, intent(in) :: case_index
      character(len=*), intent(in) :: basis

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: loc(:, :), cen(:, :), s(:, :), work(:, :), gram(:, :)
      real(dp), allocatable :: density(:, :)
      real(dp) :: c(3, 3)
      real(dp) :: l_before, l_after, worst_orth, e_before, e_after
      integer :: unit, i, j, k, n_occ, sweeps
      character(len=48) :: name

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])

      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
      if (err%has_error()) then
         write (*, "(A,A)") "[localize] basis failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call run_libcint_rhf(mol, 10, 100, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err, &
                           in_core=.true.)
      if (err%has_error() .or. .not. scf%converged) then
         write (*, "(A)") "[localize] SCF failed"
         failures = failures + 1
         call mol%destroy()
         return
      end if
      n_occ = scf%n_occupied
      e_before = scf%energy

      ! The functional before localizing, from the canonical orbitals.
      call boys_localize(mol, scf%orbitals, n_occ, loc, cen, err, max_sweeps=0, &
                         functional=l_before)
      if (.not. err%has_error()) then
         deallocate (loc, cen)
         call boys_localize(mol, scf%orbitals, n_occ, loc, cen, err, &
                            sweeps_taken=sweeps, functional=l_after)
      end if
      if (err%has_error()) then
         write (*, "(A,A)") "[localize] failed: ", err%get_message()
         failures = failures + 1
         call mol%destroy()
         return
      end if

      ! Orthonormality: C^T S C must be the identity.
      call mol%overlap(s)
      allocate (work(mol%nao, n_occ), gram(n_occ, n_occ))
      work = matmul(s, loc)
      gram = matmul(transpose(loc), work)
      worst_orth = 0.0_dp
      do j = 1, n_occ
         do i = 1, n_occ
            if (i == j) then
               worst_orth = max(worst_orth, abs(gram(i, j) - 1.0_dp))
            else
               worst_orth = max(worst_orth, abs(gram(i, j)))
            end if
         end do
      end do

      ! The density from the localized orbitals must equal the SCF one, which is
      ! the statement that the energy did not move -- and is sharper than
      ! recomputing an energy, since it compares the object the energy is built
      ! from rather than a scalar that could agree by cancellation.
      allocate (density(mol%nao, mol%nao))
      density = 2.0_dp*matmul(loc, transpose(loc))
      e_after = maxval(abs(density - scf%density))

      write (name, "(A,I0,A)") "/tmp/mqc_localize_", case_index, ".txt"
      open (newunit=unit, file=trim(name), status="replace", action="write")
      write (unit, "(A)") trim(basis)
      write (unit, "(I0,1X,I0)") mol%nao, n_occ
      write (unit, "(es25.16e3)") l_before
      write (unit, "(es25.16e3)") l_after
      write (unit, "(3es25.16e3)") cen
      write (unit, "(es25.16e3)") loc
      close (unit)

      write (*, "(A,A,A,I0,A,I0,A,F12.6,A,F12.6,A,ES9.2,A,ES9.2)") &
         "  ", trim(basis), "  n_occ ", n_occ, "  sweeps ", sweeps, &
         "   L ", l_before, " -> ", l_after, &
         "   orth ", worst_orth, "   density ", e_after

      if (l_after < l_before - 1.0e-12_dp) then
         write (*, "(A)") "    FAIL the Boys functional went down"
         failures = failures + 1
      end if
      if (worst_orth > 1.0e-10_dp) then
         write (*, "(A)") "    FAIL the localized orbitals are not orthonormal"
         failures = failures + 1
      end if
      if (e_after > 1.0e-10_dp) then
         write (*, "(A)") "    FAIL localization changed the density"
         failures = failures + 1
      end if

      call mol%destroy()
      deallocate (loc, cen, s, work, gram, density)
   end subroutine one_case

end program check_localize
