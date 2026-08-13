!! The dipole-quadrupole response summed over localized orbitals
!!
!!     cmake -B build -DMQC_ENABLE_LIBCINT=ON
!!     ./build/check_dipquad_sumrule
!!     python3 validation/check_dipquad_sumrule.py
!!
!! **Why a sum.** Every comparison of the `DIPOLE-QUADRUPOLE` block against GAMESS
!! so far has been per localized orbital, and a per-orbital comparison cannot
!! distinguish two very different failures: computing the wrong *quantity*, and
!! computing the right quantity but projecting it onto the orbitals differently.
!! Summing over the orbitals removes the second possibility entirely.
!!
!! That is not a hope, it is an identity, and it is recorded in
!! `mqc_libcint_cphf`: per orbital the two orderings of measure and respond give
!! `h^A P M^-1 h^B` and `h^A M^-1 P h^B`, which differ because the projector `P`
!! does not commute with `M^-1`; summed over every orbital `P` is the identity and
!! they coincide. So the sum is projection independent, ordering independent, and
!! localization independent -- three of the free choices that made the per-orbital
!! comparison ambiguous.
!!
!! **What each outcome would mean.** GAMESS's written values are not the tensor
!! its `LDQPOL` computes; they are `DQSHIFT`'s translation of it, which mixes the
!! dipole-dipole polarizability in. That shift is known exactly and is undone in
!! the Python from GAMESS's own dipole-dipole tensors and centroids, so both sides
!! of the comparison are pre-shift totals. Then:
!!
!!   * if the sums agree, the quantity is right and the difference is entirely in
!!     the per-orbital decomposition -- a much smaller and better-localized
!!     problem than the one on record;
!!   * if the sums disagree, the difference is upstream of any localization, and no
!!     amount of work on the projection would ever have closed it.
!!
!! Either way it is worth knowing, and neither is knowable from the per-orbital
!! numbers.
!!
!! The operators are GAMESS's: the three Cartesian dipoles measure, and the six
!! unique traceless Buckingham quadrupole components drive, in its `XX YY ZZ XY XZ
!! YZ` order, built from the raw second moments the way `prpel.src:5625` does.
!! Everything is expanded about the centre of mass, which is what `DQSHIFT` shifts
!! *from*.
program check_dipquad_sumrule
   use pic_types, only: dp
   use mqc_elements, only: element_mass
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_multipole, only: multipole_matrices
   use mqc_libcint_cphf, only: distributed_dynamic_cross, &
                               casimir_polder_frequencies, N_CASIMIR_POLDER
   use mqc_error, only: error_t
   implicit none

   !> GAMESS's Bohr, since this compares against a GAMESS potential.
   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp

   !> libcint's full-Cartesian quadrupole slots, which run xx,xy,xz,yx,...,zz.
   integer, parameter :: QXX = 1, QXY = 2, QXZ = 3, QYY = 5, QYZ = 6, QZZ = 9

   type(libcint_molecule_t) :: mol
   type(rhf_result_t) :: scf
   type(error_t) :: err
   real(dp), allocatable :: dip(:, :, :), quad(:, :, :), theta(:, :, :)
   real(dp), allocatable :: buckingham(:, :, :)
   real(dp), allocatable :: tensors(:, :, :, :), centroids(:, :)
   real(dp), allocatable :: dq_dr(:, :, :, :), dq_rd(:, :, :, :), dd(:, :, :, :)
   real(dp), allocatable :: qq(:, :, :, :), qq_raw(:, :, :, :)
   real(dp), allocatable :: total(:, :, :), buck_total(:, :, :)
   real(dp) :: c(3, 3), com(3), mass_total
   real(dp) :: nu(N_CASIMIR_POLDER)
   integer :: z(3)
   integer :: i, k, m, f, unit, n_lmo
   character(len=2) :: symbols(3)

   z = [8, 1, 1]
   symbols = ["O ", "H ", "H "]
   c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])

   ! The centre of mass, in the same masses GAMESS prints in COORDINATES.
   com = 0.0_dp
   mass_total = 0.0_dp
   do i = 1, 3
      com = com + element_mass(z(i))*c(:, i)
      mass_total = mass_total + element_mass(z(i))
   end do
   com = com/mass_total

   call build_libcint_molecule(z, symbols, c, "6-31g*", mol, err)
   if (err%has_error()) then
      write (*, "(A,A)") "[dqsum] basis failed: ", err%get_message()
      error stop 1
   end if
   call run_libcint_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
   if (err%has_error() .or. .not. scf%converged) then
      write (*, "(A)") "[dqsum] SCF failed"
      error stop 1
   end if

   call multipole_matrices(mol, com, 1, dip, err)
   if (err%has_error()) then
      write (*, "(A,A)") "[dqsum] dipole integrals failed: ", err%get_message()
      error stop 1
   end if
   call multipole_matrices(mol, com, 2, quad, err)
   if (err%has_error()) then
      write (*, "(A,A)") "[dqsum] quadrupole integrals failed: ", err%get_message()
      error stop 1
   end if

   ! Both candidate driving operators, as the *full* nine Cartesian components so
   ! that no expansion of six unique values into nine slots has to be guessed at.
   ! The traceless Buckingham form is what GAMESS builds internally
   ! (`prpel.src:5625`); the raw second moment is what the recovered tensor's
   ! non-zero trace points to. Dumping both is what makes the comparison a
   ! measurement rather than a choice.
   allocate (theta(mol%nao, mol%nao, 9))
   theta = quad
   allocate (buckingham(mol%nao, mol%nao, 9))
   buckingham(:, :, QXX) = 0.5_dp*(2.0_dp*quad(:, :, QXX) - quad(:, :, QYY) &
                                   - quad(:, :, QZZ))
   buckingham(:, :, QYY) = 0.5_dp*(2.0_dp*quad(:, :, QYY) - quad(:, :, QXX) &
                                   - quad(:, :, QZZ))
   buckingham(:, :, QZZ) = 0.5_dp*(2.0_dp*quad(:, :, QZZ) - quad(:, :, QXX) &
                                   - quad(:, :, QYY))
   buckingham(:, :, QXY) = 1.5_dp*quad(:, :, QXY)
   buckingham(:, :, QXZ) = 1.5_dp*quad(:, :, QXZ)
   buckingham(:, :, QYZ) = 1.5_dp*quad(:, :, QYZ)
   buckingham(:, :, 4) = buckingham(:, :, QXY)
   buckingham(:, :, 7) = buckingham(:, :, QXZ)
   buckingham(:, :, 8) = buckingham(:, :, QYZ)

   nu = casimir_polder_frequencies()
   allocate (total(3, 9, N_CASIMIR_POLDER), buck_total(3, 9, N_CASIMIR_POLDER))
   call summed_response(dip, theta, total)
   call summed_response(dip, buckingham, buck_total)

   ! Per orbital as well, in *both* operator orderings. The sum is insensitive to
   ! which operator drives and which measures; a single orbital is not, because
   ! the projector onto the localized set does not commute with the response
   ! operator. So both are dumped and the comparison decides.
   call distributed_dynamic_cross(mol, scf%orbitals, scf%orbital_energies, &
                                  scf%n_occupied, nu, dip, buckingham, dq_dr, &
                                  centroids, err, n_core=1)
   if (err%has_error()) error stop "dip measures"
   call distributed_dynamic_cross(mol, scf%orbitals, scf%orbital_energies, &
                                  scf%n_occupied, nu, buckingham, dip, dq_rd, &
                                  centroids, err, n_core=1)
   if (err%has_error()) error stop "quad measures"
   ! And the dipole-dipole tensors the shift mixes in -- ours, not GAMESS's, so
   ! that what is tested is a recipe we could actually run.
   call distributed_dynamic_cross(mol, scf%orbitals, scf%orbital_energies, &
                                  scf%n_occupied, nu, dip, dip, dd, centroids, &
                                  err, n_core=1)
   if (err%has_error()) error stop "dipole-dipole"
   ! And the quadrupole-quadrupole response, which `LMOQQPOL` carries and whose
   ! write-time shift `QQSHIFT` takes the two tensors above as inputs.
   call distributed_dynamic_cross(mol, scf%orbitals, scf%orbital_energies, &
                                  scf%n_occupied, nu, buckingham, buckingham, qq, &
                                  centroids, err, n_core=1)
   if (err%has_error()) error stop "quadrupole-quadrupole"

   open (newunit=unit, file="/tmp/mqc_dqsum.txt", status="replace", action="write")
   write (unit, "(I0,1X,I0)") N_CASIMIR_POLDER, n_lmo
   write (unit, "(es25.16e3)") nu
   write (unit, "(es25.16e3)") total
   write (unit, "(es25.16e3)") buck_total
   write (unit, "(es25.16e3)") com
   write (unit, "(es25.16e3)") centroids
   write (unit, "(es25.16e3)") dq_dr
   write (unit, "(es25.16e3)") dq_rd
   write (unit, "(es25.16e3)") dd
   write (unit, "(es25.16e3)") qq
   ! And the same with the *raw* second moment driving and measuring. The trace of
   ! GAMESS's LMOQQPOL block cannot be removed by any combination of QQSHIFT's
   ! terms, so its quadrupole is not the traceless Buckingham form; this is the
   ! alternative.
   call distributed_dynamic_cross(mol, scf%orbitals, scf%orbital_energies, &
                                  scf%n_occupied, nu, theta, theta, qq_raw, &
                                  centroids, err, n_core=1)
   if (err%has_error()) error stop "raw quadrupole-quadrupole"
   write (unit, "(es25.16e3)") qq_raw
   close (unit)

   write (*, "(A,I0,A,I0,A)") "  ", n_lmo, " localized orbitals, ", &
      N_CASIMIR_POLDER, " frequencies"
   write (*, "(A,3F14.8)") "  centre of mass (Bohr) ", com
   write (*, "(A)") "  summed response at the lowest frequency, raw second moment"
   write (*, "(A)") "  driving; rows XYZ (measure) by the nine quadrupole slots:"
   do m = 1, 3
      write (*, "(A,9F12.6)") "        ", (total(m, i, 1), i=1, 9)
   end do
   write (*, "(A)") ""
   write (*, "(A)") "[dqsum] wrote /tmp/mqc_dqsum.txt; "// &
      "run validation/check_dipquad_sumrule.py"

   call mol%destroy()

contains

   subroutine summed_response(measure, respond, summed)
      !! One response calculation, summed over the localized orbitals
      real(dp), intent(in) :: measure(:, :, :), respond(:, :, :)
      real(dp), intent(out) :: summed(:, :, :)

      integer :: kk, ff

      call distributed_dynamic_cross(mol, scf%orbitals, scf%orbital_energies, &
                                     scf%n_occupied, nu, measure, respond, &
                                     tensors, centroids, err, n_core=1)
      if (err%has_error()) then
         write (*, "(A,A)") "[dqsum] response failed: ", err%get_message()
         error stop 1
      end if
      n_lmo = size(tensors, 3)
      summed = 0.0_dp
      do ff = 1, N_CASIMIR_POLDER
         do kk = 1, n_lmo
            summed(:, :, ff) = summed(:, :, ff) + tensors(:, :, kk, ff)
         end do
      end do
   end subroutine summed_response

end program check_dipquad_sumrule
