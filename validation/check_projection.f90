!! Manual check of the projection data a potential carries
!!
!!     cmake -B build -DMQC_ENABLE_LIBCINT=ON
!!     ./build/check_projection
!!     python3 validation/check_projection.py
!!
!! The exchange-repulsion half of a potential is not computed, it is *dumped*:
!! GAMESS emits the localized orbitals, the Fock matrix in their basis, and the
!! basis set itself, and the energy is formed at use time from inter-fragment
!! overlaps. So this milestone is serialization, and its risks are serialization
!! risks.
!!
!! **The Fock matrix is the part that needs no basis ordering, so it comes first.**
!! In the LMO basis it is `W^T diag(eps) W`, with `W = C_occ^T S C_loc` the map
!! from canonical occupied orbitals to localized ones -- no AO Fock matrix
!! required, and nothing that depends on how libcint or GAMESS happens to order
!! basis functions. Only the *ordering of the LMOs* enters, and that is already
!! known to agree: the polarizability comparison pairs our LMOs to GAMESS's CT1..4
!! in order, by centroid.
!!
!! The orbital coefficients and the basis are a different matter -- those are
!! indexed by AO, and libcint's ordering is not GAMESS's -- so they are left to a
!! later step rather than emitted wrongly here.
program check_projection
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_localize, only: boys_localize
   use mqc_libcint_integrals, only: shell_dim
   use libcint_fortran, only: LIBCINT_ATOM_OF, LIBCINT_ANG_OF
   use mqc_error, only: error_t
   implicit none

   !> GAMESS's Bohr, as in check_dma -- this compares against a GAMESS potential.
   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp
   integer :: failures

   failures = 0
   call one_case(1, "6-31g*", 1)
   ! A second case carrying f functions, read in Cartesian form the way GAMESS reads
   ! a Pople set, so the ordering map for l = 3 can be derived against its own output.
   ! In a general orientation: planar water puts an exact zero in every function with
   ! an odd power of y, and a slot that is zero on both sides fits any scale factor,
   ! so the components it would determine are the ones left undetermined.
   call one_case(2, "6-31g(2df,p)", 1, force_cartesian=.true., rotated=.true.)

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[projection] wrote every case; run validation/check_projection.py"
   else
      write (*, "(A,I0,A)") "[projection] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine one_case(case_index, basis, n_core, force_cartesian, rotated)
      integer, intent(in) :: case_index
      character(len=*), intent(in) :: basis
      integer, intent(in) :: n_core
      logical, intent(in), optional :: force_cartesian
      !> Use the same water turned into a frame with no zero coordinate.
      logical, intent(in), optional :: rotated

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: loc(:, :), cen(:, :), s(:, :), sc(:, :), w(:, :)
      real(dp), allocatable :: fock_lmo(:, :), scaled(:, :)
      real(dp) :: c(3, 3)
      real(dp) :: core_leak
      integer :: z(3)
      integer :: unit, i, j, k, n_lmo, n_occ
      character(len=2) :: symbols(3)
      character(len=64) :: name

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      if (present(rotated)) then
         if (rotated) &
            c = reshape([0.0000000000_dp, 0.0000000000_dp, 0.0000000000_dp, &
                         0.2493568104_dp, 0.5884773575_dp, 0.7141751472_dp, &
                         -0.0661597287_dp, 0.5685269177_dp, -0.7674361694_dp], [3, 3])*ANG
      end if

      call build_libcint_molecule(z, symbols, c, basis, mol, err, &
                                  force_cartesian=force_cartesian)
      if (err%has_error()) then
         write (*, "(A,A)") "[projection] basis failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call run_libcint_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         write (*, "(A)") "[projection] SCF failed"
         failures = failures + 1
         call mol%destroy()
         return
      end if
      n_occ = scf%n_occupied
      n_lmo = n_occ - n_core

      call boys_localize(mol, scf%orbitals(:, n_core + 1:n_occ), n_lmo, loc, cen, err)
      if (err%has_error()) then
         write (*, "(A,A)") "[projection] localization failed: ", err%get_message()
         failures = failures + 1
         call mol%destroy()
         return
      end if

      ! W = C_occ^T S C_loc, over *all* occupied orbitals. The Fock operator is
      ! diagonal in the canonical basis, so the LMO Fock is W^T diag(eps) W and
      ! needs no AO Fock matrix at all.
      call mol%overlap(s)
      allocate (sc(mol%nao, n_lmo), w(n_occ, n_lmo))
      call pic_gemm(s, loc, sc)
      call pic_gemm(scf%orbitals(:, 1:n_occ), sc, w, transa="T")

      ! How much of the excluded core the localized set still touches. It should be
      ! zero -- the valence LMOs are built from valence orbitals alone -- and if it
      ! were not, the Fock matrix below would be missing a contribution.
      core_leak = 0.0_dp
      do k = 1, n_core
         core_leak = max(core_leak, maxval(abs(w(k, :))))
      end do

      allocate (scaled(n_occ, n_lmo), fock_lmo(n_lmo, n_lmo))
      do j = 1, n_lmo
         do k = 1, n_occ
            scaled(k, j) = scf%orbital_energies(k)*w(k, j)
         end do
      end do
      call pic_gemm(w, scaled, fock_lmo, transa="T")

      write (name, "(A,I0,A)") "/tmp/mqc_projection_", case_index, ".txt"
      open (newunit=unit, file=trim(name), status="replace", action="write")
      write (unit, "(A)") trim(basis)
      write (unit, "(I0,1X,I0,1X,I0)") mol%nao, n_occ, n_lmo
      write (unit, "(es25.16e3)") fock_lmo
      write (unit, "(3es25.16e3)") cen
      ! The localized coefficients, in libcint's AO order. The Python works out
      ! the permutation onto GAMESS's, which is what the projection sections need.
      write (unit, "(es25.16e3)") loc
      ! The shell layout, so the consumer can apply the Cartesian-d permutation
      ! wherever a d shell happens to sit rather than at a hard-coded offset.
      ! The canonical occupied orbital energies, which is all CTFOK turns out to be.
      write (unit, "(es25.16e3)") scf%orbital_energies(1:n_occ)
      write (unit, "(I0)") mol%nbas
      do k = 1, mol%nbas
         write (unit, "(4(I0,1X))") mol%bas(LIBCINT_ATOM_OF, k) + 1, &
            mol%bas(LIBCINT_ANG_OF, k), mol%shell_offset(k), &
            shell_dim(mol%cartesian, k - 1, mol%bas)
      end do
      ! All the canonical MOs, in libcint's AO order, appended last so the readers
      ! that parse this file positionally are unaffected. CTVEC is 5 occupied plus
      ! 2 valence virtuals, so the canonical set is what to test it against.
      write (unit, "(es25.16e3)") scf%orbitals
      close (unit)

      write (*, "(A,A8,A,I0,A,I0,A,ES9.2)") "  ", trim(basis), "  LMOs ", n_lmo, &
         "  of ", n_occ, " occupied   core leak ", core_leak
      write (*, "(A)") "        LMO Fock diagonal (Hartree):"
      write (*, "(A,4F14.8)") "        ", (fock_lmo(i, i), i=1, n_lmo)

      if (core_leak > 1.0e-10_dp) then
         write (*, "(A)") "    FAIL the localized orbitals overlap the excluded core"
         failures = failures + 1
      end if
      ! Symmetric by construction; asserted because a transposed W would not be.
      if (maxval(abs(fock_lmo - transpose(fock_lmo))) > 1.0e-12_dp) then
         write (*, "(A)") "    FAIL the LMO Fock matrix is not symmetric"
         failures = failures + 1
      end if

      call mol%destroy()
      deallocate (loc, cen, s, sc, w, scaled, fock_lmo)
   end subroutine one_case

end program check_projection
