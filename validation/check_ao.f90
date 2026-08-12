!! Manual check that basis functions evaluate to what libcint's integrals assume
!!
!!     cmake -B build -DMQC_ENABLE_LIBCINT=ON && ./build/check_ao
!!
!! Writes chi_mu(r) at a fixed set of points to a text file, for
!! `validation/check_ao.py` to compare against `pyscf.dft.numint.eval_ao` on the
!! same points and the same basis data. Two codes, one quantity, elementwise.
!!
!! The comparison is worth this much machinery because there is no cheaper way to
!! be sure. libcint does not evaluate basis functions, so the values here come
!! from our own code, and a wrong Cartesian ordering, a doubled normalisation or a
!! transposed spherical transform all produce a basis that is internally
!! consistent and different from the one the integrals were built in. Every one of
!! those converges.
!!
!! The second check -- numerical against analytic overlap -- lives in
!! `test/test_mqc_libcint_ao.f90`, because it needs no external reference.
program check_ao
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_ao, only: eval_ao_block
   use mqc_error, only: error_t
   implicit none

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   integer, parameter :: N_POINTS = 7
   integer :: failures

   failures = 0
   call one_basis("sto-3g")
   call one_basis("cc-pvdz")
   call one_basis("6-31g_st__st_")
   call one_basis("cc-pvtz")

   if (failures == 0) then
      write (*, "(A)") "[ao] wrote every basis; run validation/check_ao.py to compare"
   else
      write (*, "(A,I0,A)") "[ao] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine one_basis(basis)
      !! Dump chi at the probe points for one basis set
      character(len=*), intent(in) :: basis

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: ao(:, :)
      real(dp) :: c(3, 3), pts(3, N_POINTS)
      integer :: unit, ig, mu

      ! The geometry every CPU water case uses, so a failure here and a failure
      ! in the energy suite are about the same molecule.
      c = reshape([0.0_dp, 0.00000000009155_dp*ANG, 0.10077199490609_dp*ANG, &
                   0.0_dp, 0.77250895271063_dp*ANG, -0.46780199741728_dp*ANG, &
                   0.0_dp, -0.77250895280218_dp*ANG, -0.46780199748881_dp*ANG], [3, 3])

      ! Points chosen to exercise the awkward cases rather than the easy ones:
      ! one exactly on a nucleus, where every x**0 must still be one; one far out
      ! where the radial part underflows; and several off-axis so that no
      ! Cartesian component can be zero by symmetry and hide a transposition.
      pts = reshape([ &
                    0.0_dp, 0.00000000009155_dp*ANG, 0.10077199490609_dp*ANG, &
                    0.31_dp, 0.17_dp, -0.23_dp, &
                    -0.44_dp, 0.62_dp, 0.11_dp, &
                    1.30_dp, -0.85_dp, 0.47_dp, &
                    0.05_dp, -0.02_dp, 0.9_dp, &
                    -2.1_dp, 1.7_dp, -1.3_dp, &
                    8.0_dp, -6.0_dp, 5.0_dp], [3, N_POINTS])

      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
      if (err%has_error()) then
         write (*, "(A,A,A,A)") "[ao] ", basis, " basis failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call eval_ao_block(mol, pts, ao, err)
      if (err%has_error()) then
         write (*, "(A,A,A,A)") "[ao] ", basis, " eval failed: ", err%get_message()
         failures = failures + 1
         call mol%destroy()
         return
      end if

      open (newunit=unit, file="/tmp/mqc_ao_"//basis//".txt", status="replace", action="write")
      write (unit, "(I0,1x,I0,1x,L1)") N_POINTS, mol%nao, mol%cartesian
      do ig = 1, N_POINTS
         write (unit, "(3(es25.16e3,1x))") pts(:, ig)
      end do
      do ig = 1, N_POINTS
         do mu = 1, mol%nao
            write (unit, "(es25.16e3)") ao(ig, mu)
         end do
      end do
      close (unit)

      write (*, "(A,A,A,I0,A,L1)") "  wrote ", basis, ": nao=", mol%nao, " cartesian=", &
         mol%cartesian
      call mol%destroy()
   end subroutine one_basis

end program check_ao
