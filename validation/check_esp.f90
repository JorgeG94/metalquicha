!! Manual check of the electrostatic potential integrals
!!
!!     cmake -B build -DMQC_ENABLE_LIBCINT=ON
!!     ./build/check_esp
!!     python3 validation/check_esp.py
!!
!! `<chi_u| 1/|r - R_g| |chi_v>` at a list of points. Contracted with a density this
!! is the electronic potential the molecule makes, and it is what a charge-
!! penetration screening fit is fitted to: GAMESS optimizes its damping exponents by
!! matching a damped classical multipole potential to the quantum one on a grid, so
!! this integral is the target. The CPU path had no way to compute it -- only the
!! cuEST backend did, through its PCM plan -- so this is the enabling piece for that
!! milestone, and for a CPU polarizable-continuum model later.
!!
!! **Checked elementwise, not on a contracted potential.** libcint runs the grid
!! index fastest and the shell pair outside it, the opposite nesting from the
!! multipole blocks, and a transposed unpack gives a matrix that is plausible at one
!! point and wrong at the others. A single contracted number could hide that; every
!! element at every point cannot.
program check_esp
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_esp, only: esp_matrices
   use mqc_error, only: error_t
   implicit none
   real(dp), parameter :: ANG = 1.8897261254578281_dp

   !> Probe points. Deliberately off-centre, asymmetric, and none on an atom or a
   !> symmetry plane -- a point on either would make several integrals vanish and
   !> hide a transposed unpack.
   integer, parameter :: N_PROBE = 5
   type(libcint_molecule_t) :: mol
   type(rhf_result_t) :: scf
   type(error_t) :: err
   real(dp), allocatable :: v(:, :, :)
   real(dp) :: c(3, 3), g(3, N_PROBE)
   integer :: unit, k
   c = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
   ! Points off-centre and asymmetric, none on an atom or a symmetry plane.
   g = reshape([1.3_dp, -0.7_dp, 2.1_dp, -2.0_dp, 1.1_dp, 0.4_dp, &
                0.5_dp, 3.0_dp, -1.2_dp, 4.0_dp, 0.2_dp, 0.9_dp, &
                -1.1_dp, -2.2_dp, -3.3_dp], [3, N_PROBE])
   call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, "cc-pvdz", mol, err)
   call run_libcint_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
   call esp_matrices(mol, g, v, err)
   if (err%has_error()) then
      write (*, *) "FAILED: ", err%get_message()
      stop 1
   end if
   write (*, "(a,3(i0,1x))") " shape: ", size(v, 1), size(v, 2), size(v, 3)
   write (*, "(a)") " electronic potential at each point (a.u.):"
   do k = 1, size(g, 2)
      write (*, "(a,3f8.3,a,f16.10)") "   (", g(:, k), " )  ", -sum(scf%density*v(:, :, k))
   end do
   open (newunit=unit, file="/tmp/mqc_esp.txt", status="replace", action="write")
   write (unit, "(2(i0,1x))") mol%nao, size(g, 2)
   write (unit, "(3es25.16e3)") g
   write (unit, "(es25.16e3)") v
   close (unit)
   call mol%destroy()
end program check_esp
