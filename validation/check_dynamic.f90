!! Manual check of the dynamic polarizabilities
!!
!!     cmake -B build -DMQC_ENABLE_LIBCINT=ON
!!     ./build/check_dynamic
!!     python3 validation/check_dynamic.py
!!
!! Milestone 6: the twelve imaginary frequencies a potential is tabulated at, and
!! what fragment dispersion is built from -- the Casimir-Polder integral over these
!! is the C6 between two fragments.
!!
!! **Two checks here need no reference.** The `nu = 0` member must reproduce the
!! static polarizability exactly, since the frequency-dependent equations reduce to
!! the static ones there, and the distributed set at `nu = 0` must likewise
!! reproduce the static distributed tensors orbital by orbital. Both are exact
!! identities rather than approximations, and between them they check that the
!! dynamic solver and the two projections agree with the static machinery already
!! validated against GAMESS.
!!
!! The frequencies themselves are worth checking too: `casimir_polder_frequencies`
!! claims to reproduce the twelve values GAMESS stamps its dynamic blocks with, and
!! the Python confirms it against the reference rather than trusting the formula.
program check_dynamic
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_cphf, only: static_polarizability, dynamic_polarizability, &
                               N_CASIMIR_POLDER, &
                               casimir_polder_frequencies, &
                               distributed_dynamic_polarizability
   use mqc_error, only: error_t
   implicit none
   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp
   type(libcint_molecule_t) :: mol
   type(rhf_result_t) :: scf
   type(error_t) :: err
   real(dp), allocatable :: a(:, :, :)
   real(dp) :: c(3, 3), stat(3, 3)
   real(dp) :: nu(N_CASIMIR_POLDER), probe(N_CASIMIR_POLDER + 1)
   integer :: i
   c = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
   call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g*", mol, err)
   call run_libcint_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
   call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                              scf%n_occupied, stat, err)
   nu = casimir_polder_frequencies()
   probe(1) = 0.0_dp
   probe(2:N_CASIMIR_POLDER + 1) = nu
   call dynamic_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                               scf%n_occupied, probe, a, err)
   if (err%has_error()) then
      write (*, *) "FAILED: ", err%get_message()
      stop 1
   end if
   write (*, "(a,f14.8)") " static alpha_iso from static_polarizability: ", &
      (stat(1, 1) + stat(2, 2) + stat(3, 3))/3.0_dp
   write (*, "(a,f14.8)") " static alpha_iso from the nu=0 dynamic solve:  ", &
      (a(1, 1, 1) + a(2, 2, 1) + a(3, 3, 1))/3.0_dp
   write (*, "(a,es10.2)") " worst elementwise difference: ", &
      maxval(abs(a(:, :, 1) - stat))
   write (*, "(a)") ""
   write (*, "(a)") " nu           alpha_iso(i nu)   asymmetry"
   do i = 1, N_CASIMIR_POLDER
      write (*, "(f11.6,f16.8,es12.2)") nu(i), &
         (a(1, 1, i + 1) + a(2, 2, i + 1) + a(3, 3, i + 1))/3.0_dp, &
         maxval(abs(a(:, :, i + 1) - transpose(a(:, :, i + 1))))
   end do
   block
      real(dp), allocatable :: t(:, :, :, :), cen(:, :)
      integer :: j, unit
      call distributed_dynamic_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                              scf%n_occupied, probe, t, cen, err, &
                                              n_core=1)
      if (err%has_error()) then
         write (*, *) "distributed failed: ", err%get_message()
         stop 1
      end if
      write (*, "(a)") ""
      write (*, "(a,i0,a,i0)") " distributed: LMOs ", size(t, 3), "  frequencies ", size(t, 4)
      open (newunit=unit, file="/tmp/mqc_dyn.txt", status="replace", action="write")
      write (unit, "(3(i0,1x))") size(t, 3), size(t, 4), 0
      write (unit, "(es25.16e3)") probe
      write (unit, "(3es25.16e3)") cen
      write (unit, "(es25.16e3)") t
      close (unit)
      write (*, "(a)") " nu          sum_iso over LMOs"
      do i = 1, size(t, 4)
         write (*, "(f11.6,f16.8)") probe(i), &
            (sum(t(1, 1, :, i)) + sum(t(2, 2, :, i)) + sum(t(3, 3, :, i)))/3.0_dp
      end do
      deallocate (t, cen)
   end block
   call mol%destroy()
end program check_dynamic
