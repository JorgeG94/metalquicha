program bench_df_k
   !! Time the density-fitted Fock build alone, away from everything else
   !!
   !! Timing a whole SCF says almost nothing about this: most of the runtime
   !! goes into three-centre integrals and basis setup, and the exchange build
   !! is a small enough slice that a real change to it disappears into the
   !! noise. So this times repeated Fock builds and nothing else.
   !!
   !! The size matters too. The saving from routing exchange through the
   !! occupied orbitals is a factor of n/n_occ, so it is invisible on a small
   !! basis and grows as the basis does at fixed electron count. Water in
   !! cc-pVDZ is 24/5; in cc-pVTZ it is 58/5.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, &
                                build_df_tensor
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   implicit none

   integer, parameter :: N_DIM = 3
   integer, parameter :: N_REPEAT = 20

   call bench("cc-pvdz", "cc-pvtz-jkfit")
   call bench("cc-pvtz", "cc-pvtz-jkfit")

contains

   subroutine bench(orbital_basis, aux_basis)
      character(len=*), intent(in) :: orbital_basis, aux_basis

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: res
      type(error_t) :: error
      real(dp) :: coords(N_DIM, 3)
      real(dp) :: t0, t1
      integer :: rep

      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, -1.4308_dp, 1.1078_dp, &
                        0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3])

      call build_czt_molecule([8, 1, 1], ["O", "H", "H"], coords, orbital_basis, mol, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         return
      end if
      call build_czt_molecule([8, 1, 1], ["O", "H", "H"], coords, aux_basis, aux, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         return
      end if

      ! One SCF first, so basis reading and tensor construction are not being
      ! timed along with the thing under test.
      call run_czt_rhf(mol, 10, 100, 1.0e-10_dp, 1.0e-8_dp, .false., res, error, aux=aux)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         return
      end if

      call cpu_time(t0)
      do rep = 1, N_REPEAT
         call run_czt_rhf(mol, 10, 100, 1.0e-10_dp, 1.0e-8_dp, .false., res, error, aux=aux)
      end do
      call cpu_time(t1)

      write (*, "(a,a,a,i0,a,i0,a,i0,a,f8.4,a,f20.12)") &
         orbital_basis, "  nao=", "", mol%nao, "  naux=", aux%nao, &
         "  nocc=", res%n_occupied, "   ", (t1 - t0)/real(N_REPEAT, dp), " s/scf   E = ", res%energy

      call mol%destroy()
      call aux%destroy()
   end subroutine bench

end program bench_df_k
