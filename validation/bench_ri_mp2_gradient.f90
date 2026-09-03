!! Where the time goes in the RI-MP2 gradient
program bench_ri_mp2_gradient
   !! Warm and on a wall clock, for the reasons `bench_mp2_gradient` gives: a
   !! cold call measures lazy symbol resolution and first-touch allocation, and
   !! `cpu_time` sums over threads so a threaded routine appears to slow down as
   !! it speeds up.
   !!
   !! Two costs are mixed together in every number below, and it is worth
   !! knowing which is which. The fitted work -- the amplitudes, the three-index
   !! density and the Lagrangian -- is `n_occ^2 n_vir^2 n_aux`. The reference
   !! work -- the exact integrals for the Z-vector and the differentiated
   !! quartets -- is `n_ao^4`, and is shared with the conventional gradient next
   !! door.
   !!
   !! **The reference work dominates at every size in this file**, measured
   !! rather than assumed: on water/cc-pVTZ the Z-vector solve alone is a third
   !! of the wall clock and the fitted steps together are about a tenth. That is
   !! what these molecules are, not what the method is -- they have twenty
   !! electrons and a large basis each, and growing the basis at a fixed
   !! electron count is precisely the regime where `n_ao^4` outruns
   !! `n_occ^2 n_vir^2 n_aux`. A case with more electrons rather than more
   !! functions per electron would invert it, and none is cheap enough to sit
   !! in a benchmark that has to finish.
   use pic_types, only: dp, int64
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_ri_mp2_gradient, only: czt_ri_mp2_gradient
   use omp_lib, only: omp_get_max_threads
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none

   integer, parameter :: N_DIM = 3

   write (*, "(a,i0,a)") "== RI-MP2 gradient, ", omp_get_max_threads(), " threads"

   call bench_case("H2O / cc-pvdz", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "cc-pvdz", "cc-pvdz-rifit", 10)

   call bench_case("H2O / cc-pvtz", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "cc-pvtz", "cc-pvtz-rifit", 10)

   call bench_case("(H2O)2 / cc-pvdz", [8, 1, 1, 8, 1, 1], &
                   ["O", "H", "H", "O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp, &
                            0.0_dp, 0.0_dp, 5.6_dp, &
                            0.0_dp, -1.4308_dp, 6.7078_dp, &
                            0.0_dp, 1.4308_dp, 6.7078_dp], [N_DIM, 6]), &
                   "cc-pvdz", "cc-pvdz-rifit", 20)

   ! The largest case here, and it does *not* make the fitted work dominant --
   ! measured, not assumed. Growing the basis at a fixed electron count grows
   ! `n_ao^4` faster than `n_occ^2 n_vir^2 n_aux`, so the reference terms take a
   ! larger share here than at cc-pVDZ, not a smaller one. What would move the
   ! balance is more electrons rather than more functions per electron, and no
   ! case in this file has them.
   call bench_case("(H2O)2 / cc-pvtz", [8, 1, 1, 8, 1, 1], &
                   ["O", "H", "H", "O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp, &
                            0.0_dp, 0.0_dp, 5.6_dp, &
                            0.0_dp, -1.4308_dp, 6.7078_dp, &
                            0.0_dp, 1.4308_dp, 6.7078_dp], [N_DIM, 6]), &
                   "cc-pvtz", "cc-pvtz-rifit", 20)

contains

   subroutine bench_case(label, numbers, symbols, coords, basis, aux_basis, nelec)
      character(len=*), intent(in) :: label
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: nelec

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: error
      real(dp), allocatable :: gradient(:, :)
      integer(int64) :: t0, t1, rate
      real(dp) :: seconds

      call build_czt_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         return
      end if
      call build_czt_molecule(numbers, symbols, coords, aux_basis, aux, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         return
      end if
      call run_czt_rhf(mol, nelec, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         return
      end if

      ! Warm-up, discarded.
      call czt_ri_mp2_gradient(mol, aux, scf%orbitals, scf%orbital_energies, &
                               nelec/2, gradient, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         return
      end if
      deallocate (gradient)

      call system_clock(t0, rate)
      call czt_ri_mp2_gradient(mol, aux, scf%orbitals, scf%orbital_energies, &
                               nelec/2, gradient, error)
      call system_clock(t1)
      seconds = real(t1 - t0, dp)/real(rate, dp)

      write (*, "(a,a,a,i0,a,i0,a,f9.3,a,es12.4)") "  ", label, "  nao=", mol%nao, &
         "  naux=", aux%nao, "   ", seconds, " s   |g|=", sqrt(sum(gradient**2))
      flush (output_unit)
      deallocate (gradient)
      call mol%destroy()
      call aux%destroy()
   end subroutine bench_case

end program bench_ri_mp2_gradient
