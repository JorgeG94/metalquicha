!! Where the time goes in the MP2 gradient
program bench_mp2_gradient
   !! **Warm, and on a wall clock.** The first call into this path resolves
   !! symbols lazily and allocates for the first time, which on a small molecule
   !! is most of what a cold timing measures -- an earlier benchmark in this
   !! repository reported an 11x gap that was entirely that. And `cpu_time` sums
   !! over threads, so a threaded routine appears to get *slower* as it speeds
   !! up; `system_clock` is what answers the question being asked.
   use pic_types, only: dp, int64
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mp2_gradient, only: czt_mp2_gradient
   use omp_lib, only: omp_get_max_threads
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none

   integer, parameter :: N_DIM = 3
   real(dp), parameter :: BLOCK_BYTES = 5.0e8_dp
      !! What the blocked path would use of its own accord. Passed explicitly
      !! because these molecules are far too small to take that path otherwise.
   character(len=32) :: mode
      !! "dense" measures only that path, so a run cannot be perturbed by the
      !! allocation and release of the other one's arrays.

   mode = ""
   call get_command_argument(1, mode)
   write (*, "(a,i0,a)") "== MP2 gradient, ", omp_get_max_threads(), " threads"

   call bench_case("H2O / cc-pvdz", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "cc-pvdz", 10)

   call bench_case("H2O / cc-pvtz", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "cc-pvtz", 10)

   call bench_case("(H2O)2 / cc-pvdz", [8, 1, 1, 8, 1, 1], &
                   ["O", "H", "H", "O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp, &
                            0.0_dp, 0.0_dp, 5.6_dp, &
                            0.0_dp, -1.4308_dp, 6.7078_dp, &
                            0.0_dp, 1.4308_dp, 6.7078_dp], [N_DIM, 6]), &
                   "cc-pvdz", 20)

contains

   subroutine bench_case(label, numbers, symbols, coords, basis, nelec)
      character(len=*), intent(in) :: label
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: error
      real(dp), allocatable :: gradient(:, :)
      integer(int64) :: t0, t1, rate
      real(dp) :: seconds, blocked_seconds

      call build_czt_molecule(numbers, symbols, coords, basis, mol, error)
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
      call czt_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, nelec/2, &
                            gradient, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         return
      end if
      deallocate (gradient)

      call system_clock(t0, rate)
      call czt_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, nelec/2, &
                            gradient, error)
      call system_clock(t1)
      seconds = real(t1 - t0, dp)/real(rate, dp)

      blocked_seconds = 0.0_dp
      if (trim(mode) == "dense") then
         write (*, "(a,a,a,i0,a,f9.3,a,es12.4)") "  ", label, "  nao=", mol%nao, &
            "   dense ", seconds, " s   |g|=", sqrt(sum(gradient**2))
         flush (output_unit)
         call mol%destroy()
         return
      end if
      deallocate (gradient)

      ! And the same case through the blocked path, which these are all far too
      ! small to take on their own. What it measures is the overhead of never
      ! storing anything: the integrals are rebuilt per block and only the ket
      ! pair's symmetry survives.
      call czt_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, nelec/2, &
                            gradient, error, force_blocked=.true.)
      deallocate (gradient)
      call system_clock(t0, rate)
      call czt_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, nelec/2, &
                            gradient, error, force_blocked=.true.)
      call system_clock(t1)
      blocked_seconds = real(t1 - t0, dp)/real(rate, dp)

      write (*, "(a,a,a,i0,a,f9.3,a,f9.3,a,es12.4)") "  ", label, "  nao=", mol%nao, &
         "   dense ", seconds, " s   blocked ", blocked_seconds, &
         " s   |g|=", sqrt(sum(gradient**2))
      flush (output_unit)
      deallocate (gradient)
      call mol%destroy()
   end subroutine bench_case

end program bench_mp2_gradient
