!! Baseline: what a Fock build costs, so optimisation has a target
!!
!!     ./build/bench_fock
!!
!! The MAKEFP response solves are thousands of Fock builds, so how fast the whole
!! potential can possibly go is set by two numbers: how quickly libcint turns shell
!! quartets into integrals, and how quickly a stored tensor can be contracted. This
!! measures both, per basis, with nothing else in the way.
!!
!! **What each number is for.**
!!
!!   * `direct, 1 density` -- integrals recomputed and contracted once. Divided into
!!     quartets per second, this is libcint's throughput on this machine and the hard
!!     floor for any integral-direct scheme.
!!   * `direct, N densities` -- the same integrals contracted into several Fock
!!     matrices in one pass. The per-density cost falling with `N` is the amortisation
!!     the batched response solver relies on; where it stops falling is where the
!!     contraction rather than the integrals dominates.
!!   * `stored` -- one contraction of the four-index tensor, and separately the cost
!!     of computing it. Their ratio says how many builds a run must do before storing
!!     is the cheaper choice.
!!
!! Reported as ops per second as well as seconds, because seconds alone say nothing
!! about whether a build is near what the hardware allows. The stored contraction is
!! `2 n^4` multiply-adds, so `4 n^4` flops; against a machine peak in the hundreds of
!! GFlop/s a naive triple-nested loop lands one to two orders below, and that gap is
!! the case for BLAS.
program bench_fock
   use pic_types, only: dp, int64
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf, build_fock
   use mqc_libcint_direct, only: schwarz_bounds, build_fock_direct, &
                                 build_fock_direct_many, direct_stats_t
   use omp_lib, only: omp_get_max_threads
   implicit none

   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp
   character(len=16), parameter :: BASES(4) = [character(len=16) :: &
                                               "6-31g*", "cc-pvdz", "cc-pvtz", "cc-pvqz"]
   !> Batch sizes to show the amortisation. 108 is what the quadrupole-quadrupole
   !> block of a potential actually asks for: nine operators at twelve frequencies.
   integer, parameter :: BATCHES(4) = [1, 12, 36, 108]

   integer :: ib
   real(dp) :: c(3, 3)
   integer :: z(3)
   character(len=2) :: symbols(3)

   z = [8, 1, 1]
   symbols = ["O ", "H ", "H "]
   c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])

   write (*, "(A,I0,A)") "  water, ", omp_get_max_threads(), " threads"
   write (*, "(A)") ""
   do ib = 1, size(BASES)
      call one_basis(trim(BASES(ib)))
   end do
   write (*, "(A)") ""
   write (*, "(A)") "[bench] baseline recorded"

contains

   subroutine one_basis(basis)
      character(len=*), intent(in) :: basis

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: bounds(:, :), h(:, :), fock(:, :)
      real(dp), allocatable :: dens(:, :, :), focks(:, :, :), eri(:, :, :, :)
      real(dp) :: t0, t1, per_build, quartets_per_s, ops
      integer :: n, ib, nb, rep, reps
      logical :: stored_fits

      call build_libcint_molecule(z, symbols, c, basis, mol, err)
      if (err%has_error()) then
         write (*, "(A,A,A)") "  ", basis, ": basis failed"
         return
      end if
      n = mol%nao

      ! A converged density, so screening sees realistic magnitudes. A random one
      ! would screen differently and flatter the throughput.
      call run_libcint_rhf(mol, 10, 200, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         write (*, "(A,A,A)") "  ", basis, ": SCF failed"
         call mol%destroy()
         return
      end if
      call schwarz_bounds(mol, bounds, err)
      if (err%has_error()) return

      allocate (h(n, n), fock(n, n))
      h = 0.0_dp

      write (*, "(A,A,A,I0,A)") "  ", basis, ", ", n, " orbitals"

      ! --- one density, integral direct -------------------------------------------
      reps = max(1, min(20, 2000000/max(1, n*n)))
      call cpu_seconds(t0)
      do rep = 1, reps
         call build_fock_direct(mol, h, scf%density, bounds, fock, stats, err)
      end do
      call cpu_seconds(t1)
      per_build = (t1 - t0)/real(reps, dp)
      quartets_per_s = real(stats%quartets_computed, dp)/max(per_build, 1.0e-12_dp)
      write (*, "(A,ES11.3,A,ES11.3,A,F6.1,A)") &
         "      direct, 1 density   ", per_build, " s   ", quartets_per_s, &
         " quartets/s   ", 100.0_dp*stats%screened_fraction(), "% screened"

      ! --- many densities, one pass ----------------------------------------------
      do ib = 1, size(BATCHES)
         nb = BATCHES(ib)
         allocate (dens(n, n, nb))
         do rep = 1, nb
            dens(:, :, rep) = scf%density
         end do
         reps = max(1, min(5, 400000/max(1, n*n*nb)))
         call cpu_seconds(t0)
         do rep = 1, reps
            if (allocated(focks)) deallocate (focks)
            call build_fock_direct_many(mol, h, dens, bounds, focks, stats, err)
         end do
         call cpu_seconds(t1)
         per_build = (t1 - t0)/real(reps, dp)
         write (*, "(A,I4,A,ES11.3,A,ES11.3,A)") &
            "      direct, ", nb, " densities ", per_build, " s   ", &
            per_build/real(nb, dp), " s each"
         deallocate (dens)
         if (allocated(focks)) deallocate (focks)
      end do

      ! --- stored tensor ----------------------------------------------------------
      stored_fits = real(n, dp)**4*8.0_dp <= 4.0e9_dp
      if (stored_fits) then
         call cpu_seconds(t0)
         call mol%eris(eri)
         call cpu_seconds(t1)
         write (*, "(A,ES11.3,A,F8.2,A)") "      stored, computing it ", t1 - t0, &
            " s   ", real(n, dp)**4*8.0_dp/1.0e9_dp, " GB"

         reps = max(1, min(20, 20000000/max(1, n*n*n)))
         call cpu_seconds(t0)
         do rep = 1, reps
            call build_fock(h, eri, scf%density, fock)
         end do
         call cpu_seconds(t1)
         per_build = (t1 - t0)/real(reps, dp)
         ops = 4.0_dp*real(n, dp)**4
         write (*, "(A,ES11.3,A,F8.2,A)") "      stored, contracting  ", per_build, &
            " s   ", ops/max(per_build, 1.0e-12_dp)/1.0e9_dp, " GFlop/s"
         deallocate (eri)
      else
         write (*, "(A)") "      stored: tensor too large to time here"
      end if

      call mol%destroy()
      deallocate (bounds, h, fock)
   end subroutine one_basis

   subroutine cpu_seconds(t)
      !! Wall clock, which is what a user waits: `cpu_time` sums over threads.
      real(dp), intent(out) :: t
      integer(int64) :: count, rate
      call system_clock(count, rate)
      t = real(count, dp)/real(rate, dp)
   end subroutine cpu_seconds

end program bench_fock
