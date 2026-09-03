program check_fmo_mpi
   !! Does FMO give the same answer on many ranks as on one?
   !!
   !! The distributed path splits both phases across ranks: a monomer pass is
   !! independent within itself, so who computes which fragment is free, and the
   !! pairs are independent of everything once the monomers have settled. Neither
   !! split may move the answer, and that is the only thing worth asserting --
   !! a parallel result that differs from the serial one is wrong however fast
   !! it arrived.
   !!
   !! Run under `mpirun -np N`; N = 1 exercises the serial branch through the
   !! same code, which is worth doing because the two are meant to be one path.
   !!
   !! The reference is computed on this rank with no communicator at all, so the
   !! comparison is against the serial code and not against another rank that
   !! could be wrong in the same way.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger, warning_level
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init, pic_mpi_finalize
   use mqc_error, only: error_t
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   implicit none

   real(dp), parameter :: A2B = 1.8897261254578281_dp
   !! Tight: the two paths do identical arithmetic in a different order, so the
   !! only difference allowed is the reduction's own rounding.
   real(dp), parameter :: TOL = 1.0e-10_dp
   integer, parameter :: N_WATERS = 4

   type(comm_t) :: world
   type(fmo_options_t) :: opts
   type(fmo_result_t) :: serial_res, mpi_res
   type(error_t) :: error
   integer :: z(3*N_WATERS), owner(3*N_WATERS)
   character(len=2) :: sym(3*N_WATERS)
   real(dp) :: xyz(3, 3*N_WATERS), mono(3, 3)
   integer :: w, k, at, n_bad
   character(len=160) :: line

   call pic_mpi_init()
   world = comm_world()
   call logger%configure(warning_level)
   n_bad = 0

   mono = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, -0.7572_dp, 0.5865_dp, &
                   0.0_dp, 0.7572_dp, 0.5865_dp], [3, 3])
   at = 0
   do w = 1, N_WATERS
      do k = 1, 3
         at = at + 1
         z(at) = merge(8, 1, k == 1)
         sym(at) = merge("O ", "H ", k == 1)
         xyz(:, at) = mono(:, k)
         xyz(3, at) = xyz(3, at) + real(w - 1, dp)*2.9_dp
         owner(at) = w
      end do
   end do
   xyz = xyz*A2B

   opts%basis = "sto-3g"

   ! Each method separately: they distribute the same way but assemble
   ! differently, and a reduction dropped from one of them would not show up in
   ! the other.
   call one_method("exact", "fmo", n_bad)
   call one_method("ptc", "mbe", n_bad)
   call one_method("none", "mbe", n_bad)

   if (world%leader()) then
      if (n_bad == 0) then
         write (line, "(a,i0,a)") "FMO on ", world%size(), &
            " rank(s) agrees with the serial answer"
         call logger%configure(warning_level)
         print "(a)", trim(line)
      else
         print "(a)", "FMO MPI checks FAILED"
      end if
   end if

   call pic_mpi_finalize()
   if (n_bad /= 0) stop 1

contains

   subroutine one_method(esp, expansion, n_bad)
      character(len=*), intent(in) :: esp, expansion
      integer, intent(inout) :: n_bad

      character(len=160) :: text

      opts%esp = esp
      opts%expansion = expansion

      ! No communicator: the serial code path, on this rank.
      call run_fmo2(z, sym, xyz, owner, opts, serial_res, error)
      if (error%has_error()) then
         if (world%leader()) print "(a,a)", "  FAIL serial run: ", trim(error%get_message())
         call error%clear()
         n_bad = n_bad + 1
         return
      end if

      call run_fmo2(z, sym, xyz, owner, opts, mpi_res, error, comm=world)
      if (error%has_error()) then
         if (world%leader()) print "(a,a)", "  FAIL distributed run: ", trim(error%get_message())
         call error%clear()
         n_bad = n_bad + 1
         return
      end if

      if (world%leader()) then
         write (text, "(a,a6,a,a4,a,f18.10,a,es10.2)") "  ", esp, "/", expansion, &
            "  E = ", mpi_res%energy, "   difference from serial ", &
            mpi_res%energy - serial_res%energy
         print "(a)", trim(text)
      end if

      if (abs(mpi_res%energy - serial_res%energy) > TOL) then
         if (world%leader()) print "(a)", "    FAIL: the distributed answer moved"
         n_bad = n_bad + 1
      end if
   end subroutine one_method

end program check_fmo_mpi
