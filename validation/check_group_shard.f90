!! Round-trip the group shard assignment across the wire
!!
!!     mpirun -np 2 ./build/check_group_shard
!!
!! `send_group_assignment_matrix` and its receiving half move a group's slice
!! of the task table from the global coordinator to a group leader. They are
!! two ends of one contract and nothing else checks it, which is how GMBE came
!! to pass the transposed matrix for as long as the helper existed: group 0's
!! own shard is handed over by `move_alloc` and never crosses the wire, so a
!! single-node run -- every run anyone has done -- never exercised the pair at
!! all.
!!
!! Two ranks are enough. The topology that hid the bug needs two *nodes*, but
!! the helper itself only needs a sender and a receiver, which is why this can
!! run here when the multi-group path cannot.
!!
!! The matrix is deliberately **not square**. With `n_tasks == width` the
!! transpose is the same shape as the original and the buffer arithmetic comes
!! out right by accident -- the very case that would let the fault through.
program check_group_shard
   use pic_types, only: int32, int64, default_int
   use pic_mpi_lib, only: comm_world, pic_mpi_init, pic_mpi_finalize, comm_t
   use mqc_group_shard_io, only: send_group_assignment_matrix, &
                                 receive_group_assignment_matrix
   implicit none

   !> Deliberately different, and deliberately width < tasks: the shape a
   !> real shard has, and the one the buffer arithmetic gets wrong.
   integer, parameter :: N_TASKS = 7
   integer, parameter :: WIDTH = 4

   integer(int64), allocatable :: ids(:), got_ids(:)
   integer, allocatable :: matrix(:, :), got(:, :)
   integer :: i, j, my_rank, nranks, failures
   type(comm_t) :: world

   call pic_mpi_init()
   world = comm_world()
   my_rank = world%rank()
   nranks = world%size()
   failures = 0

   if (nranks < 2) then
      if (my_rank == 0) then
         write (*, "(A)") "check_group_shard needs at least 2 ranks:"
         write (*, "(A)") "    mpirun -np 2 ./build/check_group_shard"
      end if
      call pic_mpi_finalize()
      stop 0
   end if

   if (my_rank == 0) then
      allocate (ids(N_TASKS), matrix(N_TASKS, WIDTH))
      do i = 1, N_TASKS
         ids(i) = int(100 + i, int64)
         do j = 1, WIDTH
            ! Distinct per (task, column) so a transpose cannot survive undetected
            matrix(i, j) = 1000*i + j
         end do
      end do
      call send_group_assignment_matrix(world, 1, ids, matrix)
      write (*, "(A,I0,A,I0,A)") "rank 0 sent ", N_TASKS, " tasks x ", WIDTH, " wide"
   end if

   if (my_rank == 1) then
      call receive_group_assignment_matrix(world, got_ids, got)

      if (size(got_ids) /= N_TASKS) then
         write (*, "(A,I0,A,I0)") "FAIL ids: got ", size(got_ids), " want ", N_TASKS
         failures = failures + 1
      end if
      if (size(got, 1) /= N_TASKS .or. size(got, 2) /= WIDTH) then
         write (*, "(A,I0,A,I0,A,I0,A,I0)") "FAIL shape: got ", size(got, 1), "x", &
            size(got, 2), " want ", N_TASKS, "x", WIDTH
         failures = failures + 1
      else
         do i = 1, N_TASKS
            if (got_ids(i) /= int(100 + i, int64)) then
               write (*, "(A,I0)") "FAIL id at ", i
               failures = failures + 1
            end if
            do j = 1, WIDTH
               if (got(i, j) /= 1000*i + j) then
                  write (*, "(A,I0,A,I0,A,I0,A,I0)") "FAIL element (", i, ",", j, &
                     "): got ", got(i, j), " want ", 1000*i + j
                  failures = failures + 1
               end if
            end do
         end do
      end if

      if (failures == 0) then
         write (*, "(A)") "rank 1 received the shard intact: ids, shape and every element"
         write (*, "(A)") ""
         write (*, "(A)") "  ok  the (task, width) contract round-trips"
      else
         write (*, "(A,I0,A)") "", failures, " failure(s)"
      end if
   end if

   call pic_mpi_finalize()
   if (failures > 0) error stop 1
end program check_group_shard
