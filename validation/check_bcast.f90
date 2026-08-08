!! Manual check that the setup-time broadcasts reach every rank intact
!!
!!     mpirun -np 4 ./build/check_bcast
!!
!! Every rank prints one line; a clean run ends with "ok" from all of them and
!! nothing else. There is no ctest harness for MPI here, so this is run by hand.
!!
!! It is worth running at -np 1 as well. A broadcast with nobody to send to
!! must still leave rank 0's own data alone and must not block waiting on
!! itself, and that path is easy to get wrong while the multi-rank one works.
program check_bcast
   use pic_types, only: dp, int32
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init, pic_mpi_finalize
   use mqc_bcast, only: bcast_integer_array, bcast_real_array, bcast_string
   implicit none

   type(comm_t) :: comm
   integer(int32), allocatable :: ints(:)
   integer(int32), allocatable :: empty(:)
   real(dp), allocatable :: reals(:)
   character(len=:), allocatable :: text
   integer(int32) :: rank, i, fails

   call pic_mpi_init()
   comm = comm_world()
   rank = comm%rank()
   fails = 0

   ! ---- an integer array: the term lists and partitions travel this way -----
   if (rank == 0) then
      allocate (ints(1000))
      do i = 1, 1000
         ints(i) = i*7
      end do
   end if
   call bcast_integer_array(comm, ints, 0_int32)
   if (.not. allocated(ints)) then
      fails = fails + 1
   else if (size(ints) /= 1000) then
      fails = fails + 1
   else if (any(ints /= [(i*7, i=1, 1000)])) then
      fails = fails + 1
   end if

   ! ---- an empty one, which is a legitimate result of a strict screen -------
   if (rank == 0) allocate (empty(0))
   call bcast_integer_array(comm, empty, 0_int32)
   if (.not. allocated(empty)) fails = fails + 1
   if (allocated(empty)) then
      if (size(empty) /= 0) fails = fails + 1
   end if

   ! ---- a real array: coordinates -------------------------------------------
   if (rank == 0) then
      allocate (reals(300))
      do i = 1, 300
         reals(i) = real(i, dp)*0.5_dp
      end do
   end if
   call bcast_real_array(comm, reals, 0_int32)
   if (.not. allocated(reals)) then
      fails = fails + 1
   else if (size(reals) /= 300) then
      fails = fails + 1
   else if (maxval(abs(reals - [(real(i, dp)*0.5_dp, i=1, 300)])) > 0.0_dp) then
      ! Exact: a broadcast copies bits, it does not compute.
      fails = fails + 1
   end if

   ! ---- a string: the config document ---------------------------------------
   if (rank == 0) text = '{"schema": {"name": "t", "version": "1.0"}}'
   call bcast_string(comm, text, 0_int32)
   if (.not. allocated(text)) then
      fails = fails + 1
   else if (text /= '{"schema": {"name": "t", "version": "1.0"}}') then
      fails = fails + 1
   end if

   if (fails == 0) then
      write (*, "(A,I0,A)") "[rank ", rank, "] ok"
   else
      write (*, "(A,I0,A,I0,A)") "[rank ", rank, "] FAILED ", fails, " checks"
   end if

   call comm%finalize()
   call pic_mpi_finalize()
   if (fails /= 0) error stop 1
end program check_bcast
