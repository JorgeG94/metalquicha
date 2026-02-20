!! Helpers for sending/receiving group shard assignments
module mqc_group_shard_io
   use pic_types, only: int32, int64
   use pic_mpi_lib, only: comm_t, isend, irecv, recv, wait, request_t, MPI_Status
   use mqc_mpi_tags, only: TAG_GROUP_ASSIGN, TAG_GROUP_POLYMERS
   implicit none
   private

   public :: send_group_assignment_matrix, receive_group_assignment_matrix

contains

   subroutine send_group_assignment_matrix(world_comm, dest_rank, ids, matrix)
      type(comm_t), intent(in) :: world_comm
      integer, intent(in) :: dest_rank
      integer(int64), intent(in) :: ids(:)
      integer, intent(in) :: matrix(:, :)

      integer(int64) :: n_rows
      integer(int32) :: n_cols
      integer, allocatable :: buf(:)
      type(request_t) :: req(4)

      n_rows = size(ids, kind=int64)
      n_cols = size(matrix, 2)

      call isend(world_comm, n_rows, dest_rank, TAG_GROUP_ASSIGN, req(1))
      call isend(world_comm, ids, dest_rank, TAG_GROUP_ASSIGN, req(2))
      call isend(world_comm, n_cols, dest_rank, TAG_GROUP_POLYMERS, req(3))

      if (n_rows > 0_int64 .and. n_cols > 0) then
         allocate (buf(n_rows*n_cols))
         buf = reshape(matrix, [n_rows*n_cols])
      else
         allocate (buf(0))
      end if
      call isend(world_comm, buf, dest_rank, TAG_GROUP_POLYMERS, req(4))

      call wait(req(1))
      call wait(req(2))
      call wait(req(3))
      call wait(req(4))
      deallocate (buf)
   end subroutine send_group_assignment_matrix

   subroutine receive_group_assignment_matrix(world_comm, ids, matrix)
      type(comm_t), intent(in) :: world_comm
      integer(int64), allocatable, intent(out) :: ids(:)
      integer, allocatable, intent(out) :: matrix(:, :)

      integer(int64) :: n_rows
      integer(int32) :: n_cols
      integer, allocatable :: buf(:)
      type(MPI_Status) :: status
      type(request_t) :: req

      call irecv(world_comm, n_rows, 0, TAG_GROUP_ASSIGN, req)
      call wait(req)
      allocate (ids(n_rows))
      call recv(world_comm, ids, 0, TAG_GROUP_ASSIGN, status)

      call irecv(world_comm, n_cols, 0, TAG_GROUP_POLYMERS, req)
      call wait(req)
      allocate (matrix(int(n_rows), n_cols))

      if (n_rows > 0_int64 .and. n_cols > 0) then
         allocate (buf(n_rows*n_cols))
         call recv(world_comm, buf, 0, TAG_GROUP_POLYMERS, status)
         matrix = reshape(buf, [int(n_rows), n_cols])
         deallocate (buf)
      else
         allocate (buf(0))
         call recv(world_comm, buf, 0, TAG_GROUP_POLYMERS, status)
         deallocate (buf)
      end if
   end subroutine receive_group_assignment_matrix

end module mqc_group_shard_io
