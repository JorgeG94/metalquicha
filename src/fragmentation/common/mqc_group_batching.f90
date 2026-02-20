!! Shared batching helpers for group-global coordinator flows
module mqc_group_batching
   use pic_types, only: int32, int64
   use pic_timer, only: timer_type
   use pic_mpi_lib, only: comm_t, isend, irecv, recv, wait, iprobe, MPI_Status, request_t, MPI_ANY_SOURCE, abort_comm
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_mpi_tags, only: TAG_WORKER_SCALAR_RESULT, TAG_NODE_SCALAR_RESULT, TAG_GROUP_RESULT, TAG_GROUP_DONE
   use mqc_result_types, only: calculation_result_t, result_isend, result_irecv
   implicit none
   private

   public :: append_result_to_batch
   public :: flush_group_results
   public :: handle_local_worker_results_to_batch
   public :: handle_node_results_to_batch
   public :: handle_group_results

contains

   !! Append a completed fragment result to the current outbound batch.
   subroutine append_result_to_batch(item_idx, result, batch_count, batch_ids, batch_results)
      integer(int64), intent(in) :: item_idx
      type(calculation_result_t), intent(in) :: result
      integer(int32), intent(inout) :: batch_count
      integer(int64), intent(inout) :: batch_ids(:)
      type(calculation_result_t), intent(inout) :: batch_results(:)

      batch_count = batch_count + 1
      batch_ids(batch_count) = item_idx
      batch_results(batch_count) = result
   end subroutine append_result_to_batch

   !! Send all currently batched results to rank 0 and reset the batch.
   subroutine flush_group_results(world_comm, batch_count, batch_ids, batch_results)
      type(comm_t), intent(in) :: world_comm
      integer(int32), intent(inout) :: batch_count
      integer(int64), intent(inout) :: batch_ids(:)
      type(calculation_result_t), intent(inout) :: batch_results(:)

      type(request_t) :: req
      integer :: i

      if (batch_count <= 0) return

      call isend(world_comm, batch_count, 0, TAG_GROUP_RESULT, req)
      call wait(req)
      call isend(world_comm, batch_ids(1:batch_count), 0, TAG_GROUP_RESULT, req)
      call wait(req)
      do i = 1, batch_count
         call result_isend(batch_results(i), world_comm, 0, TAG_GROUP_RESULT, req)
         call wait(req)
         call batch_results(i)%destroy()
      end do
      batch_count = 0
   end subroutine flush_group_results

   !! Drain pending local worker results and append them to the outbound batch.
   subroutine handle_local_worker_results_to_batch(node_comm, world_comm, worker_map, batch_count, batch_ids, batch_results, &
                                                   results_received)
      type(comm_t), intent(in) :: node_comm
      type(comm_t), intent(in) :: world_comm
      integer(int64), intent(inout) :: worker_map(:)
      integer(int32), intent(inout) :: batch_count
      integer(int64), intent(inout) :: batch_ids(:)
      type(calculation_result_t), intent(inout) :: batch_results(:)
      integer(int64), intent(inout), optional :: results_received

      type(MPI_Status) :: local_status
      logical :: has_pending
      integer :: worker_source
      type(request_t) :: req
      type(calculation_result_t) :: worker_result
      integer(int64) :: item_idx

      if (node_comm%size() <= 1) return

      do
         call iprobe(node_comm, MPI_ANY_SOURCE, TAG_WORKER_SCALAR_RESULT, has_pending, local_status)
         if (.not. has_pending) exit

         worker_source = local_status%MPI_SOURCE

         if (worker_map(worker_source) == 0) then
            call logger%error("Received result from worker "//to_char(worker_source)// &
                              " but no item was assigned!")
            call abort_comm(world_comm, 1)
         end if

         call result_irecv(worker_result, node_comm, worker_source, TAG_WORKER_SCALAR_RESULT, req)
         call wait(req)

         if (worker_result%has_error) then
            call logger%error("Item "//to_char(worker_map(worker_source))// &
                              " calculation failed: "// &
                              worker_result%error%get_message())
            call abort_comm(world_comm, 1)
         end if

         item_idx = worker_map(worker_source)
         worker_map(worker_source) = 0

         if (batch_count >= size(batch_ids)) then
            call flush_group_results(world_comm, batch_count, batch_ids, batch_results)
         end if

         call append_result_to_batch(item_idx, worker_result, batch_count, batch_ids, batch_results)
         if (present(results_received)) results_received = results_received + 1_int64
         if (batch_count >= size(batch_ids)) then
            call flush_group_results(world_comm, batch_count, batch_ids, batch_results)
         end if
         call worker_result%destroy()
      end do
   end subroutine handle_local_worker_results_to_batch

   !! Drain pending node-level results and append them to the outbound batch.
   subroutine handle_node_results_to_batch(world_comm, batch_count, batch_ids, batch_results, results_received)
      type(comm_t), intent(in) :: world_comm
      integer(int32), intent(inout) :: batch_count
      integer(int64), intent(inout) :: batch_ids(:)
      type(calculation_result_t), intent(inout) :: batch_results(:)
      integer(int64), intent(inout), optional :: results_received

      integer(int64) :: item_idx
      type(MPI_Status) :: status
      logical :: has_pending
      type(request_t) :: req
      type(calculation_result_t) :: node_result

      do
         call iprobe(world_comm, MPI_ANY_SOURCE, TAG_NODE_SCALAR_RESULT, has_pending, status)
         if (.not. has_pending) exit

         call irecv(world_comm, item_idx, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT, req)
         call wait(req)
         call result_irecv(node_result, world_comm, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT, req)
         call wait(req)

         if (node_result%has_error) then
            call logger%error("Item "//to_char(item_idx)//" calculation failed: "// &
                              node_result%error%get_message())
            call abort_comm(world_comm, 1)
         end if

         if (batch_count >= size(batch_ids)) then
            call flush_group_results(world_comm, batch_count, batch_ids, batch_results)
         end if

         call append_result_to_batch(item_idx, node_result, batch_count, batch_ids, batch_results)
         if (present(results_received)) results_received = results_received + 1_int64
         if (batch_count >= size(batch_ids)) then
            call flush_group_results(world_comm, batch_count, batch_ids, batch_results)
         end if
         call node_result%destroy()
      end do
   end subroutine handle_node_results_to_batch

   !! Receive grouped result batches on rank 0 and update global progress counters.
   subroutine handle_group_results(world_comm, results, results_received, total_items, coord_timer, group_done_count, label)
      type(comm_t), intent(in) :: world_comm
      type(calculation_result_t), intent(inout) :: results(:)
      integer(int64), intent(inout) :: results_received
      integer(int64), intent(in) :: total_items
      type(timer_type), intent(in) :: coord_timer
      integer, intent(inout) :: group_done_count
      character(len=*), intent(in), optional :: label

      integer(int32) :: batch_count
      integer(int64), allocatable :: batch_ids(:)
      type(MPI_Status) :: status
      logical :: has_pending
      type(request_t) :: req
      integer :: i, dummy_msg
      character(len=32) :: item_label

      if (present(label)) then
         item_label = label
      else
         item_label = "item"
      end if

      do
         call iprobe(world_comm, MPI_ANY_SOURCE, TAG_GROUP_RESULT, has_pending, status)
         if (.not. has_pending) exit

         call irecv(world_comm, batch_count, status%MPI_SOURCE, TAG_GROUP_RESULT, req)
         call wait(req)
         if (batch_count <= 0) cycle

         allocate (batch_ids(batch_count))
         call recv(world_comm, batch_ids, status%MPI_SOURCE, TAG_GROUP_RESULT, status)
         do i = 1, batch_count
            call result_irecv(results(batch_ids(i)), world_comm, status%MPI_SOURCE, TAG_GROUP_RESULT, req)
            call wait(req)

            if (results(batch_ids(i))%has_error) then
               call logger%error(trim(item_label)//" "//to_char(batch_ids(i))//" calculation failed: "// &
                                 results(batch_ids(i))%error%get_message())
               call abort_comm(world_comm, 1)
            end if

            results_received = results_received + 1
            if (mod(results_received, max(1_int64, total_items/10_int64)) == 0 .or. &
                results_received == total_items) then
               call logger%info("  Processed "//to_char(results_received)//"/"// &
                                to_char(total_items)//" "//trim(item_label)//"s ["// &
                                to_char(coord_timer%get_elapsed_time())//" s]")
            end if
         end do
         deallocate (batch_ids)
      end do

      do
         call iprobe(world_comm, MPI_ANY_SOURCE, TAG_GROUP_DONE, has_pending, status)
         if (.not. has_pending) exit
         call irecv(world_comm, dummy_msg, status%MPI_SOURCE, TAG_GROUP_DONE, req)
         call wait(req)
         group_done_count = group_done_count + 1
      end do
   end subroutine handle_group_results

end module mqc_group_batching
