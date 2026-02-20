!! Simple int64 queue utilities shared across fragmentation schemes
module mqc_work_queue
   use pic_types, only: int64
   implicit none
   private

   public :: queue_t
   public :: queue_init_from_list, queue_pop, queue_is_empty, queue_destroy

   type :: queue_t
      integer(int64), allocatable :: ids(:)
      integer(int64) :: head = 1
      integer(int64) :: count = 0
   end type queue_t

contains

   subroutine queue_init_from_list(queue, ids)
      type(queue_t), intent(out) :: queue
      integer(int64), intent(in) :: ids(:)

      queue%count = size(ids, kind=int64)
      if (queue%count > 0) then
         allocate (queue%ids(queue%count))
         queue%ids = ids
      end if
      queue%head = 1_int64
   end subroutine queue_init_from_list

   subroutine queue_pop(queue, item_idx, has_item)
      type(queue_t), intent(inout) :: queue
      integer(int64), intent(out) :: item_idx
      logical, intent(out) :: has_item

      if (queue%head > queue%count) then
         item_idx = -1_int64
         has_item = .false.
         return
      end if

      item_idx = queue%ids(queue%head)
      queue%head = queue%head + 1_int64
      has_item = .true.
   end subroutine queue_pop

   pure logical function queue_is_empty(queue)
      type(queue_t), intent(in) :: queue
      queue_is_empty = (queue%head > queue%count)
   end function queue_is_empty

   subroutine queue_destroy(queue)
      type(queue_t), intent(inout) :: queue
      if (allocated(queue%ids)) deallocate (queue%ids)
      queue%head = 1_int64
      queue%count = 0_int64
   end subroutine queue_destroy

end module mqc_work_queue
