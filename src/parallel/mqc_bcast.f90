!! Broadcasts of the shapes pic-mpi does not collectively provide
module mqc_bcast
   !! Integer arrays and character strings, from rank 0 to everyone.
   !!
   !! **Linear, not a tree.** Rank 0 sends to every other rank in turn, which
   !! is O(P) messages where a tree would be O(log P). Setup-time only; the
   !! wrong shape if it ever moves inside a loop.
   !!
   !! Every routine is collective: all ranks call it, rank 0 supplies the data
   !! and the others receive it. A rank that skips one leaves the others
   !! waiting on a message that never comes, which is a hang rather than an
   !! error -- so these are called from paths where every rank arrives, not
   !! from inside a rank-guarded branch.
   use pic_types, only: dp, int32
   use pic_mpi_lib, only: comm_t, send, recv, bcast, MPI_Status
   use mqc_mpi_tags, only: TAG_BCAST_PAYLOAD
   implicit none
   private

   public :: bcast_integer_array   !! Allocatable integer array to every rank
   public :: bcast_real_array      !! Allocatable double array to every rank
   public :: bcast_string          !! Deferred-length character to every rank

contains

   subroutine bcast_integer_array(comm, data, root)
      !! Send an integer array from `root` to every other rank
      !!
      !! The length travels with the message, and `data` is reallocated on the
      !! receiving ranks to whatever arrived. An empty array arrives as
      !! `size 0`, allocated.
      type(comm_t), intent(in) :: comm
      integer(int32), allocatable, intent(inout) :: data(:)
      integer(int32), intent(in) :: root

      integer(int32) :: n, irank
      type(MPI_Status) :: status

      ! The count goes first as a scalar, so a receiving rank knows whether
      ! there is a payload at all before it tries to take one.
      n = 0
      if (comm%rank() == root .and. allocated(data)) n = int(size(data), int32)
      call bcast(comm, n, 1_int32, root)

      if (n == 0) then
         if (comm%rank() /= root) then
            if (allocated(data)) deallocate (data)
            allocate (data(0))
         end if
         return
      end if

      if (comm%rank() == root) then
         do irank = 0, comm%size() - 1
            if (irank == root) cycle
            call send(comm, data, irank, TAG_BCAST_PAYLOAD)
         end do
      else
         if (allocated(data)) deallocate (data)
         call recv(comm, data, root, TAG_BCAST_PAYLOAD, status)
      end if
   end subroutine bcast_integer_array

   subroutine bcast_real_array(comm, data, root)
      !! Send a double array from `root` to every other rank
      !!
      !! The length travels with the data, so the receiving ranks need not have
      !! sized their buffer in advance.
      type(comm_t), intent(in) :: comm
      real(dp), allocatable, intent(inout) :: data(:)
      integer(int32), intent(in) :: root

      integer(int32) :: n

      n = 0
      if (comm%rank() == root .and. allocated(data)) n = int(size(data), int32)
      call bcast(comm, n, 1_int32, root)

      if (comm%rank() /= root) then
         if (allocated(data)) deallocate (data)
         allocate (data(n))
      end if
      if (n > 0) call bcast(comm, data, n, root)
   end subroutine bcast_real_array

   subroutine bcast_string(comm, text, root)
      !! Send a deferred-length string from `root` to every other rank
      !!
      !! pic moves no characters, so the text travels as one integer per
      !! character. The wrong way to move anything large.
      type(comm_t), intent(in) :: comm
      character(len=:), allocatable, intent(inout) :: text
      integer(int32), intent(in) :: root

      integer(int32), allocatable :: codes(:)
      integer(int32) :: i, n

      if (comm%rank() == root) then
         n = 0
         if (allocated(text)) n = int(len(text), int32)
         allocate (codes(n))
         do i = 1, n
            codes(i) = int(iachar(text(i:i)), int32)
         end do
      else
         allocate (codes(0))
      end if

      call bcast_integer_array(comm, codes, root)

      if (comm%rank() /= root) then
         n = int(size(codes), int32)
         if (allocated(text)) deallocate (text)
         allocate (character(len=n) :: text)
         do i = 1, n
            text(i:i) = achar(codes(i))
         end do
      end if
      deallocate (codes)
   end subroutine bcast_string

end module mqc_bcast
