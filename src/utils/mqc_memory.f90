module mqc_memory
   !! What the machine has room for, asked of the kernel rather than assumed
   !!
   !! Every memory decision in the program used to be a constant sized for a
   !! laptop, and on a node with five hundred gigabytes those constants sent a
   !! 545-function MakeFP down a fourteen-hour column build when the transform
   !! it refused would have taken a hundred gigabytes and a minute. A budget
   !! read from the machine is right on both.
   use pic_types, only: dp
   use mqc_program_limits, only: MAX_LINE_LENGTH
   implicit none
   private

   public :: available_memory_bytes

contains

   function available_memory_bytes() result(bytes)
      !! MemAvailable from /proc/meminfo, or zero where that does not exist
      !!
      !! MemAvailable rather than MemFree: free memory on a warm machine is
      !! almost nothing, because the kernel has spent it on page cache it will
      !! hand back on demand. Zero is "unknown", and every caller has a blind
      !! default for it; Linux is the only platform that answers. Not `pure`:
      !! it reads the machine.
      real(dp) :: bytes
      integer :: unit, stat
      character(len=MAX_LINE_LENGTH) :: line
      real(dp) :: kb

      bytes = 0.0_dp
      open (newunit=unit, file="/proc/meminfo", status="old", action="read", iostat=stat)
      if (stat /= 0) return
      do
         read (unit, "(a)", iostat=stat) line
         if (stat /= 0) exit
         if (line(1:13) == "MemAvailable:") then
            read (line(14:), *, iostat=stat) kb
            if (stat == 0) bytes = kb*1024.0_dp
            exit
         end if
      end do
      close (unit)
   end function available_memory_bytes

end module mqc_memory
