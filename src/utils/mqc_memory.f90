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
   public :: set_memory_budget
   public :: memory_budget

   real(dp), save :: budget_override = -1.0_dp
      !! Bytes the deck said this run may plan on, `system.memory_gb`; negative
      !! is unset, and the machine decides.

contains

   subroutine set_memory_budget(gigabytes)
      !! Fix the budget from the deck; a non-positive value hands it back to the machine
      real(dp), intent(in) :: gigabytes
      if (gigabytes > 0.0_dp) then
         budget_override = gigabytes*1.0e9_dp
      else
         budget_override = -1.0_dp
      end if
   end subroutine set_memory_budget

   function memory_budget(blind, share) result(budget)
      !! Bytes one memory decision may plan on
      !!
      !! The deck's figure when it gave one, otherwise `share` of what the
      !! machine reports available, otherwise `blind` -- the caller's own
      !! constant for a machine that reports nothing. The deck's figure is
      !! taken whole: it is what the user has decided this run may have, and
      !! sharing it again would second-guess that.
      real(dp), intent(in) :: blind
      real(dp), intent(in) :: share
      real(dp) :: budget

      real(dp) :: available

      if (budget_override > 0.0_dp) then
         budget = budget_override
         return
      end if
      available = available_memory_bytes()
      if (available > 0.0_dp) then
         budget = share*available
      else
         budget = blind
      end if
   end function memory_budget

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
