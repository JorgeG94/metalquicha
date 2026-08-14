!! Stage timing, and a report that says where a method spent its wall clock
module mqc_timing
   !! One accumulator per named stage, and one table at the end.
   !!
   !! **Why this exists.** A direct-SCF cholesterol run spends 88% of its wall
   !! clock in Fock builds and 0.4% in diagonalisation. Nothing in the program
   !! said so until it was measured, and "which part is slow" is not a question
   !! that should need a rebuild to answer. Every method that takes long enough
   !! to notice should be able to answer it for itself.
   !!
   !! **Wall, not CPU.** `pic_timer` reads `omp_get_wtime` when PIC is built
   !! threaded and `system_clock` otherwise -- both are elapsed real time. CPU
   !! seconds would report thread count multiplied by duration, which is the
   !! wrong number for every threaded stage here.
   !!
   !! **The shape is lap timing.** One clock runs from `start` to `finish`, and
   !! each `lap` attributes everything since the previous lap to a named stage:
   !!
   !! ```fortran
   !! type(timing_report_t) :: clk
   !! call clk%start()
   !! call build_integrals(...)
   !! call clk%lap("integrals")
   !! do iter = 1, n
   !!    call amplitude_equations(...)
   !!    call clk%lap("amplitudes")     ! accumulates, one entry not n
   !! end do
   !! call clk%finish()
   !! call clk%report("CCSD", verbose)
   !! ```
   !!
   !! Stages are created on first mention, so a caller adds one by naming it and
   !! nothing has to be declared up front. Laps accumulate, so a stage inside a
   !! loop reports its total and its call count rather than the last pass.
   !!
   !! **Nothing is silently unaccounted.** `report` prints the residual between
   !! the total and the sum of the stages as `other`. A missing `lap` shows up
   !! there rather than vanishing, which is the property that makes the table
   !! worth believing.
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use pic_logger, only: logger => global_logger
   use mqc_program_limits, only: MAX_LINE_LENGTH
   implicit none
   private

   public :: timing_report_t
   public :: MAX_TIMED_STAGES, MAX_STAGE_NAME

   !> Most stages one method may name. Exceeding it drops the extras into
   !> `other` rather than failing: a timer must never be the reason a
   !> calculation stops.
   integer, parameter :: MAX_TIMED_STAGES = 24
   !> Longest stage label kept. Longer names are truncated for the table.
   integer, parameter :: MAX_STAGE_NAME = 28
   !> Rule width, matching the row format: 4 indent + name + 12 + 10 + 14.
   integer, parameter :: TABLE_WIDTH = 4 + MAX_STAGE_NAME + 12 + 10 + 14

   type :: stage_t
      !! One named accumulator
      character(len=MAX_STAGE_NAME) :: name = ""
      real(dp) :: seconds = 0.0_dp
      integer :: calls = 0
   end type stage_t

   type :: timing_report_t
      !! Lap timer over named stages, plus the report
      private
      type(stage_t) :: stage(MAX_TIMED_STAGES)
      integer :: n_stages = 0
      type(timer_type) :: lap_clock       !! Reset at every lap
      type(timer_type) :: total_clock     !! Runs start to finish
      real(dp) :: total_seconds = 0.0_dp
      logical :: running = .false.
      logical :: overflowed = .false.
   contains
      procedure :: start => report_start
      procedure :: lap => report_lap
      procedure :: finish => report_finish
      procedure :: report => report_emit
      procedure :: seconds_of => report_seconds_of
      procedure :: total => report_total
   end type timing_report_t

contains

   subroutine report_start(self)
      !! Begin timing. Clears any previous stages.
      class(timing_report_t), intent(inout) :: self

      self%n_stages = 0
      self%total_seconds = 0.0_dp
      self%overflowed = .false.
      self%stage = stage_t()
      call self%total_clock%start()
      call self%lap_clock%start()
      self%running = .true.
   end subroutine report_start

   subroutine report_lap(self, name)
      !! Attribute the time since the previous lap to `name`
      !!
      !! Creates the stage on first mention and accumulates thereafter, so a
      !! call inside a loop reports the total for the loop and how many passes
      !! it took.
      class(timing_report_t), intent(inout) :: self
      character(len=*), intent(in) :: name

      real(dp) :: elapsed
      integer :: i, slot

      ! A lap before start would read an unstarted clock, which pic_timer stops
      ! the program for. Silently doing nothing is the right failure here.
      if (.not. self%running) return

      call self%lap_clock%stop()
      elapsed = self%lap_clock%get_elapsed_time()
      call self%lap_clock%start()

      slot = 0
      do i = 1, self%n_stages
         if (self%stage(i)%name == name) then
            slot = i
            exit
         end if
      end do

      if (slot == 0) then
         if (self%n_stages >= MAX_TIMED_STAGES) then
            ! Out of slots. The time is not lost -- it stays inside the total and
            ! surfaces as `other` -- but the report says the table is incomplete.
            self%overflowed = .true.
            return
         end if
         self%n_stages = self%n_stages + 1
         slot = self%n_stages
         self%stage(slot)%name = name
      end if

      self%stage(slot)%seconds = self%stage(slot)%seconds + elapsed
      self%stage(slot)%calls = self%stage(slot)%calls + 1
   end subroutine report_lap

   subroutine report_finish(self)
      !! Stop the total clock. Stages are whatever the laps recorded.
      class(timing_report_t), intent(inout) :: self

      if (.not. self%running) return
      call self%lap_clock%stop()
      call self%total_clock%stop()
      self%total_seconds = self%total_clock%get_elapsed_time()
      self%running = .false.
   end subroutine report_finish

   pure function report_total(self) result(seconds)
      !! Wall seconds between `start` and `finish`
      class(timing_report_t), intent(in) :: self
      real(dp) :: seconds

      seconds = self%total_seconds
   end function report_total

   pure function report_seconds_of(self, name) result(seconds)
      !! Wall seconds accumulated against one stage, or zero if never lapped
      !!
      !! For a caller that wants to assert on a stage rather than print it --
      !! a regression test on where the time goes, say.
      class(timing_report_t), intent(in) :: self
      character(len=*), intent(in) :: name
      real(dp) :: seconds

      integer :: i

      seconds = 0.0_dp
      do i = 1, self%n_stages
         if (self%stage(i)%name == name) then
            seconds = self%stage(i)%seconds
            return
         end if
      end do
   end function report_seconds_of

   subroutine report_emit(self, title, verbose)
      !! Print the table
      !!
      !! At `performance` level rather than `info`: these are the numbers a
      !! tuning run is for, and keeping them on their own level means they can
      !! be turned up or down without moving everything else with them.
      !!
      !! `verbose` is separate from the log level on purpose. The atomic guess
      !! runs a full SCF per element; those are runs nobody wants a table for.
      !! The level decides how loud a reporting run is, `verbose` decides
      !! whether this particular run reports at all.
      class(timing_report_t), intent(in) :: self
      character(len=*), intent(in) :: title
      logical, intent(in), optional :: verbose

      character(len=MAX_LINE_LENGTH) :: line
      real(dp) :: accounted, other
      integer :: i
      logical :: talk

      talk = .true.
      if (present(verbose)) talk = verbose
      if (.not. talk) return
      if (self%n_stages == 0 .and. self%total_seconds <= 0.0_dp) return

      accounted = 0.0_dp
      do i = 1, self%n_stages
         accounted = accounted + self%stage(i)%seconds
      end do
      other = self%total_seconds - accounted

      call logger%performance("")
      call logger%performance("  "//trim(title)//" timings")
      call logger%performance("  "//repeat("-", TABLE_WIDTH))
      ! Column widths match the row format below: name, then 12 for the seconds
      ! (f10.2 plus " s"), 10 for the call count, 14 for the per-call figure.
      write (line, "(a,a,a12,a10,a14)") "    ", pad_name("stage"), "total", "calls", "per call"
      call logger%performance(trim(line))
      call logger%performance("  "//repeat("-", TABLE_WIDTH))
      do i = 1, self%n_stages
         write (line, "(a,a,f10.2,a,i10,f12.4,a)") "    ", &
            pad_name(self%stage(i)%name), self%stage(i)%seconds, " s", &
            self%stage(i)%calls, &
            self%stage(i)%seconds/real(max(self%stage(i)%calls, 1), dp), " s"
         call logger%performance(trim(line))
      end do
      ! Always shown, including when it is zero: the parts visibly summing to the
      ! whole is what makes a missing lap detectable rather than invisible.
      write (line, "(a,a,f10.2,a)") "    ", pad_name("other"), other, " s"
      call logger%performance(trim(line))
      call logger%performance("  "//repeat("-", TABLE_WIDTH))
      write (line, "(a,a,f10.2,a)") "    ", pad_name("total"), self%total_seconds, " s"
      call logger%performance(trim(line))
      if (self%overflowed) then
         call logger%warning("  timing: more than "//trim(int_str(MAX_TIMED_STAGES))// &
                             " stages were named; the excess is inside 'other'")
      end if
      call logger%performance("")
   end subroutine report_emit

   pure function pad_name(name) result(padded)
      !! Stage label in a fixed column, truncated if it will not fit
      character(len=*), intent(in) :: name
      character(len=MAX_STAGE_NAME) :: padded

      padded = name
   end function pad_name

   pure function int_str(value) result(text)
      !! Small integer to text, so the warning above needs no io module
      integer, intent(in) :: value
      character(len=12) :: text

      write (text, "(i0)") value
   end function int_str

end module mqc_timing
