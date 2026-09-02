!! Stage timing, and a report that says where a method spent its wall clock
module mqc_timing
   !! One accumulator per named stage, and one table at the end.
   !!
   !! **Wall, not CPU.** `pic_timer` reads `omp_get_wtime` when PIC is built
   !! threaded and `system_clock` otherwise; both are elapsed real time.
   !!
   !! **The shape is lap timing.** One clock runs from `start` to `finish`, and
   !! each `lap` attributes everything since the previous lap to a named stage.
   !! `begin` is the announcement and `lap` the measurement: `begin` names the
   !! stage on the way in, and the matching `lap` needs no argument.
   !!
   !! Stages are created on first mention, so a caller adds one by naming it and
   !! nothing has to be declared up front. Laps accumulate, so a stage inside a
   !! loop reports its total and its call count rather than the last pass.
   !!
   !! `report` prints the residual between the total and the sum of the stages
   !! as `other`, so a missing `lap` shows up there rather than vanishing.
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use pic_logger, only: logger => global_logger
   use mqc_program_limits, only: MAX_LINE_LENGTH
   implicit none
   private

   public :: timing_report_t
   public :: MAX_TIMED_STAGES, MAX_STAGE_NAME

   integer, parameter :: MAX_TIMED_STAGES = 24
      !! Most stages one method may name. Exceeding it drops the extras into
      !! `other` rather than failing.
   integer, parameter :: MAX_STAGE_NAME = 28
      !! Longest stage label kept. Longer names are truncated for the table.
   integer, parameter :: TABLE_WIDTH = 4 + MAX_STAGE_NAME + 12 + 10 + 14
      !! Rule width, matching the row format: 4 indent + name + 12 + 10 + 14.

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
      character(len=MAX_STAGE_NAME) :: pending = ""
         !! Stage opened by `begin` and not yet closed by `lap`
   contains
      procedure :: start => report_start
      procedure :: begin => report_begin
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
      self%pending = ""
      call self%total_clock%start()
      call self%lap_clock%start()
      self%running = .true.
   end subroutine report_start

   subroutine report_begin(self, name)
      !! Say that `name` is starting, and let the next `lap` close it
      !!
      !! Logged at `performance` level, the same level the table uses, which is
      !! numerically below `info` and so shows by default.
      !!
      !! Announcing does not touch the lap clock. `begin` is free to be called
      !! anywhere; only `lap` divides the timeline.
      ! TODO(mqc): the `verbose` argument to `report` suppresses the table but
      ! nothing suppresses these, so a run that asked for no timing still logs
      ! a line per stage.
      class(timing_report_t), intent(inout) :: self
      character(len=*), intent(in) :: name

      self%pending = name
      call logger%performance("    "//trim(name)//" ...")
   end subroutine report_begin

   subroutine report_lap(self, name)
      !! Attribute the time since the previous lap to `name`
      !!
      !! Creates the stage on first mention and accumulates thereafter, so a
      !! call inside a loop reports the total for the loop and how many passes
      !! it took.
      class(timing_report_t), intent(inout) :: self
      character(len=*), intent(in), optional :: name
         !! Omitted closes whatever `begin` last opened. A lap with neither a
         !! name nor a pending stage does not divide the timeline at all -- the
         !! clock keeps running and the time lands in the next stage that is
         !! named.

      real(dp) :: elapsed
      integer :: i, slot
      character(len=MAX_STAGE_NAME) :: label

      ! A lap before start would read an unstarted clock, which pic_timer stops
      ! the program for.
      if (.not. self%running) return

      if (present(name)) then
         label = name
      else
         label = self%pending
      end if
      self%pending = ""
      if (label == "") return

      call self%lap_clock%stop()
      elapsed = self%lap_clock%get_elapsed_time()
      call self%lap_clock%start()

      slot = 0
      do i = 1, self%n_stages
         if (self%stage(i)%name == label) then
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
         self%stage(slot)%name = label
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
      !! At `performance` level rather than `info`. `verbose` is separate from
      !! the log level: the level decides how loud a reporting run is,
      !! `verbose` whether this particular run reports at all.
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
      ! Always shown, including when it is zero, so the parts visibly sum to
      ! the whole and a missing lap is detectable.
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
      !! Small integer to text
      integer, intent(in) :: value
      character(len=12) :: text

      write (text, "(i0)") value
   end function int_str

end module mqc_timing
