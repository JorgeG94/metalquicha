!! Calculation settings section parsers for MQC config files
!! Handles: hessian, aimd, fragmentation sections
submodule(mqc_config_parser) mqc_config_parser_calc_settings
   implicit none

contains

   module subroutine parse_hessian_section(unit, config, error)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      integer :: io_stat, eq_pos
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %hessian section")
            return
         end if

         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle

         if (trim(strip_comment(line)) == 'end') exit

         eq_pos = index(line, '=')
         if (eq_pos == 0) cycle

         key = adjustl(line(1:eq_pos - 1))
         value = adjustl(line(eq_pos + 1:))

         select case (trim(key))
         case ('finite_difference_displacement', 'displacement')
            read (value, *, iostat=io_stat) config%hessian_displacement
         case ('temperature')
            read (value, *, iostat=io_stat) config%hessian_temperature
            if (io_stat /= 0) then
               call error%set(ERROR_PARSE, "Invalid temperature value: "//trim(value))
               return
            end if
         case ('pressure')
            read (value, *, iostat=io_stat) config%hessian_pressure
            if (io_stat /= 0) then
               call error%set(ERROR_PARSE, "Invalid pressure value: "//trim(value))
               return
            end if
         case default
            call error%set(ERROR_PARSE, "Unknown key in %hessian section: "//trim(key))
            return
         end select
      end do

   end subroutine parse_hessian_section

   module subroutine parse_aimd_section(unit, config, error)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      integer :: io_stat, eq_pos
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %aimd section")
            return
         end if

         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle

         if (trim(strip_comment(line)) == 'end') exit

         eq_pos = index(line, '=')
         if (eq_pos == 0) cycle

         key = adjustl(line(1:eq_pos - 1))
         value = adjustl(line(eq_pos + 1:))

         select case (trim(key))
         case ('dt', 'timestep')
            read (value, *, iostat=io_stat) config%aimd_dt
         case ('nsteps', 'steps')
            read (value, *, iostat=io_stat) config%aimd_nsteps
         case ('initial_temperature', 'temperature')
            read (value, *, iostat=io_stat) config%aimd_initial_temperature
         case ('output_frequency', 'output_freq')
            read (value, *, iostat=io_stat) config%aimd_output_frequency
         case default
            call error%set(ERROR_PARSE, "Unknown key in %aimd section: "//trim(key))
            return
         end select
      end do

   end subroutine parse_aimd_section

   module subroutine parse_fragmentation_section(unit, config, error)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      integer :: io_stat, eq_pos
      logical :: in_cutoffs
      in_cutoffs = .false.

      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %fragmentation section")
            return
         end if

         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle

         if (trim(strip_comment(line)) == 'end') then
            if (in_cutoffs) then
               ! Validate cutoffs before leaving the cutoffs section
               call validate_cutoffs(config, error)
               if (error%has_error()) then
                  call error%add_context("mqc_config_parser:parse_fragmentation_section")
                  return
               end if
               in_cutoffs = .false.
               cycle
            else
               exit
            end if
         end if

         if (trim(line) == '%cutoffs') then
            in_cutoffs = .true.
            cycle
         end if

         eq_pos = index(line, '=')
         if (eq_pos == 0) cycle

         key = adjustl(line(1:eq_pos - 1))
         value = adjustl(line(eq_pos + 1:))

         if (in_cutoffs) then
            ! Parse cutoffs: numeric keys like "2", "3", "4", etc.
            ! representing n-mer level (2=dimer, 3=trimer, etc.)
            block
               integer :: nmer_level
               real(dp) :: cutoff_value

               ! Try to read the key as an integer (n-mer level)
               read (key, *, iostat=io_stat) nmer_level
               if (io_stat /= 0) then
                  call error%set(ERROR_PARSE, "Invalid n-mer level in cutoffs (expected integer): "//trim(key))
                  return
               end if

               ! Validate n-mer level
               if (nmer_level < 2) then
                  call error%set(ERROR_PARSE, "N-mer level must be >= 2 in cutoffs")
                  return
               end if

               if (nmer_level > 10) then
                  call error%set(ERROR_PARSE, "N-mer level too large in cutoffs (max 10 for decamer)")
                  return
               end if

               ! Read the cutoff value
               read (value, *, iostat=io_stat) cutoff_value
               if (io_stat /= 0) then
                  call error%set(ERROR_PARSE, "Invalid cutoff value: "//trim(value))
                  return
               end if

               ! Allocate array if not yet allocated (up to decamer = 10)
               if (.not. allocated(config%fragment_cutoffs)) then
                  allocate (config%fragment_cutoffs(10))
                  config%fragment_cutoffs = -1.0_dp  ! Initialize with sentinel value
               end if

               ! Store the cutoff value at the appropriate index
               config%fragment_cutoffs(nmer_level) = cutoff_value
            end block
         else
            select case (trim(key))
            case ('method')
               config%frag_method = trim(value)
            case ('level')
               read (value, *, iostat=io_stat) config%frag_level
               if (io_stat == 0) then
                  if (config%frag_level < 0) then
                     call error%set(ERROR_VALIDATION, "Fragmentation level must be >= 0 (0 = unfragmented)")
                     return
                  end if
                  if (config%frag_level > 10) then
                     call error%set(ERROR_VALIDATION, &
                                    "Fragmentation level must be <= 10 (decamers). Higher levels not supported.")
                     return
                  end if
               end if
            case ('allow_overlapping_fragments')
               config%allow_overlapping_fragments = (trim(value) == 'true')
            case ('max_intersection_level')
               read (value, *, iostat=io_stat) config%max_intersection_level
               if (io_stat == 0) then
                  if (config%max_intersection_level < 1) then
                     call error%set(ERROR_VALIDATION, "max_intersection_level must be >= 1")
                     return
                  end if
                  if (config%max_intersection_level > 10) then
                     call error%set(ERROR_VALIDATION, &
                                    "max_intersection_level must be <= 10 (decamers). Higher levels not supported.")
                     return
                  end if
               end if
            case ('embedding')
               config%embedding = trim(value)
            case ('cutoff_method')
               config%cutoff_method = trim(value)
            case ('distance_metric')
               config%distance_metric = trim(value)
            case default
               call error%set(ERROR_PARSE, "Unknown key in %fragmentation section: "//trim(key))
               return
            end select
         end if
      end do

   end subroutine parse_fragmentation_section

   subroutine validate_cutoffs(config, error)
      !! Validate that fragment cutoffs are monotonically decreasing
      !! For n-mer level N, cutoff(N) must be <= cutoff(N-1)
      type(mqc_config_t), intent(in) :: config
      type(error_t), intent(out) :: error

      integer :: i, level_low, level_high
      real(dp) :: cutoff_low, cutoff_high
      character(len=256) :: msg

      if (.not. allocated(config%fragment_cutoffs)) return

      ! Check monotonicity for consecutive levels with defined cutoffs
      do i = 2, size(config%fragment_cutoffs)
         level_low = i - 1
         level_high = i

         cutoff_low = config%fragment_cutoffs(level_low)
         cutoff_high = config%fragment_cutoffs(level_high)

         ! Skip if either cutoff is not defined (negative or zero sentinel value)
         if (cutoff_low <= 0.0_dp .or. cutoff_high <= 0.0_dp) cycle

         ! Validate monotonic decreasing
         if (cutoff_high > cutoff_low) then
            write (msg, '(a,i0,a,f0.2,a,i0,a,f0.2,a)') &
               "Fragment cutoffs must be monotonically decreasing: ", &
               level_high, "-mer cutoff (", cutoff_high, ") cannot be larger than ", &
               level_low, "-mer cutoff (", cutoff_low, "). Check %cutoffs section."
            call error%set(ERROR_PARSE, trim(msg))
            return
         end if
      end do

   end subroutine validate_cutoffs

end submodule mqc_config_parser_calc_settings
