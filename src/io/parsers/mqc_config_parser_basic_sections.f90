!! Basic section parsers for MQC config files
!! Handles: schema, model, driver, scf, xtb, system sections
submodule(mqc_config_parser) mqc_config_parser_basic_sections
   implicit none

contains

   module subroutine parse_schema_section(unit, config, error)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      integer :: io_stat, eq_pos

      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_PARSE, "Unexpected end of file in %schema section")
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
         case ('name')
            config%schema_name = trim(value)
         case ('version')
            config%schema_version = trim(value)
         case ('index_base')
            read (value, *, iostat=io_stat) config%index_base
         case ('units')
            config%units = trim(value)
         case default
            call error%set(ERROR_PARSE, "Unknown key in %schema section: "//trim(key))
            return
         end select
      end do

   end subroutine parse_schema_section

   module subroutine parse_model_section(unit, config, error)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      integer :: io_stat, eq_pos
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %model section")
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
         case ('method')
            ! Parse method string (e.g., "XTB-GFN1" -> "gfn1")
            config%method = parse_method_string(trim(value))
            if (config%method == METHOD_TYPE_UNKNOWN) then
               call error%set(ERROR_PARSE, "Invalid method: "//trim(value))
               return
            end if
         case ('basis')
            config%basis = trim(value)
         case ('aux_basis')
            config%aux_basis = trim(value)
         case default
            call error%set(ERROR_PARSE, "Unknown key in %model section: "//trim(key))
            return
         end select
      end do

   end subroutine parse_model_section

   module subroutine parse_driver_section(unit, config, error)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      integer :: io_stat, eq_pos
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %driver section")
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
         case ('type')
            config%calc_type = calc_type_from_string(trim(value))
            if (config%calc_type == CALC_TYPE_UNKNOWN) then
               call error%set(ERROR_PARSE, "Invalid calc_type: "//trim(value))
               return
            end if
         case default
            call error%set(ERROR_PARSE, "Unknown key in %driver section: "//trim(key))
            return
         end select
      end do

   end subroutine parse_driver_section

   module subroutine parse_scf_section(unit, config, error)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      integer :: io_stat, eq_pos
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %scf section")
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
         case ('maxiter')
            read (value, *, iostat=io_stat) config%scf_maxiter
         case ('tolerance')
            read (value, *, iostat=io_stat) config%scf_tolerance
         case default
            call error%set(ERROR_PARSE, "Unknown key in %scf section: "//trim(key))
            return
         end select
      end do

   end subroutine parse_scf_section

   module subroutine parse_xtb_section(unit, config, error)
      !! Parse %xtb section for XTB-specific settings (solvation, etc.)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      integer :: io_stat, eq_pos
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %xtb section")
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
         case ('solvent')
            config%solvent = trim(value)
         case ('solvation_model')
            config%solvation_model = trim(value)
         case ('use_cds')
            config%use_cds = (trim(value) == 'true')
         case ('use_shift')
            config%use_shift = (trim(value) == 'true')
         case ('dielectric')
            read (value, *, iostat=io_stat) config%dielectric
            if (io_stat /= 0) then
               call error%set(ERROR_PARSE, "Invalid dielectric value: "//trim(value))
               return
            end if
         case ('cpcm_nang')
            read (value, *, iostat=io_stat) config%cpcm_nang
            if (io_stat /= 0) then
               call error%set(ERROR_PARSE, "Invalid cpcm_nang value: "//trim(value))
               return
            end if
         case ('cpcm_rscale')
            read (value, *, iostat=io_stat) config%cpcm_rscale
            if (io_stat /= 0) then
               call error%set(ERROR_PARSE, "Invalid cpcm_rscale value: "//trim(value))
               return
            end if
         case default
            call error%set(ERROR_PARSE, "Unknown key in %xtb section: "//trim(key))
            return
         end select
      end do

   end subroutine parse_xtb_section

   module subroutine parse_system_section(unit, config, error)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      integer :: io_stat, eq_pos
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %system section")
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
         case ('log_level')
            config%log_level = trim(value)
         case ('skip_json_output')
            config%skip_json_output = (trim(value) == 'true' .or. trim(value) == '.true.')
         case default
            call error%set(ERROR_PARSE, "Unknown key in %system section: "//trim(key))
            return
         end select
      end do

   end subroutine parse_system_section

end submodule mqc_config_parser_basic_sections
