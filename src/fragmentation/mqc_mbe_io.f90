module mqc_mbe_io
   use pic_types, only: int32, int64, dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_physical_fragment, only: physical_fragment_t, to_angstrom
   use mqc_elements, only: element_number_to_symbol
   implicit none
   private
   public :: print_fragment_xyz, print_detailed_breakdown, print_detailed_breakdown_json

contains

   function get_body_level_name(body_level) result(level_name)
      !! Map body level (n-mer) to descriptive name
      !! Supports up to decamers (10-mers), then falls back to "N-mers" format
      integer, intent(in) :: body_level
      character(len=32) :: level_name

      select case (body_level)
      case (1)
         level_name = "monomers"
      case (2)
         level_name = "dimers"
      case (3)
         level_name = "trimers"
      case (4)
         level_name = "tetramers"
      case (5)
         level_name = "pentamers"
      case (6)
         level_name = "hexamers"
      case (7)
         level_name = "heptamers"
      case (8)
         level_name = "octamers"
      case (9)
         level_name = "nonamers"
      case (10)
         level_name = "decamers"
      case default
         ! For levels > 10, use generic format
         write (level_name, '(i0,a)') body_level, "-mers"
      end select
   end function get_body_level_name

   subroutine print_fragment_xyz(fragment_idx, phys_frag)
      !! Print fragment geometry in XYZ format
      integer, intent(in) :: fragment_idx
      type(physical_fragment_t), intent(in) :: phys_frag
      integer :: i
      character(len=2) :: symbol
      character(len=256) :: coord_line

      call logger%info("=========================================")
      call logger%info(" Fragment "//to_char(fragment_idx))
      call logger%info(" Number of atoms: "//to_char(phys_frag%n_atoms))
      call logger%info(" Coordinates in Angstroms:")
      call logger%info("-----------------------------------------")
      do i = 1, phys_frag%n_atoms
         symbol = element_number_to_symbol(phys_frag%element_numbers(i))
         ! Convert from Bohr back to Angstroms for printing
         write (coord_line, '(a2,3f15.8)') symbol, to_angstrom(phys_frag%coordinates(1:3, i))
         call logger%info(trim(coord_line))
      end do
      call logger%info("=========================================")

   end subroutine print_fragment_xyz

   subroutine print_detailed_breakdown(polymers, fragment_count, max_level, energies, delta_energies)
      !! Print detailed energy breakdown for each fragment
      !! Shows full energy and deltaE correction for all monomers, dimers, trimers, etc.
      !! Uses int64 for fragment_count to handle large fragment counts that overflow int32.
      integer, intent(in) :: polymers(:, :), max_level
      integer(int64), intent(in) :: fragment_count
      real(dp), intent(in) :: energies(:), delta_energies(:)

      integer(int64) :: i
      integer :: fragment_size, j, body_level
      character(len=512) :: fragment_str, energy_line
      integer(int64) :: count_by_level

      call logger%verbose(" ")
      call logger%verbose("============================================")
      call logger%verbose("Detailed Energy Breakdown by Fragment")
      call logger%verbose("============================================")

      ! Warn if we have very high fragmentation levels
      if (max_level > 10) then
         call logger%warning("Fragment levels exceed decamers (10-mers). Using generic N-mers notation.")
      end if

      do body_level = 1, max_level
         count_by_level = 0_int64

         do i = 1_int64, fragment_count
            fragment_size = count(polymers(i, :) > 0)
            if (fragment_size == body_level) count_by_level = count_by_level + 1_int64
         end do

         if (count_by_level > 0_int64) then
            call logger%verbose(" ")
            block
               character(len=256) :: header
               character(len=32) :: level_name
               level_name = get_body_level_name(body_level)
               write (header, '(a,a,i0,a)') trim(level_name), " (", count_by_level, " fragments):"
               ! Capitalize first letter
               if (len_trim(level_name) > 0) then
                  if (level_name(1:1) >= 'a' .and. level_name(1:1) <= 'z') then
                     header(1:1) = achar(iachar(header(1:1)) - 32)
                  end if
               end if
               call logger%verbose(trim(header))
            end block
            call logger%verbose("--------------------------------------------")

            do i = 1_int64, fragment_count
               fragment_size = count(polymers(i, :) > 0)

               if (fragment_size == body_level) then
                  fragment_str = "["
                  do j = 1, fragment_size
                     if (j > 1) then
                        write (fragment_str, '(a,a,i0)') trim(fragment_str), ",", polymers(i, j)
                     else
                        write (fragment_str, '(a,i0)') trim(fragment_str), polymers(i, j)
                     end if
                  end do
                  write (fragment_str, '(a,a)') trim(fragment_str), "]"

                  if (body_level == 1) then
                     write (energy_line, '(a,a,f20.10)') &
                        "  Fragment ", trim(adjustl(fragment_str)), energies(i)
                  else
                     write (energy_line, '(a,a,f20.10,a,f20.10)') &
                        "  Fragment ", trim(adjustl(fragment_str)), energies(i), &
                        "   deltaE: ", delta_energies(i)
                  end if
                  call logger%verbose(trim(energy_line))
               end if
            end do
         end if
      end do

      call logger%verbose(" ")
      call logger%verbose("============================================")

   end subroutine print_detailed_breakdown

   subroutine print_detailed_breakdown_json(polymers, fragment_count, max_level, &
                                            energies, delta_energies, sum_by_level, total_energy)
      !! Write detailed energy breakdown to results.json file
      !! Outputs structured JSON with all fragment energies and deltaE corrections
      !! Uses int64 for fragment_count to handle large fragment counts that overflow int32.
      integer, intent(in) :: polymers(:, :), max_level
      integer(int64), intent(in) :: fragment_count
      real(dp), intent(in) :: energies(:), delta_energies(:)
      real(dp), intent(in) :: sum_by_level(:), total_energy

      integer(int64) :: i
      integer :: fragment_size, j, body_level, unit, io_stat
      character(len=512) :: json_line
      integer(int64) :: count_by_level
      logical :: first_level, first_fragment
      character(len=32) :: level_name

      open (newunit=unit, file='results.json', status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open results.json for writing")
         return
      end if

      call logger%info("Writing JSON output to results.json")

      ! Warn if we have very high fragmentation levels
      if (max_level > 10) then
         call logger%warning("Fragment levels exceed decamers (10-mers). JSON will use generic N-mers notation.")
      end if

      write (unit, '(a)') "{"
      write (unit, '(a)') '  "mbe_breakdown": {'

      write (json_line, '(a,f20.10,a)') '    "total_energy": ', total_energy, ','
      write (unit, '(a)') trim(json_line)

      write (unit, '(a)') '    "levels": ['

      first_level = .true.
      do body_level = 1, max_level
         count_by_level = 0_int64

         do i = 1_int64, fragment_count
            fragment_size = count(polymers(i, :) > 0)
            if (fragment_size == body_level) count_by_level = count_by_level + 1_int64
         end do

         if (count_by_level > 0_int64) then
            if (.not. first_level) then
               write (unit, '(a)') '      },'
            end if
            first_level = .false.

            write (unit, '(a)') '      {'

            level_name = get_body_level_name(body_level)

            write (json_line, '(a,i0,a)') '        "body_level": ', body_level, ','
            write (unit, '(a)') trim(json_line)
            write (json_line, '(a,a,a)') '        "name": "', trim(level_name), '",'
            write (unit, '(a)') trim(json_line)
            write (json_line, '(a,i0,a)') '        "count": ', count_by_level, ','
            write (unit, '(a)') trim(json_line)
            write (json_line, '(a,f20.10,a)') '        "total_energy": ', sum_by_level(body_level), ','
            write (unit, '(a)') trim(json_line)
            write (unit, '(a)') '        "fragments": ['

            first_fragment = .true.
            do i = 1_int64, fragment_count
               fragment_size = count(polymers(i, :) > 0)

               if (fragment_size == body_level) then
                  if (.not. first_fragment) then
                     write (unit, '(a)') '          },'
                  end if
                  first_fragment = .false.

                  write (unit, '(a)') '          {'

                  json_line = '            "indices": ['
                  do j = 1, fragment_size
                     if (j > 1) then
                        write (json_line, '(a,a,i0)') trim(json_line), ', ', polymers(i, j)
                     else
                        write (json_line, '(a,i0)') trim(json_line), polymers(i, j)
                     end if
                  end do
                  write (json_line, '(a,a)') trim(json_line), '],'
                  write (unit, '(a)') trim(json_line)

                  write (json_line, '(a,f20.10)') '            "energy": ', energies(i)
                  if (body_level > 1) then
                     write (json_line, '(a,a)') trim(json_line), ','
                     write (unit, '(a)') trim(json_line)
                     write (json_line, '(a,f20.10)') '            "delta_energy": ', delta_energies(i)
                  end if
                  write (unit, '(a)') trim(json_line)
               end if
            end do

            if (.not. first_fragment) then
               write (unit, '(a)') '          }'
            end if
            write (unit, '(a)') '        ]'
         end if
      end do

      if (.not. first_level) then
         write (unit, '(a)') '      }'
      end if

      write (unit, '(a)') '    ]'
      write (unit, '(a)') '  }'
      write (unit, '(a)') '}'

      close (unit)
      call logger%info("JSON output written successfully to results.json")

   end subroutine print_detailed_breakdown_json

end module mqc_mbe_io
