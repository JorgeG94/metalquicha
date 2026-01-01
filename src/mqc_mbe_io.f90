module mqc_mbe_io
   use pic_types, only: int32, int64, dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_physical_fragment, only: physical_fragment_t, to_angstrom
   use mqc_elements, only: element_number_to_symbol
   use mqc_result_types, only: calculation_result_t
   use mqc_io_helpers, only: get_output_json_filename, get_basename
   implicit none
   private
   public :: print_fragment_xyz, print_detailed_breakdown, print_detailed_breakdown_json
   public :: print_unfragmented_json, print_gmbe_json, print_gmbe_pie_json

contains

   function get_frag_level_name(frag_level) result(level_name)
      !! Map body level (n-mer) to descriptive name
      !! Supports up to decamers (10-mers), then falls back to "N-mers" format
      integer, intent(in) :: frag_level
      character(len=32) :: level_name

      select case (frag_level)
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
         write (level_name, '(i0,a)') frag_level, "-mers"
      end select
   end function get_frag_level_name

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
      integer :: fragment_size, j, frag_level
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

      do frag_level = 1, max_level
         count_by_level = 0_int64

         do i = 1_int64, fragment_count
            fragment_size = count(polymers(i, :) > 0)
            if (fragment_size == frag_level) count_by_level = count_by_level + 1_int64
         end do

         if (count_by_level > 0_int64) then
            call logger%verbose(" ")
            block
               character(len=256) :: header
               character(len=32) :: level_name
               level_name = get_frag_level_name(frag_level)
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

               if (fragment_size == frag_level) then
                  fragment_str = "["
                  do j = 1, fragment_size
                     if (j > 1) then
                        write (fragment_str, '(a,a,i0)') trim(fragment_str), ",", polymers(i, j)
                     else
                        write (fragment_str, '(a,i0)') trim(fragment_str), polymers(i, j)
                     end if
                  end do
                  write (fragment_str, '(a,a)') trim(fragment_str), "]"

                  if (frag_level == 1) then
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
                                            energies, delta_energies, sum_by_level, total_energy, &
                                            total_gradient)
      !! Write detailed energy breakdown to results.json file
      !! Outputs structured JSON with all fragment energies and deltaE corrections
      !! Optionally includes total gradient if provided
      !! Uses int64 for fragment_count to handle large fragment counts that overflow int32.
      integer, intent(in) :: polymers(:, :), max_level
      integer(int64), intent(in) :: fragment_count
      real(dp), intent(in) :: energies(:), delta_energies(:)
      real(dp), intent(in) :: sum_by_level(:), total_energy
      real(dp), intent(in), optional :: total_gradient(:, :)  !! (3, total_atoms)

      integer(int64) :: i
      integer :: fragment_size, j, frag_level, unit, io_stat, iatom
      character(len=512) :: json_line
      integer(int64) :: count_by_level
      logical :: first_level, first_fragment
      character(len=32) :: level_name
      integer :: total_atoms
      character(len=256) :: output_file, basename

      output_file = get_output_json_filename()
      basename = get_basename()

      open (newunit=unit, file=trim(output_file), status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         return
      end if

      call logger%info("Writing JSON output to "//trim(output_file))

      ! Warn if we have very high fragmentation levels
      if (max_level > 10) then
         call logger%warning("Fragment levels exceed decamers (10-mers). JSON will use generic N-mers notation.")
      end if

      write (unit, '(a)') "{"
      write (json_line, '(a,a,a)') '  "', trim(basename), '": {'
      write (unit, '(a)') trim(json_line)

      write (json_line, '(a,f20.10,a)') '    "total_energy": ', total_energy, ','
      write (unit, '(a)') trim(json_line)

      write (unit, '(a)') '    "levels": ['

      first_level = .true.
      do frag_level = 1, max_level
         count_by_level = 0_int64

         do i = 1_int64, fragment_count
            fragment_size = count(polymers(i, :) > 0)
            if (fragment_size == frag_level) count_by_level = count_by_level + 1_int64
         end do

         if (count_by_level > 0_int64) then
            if (.not. first_level) then
               write (unit, '(a)') '      },'
            end if
            first_level = .false.

            write (unit, '(a)') '      {'

            level_name = get_frag_level_name(frag_level)

            write (json_line, '(a,i0,a)') '        "frag_level": ', frag_level, ','
            write (unit, '(a)') trim(json_line)
            write (json_line, '(a,a,a)') '        "name": "', trim(level_name), '",'
            write (unit, '(a)') trim(json_line)
            write (json_line, '(a,i0,a)') '        "count": ', count_by_level, ','
            write (unit, '(a)') trim(json_line)
            write (json_line, '(a,f20.10,a)') '        "total_energy": ', sum_by_level(frag_level), ','
            write (unit, '(a)') trim(json_line)
            write (unit, '(a)') '        "fragments": ['

            first_fragment = .true.
            do i = 1_int64, fragment_count
               fragment_size = count(polymers(i, :) > 0)

               if (fragment_size == frag_level) then
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
                  if (frag_level > 1) then
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

      ! Close mbe_breakdown section with or without comma depending on gradient presence
      if (present(total_gradient)) then
         write (unit, '(a)') '  },'
      else
         write (unit, '(a)') '  }'
      end if

      ! Add gradient section if present
      if (present(total_gradient)) then
         total_atoms = size(total_gradient, 2)

         write (unit, '(a)') '  "gradient": {'
         write (json_line, '(a,f20.10,a)') '    "norm": ', sqrt(sum(total_gradient**2)), ','
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '    "components": ['

         do iatom = 1, total_atoms
            write (json_line, '(a,3(f20.12,a))') '      [', &
               total_gradient(1, iatom), ', ', &
               total_gradient(2, iatom), ', ', &
               total_gradient(3, iatom), ']'
            if (iatom < total_atoms) then
               write (json_line, '(a,a)') trim(json_line), ','
            end if
            write (unit, '(a)') trim(json_line)
         end do

         write (unit, '(a)') '    ]'
         write (unit, '(a)') '  }'
      end if

      write (unit, '(a)') '}'

      close (unit)
      call logger%info("JSON output written successfully to "//trim(output_file))

   end subroutine print_detailed_breakdown_json

   subroutine print_unfragmented_json(result)
      !! Write unfragmented calculation results to output JSON file
      !! Outputs structured JSON with energy and optionally gradient
      type(calculation_result_t), intent(in) :: result

      integer :: unit, io_stat, iatom, total_atoms
      character(len=512) :: json_line
      character(len=256) :: output_file, basename

      output_file = get_output_json_filename()
      basename = get_basename()

      open (newunit=unit, file=trim(output_file), status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         return
      end if

      call logger%info("Writing JSON output to "//trim(output_file))

      write (unit, '(a)') "{"
      write (json_line, '(a,a,a)') '  "', trim(basename), '": {'
      write (unit, '(a)') trim(json_line)

      if (result%has_energy) then
         write (json_line, '(a,f25.15)') '    "total_energy": ', result%energy%total()
         if (result%has_gradient) then
            write (json_line, '(a,a)') trim(json_line), ','
         end if
         write (unit, '(a)') trim(json_line)
      end if

      ! Add gradient section if present
      if (result%has_gradient) then
         total_atoms = size(result%gradient, 2)

         write (unit, '(a)') '    "gradient": {'
         write (json_line, '(a,f25.15,a)') '      "norm": ', sqrt(sum(result%gradient**2)), ','
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '      "components": ['

         do iatom = 1, total_atoms
            write (json_line, '(a,3(f25.15,a))') '        [', &
               result%gradient(1, iatom), ', ', &
               result%gradient(2, iatom), ', ', &
               result%gradient(3, iatom), ']'
            if (iatom < total_atoms) then
               write (json_line, '(a,a)') trim(json_line), ','
            end if
            write (unit, '(a)') trim(json_line)
         end do

         write (unit, '(a)') '      ]'
         write (unit, '(a)') '    }'
      end if

      write (unit, '(a)') '  }'
      write (unit, '(a)') '}'

      close (unit)
      call logger%info("JSON output written successfully to "//trim(output_file))

   end subroutine print_unfragmented_json

   subroutine print_gmbe_json(n_monomers, monomer_indices, monomer_results, &
                              n_intersections, intersection_results, &
                              intersection_sets, intersection_levels, total_energy)
      !! Write GMBE calculation results to output JSON file
      !! Outputs structured JSON with monomers, intersections, and total energy
      !! Intersection parameters are optional and should be omitted when n_intersections=0
      integer, intent(in) :: n_monomers
      integer, intent(in) :: monomer_indices(:)
      type(calculation_result_t), intent(in) :: monomer_results(:)
      integer, intent(in) :: n_intersections
      type(calculation_result_t), intent(in), optional :: intersection_results(:)
      integer, intent(in), optional :: intersection_sets(:, :)  !! (n_monomers, n_intersections)
      integer, intent(in), optional :: intersection_levels(:)
      real(dp), intent(in) :: total_energy

      integer :: i, j, k, max_level, unit, io_stat
      character(len=512) :: json_line
      character(len=256) :: output_file, basename
      logical :: first_level, first_intersection
      integer :: level_count

      output_file = get_output_json_filename()
      basename = get_basename()

      open (newunit=unit, file=trim(output_file), status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         return
      end if

      call logger%info("Writing GMBE JSON output to "//trim(output_file))

      write (unit, '(a)') "{"
      write (json_line, '(a,a,a)') '  "', trim(basename), '": {'
      write (unit, '(a)') trim(json_line)

      write (json_line, '(a,f20.10,a)') '    "total_energy": ', total_energy, ','
      write (unit, '(a)') trim(json_line)

      ! Monomers section
      write (unit, '(a)') '    "monomers": {'
      write (json_line, '(a,i0,a)') '      "count": ', n_monomers, ','
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '      "fragments": ['

      do i = 1, n_monomers
         write (unit, '(a)') '        {'
         write (json_line, '(a,i0,a)') '          "index": ', monomer_indices(i), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f20.10)') '          "energy": ', monomer_results(i)%energy%total()
         write (unit, '(a)') trim(json_line)
         if (i < n_monomers) then
            write (unit, '(a)') '        },'
         else
            write (unit, '(a)') '        }'
         end if
      end do

      write (unit, '(a)') '      ]'

      ! Add comma after monomers if we have intersections
      if (n_intersections > 0 .and. present(intersection_results) .and. &
          present(intersection_sets) .and. present(intersection_levels)) then
         write (unit, '(a)') '    },'
      else
         write (unit, '(a)') '    }'
      end if

      ! Intersections section
      if (n_intersections > 0 .and. present(intersection_results) .and. &
          present(intersection_sets) .and. present(intersection_levels)) then
         max_level = maxval(intersection_levels)

         write (unit, '(a)') '    "intersections": {'
         write (json_line, '(a,i0,a)') '      "total_count": ', n_intersections, ','
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '      "levels": ['

         first_level = .true.
         do k = 2, max_level
            ! Count intersections at this level
            level_count = 0
            do i = 1, n_intersections
               if (intersection_levels(i) == k) level_count = level_count + 1
            end do

            if (level_count > 0) then
               if (.not. first_level) then
                  write (unit, '(a)') '        },'
               end if
               first_level = .false.

               write (unit, '(a)') '        {'
               write (json_line, '(a,i0,a)') '          "level": ', k, ','
               write (unit, '(a)') trim(json_line)
               write (json_line, '(a,i0,a)') '          "count": ', level_count, ','
               write (unit, '(a)') trim(json_line)
               write (unit, '(a)') '          "fragments": ['

               first_intersection = .true.
               do i = 1, n_intersections
                  if (intersection_levels(i) == k) then
                     if (.not. first_intersection) then
                        write (unit, '(a)') '            },'
                     end if
                     first_intersection = .false.

                     write (unit, '(a)') '            {'

                     ! Write indices
                     json_line = '              "indices": ['
                     do j = 1, n_monomers
                        if (intersection_sets(j, i) > 0) then
                           if (j > 1 .and. intersection_sets(j - 1, i) > 0) then
                              write (json_line, '(a,a,i0)') trim(json_line), ', ', intersection_sets(j, i)
                           else
                              write (json_line, '(a,i0)') trim(json_line), intersection_sets(j, i)
                           end if
                        end if
                     end do
                     write (json_line, '(a,a)') trim(json_line), '],'
                     write (unit, '(a)') trim(json_line)

                     write (json_line, '(a,f20.10)') '              "energy": ', intersection_results(i)%energy%total()
                     write (unit, '(a)') trim(json_line)
                  end if
               end do

               if (.not. first_intersection) then
                  write (unit, '(a)') '            }'
               end if
               write (unit, '(a)') '          ]'
            end if
         end do

         if (.not. first_level) then
            write (unit, '(a)') '        }'
         end if

         write (unit, '(a)') '      ]'
         write (unit, '(a)') '    }'
      end if

      write (unit, '(a)') '  }'
      write (unit, '(a)') '}'

      close (unit)
      call logger%info("GMBE JSON output written successfully to "//trim(output_file))

   end subroutine print_gmbe_json

   subroutine print_gmbe_pie_json(pie_atom_sets, pie_coefficients, pie_energies, n_pie_terms, total_energy)
      !! Write GMBE PIE calculation results to output JSON file
      !! Outputs structured JSON with PIE terms (atom sets with coefficients and energies)
      integer, intent(in) :: pie_atom_sets(:, :)  !! Unique atom sets (max_atoms, n_pie_terms)
      integer, intent(in) :: pie_coefficients(:)  !! PIE coefficient for each term
      real(dp), intent(in) :: pie_energies(:)  !! Raw energy for each term
      integer, intent(in) :: n_pie_terms
      real(dp), intent(in) :: total_energy

      integer :: i, j, max_atoms, n_atoms
      integer :: unit, io_stat
      integer :: n_nonzero_terms
      logical :: first_term
      character(len=512) :: json_line
      character(len=256) :: output_file, basename

      output_file = get_output_json_filename()
      basename = get_basename()

      open (newunit=unit, file=trim(output_file), status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         return
      end if

      call logger%info("Writing GMBE PIE JSON output to "//trim(output_file))

      write (unit, '(a)') "{"
      write (json_line, '(a,a,a)') '  "', trim(basename), '": {'
      write (unit, '(a)') trim(json_line)

      write (json_line, '(a,f20.10,a)') '    "total_energy": ', total_energy, ','
      write (unit, '(a)') trim(json_line)

      ! PIE terms section
      ! First count non-zero coefficient terms
      n_nonzero_terms = 0
      do i = 1, n_pie_terms
         if (pie_coefficients(i) /= 0) n_nonzero_terms = n_nonzero_terms + 1
      end do

      write (unit, '(a)') '    "pie_terms": {'
      write (json_line, '(a,i0,a)') '      "count": ', n_nonzero_terms, ','
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '      "terms": ['

      max_atoms = size(pie_atom_sets, 1)
      first_term = .true.

      do i = 1, n_pie_terms
         ! Skip terms with zero coefficient
         if (pie_coefficients(i) == 0) cycle

         if (.not. first_term) write (unit, '(a)') '        },'
         first_term = .false.

         write (unit, '(a)') '        {'

         ! Extract atom list size
         n_atoms = 0
         do while (n_atoms < max_atoms .and. pie_atom_sets(n_atoms + 1, i) >= 0)
            n_atoms = n_atoms + 1
         end do

         ! Write atom indices
         json_line = '          "atom_indices": ['
         do j = 1, n_atoms
            if (j > 1) then
               write (json_line, '(a,a,i0)') trim(json_line), ', ', pie_atom_sets(j, i)
            else
               write (json_line, '(a,i0)') trim(json_line), pie_atom_sets(j, i)
            end if
         end do
         write (json_line, '(a,a)') trim(json_line), '],'
         write (unit, '(a)') trim(json_line)

         write (json_line, '(a,i0,a)') '          "coefficient": ', pie_coefficients(i), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f20.10,a)') '          "energy": ', pie_energies(i), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f20.10)') '          "weighted_energy": ', &
            real(pie_coefficients(i), dp)*pie_energies(i)
         write (unit, '(a)') trim(json_line)
      end do

      if (.not. first_term) write (unit, '(a)') '        }'

      write (unit, '(a)') '      ]'
      write (unit, '(a)') '    }'
      write (unit, '(a)') '  }'
      write (unit, '(a)') '}'

      close (unit)
      call logger%info("GMBE PIE JSON output written successfully to "//trim(output_file))

   end subroutine print_gmbe_pie_json

end module mqc_mbe_io
