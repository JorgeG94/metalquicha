module mqc_mbe_io
   use pic_types, only: int32, int64, dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_physical_fragment, only: physical_fragment_t, to_angstrom
   use mqc_elements, only: element_number_to_symbol
   use mqc_result_types, only: calculation_result_t, mbe_result_t
   use mqc_io_helpers, only: get_output_json_filename, get_basename
   use mqc_thermochemistry, only: thermochemistry_result_t
   use mqc_physical_constants, only: HARTREE_TO_CALMOL, R_CALMOLK
   implicit none
   private
   public :: print_fragment_xyz, print_detailed_breakdown, print_detailed_breakdown_json
   public :: print_unfragmented_json, print_gmbe_json, print_gmbe_pie_json
   public :: print_vibrational_json, print_vibrational_json_mbe

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
      integer(int64), intent(in) :: fragment_idx
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
                                            energies, delta_energies, sum_by_level, &
                                            mbe_result, results)
      !! Write detailed energy breakdown to results.json file
      !! Outputs structured JSON with all fragment energies and deltaE corrections
      !! Includes total gradient, Hessian, and dipole based on mbe_result%has_* flags
      !! Uses int64 for fragment_count to handle large fragment counts that overflow int32.
      use mqc_result_types, only: calculation_result_t, mbe_result_t
      integer, intent(in) :: polymers(:, :), max_level
      integer(int64), intent(in) :: fragment_count
      real(dp), intent(in) :: energies(:), delta_energies(:)
      real(dp), intent(in) :: sum_by_level(:)
      type(mbe_result_t), intent(in) :: mbe_result  !! MBE aggregated results
      type(calculation_result_t), intent(in), optional :: results(:)  !! Fragment results with distance info

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

      write (json_line, '(a,f20.10,a)') '    "total_energy": ', mbe_result%total_energy, ','
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

                  write (json_line, '(a,f20.10,a)') '            "energy": ', energies(i), ','
                  write (unit, '(a)') trim(json_line)

                  ! Add distance field if available
                  if (present(results)) then
                     write (json_line, '(a,f20.10)') '            "distance": ', results(i)%distance
                     if (frag_level > 1) then
                        write (json_line, '(a,a)') trim(json_line), ','
                        write (unit, '(a)') trim(json_line)
                        write (json_line, '(a,f20.10)') '            "delta_energy": ', delta_energies(i)
                     end if
                  else
                     if (frag_level > 1) then
                        write (json_line, '(a,f20.10)') '            "delta_energy": ', delta_energies(i)
                     end if
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

      ! Close levels array (with comma if we have more fields)
      if (mbe_result%has_gradient .or. mbe_result%has_hessian .or. mbe_result%has_dipole) then
         write (unit, '(a)') '    ],'
      else
         write (unit, '(a)') '    ]'
      end if

      ! Add dipole if computed (inside basename object)
      if (mbe_result%has_dipole) then
         write (unit, '(a)') '    "dipole": {'
         write (json_line, '(a,f25.15,a)') '      "x": ', mbe_result%dipole(1), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f25.15,a)') '      "y": ', mbe_result%dipole(2), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f25.15,a)') '      "z": ', mbe_result%dipole(3), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f25.15)') '      "magnitude_debye": ', norm2(mbe_result%dipole)*2.541746_dp
         write (unit, '(a)') trim(json_line)
         if (mbe_result%has_gradient .or. mbe_result%has_hessian) then
            write (unit, '(a)') '    },'
         else
            write (unit, '(a)') '    }'
         end if
      end if

      ! Add gradient norm if computed (inside basename object)
      if (mbe_result%has_gradient) then
         write (json_line, '(a,f20.10)') '    "gradient_norm": ', sqrt(sum(mbe_result%gradient**2))
         if (mbe_result%has_hessian) then
            write (json_line, '(a,a)') trim(json_line), ','
         end if
         write (unit, '(a)') trim(json_line)
      end if

      ! Add Hessian Frobenius norm if computed (inside basename object)
      if (mbe_result%has_hessian) then
         write (json_line, '(a,f20.10)') '    "hessian_frobenius_norm": ', sqrt(sum(mbe_result%hessian**2))
         write (unit, '(a)') trim(json_line)
      end if

      ! Close basename object
      write (unit, '(a)') '  }'
      ! Close outer object
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
         if (result%has_dipole .or. result%has_gradient .or. result%has_hessian) then
            write (json_line, '(a,a)') trim(json_line), ','
         end if
         write (unit, '(a)') trim(json_line)
      end if

      ! Add dipole if present
      if (result%has_dipole) then
         write (unit, '(a)') '    "dipole": {'
         write (json_line, '(a,f25.15,a)') '      "x": ', result%dipole(1), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f25.15,a)') '      "y": ', result%dipole(2), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f25.15,a)') '      "z": ', result%dipole(3), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f25.15)') '      "magnitude_debye": ', norm2(result%dipole)*2.541746_dp
         write (unit, '(a)') trim(json_line)
         if (result%has_gradient .or. result%has_hessian) then
            write (unit, '(a)') '    },'
         else
            write (unit, '(a)') '    }'
         end if
      end if

      ! Add gradient norm if present
      if (result%has_gradient) then
         write (json_line, '(a,f25.15)') '    "gradient_norm": ', sqrt(sum(result%gradient**2))
         if (result%has_hessian) then
            write (json_line, '(a,a)') trim(json_line), ','
         end if
         write (unit, '(a)') trim(json_line)
      end if

      ! Add Hessian Frobenius norm if present
      if (result%has_hessian) then
         write (json_line, '(a,f25.15)') '    "hessian_frobenius_norm": ', sqrt(sum(result%hessian**2))
         write (unit, '(a)') trim(json_line)
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

   subroutine print_gmbe_pie_json(pie_atom_sets, pie_coefficients, pie_energies, n_pie_terms, total_energy, &
                                  total_gradient, total_hessian)
      !! Write GMBE PIE calculation results to output JSON file
      !! Outputs structured JSON with PIE terms (atom sets with coefficients and energies)
      !! Optionally includes total gradient and Hessian norms
      integer, intent(in) :: pie_atom_sets(:, :)  !! Unique atom sets (max_atoms, n_pie_terms)
      integer, intent(in) :: pie_coefficients(:)  !! PIE coefficient for each term
      real(dp), intent(in) :: pie_energies(:)  !! Raw energy for each term
      integer(int64), intent(in) :: n_pie_terms
      real(dp), intent(in) :: total_energy
      real(dp), intent(in), optional :: total_gradient(:, :)  !! (3, total_atoms)
      real(dp), intent(in), optional :: total_hessian(:, :)  !! (3*total_atoms, 3*total_atoms)

      integer :: j, max_atoms, n_atoms
      integer(int64) :: i, n_nonzero_terms
      integer :: unit, io_stat
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

      ! Add gradient if present
      if (present(total_gradient)) then
         write (json_line, '(a,f20.10,a)') '    "gradient_norm": ', sqrt(sum(total_gradient**2)), ','
         write (unit, '(a)') trim(json_line)
      end if

      ! Add Hessian if present
      if (present(total_hessian)) then
         write (json_line, '(a,f20.10,a)') '    "hessian_frobenius_norm": ', sqrt(sum(total_hessian**2)), ','
         write (unit, '(a)') trim(json_line)
      end if

      ! PIE terms section
      ! First count non-zero coefficient terms
      n_nonzero_terms = 0_int64
      do i = 1_int64, n_pie_terms
         if (pie_coefficients(i) /= 0) n_nonzero_terms = n_nonzero_terms + 1
      end do

      write (unit, '(a)') '    "pie_terms": {'
      write (json_line, '(a,i0,a)') '      "count": ', n_nonzero_terms, ','
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '      "terms": ['

      max_atoms = size(pie_atom_sets, 1)
      first_term = .true.

      do i = 1_int64, n_pie_terms
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

   subroutine print_vibrational_json(result, frequencies, reduced_masses, force_constants, &
                                     thermo, ir_intensities)
      !! Write vibrational analysis and thermochemistry results to JSON file
      !! This extends the existing results.json with vibrational data
      type(calculation_result_t), intent(in) :: result
      real(dp), intent(in) :: frequencies(:)
      real(dp), intent(in) :: reduced_masses(:)
      real(dp), intent(in) :: force_constants(:)
      type(thermochemistry_result_t), intent(in) :: thermo
      real(dp), intent(in), optional :: ir_intensities(:)

      integer :: unit, io_stat, i, n_modes
      character(len=512) :: json_line
      character(len=256) :: output_file, basename
      logical :: has_ir

      output_file = get_output_json_filename()
      basename = get_basename()
      n_modes = size(frequencies)
      has_ir = present(ir_intensities)

      open (newunit=unit, file=trim(output_file), status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         return
      end if

      call logger%info("Writing vibrational/thermochemistry JSON to "//trim(output_file))

      write (unit, '(a)') "{"
      write (json_line, '(a,a,a)') '  "', trim(basename), '": {'
      write (unit, '(a)') trim(json_line)

      ! Total energy
      if (result%has_energy) then
         write (json_line, '(a,f25.15,a)') '    "total_energy": ', result%energy%total(), ','
         write (unit, '(a)') trim(json_line)
      end if

      ! Dipole
      if (result%has_dipole) then
         write (unit, '(a)') '    "dipole": {'
         write (json_line, '(a,f25.15,a)') '      "x": ', result%dipole(1), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f25.15,a)') '      "y": ', result%dipole(2), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f25.15,a)') '      "z": ', result%dipole(3), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f25.15)') '      "magnitude_debye": ', norm2(result%dipole)*2.541746_dp
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '    },'
      end if

      ! Gradient norm
      if (result%has_gradient) then
         write (json_line, '(a,f25.15,a)') '    "gradient_norm": ', sqrt(sum(result%gradient**2)), ','
         write (unit, '(a)') trim(json_line)
      end if

      ! Hessian norm
      if (result%has_hessian) then
         write (json_line, '(a,f25.15,a)') '    "hessian_frobenius_norm": ', sqrt(sum(result%hessian**2)), ','
         write (unit, '(a)') trim(json_line)
      end if

      ! Vibrational analysis section
      write (unit, '(a)') '    "vibrational_analysis": {'
      write (json_line, '(a,i0,a)') '      "n_modes": ', n_modes, ','
      write (unit, '(a)') trim(json_line)

      ! Frequencies array
      write (unit, '(a)', advance='no') '      "frequencies_cm1": ['
      do i = 1, n_modes
         if (i > 1) write (unit, '(a)', advance='no') ', '
         write (json_line, '(ES15.6)') frequencies(i)
         write (unit, '(a)', advance='no') trim(adjustl(json_line))
      end do
      write (unit, '(a)') '],'

      ! Reduced masses array
      write (unit, '(a)', advance='no') '      "reduced_masses_amu": ['
      do i = 1, n_modes
         if (i > 1) write (unit, '(a)', advance='no') ', '
         write (json_line, '(ES15.6)') reduced_masses(i)
         write (unit, '(a)', advance='no') trim(adjustl(json_line))
      end do
      write (unit, '(a)') '],'

      ! Force constants array
      write (unit, '(a)', advance='no') '      "force_constants_mdyne_ang": ['
      do i = 1, n_modes
         if (i > 1) write (unit, '(a)', advance='no') ', '
         write (json_line, '(ES15.6)') force_constants(i)
         write (unit, '(a)', advance='no') trim(adjustl(json_line))
      end do
      if (has_ir) then
         write (unit, '(a)') '],'
      else
         write (unit, '(a)') ']'
      end if

      ! IR intensities array (optional)
      if (has_ir) then
         write (unit, '(a)', advance='no') '      "ir_intensities_km_mol": ['
         do i = 1, n_modes
            if (i > 1) write (unit, '(a)', advance='no') ', '
            write (json_line, '(ES15.6)') ir_intensities(i)
            write (unit, '(a)', advance='no') trim(adjustl(json_line))
         end do
         write (unit, '(a)') ']'
      end if

      write (unit, '(a)') '    },'

      ! Thermochemistry section
      write (unit, '(a)') '    "thermochemistry": {'

      ! Conditions
      write (json_line, '(a,f10.4,a)') '      "temperature_K": ', thermo%temperature, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f10.4,a)') '      "pressure_atm": ', thermo%pressure, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6,a)') '      "molecular_mass_amu": ', thermo%total_mass, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,i0,a)') '      "symmetry_number": ', thermo%symmetry_number, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,i0,a)') '      "spin_multiplicity": ', thermo%spin_multiplicity, ','
      write (unit, '(a)') trim(json_line)
      if (thermo%is_linear) then
         write (unit, '(a)') '      "is_linear": true,'
      else
         write (unit, '(a)') '      "is_linear": false,'
      end if
      write (json_line, '(a,i0,a)') '      "n_real_frequencies": ', thermo%n_real_freqs, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,i0,a)') '      "n_imaginary_frequencies": ', thermo%n_imag_freqs, ','
      write (unit, '(a)') trim(json_line)

      ! Moments of inertia and rotational constants
      write (unit, '(a)') '      "moments_of_inertia_amu_ang2": {'
      write (json_line, '(a,f15.8,a)') '        "Ia": ', thermo%moments(1), ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f15.8,a)') '        "Ib": ', thermo%moments(2), ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f15.8)') '        "Ic": ', thermo%moments(3)
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '      },'

      write (unit, '(a)') '      "rotational_constants_GHz": {'
      write (json_line, '(a,f15.8,a)') '        "A": ', thermo%rot_const(1), ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f15.8,a)') '        "B": ', thermo%rot_const(2), ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f15.8)') '        "C": ', thermo%rot_const(3)
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '      },'

      ! Partition functions
      write (unit, '(a)') '      "partition_functions": {'
      write (json_line, '(a,ES15.6,a)') '        "translational": ', thermo%q_trans, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f15.6,a)') '        "rotational": ', thermo%q_rot, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f15.6)') '        "vibrational": ', thermo%q_vib
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '      },'

      ! Thermodynamic contributions
      write (unit, '(a)') '      "contributions": {'
      write (unit, '(a)') '        "translational": {'
      write (json_line, '(a,f18.12,a)') '          "energy_hartree": ', thermo%E_trans, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6,a)') '          "entropy_cal_mol_K": ', thermo%S_trans, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6)') '          "Cv_cal_mol_K": ', thermo%Cv_trans
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '        },'

      write (unit, '(a)') '        "rotational": {'
      write (json_line, '(a,f18.12,a)') '          "energy_hartree": ', thermo%E_rot, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6,a)') '          "entropy_cal_mol_K": ', thermo%S_rot, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6)') '          "Cv_cal_mol_K": ', thermo%Cv_rot
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '        },'

      write (unit, '(a)') '        "vibrational": {'
      write (json_line, '(a,f18.12,a)') '          "energy_hartree": ', thermo%E_vib, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6,a)') '          "entropy_cal_mol_K": ', thermo%S_vib, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6)') '          "Cv_cal_mol_K": ', thermo%Cv_vib
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '        },'

      write (unit, '(a)') '        "electronic": {'
      write (json_line, '(a,f18.12,a)') '          "energy_hartree": ', thermo%E_elec, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6)') '          "entropy_cal_mol_K": ', thermo%S_elec
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '        }'
      write (unit, '(a)') '      },'

      ! Contribution table (matching stdout format)
      block
         real(dp) :: pV_cal, H_vib_cal, H_rot_cal, H_trans_cal, H_total_cal
         real(dp) :: Cv_total, S_total, S_total_J
         real(dp) :: H_int_cal, Cv_int, S_int, S_int_J
         real(dp) :: Cp_trans

         pV_cal = R_CALMOLK*thermo%temperature
         H_vib_cal = thermo%E_vib*HARTREE_TO_CALMOL
         H_rot_cal = thermo%E_rot*HARTREE_TO_CALMOL
         H_trans_cal = thermo%E_trans*HARTREE_TO_CALMOL + pV_cal
         H_total_cal = H_vib_cal + H_rot_cal + H_trans_cal
         H_int_cal = H_vib_cal + H_rot_cal
         Cp_trans = thermo%Cv_trans + R_CALMOLK
         Cv_int = thermo%Cv_vib + thermo%Cv_rot
         Cv_total = Cp_trans + thermo%Cv_rot + thermo%Cv_vib
         S_int = thermo%S_vib + thermo%S_rot
         S_int_J = S_int*4.184_dp
         S_total = thermo%S_trans + thermo%S_rot + thermo%S_vib + thermo%S_elec
         S_total_J = S_total*4.184_dp

         write (unit, '(a)') '      "contribution_table": {'
         ! VIB row
         write (unit, '(a)') '        "VIB": {'
         write (json_line, '(a,f12.4,a)') '          "H_cal_mol": ', H_vib_cal, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "Cp_cal_mol_K": ', thermo%Cv_vib, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "S_cal_mol_K": ', thermo%S_vib, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4)') '          "S_J_mol_K": ', thermo%S_vib*4.184_dp
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '        },'

         ! ROT row
         write (unit, '(a)') '        "ROT": {'
         write (json_line, '(a,f12.4,a)') '          "H_cal_mol": ', H_rot_cal, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "Cp_cal_mol_K": ', thermo%Cv_rot, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "S_cal_mol_K": ', thermo%S_rot, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4)') '          "S_J_mol_K": ', thermo%S_rot*4.184_dp
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '        },'

         ! INT row (internal = VIB + ROT)
         write (unit, '(a)') '        "INT": {'
         write (json_line, '(a,f12.4,a)') '          "H_cal_mol": ', H_int_cal, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "Cp_cal_mol_K": ', Cv_int, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "S_cal_mol_K": ', S_int, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4)') '          "S_J_mol_K": ', S_int_J
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '        },'

         ! TR row (translational)
         write (unit, '(a)') '        "TR": {'
         write (json_line, '(a,f12.4,a)') '          "H_cal_mol": ', H_trans_cal, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "Cp_cal_mol_K": ', Cp_trans, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "S_cal_mol_K": ', thermo%S_trans, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4)') '          "S_J_mol_K": ', thermo%S_trans*4.184_dp
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '        },'

         ! TOT row (total)
         write (unit, '(a)') '        "TOT": {'
         write (json_line, '(a,f12.4,a)') '          "H_cal_mol": ', H_total_cal, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "Cp_cal_mol_K": ', Cv_total, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "S_cal_mol_K": ', S_total, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4)') '          "S_J_mol_K": ', S_total_J
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '        }'
         write (unit, '(a)') '      },'
      end block

      ! Zero-point energy
      write (json_line, '(a,f18.12,a)') '      "zero_point_energy_hartree": ', thermo%zpe_hartree, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6,a)') '      "zero_point_energy_kcal_mol": ', thermo%zpe_kcalmol, ','
      write (unit, '(a)') trim(json_line)

      ! Thermal corrections
      write (unit, '(a)') '      "thermal_corrections_hartree": {'
      write (json_line, '(a,f18.12,a)') '        "to_energy": ', thermo%thermal_correction_energy, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f18.12,a)') '        "to_enthalpy": ', thermo%thermal_correction_enthalpy, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f18.12)') '        "to_gibbs": ', thermo%thermal_correction_gibbs
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '      },'

      ! Total energies (with electronic energy)
      write (unit, '(a)') '      "total_energies_hartree": {'
      write (json_line, '(a,f20.12,a)') '        "electronic": ', result%energy%total(), ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f20.12,a)') '        "electronic_plus_zpe": ', &
         result%energy%total() + thermo%zpe_hartree, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f20.12,a)') '        "electronic_plus_thermal_E": ', &
         result%energy%total() + thermo%thermal_correction_energy, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f20.12,a)') '        "electronic_plus_thermal_H": ', &
         result%energy%total() + thermo%thermal_correction_enthalpy, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f20.12)') '        "electronic_plus_thermal_G": ', &
         result%energy%total() + thermo%thermal_correction_gibbs
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '      }'

      write (unit, '(a)') '    }'

      ! Close main object
      write (unit, '(a)') '  }'
      write (unit, '(a)') '}'

      close (unit)
      call logger%info("Vibrational/thermochemistry JSON written successfully to "//trim(output_file))

   end subroutine print_vibrational_json

   subroutine print_vibrational_json_mbe(mbe_result, frequencies, reduced_masses, force_constants, &
                                         thermo, ir_intensities)
      !! Write vibrational analysis and thermochemistry results to JSON file for MBE calculations
      !! This is analogous to print_vibrational_json but takes mbe_result_t instead of calculation_result_t
      type(mbe_result_t), intent(in) :: mbe_result
      real(dp), intent(in) :: frequencies(:)
      real(dp), intent(in) :: reduced_masses(:)
      real(dp), intent(in) :: force_constants(:)
      type(thermochemistry_result_t), intent(in) :: thermo
      real(dp), intent(in), optional :: ir_intensities(:)

      integer :: unit, io_stat, i, n_modes
      character(len=512) :: json_line
      character(len=256) :: output_file, basename
      logical :: has_ir

      output_file = get_output_json_filename()
      basename = get_basename()
      n_modes = size(frequencies)
      has_ir = present(ir_intensities)

      open (newunit=unit, file=trim(output_file), status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         return
      end if

      call logger%info("Writing vibrational/thermochemistry JSON to "//trim(output_file))

      write (unit, '(a)') "{"
      write (json_line, '(a,a,a)') '  "', trim(basename), '": {'
      write (unit, '(a)') trim(json_line)

      ! Total energy
      if (mbe_result%has_energy) then
         write (json_line, '(a,f25.15,a)') '    "total_energy": ', mbe_result%total_energy, ','
         write (unit, '(a)') trim(json_line)
      end if

      ! Dipole
      if (mbe_result%has_dipole) then
         write (unit, '(a)') '    "dipole": {'
         write (json_line, '(a,f25.15,a)') '      "x": ', mbe_result%dipole(1), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f25.15,a)') '      "y": ', mbe_result%dipole(2), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f25.15,a)') '      "z": ', mbe_result%dipole(3), ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f25.15)') '      "magnitude_debye": ', norm2(mbe_result%dipole)*2.541746_dp
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '    },'
      end if

      ! Gradient norm
      if (mbe_result%has_gradient) then
         write (json_line, '(a,f25.15,a)') '    "gradient_norm": ', sqrt(sum(mbe_result%gradient**2)), ','
         write (unit, '(a)') trim(json_line)
      end if

      ! Hessian norm
      if (mbe_result%has_hessian) then
         write (json_line, '(a,f25.15,a)') '    "hessian_frobenius_norm": ', sqrt(sum(mbe_result%hessian**2)), ','
         write (unit, '(a)') trim(json_line)
      end if

      ! Vibrational analysis section
      write (unit, '(a)') '    "vibrational_analysis": {'
      write (json_line, '(a,i0,a)') '      "n_modes": ', n_modes, ','
      write (unit, '(a)') trim(json_line)

      ! Frequencies array
      write (unit, '(a)', advance='no') '      "frequencies_cm1": ['
      do i = 1, n_modes
         if (i > 1) write (unit, '(a)', advance='no') ', '
         write (json_line, '(ES15.6)') frequencies(i)
         write (unit, '(a)', advance='no') trim(adjustl(json_line))
      end do
      write (unit, '(a)') '],'

      ! Reduced masses array
      write (unit, '(a)', advance='no') '      "reduced_masses_amu": ['
      do i = 1, n_modes
         if (i > 1) write (unit, '(a)', advance='no') ', '
         write (json_line, '(ES15.6)') reduced_masses(i)
         write (unit, '(a)', advance='no') trim(adjustl(json_line))
      end do
      write (unit, '(a)') '],'

      ! Force constants array
      write (unit, '(a)', advance='no') '      "force_constants_mdyne_ang": ['
      do i = 1, n_modes
         if (i > 1) write (unit, '(a)', advance='no') ', '
         write (json_line, '(ES15.6)') force_constants(i)
         write (unit, '(a)', advance='no') trim(adjustl(json_line))
      end do
      if (has_ir) then
         write (unit, '(a)') '],'
      else
         write (unit, '(a)') ']'
      end if

      ! IR intensities array (optional)
      if (has_ir) then
         write (unit, '(a)', advance='no') '      "ir_intensities_km_mol": ['
         do i = 1, n_modes
            if (i > 1) write (unit, '(a)', advance='no') ', '
            write (json_line, '(ES15.6)') ir_intensities(i)
            write (unit, '(a)', advance='no') trim(adjustl(json_line))
         end do
         write (unit, '(a)') ']'
      end if

      write (unit, '(a)') '    },'

      ! Thermochemistry section
      write (unit, '(a)') '    "thermochemistry": {'

      ! Conditions
      write (json_line, '(a,f10.4,a)') '      "temperature_K": ', thermo%temperature, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f10.4,a)') '      "pressure_atm": ', thermo%pressure, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6,a)') '      "molecular_mass_amu": ', thermo%total_mass, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,i0,a)') '      "symmetry_number": ', thermo%symmetry_number, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,i0,a)') '      "spin_multiplicity": ', thermo%spin_multiplicity, ','
      write (unit, '(a)') trim(json_line)
      if (thermo%is_linear) then
         write (unit, '(a)') '      "is_linear": true,'
      else
         write (unit, '(a)') '      "is_linear": false,'
      end if
      write (json_line, '(a,i0,a)') '      "n_real_frequencies": ', thermo%n_real_freqs, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,i0,a)') '      "n_imaginary_frequencies": ', thermo%n_imag_freqs, ','
      write (unit, '(a)') trim(json_line)

      ! Moments of inertia and rotational constants
      write (unit, '(a)') '      "moments_of_inertia_amu_ang2": {'
      write (json_line, '(a,f15.8,a)') '        "Ia": ', thermo%moments(1), ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f15.8,a)') '        "Ib": ', thermo%moments(2), ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f15.8)') '        "Ic": ', thermo%moments(3)
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '      },'

      write (unit, '(a)') '      "rotational_constants_GHz": {'
      write (json_line, '(a,f15.8,a)') '        "A": ', thermo%rot_const(1), ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f15.8,a)') '        "B": ', thermo%rot_const(2), ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f15.8)') '        "C": ', thermo%rot_const(3)
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '      },'

      ! Partition functions
      write (unit, '(a)') '      "partition_functions": {'
      write (json_line, '(a,ES15.6,a)') '        "translational": ', thermo%q_trans, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f15.6,a)') '        "rotational": ', thermo%q_rot, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f15.6)') '        "vibrational": ', thermo%q_vib
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '      },'

      ! Thermodynamic contributions
      write (unit, '(a)') '      "contributions": {'
      write (unit, '(a)') '        "translational": {'
      write (json_line, '(a,f18.12,a)') '          "energy_hartree": ', thermo%E_trans, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6,a)') '          "entropy_cal_mol_K": ', thermo%S_trans, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6)') '          "Cv_cal_mol_K": ', thermo%Cv_trans
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '        },'

      write (unit, '(a)') '        "rotational": {'
      write (json_line, '(a,f18.12,a)') '          "energy_hartree": ', thermo%E_rot, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6,a)') '          "entropy_cal_mol_K": ', thermo%S_rot, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6)') '          "Cv_cal_mol_K": ', thermo%Cv_rot
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '        },'

      write (unit, '(a)') '        "vibrational": {'
      write (json_line, '(a,f18.12,a)') '          "energy_hartree": ', thermo%E_vib, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6,a)') '          "entropy_cal_mol_K": ', thermo%S_vib, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6)') '          "Cv_cal_mol_K": ', thermo%Cv_vib
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '        },'

      write (unit, '(a)') '        "electronic": {'
      write (json_line, '(a,f18.12,a)') '          "energy_hartree": ', thermo%E_elec, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6)') '          "entropy_cal_mol_K": ', thermo%S_elec
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '        }'
      write (unit, '(a)') '      },'

      ! Contribution table (matching stdout format)
      block
         real(dp) :: pV_cal, H_vib_cal, H_rot_cal, H_trans_cal, H_total_cal
         real(dp) :: Cv_total, S_total, S_total_J
         real(dp) :: H_int_cal, Cv_int, S_int, S_int_J
         real(dp) :: Cp_trans

         pV_cal = R_CALMOLK*thermo%temperature
         H_vib_cal = thermo%E_vib*HARTREE_TO_CALMOL
         H_rot_cal = thermo%E_rot*HARTREE_TO_CALMOL
         H_trans_cal = thermo%E_trans*HARTREE_TO_CALMOL + pV_cal
         H_total_cal = H_vib_cal + H_rot_cal + H_trans_cal
         H_int_cal = H_vib_cal + H_rot_cal
         Cp_trans = thermo%Cv_trans + R_CALMOLK
         Cv_int = thermo%Cv_vib + thermo%Cv_rot
         Cv_total = Cp_trans + thermo%Cv_rot + thermo%Cv_vib
         S_int = thermo%S_vib + thermo%S_rot
         S_int_J = S_int*4.184_dp
         S_total = thermo%S_trans + thermo%S_rot + thermo%S_vib + thermo%S_elec
         S_total_J = S_total*4.184_dp

         write (unit, '(a)') '      "contribution_table": {'
         ! VIB row
         write (unit, '(a)') '        "VIB": {'
         write (json_line, '(a,f12.4,a)') '          "H_cal_mol": ', H_vib_cal, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "Cp_cal_mol_K": ', thermo%Cv_vib, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "S_cal_mol_K": ', thermo%S_vib, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4)') '          "S_J_mol_K": ', thermo%S_vib*4.184_dp
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '        },'

         ! ROT row
         write (unit, '(a)') '        "ROT": {'
         write (json_line, '(a,f12.4,a)') '          "H_cal_mol": ', H_rot_cal, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "Cp_cal_mol_K": ', thermo%Cv_rot, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "S_cal_mol_K": ', thermo%S_rot, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4)') '          "S_J_mol_K": ', thermo%S_rot*4.184_dp
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '        },'

         ! INT row (internal = VIB + ROT)
         write (unit, '(a)') '        "INT": {'
         write (json_line, '(a,f12.4,a)') '          "H_cal_mol": ', H_int_cal, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "Cp_cal_mol_K": ', Cv_int, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "S_cal_mol_K": ', S_int, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4)') '          "S_J_mol_K": ', S_int_J
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '        },'

         ! TR row (translational)
         write (unit, '(a)') '        "TR": {'
         write (json_line, '(a,f12.4,a)') '          "H_cal_mol": ', H_trans_cal, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "Cp_cal_mol_K": ', Cp_trans, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "S_cal_mol_K": ', thermo%S_trans, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4)') '          "S_J_mol_K": ', thermo%S_trans*4.184_dp
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '        },'

         ! TOT row (total)
         write (unit, '(a)') '        "TOT": {'
         write (json_line, '(a,f12.4,a)') '          "H_cal_mol": ', H_total_cal, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "Cp_cal_mol_K": ', Cv_total, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4,a)') '          "S_cal_mol_K": ', S_total, ','
         write (unit, '(a)') trim(json_line)
         write (json_line, '(a,f12.4)') '          "S_J_mol_K": ', S_total_J
         write (unit, '(a)') trim(json_line)
         write (unit, '(a)') '        }'
         write (unit, '(a)') '      },'
      end block

      ! Zero-point energy
      write (json_line, '(a,f18.12,a)') '      "zero_point_energy_hartree": ', thermo%zpe_hartree, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f12.6,a)') '      "zero_point_energy_kcal_mol": ', thermo%zpe_kcalmol, ','
      write (unit, '(a)') trim(json_line)

      ! Thermal corrections
      write (unit, '(a)') '      "thermal_corrections_hartree": {'
      write (json_line, '(a,f18.12,a)') '        "to_energy": ', thermo%thermal_correction_energy, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f18.12,a)') '        "to_enthalpy": ', thermo%thermal_correction_enthalpy, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f18.12)') '        "to_gibbs": ', thermo%thermal_correction_gibbs
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '      },'

      ! Total energies (with electronic energy)
      write (unit, '(a)') '      "total_energies_hartree": {'
      write (json_line, '(a,f20.12,a)') '        "electronic": ', mbe_result%total_energy, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f20.12,a)') '        "electronic_plus_zpe": ', &
         mbe_result%total_energy + thermo%zpe_hartree, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f20.12,a)') '        "electronic_plus_thermal_E": ', &
         mbe_result%total_energy + thermo%thermal_correction_energy, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f20.12,a)') '        "electronic_plus_thermal_H": ', &
         mbe_result%total_energy + thermo%thermal_correction_enthalpy, ','
      write (unit, '(a)') trim(json_line)
      write (json_line, '(a,f20.12)') '        "electronic_plus_thermal_G": ', &
         mbe_result%total_energy + thermo%thermal_correction_gibbs
      write (unit, '(a)') trim(json_line)
      write (unit, '(a)') '      }'

      write (unit, '(a)') '    }'

      ! Close main object
      write (unit, '(a)') '  }'
      write (unit, '(a)') '}'

      close (unit)
      call logger%info("Vibrational/thermochemistry JSON written successfully to "//trim(output_file))

   end subroutine print_vibrational_json_mbe

end module mqc_mbe_io
