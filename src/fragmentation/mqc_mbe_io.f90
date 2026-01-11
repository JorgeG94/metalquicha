module mqc_mbe_io
   use pic_types, only: int32, int64, dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_physical_fragment, only: physical_fragment_t, to_angstrom
   use mqc_elements, only: element_number_to_symbol
   use mqc_result_types, only: calculation_result_t, mbe_result_t
   use mqc_io_helpers, only: get_output_json_filename, get_basename
   use mqc_thermochemistry, only: thermochemistry_result_t
   use mqc_physical_constants, only: HARTREE_TO_CALMOL, R_CALMOLK, AU_TO_DEBYE, CAL_TO_J
   use json_module, only: json_core, json_value
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
      !! Uses json-fortran for clean, maintainable JSON output
      use mqc_result_types, only: calculation_result_t, mbe_result_t
      integer, intent(in) :: polymers(:, :), max_level
      integer(int64), intent(in) :: fragment_count
      real(dp), intent(in) :: energies(:), delta_energies(:)
      real(dp), intent(in) :: sum_by_level(:)
      type(mbe_result_t), intent(in) :: mbe_result
      type(calculation_result_t), intent(in), optional :: results(:)

      type(json_core) :: json
      type(json_value), pointer :: root, main_obj, levels_arr, level_obj, frags_arr, frag_obj
      type(json_value), pointer :: dipole_obj
      integer(int64) :: i, count_by_level
      integer :: fragment_size, j, frag_level, iunit, io_stat
      integer, allocatable :: indices(:)
      character(len=32) :: level_name
      character(len=256) :: output_file, basename

      output_file = get_output_json_filename()
      basename = get_basename()

      if (max_level > 10) then
         call logger%warning("Fragment levels exceed decamers (10-mers). JSON will use generic N-mers notation.")
      end if

      call json%initialize()
      call json%create_object(root, '')
      call json%create_object(main_obj, trim(basename))
      call json%add(root, main_obj)

      call json%add(main_obj, 'total_energy', mbe_result%total_energy)

      ! Build levels array
      call json%create_array(levels_arr, 'levels')
      call json%add(main_obj, levels_arr)

      do frag_level = 1, max_level
         count_by_level = 0_int64
         do i = 1_int64, fragment_count
            fragment_size = count(polymers(i, :) > 0)
            if (fragment_size == frag_level) count_by_level = count_by_level + 1_int64
         end do

         if (count_by_level > 0_int64) then
            call json%create_object(level_obj, '')
            call json%add(levels_arr, level_obj)

            level_name = get_frag_level_name(frag_level)
            call json%add(level_obj, 'frag_level', frag_level)
            call json%add(level_obj, 'name', trim(level_name))
            call json%add(level_obj, 'count', int(count_by_level))
            call json%add(level_obj, 'total_energy', sum_by_level(frag_level))

            call json%create_array(frags_arr, 'fragments')
            call json%add(level_obj, frags_arr)

            do i = 1_int64, fragment_count
               fragment_size = count(polymers(i, :) > 0)
               if (fragment_size == frag_level) then
                  call json%create_object(frag_obj, '')
                  call json%add(frags_arr, frag_obj)

                  allocate (indices(fragment_size))
                  indices = polymers(i, 1:fragment_size)
                  call json%add(frag_obj, 'indices', indices)
                  deallocate (indices)

                  call json%add(frag_obj, 'energy', energies(i))

                  if (present(results)) then
                     call json%add(frag_obj, 'distance', results(i)%distance)
                  end if
                  if (frag_level > 1) then
                     call json%add(frag_obj, 'delta_energy', delta_energies(i))
                  end if
               end if
            end do
         end if
      end do

      ! Add dipole if computed
      if (mbe_result%has_dipole) then
         call json%create_object(dipole_obj, 'dipole')
         call json%add(main_obj, dipole_obj)
         call json%add(dipole_obj, 'x', mbe_result%dipole(1))
         call json%add(dipole_obj, 'y', mbe_result%dipole(2))
         call json%add(dipole_obj, 'z', mbe_result%dipole(3))
         call json%add(dipole_obj, 'magnitude_debye', norm2(mbe_result%dipole)*AU_TO_DEBYE)
      end if

      if (mbe_result%has_gradient) then
         call json%add(main_obj, 'gradient_norm', sqrt(sum(mbe_result%gradient**2)))
      end if

      if (mbe_result%has_hessian) then
         call json%add(main_obj, 'hessian_frobenius_norm', sqrt(sum(mbe_result%hessian**2)))
      end if

      ! Write to file
      call logger%info("Writing JSON output to "//trim(output_file))
      open (newunit=iunit, file=trim(output_file), status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         call json%destroy(root)
         return
      end if

      call json%print(root, iunit)
      close (iunit)
      call json%destroy(root)
      call logger%info("JSON output written successfully to "//trim(output_file))

   end subroutine print_detailed_breakdown_json

   subroutine print_unfragmented_json(result)
      !! Write unfragmented calculation results to output JSON file
      !! Uses json-fortran for clean, maintainable JSON output
      type(calculation_result_t), intent(in) :: result

      type(json_core) :: json
      type(json_value), pointer :: root, main_obj, dipole_obj
      integer :: iunit, io_stat
      character(len=256) :: output_file, basename

      output_file = get_output_json_filename()
      basename = get_basename()

      call json%initialize()
      call json%create_object(root, '')
      call json%create_object(main_obj, trim(basename))
      call json%add(root, main_obj)

      if (result%has_energy) call json%add(main_obj, 'total_energy', result%energy%total())

      if (result%has_dipole) then
         call json%create_object(dipole_obj, 'dipole')
         call json%add(main_obj, dipole_obj)
         call json%add(dipole_obj, 'x', result%dipole(1))
         call json%add(dipole_obj, 'y', result%dipole(2))
         call json%add(dipole_obj, 'z', result%dipole(3))
         call json%add(dipole_obj, 'magnitude_debye', norm2(result%dipole)*AU_TO_DEBYE)
      end if

      if (result%has_gradient) call json%add(main_obj, 'gradient_norm', sqrt(sum(result%gradient**2)))
      if (result%has_hessian) call json%add(main_obj, 'hessian_frobenius_norm', sqrt(sum(result%hessian**2)))

      call logger%info("Writing JSON output to "//trim(output_file))
      open (newunit=iunit, file=trim(output_file), status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         call json%destroy(root)
         return
      end if

      call json%print(root, iunit)
      close (iunit)
      call json%destroy(root)
      call logger%info("JSON output written successfully to "//trim(output_file))

   end subroutine print_unfragmented_json

   subroutine print_gmbe_json(n_monomers, monomer_indices, monomer_results, &
                              n_intersections, intersection_results, &
                              intersection_sets, intersection_levels, total_energy)
      !! Write GMBE calculation results to output JSON file
      !! Uses json-fortran for clean, maintainable JSON output
      integer, intent(in) :: n_monomers
      integer, intent(in) :: monomer_indices(:)
      type(calculation_result_t), intent(in) :: monomer_results(:)
      integer, intent(in) :: n_intersections
      type(calculation_result_t), intent(in), optional :: intersection_results(:)
      integer, intent(in), optional :: intersection_sets(:, :)
      integer, intent(in), optional :: intersection_levels(:)
      real(dp), intent(in) :: total_energy

      type(json_core) :: json
      type(json_value), pointer :: root, main_obj, monomers_obj, frags_arr, frag_obj
      type(json_value), pointer :: intersect_obj, levels_arr, level_obj, int_frags_arr, int_frag_obj
      integer :: i, j, k, max_level, iunit, io_stat, level_count, n_indices
      integer, allocatable :: indices(:)
      character(len=256) :: output_file, basename

      output_file = get_output_json_filename()
      basename = get_basename()

      call json%initialize()
      call json%create_object(root, '')
      call json%create_object(main_obj, trim(basename))
      call json%add(root, main_obj)

      call json%add(main_obj, 'total_energy', total_energy)

      ! Monomers section
      call json%create_object(monomers_obj, 'monomers')
      call json%add(main_obj, monomers_obj)
      call json%add(monomers_obj, 'count', n_monomers)

      call json%create_array(frags_arr, 'fragments')
      call json%add(monomers_obj, frags_arr)

      do i = 1, n_monomers
         call json%create_object(frag_obj, '')
         call json%add(frags_arr, frag_obj)
         call json%add(frag_obj, 'index', monomer_indices(i))
         call json%add(frag_obj, 'energy', monomer_results(i)%energy%total())
      end do

      ! Intersections section
      if (n_intersections > 0 .and. present(intersection_results) .and. &
          present(intersection_sets) .and. present(intersection_levels)) then
         max_level = maxval(intersection_levels)

         call json%create_object(intersect_obj, 'intersections')
         call json%add(main_obj, intersect_obj)
         call json%add(intersect_obj, 'total_count', n_intersections)

         call json%create_array(levels_arr, 'levels')
         call json%add(intersect_obj, levels_arr)

         do k = 2, max_level
            level_count = count(intersection_levels == k)

            if (level_count > 0) then
               call json%create_object(level_obj, '')
               call json%add(levels_arr, level_obj)
               call json%add(level_obj, 'level', k)
               call json%add(level_obj, 'count', level_count)

               call json%create_array(int_frags_arr, 'fragments')
               call json%add(level_obj, int_frags_arr)

               do i = 1, n_intersections
                  if (intersection_levels(i) == k) then
                     call json%create_object(int_frag_obj, '')
                     call json%add(int_frags_arr, int_frag_obj)

                     ! Count and extract non-zero indices
                     n_indices = count(intersection_sets(:, i) > 0)
                     allocate (indices(n_indices))
                     n_indices = 0
                     do j = 1, size(intersection_sets, 1)
                        if (intersection_sets(j, i) > 0) then
                           n_indices = n_indices + 1
                           indices(n_indices) = intersection_sets(j, i)
                        end if
                     end do
                     call json%add(int_frag_obj, 'indices', indices)
                     deallocate (indices)

                     call json%add(int_frag_obj, 'energy', intersection_results(i)%energy%total())
                  end if
               end do
            end if
         end do
      end if

      ! Write to file
      call logger%info("Writing GMBE JSON output to "//trim(output_file))
      open (newunit=iunit, file=trim(output_file), status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         call json%destroy(root)
         return
      end if

      call json%print(root, iunit)
      close (iunit)
      call json%destroy(root)
      call logger%info("GMBE JSON output written successfully to "//trim(output_file))

   end subroutine print_gmbe_json

   subroutine print_gmbe_pie_json(pie_atom_sets, pie_coefficients, pie_energies, n_pie_terms, total_energy, &
                                  total_gradient, total_hessian)
      !! Write GMBE PIE calculation results to output JSON file
      !! Uses json-fortran for clean, maintainable JSON output
      integer, intent(in) :: pie_atom_sets(:, :)
      integer, intent(in) :: pie_coefficients(:)
      real(dp), intent(in) :: pie_energies(:)
      integer(int64), intent(in) :: n_pie_terms
      real(dp), intent(in) :: total_energy
      real(dp), intent(in), optional :: total_gradient(:, :)
      real(dp), intent(in), optional :: total_hessian(:, :)

      type(json_core) :: json
      type(json_value), pointer :: root, main_obj, pie_obj, terms_arr, term_obj
      integer :: j, max_atoms, n_atoms, iunit, io_stat
      integer(int64) :: i, n_nonzero_terms
      integer, allocatable :: atom_indices(:)
      character(len=256) :: output_file, basename

      output_file = get_output_json_filename()
      basename = get_basename()

      call json%initialize()
      call json%create_object(root, '')
      call json%create_object(main_obj, trim(basename))
      call json%add(root, main_obj)

      call json%add(main_obj, 'total_energy', total_energy)
      if (present(total_gradient)) call json%add(main_obj, 'gradient_norm', sqrt(sum(total_gradient**2)))
      if (present(total_hessian)) call json%add(main_obj, 'hessian_frobenius_norm', sqrt(sum(total_hessian**2)))

      ! Count non-zero coefficient terms
      n_nonzero_terms = count(pie_coefficients(1:n_pie_terms) /= 0)

      ! PIE terms section
      call json%create_object(pie_obj, 'pie_terms')
      call json%add(main_obj, pie_obj)
      call json%add(pie_obj, 'count', int(n_nonzero_terms))

      call json%create_array(terms_arr, 'terms')
      call json%add(pie_obj, terms_arr)

      max_atoms = size(pie_atom_sets, 1)

      do i = 1_int64, n_pie_terms
         if (pie_coefficients(i) == 0) cycle

         call json%create_object(term_obj, '')
         call json%add(terms_arr, term_obj)

         ! Extract atom list size (atoms until negative sentinel)
         n_atoms = 0
         do while (n_atoms < max_atoms .and. pie_atom_sets(n_atoms + 1, i) >= 0)
            n_atoms = n_atoms + 1
         end do

         allocate (atom_indices(n_atoms))
         atom_indices = pie_atom_sets(1:n_atoms, i)
         call json%add(term_obj, 'atom_indices', atom_indices)
         deallocate (atom_indices)

         call json%add(term_obj, 'coefficient', pie_coefficients(i))
         call json%add(term_obj, 'energy', pie_energies(i))
         call json%add(term_obj, 'weighted_energy', real(pie_coefficients(i), dp)*pie_energies(i))
      end do

      ! Write to file
      call logger%info("Writing GMBE PIE JSON output to "//trim(output_file))
      open (newunit=iunit, file=trim(output_file), status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         call json%destroy(root)
         return
      end if

      call json%print(root, iunit)
      close (iunit)
      call json%destroy(root)
      call logger%info("GMBE PIE JSON output written successfully to "//trim(output_file))

   end subroutine print_gmbe_pie_json

   subroutine write_vibrational_json_impl(total_energy, has_energy, &
                                          dipole, has_dipole, &
                                          gradient_norm, has_gradient, &
                                          hessian_norm, has_hessian, &
                                          frequencies, reduced_masses, force_constants, &
                                          thermo, ir_intensities)
      !! Private implementation for writing vibrational/thermochemistry JSON
      !! Uses json-fortran library for clean, maintainable JSON output
      real(dp), intent(in) :: total_energy
      logical, intent(in) :: has_energy
      real(dp), intent(in), optional :: dipole(:)
      logical, intent(in) :: has_dipole
      real(dp), intent(in), optional :: gradient_norm
      logical, intent(in) :: has_gradient
      real(dp), intent(in), optional :: hessian_norm
      logical, intent(in) :: has_hessian
      real(dp), intent(in) :: frequencies(:)
      real(dp), intent(in) :: reduced_masses(:)
      real(dp), intent(in) :: force_constants(:)
      type(thermochemistry_result_t), intent(in) :: thermo
      real(dp), intent(in), optional :: ir_intensities(:)

      type(json_core) :: json
      type(json_value), pointer :: root, main_obj, dipole_obj, vib_obj, thermo_obj
      type(json_value), pointer :: moi_obj, rot_obj, pf_obj, contrib_obj, table_obj
      type(json_value), pointer :: trans_obj, rot_contrib, vib_contrib, elec_obj
      type(json_value), pointer :: vib_row, rot_row, int_row, tr_row, tot_row
      type(json_value), pointer :: thermal_obj, total_e_obj
      integer :: io_stat, iunit
      character(len=256) :: output_file, basename
      real(dp) :: pV_cal, H_vib_cal, H_rot_cal, H_trans_cal, H_total_cal
      real(dp) :: Cv_total, S_total, S_total_J, H_int_cal, Cv_int, S_int, S_int_J, Cp_trans

      output_file = get_output_json_filename()
      basename = get_basename()

      call json%initialize()

      ! Create root object
      call json%create_object(root, '')

      ! Create main object with basename as key
      call json%create_object(main_obj, trim(basename))
      call json%add(root, main_obj)

      ! Total energy
      if (has_energy) call json%add(main_obj, 'total_energy', total_energy)

      ! Dipole
      if (has_dipole .and. present(dipole)) then
         call json%create_object(dipole_obj, 'dipole')
         call json%add(main_obj, dipole_obj)
         call json%add(dipole_obj, 'x', dipole(1))
         call json%add(dipole_obj, 'y', dipole(2))
         call json%add(dipole_obj, 'z', dipole(3))
         call json%add(dipole_obj, 'magnitude_debye', norm2(dipole)*AU_TO_DEBYE)
      end if

      ! Gradient and Hessian norms
      if (has_gradient .and. present(gradient_norm)) call json%add(main_obj, 'gradient_norm', gradient_norm)
      if (has_hessian .and. present(hessian_norm)) call json%add(main_obj, 'hessian_frobenius_norm', hessian_norm)

      ! Vibrational analysis section
      call json%create_object(vib_obj, 'vibrational_analysis')
      call json%add(main_obj, vib_obj)
      call json%add(vib_obj, 'n_modes', size(frequencies))
      call json%add(vib_obj, 'frequencies_cm1', frequencies)
      call json%add(vib_obj, 'reduced_masses_amu', reduced_masses)
      call json%add(vib_obj, 'force_constants_mdyne_ang', force_constants)
      if (present(ir_intensities)) call json%add(vib_obj, 'ir_intensities_km_mol', ir_intensities)

      ! Thermochemistry section
      call json%create_object(thermo_obj, 'thermochemistry')
      call json%add(main_obj, thermo_obj)

      ! Conditions
      call json%add(thermo_obj, 'temperature_K', thermo%temperature)
      call json%add(thermo_obj, 'pressure_atm', thermo%pressure)
      call json%add(thermo_obj, 'molecular_mass_amu', thermo%total_mass)
      call json%add(thermo_obj, 'symmetry_number', thermo%symmetry_number)
      call json%add(thermo_obj, 'spin_multiplicity', thermo%spin_multiplicity)
      call json%add(thermo_obj, 'is_linear', thermo%is_linear)
      call json%add(thermo_obj, 'n_real_frequencies', thermo%n_real_freqs)
      call json%add(thermo_obj, 'n_imaginary_frequencies', thermo%n_imag_freqs)

      ! Moments of inertia
      call json%create_object(moi_obj, 'moments_of_inertia_amu_ang2')
      call json%add(thermo_obj, moi_obj)
      call json%add(moi_obj, 'Ia', thermo%moments(1))
      call json%add(moi_obj, 'Ib', thermo%moments(2))
      call json%add(moi_obj, 'Ic', thermo%moments(3))

      ! Rotational constants
      call json%create_object(rot_obj, 'rotational_constants_GHz')
      call json%add(thermo_obj, rot_obj)
      call json%add(rot_obj, 'A', thermo%rot_const(1))
      call json%add(rot_obj, 'B', thermo%rot_const(2))
      call json%add(rot_obj, 'C', thermo%rot_const(3))

      ! Partition functions
      call json%create_object(pf_obj, 'partition_functions')
      call json%add(thermo_obj, pf_obj)
      call json%add(pf_obj, 'translational', thermo%q_trans)
      call json%add(pf_obj, 'rotational', thermo%q_rot)
      call json%add(pf_obj, 'vibrational', thermo%q_vib)

      ! Thermodynamic contributions
      call json%create_object(contrib_obj, 'contributions')
      call json%add(thermo_obj, contrib_obj)

      call json%create_object(trans_obj, 'translational')
      call json%add(contrib_obj, trans_obj)
      call json%add(trans_obj, 'energy_hartree', thermo%E_trans)
      call json%add(trans_obj, 'entropy_cal_mol_K', thermo%S_trans)
      call json%add(trans_obj, 'Cv_cal_mol_K', thermo%Cv_trans)

      call json%create_object(rot_contrib, 'rotational')
      call json%add(contrib_obj, rot_contrib)
      call json%add(rot_contrib, 'energy_hartree', thermo%E_rot)
      call json%add(rot_contrib, 'entropy_cal_mol_K', thermo%S_rot)
      call json%add(rot_contrib, 'Cv_cal_mol_K', thermo%Cv_rot)

      call json%create_object(vib_contrib, 'vibrational')
      call json%add(contrib_obj, vib_contrib)
      call json%add(vib_contrib, 'energy_hartree', thermo%E_vib)
      call json%add(vib_contrib, 'entropy_cal_mol_K', thermo%S_vib)
      call json%add(vib_contrib, 'Cv_cal_mol_K', thermo%Cv_vib)

      call json%create_object(elec_obj, 'electronic')
      call json%add(contrib_obj, elec_obj)
      call json%add(elec_obj, 'energy_hartree', thermo%E_elec)
      call json%add(elec_obj, 'entropy_cal_mol_K', thermo%S_elec)

      ! Contribution table
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
      S_int_J = S_int*CAL_TO_J
      S_total = thermo%S_trans + thermo%S_rot + thermo%S_vib + thermo%S_elec
      S_total_J = S_total*CAL_TO_J

      call json%create_object(table_obj, 'contribution_table')
      call json%add(thermo_obj, table_obj)

      call json%create_object(vib_row, 'VIB')
      call json%add(table_obj, vib_row)
      call json%add(vib_row, 'H_cal_mol', H_vib_cal)
      call json%add(vib_row, 'Cp_cal_mol_K', thermo%Cv_vib)
      call json%add(vib_row, 'S_cal_mol_K', thermo%S_vib)
      call json%add(vib_row, 'S_J_mol_K', thermo%S_vib*CAL_TO_J)

      call json%create_object(rot_row, 'ROT')
      call json%add(table_obj, rot_row)
      call json%add(rot_row, 'H_cal_mol', H_rot_cal)
      call json%add(rot_row, 'Cp_cal_mol_K', thermo%Cv_rot)
      call json%add(rot_row, 'S_cal_mol_K', thermo%S_rot)
      call json%add(rot_row, 'S_J_mol_K', thermo%S_rot*CAL_TO_J)

      call json%create_object(int_row, 'INT')
      call json%add(table_obj, int_row)
      call json%add(int_row, 'H_cal_mol', H_int_cal)
      call json%add(int_row, 'Cp_cal_mol_K', Cv_int)
      call json%add(int_row, 'S_cal_mol_K', S_int)
      call json%add(int_row, 'S_J_mol_K', S_int_J)

      call json%create_object(tr_row, 'TR')
      call json%add(table_obj, tr_row)
      call json%add(tr_row, 'H_cal_mol', H_trans_cal)
      call json%add(tr_row, 'Cp_cal_mol_K', Cp_trans)
      call json%add(tr_row, 'S_cal_mol_K', thermo%S_trans)
      call json%add(tr_row, 'S_J_mol_K', thermo%S_trans*CAL_TO_J)

      call json%create_object(tot_row, 'TOT')
      call json%add(table_obj, tot_row)
      call json%add(tot_row, 'H_cal_mol', H_total_cal)
      call json%add(tot_row, 'Cp_cal_mol_K', Cv_total)
      call json%add(tot_row, 'S_cal_mol_K', S_total)
      call json%add(tot_row, 'S_J_mol_K', S_total_J)

      ! Zero-point energy
      call json%add(thermo_obj, 'zero_point_energy_hartree', thermo%zpe_hartree)
      call json%add(thermo_obj, 'zero_point_energy_kcal_mol', thermo%zpe_kcalmol)

      ! Thermal corrections
      call json%create_object(thermal_obj, 'thermal_corrections_hartree')
      call json%add(thermo_obj, thermal_obj)
      call json%add(thermal_obj, 'to_energy', thermo%thermal_correction_energy)
      call json%add(thermal_obj, 'to_enthalpy', thermo%thermal_correction_enthalpy)
      call json%add(thermal_obj, 'to_gibbs', thermo%thermal_correction_gibbs)

      ! Total energies
      call json%create_object(total_e_obj, 'total_energies_hartree')
      call json%add(thermo_obj, total_e_obj)
      call json%add(total_e_obj, 'electronic', total_energy)
      call json%add(total_e_obj, 'electronic_plus_zpe', total_energy + thermo%zpe_hartree)
      call json%add(total_e_obj, 'electronic_plus_thermal_E', total_energy + thermo%thermal_correction_energy)
      call json%add(total_e_obj, 'electronic_plus_thermal_H', total_energy + thermo%thermal_correction_enthalpy)
      call json%add(total_e_obj, 'electronic_plus_thermal_G', total_energy + thermo%thermal_correction_gibbs)

      ! Write to file
      call logger%info("Writing vibrational/thermochemistry JSON to "//trim(output_file))
      open (newunit=iunit, file=trim(output_file), status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         call json%destroy(root)
         return
      end if

      call json%print(root, iunit)
      close (iunit)

      call json%destroy(root)
      call logger%info("Vibrational/thermochemistry JSON written successfully to "//trim(output_file))

   end subroutine write_vibrational_json_impl

   subroutine print_vibrational_json(result, frequencies, reduced_masses, force_constants, &
                                     thermo, ir_intensities)
      !! Write vibrational analysis and thermochemistry results to JSON file
      !! Wrapper that extracts values from calculation_result_t and calls the implementation
      type(calculation_result_t), intent(in) :: result
      real(dp), intent(in) :: frequencies(:)
      real(dp), intent(in) :: reduced_masses(:)
      real(dp), intent(in) :: force_constants(:)
      type(thermochemistry_result_t), intent(in) :: thermo
      real(dp), intent(in), optional :: ir_intensities(:)

      real(dp) :: energy_val, grad_norm, hess_norm

      ! Extract energy
      if (result%has_energy) then
         energy_val = result%energy%total()
      else
         energy_val = 0.0_dp
      end if

      ! Compute norms if available
      if (result%has_gradient) then
         grad_norm = sqrt(sum(result%gradient**2))
      end if
      if (result%has_hessian) then
         hess_norm = sqrt(sum(result%hessian**2))
      end if

      ! Call implementation with extracted values
      if (present(ir_intensities)) then
         call write_vibrational_json_impl(energy_val, result%has_energy, &
                                          result%dipole, result%has_dipole, &
                                          grad_norm, result%has_gradient, &
                                          hess_norm, result%has_hessian, &
                                          frequencies, reduced_masses, force_constants, &
                                          thermo, ir_intensities)
      else
         call write_vibrational_json_impl(energy_val, result%has_energy, &
                                          result%dipole, result%has_dipole, &
                                          grad_norm, result%has_gradient, &
                                          hess_norm, result%has_hessian, &
                                          frequencies, reduced_masses, force_constants, &
                                          thermo)
      end if

   end subroutine print_vibrational_json

   subroutine print_vibrational_json_mbe(mbe_result, frequencies, reduced_masses, force_constants, &
                                         thermo, ir_intensities)
      !! Write vibrational analysis and thermochemistry results to JSON file for MBE calculations
      !! Wrapper that extracts values from mbe_result_t and calls the implementation
      type(mbe_result_t), intent(in) :: mbe_result
      real(dp), intent(in) :: frequencies(:)
      real(dp), intent(in) :: reduced_masses(:)
      real(dp), intent(in) :: force_constants(:)
      type(thermochemistry_result_t), intent(in) :: thermo
      real(dp), intent(in), optional :: ir_intensities(:)

      real(dp) :: grad_norm, hess_norm

      ! Compute norms if available
      if (mbe_result%has_gradient) then
         grad_norm = sqrt(sum(mbe_result%gradient**2))
      end if
      if (mbe_result%has_hessian) then
         hess_norm = sqrt(sum(mbe_result%hessian**2))
      end if

      ! Call implementation with extracted values
      ! mbe_result_t stores energy directly as total_energy
      if (present(ir_intensities)) then
         call write_vibrational_json_impl(mbe_result%total_energy, mbe_result%has_energy, &
                                          mbe_result%dipole, mbe_result%has_dipole, &
                                          grad_norm, mbe_result%has_gradient, &
                                          hess_norm, mbe_result%has_hessian, &
                                          frequencies, reduced_masses, force_constants, &
                                          thermo, ir_intensities)
      else
         call write_vibrational_json_impl(mbe_result%total_energy, mbe_result%has_energy, &
                                          mbe_result%dipole, mbe_result%has_dipole, &
                                          grad_norm, mbe_result%has_gradient, &
                                          hess_norm, mbe_result%has_hessian, &
                                          frequencies, reduced_masses, force_constants, &
                                          thermo)
      end if

   end subroutine print_vibrational_json_mbe

end module mqc_mbe_io
