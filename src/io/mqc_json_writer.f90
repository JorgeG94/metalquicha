!! Centralized JSON output writer
!! Single entry point for all JSON output in metalquicha
module mqc_json_writer
   use pic_types, only: int64, dp
   use pic_logger, only: logger => global_logger
   use mqc_json_output_types, only: json_output_data_t, &
                                    OUTPUT_MODE_UNFRAGMENTED, OUTPUT_MODE_MBE, OUTPUT_MODE_GMBE_PIE
   use mqc_io_helpers, only: get_output_json_filename, get_basename
   use mqc_physical_constants, only: HARTREE_TO_CALMOL, R_CALMOLK, AU_TO_DEBYE, CAL_TO_J
   use mqc_program_limits, only: JSON_REAL_FORMAT
   use mqc_mbe_io, only: get_frag_level_name
   use json_module, only: json_core, json_value
   implicit none
   private

   public :: write_json_output  !! Single entry point for all JSON output

contains

   subroutine write_json_output(output_data)
      !! THE single entry point for all JSON output
      !!
      !! Dispatches to the appropriate JSON writer based on output_mode.
      !! This is the ONLY place in the codebase where JSON files are written.
      type(json_output_data_t), intent(in) :: output_data

      select case (output_data%output_mode)
      case (OUTPUT_MODE_UNFRAGMENTED)
         if (output_data%has_vibrational) then
            call write_vibrational_json_impl(output_data)
         else
            call write_unfragmented_json_impl(output_data)
         end if

      case (OUTPUT_MODE_MBE)
         if (output_data%has_vibrational) then
            call write_vibrational_json_impl(output_data)
         else
            call write_mbe_breakdown_json_impl(output_data)
         end if

      case (OUTPUT_MODE_GMBE_PIE)
         if (output_data%has_vibrational) then
            call write_vibrational_json_impl(output_data)
         else
            call write_gmbe_pie_json_impl(output_data)
         end if

      case default
         call logger%error("Unknown output mode in write_json_output")
      end select

   end subroutine write_json_output

   subroutine write_unfragmented_json_impl(data)
      !! Write unfragmented calculation results to output JSON file
      type(json_output_data_t), intent(in) :: data

      type(json_core) :: json
      type(json_value), pointer :: root, main_obj, dipole_obj
      integer :: iunit, io_stat
      character(len=256) :: output_file, basename

      output_file = get_output_json_filename()
      basename = get_basename()

      call json%initialize(real_format=JSON_REAL_FORMAT)
      call json%create_object(root, '')
      call json%create_object(main_obj, trim(basename))
      call json%add(root, main_obj)

      if (data%has_energy) call json%add(main_obj, 'total_energy', data%total_energy)

      if (data%has_dipole .and. allocated(data%dipole)) then
         call json%create_object(dipole_obj, 'dipole')
         call json%add(main_obj, dipole_obj)
         call json%add(dipole_obj, 'x', data%dipole(1))
         call json%add(dipole_obj, 'y', data%dipole(2))
         call json%add(dipole_obj, 'z', data%dipole(3))
         call json%add(dipole_obj, 'magnitude_debye', norm2(data%dipole)*AU_TO_DEBYE)
      end if

      if (data%has_gradient .and. allocated(data%gradient)) then
         call json%add(main_obj, 'gradient_norm', sqrt(sum(data%gradient**2)))
      end if
      if (data%has_hessian .and. allocated(data%hessian)) then
         call json%add(main_obj, 'hessian_frobenius_norm', sqrt(sum(data%hessian**2)))
      end if

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

   end subroutine write_unfragmented_json_impl

   subroutine write_mbe_breakdown_json_impl(data)
      !! Write detailed MBE energy breakdown to JSON file
      type(json_output_data_t), intent(in) :: data

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

      if (data%max_level > 10) then
         call logger%warning("Fragment levels exceed decamers (10-mers). JSON will use generic N-mers notation.")
      end if

      call json%initialize(real_format=JSON_REAL_FORMAT)
      call json%create_object(root, '')
      call json%create_object(main_obj, trim(basename))
      call json%add(root, main_obj)

      call json%add(main_obj, 'total_energy', data%total_energy)

      ! Build levels array
      call json%create_array(levels_arr, 'levels')
      call json%add(main_obj, levels_arr)

      do frag_level = 1, data%max_level
         count_by_level = 0_int64
         do i = 1_int64, data%fragment_count
            fragment_size = count(data%polymers(i, :) > 0)
            if (fragment_size == frag_level) count_by_level = count_by_level + 1_int64
         end do

         if (count_by_level > 0_int64) then
            call json%create_object(level_obj, '')
            call json%add(levels_arr, level_obj)

            level_name = get_frag_level_name(frag_level)
            call json%add(level_obj, 'frag_level', frag_level)
            call json%add(level_obj, 'name', trim(level_name))
            call json%add(level_obj, 'count', int(count_by_level))
            if (allocated(data%sum_by_level)) then
               call json%add(level_obj, 'total_energy', data%sum_by_level(frag_level))
            end if

            call json%create_array(frags_arr, 'fragments')
            call json%add(level_obj, frags_arr)

            do i = 1_int64, data%fragment_count
               fragment_size = count(data%polymers(i, :) > 0)
               if (fragment_size == frag_level) then
                  call json%create_object(frag_obj, '')
                  call json%add(frags_arr, frag_obj)

                  allocate (indices(fragment_size))
                  indices = data%polymers(i, 1:fragment_size)
                  call json%add(frag_obj, 'indices', indices)
                  deallocate (indices)

                  if (allocated(data%fragment_energies)) then
                     call json%add(frag_obj, 'energy', data%fragment_energies(i))
                  end if

                  if (allocated(data%fragment_distances)) then
                     call json%add(frag_obj, 'distance', data%fragment_distances(i))
                  end if

                  if (frag_level > 1 .and. allocated(data%delta_energies)) then
                     call json%add(frag_obj, 'delta_energy', data%delta_energies(i))
                  end if
               end if
            end do
         end if
      end do

      ! Add dipole if computed
      if (data%has_dipole .and. allocated(data%dipole)) then
         call json%create_object(dipole_obj, 'dipole')
         call json%add(main_obj, dipole_obj)
         call json%add(dipole_obj, 'x', data%dipole(1))
         call json%add(dipole_obj, 'y', data%dipole(2))
         call json%add(dipole_obj, 'z', data%dipole(3))
         call json%add(dipole_obj, 'magnitude_debye', norm2(data%dipole)*AU_TO_DEBYE)
      end if

      if (data%has_gradient .and. allocated(data%gradient)) then
         call json%add(main_obj, 'gradient_norm', sqrt(sum(data%gradient**2)))
      end if

      if (data%has_hessian .and. allocated(data%hessian)) then
         call json%add(main_obj, 'hessian_frobenius_norm', sqrt(sum(data%hessian**2)))
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

   end subroutine write_mbe_breakdown_json_impl

   subroutine write_gmbe_pie_json_impl(data)
      !! Write GMBE PIE calculation results to output JSON file
      type(json_output_data_t), intent(in) :: data

      type(json_core) :: json
      type(json_value), pointer :: root, main_obj, pie_obj, terms_arr, term_obj
      integer :: j, max_atoms, n_atoms, iunit, io_stat
      integer(int64) :: i, n_nonzero_terms
      integer, allocatable :: atom_indices(:)
      character(len=256) :: output_file, basename

      output_file = get_output_json_filename()
      basename = get_basename()

      call json%initialize(real_format=JSON_REAL_FORMAT)
      call json%create_object(root, '')
      call json%create_object(main_obj, trim(basename))
      call json%add(root, main_obj)

      call json%add(main_obj, 'total_energy', data%total_energy)

      if (data%has_gradient .and. allocated(data%gradient)) then
         call json%add(main_obj, 'gradient_norm', sqrt(sum(data%gradient**2)))
      end if
      if (data%has_hessian .and. allocated(data%hessian)) then
         call json%add(main_obj, 'hessian_frobenius_norm', sqrt(sum(data%hessian**2)))
      end if

      ! Count non-zero coefficient terms
      if (allocated(data%pie_coefficients)) then
         n_nonzero_terms = count(data%pie_coefficients(1:data%n_pie_terms) /= 0)
      else
         n_nonzero_terms = 0
      end if

      ! PIE terms section
      call json%create_object(pie_obj, 'pie_terms')
      call json%add(main_obj, pie_obj)
      call json%add(pie_obj, 'count', int(n_nonzero_terms))

      call json%create_array(terms_arr, 'terms')
      call json%add(pie_obj, terms_arr)

      if (allocated(data%pie_atom_sets) .and. allocated(data%pie_coefficients) .and. &
          allocated(data%pie_energies)) then
         max_atoms = size(data%pie_atom_sets, 1)

         do i = 1_int64, data%n_pie_terms
            if (data%pie_coefficients(i) == 0) cycle

            call json%create_object(term_obj, '')
            call json%add(terms_arr, term_obj)

            ! Extract atom list size (atoms until negative sentinel)
            n_atoms = 0
            do while (n_atoms < max_atoms .and. data%pie_atom_sets(n_atoms + 1, i) >= 0)
               n_atoms = n_atoms + 1
            end do

            allocate (atom_indices(n_atoms))
            atom_indices = data%pie_atom_sets(1:n_atoms, i)
            call json%add(term_obj, 'atom_indices', atom_indices)
            deallocate (atom_indices)

            call json%add(term_obj, 'coefficient', data%pie_coefficients(i))
            call json%add(term_obj, 'energy', data%pie_energies(i))
            call json%add(term_obj, 'weighted_energy', real(data%pie_coefficients(i), dp)*data%pie_energies(i))
         end do
      end if

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

   end subroutine write_gmbe_pie_json_impl

   subroutine write_vibrational_json_impl(data)
      !! Write vibrational analysis and thermochemistry results to JSON file
      type(json_output_data_t), intent(in) :: data

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
      real(dp) :: grad_norm, hess_norm

      output_file = get_output_json_filename()
      basename = get_basename()

      call json%initialize(real_format=JSON_REAL_FORMAT)

      ! Create root object
      call json%create_object(root, '')

      ! Create main object with basename as key
      call json%create_object(main_obj, trim(basename))
      call json%add(root, main_obj)

      ! Total energy
      if (data%has_energy) call json%add(main_obj, 'total_energy', data%total_energy)

      ! Dipole
      if (data%has_dipole .and. allocated(data%dipole)) then
         call json%create_object(dipole_obj, 'dipole')
         call json%add(main_obj, dipole_obj)
         call json%add(dipole_obj, 'x', data%dipole(1))
         call json%add(dipole_obj, 'y', data%dipole(2))
         call json%add(dipole_obj, 'z', data%dipole(3))
         call json%add(dipole_obj, 'magnitude_debye', norm2(data%dipole)*AU_TO_DEBYE)
      end if

      ! Gradient and Hessian norms
      if (data%has_gradient .and. allocated(data%gradient)) then
         grad_norm = sqrt(sum(data%gradient**2))
         call json%add(main_obj, 'gradient_norm', grad_norm)
      end if
      if (data%has_hessian .and. allocated(data%hessian)) then
         hess_norm = sqrt(sum(data%hessian**2))
         call json%add(main_obj, 'hessian_frobenius_norm', hess_norm)
      end if

      ! Vibrational analysis section
      call json%create_object(vib_obj, 'vibrational_analysis')
      call json%add(main_obj, vib_obj)
      if (allocated(data%frequencies)) then
         call json%add(vib_obj, 'n_modes', size(data%frequencies))
         call json%add(vib_obj, 'frequencies_cm1', data%frequencies)
      end if
      if (allocated(data%reduced_masses)) then
         call json%add(vib_obj, 'reduced_masses_amu', data%reduced_masses)
      end if
      if (allocated(data%force_constants)) then
         call json%add(vib_obj, 'force_constants_mdyne_ang', data%force_constants)
      end if
      if (data%has_ir_intensities .and. allocated(data%ir_intensities)) then
         call json%add(vib_obj, 'ir_intensities_km_mol', data%ir_intensities)
      end if

      ! Thermochemistry section
      call json%create_object(thermo_obj, 'thermochemistry')
      call json%add(main_obj, thermo_obj)

      ! Conditions
      call json%add(thermo_obj, 'temperature_K', data%thermo%temperature)
      call json%add(thermo_obj, 'pressure_atm', data%thermo%pressure)
      call json%add(thermo_obj, 'molecular_mass_amu', data%thermo%total_mass)
      call json%add(thermo_obj, 'symmetry_number', data%thermo%symmetry_number)
      call json%add(thermo_obj, 'spin_multiplicity', data%thermo%spin_multiplicity)
      call json%add(thermo_obj, 'is_linear', data%thermo%is_linear)
      call json%add(thermo_obj, 'n_real_frequencies', data%thermo%n_real_freqs)
      call json%add(thermo_obj, 'n_imaginary_frequencies', data%thermo%n_imag_freqs)

      ! Moments of inertia
      call json%create_object(moi_obj, 'moments_of_inertia_amu_ang2')
      call json%add(thermo_obj, moi_obj)
      call json%add(moi_obj, 'Ia', data%thermo%moments(1))
      call json%add(moi_obj, 'Ib', data%thermo%moments(2))
      call json%add(moi_obj, 'Ic', data%thermo%moments(3))

      ! Rotational constants
      call json%create_object(rot_obj, 'rotational_constants_GHz')
      call json%add(thermo_obj, rot_obj)
      call json%add(rot_obj, 'A', data%thermo%rot_const(1))
      call json%add(rot_obj, 'B', data%thermo%rot_const(2))
      call json%add(rot_obj, 'C', data%thermo%rot_const(3))

      ! Partition functions
      call json%create_object(pf_obj, 'partition_functions')
      call json%add(thermo_obj, pf_obj)
      call json%add(pf_obj, 'translational', data%thermo%q_trans)
      call json%add(pf_obj, 'rotational', data%thermo%q_rot)
      call json%add(pf_obj, 'vibrational', data%thermo%q_vib)

      ! Thermodynamic contributions
      call json%create_object(contrib_obj, 'contributions')
      call json%add(thermo_obj, contrib_obj)

      call json%create_object(trans_obj, 'translational')
      call json%add(contrib_obj, trans_obj)
      call json%add(trans_obj, 'energy_hartree', data%thermo%E_trans)
      call json%add(trans_obj, 'entropy_cal_mol_K', data%thermo%S_trans)
      call json%add(trans_obj, 'Cv_cal_mol_K', data%thermo%Cv_trans)

      call json%create_object(rot_contrib, 'rotational')
      call json%add(contrib_obj, rot_contrib)
      call json%add(rot_contrib, 'energy_hartree', data%thermo%E_rot)
      call json%add(rot_contrib, 'entropy_cal_mol_K', data%thermo%S_rot)
      call json%add(rot_contrib, 'Cv_cal_mol_K', data%thermo%Cv_rot)

      call json%create_object(vib_contrib, 'vibrational')
      call json%add(contrib_obj, vib_contrib)
      call json%add(vib_contrib, 'energy_hartree', data%thermo%E_vib)
      call json%add(vib_contrib, 'entropy_cal_mol_K', data%thermo%S_vib)
      call json%add(vib_contrib, 'Cv_cal_mol_K', data%thermo%Cv_vib)

      call json%create_object(elec_obj, 'electronic')
      call json%add(contrib_obj, elec_obj)
      call json%add(elec_obj, 'energy_hartree', data%thermo%E_elec)
      call json%add(elec_obj, 'entropy_cal_mol_K', data%thermo%S_elec)

      ! Contribution table
      pV_cal = R_CALMOLK*data%thermo%temperature
      H_vib_cal = data%thermo%E_vib*HARTREE_TO_CALMOL
      H_rot_cal = data%thermo%E_rot*HARTREE_TO_CALMOL
      H_trans_cal = data%thermo%E_trans*HARTREE_TO_CALMOL + pV_cal
      H_total_cal = H_vib_cal + H_rot_cal + H_trans_cal
      H_int_cal = H_vib_cal + H_rot_cal
      Cp_trans = data%thermo%Cv_trans + R_CALMOLK
      Cv_int = data%thermo%Cv_vib + data%thermo%Cv_rot
      Cv_total = Cp_trans + data%thermo%Cv_rot + data%thermo%Cv_vib
      S_int = data%thermo%S_vib + data%thermo%S_rot
      S_int_J = S_int*CAL_TO_J
      S_total = data%thermo%S_trans + data%thermo%S_rot + data%thermo%S_vib + data%thermo%S_elec
      S_total_J = S_total*CAL_TO_J

      call json%create_object(table_obj, 'contribution_table')
      call json%add(thermo_obj, table_obj)

      call json%create_object(vib_row, 'VIB')
      call json%add(table_obj, vib_row)
      call json%add(vib_row, 'H_cal_mol', H_vib_cal)
      call json%add(vib_row, 'Cp_cal_mol_K', data%thermo%Cv_vib)
      call json%add(vib_row, 'S_cal_mol_K', data%thermo%S_vib)
      call json%add(vib_row, 'S_J_mol_K', data%thermo%S_vib*CAL_TO_J)

      call json%create_object(rot_row, 'ROT')
      call json%add(table_obj, rot_row)
      call json%add(rot_row, 'H_cal_mol', H_rot_cal)
      call json%add(rot_row, 'Cp_cal_mol_K', data%thermo%Cv_rot)
      call json%add(rot_row, 'S_cal_mol_K', data%thermo%S_rot)
      call json%add(rot_row, 'S_J_mol_K', data%thermo%S_rot*CAL_TO_J)

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
      call json%add(tr_row, 'S_cal_mol_K', data%thermo%S_trans)
      call json%add(tr_row, 'S_J_mol_K', data%thermo%S_trans*CAL_TO_J)

      call json%create_object(tot_row, 'TOT')
      call json%add(table_obj, tot_row)
      call json%add(tot_row, 'H_cal_mol', H_total_cal)
      call json%add(tot_row, 'Cp_cal_mol_K', Cv_total)
      call json%add(tot_row, 'S_cal_mol_K', S_total)
      call json%add(tot_row, 'S_J_mol_K', S_total_J)

      ! Zero-point energy
      call json%add(thermo_obj, 'zero_point_energy_hartree', data%thermo%zpe_hartree)
      call json%add(thermo_obj, 'zero_point_energy_kcal_mol', data%thermo%zpe_kcalmol)

      ! Thermal corrections
      call json%create_object(thermal_obj, 'thermal_corrections_hartree')
      call json%add(thermo_obj, thermal_obj)
      call json%add(thermal_obj, 'to_energy', data%thermo%thermal_correction_energy)
      call json%add(thermal_obj, 'to_enthalpy', data%thermo%thermal_correction_enthalpy)
      call json%add(thermal_obj, 'to_gibbs', data%thermo%thermal_correction_gibbs)

      ! Total energies
      call json%create_object(total_e_obj, 'total_energies_hartree')
      call json%add(thermo_obj, total_e_obj)
      call json%add(total_e_obj, 'electronic', data%total_energy)
      call json%add(total_e_obj, 'electronic_plus_zpe', data%total_energy + data%thermo%zpe_hartree)
      call json%add(total_e_obj, 'electronic_plus_thermal_E', data%total_energy + data%thermo%thermal_correction_energy)
      call json%add(total_e_obj, 'electronic_plus_thermal_H', data%total_energy + data%thermo%thermal_correction_enthalpy)
      call json%add(total_e_obj, 'electronic_plus_thermal_G', data%total_energy + data%thermo%thermal_correction_gibbs)

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

end module mqc_json_writer
