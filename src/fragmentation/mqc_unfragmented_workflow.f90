submodule(mqc_mbe_fragment_distribution_scheme) mqc_unfragmented_workflow
   implicit none
contains
   module subroutine unfragmented_calculation(sys_geom, method, calc_type, bonds, result_out, &
                                              temperature, pressure, json_data)
      !! Run unfragmented calculation on the entire system (nlevel=0)
      !! This is a simple single-process calculation without MPI distribution
      !! If result_out is present, returns result instead of writing JSON and destroying it
      !! If json_data is present, populates it for centralized JSON output
      use mqc_error, only: error_t
      use mqc_vibrational_analysis, only: compute_vibrational_frequencies, &
                                          compute_vibrational_analysis, print_vibrational_analysis
      use mqc_thermochemistry, only: thermochemistry_result_t, compute_thermochemistry
      use mqc_json_output_types, only: json_output_data_t, OUTPUT_MODE_UNFRAGMENTED
      type(system_geometry_t), intent(in), optional :: sys_geom
      integer(int32), intent(in) :: method
      integer(int32), intent(in) :: calc_type
      type(bond_t), intent(in), optional :: bonds(:)
      type(calculation_result_t), intent(out), optional :: result_out
      real(dp), intent(in), optional :: temperature  !! Temperature for thermochemistry (K)
      real(dp), intent(in), optional :: pressure     !! Pressure for thermochemistry (atm)
      type(json_output_data_t), intent(out), optional :: json_data

      type(calculation_result_t) :: result
      integer :: total_atoms
      type(physical_fragment_t) :: full_system
      type(error_t) :: error
      integer :: i

      if (.not. present(sys_geom)) then
         call logger%error("sys_geom required for unfragmented calculation")
         error stop "Missing geometry in unfragmented_calculation"
      end if

      total_atoms = sys_geom%total_atoms

      call logger%info("============================================")
      call logger%info("Running unfragmented calculation")
      call logger%info("  Total atoms: "//to_char(total_atoms))
      call logger%info("============================================")

      ! Build the full system as a single fragment
      ! For overlapping fragments, we use the full system directly (not concatenating fragments)
      full_system%n_atoms = total_atoms
      full_system%n_caps = 0
      allocate (full_system%element_numbers(total_atoms))
      allocate (full_system%coordinates(3, total_atoms))

      ! Copy all atoms from system geometry
      full_system%element_numbers = sys_geom%element_numbers
      full_system%coordinates = sys_geom%coordinates

      ! Set charge and multiplicity from system
      full_system%charge = sys_geom%charge
      full_system%multiplicity = sys_geom%multiplicity
      call full_system%compute_nelec()

      ! Validate geometry (check for spatially overlapping atoms)
      call check_duplicate_atoms(full_system, error)
      if (error%has_error()) then
         call logger%error(error%get_full_trace())
         error stop "Overlapping atoms in unfragmented system"
      end if

      ! Process the full system
      call do_fragment_work(0_int64, result, method, phys_frag=full_system, calc_type=calc_type)

      ! Check for calculation errors
      if (result%has_error) then
         call logger%error("Unfragmented calculation failed: "//result%error%get_message())
         if (present(result_out)) then
            result_out = result
            return
         else
            error stop "Unfragmented calculation failed"
         end if
      end if

      call logger%info("============================================")
      call logger%info("Unfragmented calculation completed")
      block
         character(len=2048) :: result_line  ! Large buffer for Hessian matrix rows
         integer :: current_log_level, iatom, i, j
         real(dp) :: hess_norm

         write (result_line, '(a,f25.15)') "  Final energy: ", result%energy%total()
         call logger%info(trim(result_line))

         if (result%has_dipole) then
            write (result_line, '(a,3f15.8)') "  Dipole (e*Bohr): ", result%dipole
            call logger%info(trim(result_line))
            write (result_line, '(a,f15.8)') "  Dipole magnitude (Debye): ", norm2(result%dipole)*2.541746_dp
            call logger%info(trim(result_line))
         end if

         if (result%has_gradient) then
            write (result_line, '(a,f25.15)') "  Gradient norm: ", sqrt(sum(result%gradient**2))
            call logger%info(trim(result_line))

            ! Print full gradient if verbose and system is small
            call logger%configuration(level=current_log_level)
            if (current_log_level >= verbose_level .and. total_atoms < 100) then
               call logger%info(" ")
               call logger%info("Gradient (Hartree/Bohr):")
               do iatom = 1, total_atoms
                  write (result_line, '(a,i5,a,3f20.12)') "  Atom ", iatom, ": ", &
                     result%gradient(1, iatom), result%gradient(2, iatom), result%gradient(3, iatom)
                  call logger%info(trim(result_line))
               end do
               call logger%info(" ")
            end if
         end if

         if (result%has_hessian) then
            ! Compute Frobenius norm of Hessian
            hess_norm = sqrt(sum(result%hessian**2))
            write (result_line, '(a,f25.15)') "  Hessian Frobenius norm: ", hess_norm
            call logger%info(trim(result_line))

            ! Print full Hessian if verbose and system is small
            call logger%configuration(level=current_log_level)
            if (current_log_level >= verbose_level .and. total_atoms < 20) then
               call logger%info(" ")
               call logger%info("Hessian matrix (Hartree/Bohr^2):")
               do i = 1, 3*total_atoms
                  write (result_line, '(a,i5,a,999f15.8)') "  Row ", i, ": ", (result%hessian(i, j), j=1, 3*total_atoms)
                  call logger%info(trim(result_line))
               end do
               call logger%info(" ")
            end if

            ! Compute and print vibrational analysis
            block
               real(dp), allocatable :: frequencies(:), eigenvalues(:), projected_hessian(:, :)
               real(dp), allocatable :: reduced_masses(:), force_constants(:)
               real(dp), allocatable :: cart_disp(:, :), fc_mdyne(:), ir_intensities(:)
               integer :: ii, jj

               ! First get projected Hessian for verbose output
               call logger%info("  Computing vibrational analysis (projecting trans/rot modes)...")
               call compute_vibrational_frequencies(result%hessian, sys_geom%element_numbers, frequencies, eigenvalues, &
                                                    coordinates=sys_geom%coordinates, project_trans_rot=.true., &
                                                    projected_hessian_out=projected_hessian)

               ! Print projected mass-weighted Hessian if verbose and small system
               if (current_log_level >= verbose_level .and. total_atoms < 20) then
                  if (allocated(projected_hessian)) then
                     call logger%info(" ")
                     call logger%info("Mass-weighted Hessian after trans/rot projection (a.u.):")
                     do ii = 1, 3*total_atoms
                        write (result_line, '(a,i5,a,999f15.8)') "  Row ", ii, ": ", &
                           (projected_hessian(ii, jj), jj=1, 3*total_atoms)
                        call logger%info(trim(result_line))
                     end do
                     call logger%info(" ")
                  end if
               end if

               ! Compute full vibrational analysis and print (with IR intensities if available)
               if (result%has_dipole_derivatives) then
                  call compute_vibrational_analysis(result%hessian, sys_geom%element_numbers, frequencies, &
                                                    reduced_masses, force_constants, cart_disp, &
                                                    coordinates=sys_geom%coordinates, &
                                                    project_trans_rot=.true., &
                                                    force_constants_mdyne=fc_mdyne, &
                                                    dipole_derivatives=result%dipole_derivatives, &
                                                    ir_intensities=ir_intensities)
               else
                  call compute_vibrational_analysis(result%hessian, sys_geom%element_numbers, frequencies, &
                                                    reduced_masses, force_constants, cart_disp, &
                                                    coordinates=sys_geom%coordinates, &
                                                    project_trans_rot=.true., &
                                                    force_constants_mdyne=fc_mdyne)
               end if

               if (allocated(frequencies)) then
                  ! Compute thermochemistry for JSON output
                  block
                     type(thermochemistry_result_t) :: thermo_result
                     integer :: n_modes, n_at

                     n_at = size(sys_geom%element_numbers)
                     n_modes = size(frequencies)

                     call compute_thermochemistry(sys_geom%coordinates, sys_geom%element_numbers, &
                                                  frequencies, n_at, n_modes, thermo_result, &
                                                  temperature=temperature, pressure=pressure)

                     ! Print vibrational analysis to log
                     if (allocated(ir_intensities)) then
                        call print_vibrational_analysis(frequencies, reduced_masses, force_constants, &
                                                        cart_disp, sys_geom%element_numbers, &
                                                        force_constants_mdyne=fc_mdyne, &
                                                        ir_intensities=ir_intensities, &
                                                        coordinates=sys_geom%coordinates, &
                                                        electronic_energy=result%energy%total(), &
                                                        temperature=temperature, pressure=pressure)
                     else
                        call print_vibrational_analysis(frequencies, reduced_masses, force_constants, &
                                                        cart_disp, sys_geom%element_numbers, &
                                                        force_constants_mdyne=fc_mdyne, &
                                                        coordinates=sys_geom%coordinates, &
                                                        electronic_energy=result%energy%total(), &
                                                        temperature=temperature, pressure=pressure)
                     end if

                     ! Populate json_data if present (for centralized JSON output)
                     if (present(json_data)) then
                        json_data%output_mode = OUTPUT_MODE_UNFRAGMENTED
                        json_data%total_energy = result%energy%total()
                        json_data%has_energy = result%has_energy
                        json_data%has_vibrational = .true.

                        ! Copy vibrational data
                        allocate (json_data%frequencies(n_modes))
                        allocate (json_data%reduced_masses(n_modes))
                        allocate (json_data%force_constants(n_modes))
                        json_data%frequencies = frequencies
                        json_data%reduced_masses = reduced_masses
                        json_data%force_constants = fc_mdyne
                        json_data%thermo = thermo_result

                        if (allocated(ir_intensities)) then
                           allocate (json_data%ir_intensities(n_modes))
                           json_data%ir_intensities = ir_intensities
                           json_data%has_ir_intensities = .true.
                        end if

                        ! Copy dipole if available
                        if (result%has_dipole) then
                           allocate (json_data%dipole(3))
                           json_data%dipole = result%dipole
                           json_data%has_dipole = .true.
                        end if

                        ! Copy gradient if available
                        if (result%has_gradient) then
                           allocate (json_data%gradient(3, total_atoms))
                           json_data%gradient = result%gradient
                           json_data%has_gradient = .true.
                        end if

                        ! Copy hessian if available
                        if (result%has_hessian) then
                           allocate (json_data%hessian(3*total_atoms, 3*total_atoms))
                           json_data%hessian = result%hessian
                           json_data%has_hessian = .true.
                        end if
                     end if

                     if (allocated(ir_intensities)) deallocate (ir_intensities)
                  end block
                  deallocate (frequencies, reduced_masses, force_constants, cart_disp, fc_mdyne)
               end if

               if (allocated(eigenvalues)) deallocate (eigenvalues)
               if (allocated(projected_hessian)) deallocate (projected_hessian)
            end block
         end if
      end block
      call logger%info("============================================")

      ! Return result to caller or handle json_data
      if (present(result_out)) then
         ! Transfer result to output (for dynamics/optimization)
         result_out = result
      else
         ! Populate json_data for non-Hessian case if present
         ! (Hessian case already handled above in the vibrational block)
         if (present(json_data) .and. .not. result%has_hessian) then
            json_data%output_mode = OUTPUT_MODE_UNFRAGMENTED
            json_data%total_energy = result%energy%total()
            json_data%has_energy = result%has_energy

            if (result%has_dipole) then
               allocate (json_data%dipole(3))
               json_data%dipole = result%dipole
               json_data%has_dipole = .true.
            end if

            if (result%has_gradient) then
               allocate (json_data%gradient(3, total_atoms))
               json_data%gradient = result%gradient
               json_data%has_gradient = .true.
            end if
         end if
         call result%destroy()
      end if

   end subroutine unfragmented_calculation

end submodule mqc_unfragmented_workflow
