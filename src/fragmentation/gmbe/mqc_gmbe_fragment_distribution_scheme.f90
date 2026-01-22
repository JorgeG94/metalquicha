!! Generalized Many-Body Expansion (GMBE) fragment distribution module
module mqc_gmbe_fragment_distribution_scheme
   !! Implements fragment distribution schemes for GMBE calculations with overlapping fragments
   !! Handles both serial and MPI-parallelized distribution of monomers and intersection fragments
   use pic_types, only: int32, int64, dp
   use pic_timer, only: timer_type
   use mqc_calc_types, only: CALC_TYPE_GRADIENT
   use pic_mpi_lib, only: comm_t, send, recv, isend, irecv, &
                          wait, iprobe, MPI_Status, request_t, MPI_ANY_SOURCE, MPI_ANY_TAG, abort_comm
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_mpi_tags, only: TAG_WORKER_REQUEST, TAG_WORKER_FRAGMENT, TAG_WORKER_FINISH, &
                           TAG_WORKER_SCALAR_RESULT, &
                           TAG_NODE_REQUEST, TAG_NODE_FRAGMENT, TAG_NODE_FINISH, &
                           TAG_NODE_SCALAR_RESULT
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, build_fragment_from_indices, &
                                    build_fragment_from_atom_list
   use mqc_config_parser, only: bond_t
   use mqc_result_types, only: calculation_result_t, result_send, result_isend, &
                               result_recv, result_irecv, mbe_result_t
   use mqc_mbe_fragment_distribution_scheme, only: do_fragment_work
   use mqc_method_config, only: method_config_t
   use mqc_json_output_types, only: json_output_data_t, OUTPUT_MODE_GMBE_PIE
   use mqc_vibrational_analysis, only: compute_vibrational_analysis, print_vibrational_analysis
   use mqc_thermochemistry, only: thermochemistry_result_t, compute_thermochemistry
   use mqc_calculation_defaults, only: FRAGMENT_TYPE_ATOMS
   implicit none
   ! Error handling imported where needed
   private

   ! Public interface
   public :: serial_gmbe_pie_processor  !! PIE-based serial processor
   public :: gmbe_pie_coordinator  !! PIE-based MPI coordinator

contains

   subroutine serial_gmbe_pie_processor(pie_atom_sets, pie_coefficients, n_pie_terms, &
                                        sys_geom, method_config, calc_type, bonds, json_data)
      !! Serial GMBE processor using PIE coefficients
      !! Evaluates each unique atom set once and sums with PIE coefficients
      !! Supports energy-only, energy+gradient, and energy+gradient+Hessian calculations
      !! If json_data is present, populates it for centralized JSON output
      use mqc_calc_types, only: CALC_TYPE_GRADIENT, CALC_TYPE_HESSIAN, CALC_TYPE_ENERGY, calc_type_to_string
      use mqc_physical_fragment, only: redistribute_cap_gradients, redistribute_cap_hessian, &
                                       redistribute_cap_dipole_derivatives
      use mqc_error, only: error_t
      use pic_logger, only: info_level
      integer, intent(in) :: pie_atom_sets(:, :)  !! Unique atom sets (max_atoms, n_pie_terms)
      integer, intent(in) :: pie_coefficients(:)  !! PIE coefficient for each term
      integer(int64), intent(in) :: n_pie_terms
      type(system_geometry_t), intent(in) :: sys_geom
      type(method_config_t), intent(in) :: method_config  !! Method configuration
      integer(int32), intent(in) :: calc_type
      type(bond_t), intent(in), optional :: bonds(:)
      type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data

      type(physical_fragment_t) :: phys_frag
      type(calculation_result_t), allocatable :: results(:)
      type(error_t) :: error
      integer :: n_atoms, max_atoms, iatom, current_log_level, hess_dim
      integer(int64) :: term_idx
      integer, allocatable :: atom_list(:)
      real(dp) :: total_energy, term_energy
      real(dp), allocatable :: pie_energies(:)  !! Store individual energies for JSON output
      real(dp), allocatable :: total_gradient(:, :)  !! Total gradient (3, total_atoms)
      real(dp), allocatable :: term_gradient(:, :)  !! Temporary gradient for each term
      real(dp), allocatable :: total_hessian(:, :)  !! Total Hessian (3*total_atoms, 3*total_atoms)
      real(dp), allocatable :: term_hessian(:, :)  !! Temporary Hessian for each term
      real(dp), allocatable :: total_dipole_derivs(:, :)  !! Total dipole derivatives (3, 3*total_atoms)
      real(dp), allocatable :: term_dipole_derivs(:, :)   !! Temporary dipole derivatives for each term
      real(dp), allocatable :: ir_intensities(:)          !! IR intensities in km/mol
      integer :: coeff

      if (int(size(pie_atom_sets, 2), int64) < n_pie_terms .or. &
          int(size(pie_coefficients), int64) < n_pie_terms) then
         call logger%error("PIE term arrays are smaller than n_pie_terms")
         error stop "Invalid PIE term array sizes"
      end if

      call logger%info("Processing "//to_char(n_pie_terms)//" unique PIE terms...")
      call logger%info("  Calculation type: "//calc_type_to_string(calc_type))

      total_energy = 0.0_dp
      max_atoms = size(pie_atom_sets, 1)
      allocate (pie_energies(n_pie_terms))
      allocate (results(n_pie_terms))

      ! Allocate gradient and Hessian arrays if needed
      if (calc_type == CALC_TYPE_GRADIENT .or. calc_type == CALC_TYPE_HESSIAN) then
         allocate (total_gradient(3, sys_geom%total_atoms))
         allocate (term_gradient(3, sys_geom%total_atoms))
         total_gradient = 0.0_dp
      end if

      if (calc_type == CALC_TYPE_HESSIAN) then
         hess_dim = 3*sys_geom%total_atoms
         allocate (total_hessian(hess_dim, hess_dim))
         allocate (term_hessian(hess_dim, hess_dim))
         total_hessian = 0.0_dp
         ! Allocate dipole derivative arrays for IR intensities
         allocate (total_dipole_derivs(3, hess_dim))
         allocate (term_dipole_derivs(3, hess_dim))
         total_dipole_derivs = 0.0_dp
      end if

      do term_idx = 1_int64, n_pie_terms
         coeff = pie_coefficients(term_idx)

         ! Skip terms with zero coefficient (shouldn't happen, but safety check)
         if (coeff == 0) then
            pie_energies(term_idx) = 0.0_dp  ! Mark as skipped
            cycle
         end if

         ! Extract atom list for this term
         n_atoms = 0
         do while (n_atoms < max_atoms .and. pie_atom_sets(n_atoms + 1, term_idx) >= 0)
            n_atoms = n_atoms + 1
         end do

         if (n_atoms == 0) then
            pie_energies(term_idx) = 0.0_dp  ! Mark as skipped
            cycle
         end if

         allocate (atom_list(n_atoms))
         atom_list = pie_atom_sets(1:n_atoms, term_idx)

         ! Build fragment from atom list
         call build_fragment_from_atom_list(sys_geom, atom_list, n_atoms, phys_frag, error, bonds)
         if (error%has_error()) then
            call logger%error(error%get_full_trace())
            error stop "Failed to build intersection fragment"
         end if

         ! Compute energy (and gradient if requested)
         call do_fragment_work(term_idx, results(term_idx), method_config, phys_frag, calc_type)

         ! Check for calculation errors
         if (results(term_idx)%has_error) then
            call logger%error("PIE term "//to_char(term_idx)//" calculation failed: "// &
                              results(term_idx)%error%get_message())
            error stop "PIE term calculation failed in serial processing"
         end if

         term_energy = results(term_idx)%energy%total()

         ! Store energy for JSON output
         pie_energies(term_idx) = term_energy

         ! Accumulate with PIE coefficient
         total_energy = total_energy + real(coeff, dp)*term_energy

         ! Accumulate gradient if present
         if ((calc_type == CALC_TYPE_GRADIENT .or. calc_type == CALC_TYPE_HESSIAN) .and. &
             results(term_idx)%has_gradient) then
            ! Map fragment gradient to system coordinates with proper cap handling
            term_gradient = 0.0_dp
            call redistribute_cap_gradients(phys_frag, results(term_idx)%gradient, term_gradient)

            ! Accumulate with PIE coefficient
            total_gradient = total_gradient + real(coeff, dp)*term_gradient
         end if

         ! Accumulate Hessian if present
         if (calc_type == CALC_TYPE_HESSIAN .and. results(term_idx)%has_hessian) then
            ! Map fragment Hessian to system coordinates with proper cap handling
            term_hessian = 0.0_dp
            call redistribute_cap_hessian(phys_frag, results(term_idx)%hessian, term_hessian)

            ! Accumulate with PIE coefficient
            total_hessian = total_hessian + real(coeff, dp)*term_hessian

            ! Accumulate dipole derivatives if present (for IR intensities)
            if (results(term_idx)%has_dipole_derivatives) then
               term_dipole_derivs = 0.0_dp
               call redistribute_cap_dipole_derivatives(phys_frag, results(term_idx)%dipole_derivatives, &
                                                        term_dipole_derivs)
               total_dipole_derivs = total_dipole_derivs + real(coeff, dp)*term_dipole_derivs
            end if
         end if

         call logger%verbose("PIE term "//to_char(term_idx)//"/"//to_char(n_pie_terms)// &
                             ": "//to_char(n_atoms)//" atoms, coeff="//to_char(coeff)// &
                             ", E="//to_char(term_energy))

         deallocate (atom_list)
         call phys_frag%destroy()
      end do

      call logger%info(" ")
      call logger%info("GMBE PIE calculation completed successfully")
      call logger%info("Final GMBE energy: "//to_char(total_energy)//" Hartree")

      ! Print gradient info if computed
      if (calc_type == CALC_TYPE_GRADIENT .or. calc_type == CALC_TYPE_HESSIAN) then
         call logger%info("GMBE PIE gradient computation completed")
         call logger%info("  Total gradient norm: "//to_char(sqrt(sum(total_gradient**2))))

         ! Print detailed gradient if info level and small system
         call logger%configuration(level=current_log_level)
         if (current_log_level >= info_level .and. sys_geom%total_atoms < 100) then
            call logger%info(" ")
            call logger%info("Total GMBE PIE Gradient (Hartree/Bohr):")
            do iatom = 1, sys_geom%total_atoms
               block
                  character(len=256) :: grad_line
                  write (grad_line, '(a,i5,a,3f20.12)') "  Atom ", iatom, ": ", &
                     total_gradient(1, iatom), total_gradient(2, iatom), total_gradient(3, iatom)
                  call logger%info(trim(grad_line))
               end block
            end do
            call logger%info(" ")
         end if
      end if

      ! Print Hessian info if computed
      if (calc_type == CALC_TYPE_HESSIAN) then
         call logger%info("GMBE PIE Hessian computation completed")
         call logger%info("  Total Hessian Frobenius norm: "//to_char(sqrt(sum(total_hessian**2))))

         ! Compute and print full vibrational analysis with thermochemistry
         block
            real(dp), allocatable :: frequencies(:), reduced_masses(:), force_constants(:)
            real(dp), allocatable :: cart_disp(:, :), fc_mdyne(:)
            type(thermochemistry_result_t) :: thermo_result
            type(mbe_result_t) :: gmbe_result
            integer :: n_at, n_modes

            call logger%info("  Computing vibrational analysis (projecting trans/rot modes)...")
            call compute_vibrational_analysis(total_hessian, sys_geom%element_numbers, frequencies, &
                                              reduced_masses, force_constants, cart_disp, &
                                              coordinates=sys_geom%coordinates, &
                                              project_trans_rot=.true., &
                                              force_constants_mdyne=fc_mdyne, &
                                              dipole_derivatives=total_dipole_derivs, &
                                              ir_intensities=ir_intensities)

            if (allocated(frequencies)) then
               ! Compute thermochemistry
               n_at = size(sys_geom%element_numbers)
               n_modes = size(frequencies)
               call compute_thermochemistry(sys_geom%coordinates, sys_geom%element_numbers, &
                                            frequencies, n_at, n_modes, thermo_result)

               ! Print vibrational analysis to log
               call print_vibrational_analysis(frequencies, reduced_masses, force_constants, &
                                               cart_disp, sys_geom%element_numbers, &
                                               force_constants_mdyne=fc_mdyne, &
                                               ir_intensities=ir_intensities, &
                                               coordinates=sys_geom%coordinates, &
                                               electronic_energy=total_energy)

               ! Build temporary mbe_result for JSON output
               gmbe_result%total_energy = total_energy
               gmbe_result%has_energy = .true.
               gmbe_result%has_hessian = .true.
               if (allocated(total_gradient)) then
                  gmbe_result%has_gradient = .true.
                  allocate (gmbe_result%gradient, source=total_gradient)
               end if
               allocate (gmbe_result%hessian, source=total_hessian)

               ! Populate json_data for vibrational output if present
               if (present(json_data)) then
                  json_data%output_mode = OUTPUT_MODE_GMBE_PIE
                  json_data%total_energy = total_energy
                  json_data%has_energy = .true.
                  json_data%has_vibrational = .true.

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

                  if (allocated(total_gradient)) then
                     allocate (json_data%gradient, source=total_gradient)
                     json_data%has_gradient = .true.
                  end if

                  allocate (json_data%hessian, source=total_hessian)
                  json_data%has_hessian = .true.
               end if

               if (allocated(ir_intensities)) deallocate (ir_intensities)
               call gmbe_result%destroy()
               deallocate (frequencies, reduced_masses, force_constants, cart_disp, fc_mdyne)
            end if
         end block
      end if

      call logger%info(" ")

      ! Populate json_data for non-Hessian case if present
      if (present(json_data) .and. calc_type /= CALC_TYPE_HESSIAN) then
         json_data%output_mode = OUTPUT_MODE_GMBE_PIE
         json_data%total_energy = total_energy
         json_data%has_energy = .true.
         json_data%n_pie_terms = n_pie_terms

         ! Copy PIE data
         allocate (json_data%pie_atom_sets, source=pie_atom_sets(:, 1:n_pie_terms))
         allocate (json_data%pie_coefficients(n_pie_terms))
         json_data%pie_coefficients = pie_coefficients(1:n_pie_terms)
         allocate (json_data%pie_energies(n_pie_terms))
         json_data%pie_energies = pie_energies

         if (allocated(total_gradient)) then
            allocate (json_data%gradient, source=total_gradient)
            json_data%has_gradient = .true.
         end if
         if (allocated(total_hessian)) then
            allocate (json_data%hessian, source=total_hessian)
            json_data%has_hessian = .true.
         end if
      end if

      deallocate (pie_energies, results)
      if (allocated(total_gradient)) deallocate (total_gradient)
      if (allocated(term_gradient)) deallocate (term_gradient)
      if (allocated(total_hessian)) deallocate (total_hessian)
      if (allocated(term_hessian)) deallocate (term_hessian)
      if (allocated(total_dipole_derivs)) deallocate (total_dipole_derivs)
      if (allocated(term_dipole_derivs)) deallocate (term_dipole_derivs)

   end subroutine serial_gmbe_pie_processor

   subroutine gmbe_pie_coordinator(resources, pie_atom_sets, pie_coefficients, n_pie_terms, &
                                   node_leader_ranks, num_nodes, sys_geom, method_config, calc_type, bonds, json_data)
      !! MPI coordinator for PIE-based GMBE calculations
      !! Distributes PIE terms across MPI ranks and accumulates results
      !! If json_data is present, populates it for centralized JSON output
      use mqc_calc_types, only: CALC_TYPE_GRADIENT, CALC_TYPE_HESSIAN
      use mqc_physical_fragment, only: redistribute_cap_gradients, redistribute_cap_hessian, &
                                       redistribute_cap_dipole_derivatives
      use mqc_resources, only: resources_t

      type(resources_t), intent(in) :: resources
      integer, intent(in) :: pie_atom_sets(:, :)  !! Unique atom sets (max_atoms, n_pie_terms)
      integer, intent(in) :: pie_coefficients(:)  !! PIE coefficient for each term
      integer(int64), intent(in) :: n_pie_terms
      integer, intent(in) :: node_leader_ranks(:), num_nodes
      type(system_geometry_t), intent(in) :: sys_geom
      type(method_config_t), intent(in) :: method_config  !! Method configuration
      integer(int32), intent(in) :: calc_type
      type(bond_t), intent(in), optional :: bonds(:)
      type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data

      type(timer_type) :: coord_timer
      integer(int64) :: current_term_idx, results_received, term_idx
      integer :: finished_nodes
      integer :: request_source, dummy_msg
      type(MPI_Status) :: status, local_status
      logical :: handling_local_workers, has_pending
      integer :: local_finished_workers, local_dummy

      ! Storage for results
      type(calculation_result_t), allocatable :: results(:)
      integer(int64) :: worker_term_map(resources%mpi_comms%node_comm%size())
      integer :: worker_source
      real(dp) :: total_energy
      real(dp), allocatable :: total_gradient(:, :)
      real(dp), allocatable :: total_hessian(:, :)
      real(dp), allocatable :: total_dipole_derivs(:, :)  !! Total dipole derivatives (3, 3*total_atoms)
      real(dp), allocatable :: ir_intensities(:)          !! IR intensities in km/mol
      integer :: hess_dim

      ! MPI request handles
      type(request_t) :: req

      if (int(size(pie_atom_sets, 2), int64) < n_pie_terms .or. &
          int(size(pie_coefficients), int64) < n_pie_terms) then
         call logger%error("PIE term arrays are smaller than n_pie_terms")
         call abort_comm(resources%mpi_comms%world_comm, 1)
      end if

      current_term_idx = n_pie_terms
      finished_nodes = 0
      local_finished_workers = 0
      handling_local_workers = (resources%mpi_comms%node_comm%size() > 1)
      results_received = 0_int64
      worker_term_map = 0

      allocate (results(n_pie_terms))

      call logger%verbose("GMBE PIE coordinator starting with "//to_char(n_pie_terms)// &
                          " PIE terms for "//to_char(num_nodes)//" nodes")

      call coord_timer%start()
      do while (finished_nodes < num_nodes)

         ! PRIORITY 1: Check for incoming results from local workers
         if (handling_local_workers) then
            do
           call iprobe(resources%mpi_comms%node_comm, MPI_ANY_SOURCE, TAG_WORKER_SCALAR_RESULT, has_pending, local_status)
               if (.not. has_pending) exit

               worker_source = local_status%MPI_SOURCE
               if (worker_term_map(worker_source) == 0) then
                  call logger%error("Received result from worker "//to_char(worker_source)// &
                                    " but no term was assigned!")
                  call abort_comm(resources%mpi_comms%world_comm, 1)
               end if

               call result_irecv(results(worker_term_map(worker_source)), resources%mpi_comms%node_comm, worker_source, &
                                 TAG_WORKER_SCALAR_RESULT, req)
               call wait(req)

               ! Check for calculation errors from worker
               if (results(worker_term_map(worker_source))%has_error) then
                  call logger%error("PIE term "//to_char(worker_term_map(worker_source))// &
                                    " calculation failed: "// &
                                    results(worker_term_map(worker_source))%error%get_message())
                  call abort_comm(resources%mpi_comms%world_comm, 1)
               end if

               worker_term_map(worker_source) = 0
               results_received = results_received + 1
               if (mod(results_received, max(1_int64, n_pie_terms/10_int64)) == 0 .or. &
                   results_received == n_pie_terms) then
                  call logger%info("  Processed "//to_char(results_received)//"/"// &
                                   to_char(n_pie_terms)//" PIE terms ["// &
                                   to_char(coord_timer%get_elapsed_time())//" s]")
               end if
            end do
         end if

         ! PRIORITY 1b: Check for incoming results from remote node coordinators
         do
            call iprobe(resources%mpi_comms%world_comm, MPI_ANY_SOURCE, TAG_NODE_SCALAR_RESULT, has_pending, status)
            if (.not. has_pending) exit

            call irecv(resources%mpi_comms%world_comm, term_idx, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT, req)
            call wait(req)
      call result_irecv(results(term_idx), resources%mpi_comms%world_comm, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT, req)
            call wait(req)

            ! Check for calculation errors from node coordinator
            if (results(term_idx)%has_error) then
               call logger%error("PIE term "//to_char(term_idx)//" calculation failed: "// &
                                 results(term_idx)%error%get_message())
               call abort_comm(resources%mpi_comms%world_comm, 1)
            end if

            results_received = results_received + 1
            if (mod(results_received, max(1_int64, n_pie_terms/10_int64)) == 0 .or. &
                results_received == n_pie_terms) then
               call logger%info("  Processed "//to_char(results_received)//"/"// &
                                to_char(n_pie_terms)//" PIE terms ["// &
                                to_char(coord_timer%get_elapsed_time())//" s]")
            end if
         end do

         ! PRIORITY 2: Remote node coordinator requests
         call iprobe(resources%mpi_comms%world_comm, MPI_ANY_SOURCE, TAG_NODE_REQUEST, has_pending, status)
         if (has_pending) then
            call irecv(resources%mpi_comms%world_comm, dummy_msg, status%MPI_SOURCE, TAG_NODE_REQUEST, req)
            call wait(req)
            request_source = status%MPI_SOURCE

            if (current_term_idx >= 1) then
               call send_pie_term_to_node(resources%mpi_comms%world_comm, current_term_idx, pie_atom_sets, request_source)
               current_term_idx = current_term_idx - 1
            else
               call isend(resources%mpi_comms%world_comm, -1, request_source, TAG_NODE_FINISH, req)
               call wait(req)
               finished_nodes = finished_nodes + 1
            end if
         end if

         ! PRIORITY 3: Local workers (shared memory) - send new work
         if (handling_local_workers .and. local_finished_workers < resources%mpi_comms%node_comm%size() - 1) then
            call iprobe(resources%mpi_comms%node_comm, MPI_ANY_SOURCE, TAG_WORKER_REQUEST, has_pending, local_status)
            if (has_pending) then
               if (worker_term_map(local_status%MPI_SOURCE) == 0) then
                  call irecv(resources%mpi_comms%node_comm, local_dummy, local_status%MPI_SOURCE, TAG_WORKER_REQUEST, req)
                  call wait(req)

                  if (current_term_idx >= 1) then
                     call send_pie_term_to_worker(resources%mpi_comms%node_comm, &
                                                  current_term_idx, pie_atom_sets, local_status%MPI_SOURCE)
                     worker_term_map(local_status%MPI_SOURCE) = current_term_idx
                     current_term_idx = current_term_idx - 1
                  else
                     call isend(resources%mpi_comms%node_comm, -1, local_status%MPI_SOURCE, TAG_WORKER_FINISH, req)
                     call wait(req)
                     local_finished_workers = local_finished_workers + 1
                  end if
               end if
            end if
         end if

         ! Finalize local worker completion
         if (handling_local_workers .and. local_finished_workers >= resources%mpi_comms%node_comm%size() - 1 &
             .and. results_received >= n_pie_terms) then
            handling_local_workers = .false.
            finished_nodes = finished_nodes + 1
         end if
      end do

      call logger%verbose("GMBE PIE coordinator finished all terms")
      call coord_timer%stop()
      call logger%info("Time to evaluate all PIE terms "//to_char(coord_timer%get_elapsed_time())//" s")

      ! Accumulate results with PIE coefficients
      call logger%info(" ")
      call logger%info("Computing GMBE PIE energy...")
      call coord_timer%start()

      total_energy = 0.0_dp
      do term_idx = 1_int64, n_pie_terms
         total_energy = total_energy + real(pie_coefficients(term_idx), dp)*results(term_idx)%energy%total()
      end do

      ! Handle gradients if computed
      if (calc_type == CALC_TYPE_GRADIENT) then
         allocate (total_gradient(3, sys_geom%total_atoms))
         total_gradient = 0.0_dp

         do term_idx = 1_int64, n_pie_terms
            if (results(term_idx)%has_gradient) then
               ! Map fragment gradient to system coordinates
               block
                  use mqc_error, only: error_t
                  real(dp), allocatable :: term_gradient(:, :)
                  type(physical_fragment_t) :: phys_frag
                  type(error_t) :: error
                  integer :: n_atoms, max_atoms
                  integer, allocatable :: atom_list(:)

                  allocate (term_gradient(3, sys_geom%total_atoms))
                  term_gradient = 0.0_dp

                  ! Extract atom list for this term
                  max_atoms = size(pie_atom_sets, 1)
                  n_atoms = 0
                  do while (n_atoms < max_atoms .and. pie_atom_sets(n_atoms + 1, term_idx) >= 0)
                     n_atoms = n_atoms + 1
                  end do

                  if (n_atoms > 0) then
                     allocate (atom_list(n_atoms))
                     atom_list = pie_atom_sets(1:n_atoms, term_idx)

                     ! Build fragment to get proper mapping
                     call build_fragment_from_atom_list(sys_geom, atom_list, n_atoms, phys_frag, error, bonds)
                     call redistribute_cap_gradients(phys_frag, results(term_idx)%gradient, term_gradient)
                     call phys_frag%destroy()
                     deallocate (atom_list)
                  end if

                  ! Accumulate with PIE coefficient
                  total_gradient = total_gradient + real(pie_coefficients(term_idx), dp)*term_gradient
                  deallocate (term_gradient)
               end block
            end if
         end do

         ! Print gradient information
         call logger%info("GMBE PIE gradient computation completed")
         call logger%info("  Total gradient norm: "//to_char(sqrt(sum(total_gradient**2))))

         ! Print detailed gradient if info level and small system
         block
            use pic_logger, only: info_level
            integer :: iatom, current_log_level
            call logger%configuration(level=current_log_level)
            if (current_log_level >= info_level .and. sys_geom%total_atoms < 100) then
               call logger%info(" ")
               call logger%info("Total GMBE PIE Gradient (Hartree/Bohr):")
               do iatom = 1, sys_geom%total_atoms
                  block
                     character(len=256) :: grad_line
                     write (grad_line, '(a,i5,a,3f20.12)') "  Atom ", iatom, ": ", &
                        total_gradient(1, iatom), total_gradient(2, iatom), total_gradient(3, iatom)
                     call logger%info(trim(grad_line))
                  end block
               end do
               call logger%info(" ")
            end if
         end block

         deallocate (total_gradient)
      end if

      ! Handle Hessians if computed
      if (calc_type == CALC_TYPE_HESSIAN) then
         hess_dim = 3*sys_geom%total_atoms
         allocate (total_hessian(hess_dim, hess_dim))
         total_hessian = 0.0_dp

         ! Also allocate gradient for Hessian calculations
         if (.not. allocated(total_gradient)) then
            allocate (total_gradient(3, sys_geom%total_atoms))
            total_gradient = 0.0_dp
         end if

         ! Allocate dipole derivative arrays for IR intensities
         allocate (total_dipole_derivs(3, hess_dim))
         total_dipole_derivs = 0.0_dp

         do term_idx = 1_int64, n_pie_terms
            if (results(term_idx)%has_hessian .or. results(term_idx)%has_gradient) then
               block
                  use mqc_error, only: error_t
                  real(dp), allocatable :: term_gradient(:, :), term_hessian(:, :), term_dipole_derivs(:, :)
                  type(physical_fragment_t) :: phys_frag
                  type(error_t) :: error
                  integer :: n_atoms, max_atoms
                  integer, allocatable :: atom_list(:)

                  ! Extract atom list for this term
                  max_atoms = size(pie_atom_sets, 1)
                  n_atoms = 0
                  do while (n_atoms < max_atoms .and. pie_atom_sets(n_atoms + 1, term_idx) >= 0)
                     n_atoms = n_atoms + 1
                  end do

                  if (n_atoms > 0) then
                     allocate (atom_list(n_atoms))
                     atom_list = pie_atom_sets(1:n_atoms, term_idx)

                     ! Build fragment to get proper mapping
                     call build_fragment_from_atom_list(sys_geom, atom_list, n_atoms, phys_frag, error, bonds)

                     ! Redistribute gradient if present
                     if (results(term_idx)%has_gradient) then
                        allocate (term_gradient(3, sys_geom%total_atoms))
                        term_gradient = 0.0_dp
                        call redistribute_cap_gradients(phys_frag, results(term_idx)%gradient, term_gradient)
                        total_gradient = total_gradient + real(pie_coefficients(term_idx), dp)*term_gradient
                        deallocate (term_gradient)
                     end if

                     ! Redistribute Hessian if present
                     if (results(term_idx)%has_hessian) then
                        allocate (term_hessian(hess_dim, hess_dim))
                        term_hessian = 0.0_dp
                        call redistribute_cap_hessian(phys_frag, results(term_idx)%hessian, term_hessian)
                        total_hessian = total_hessian + real(pie_coefficients(term_idx), dp)*term_hessian
                        deallocate (term_hessian)

                        ! Accumulate dipole derivatives if present (for IR intensities)
                        if (results(term_idx)%has_dipole_derivatives) then
                           allocate (term_dipole_derivs(3, hess_dim))
                           term_dipole_derivs = 0.0_dp
                           call redistribute_cap_dipole_derivatives(phys_frag, &
                                                                    results(term_idx)%dipole_derivatives, &
                                                                    term_dipole_derivs)
                           total_dipole_derivs = total_dipole_derivs + &
                                                 real(pie_coefficients(term_idx), dp)*term_dipole_derivs
                           deallocate (term_dipole_derivs)
                        end if
                     end if

                     call phys_frag%destroy()
                     deallocate (atom_list)
                  end if
               end block
            end if
         end do

         ! Print gradient information
         call logger%info("GMBE PIE gradient computation completed")
         call logger%info("  Total gradient norm: "//to_char(sqrt(sum(total_gradient**2))))

         ! Print Hessian information
         call logger%info("GMBE PIE Hessian computation completed")
         call logger%info("  Total Hessian Frobenius norm: "//to_char(sqrt(sum(total_hessian**2))))

         ! Compute and print full vibrational analysis with thermochemistry
         block
            real(dp), allocatable :: frequencies(:), reduced_masses(:), force_constants(:)
            real(dp), allocatable :: cart_disp(:, :), fc_mdyne(:)
            type(thermochemistry_result_t) :: thermo_result
            type(mbe_result_t) :: gmbe_result
            integer :: n_at, n_modes

            call logger%info("  Computing vibrational analysis (projecting trans/rot modes)...")
            call compute_vibrational_analysis(total_hessian, sys_geom%element_numbers, frequencies, &
                                              reduced_masses, force_constants, cart_disp, &
                                              coordinates=sys_geom%coordinates, &
                                              project_trans_rot=.true., &
                                              force_constants_mdyne=fc_mdyne, &
                                              dipole_derivatives=total_dipole_derivs, &
                                              ir_intensities=ir_intensities)

            if (allocated(frequencies)) then
               ! Compute thermochemistry
               n_at = size(sys_geom%element_numbers)
               n_modes = size(frequencies)
               call compute_thermochemistry(sys_geom%coordinates, sys_geom%element_numbers, &
                                            frequencies, n_at, n_modes, thermo_result)

               ! Print vibrational analysis to log
               call print_vibrational_analysis(frequencies, reduced_masses, force_constants, &
                                               cart_disp, sys_geom%element_numbers, &
                                               force_constants_mdyne=fc_mdyne, &
                                               ir_intensities=ir_intensities, &
                                               coordinates=sys_geom%coordinates, &
                                               electronic_energy=total_energy)

               ! Build temporary mbe_result for JSON output
               gmbe_result%total_energy = total_energy
               gmbe_result%has_energy = .true.
               gmbe_result%has_hessian = .true.
               if (allocated(total_gradient)) then
                  gmbe_result%has_gradient = .true.
                  allocate (gmbe_result%gradient, source=total_gradient)
               end if
               allocate (gmbe_result%hessian, source=total_hessian)

               ! Populate json_data for vibrational output if present
               if (present(json_data)) then
                  json_data%output_mode = OUTPUT_MODE_GMBE_PIE
                  json_data%total_energy = total_energy
                  json_data%has_energy = .true.
                  json_data%has_vibrational = .true.

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

                  if (allocated(total_gradient)) then
                     allocate (json_data%gradient, source=total_gradient)
                     json_data%has_gradient = .true.
                  end if

                  allocate (json_data%hessian, source=total_hessian)
                  json_data%has_hessian = .true.
               end if

               if (allocated(ir_intensities)) deallocate (ir_intensities)
               call gmbe_result%destroy()
               deallocate (frequencies, reduced_masses, force_constants, cart_disp, fc_mdyne)
            end if
         end block
         if (allocated(total_dipole_derivs)) deallocate (total_dipole_derivs)
      end if

      call coord_timer%stop()
      call logger%info("Time to compute GMBE PIE "//to_char(coord_timer%get_elapsed_time())//" s")
      call logger%info(" ")
      call logger%info("GMBE PIE calculation completed successfully")
      call logger%info("Final GMBE energy: "//to_char(total_energy)//" Hartree")
      call logger%info(" ")

      ! Populate json_data for non-Hessian case if present
      if (present(json_data) .and. calc_type /= CALC_TYPE_HESSIAN) then
         block
            real(dp), allocatable :: pie_energies(:)
            allocate (pie_energies(n_pie_terms))
            do term_idx = 1_int64, n_pie_terms
               pie_energies(term_idx) = results(term_idx)%energy%total()
            end do

            json_data%output_mode = OUTPUT_MODE_GMBE_PIE
            json_data%total_energy = total_energy
            json_data%has_energy = .true.
            json_data%n_pie_terms = n_pie_terms

            allocate (json_data%pie_atom_sets, source=pie_atom_sets(:, 1:n_pie_terms))
            allocate (json_data%pie_coefficients(n_pie_terms))
            json_data%pie_coefficients = pie_coefficients(1:n_pie_terms)
            allocate (json_data%pie_energies(n_pie_terms))
            json_data%pie_energies = pie_energies

            if (allocated(total_gradient)) then
               allocate (json_data%gradient, source=total_gradient)
               json_data%has_gradient = .true.
            end if
            if (allocated(total_hessian)) then
               allocate (json_data%hessian, source=total_hessian)
               json_data%has_hessian = .true.
            end if

            deallocate (pie_energies)
         end block
      end if

      deallocate (results)
      if (allocated(total_gradient)) deallocate (total_gradient)
      if (allocated(total_hessian)) deallocate (total_hessian)

   end subroutine gmbe_pie_coordinator

   subroutine send_pie_term_to_node(world_comm, term_idx, pie_atom_sets, dest_rank)
      !! Send PIE term (atom list) to remote node coordinator
      type(comm_t), intent(in) :: world_comm
      integer(int64), intent(in) :: term_idx
      integer, intent(in) :: pie_atom_sets(:, :)
      integer, intent(in) :: dest_rank

      integer :: n_atoms, max_atoms
      integer, allocatable :: atom_list(:)
      integer(int32) :: fragment_type
      type(request_t) :: req(4)

      ! PIE terms always use atom lists
      fragment_type = FRAGMENT_TYPE_ATOMS

      ! Extract atom list for this term
      max_atoms = size(pie_atom_sets, 1)
      n_atoms = 0
      do while (n_atoms < max_atoms .and. pie_atom_sets(n_atoms + 1, term_idx) >= 0)
         n_atoms = n_atoms + 1
      end do

      allocate (atom_list(n_atoms))
      atom_list = pie_atom_sets(1:n_atoms, term_idx)

      call isend(world_comm, term_idx, dest_rank, TAG_NODE_FRAGMENT, req(1))
      call isend(world_comm, fragment_type, dest_rank, TAG_NODE_FRAGMENT, req(2))
      call isend(world_comm, n_atoms, dest_rank, TAG_NODE_FRAGMENT, req(3))
      call isend(world_comm, atom_list, dest_rank, TAG_NODE_FRAGMENT, req(4))

      call wait(req(1))
      call wait(req(2))
      call wait(req(3))
      call wait(req(4))

      deallocate (atom_list)
   end subroutine send_pie_term_to_node

   subroutine send_pie_term_to_worker(node_comm, term_idx, pie_atom_sets, dest_rank)
      !! Send PIE term (atom list) to local worker
      type(comm_t), intent(in) :: node_comm
      integer(int64), intent(in) :: term_idx
      integer, intent(in) :: pie_atom_sets(:, :)
      integer, intent(in) :: dest_rank

      integer :: n_atoms, max_atoms
      integer, allocatable :: atom_list(:)
      integer(int32) :: fragment_type
      type(request_t) :: req(4)

      ! PIE terms always use atom lists
      fragment_type = FRAGMENT_TYPE_ATOMS

      ! Extract atom list for this term
      max_atoms = size(pie_atom_sets, 1)
      n_atoms = 0
      do while (n_atoms < max_atoms .and. pie_atom_sets(n_atoms + 1, term_idx) >= 0)
         n_atoms = n_atoms + 1
      end do

      allocate (atom_list(n_atoms))
      atom_list = pie_atom_sets(1:n_atoms, term_idx)

      call isend(node_comm, term_idx, dest_rank, TAG_WORKER_FRAGMENT, req(1))
      call isend(node_comm, fragment_type, dest_rank, TAG_WORKER_FRAGMENT, req(2))
      call isend(node_comm, n_atoms, dest_rank, TAG_WORKER_FRAGMENT, req(3))
      call isend(node_comm, atom_list, dest_rank, TAG_WORKER_FRAGMENT, req(4))

      call wait(req(1))
      call wait(req(2))
      call wait(req(3))
      call wait(req(4))

      deallocate (atom_list)
   end subroutine send_pie_term_to_worker

end module mqc_gmbe_fragment_distribution_scheme
