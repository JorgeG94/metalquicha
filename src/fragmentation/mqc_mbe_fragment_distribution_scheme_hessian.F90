submodule(mqc_mbe_fragment_distribution_scheme) mqc_hessian_distribution_scheme
   implicit none

contains

   module subroutine distributed_unfragmented_hessian(world_comm, sys_geom, method, driver_config)
      !! Compute Hessian for unfragmented system using MPI distribution
      !!
      !! Uses a dynamic work queue approach: workers request displacement indices
      !! from rank 0, compute gradients, and send results back. This provides
      !! better load balancing than static work distribution.
      use mqc_finite_differences, only: generate_perturbed_geometries, displaced_geometry_t, &
                                        finite_diff_hessian_from_gradients, DEFAULT_DISPLACEMENT, &
                                        copy_and_displace_geometry
      use mqc_config_adapter, only: driver_config_t
#ifndef MQC_WITHOUT_TBLITE
      use mqc_method_xtb, only: xtb_method_t
#endif

      type(comm_t), intent(in) :: world_comm
      type(system_geometry_t), intent(in) :: sys_geom
      integer(int32), intent(in) :: method
      type(driver_config_t), intent(in), optional :: driver_config  !! Driver configuration

      integer :: my_rank, n_ranks
      real(dp) :: displacement

      my_rank = world_comm%rank()
      n_ranks = world_comm%size()

      ! Use provided displacement or default
      if (present(driver_config)) then
         displacement = driver_config%hessian%displacement
      else
         displacement = DEFAULT_DISPLACEMENT
      end if

      if (my_rank == 0) then
         ! Rank 0 is the coordinator
         call hessian_coordinator(world_comm, sys_geom, method, displacement)
      else
         ! Other ranks are workers
         call hessian_worker(world_comm, sys_geom, method, displacement)
      end if

      ! Synchronize all ranks before returning
      call world_comm%barrier()

   end subroutine distributed_unfragmented_hessian

   module subroutine hessian_coordinator(world_comm, sys_geom, method, displacement)
      !! Coordinator for distributed Hessian calculation
      !! Distributes displacement work and collects gradient results
      use mqc_finite_differences, only: finite_diff_hessian_from_gradients
#ifndef MQC_WITHOUT_TBLITE
      use mqc_method_xtb, only: xtb_method_t
#endif

      type(comm_t), intent(in) :: world_comm
      type(system_geometry_t), intent(in) :: sys_geom
      integer(int32), intent(in) :: method
      real(dp), intent(in) :: displacement  !! Finite difference displacement (Bohr)

      type(physical_fragment_t) :: full_system
      type(timer_type) :: coord_timer
      real(dp), allocatable :: forward_gradients(:, :, :)  ! (n_displacements, 3, n_atoms)
      real(dp), allocatable :: backward_gradients(:, :, :)  ! (n_displacements, 3, n_atoms)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: grad_buffer(:, :)
      type(calculation_result_t) :: result
      integer :: n_atoms, n_displacements, n_ranks
      integer :: current_disp, finished_workers, dummy_msg, worker_rank
      integer :: disp_idx, gradient_type  ! gradient_type: 1=forward, 2=backward
      type(MPI_Status) :: status
      logical :: has_pending
      type(request_t) :: req
      integer :: current_log_level
      logical :: is_verbose
      character(len=256) :: result_line
      real(dp) :: hess_norm
      integer :: i, j
#ifndef MQC_WITHOUT_TBLITE
      type(xtb_method_t) :: xtb_calc
#endif

      n_ranks = world_comm%size()
      n_atoms = sys_geom%total_atoms
      n_displacements = 3*n_atoms

      call logger%configuration(level=current_log_level)
      is_verbose = (current_log_level >= verbose_level)

      call logger%info("============================================")
      call logger%info("Distributed unfragmented Hessian calculation")
      call logger%info("  Total atoms: "//to_char(n_atoms))
      call logger%info("  Gradient calculations needed: "//to_char(2*n_displacements))
      call logger%info("  MPI ranks: "//to_char(n_ranks))
      call logger%info("  Work distribution: Dynamic queue")
      call logger%info("============================================")

      ! Build full system geometry
      full_system%n_atoms = n_atoms
      full_system%n_caps = 0
      allocate (full_system%element_numbers(n_atoms))
      allocate (full_system%coordinates(3, n_atoms))
      full_system%element_numbers = sys_geom%element_numbers
      full_system%coordinates = sys_geom%coordinates
      full_system%charge = sys_geom%charge
      full_system%multiplicity = sys_geom%multiplicity
      call full_system%compute_nelec()

      ! Allocate storage for all gradients
      allocate (forward_gradients(n_displacements, 3, n_atoms))
      allocate (backward_gradients(n_displacements, 3, n_atoms))
      allocate (grad_buffer(3, n_atoms))

      current_disp = 1
      finished_workers = 0

      ! Process work requests and collect results
      call coord_timer%start()
      do while (finished_workers < n_ranks - 1)

         ! Check for incoming gradient results
         call iprobe(world_comm, MPI_ANY_SOURCE, TAG_WORKER_SCALAR_RESULT, has_pending, status)
         if (has_pending) then
            worker_rank = status%MPI_SOURCE

            ! Receive: displacement index, gradient type (1=forward, 2=backward), gradient data
            call irecv(world_comm, disp_idx, worker_rank, TAG_WORKER_SCALAR_RESULT, req)
            call wait(req)
            call irecv(world_comm, gradient_type, worker_rank, TAG_WORKER_SCALAR_RESULT, req)
            call wait(req)
            call recv(world_comm, grad_buffer, worker_rank, TAG_WORKER_SCALAR_RESULT, status)

            ! Store gradient in appropriate array
            if (gradient_type == 1) then
               forward_gradients(disp_idx, :, :) = grad_buffer
            else
               backward_gradients(disp_idx, :, :) = grad_buffer
            end if

            ! Log progress every 10% or at completion (count both forward and backward)
            if (gradient_type == 2) then  ! Only log after backward gradient to count complete displacements
               if (mod(disp_idx, max(1, n_displacements/10)) == 0 .or. disp_idx == n_displacements) then
                  call logger%info("  Completed "//to_char(disp_idx)//"/"//to_char(n_displacements)// &
                                   " displacement pairs in "//to_char(coord_timer%get_elapsed_time())//" s")
               end if
            end if
         end if

         ! Check for work requests from workers
         call iprobe(world_comm, MPI_ANY_SOURCE, TAG_WORKER_REQUEST, has_pending, status)
         if (has_pending) then
            worker_rank = status%MPI_SOURCE
            call irecv(world_comm, dummy_msg, worker_rank, TAG_WORKER_REQUEST, req)
            call wait(req)

            if (current_disp <= n_displacements) then
               ! Send next displacement index to worker
               call isend(world_comm, current_disp, worker_rank, TAG_WORKER_FRAGMENT, req)
               call wait(req)
               current_disp = current_disp + 1
            else
               ! No more work - tell worker to finish
               call isend(world_comm, -1, worker_rank, TAG_WORKER_FINISH, req)
               call wait(req)
               finished_workers = finished_workers + 1
            end if
         end if
      end do

      deallocate (grad_buffer)
      call coord_timer%stop()
      call logger%info("All gradient calculations completed in "//to_char(coord_timer%get_elapsed_time())//" s")

      ! Assemble Hessian from finite differences
      call logger%info("  Assembling Hessian matrix...")
      call coord_timer%start()
      call finite_diff_hessian_from_gradients(full_system, forward_gradients, backward_gradients, &
                                              displacement, hessian)
      call coord_timer%stop()
      call logger%info("Hessian assembly completed in "//to_char(coord_timer%get_elapsed_time())//" s")

      ! Compute energy and gradient at reference geometry
      call logger%info("  Computing reference energy and gradient...")
#ifndef MQC_WITHOUT_TBLITE
      xtb_calc%variant = method_type_to_string(method)
      xtb_calc%verbose = is_verbose
      call xtb_calc%calc_gradient(full_system, result)
#endif

      ! Store Hessian in result
      if (allocated(result%hessian)) deallocate (result%hessian)
      allocate (result%hessian(size(hessian, 1), size(hessian, 2)))
      result%hessian = hessian
      result%has_hessian = .true.

      ! Print results
      call logger%info("============================================")
      call logger%info("Distributed Hessian calculation completed")
      write (result_line, '(a,f25.15)') "  Final energy: ", result%energy%total()
      call logger%info(trim(result_line))

      if (result%has_gradient) then
         write (result_line, '(a,f25.15)') "  Gradient norm: ", sqrt(sum(result%gradient**2))
         call logger%info(trim(result_line))
      end if

      if (result%has_hessian) then
         hess_norm = sqrt(sum(result%hessian**2))
         write (result_line, '(a,f25.15)') "  Hessian Frobenius norm: ", hess_norm
         call logger%info(trim(result_line))

         if (is_verbose .and. n_atoms < 20) then
            call logger%info(" ")
            call logger%info("Hessian matrix (Hartree/Bohr^2):")
            do i = 1, 3*n_atoms
               write (result_line, '(a,i5,a,999f15.8)') "  Row ", i, ": ", (result%hessian(i, j), j=1, 3*n_atoms)
               call logger%info(trim(result_line))
            end do
            call logger%info(" ")
         end if
      end if
      call logger%info("============================================")

      ! Output JSON
      call print_unfragmented_json(result)

      ! Cleanup
      call result%destroy()
      deallocate (forward_gradients, backward_gradients)
      if (allocated(hessian)) deallocate (hessian)

   end subroutine hessian_coordinator

   module subroutine hessian_worker(world_comm, sys_geom, method, displacement)
      !! Worker for distributed Hessian calculation
      !! Requests displacement indices, computes gradients, and sends results back
      use mqc_finite_differences, only: copy_and_displace_geometry
#ifndef MQC_WITHOUT_TBLITE
      use mqc_method_xtb, only: xtb_method_t
#endif

      type(comm_t), intent(in) :: world_comm
      type(system_geometry_t), intent(in) :: sys_geom
      integer(int32), intent(in) :: method
      real(dp), intent(in) :: displacement  !! Finite difference displacement (Bohr)

      type(physical_fragment_t) :: full_system, displaced_geom
      type(calculation_result_t) :: grad_result
      integer :: n_atoms, disp_idx, atom_idx, coord, gradient_type, dummy_msg
      type(MPI_Status) :: status
      type(request_t) :: req
#ifndef MQC_WITHOUT_TBLITE
      type(xtb_method_t) :: xtb_calc
#endif

      n_atoms = sys_geom%total_atoms

      ! Build full system geometry
      full_system%n_atoms = n_atoms
      full_system%n_caps = 0
      allocate (full_system%element_numbers(n_atoms))
      allocate (full_system%coordinates(3, n_atoms))
      full_system%element_numbers = sys_geom%element_numbers
      full_system%coordinates = sys_geom%coordinates
      full_system%charge = sys_geom%charge
      full_system%multiplicity = sys_geom%multiplicity
      call full_system%compute_nelec()

#ifndef MQC_WITHOUT_TBLITE
      ! Setup XTB method
      xtb_calc%variant = method_type_to_string(method)
      xtb_calc%verbose = .false.

      dummy_msg = 0
      do
         ! Request work from coordinator
         call isend(world_comm, dummy_msg, 0, TAG_WORKER_REQUEST, req)
         call wait(req)
         call irecv(world_comm, disp_idx, 0, MPI_ANY_TAG, req)
         call wait(req, status)

         if (status%MPI_TAG == TAG_WORKER_FINISH) exit

         ! Compute displacement index to atom and coordinate
         atom_idx = (disp_idx - 1)/3 + 1
         coord = mod(disp_idx - 1, 3) + 1

         ! Compute FORWARD gradient
         call copy_and_displace_geometry(full_system, atom_idx, coord, displacement, displaced_geom)
         call xtb_calc%calc_gradient(displaced_geom, grad_result)

         if (.not. grad_result%has_gradient) then
            call logger%error("Worker failed gradient for displacement "//to_char(disp_idx))
            call abort_comm(world_comm, 1)
         end if

         ! Send: displacement index, gradient type (1=forward), gradient data
         gradient_type = 1
         call isend(world_comm, disp_idx, 0, TAG_WORKER_SCALAR_RESULT, req)
         call wait(req)
         call isend(world_comm, gradient_type, 0, TAG_WORKER_SCALAR_RESULT, req)
         call wait(req)
         call send(world_comm, grad_result%gradient, 0, TAG_WORKER_SCALAR_RESULT)

         call grad_result%destroy()
         call displaced_geom%destroy()

         ! Compute BACKWARD gradient
         call copy_and_displace_geometry(full_system, atom_idx, coord, -displacement, displaced_geom)
         call xtb_calc%calc_gradient(displaced_geom, grad_result)

         if (.not. grad_result%has_gradient) then
            call logger%error("Worker failed gradient for displacement "//to_char(disp_idx))
            call abort_comm(world_comm, 1)
         end if

         ! Send: displacement index, gradient type (2=backward), gradient data
         gradient_type = 2
         call isend(world_comm, disp_idx, 0, TAG_WORKER_SCALAR_RESULT, req)
         call wait(req)
         call isend(world_comm, gradient_type, 0, TAG_WORKER_SCALAR_RESULT, req)
         call wait(req)
         call send(world_comm, grad_result%gradient, 0, TAG_WORKER_SCALAR_RESULT)

         call grad_result%destroy()
         call displaced_geom%destroy()
      end do
#else
      call logger%error("XTB method requested but tblite support not compiled in")
      call abort_comm(world_comm, 1)
#endif

   end subroutine hessian_worker
end submodule mqc_hessian_distribution_scheme
