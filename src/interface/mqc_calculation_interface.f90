!! External calculation interface for geometry optimization, AIMD, and Monte Carlo
module mqc_calculation_interface
   !! Energies and forces for a geometry, in the shape an optimizer, an MD
   !! integrator or an MC sampler wants them.
   use pic_types, only: int32, dp
   use pic_mpi_lib, only: comm_t, bcast, abort_comm
   use pic_logger, only: logger => global_logger
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_config_types, only: bond_t
   use mqc_result_types, only: calculation_result_t
   use mqc_calc_types, only: CALC_TYPE_ENERGY, CALC_TYPE_GRADIENT, CALC_TYPE_HESSIAN
   use mqc_resources, only: resources_t

   implicit none
   private

   public :: compute_energy_and_forces
   public :: sync_geometry_to_workers

contains

   subroutine sync_geometry_to_workers(sys_geom, comm)
      !! Where a broadcast of updated coordinates to the workers would go
      !!
      !! Does nothing today.
      ! TODO(mqc): a no-op that touches neither of its arguments, called on
      ! every energy of an optimization or dynamics loop. The ranks agree on a
      ! geometry only if `run_calculation` distributes it, so a caller reading
      ! the name of this routine is told something that is not being done.
      type(system_geometry_t), intent(inout) :: sys_geom
      type(comm_t), intent(in) :: comm

   end subroutine sync_geometry_to_workers

   subroutine compute_energy_and_forces(sys_geom, driver_config, resources, &
                                        energy, gradient, hessian, bonds, write_output)
      !! Compute energy and forces for the current geometry
      !!
      !! Called on all ranks. The master rank provides the updated geometry and
      !! is the only one where `energy`, `gradient` and `hessian` are set.
      use mqc_driver, only: run_calculation
      use mqc_config_adapter, only: driver_config_t
      type(system_geometry_t), intent(inout) :: sys_geom
      type(driver_config_t), intent(in) :: driver_config
      type(resources_t), intent(in) :: resources
      real(dp), intent(out) :: energy
      real(dp), intent(out), optional :: gradient(:, :)  !! (3, total_atoms)
      real(dp), intent(out), optional :: hessian(:, :)  !! (3*total_atoms, 3*total_atoms)
      type(bond_t), intent(in), optional :: bonds(:)
      logical, intent(in), optional :: write_output
         !! Whether the driver writes its JSON output file. Defaults to what
         !! `run_calculation` defaults to, which is yes. A caller driving this
         !! in a loop wants no: one write per call, each overwriting the last.

      type(calculation_result_t) :: result
      logical :: need_gradient, need_hessian

      ! Determine what we need based on what's requested
      need_gradient = present(gradient)
      need_hessian = present(hessian)

      ! Synchronize geometry from master to all ranks
      ! (Master may have updated coordinates for optimization/dynamics)
      call sync_geometry_to_workers(sys_geom, resources%mpi_comms%world_comm)

      ! Call the main calculation driver
      ! This handles both fragmented and unfragmented cases
      call run_calculation(resources, driver_config, sys_geom, bonds, result, &
                           write_output=write_output)

      ! Extract results (only valid on master rank)
      if (resources%mpi_comms%world_comm%rank() == 0) then
         energy = result%energy%total()

         if (need_gradient .and. result%has_gradient) then
            gradient = result%gradient
         else if (need_gradient) then
            call logger%error("Gradient requested but not computed!")
            call abort_comm(resources%mpi_comms%world_comm, 1)
         end if

         if (need_hessian .and. result%has_hessian) then
            hessian = result%hessian
         else if (need_hessian) then
            call logger%error("Hessian requested but not computed!")
            call abort_comm(resources%mpi_comms%world_comm, 1)
         end if

         ! Clean up
         call result%destroy()
      end if

   end subroutine compute_energy_and_forces

end module mqc_calculation_interface
