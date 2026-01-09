!! External calculation interface for geometry optimization, AIMD, and Monte Carlo
module mqc_calculation_interface
   !! Provides a clean interface for computing energies and forces
   !! that can be used by optimization algorithms, MD integrators, and MC samplers
   use pic_types, only: int32, dp
   use pic_mpi_lib, only: comm_t, bcast, abort_comm
   use pic_logger, only: logger => global_logger
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_config_parser, only: bond_t
   use mqc_result_types, only: calculation_result_t
   use mqc_calc_types, only: CALC_TYPE_ENERGY, CALC_TYPE_GRADIENT, CALC_TYPE_HESSIAN

   implicit none
   private

   public :: compute_energy_and_forces
   public :: sync_geometry_to_workers

contains

   subroutine sync_geometry_to_workers(sys_geom, comm)
      !! Synchronize geometry coordinates from master rank to all worker ranks
      !! This is needed when master rank updates coordinates for optimization/dynamics
      !!
      !! TODO: Implement explicit broadcast if needed. Currently, the fragmented
      !! calculation infrastructure may already handle geometry distribution.
      !! For unfragmented calculations on master rank only, this is not needed.
      type(system_geometry_t), intent(inout) :: sys_geom
      type(comm_t), intent(in) :: comm

      ! NOTE: For now, we rely on the existing calculation infrastructure
      ! to handle geometry as needed. If explicit broadcasting is required,
      ! we can add MPI send/recv logic here later.

   end subroutine sync_geometry_to_workers

   subroutine compute_energy_and_forces(sys_geom, driver_config, world_comm, node_comm, &
                                        energy, gradient, hessian, bonds)
      !! Compute energy and forces for current geometry
      !! This is the main interface for optimization/dynamics codes
      !!
      !! Master rank provides updated geometry, all ranks compute fragments,
      !! results are returned on master rank only
      !!
      !! Usage:
      !!   1. Master rank updates sys_geom%coordinates
      !!   2. Call this subroutine (all ranks)
      !!   3. Master rank receives energy/gradient/hessian
      !!   4. Master rank updates geometry based on forces
      !!   5. Repeat from step 1
      use mqc_driver, only: run_calculation
      use mqc_config_adapter, only: driver_config_t
      type(system_geometry_t), intent(inout) :: sys_geom
      type(driver_config_t), intent(in) :: driver_config
      type(comm_t), intent(in) :: world_comm, node_comm
      real(dp), intent(out) :: energy
      real(dp), intent(out), optional :: gradient(:, :)  !! (3, total_atoms)
      real(dp), intent(out), optional :: hessian(:, :)  !! (3*total_atoms, 3*total_atoms)
      type(bond_t), intent(in), optional :: bonds(:)

      type(calculation_result_t) :: result
      logical :: need_gradient, need_hessian

      ! Determine what we need based on what's requested
      need_gradient = present(gradient)
      need_hessian = present(hessian)

      ! Synchronize geometry from master to all ranks
      ! (Master may have updated coordinates for optimization/dynamics)
      call sync_geometry_to_workers(sys_geom, world_comm)

      ! Call the main calculation driver
      ! This handles both fragmented and unfragmented cases
      call run_calculation(world_comm, node_comm, driver_config, sys_geom, bonds, result)

      ! Extract results (only valid on master rank)
      if (world_comm%rank() == 0) then
         energy = result%energy%total()

         if (need_gradient .and. result%has_gradient) then
            gradient = result%gradient
         else if (need_gradient) then
            call logger%error("Gradient requested but not computed!")
            call abort_comm(world_comm, 1)
         end if

         if (need_hessian .and. result%has_hessian) then
            hessian = result%hessian
         else if (need_hessian) then
            call logger%error("Hessian requested but not computed!")
            call abort_comm(world_comm, 1)
         end if

         ! Clean up
         call result%destroy()
      end if

   end subroutine compute_energy_and_forces

end module mqc_calculation_interface
