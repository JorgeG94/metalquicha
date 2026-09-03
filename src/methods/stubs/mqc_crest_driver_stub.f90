!! Stand-in for the conformer sampler when the build has no CREST
module mqc_crest_driver
   !! Same name, same procedures as the real backend, and both decline.
   !!
   !! A stub rather than a preprocessor conditional at the call site, matching
   !! `mqc_hdf5_checkpoint_stub` and `mqc_cuest_bridge_stub`: the driver reads
   !! the same either way, and there is exactly one place where "this build
   !! cannot do that" is decided. `crest_available` is how a caller asks in
   !! advance, so the refusal can name the build option rather than failing
   !! somewhere less useful.
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_config_adapter, only: driver_config_t
   use mqc_resources, only: resources_t
   implicit none
   private

   public :: run_conformer_search
   public :: crest_available

contains

   pure function crest_available() result(available)
      !! Whether this build can sample conformers
      logical :: available
      available = .false.
   end function crest_available

   subroutine run_conformer_search(resources, config, sys_geom, error)
      type(resources_t), intent(in) :: resources
      type(driver_config_t), intent(in) :: config
      type(system_geometry_t), intent(in) :: sys_geom
      type(error_t), intent(inout) :: error

      call error%set(ERROR_VALIDATION, "this build cannot sample conformers; "// &
                     "configure with -DMQC_ENABLE_CREST=ON. Every other driver "// &
                     "runs unchanged without it.")
      if (sys_geom%total_atoms < 0 .or. config%nlevel < 0 .or. &
          resources%num_threads < 0) return
   end subroutine run_conformer_search

end module mqc_crest_driver
