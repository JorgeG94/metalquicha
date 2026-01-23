!! Many-Body Expansion abstract base type and concrete implementations
module mqc_many_body_expansion
   !! Provides an abstract base class for all many-body expansion methods
   !! with concrete implementations for standard MBE and generalized MBE (GMBE).
   use pic_types, only: int32, int64, dp
   use mqc_method_config, only: method_config_t
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_resources, only: resources_t
   use mqc_config_adapter, only: driver_config_t
   use mqc_json_output_types, only: json_output_data_t
   implicit none
   private

   public :: many_body_expansion_t
   public :: mbe_context_t
   public :: gmbe_context_t

   !============================================================================
   ! Abstract base type for all many-body expansion methods
   !============================================================================
   type, abstract :: many_body_expansion_t
      !! Abstract base for all many-body expansion methods
      !!
      !! Encapsulates shared configuration for MBE and GMBE calculations.
      !! Concrete implementations provide specific fragment representations
      !! and expansion computation logic.

      ! Required configuration
      type(method_config_t) :: method_config
         !! Method configuration (XTB settings, etc.)
      integer(int32) :: calc_type = 0
         !! Calculation type (energy, gradient, hessian)

      ! System geometry (includes connectivity via sys_geom%bonds)
      type(system_geometry_t), allocatable :: sys_geom
         !! System geometry (coordinates, elements, fragments, bonds)

      ! MPI configuration (optional - for distributed calculations)
      type(resources_t), pointer :: resources => null()
         !! MPI communicators and hardware resources
      integer, allocatable :: node_leader_ranks(:)
         !! Ranks of processes that lead each compute node
      integer :: num_nodes = 1
         !! Number of compute nodes

      ! Driver configuration (for Hessian parameters, etc.)
      type(driver_config_t), pointer :: driver_config => null()
         !! Driver configuration with calculation-specific settings

   contains
      ! Deferred procedures - must be implemented by subclasses
      procedure(run_serial_sub), deferred :: run_serial
      procedure(run_distributed_sub), deferred :: run_distributed

      ! Common procedures
      procedure :: has_mpi => mbe_base_has_mpi
      procedure :: has_geometry => mbe_base_has_geometry
      procedure :: destroy_base => mbe_base_destroy
   end type many_body_expansion_t

   !============================================================================
   ! Abstract interfaces for deferred procedures
   !============================================================================
   abstract interface
      subroutine run_serial_sub(this, json_data)
         import :: many_body_expansion_t, json_output_data_t
         implicit none
         class(many_body_expansion_t), intent(inout) :: this
         type(json_output_data_t), intent(out), optional :: json_data
      end subroutine run_serial_sub

      subroutine run_distributed_sub(this, json_data)
         import :: many_body_expansion_t, json_output_data_t
         implicit none
         class(many_body_expansion_t), intent(inout) :: this
         type(json_output_data_t), intent(out), optional :: json_data
      end subroutine run_distributed_sub
   end interface

   !============================================================================
   ! Standard MBE for non-overlapping fragments
   !============================================================================
   type, extends(many_body_expansion_t) :: mbe_context_t
      !! Standard Many-Body Expansion for non-overlapping fragments
      !!
      !! Uses polymer representation where each fragment is defined by
      !! monomer indices. Coefficients are implicit: (-1)^(n+1) based on
      !! fragment size.

      integer, allocatable :: polymers(:, :)
         !! Fragment composition array (fragment_idx, monomer_indices)
      integer(int64) :: total_fragments = 0
         !! Total number of fragments to process
      integer :: max_level = 0
         !! Maximum MBE level (e.g., 2 for dimers, 3 for trimers)

   contains
      procedure :: run_serial => mbe_run_serial
      procedure :: run_distributed => mbe_run_distributed
      procedure :: init => mbe_init
      procedure :: destroy => mbe_destroy
   end type mbe_context_t

   !============================================================================
   ! Generalized MBE for overlapping fragments (PIE-based)
   !============================================================================
   type, extends(many_body_expansion_t) :: gmbe_context_t
      !! Generalized Many-Body Expansion for overlapping fragments
      !!
      !! Uses PIE (Principle of Inclusion-Exclusion) representation where
      !! each term is defined by atom indices with explicit coefficients
      !! from the inclusion-exclusion principle.

      integer, allocatable :: pie_atom_sets(:, :)
         !! Unique atom sets (max_atoms, n_pie_terms)
      integer, allocatable :: pie_coefficients(:)
         !! PIE coefficient for each term (+1 or -1)
      integer(int64) :: n_pie_terms = 0
         !! Number of unique PIE terms to evaluate

   contains
      procedure :: run_serial => gmbe_run_serial
      procedure :: run_distributed => gmbe_run_distributed
      procedure :: init => gmbe_init
      procedure :: destroy => gmbe_destroy
   end type gmbe_context_t

contains

   !============================================================================
   ! Base class methods
   !============================================================================

   pure logical function mbe_base_has_mpi(this)
      !! Check if MPI resources are available
      class(many_body_expansion_t), intent(in) :: this
      mbe_base_has_mpi = associated(this%resources)
   end function mbe_base_has_mpi

   pure logical function mbe_base_has_geometry(this)
      !! Check if system geometry is available
      class(many_body_expansion_t), intent(in) :: this
      mbe_base_has_geometry = allocated(this%sys_geom)
   end function mbe_base_has_geometry

   subroutine mbe_base_destroy(this)
      !! Clean up base class resources
      class(many_body_expansion_t), intent(inout) :: this

      if (allocated(this%sys_geom)) then
         call this%sys_geom%destroy()
         deallocate (this%sys_geom)
      end if
      if (allocated(this%node_leader_ranks)) deallocate (this%node_leader_ranks)

      ! Clear pointers (don't deallocate - we don't own these)
      this%resources => null()
      this%driver_config => null()

      ! Reset scalars
      this%num_nodes = 1
      this%calc_type = 0
   end subroutine mbe_base_destroy

   !============================================================================
   ! MBE (standard) methods
   !============================================================================

   subroutine mbe_init(this, method_config, calc_type)
      !! Initialize MBE context with required configuration
      class(mbe_context_t), intent(out) :: this
      type(method_config_t), intent(in) :: method_config
      integer(int32), intent(in) :: calc_type

      this%method_config = method_config
      this%calc_type = calc_type
   end subroutine mbe_init

   subroutine mbe_destroy(this)
      !! Clean up MBE context
      class(mbe_context_t), intent(inout) :: this

      ! Clean up MBE-specific data
      if (allocated(this%polymers)) deallocate (this%polymers)
      this%total_fragments = 0
      this%max_level = 0

      ! Clean up base class data
      call this%destroy_base()
   end subroutine mbe_destroy

   subroutine mbe_run_serial(this, json_data)
      !! Run serial MBE calculation
      use mqc_mbe_fragment_distribution_scheme, only: serial_fragment_processor
      use pic_logger, only: logger => global_logger

      class(mbe_context_t), intent(inout) :: this
      type(json_output_data_t), intent(out), optional :: json_data

      if (.not. this%has_geometry()) then
         call logger%error("mbe_run_serial: sys_geom required but not set")
         return
      end if

      call serial_fragment_processor(this%total_fragments, this%polymers, this%max_level, &
                                     this%sys_geom, this%method_config, this%calc_type, json_data)
   end subroutine mbe_run_serial

   subroutine mbe_run_distributed(this, json_data)
      !! Run distributed MBE calculation
      use mqc_mbe_fragment_distribution_scheme, only: global_coordinator, node_coordinator, node_worker
      use pic_logger, only: logger => global_logger
      use pic_io, only: to_char
      use omp_lib, only: omp_set_num_threads, omp_get_max_threads

      class(mbe_context_t), intent(inout) :: this
      type(json_output_data_t), intent(out), optional :: json_data

      if (.not. this%has_mpi()) then
         call logger%error("mbe_run_distributed: resources not set in context")
         return
      end if

      ! Determine role based on MPI rank
      if (this%resources%mpi_comms%world_comm%leader() .and. &
          this%resources%mpi_comms%node_comm%leader()) then
         ! Global coordinator (rank 0, node leader on node 0)
         call omp_set_num_threads(omp_get_max_threads())
         call logger%verbose("Rank 0: Acting as global coordinator")
         if (this%has_geometry()) then
            call global_coordinator(this%resources, this%total_fragments, this%polymers, &
                                    this%max_level, this%node_leader_ranks, this%num_nodes, &
                                    this%sys_geom, this%method_config, this%calc_type, json_data)
         else
            call global_coordinator(this%resources, this%total_fragments, this%polymers, &
                                    this%max_level, this%node_leader_ranks, this%num_nodes, &
                                    method_config=this%method_config, calc_type=this%calc_type, &
                                    json_data=json_data)
         end if
      else if (this%resources%mpi_comms%node_comm%leader()) then
         ! Node coordinator (node leader on other nodes)
         call logger%verbose("Rank "//to_char(this%resources%mpi_comms%world_comm%rank())// &
                             ": Acting as node coordinator")
         call node_coordinator(this%resources, this%method_config, this%calc_type)
      else
         ! Worker
         call omp_set_num_threads(1)
         call logger%verbose("Rank "//to_char(this%resources%mpi_comms%world_comm%rank())// &
                             ": Acting as worker")
         if (this%has_geometry()) then
            call node_worker(this%resources, this%sys_geom, this%method_config, this%calc_type)
         else
            call node_worker(this%resources, method_config=this%method_config, &
                             calc_type=this%calc_type)
         end if
      end if
   end subroutine mbe_run_distributed

   !============================================================================
   ! GMBE (generalized) methods
   !============================================================================

   subroutine gmbe_init(this, method_config, calc_type)
      !! Initialize GMBE context with required configuration
      class(gmbe_context_t), intent(out) :: this
      type(method_config_t), intent(in) :: method_config
      integer(int32), intent(in) :: calc_type

      this%method_config = method_config
      this%calc_type = calc_type
   end subroutine gmbe_init

   subroutine gmbe_destroy(this)
      !! Clean up GMBE context
      class(gmbe_context_t), intent(inout) :: this

      ! Clean up GMBE-specific data
      if (allocated(this%pie_atom_sets)) deallocate (this%pie_atom_sets)
      if (allocated(this%pie_coefficients)) deallocate (this%pie_coefficients)
      this%n_pie_terms = 0

      ! Clean up base class data
      call this%destroy_base()
   end subroutine gmbe_destroy

   subroutine gmbe_run_serial(this, json_data)
      !! Run serial GMBE calculation
      use mqc_gmbe_fragment_distribution_scheme, only: serial_gmbe_pie_processor
      use pic_logger, only: logger => global_logger

      class(gmbe_context_t), intent(inout) :: this
      type(json_output_data_t), intent(out), optional :: json_data

      if (.not. this%has_geometry()) then
         call logger%error("gmbe_run_serial: sys_geom required but not set")
         return
      end if

      call serial_gmbe_pie_processor(this%pie_atom_sets, this%pie_coefficients, &
                                     this%n_pie_terms, this%sys_geom, this%method_config, &
                                     this%calc_type, json_data)
   end subroutine gmbe_run_serial

   subroutine gmbe_run_distributed(this, json_data)
      !! Run distributed GMBE calculation
      use mqc_gmbe_fragment_distribution_scheme, only: gmbe_pie_coordinator
      use mqc_mbe_fragment_distribution_scheme, only: node_coordinator, node_worker
      use pic_logger, only: logger => global_logger
      use pic_io, only: to_char
      use omp_lib, only: omp_set_num_threads, omp_get_max_threads

      class(gmbe_context_t), intent(inout) :: this
      type(json_output_data_t), intent(out), optional :: json_data

      if (.not. this%has_mpi()) then
         call logger%error("gmbe_run_distributed: resources not set in context")
         return
      end if

      ! Determine role based on MPI rank
      if (this%resources%mpi_comms%world_comm%leader() .and. &
          this%resources%mpi_comms%node_comm%leader()) then
         ! Global coordinator (rank 0, node leader on node 0)
         call omp_set_num_threads(omp_get_max_threads())
         call logger%verbose("Rank 0: Acting as GMBE PIE coordinator")
         call gmbe_pie_coordinator(this%resources, this%pie_atom_sets, this%pie_coefficients, &
                                   this%n_pie_terms, this%node_leader_ranks, this%num_nodes, &
                                   this%sys_geom, this%method_config, this%calc_type, json_data)
      else if (this%resources%mpi_comms%node_comm%leader()) then
         ! Node coordinator (node leader on other nodes)
         call logger%verbose("Rank "//to_char(this%resources%mpi_comms%world_comm%rank())// &
                             ": Acting as node coordinator")
         ! Note: node_coordinator works for both MBE and GMBE
         call node_coordinator(this%resources, this%method_config, this%calc_type)
      else
         ! Worker
         call omp_set_num_threads(1)
         call logger%verbose("Rank "//to_char(this%resources%mpi_comms%world_comm%rank())// &
                             ": Acting as worker")
         ! Note: node_worker works for both MBE and GMBE (fragment_type distinguishes)
         if (this%has_geometry()) then
            call node_worker(this%resources, this%sys_geom, this%method_config, this%calc_type)
         else
            call node_worker(this%resources, method_config=this%method_config, &
                             calc_type=this%calc_type)
         end if
      end if
   end subroutine gmbe_run_distributed

end module mqc_many_body_expansion
