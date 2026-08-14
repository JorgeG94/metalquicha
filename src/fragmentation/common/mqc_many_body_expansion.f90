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
   use mqc_checkpoint, only: checkpoint_t
   implicit none
   private

   public :: many_body_expansion_t
   public :: mbe_context_t
   public :: gmbe_context_t
   public :: fmo_context_t

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
      type(checkpoint_t) :: checkpoint
         !! Fragments already computed by an earlier run, and where the ones
         !! computed by this run are appended. Inactive unless a path was
         !! configured, in which case everything below behaves as before.
         !! System geometry (coordinates, elements, fragments, bonds)

      ! MPI configuration (optional - for distributed calculations)
      type(resources_t), pointer :: resources => null()
         !! MPI communicators and hardware resources
      integer, allocatable :: node_leader_ranks(:)
         !! Ranks of processes that lead each compute node
      integer :: num_nodes = 1
         !! Number of compute nodes
      integer :: global_groups = 1
         !! Number of global groups (multi-global coordinator)
      integer :: nodes_per_group = 0
         !! Nodes per group (optional config override)
      integer, allocatable :: group_leader_ranks(:)
         !! Group leader world rank for each node leader (index matches node_leader_ranks)
      integer, allocatable :: group_ids(:)
         !! Group id for each node leader (index matches node_leader_ranks)

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

   !============================================================================
   ! Fragment molecular orbital method, and electrostatically embedded MBE
   !============================================================================
   type, extends(many_body_expansion_t) :: fmo_context_t
      !! FMO2 and EE-MBE over non-covalently bonded fragments
      !!
      !! Two phases rather than one, which is what separates this from
      !! `mbe_context_t`. The monomers are iterated to self-consistency -- each
      !! solved in the field the others make, which changes its density, which
      !! changes their field -- and only then are the pairs computed, in the
      !! field the settled monomers make. Within a phase the tasks are
      !! independent; between monomer passes there is a barrier and a density
      !! exchange.
      !!
      !! `esp` and `expansion` pick which method this is:
      !!
      !!     "exact" + "fmo"   FMO2
      !!     "ptc"   + "mbe"   electrostatically embedded MBE
      !!     "none"  + "mbe"   plain MBE
      !!
      !! The physics lives in the libcint backend and is reached through
      !! `mqc_libcint_bridge`, so a build without the backend refuses with a
      !! message naming the option rather than failing to link.
      integer, allocatable :: owner(:)
         !! Fragment index per atom, numbered from one with no gaps
      integer :: n_fragments = 0
      character(len=64) :: basis = "6-31g"
      character(len=16) :: esp = "exact"
      character(len=16) :: expansion = "fmo"
      character(len=16) :: far_field = "mulliken"
      real(dp) :: resppc = 2.0_dp
         !! Separation past which a neighbour becomes point charges. Negative
         !! disables the approximation.
      integer :: max_outer = 50
      real(dp) :: outer_tol = 1.0e-7_dp
      real(dp) :: energy = 0.0_dp
         !! What the run produced

   contains
      procedure :: run_serial => fmo_run_serial
      procedure :: run_distributed => fmo_run_distributed
      procedure :: init => fmo_init
      procedure :: destroy => fmo_destroy
   end type fmo_context_t

contains

   !============================================================================
   ! Base class methods
   !============================================================================

   pure function mbe_base_has_mpi(this) result(has_mpi)
      !! Check if MPI resources are available
      class(many_body_expansion_t), intent(in) :: this
      logical :: has_mpi
      has_mpi = associated(this%resources)
   end function mbe_base_has_mpi

   pure function mbe_base_has_geometry(this) result(has_geometry)
      !! Check if system geometry is available
      class(many_body_expansion_t), intent(in) :: this
      logical :: has_geometry
      has_geometry = allocated(this%sys_geom)
   end function mbe_base_has_geometry

   subroutine mbe_base_destroy(this)
      !! Clean up base class resources
      class(many_body_expansion_t), intent(inout) :: this

      if (allocated(this%sys_geom)) then
         call this%sys_geom%destroy()
         deallocate (this%sys_geom)
      end if
      if (allocated(this%node_leader_ranks)) deallocate (this%node_leader_ranks)
      if (allocated(this%group_leader_ranks)) deallocate (this%group_leader_ranks)
      if (allocated(this%group_ids)) deallocate (this%group_ids)

      ! Clear pointers (don't deallocate - we don't own these)
      this%resources => null()
      this%driver_config => null()

      ! Reset scalars
      this%num_nodes = 1
      this%calc_type = 0
      this%global_groups = 1
      this%nodes_per_group = 0
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

   subroutine fmo_init(this, method_config, calc_type)
      !! Initialise an FMO context
      class(fmo_context_t), intent(out) :: this
      type(method_config_t), intent(in) :: method_config
      integer(int32), intent(in) :: calc_type

      this%method_config = method_config
      this%calc_type = calc_type
   end subroutine fmo_init

   subroutine fmo_destroy(this)
      !! Clean up an FMO context
      class(fmo_context_t), intent(inout) :: this

      if (allocated(this%owner)) deallocate (this%owner)
      this%n_fragments = 0
      this%energy = 0.0_dp
      call this%destroy_base()
   end subroutine fmo_destroy

   subroutine fmo_run_serial(this, json_data)
      !! Run FMO on one rank
      use mqc_libcint_bridge, only: run_libcint_fmo
      use mqc_error, only: error_t
      use pic_logger, only: logger => global_logger
      use pic_io, only: to_char

      class(fmo_context_t), intent(inout) :: this
      type(json_output_data_t), intent(out), optional :: json_data

      type(error_t) :: error
      character(len=2), allocatable :: symbols(:)
      integer :: i

      if (.not. this%has_geometry()) then
         call logger%error("fmo_run_serial: sys_geom required but not set")
         return
      end if
      if (.not. allocated(this%owner)) then
         call logger%error("fmo_run_serial: no fragment partition set")
         return
      end if

      allocate (symbols(this%sys_geom%total_atoms))
      do i = 1, this%sys_geom%total_atoms
         symbols(i) = element_symbol_of(this%sys_geom%element_numbers(i))
      end do

      call run_libcint_fmo(this%sys_geom%element_numbers, symbols, &
                           this%sys_geom%coordinates, this%owner, &
                           trim(this%basis), trim(this%esp), trim(this%expansion), &
                           trim(this%far_field), this%resppc, this%max_outer, &
                           this%outer_tol, this%energy, error)
      if (error%has_error()) then
         call logger%error("fmo_run_serial: "//error%get_message())
         return
      end if
      call logger%info("FMO total energy: "//to_char(this%energy)//" Hartree")
      call fmo_report(this, json_data)
   end subroutine fmo_run_serial

   subroutine fmo_report(this, json_data)
      !! Hand the total to whatever writes the output file
      !!
      !! Reported as a fragmented total, which is what it is: the expansion is
      !! different from MBE's but the shape of the answer -- one energy for the
      !! whole system, assembled from fragments -- is the same, and a reader
      !! should not have to learn a third format to find it.
      use mqc_json_output_types, only: OUTPUT_MODE_MBE

      class(fmo_context_t), intent(in) :: this
      type(json_output_data_t), intent(inout), optional :: json_data

      if (.not. present(json_data)) return
      json_data%output_mode = OUTPUT_MODE_MBE
      json_data%total_energy = this%energy
      json_data%has_energy = .true.
   end subroutine fmo_report

   subroutine fmo_run_distributed(this, json_data)
      !! Run FMO across ranks
      !!
      !! The two phases distribute differently and both go through the backend,
      !! which is where the loops are: a monomer pass is a bag of independent
      !! tasks followed by a barrier and a density exchange, and the pair phase
      !! is one bag with no barrier at all. Handing the communicator down rather
      !! than distributing here keeps the fragment geometries from ever crossing
      !! a wire -- every rank assembles what it was asked for from the geometry
      !! it already holds.
      use mqc_libcint_bridge, only: run_libcint_fmo
      use mqc_error, only: error_t
      use pic_logger, only: logger => global_logger
      use pic_io, only: to_char

      class(fmo_context_t), intent(inout) :: this
      type(json_output_data_t), intent(out), optional :: json_data

      type(error_t) :: error
      character(len=2), allocatable :: symbols(:)
      integer :: i

      if (.not. this%has_mpi()) then
         call logger%error("fmo_run_distributed: resources not set in context")
         return
      end if
      if (.not. this%has_geometry()) then
         call logger%error("fmo_run_distributed: sys_geom required but not set")
         return
      end if
      if (.not. allocated(this%owner)) then
         call logger%error("fmo_run_distributed: no fragment partition set")
         return
      end if

      allocate (symbols(this%sys_geom%total_atoms))
      do i = 1, this%sys_geom%total_atoms
         symbols(i) = element_symbol_of(this%sys_geom%element_numbers(i))
      end do

      call run_libcint_fmo(this%sys_geom%element_numbers, symbols, &
                           this%sys_geom%coordinates, this%owner, &
                           trim(this%basis), trim(this%esp), trim(this%expansion), &
                           trim(this%far_field), this%resppc, this%max_outer, &
                           this%outer_tol, this%energy, error, &
                           comm=this%resources%mpi_comms%world_comm)
      if (error%has_error()) then
         call logger%error("fmo_run_distributed: "//error%get_message())
         return
      end if
      if (this%resources%mpi_comms%world_comm%leader()) then
         call logger%info("FMO total energy: "//to_char(this%energy)//" Hartree")
      end if
      ! Every rank has the same total, but only the leader writes a file.
      if (this%resources%mpi_comms%world_comm%leader()) call fmo_report(this, json_data)
   end subroutine fmo_run_distributed

   function element_symbol_of(z) result(sym)
      !! Two-character element symbol, for handing a geometry to the backend
      use mqc_elements, only: element_number_to_symbol
      integer, intent(in) :: z
      character(len=2) :: sym

      sym = element_number_to_symbol(z)
   end function element_symbol_of

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
                                     this%sys_geom, this%method_config, this%calc_type, json_data, &
                                     this%checkpoint)
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
         call logger%verbose("Rank 0: global coordinator, "// &
                             to_char(omp_get_max_threads())//" thread(s)")
         call global_coordinator(this, json_data)
      else if (this%resources%mpi_comms%node_comm%leader()) then
         ! Node coordinator (node leader on other nodes)
         call logger%verbose("Rank "//to_char(this%resources%mpi_comms%world_comm%rank())// &
                             ": Acting as node coordinator")
         call node_coordinator(this)
      else
         ! Worker
         call set_worker_threads(this, "worker")
         call node_worker(this)
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
      use mqc_gmbe_fragment_distribution_scheme, only: gmbe_pie_coordinator, gmbe_group_global_coordinator
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
         call logger%verbose("Rank 0: GMBE PIE coordinator, "// &
                             to_char(omp_get_max_threads())//" thread(s)")
         call gmbe_pie_coordinator(this%resources, this%pie_atom_sets, this%pie_coefficients, &
                                   this%n_pie_terms, this%node_leader_ranks, this%num_nodes, &
                                   this%group_leader_ranks, this%group_ids, this%global_groups, &
                                   this%sys_geom, this%method_config, this%calc_type, json_data)
      else if (this%resources%mpi_comms%node_comm%leader()) then
         ! Node coordinator (node leader on other nodes)
         call logger%verbose("Rank "//to_char(this%resources%mpi_comms%world_comm%rank())// &
                             ": Acting as node coordinator")
         block
            integer :: i, group_leader_rank
            group_leader_rank = 0
            if (allocated(this%node_leader_ranks) .and. allocated(this%group_leader_ranks)) then
               do i = 1, size(this%node_leader_ranks)
                  if (this%node_leader_ranks(i) == this%resources%mpi_comms%world_comm%rank()) then
                     group_leader_rank = this%group_leader_ranks(i)
                     exit
                  end if
               end do
            end if
            if (group_leader_rank == this%resources%mpi_comms%world_comm%rank()) then
               call gmbe_group_global_coordinator(this%resources, this%node_leader_ranks, this%group_ids)
            else
               ! Note: node_coordinator works for both MBE and GMBE
               call node_coordinator(this)
            end if
         end block
      else
         ! Worker
         call set_worker_threads(this, "worker")
         ! Note: node_worker works for both MBE and GMBE (fragment_type distinguishes)
         call node_worker(this)
      end if
   end subroutine gmbe_run_distributed

   subroutine set_worker_threads(this, role)
      !! Give a worker rank its threads, and say what it got
      !!
      !! This used to be an unconditional `omp_set_num_threads(1)`, which made
      !! "one fragment per rank, several threads per rank" impossible: the ranks
      !! that run the chemistry were clamped to one thread whatever the launcher
      !! asked for, so raising OMP_NUM_THREADS did nothing and the cause was
      !! nowhere near the command line.
      !!
      !! The clamp exists for tblite, which corrupts a result when threaded
      !! rather than failing, so it is kept for exactly that -- and the ranks
      !! running anything else now use the threads they were given. A libcint
      !! Fock build threads its own quartet loop, which is where they go.
      !!
      !! Every rank reports its role and thread count because the failure this
      !! replaces was invisible: nothing in the output distinguished "four
      !! threads, no speedup" from "one thread all along".
      use omp_lib, only: omp_set_num_threads, omp_get_max_threads
      use mqc_method_types, only: needs_serial_execution
      use pic_logger, only: logger => global_logger
      use pic_io, only: to_char
      class(many_body_expansion_t), intent(inout) :: this
      character(len=*), intent(in) :: role

      if (needs_serial_execution(this%method_config%method_type)) then
         call omp_set_num_threads(1)
      end if
      call logger%verbose("Rank "//to_char(this%resources%mpi_comms%world_comm%rank())// &
                          ": "//role//", "//to_char(omp_get_max_threads())//" thread(s)")
   end subroutine set_worker_threads

end module mqc_many_body_expansion
