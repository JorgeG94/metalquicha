!! Main calculation driver module for metalquicha
module mqc_driver
   !! Handles both fragmented (many-body expansion) and unfragmented calculations
   !! with MPI parallelization and node-based work distribution.
   use pic_types, only: int32, int64, dp
   use pic_mpi_lib, only: comm_t, abort_comm, bcast, allgather
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   use mqc_mbe_fragment_distribution_scheme, only: global_coordinator, node_coordinator, node_worker, unfragmented_calculation, &
                                             serial_fragment_processor, do_fragment_work, distributed_unfragmented_hessian
   use mqc_gmbe_fragment_distribution_scheme, only: serial_gmbe_pie_processor, gmbe_pie_coordinator
   use mqc_frag_utils, only: get_nfrags, create_monomer_list, generate_fragment_list, generate_intersections, &
                             gmbe_enumerate_pie_terms, binomial, combine, apply_distance_screening, sort_fragments_by_size
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, &
                                    build_fragment_from_indices, build_fragment_from_atom_list
   use mqc_config_adapter, only: driver_config_t, config_to_driver, config_to_system_geometry
   use mqc_method_types, only: method_type_to_string
   use mqc_calc_types, only: calc_type_to_string, CALC_TYPE_GRADIENT, CALC_TYPE_HESSIAN
   use mqc_config_parser, only: bond_t, mqc_config_t
   use mqc_mbe, only: compute_gmbe_energy
   use mqc_result_types, only: calculation_result_t
   use mqc_error, only: error_t
   use mqc_io_helpers, only: set_molecule_suffix, get_output_json_filename
   use mqc_json, only: merge_multi_molecule_json
   implicit none
   private

   public :: run_calculation  !! Main entry point for all calculations
   public :: run_multi_molecule_calculations  !! Multi-molecule calculation dispatcher

contains

   subroutine run_calculation(world_comm, node_comm, config, sys_geom, bonds, result_out)
      !! Main calculation dispatcher - routes to fragmented or unfragmented calculation
      !!
      !! Determines calculation type based on nlevel and dispatches to appropriate
      !! calculation routine with proper MPI setup and validation.
      !! If result_out is present, returns result instead of writing JSON (for dynamics/optimization)
      type(comm_t), intent(in) :: world_comm  !! Global MPI communicator
      type(comm_t), intent(in) :: node_comm   !! Node-local MPI communicator
      type(driver_config_t), intent(in) :: config  !! Driver configuration
      type(system_geometry_t), intent(in) :: sys_geom  !! System geometry and fragment info
      type(bond_t), intent(in), optional :: bonds(:)  !! Bond connectivity information
      type(calculation_result_t), intent(out), optional :: result_out  !! Optional result output

      ! Local variables
      integer :: max_level   !! Maximum fragment level (nlevel from config)
      integer :: i  !! Loop counter

      ! Set max_level from config
      max_level = config%nlevel

      if (world_comm%rank() == 0) then
         call logger%info("============================================")
         call logger%info("Loaded geometry:")
         call logger%info("  Total monomers: "//to_char(sys_geom%n_monomers))
         call logger%info("  Atoms per monomer: "//to_char(sys_geom%atoms_per_monomer))
         call logger%info("  Total atoms: "//to_char(sys_geom%total_atoms))
         call logger%info("  Fragment level: "//to_char(max_level))
         call logger%info("============================================")
      end if

      ! Warn if overlapping fragments flag is set but nlevel=0
      if (config%allow_overlapping_fragments .and. max_level == 0) then
         if (world_comm%rank() == 0) then
            call logger%warning("allow_overlapping_fragments is set to true, but nlevel=0")
            call logger%warning("Running unfragmented calculation - overlapping fragments flag will be ignored")
         end if
      end if

      ! GMBE (overlapping fragments) with inclusion-exclusion principle
      ! GMBE(1): Base fragments are monomers
      ! GMBE(N): Base fragments are N-mers (e.g., dimers for N=2)
      ! Algorithm: Generate primaries, use DFS to enumerate overlapping cliques,
      ! accumulate PIE coefficients per unique atom set, evaluate each once

      if (max_level == 0) then
         call omp_set_num_threads(1)
       call run_unfragmented_calculation(world_comm, sys_geom, config%method, config%calc_type, bonds, config, result_out)
      else
         call run_fragmented_calculation(world_comm, node_comm, config%method, config%calc_type, sys_geom, max_level, &
                                         config%allow_overlapping_fragments, &
                                         config%max_intersection_level, bonds, config)
      end if

   end subroutine run_calculation

   subroutine run_unfragmented_calculation(world_comm, sys_geom, method, calc_type, bonds, driver_config, result_out)
      !! Handle unfragmented calculation (nlevel=0)
      !!
      !! For single-molecule mode: Only rank 0 runs (validates single rank)
      !! For multi-molecule mode: ALL ranks can run (each with their own molecule)
      !! For Hessian calculations with multiple ranks: Uses distributed parallelization
      !! If result_out is present, returns result instead of writing JSON
      type(comm_t), intent(in) :: world_comm  !! Global MPI communicator
      type(system_geometry_t), intent(in) :: sys_geom  !! Complete system geometry
      integer(int32), intent(in) :: method  !! Quantum chemistry method
      integer(int32), intent(in) :: calc_type  !! Calculation type
      type(bond_t), intent(in), optional :: bonds(:)  !! Bond connectivity information
      type(driver_config_t), intent(in), optional :: driver_config  !! Driver configuration
      type(calculation_result_t), intent(out), optional :: result_out  !! Optional result output

      ! For Hessian calculations with multiple ranks, use distributed approach
      if (calc_type == CALC_TYPE_HESSIAN .and. world_comm%size() > 1) then
         if (world_comm%rank() == 0) then
            call logger%info(" ")
            call logger%info("Running distributed unfragmented Hessian calculation")
            call logger%info("  MPI ranks: "//to_char(world_comm%size()))
            call logger%info(" ")
         end if
         call distributed_unfragmented_hessian(world_comm, sys_geom, method, driver_config)
         return
      end if

      ! Check if this is multi-molecule mode or single-molecule mode
      ! In multi-molecule mode, each rank processes its own molecule
      ! In single-molecule mode, only rank 0 should work
      if (world_comm%size() == 1 .or. world_comm%rank() == 0) then
         ! Either single-rank calculation, or rank 0 in multi-rank setup
         call logger%info(" ")
         call logger%info("Running unfragmented calculation")
         call logger%info("  Calculation type: "//calc_type_to_string(calc_type))
         call logger%info(" ")
         call unfragmented_calculation(sys_geom, method, calc_type, bonds, result_out)
      else if (sys_geom%total_atoms > 0) then
         ! Multi-molecule mode: non-zero rank with a molecule
         call logger%verbose("Rank "//to_char(world_comm%rank())//": Running unfragmented calculation")
         call unfragmented_calculation(sys_geom, method, calc_type, bonds, result_out)
      end if

   end subroutine run_unfragmented_calculation

   subroutine run_fragmented_calculation(world_comm, node_comm, method, calc_type, sys_geom, max_level, &
                                         allow_overlapping_fragments, max_intersection_level, bonds, driver_config)
      !! Handle fragmented calculation (nlevel > 0)
      !!
      !! Generates fragments, distributes work across MPI processes organized in nodes,
      !! and coordinates many-body expansion calculation using hierarchical parallelism.
      !! If allow_overlapping_fragments=true, uses GMBE with intersection correction.
      type(comm_t), intent(in) :: world_comm  !! Global MPI communicator
      type(comm_t), intent(in) :: node_comm   !! Node-local MPI communicator
      integer(int32), intent(in) :: method  !! Quantum chemistry method
      integer(int32), intent(in) :: calc_type  !! Calculation type
      type(system_geometry_t), intent(in) :: sys_geom  !! System geometry and fragment info
      integer, intent(in) :: max_level    !! Maximum fragment level for MBE
      logical, intent(in) :: allow_overlapping_fragments  !! Use GMBE for overlapping fragments
      integer, intent(in) :: max_intersection_level  !! Maximum k-way intersection depth for GMBE
      type(bond_t), intent(in), optional :: bonds(:)  !! Bond connectivity information
      type(driver_config_t), intent(in) :: driver_config  !! Driver configuration with cutoffs

      integer(int64) :: total_fragments  !! Total number of fragments generated (int64 to handle large systems)
      integer, allocatable :: polymers(:, :)  !! Fragment composition array (fragment, monomer_indices)
      integer :: num_nodes   !! Number of compute nodes
      integer :: i, j        !! Loop counters
      integer, allocatable :: node_leader_ranks(:)  !! Ranks of processes that lead each node
      integer, allocatable :: monomers(:)     !! Temporary monomer list for fragment generation
      integer(int64) :: n_expected_frags  !! Expected number of fragments based on combinatorics (int64 to handle large systems)
      integer(int64) :: n_rows      !! Number of rows needed for polymers array (int64 to handle large systems)
      integer :: global_node_rank  !! Global rank if this process leads a node, -1 otherwise
      integer, allocatable :: all_node_leader_ranks(:)  !! Node leader status for all ranks

      ! GMBE-specific variables (old approach - kept for compatibility)
      integer, allocatable :: intersections(:, :)  !! Intersection atom lists (max_atoms, n_intersections)
      integer, allocatable :: intersection_sets(:, :)  !! k-tuples for each intersection (n_monomers, n_intersections)
      integer, allocatable :: intersection_levels(:)  !! Level k of each intersection (n_intersections)
      integer :: n_intersections, n_monomers  !! Counts for GMBE

      ! GMBE PIE-based variables (new approach)
      integer :: n_primaries  !! Number of primary polymers
      integer(int64) :: n_primaries_i64  !! For binomial calculation
      integer, allocatable :: pie_atom_sets(:, :)  !! Unique atom sets (max_atoms, n_pie_terms)
      integer, allocatable :: pie_coefficients(:)  !! PIE coefficient for each term
      integer(int64) :: n_pie_terms  !! Number of unique PIE terms

      ! Generate fragments
      if (world_comm%rank() == 0) then
         if (allow_overlapping_fragments) then
            ! GMBE mode: PIE-based inclusion-exclusion
            ! GMBE(1): primaries are monomers
            ! GMBE(N): primaries are N-mers (e.g., dimers for N=2)

            ! Generate primaries
            if (max_level == 1) then
               ! GMBE(1): primaries are base monomers
               n_primaries = sys_geom%n_monomers
               allocate (polymers(n_primaries, 1))
               do i = 1, n_primaries
                  polymers(i, 1) = i
               end do
            else
               ! GMBE(N): primaries are all C(M, N) N-tuples
               n_primaries_i64 = binomial(sys_geom%n_monomers, max_level)
               n_primaries = int(n_primaries_i64)
               allocate (monomers(sys_geom%n_monomers))
               allocate (polymers(n_primaries, max_level))
               polymers = 0

               call create_monomer_list(monomers)
               total_fragments = 0_int64
               call combine(monomers, sys_geom%n_monomers, max_level, polymers, total_fragments)
               n_primaries = int(total_fragments)
               deallocate (monomers)

               ! Apply distance-based screening to primaries if cutoffs are provided
               if (max_level > 1) then
                  ! Only screen if primaries are n-mers (not for GMBE(1) where primaries are monomers)
                  total_fragments = int(n_primaries, int64)
                  call apply_distance_screening(polymers, total_fragments, sys_geom, driver_config, max_level)
                  n_primaries = int(total_fragments)
               end if

               ! Sort primaries by size (largest first)
               ! TODO: Currently disabled - see comment in MBE section above
               ! total_fragments = int(n_primaries, int64)
               call sort_fragments_by_size(polymers, total_fragments, max_level)
            end if

            call logger%info("Generated "//to_char(n_primaries)//" primary "//to_char(max_level)//"-mers for GMBE("// &
                             to_char(max_level)//")")

            ! Use DFS to enumerate PIE terms with coefficients
            call gmbe_enumerate_pie_terms(sys_geom, polymers, n_primaries, max_level, max_intersection_level, &
                                          pie_atom_sets, pie_coefficients, n_pie_terms)

            call logger%info("GMBE PIE enumeration complete: "//to_char(n_pie_terms)//" unique subsystems to evaluate")

            ! For now: total_fragments = n_pie_terms (each PIE term is a subsystem to evaluate)
            total_fragments = n_pie_terms
         else
            ! Standard MBE mode
            ! Calculate expected number of fragments
            n_expected_frags = get_nfrags(sys_geom%n_monomers, max_level)
            n_rows = n_expected_frags

            ! Allocate monomer list and polymers array
            allocate (monomers(sys_geom%n_monomers))
            allocate (polymers(n_rows, max_level))
            polymers = 0

            ! Create monomer list [1, 2, 3, ..., n_monomers]
            call create_monomer_list(monomers)

            ! Generate all fragments (includes monomers in polymers array)
            total_fragments = 0_int64

            ! First add monomers
            do i = 1, sys_geom%n_monomers
               total_fragments = total_fragments + 1_int64
               polymers(total_fragments, 1) = i
            end do

            ! Then add n-mers for n >= 2
            call generate_fragment_list(monomers, max_level, polymers, total_fragments)

            deallocate (monomers)

            ! Apply distance-based screening if cutoffs are provided
            call apply_distance_screening(polymers, total_fragments, sys_geom, driver_config, max_level)

            ! Sort fragments by size (largest first) for better load balancing
            ! TODO: Currently disabled - MBE assembly is now order-independent (uses nested loops),
            ! but sorting still causes "Subset not found" errors in real validation cases.
            ! Unit tests pass with arbitrary order, so there may be an issue with the hash table
            ! or fragment generation in production code. Needs investigation.
            call sort_fragments_by_size(polymers, total_fragments, max_level)

            call logger%info("Generated fragments:")
            call logger%info("  Total fragments: "//to_char(total_fragments))
            call logger%info("  Max level: "//to_char(max_level))
         end if
      end if

      ! Broadcast total_fragments to all ranks
      call bcast(world_comm, total_fragments, 1, 0)

      ! Determine node leaders
      global_node_rank = -1
      if (node_comm%rank() == 0) global_node_rank = world_comm%rank()

      allocate (all_node_leader_ranks(world_comm%size()))
      call allgather(world_comm, global_node_rank, all_node_leader_ranks)

      num_nodes = count(all_node_leader_ranks /= -1)

      if (world_comm%rank() == 0) then
         call logger%info("Running with "//to_char(num_nodes)//" node(s)")
      end if

      allocate (node_leader_ranks(num_nodes))
      i = 0
      do j = 1, world_comm%size()
         if (all_node_leader_ranks(j) /= -1) then
            i = i + 1
            node_leader_ranks(i) = all_node_leader_ranks(j)
         end if
      end do
      deallocate (all_node_leader_ranks)

      ! Execute appropriate role
      if (world_comm%size() == 1) then
         ! Single rank: process fragments serially
         call logger%info("Running in serial mode (single MPI rank)")
         if (allow_overlapping_fragments) then
            ! GMBE serial processing with PIE coefficients
          call serial_gmbe_pie_processor(pie_atom_sets, pie_coefficients, n_pie_terms, sys_geom, method, calc_type, bonds)
         else
            ! Standard MBE serial processing
            call serial_fragment_processor(total_fragments, polymers, max_level, sys_geom, method, calc_type, bonds)
         end if
      else if (world_comm%leader() .and. node_comm%leader()) then
         ! Global coordinator (rank 0, node leader on node 0)
         call omp_set_num_threads(omp_get_max_threads())
         call logger%verbose("Rank 0: Acting as global coordinator")
         if (allow_overlapping_fragments) then
            ! GMBE MPI processing - PIE-based approach
            call gmbe_pie_coordinator(world_comm, node_comm, pie_atom_sets, pie_coefficients, n_pie_terms, &
                                      node_leader_ranks, num_nodes, sys_geom, method, calc_type, bonds)
         else
            ! Standard MBE MPI processing
            call global_coordinator(world_comm, node_comm, total_fragments, polymers, max_level, &
                                    node_leader_ranks, num_nodes, sys_geom, calc_type, bonds)
         end if
      else if (node_comm%leader()) then
         ! Node coordinator (node leader on other nodes)
         call logger%verbose("Rank "//to_char(world_comm%rank())//": Acting as node coordinator")
         ! Node coordinator works for both MBE and GMBE (receives fragments from global coordinator)
         call node_coordinator(world_comm, node_comm, calc_type)
      else
         ! Worker
         call omp_set_num_threads(1)
         call logger%verbose("Rank "//to_char(world_comm%rank())//": Acting as worker")
         ! Worker processes work for both MBE and GMBE (fragment_type distinguishes them)
         call node_worker(world_comm, node_comm, sys_geom, method, calc_type, bonds)
      end if

      ! Cleanup
      if (world_comm%rank() == 0) then
         if (allocated(polymers)) deallocate (polymers)
         if (allocated(node_leader_ranks)) deallocate (node_leader_ranks)
         if (allocated(intersections)) deallocate (intersections)
         if (allocated(intersection_sets)) deallocate (intersection_sets)
         if (allocated(intersection_levels)) deallocate (intersection_levels)
         if (allocated(pie_atom_sets)) deallocate (pie_atom_sets)
         if (allocated(pie_coefficients)) deallocate (pie_coefficients)
      end if

   end subroutine run_fragmented_calculation

   subroutine run_multi_molecule_calculations(world_comm, node_comm, mqc_config)
      !! Run calculations for multiple molecules with MPI parallelization
      !! Each molecule is independent, so assign one molecule per rank
      use mqc_config_parser, only: mqc_config_t
      use mqc_config_adapter, only: config_to_system_geometry
      use mqc_error, only: error_t
      use mqc_io_helpers, only: set_molecule_suffix, get_output_json_filename
      use mqc_json, only: merge_multi_molecule_json

      type(comm_t), intent(in) :: world_comm
      type(comm_t), intent(in) :: node_comm
      type(mqc_config_t), intent(in) :: mqc_config

      type(driver_config_t) :: config
      type(system_geometry_t) :: sys_geom
      type(comm_t) :: mol_comm, mol_node_comm
      type(error_t) :: error
      integer :: imol, my_rank, num_ranks, color
      integer :: molecules_processed
      character(len=:), allocatable :: mol_name
      logical :: has_fragmented_molecules
      character(len=256), allocatable :: individual_json_files(:)

      my_rank = world_comm%rank()
      num_ranks = world_comm%size()

      ! Allocate array to track individual JSON files for merging
      allocate (individual_json_files(mqc_config%nmol))

      ! Check if any molecules have fragments (nlevel > 0)
      has_fragmented_molecules = .false.
      do imol = 1, mqc_config%nmol
         if (mqc_config%molecules(imol)%nfrag > 0) then
            has_fragmented_molecules = .true.
            exit
         end if
      end do

      if (my_rank == 0) then
         call logger%info(" ")
         call logger%info("============================================")
         call logger%info("Multi-molecule mode: "//to_char(mqc_config%nmol)//" molecules")
         call logger%info("MPI ranks: "//to_char(num_ranks))
         if (has_fragmented_molecules) then
            call logger%info("Mode: Sequential execution (fragmented molecules detected)")
            call logger%info("  Each molecule will use all "//to_char(num_ranks)//" rank(s) for its calculation")
         else if (num_ranks == 1) then
            call logger%info("Mode: Sequential execution (single rank)")
         else if (num_ranks > mqc_config%nmol) then
            call logger%info("Mode: Parallel execution (one molecule per rank)")
            call logger%info("Note: More ranks than molecules - ranks "//to_char(mqc_config%nmol)// &
                             " to "//to_char(num_ranks - 1)//" will be idle")
         else
            call logger%info("Mode: Parallel execution (one molecule per rank)")
         end if
         call logger%info("============================================")
         call logger%info(" ")
      end if

      ! Determine execution mode:
      ! 1. Sequential: Single rank OR fragmented molecules (each molecule needs all ranks)
      ! 2. Parallel: Multiple ranks AND unfragmented molecules (distribute molecules across ranks)
      molecules_processed = 0

      if (num_ranks == 1 .or. has_fragmented_molecules) then
         ! Sequential mode: process all molecules one after another
         ! Each molecule uses all available ranks for its calculation
         do imol = 1, mqc_config%nmol
            ! Determine molecule name for logging
            if (allocated(mqc_config%molecules(imol)%name)) then
               mol_name = mqc_config%molecules(imol)%name
            else
               mol_name = "molecule_"//to_char(imol)
            end if

            if (my_rank == 0) then
               call logger%info(" ")
               call logger%info("--------------------------------------------")
               call logger%info("Processing molecule "//to_char(imol)//"/"//to_char(mqc_config%nmol)//": "//mol_name)
               call logger%info("--------------------------------------------")
            end if

            ! Convert to driver configuration for this molecule
            call config_to_driver(mqc_config, config, molecule_index=imol)

            ! Convert geometry for this molecule
            call config_to_system_geometry(mqc_config, sys_geom, error, molecule_index=imol)
            if (error%has_error()) then
               call error%add_context("mqc_driver:run_multi_molecule_calculation")
               if (my_rank == 0) then
                  call logger%error("Error converting geometry for "//mol_name//": "//error%get_full_trace())
               end if
               call abort_comm(world_comm, 1)
            end if

            ! Set output filename suffix for this molecule
            call set_molecule_suffix("_"//trim(mol_name))

            ! Run calculation for this molecule
            call run_calculation(world_comm, node_comm, config, sys_geom, mqc_config%molecules(imol)%bonds)

            ! Track the JSON filename for later merging
            individual_json_files(imol) = get_output_json_filename()

            ! Clean up for this molecule
            call sys_geom%destroy()

            if (my_rank == 0) then
               call logger%info("Completed molecule "//to_char(imol)//"/"//to_char(mqc_config%nmol)//": "//mol_name)
            end if
            molecules_processed = molecules_processed + 1
         end do
      else
         ! Multiple ranks: distribute molecules across ranks in round-robin fashion
         molecules_processed = 0
         do imol = 1, mqc_config%nmol
            ! This rank processes molecules where (imol - 1) mod num_ranks == my_rank
            if (mod(imol - 1, num_ranks) == my_rank) then
               ! Determine molecule name for logging
               if (allocated(mqc_config%molecules(imol)%name)) then
                  mol_name = mqc_config%molecules(imol)%name
               else
                  mol_name = "molecule_"//to_char(imol)
               end if

               call logger%info(" ")
               call logger%info("--------------------------------------------")
               call logger%info("Rank "//to_char(my_rank)//": Processing molecule "//to_char(imol)// &
                                "/"//to_char(mqc_config%nmol)//": "//mol_name)
               call logger%info("--------------------------------------------")

               ! Convert to driver configuration for this molecule
               call config_to_driver(mqc_config, config, molecule_index=imol)

               ! Convert geometry for this molecule
               call config_to_system_geometry(mqc_config, sys_geom, error, molecule_index=imol)
               if (error%has_error()) then
                  call error%add_context("mqc_driver:run_multi_molecule_calculation")
  call logger%error("Rank "//to_char(my_rank)//": Error converting geometry for "//mol_name//": "//error%get_full_trace())
                  call abort_comm(world_comm, 1)
               end if

               ! Set output filename suffix for this molecule
               call set_molecule_suffix("_"//trim(mol_name))

               ! Run calculation for this molecule
               call run_calculation(world_comm, node_comm, config, sys_geom, mqc_config%molecules(imol)%bonds)

               ! Track the JSON filename for later merging
               individual_json_files(imol) = get_output_json_filename()

               ! Clean up for this molecule
               call sys_geom%destroy()

               call logger%info("Rank "//to_char(my_rank)//": Completed molecule "//to_char(imol)// &
                                "/"//to_char(mqc_config%nmol)//": "//mol_name)

               molecules_processed = molecules_processed + 1
            end if
         end do

         if (molecules_processed == 0) then
            ! Idle rank - no molecules assigned
            call logger%verbose("Rank "//to_char(my_rank)//": No molecules assigned (idle)")
         end if
      end if

      ! Synchronize all ranks
      call world_comm%barrier()

      ! In parallel execution, rank 0 needs to reconstruct all JSON filenames for merging
      ! since each rank only populated its own entry
      if (my_rank == 0 .and. num_ranks > 1 .and. .not. has_fragmented_molecules) then
         ! Rank 0 constructs filenames for all molecules
         do imol = 1, mqc_config%nmol
            ! Get molecule name
            if (allocated(mqc_config%molecules(imol)%name)) then
               mol_name = mqc_config%molecules(imol)%name
            else
               mol_name = "molecule_"//to_char(imol)
            end if
            ! Construct JSON filename pattern: output_<basename>_<molname>.json
            ! This mirrors what get_output_json_filename() returns after set_molecule_suffix()
            call set_molecule_suffix("_"//trim(mol_name))
            individual_json_files(imol) = get_output_json_filename()
         end do
      end if

      ! Merge individual JSON files into one combined file (rank 0 only)
      if (my_rank == 0) then
         call merge_multi_molecule_json(individual_json_files, mqc_config%nmol)
      end if

      if (my_rank == 0) then
         call logger%info(" ")
         call logger%info("============================================")
         call logger%info("All "//to_char(mqc_config%nmol)//" molecules completed")
         if (has_fragmented_molecules) then
            call logger%info("Execution: Sequential (each molecule used all ranks)")
         else if (num_ranks == 1) then
            call logger%info("Execution: Sequential (single rank)")
         else if (num_ranks > mqc_config%nmol) then
           call logger%info("Execution: Parallel (active ranks: "//to_char(mqc_config%nmol)//"/"//to_char(num_ranks)//")")
         else
            call logger%info("Execution: Parallel (all ranks active)")
         end if
         call logger%info("============================================")
      end if

   end subroutine run_multi_molecule_calculations

end module mqc_driver
