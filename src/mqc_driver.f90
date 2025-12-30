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
                                                   serial_fragment_processor, do_fragment_work
   use mqc_frag_utils, only: get_nfrags, create_monomer_list, generate_fragment_list, generate_intersections
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, &
                                    build_fragment_from_indices, build_fragment_from_atom_list
   use mqc_config_adapter, only: driver_config_t, config_to_driver, config_to_system_geometry
   use mqc_method_types, only: method_type_to_string
   use mqc_calc_types, only: calc_type_to_string, CALC_TYPE_GRADIENT
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

   subroutine run_calculation(world_comm, node_comm, config, sys_geom, bonds)
      !! Main calculation dispatcher - routes to fragmented or unfragmented calculation
      !!
      !! Determines calculation type based on nlevel and dispatches to appropriate
      !! calculation routine with proper MPI setup and validation.
      type(comm_t), intent(in) :: world_comm  !! Global MPI communicator
      type(comm_t), intent(in) :: node_comm   !! Node-local MPI communicator
      type(driver_config_t), intent(in) :: config  !! Driver configuration
      type(system_geometry_t), intent(in) :: sys_geom  !! System geometry and fragment info
      type(bond_t), intent(in), optional :: bonds(:)  !! Bond connectivity information

      ! Local variables
      integer :: max_level   !! Maximum fragment level (nlevel from config)
      integer :: matrix_size  !! Size of gradient matrix (natoms*3), tmp
      integer :: i  !! Loop counter
      logical :: has_broken_bonds  !! Flag for broken bonds (H-capping will occur)

      ! Set max_level from config
      ! TODO JORGE: change to max fragmentation level
      max_level = config%nlevel

      ! Set matrix_size based on atoms per monomer (natoms * 3 for gradient)
      ! TODO JORGE: this is temporary, until we define a result_struct that will initialize itself
      matrix_size = sys_geom%atoms_per_monomer*3

      if (world_comm%rank() == 0) then
         call logger%info("============================================")
         call logger%info("Loaded geometry:")
         call logger%info("  Total monomers: "//to_char(sys_geom%n_monomers))
         call logger%info("  Atoms per monomer: "//to_char(sys_geom%atoms_per_monomer))
         call logger%info("  Total atoms: "//to_char(sys_geom%total_atoms))
         call logger%info("  Fragment level: "//to_char(max_level))
         call logger%info("  Matrix size (natoms*3): "//to_char(matrix_size))
         call logger%info("============================================")
      end if

      ! Warn if overlapping fragments flag is set but nlevel=0
      if (config%allow_overlapping_fragments .and. max_level == 0) then
         if (world_comm%rank() == 0) then
            call logger%warning("allow_overlapping_fragments is set to true, but nlevel=0")
            call logger%warning("Running unfragmented calculation - overlapping fragments flag will be ignored")
         end if
      end if

      ! Validate GMBE (overlapping fragments) settings
      if (config%allow_overlapping_fragments .and. max_level > 1) then
         if (world_comm%rank() == 0) then
            call logger%error(" ")
            call logger%error("ERROR: Overlapping fragments (GMBE) only supported with nlevel=1")
            call logger%error(" ")
            call logger%error("Current settings:")
            call logger%error("  nlevel = "//to_char(max_level))
            call logger%error("  allow_overlapping_fragments = true")
            call logger%error(" ")
            call logger%error("GMBE (Generalized Many-Body Expansion) with overlapping fragments")
            call logger%error("is currently restricted to monomer-level calculations (nlevel=1).")
            call logger%error(" ")
            call logger%error("Solution: Set nlevel = 1 in your input file")
            call logger%error(" ")
         end if
         call abort_comm(world_comm, 1)
      end if

      ! Validate gradient calculations with overlapping fragments
      if (config%allow_overlapping_fragments .and. max_level > 0 .and. config%calc_type == CALC_TYPE_GRADIENT) then
         if (world_comm%rank() == 0) then
            call logger%error(" ")
            call logger%error("ERROR: Gradient calculations with overlapping fragments are not supported")
            call logger%error(" ")
            call logger%error("Current settings:")
            call logger%error("  nlevel = "//to_char(max_level))
            call logger%error("  allow_overlapping_fragments = true")
            call logger%error("  calc_type = Gradient")
            call logger%error(" ")
            call logger%error("GMBE gradients require implementing gradient intersection corrections,")
            call logger%error("which is not currently supported.")
            call logger%error(" ")
            call logger%error("Solutions:")
            call logger%error("  1. Use calc_type = Energy instead of Gradient")
            call logger%error("  2. Set allow_overlapping_fragments = false")
            call logger%error("  3. Use unfragmented calculation (nlevel = 0)")
            call logger%error(" ")
         end if
         call abort_comm(world_comm, 1)
      end if

      ! Validate gradient calculations with H-capping
      ! Gradients with H-caps are not physically meaningful since H-caps are artificial atoms
      if (max_level > 0 .and. config%calc_type == CALC_TYPE_GRADIENT) then
         if (present(bonds)) then
            has_broken_bonds = .false.
            do i = 1, size(bonds)
               if (bonds(i)%is_broken) then
                  has_broken_bonds = .true.
                  exit
               end if
            end do

            if (has_broken_bonds) then
               if (world_comm%rank() == 0) then
                  call logger%error(" ")
                  call logger%error("ERROR: Gradient calculations with hydrogen capping are not supported")
                  call logger%error(" ")
                  call logger%error("Hydrogen caps are artificial atoms added where bonds are broken.")
                  call logger%error("Their gradients are not physically meaningful and would give")
                  call logger%error("incorrect results for geometry optimization or force calculations.")
                  call logger%error(" ")
                  call logger%error("Solutions:")
                  call logger%error("  1. Use calc_type = Energy instead of Gradient")
                  call logger%error("  2. Remove broken bonds from connectivity (no H-capping)")
                  call logger%error("  3. Use unfragmented calculation (nlevel = 0)")
                  call logger%error(" ")
               end if
               call abort_comm(world_comm, 1)
            end if
         end if
      end if

      if (max_level == 0) then
         call omp_set_num_threads(1)
         call run_unfragmented_calculation(world_comm, sys_geom, config%method, config%calc_type, bonds)
      else
         call run_fragmented_calculation(world_comm, node_comm, config%method, config%calc_type, sys_geom, max_level, &
                                         matrix_size, config%allow_overlapping_fragments, bonds)
      end if

   end subroutine run_calculation

   subroutine run_unfragmented_calculation(world_comm, sys_geom, method, calc_type, bonds)
      !! Handle unfragmented calculation (nlevel=0)
      !!
      !! For single-molecule mode: Only rank 0 runs (validates single rank)
      !! For multi-molecule mode: ALL ranks can run (each with their own molecule)
      type(comm_t), intent(in) :: world_comm  !! Global MPI communicator
      type(system_geometry_t), intent(in) :: sys_geom  !! Complete system geometry
      integer(int32), intent(in) :: method  !! Quantum chemistry method
      integer(int32), intent(in) :: calc_type  !! Calculation type
      type(bond_t), intent(in), optional :: bonds(:)  !! Bond connectivity information

      ! Check if this is multi-molecule mode or single-molecule mode
      ! In multi-molecule mode, each rank processes its own molecule
      ! In single-molecule mode, only rank 0 should work
      if (world_comm%size() == 1 .or. world_comm%rank() == 0) then
         ! Either single-rank calculation, or rank 0 in multi-rank setup
         call logger%info(" ")
         call logger%info("Running unfragmented calculation")
         call logger%info("  Calculation type: "//calc_type_to_string(calc_type))
         call logger%info(" ")
         call unfragmented_calculation(sys_geom, method, calc_type, bonds)
      else if (sys_geom%total_atoms > 0) then
         ! Multi-molecule mode: non-zero rank with a molecule
         call logger%verbose("Rank "//to_char(world_comm%rank())//": Running unfragmented calculation")
         call unfragmented_calculation(sys_geom, method, calc_type, bonds)
      end if

   end subroutine run_unfragmented_calculation

   subroutine run_fragmented_calculation(world_comm, node_comm, method, calc_type, sys_geom, max_level, matrix_size, &
                                         allow_overlapping_fragments, bonds)
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
      integer, intent(in) :: matrix_size  !! Size of gradient matrix (natoms*3)
      logical, intent(in) :: allow_overlapping_fragments  !! Use GMBE for overlapping fragments
      type(bond_t), intent(in), optional :: bonds(:)  !! Bond connectivity information

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

      ! GMBE-specific variables
      integer, allocatable :: intersections(:, :)  !! Intersection atom lists (max_atoms, n_intersections)
      integer, allocatable :: intersection_pairs(:, :)  !! Monomer pairs for each intersection (2, n_intersections)
      integer :: n_intersections, n_monomers  !! Counts for GMBE

      ! Generate fragments
      if (world_comm%rank() == 0) then
         if (allow_overlapping_fragments) then
            ! GMBE mode: generate monomers + intersections
            n_monomers = sys_geom%n_monomers
            allocate (monomers(n_monomers))
            allocate (polymers(n_monomers, 1))  ! Only monomers
            polymers = 0

            ! Create monomer list
            call create_monomer_list(monomers)

            ! Add monomers to polymers
            do i = 1, n_monomers
               polymers(i, 1) = i
            end do

            ! Generate intersections
            call generate_intersections(sys_geom, monomers, polymers, n_monomers, &
                                        intersections, intersection_pairs, n_intersections)

            ! Total fragments = monomers + intersections
            total_fragments = int(n_monomers, int64) + int(n_intersections, int64)

            call logger%info("Generated GMBE fragments:")
            call logger%info("  Monomers: "//to_char(n_monomers))
            call logger%info("  Intersections: "//to_char(n_intersections))
            call logger%info("  Total fragments: "//to_char(total_fragments))

            deallocate (monomers)
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
            ! GMBE serial processing
            call serial_gmbe_processor(n_monomers, polymers, intersections, intersection_pairs, n_intersections, &
                                       sys_geom, method, calc_type, bonds)
         else
            ! Standard MBE serial processing
            call serial_fragment_processor(total_fragments, polymers, max_level, sys_geom, method, calc_type, bonds)
         end if
      else if (world_comm%leader() .and. node_comm%leader()) then
         ! Global coordinator (rank 0, node leader on node 0)
         call omp_set_num_threads(omp_get_max_threads())
         call logger%verbose("Rank 0: Acting as global coordinator")
         call global_coordinator(world_comm, node_comm, total_fragments, polymers, max_level, &
                                 node_leader_ranks, num_nodes, sys_geom, calc_type)
      else if (node_comm%leader()) then
         ! Node coordinator (node leader on other nodes)
         call logger%verbose("Rank "//to_char(world_comm%rank())//": Acting as node coordinator")
         call node_coordinator(world_comm, node_comm, calc_type)
      else
         ! Worker
         call omp_set_num_threads(1)
         call logger%verbose("Rank "//to_char(world_comm%rank())//": Acting as worker")
         call node_worker(world_comm, node_comm, sys_geom, method, calc_type, bonds)
      end if

      ! Cleanup
      if (world_comm%rank() == 0) then
         if (allocated(polymers)) deallocate (polymers)
         if (allocated(node_leader_ranks)) deallocate (node_leader_ranks)
         if (allocated(intersections)) deallocate (intersections)
         if (allocated(intersection_pairs)) deallocate (intersection_pairs)
      end if

   end subroutine run_fragmented_calculation

   subroutine serial_gmbe_processor(n_monomers, polymers, intersections, intersection_pairs, n_intersections, &
                                    sys_geom, method, calc_type, bonds)
      !! Serial GMBE processor: builds all fragments (monomers + intersections) and computes GMBE energy
      integer, intent(in) :: n_monomers  !! Number of monomers
      integer, intent(in) :: polymers(:, :)  !! Monomer indices (n_monomers, 1)
      integer, intent(in) :: intersections(:, :)  !! Intersection atom lists
      integer, intent(in) :: intersection_pairs(:, :)  !! Pairs that created intersections
      integer, intent(in) :: n_intersections  !! Number of intersection fragments
      type(system_geometry_t), intent(in) :: sys_geom
      integer(int32), intent(in) :: method, calc_type
      type(bond_t), intent(in), optional :: bonds(:)

      type(calculation_result_t), allocatable :: all_results(:)
      type(physical_fragment_t) :: phys_frag
      integer :: i, intersection_size
      integer, allocatable :: monomer_idx(:), monomer_indices(:)
      real(dp) :: total_energy

      ! Total results = monomers + intersections
      allocate (all_results(n_monomers + n_intersections))

      ! Build and calculate monomers
      call logger%info("Processing "//to_char(n_monomers)//" monomer fragments...")
      do i = 1, n_monomers
         allocate (monomer_idx(1))
         monomer_idx(1) = polymers(i, 1)
         ! Use build_fragment_from_indices for monomers to preserve fragment charges/multiplicities
         call build_fragment_from_indices(sys_geom, monomer_idx, phys_frag, bonds)
         call do_fragment_work(i, all_results(i), method, phys_frag, calc_type)
         call phys_frag%destroy()
         deallocate (monomer_idx)
      end do

      ! Build and calculate intersections
      if (n_intersections > 0) then
         call logger%info("Processing "//to_char(n_intersections)//" intersection fragments...")
         do i = 1, n_intersections
            ! Find size of this intersection (first zero marks end)
            intersection_size = 0
            do while (intersection_size < size(intersections, 1))
               if (intersections(intersection_size + 1, i) == 0) exit
               intersection_size = intersection_size + 1
            end do

            ! Build intersection fragment from atom list
            call build_fragment_from_atom_list(sys_geom, intersections(1:intersection_size, i), &
                                               intersection_size, phys_frag, bonds)
            call do_fragment_work(n_monomers + i, all_results(n_monomers + i), method, phys_frag, calc_type)
            call phys_frag%destroy()
         end do
      end if

      ! Compute GMBE energy
      allocate (monomer_indices(n_monomers))
      do i = 1, n_monomers
         monomer_indices(i) = polymers(i, 1)
      end do

      call compute_gmbe_energy(monomer_indices, n_monomers, all_results(1:n_monomers), &
                               n_intersections, all_results(n_monomers + 1:n_monomers + n_intersections), &
                               intersection_pairs, total_energy)

      call logger%info(" ")
      call logger%info("GMBE calculation completed successfully")
      call logger%info("Final GMBE energy: "//to_char(total_energy)//" Hartree")
      call logger%info(" ")

      deallocate (all_results, monomer_indices)

   end subroutine serial_gmbe_processor

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
               if (my_rank == 0) then
                  call logger%error("Error converting geometry for "//mol_name//": "//error%get_message())
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
         ! Multiple ranks: distribute molecules across ranks (one per rank)
         if (my_rank < mqc_config%nmol) then
            imol = my_rank + 1
            ! This rank has a molecule to process
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
     call logger%error("Rank "//to_char(my_rank)//": Error converting geometry for "//mol_name//": "//error%get_message())
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

            molecules_processed = 1
         else
            ! Idle rank - no molecule assigned
            call logger%verbose("Rank "//to_char(my_rank)//": No molecule assigned (idle)")
         end if
      end if

      ! Synchronize all ranks
      call world_comm%barrier()

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
