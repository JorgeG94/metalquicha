!! Main calculation driver module for metalquicha
module mqc_driver
   !! Handles both fragmented (many-body expansion) and unfragmented calculations
   !! with MPI parallelization and node-based work distribution.
   use pic_types, only: int32, int64, dp, default_int
   use pic_mpi_lib, only: comm_t, abort_comm, bcast, allgather
   use mqc_resources, only: resources_t
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   use mqc_method_types, only: needs_serial_execution
   use mqc_mbe_fragment_distribution_scheme, only: unfragmented_calculation, distributed_unfragmented_hessian
   use mqc_many_body_expansion, only: many_body_expansion_t, mbe_context_t, gmbe_context_t, &
                                      fmo_context_t
   use mqc_method_config, only: method_config_t
   ! GMBE functions are now called via type-bound procedures in gmbe_context_t
   use mqc_validate, only: validate_system, validate_terms
   use mqc_fraglist, only: fraglist_t
   use mqc_frag_utils, only: get_nfrags, create_monomer_list, generate_fragment_list, generate_intersections, &
                             gmbe_enumerate_pie_terms, binomial, combine, apply_distance_screening, &
                             sort_fragments_by_size, generate_mbe_term_list
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, &
                                    build_fragment_from_indices, build_fragment_from_atom_list
   use mqc_config_adapter, only: driver_config_t, config_to_driver, config_to_system_geometry, &
                                 check_counterpoise_support
   use mqc_method_types, only: method_type_to_string
   use mqc_calc_types, only: calc_type_to_string, CALC_TYPE_ENERGY, CALC_TYPE_GRADIENT, &
                             CALC_TYPE_OPTIMIZE, &
                             CALC_TYPE_HESSIAN, CALC_TYPE_MAKEFP
   use mqc_config_types, only: bond_t, mqc_config_t
   use mqc_mbe, only: compute_gmbe
   use mqc_result_types, only: calculation_result_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_fingerprint, only: calculation_fingerprint
   use mqc_checkpoint, only: checkpoint_t
   use mqc_io_helpers, only: set_molecule_suffix, get_output_json_filename
   use mqc_json, only: merge_multi_molecule_json
   use mqc_json_output_types, only: json_output_data_t, OUTPUT_MODE_NONE
   use mqc_json_writer, only: write_json_output
   implicit none
   private

   public :: run_calculation  !! Main entry point for all calculations
   public :: run_multi_molecule_calculations  !! Multi-molecule calculation dispatcher

contains

   subroutine run_calculation(resources, config, sys_geom, bonds, result_out, all_ranks_write_json, &
                              supplied_terms, n_supplied_terms, write_output)
      !! Main calculation dispatcher - routes to fragmented or unfragmented calculation
      !!
      !! Determines calculation type based on nlevel and dispatches to appropriate
      !! calculation routine with proper MPI setup and validation.
      !! `result_out` and `write_output` are independent. A caller can take the
      !! energy back, write the files, or both -- the common case for a driven
      !! run being both, so a script can print or branch on a number without
      !! re-reading the file this just wrote. They used to be one switch, which
      !! meant asking for the result silently suppressed the output.
      !!
      !! `supplied_terms` hands over a term list instead of generating one, which
      !! is how a caller applies a criterion this code knows nothing about --
      !! an energy threshold from a previous run's breakdown, the terms a dead
      !! run never reached, a hand-picked set. Distance screening and sorting
      !! are skipped with it: the list is taken as given, in the order given.
      !!
      !! Only rank 0 needs it. The coordinator hands individual fragment
      !! definitions to workers as tasks and they build them from `sys_geom`,
      !! so the list itself never leaves rank 0 -- but every rank must already
      !! have `sys_geom`, which the normal program gets by having every rank
      !! parse the input and a session gets by broadcasting it.
      !!
      !! **The list must be closed under subsets.** An n-body term's delta is
      !! its energy less every proper subset's delta, so keeping a trimer whose
      !! dimers were screened away fails the lookup rather than approximating
      !! anything. `fraglist_t%close_subsets` exists for this.
      use mqc_method_types, only: METHOD_TYPE_EFP2, METHOD_TYPE_SAPT0
      type(resources_t), intent(in) :: resources  !! Resources container (MPI comms, etc.)
      type(driver_config_t), intent(in) :: config  !! Driver configuration
      type(system_geometry_t), intent(in) :: sys_geom  !! System geometry and fragment info
      type(bond_t), intent(in), optional :: bonds(:)  !! Bond connectivity information
      type(calculation_result_t), intent(out), optional :: result_out  !! Optional result output
      logical, intent(in), optional :: all_ranks_write_json  !! If true, all ranks write JSON (for multi-molecule)
      logical, intent(in), optional :: write_output
         !! Write the JSON summary and fragment breakdown. Default true.
         !! `skip_json_output` in the config still overrides it.
      integer, intent(in), optional :: supplied_terms(:, :)
         !! (n_terms, max_level), 1-based monomer indices, zero-padded. Rank 0 only.
      integer(int64), intent(in), optional :: n_supplied_terms

      ! Local variables
      integer :: max_level   !! Maximum fragment level (nlevel from config)
      integer :: i  !! Loop counter
      type(json_output_data_t) :: json_data  !! Cached output data for centralized JSON writing
      logical :: should_write_json  !! Whether this rank should write JSON
      logical :: wants_output       !! Whether the caller asked for files at all

      ! Set max_level from config
      max_level = config%nlevel

      ! Log method-specific settings (rank 0 only)
      if (resources%mpi_comms%world_comm%rank() == 0) then
         call config%method_config%log_settings()
      end if

      if (resources%mpi_comms%world_comm%rank() == 0 .and. max_level > 0) then
         call logger%info("============================================")
         call logger%info("Loaded geometry:")
         call logger%info("  Total monomers: "//to_char(sys_geom%n_monomers))
         ! Zero is the "variable-sized fragments" sentinel, not a count -- say so
         ! rather than printing "0 atoms per monomer", which reads as nonsense.
         if (sys_geom%atoms_per_monomer > 0) then
            call logger%info("  Atoms per monomer: "//to_char(sys_geom%atoms_per_monomer))
         else
            call logger%info("  Atoms per monomer: variable")
         end if
         call logger%info("  Fragment level: "//to_char(max_level))
         call logger%info("  Total atoms: "//to_char(sys_geom%total_atoms))
         call logger%info("============================================")
      end if

      ! An optimization is a loop over calculations and is driven from above
      ! this routine, so reaching here with one means a caller took a path that
      ! does not know about it -- the multi-molecule loop, a session, the C
      ! API. Refused by name rather than left to fall through: `calc_type`
      ! reaches the method layer as an unhandled value, and what that produced
      ! was a backtrace out of `unfragmented_calculation` rather than anything
      ! a user could act on.
      if (config%calc_type == CALC_TYPE_OPTIMIZE) then
         if (resources%mpi_comms%world_comm%rank() == 0) then
            call logger%error('driver "Optimize" is not available on this path.')
            call logger%error("  It is driven above run_calculation, so it works for a "// &
                              "single molecule")
            call logger%error("  from an input deck, and not yet from a multi-molecule "// &
                              "deck, a session")
            call logger%error("  or the C API.")
         end if
         call abort_comm(resources%mpi_comms%world_comm, 1)
      end if

      ! Resolved before the branches below rather than beside the JSON write at
      ! the end, because those branches write their own summary and return
      ! without ever reaching it. Left where it was, a SAPT run asked for no
      ! files wrote one anyway.
      wants_output = .true.
      if (present(write_output)) wants_output = write_output

      ! MAKEFP writes a file and returns no energy, so it takes neither the
      ! fragmented nor the unfragmented path: a fragment potential is built from
      ! the whole system by definition, and there is no result_t to fill.
      if (config%calc_type == CALC_TYPE_MAKEFP) then
         call run_makefp(config, sys_geom, resources%mpi_comms%world_comm%rank(), &
                         result_out)
         return
      end if

      ! EFP takes neither path either, for the opposite reason to MAKEFP: the
      ! fragments already carry their wavefunctions, so there is no SCF to fragment
      ! and nothing for a many-body expansion to expand. What is evaluated is the
      ! interaction between potentials that already exist.
      if (config%method_config%method_type == METHOD_TYPE_EFP2) then
         call run_efp(config, sys_geom, resources%mpi_comms%world_comm%rank(), &
                      wants_output, result_out)
         return
      end if

      ! SAPT takes neither path either: it returns the interaction between two
      ! monomers rather than the energy of one system, so there is nothing for a
      ! many-body expansion to expand and no single wavefunction to fragment.
      if (config%method_config%method_type == METHOD_TYPE_SAPT0) then
         call run_sapt(config, sys_geom, resources%mpi_comms%world_comm%rank(), &
                       wants_output, result_out)
         return
      end if

      ! Warn if overlapping fragments flag is set but nlevel=0
      if (config%allow_overlapping_fragments .and. max_level == 0) then
         if (resources%mpi_comms%world_comm%rank() == 0) then
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
         ! One thread, for the methods that need it rather than for all of them.
         !
         ! tblite is the reason this pin exists: run threaded it corrupts a
         ! result instead of failing, so xTB is clamped here and stays clamped.
         ! The libcint Hartree-Fock path is a different case -- libcint keeps no
         ! state across calls, and its Fock build threads its own quartet loop --
         ! so clamping it does not make it safer, only slower, and on this box
         ! by a factor of the core count.
         !
         ! The clamp is deliberately not restored afterwards. `omp_set_num_threads(1)`
         ! makes `omp_get_max_threads()` report 1, so the value has to be saved
         ! before the call to be recoverable; an unfragmented run ends here, and
         ! the one path that does need the count back saves it first. See
         ! mqc_serial_fragment_processor.
         if (needs_serial_execution(config%method_config%method_type)) then
            call omp_set_num_threads(1)
         end if
         if (present(result_out)) then
            call run_unfragmented_calculation(resources%mpi_comms%world_comm, sys_geom, config, &
                                              result_out, json_data=json_data)
         else
            call run_unfragmented_calculation(resources%mpi_comms%world_comm, sys_geom, config, &
                                              json_data=json_data)
         end if
      else
         ! json_data is collected whether or not it will be written, because it
         ! is also where a fragmented result comes from -- the expansion has no
         ! other route back to a calculation_result_t.
         call run_fragmented_calculation(resources, config, sys_geom, bonds, json_data, &
                                         supplied_terms=supplied_terms, &
                                         n_supplied_terms=n_supplied_terms)
         if (present(result_out)) call result_from_json(json_data, result_out)
      end if

      ! Stamped whether or not it is written, because the fingerprint is part
      ! of the result rather than part of the reporting -- a caller taking the
      ! energy back without files still needs to know what produced it.
      json_data%fingerprint = calculation_fingerprint(sys_geom, config%method_config, &
                                                      config%calc_type)

      ! Centralized JSON output (rank 0 only by default, or all ranks if all_ranks_write_json is set)
      if (wants_output) then
         ! Check if JSON output should be skipped
         if (config%skip_json_output) then
            if (resources%mpi_comms%world_comm%rank() == 0) then
               call logger%info("Skipping JSON output (skip_json_output = true)")
            end if
         else
            ! Determine if this rank should write JSON
            should_write_json = (resources%mpi_comms%world_comm%rank() == 0)
            if (present(all_ranks_write_json)) then
               if (all_ranks_write_json) should_write_json = .true.
            end if

            if (should_write_json) then
               if (json_data%output_mode /= OUTPUT_MODE_NONE) then
                  json_data%fragment_breakdown = config%fragment_breakdown
                  call write_json_output(json_data)
                  call json_data%destroy()
               end if
            end if
         end if
      end if

   end subroutine run_calculation

   subroutine result_from_json(json_data, result_out)
      !! Take a fragmented run's headline numbers back to the caller
      !!
      !! The expansion assembles into `json_output_data_t` and has no other
      !! route to a `calculation_result_t`, so this reads across rather than
      !! recomputing. Only what a driving script would branch on: the total,
      !! and the gradient and dipole when the run produced them. The
      !! per-fragment detail stays in the CSV, where twenty million rows
      !! belong.
      type(json_output_data_t), intent(in) :: json_data
      type(calculation_result_t), intent(inout) :: result_out

      ! energy_t computes its total from components; an expansion total has
      ! no correlation breakdown to put in the others, so it lands in scf --
      ! which is what `total()` will then report.
      result_out%energy%scf = json_data%total_energy
      result_out%has_energy = json_data%has_energy
      if (allocated(json_data%gradient)) then
         result_out%gradient = json_data%gradient
         result_out%has_gradient = .true.
      end if
      if (allocated(json_data%dipole)) result_out%dipole = json_data%dipole
   end subroutine result_from_json

   subroutine run_unfragmented_calculation(world_comm, sys_geom, config, result_out, json_data)
      !! Handle unfragmented calculation (nlevel=0)
      !!
      !! For single-molecule mode: Only rank 0 runs (validates single rank)
      !! For multi-molecule mode: ALL ranks can run (each with their own molecule)
      !! For Hessian calculations with multiple ranks: Uses distributed parallelization
      !! If result_out is present, returns result instead of writing JSON
      type(comm_t), intent(in) :: world_comm  !! Global MPI communicator
      type(system_geometry_t), intent(in) :: sys_geom  !! Complete system geometry
      type(driver_config_t), intent(in) :: config  !! Driver configuration (includes method_config, calc_type, etc.)
      type(calculation_result_t), intent(out), optional :: result_out  !! Optional result output
      type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data

      ! For Hessian calculations with multiple ranks, use distributed approach
      if (config%calc_type == CALC_TYPE_HESSIAN .and. world_comm%size() > 1) then
         if (world_comm%rank() == 0) then
            call logger%info(" ")
            call logger%info("Running distributed unfragmented Hessian calculation")
            call logger%info("  MPI ranks: "//to_char(world_comm%size()))
            call logger%info(" ")
         end if
         call distributed_unfragmented_hessian(world_comm, sys_geom, config, json_data)
         return
      end if

      ! Check if this is multi-molecule mode or single-molecule mode
      ! In multi-molecule mode, each rank processes its own molecule
      ! In single-molecule mode, only rank 0 should work
      if (world_comm%size() == 1 .or. world_comm%rank() == 0) then
         ! Either single-rank calculation, or rank 0 in multi-rank setup
         call logger%info(" ")
         call logger%info("Running unfragmented calculation")
         call logger%info("  Calculation type: "//calc_type_to_string(config%calc_type))
         call logger%info(" ")
         call unfragmented_calculation(sys_geom, config, result_out, json_data)
      else if (sys_geom%total_atoms > 0) then
         ! Multi-molecule mode: non-zero rank with a molecule
         call logger%verbose("Rank "//to_char(world_comm%rank())//": Running unfragmented calculation")
         call unfragmented_calculation(sys_geom, config, result_out, json_data)
      end if

   end subroutine run_unfragmented_calculation

   subroutine run_fragmented_calculation(resources, config, sys_geom, bonds, json_data, &
                                         supplied_terms, n_supplied_terms)
      !! Handle fragmented calculation (nlevel > 0)
      !!
      !! Generates fragments, distributes work across MPI processes organized in nodes,
      !! and coordinates many-body expansion calculation using hierarchical parallelism.
      !! If allow_overlapping_fragments=true, uses GMBE with intersection correction.
      type(resources_t), intent(in), target :: resources  !! Resources container (MPI comms, etc.)
      type(driver_config_t), intent(in) :: config  !! Driver configuration (includes method_config, calc_type, etc.)
      type(system_geometry_t), intent(in) :: sys_geom  !! System geometry and fragment info
      type(bond_t), intent(in), optional :: bonds(:)  !! Bond connectivity information
      type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data
      integer, intent(in), optional :: supplied_terms(:, :)
         !! (n_terms, max_level), 1-based monomer indices, zero-padded. Rank 0 only.
      integer(int64), intent(in), optional :: n_supplied_terms

      ! Local variables extracted from config for readability
      integer :: max_level    !! Maximum fragment level for MBE
      logical :: allow_overlapping_fragments  !! Use GMBE for overlapping fragments
      integer :: max_intersection_level  !! Maximum k-way intersection depth for GMBE

      integer(int64) :: total_fragments  !! Total number of fragments generated (int64 to handle large systems)
      integer(default_int) :: supplied_width  !! Columns the caller actually provided
      type(error_t) :: checkpoint_error
      integer, allocatable :: polymers(:, :)  !! Fragment composition array (fragment, monomer_indices)
      integer :: num_nodes   !! Number of compute nodes
      integer :: i, j        !! Loop counters
      integer, allocatable :: node_leader_ranks(:)  !! Ranks of processes that lead each node
      integer :: global_groups  !! Number of global groups for multi-global coordinator
      integer :: nodes_per_group  !! Nodes per group (optional override)
      integer :: group_size  !! Nodes per group after resolving config
      integer :: group_id, group_start
      integer, allocatable :: group_leader_ranks(:)  !! Group leader rank for each node leader
      integer, allocatable :: group_ids(:)  !! Group id for each node leader
      integer, allocatable :: monomers(:)     !! Temporary monomer list for fragment generation
      integer(int64) :: n_expected_frags  !! Expected number of fragments based on combinatorics (int64 to handle large systems)
      integer(int64) :: n_rows      !! Number of rows needed for polymers array (int64 to handle large systems)
      integer :: global_node_rank  !! Global rank if this process leads a node, -1 otherwise
      integer, allocatable :: all_node_leader_ranks(:)  !! Node leader status for all ranks

      ! Polymorphic expansion context for unified MBE/GMBE handling
      class(many_body_expansion_t), allocatable :: expansion

      ! GMBE PIE-based variables
      integer :: n_primaries  !! Number of primary polymers
      integer(int64) :: n_primaries_i64  !! For binomial calculation
      integer, allocatable :: pie_atom_sets(:, :)  !! Unique atom sets (max_atoms, n_pie_terms)
      integer, allocatable :: pie_coefficients(:)  !! PIE coefficient for each term
      integer(int64) :: n_pie_terms  !! Number of unique PIE terms
      type(error_t) :: pie_error  !! Error from PIE enumeration
      type(error_t) :: validation_error  !! Error from semantic validation
      type(fraglist_t) :: supplied_check  !! Supplied terms, for validation only

      ! Extract values from config for readability
      max_level = config%nlevel
      allow_overlapping_fragments = config%allow_overlapping_fragments
      max_intersection_level = config%max_intersection_level

      ! Every input path arrives here, so this is where the system is checked:
      ! the JSON reader, the C interface and a supplied term list alike.
      if (resources%mpi_comms%world_comm%rank() == 0) then
         call validate_system(sys_geom,.not. config%unchecked_input, validation_error, &
                              check_bonds=allocated(sys_geom%bonds))
         if (validation_error%has_error()) then
            call logger%error("invalid system: "//validation_error%get_message())
            call abort_comm(resources%mpi_comms%world_comm, 1)
         end if

         ! Ahead of the branch, because each of the three expansions below
         ! ignores counterpoise in its own way and none of them says so.
         call validation_error%clear()
         call check_counterpoise_support(config, validation_error)
         if (validation_error%has_error()) then
            call logger%error(validation_error%get_message())
            call abort_comm(resources%mpi_comms%world_comm, 1)
         end if
      end if

      ! Generate fragments -- unless the caller brought their own
      if (present(supplied_terms) .and. present(n_supplied_terms)) then
         if (resources%mpi_comms%world_comm%rank() == 0) then
            total_fragments = n_supplied_terms
            allocate (polymers(max(total_fragments, 1_int64), max_level))
            polymers = 0
            ! The caller's array is as wide as their highest term, which need
            ! not be `max_level`: a screen that keeps no trimers hands over a
            ! two-column list for a level-3 expansion, and copying `max_level`
            ! columns from it reads off the end. That produced a monomer index
            ! of 0 the once, which the validator caught; it is undefined
            ! memory and could as easily have been a plausible index.
            supplied_width = int(size(supplied_terms, 2), default_int)
            if (supplied_width > max_level) then
               ! Wider than the expansion allows. Truncating would silently
               ! turn an n-mer into a smaller one, so refuse instead -- but
               ! only if the extra columns actually hold anything.
               if (any(supplied_terms(1:total_fragments, max_level + 1:supplied_width) /= 0)) then
                  call logger%error("invalid fragment list: a term names more than "// &
                                    to_char(max_level)//" monomers, which is the level "// &
                                    "this expansion was configured for")
                  call abort_comm(resources%mpi_comms%world_comm, 1)
               end if
               supplied_width = max_level
            end if
            if (total_fragments > 0) then
               polymers(1:total_fragments, 1:supplied_width) = &
                  supplied_terms(1:total_fragments, 1:supplied_width)
            end if
            ! A supplied list has had no screening or sorting applied to it and
            ! comes from outside, so it is checked before anything is spent on
            ! it -- above all for subset closure, which a reasonable-looking
            ! screen breaks silently.
            call supplied_check%replace(polymers, total_fragments, max_level, validation_error)
            if (.not. validation_error%has_error()) then
               call validate_terms(supplied_check, sys_geom,.not. config%unchecked_input, &
                                   validation_error)
            end if
            call supplied_check%destroy()
            if (validation_error%has_error()) then
               call logger%error("invalid fragment list: "//validation_error%get_message())
               call abort_comm(resources%mpi_comms%world_comm, 1)
            end if

            call logger%info("Using a supplied fragment list:")
            call logger%info("  Total fragments: "//to_char(total_fragments))
            call logger%info("  Max level: "//to_char(max_level))
         end if
      else if (resources%mpi_comms%world_comm%rank() == 0) then
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
                  call apply_distance_screening(polymers, total_fragments, sys_geom, config, max_level)
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
                                          pie_atom_sets, pie_coefficients, n_pie_terms, pie_error)
            if (pie_error%has_error()) then
               call logger%error("GMBE PIE enumeration failed: "//pie_error%get_message())
               call abort_comm(resources%mpi_comms%world_comm, 1)
            end if

            call logger%info("GMBE PIE enumeration complete: "//to_char(n_pie_terms)//" unique subsystems to evaluate")

            ! For now: total_fragments = n_pie_terms (each PIE term is a subsystem to evaluate)
            total_fragments = n_pie_terms
         else
            ! Standard MBE mode. Monomers, then n-mers, then screening, then
            ! the size sort -- all of it in `generate_mbe_term_list`, which a
            ! geometry optimization also calls so that the list it freezes is
            ! the same one this would have built.
            call generate_mbe_term_list(sys_geom, config, max_level, polymers, total_fragments)

            call logger%info("Generated fragments:")
            call logger%info("  Total fragments: "//to_char(total_fragments))
            call logger%info("  Max level: "//to_char(max_level))
         end if
      end if

      ! Broadcast total_fragments to all ranks
      call bcast(resources%mpi_comms%world_comm, total_fragments, 1, 0)

      ! Determine node leaders
      global_node_rank = -1
      if (resources%mpi_comms%node_comm%rank() == 0) global_node_rank = resources%mpi_comms%world_comm%rank()

      allocate (all_node_leader_ranks(resources%mpi_comms%world_comm%size()))
      call allgather(resources%mpi_comms%world_comm, global_node_rank, all_node_leader_ranks)

      num_nodes = count(all_node_leader_ranks /= -1)

      if (resources%mpi_comms%world_comm%rank() == 0) then
         call logger%info("Running with "//to_char(num_nodes)//" node(s)")
      end if

      allocate (node_leader_ranks(num_nodes))
      i = 0
      do j = 1, resources%mpi_comms%world_comm%size()
         if (all_node_leader_ranks(j) /= -1) then
            i = i + 1
            node_leader_ranks(i) = all_node_leader_ranks(j)
         end if
      end do
      deallocate (all_node_leader_ranks)

      global_groups = config%global_groups
      nodes_per_group = config%nodes_per_group

      if (global_groups > 0 .and. nodes_per_group > 0) then
         if (resources%mpi_comms%world_comm%rank() == 0) then
            call logger%error("Both global_groups and nodes_per_group are set; only one may be specified")
         end if
         call abort_comm(resources%mpi_comms%world_comm, 1)
      end if

      if (nodes_per_group > 0) then
         group_size = nodes_per_group
         global_groups = (num_nodes + group_size - 1)/group_size
      else if (global_groups > 0) then
         if (global_groups > num_nodes) global_groups = num_nodes
         group_size = (num_nodes + global_groups - 1)/global_groups
      else
         global_groups = 1
         group_size = num_nodes
      end if

      allocate (group_leader_ranks(num_nodes))
      allocate (group_ids(num_nodes))
      do i = 1, num_nodes
         group_id = min((i - 1)/group_size + 1, global_groups)
         group_ids(i) = group_id
         group_start = (group_id - 1)*group_size + 1
         if (group_start > num_nodes) group_start = num_nodes
         group_leader_ranks(i) = node_leader_ranks(group_start)
      end do

      if (resources%mpi_comms%world_comm%rank() == 0 .and. num_nodes > 1) then
         call logger%info("Multi-global groups: "//to_char(global_groups)//" (nodes_per_group="// &
                          to_char(nodes_per_group)//")")
      end if

      ! Build polymorphic expansion context
      if (config%expansion_kind == "fmo" .or. config%expansion_kind == "ee-mbe") then
         ! FMO or electrostatically embedded MBE. Both are the same machinery
         ! and differ only in what a fragment sees of its neighbours and how
         ! the pieces are added up, so one context serves both and the deck
         ! chooses by name rather than by setting two switches consistently.
         allocate (fmo_context_t :: expansion)
         select type (expansion)
         type is (fmo_context_t)
            call expansion%init(config%method_config, config%calc_type)
            allocate (expansion%sys_geom, source=sys_geom)
            if (present(bonds)) then
               if (allocated(expansion%sys_geom%bonds)) deallocate (expansion%sys_geom%bonds)
               allocate (expansion%sys_geom%bonds, source=bonds)
            end if
            call fragment_owner_map(sys_geom, expansion%owner, expansion%n_fragments)
            expansion%basis = config%method_config%basis_set
            if (config%expansion_kind == "fmo") then
               expansion%esp = "exact"
               expansion%expansion = "fmo"
            else
               expansion%esp = "ptc"
               expansion%expansion = "mbe"
            end if
            expansion%far_field = config%fmo_far_field
            ! The deck's fragmentation level means the same thing here as it
            ! does for MBE: how many fragments at a time.
            expansion%level = max_level
            expansion%resppc = config%fmo_resppc
            expansion%max_outer = config%fmo_max_outer
            expansion%outer_tol = config%fmo_tolerance
            expansion%scf_max_iter = config%fmo_scf_max_iter
            expansion%scf_energy_tol = config%fmo_scf_energy_tol
            expansion%scf_density_tol = config%fmo_scf_density_tol
            expansion%resources => resources
            expansion%node_leader_ranks = node_leader_ranks
            expansion%num_nodes = num_nodes
            expansion%global_groups = global_groups
            expansion%nodes_per_group = nodes_per_group
            expansion%group_leader_ranks = group_leader_ranks
            expansion%group_ids = group_ids
         end select
      else if (allow_overlapping_fragments) then
         ! GMBE: allocate gmbe_context_t
         allocate (gmbe_context_t :: expansion)
         select type (expansion)
         type is (gmbe_context_t)
            call expansion%init(config%method_config, config%calc_type)
            allocate (expansion%sys_geom, source=sys_geom)
            if (present(bonds)) then
               ! source=sys_geom above already copies its bonds when the caller
               ! embeds them there (the C API does), so replace rather than
               ! allocate -- allocating an allocated component aborts.
               if (allocated(expansion%sys_geom%bonds)) deallocate (expansion%sys_geom%bonds)
               allocate (expansion%sys_geom%bonds, source=bonds)
            end if
            expansion%n_pie_terms = n_pie_terms
            if (resources%mpi_comms%world_comm%rank() == 0) then
               allocate (expansion%pie_atom_sets, source=pie_atom_sets)
               allocate (expansion%pie_coefficients, source=pie_coefficients)
            end if
            expansion%resources => resources
            expansion%node_leader_ranks = node_leader_ranks
            expansion%num_nodes = num_nodes
            expansion%global_groups = global_groups
            expansion%nodes_per_group = nodes_per_group
            expansion%group_leader_ranks = group_leader_ranks
            expansion%group_ids = group_ids
         end select
      else
         ! Standard MBE: allocate mbe_context_t
         allocate (mbe_context_t :: expansion)
         select type (expansion)
         type is (mbe_context_t)
            call expansion%init(config%method_config, config%calc_type)
            allocate (expansion%sys_geom, source=sys_geom)
            if (present(bonds)) then
               ! source=sys_geom above already copies its bonds when the caller
               ! embeds them there (the C API does), so replace rather than
               ! allocate -- allocating an allocated component aborts.
               if (allocated(expansion%sys_geom%bonds)) deallocate (expansion%sys_geom%bonds)
               allocate (expansion%sys_geom%bonds, source=bonds)
            end if
            expansion%total_fragments = total_fragments
            if (resources%mpi_comms%world_comm%rank() == 0) then
               allocate (expansion%polymers, source=polymers)
            end if
            expansion%max_level = max_level
            expansion%resources => resources
            expansion%node_leader_ranks = node_leader_ranks
            expansion%num_nodes = num_nodes
            expansion%global_groups = global_groups
            expansion%nodes_per_group = nodes_per_group
            expansion%group_leader_ranks = group_leader_ranks
            expansion%group_ids = group_ids
         end select
      end if

      ! Opened on rank 0 only: it is the rank that collects every result, so
      ! it is the only one with anything to write, and N ranks appending to one
      ! file would interleave.
      ! GMBE keys its terms by atom set, not by monomer tuple, so the file
      ! this writes would not describe them. Announced rather than opened and
      ! left unused -- a checkpoint that silently records nothing is worse
      ! than no checkpoint, because it is believed.
      if (len_trim(config%checkpoint_file) > 0 .and. allow_overlapping_fragments .and. &
          resources%mpi_comms%world_comm%rank() == 0) then
         call logger%warning("Checkpointing does not support GMBE; this run will write "// &
                             "nothing and resume nothing. Its PIE terms are atom sets "// &
                             "rather than monomer tuples, which the file cannot express.")
      end if

      if (resources%mpi_comms%world_comm%rank() == 0 .and. &
          .not. allow_overlapping_fragments .and. &
          len_trim(config%checkpoint_file) > 0) then
         call expansion%checkpoint%open(trim(config%checkpoint_file), &
                                        calculation_fingerprint(sys_geom, config%method_config, &
                                                                config%calc_type), &
                                        max_level + 1, &
                                        config%calc_type == CALC_TYPE_ENERGY, &
                                        checkpoint_error)
         if (checkpoint_error%has_error()) then
            ! A checkpoint from another calculation is not a warning. Its
            ! energies would be spliced into this one and the total would come
            ! out converged and meaningless.
            call logger%error(checkpoint_error%get_message())
            call abort_comm(resources%mpi_comms%world_comm, 1)
         end if
      end if

      ! Execute calculation using polymorphic dispatch
      if (resources%mpi_comms%world_comm%size() == 1) then
         call logger%info("Running in serial mode (single MPI rank)")
         call expansion%run_serial(json_data)
      else
         call expansion%run_distributed(json_data)
      end if
      call expansion%checkpoint%close()

      ! Clean up expansion context
      select type (expansion)
      type is (mbe_context_t)
         call expansion%destroy()
      type is (gmbe_context_t)
         call expansion%destroy()
      type is (fmo_context_t)
         call expansion%destroy()
      end select
      deallocate (expansion)

      ! Cleanup
      if (resources%mpi_comms%world_comm%rank() == 0) then
         if (allocated(polymers)) deallocate (polymers)
         if (allocated(node_leader_ranks)) deallocate (node_leader_ranks)
         if (allocated(group_leader_ranks)) deallocate (group_leader_ranks)
         if (allocated(group_ids)) deallocate (group_ids)
         if (allocated(pie_atom_sets)) deallocate (pie_atom_sets)
         if (allocated(pie_coefficients)) deallocate (pie_coefficients)
      end if

   end subroutine run_fragmented_calculation

   subroutine fragment_owner_map(sys_geom, owner, n_fragments)
      !! Which fragment each atom belongs to, from the declared monomers
      !!
      !! `fragment_atoms` is 0-indexed and padded to the widest fragment, so
      !! only the first `fragment_sizes(f)` entries of a column mean anything.
      type(system_geometry_t), intent(in) :: sys_geom
      integer, allocatable, intent(out) :: owner(:)
      integer, intent(out) :: n_fragments

      integer :: f, k

      n_fragments = sys_geom%n_monomers
      allocate (owner(sys_geom%total_atoms), source=0)
      do f = 1, n_fragments
         do k = 1, sys_geom%fragment_sizes(f)
            owner(sys_geom%fragment_atoms(k, f) + 1) = f
         end do
      end do
   end subroutine fragment_owner_map

   subroutine run_multi_molecule_calculations(resources, mqc_config)
      !! Run calculations for multiple molecules with MPI parallelization
      !! Each molecule is independent, so assign one molecule per rank
      use mqc_config_types, only: mqc_config_t
      use mqc_config_adapter, only: config_to_system_geometry
      use mqc_error, only: error_t
      use mqc_io_helpers, only: set_molecule_suffix, get_output_json_filename
      use mqc_json, only: merge_multi_molecule_json

      type(resources_t), intent(in) :: resources
      type(mqc_config_t), intent(in) :: mqc_config

      type(driver_config_t) :: config
      type(system_geometry_t) :: sys_geom
      type(resources_t) :: mol_resources
      type(error_t) :: error
      integer :: imol, my_rank, num_ranks, color
      integer :: molecules_processed
      character(len=:), allocatable :: mol_name
      logical :: has_fragmented_molecules
      character(len=256), allocatable :: individual_json_files(:)

      my_rank = resources%mpi_comms%world_comm%rank()
      num_ranks = resources%mpi_comms%world_comm%size()

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
            call config_to_driver(mqc_config, config, molecule_index=imol, &
                                  node_rank=resources%mpi_comms%node_comm%rank())

            ! Convert geometry for this molecule
            call config_to_system_geometry(mqc_config, sys_geom, error, molecule_index=imol)
            if (error%has_error()) then
               call error%add_context("mqc_driver:run_multi_molecule_calculation")
               if (my_rank == 0) then
                  call logger%error("Error converting geometry for "//mol_name//": "//error%get_full_trace())
               end if
               call abort_comm(resources%mpi_comms%world_comm, 1)
            end if

            ! Set output filename suffix for this molecule
            call set_molecule_suffix("_"//trim(mol_name))

            ! Run calculation for this molecule
            call run_calculation(resources, config, sys_geom, mqc_config%molecules(imol)%bonds)

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
               call config_to_driver(mqc_config, config, molecule_index=imol, &
                                     node_rank=resources%mpi_comms%node_comm%rank())

               ! Convert geometry for this molecule
               call config_to_system_geometry(mqc_config, sys_geom, error, molecule_index=imol)
               if (error%has_error()) then
                  call error%add_context("mqc_driver:run_multi_molecule_calculation")
  call logger%error("Rank "//to_char(my_rank)//": Error converting geometry for "//mol_name//": "//error%get_full_trace())
                  call abort_comm(resources%mpi_comms%world_comm, 1)
               end if

               ! Set output filename suffix for this molecule
               call set_molecule_suffix("_"//trim(mol_name))

               ! Run calculation for this molecule (all ranks write JSON in parallel mode)
               call run_calculation(resources, config, sys_geom, mqc_config%molecules(imol)%bonds, &
                                    all_ranks_write_json=.true.)

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
      call resources%mpi_comms%world_comm%barrier()

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

   subroutine run_sapt(config, sys_geom, rank, write_output, result_out)
      !! SAPT0 between the deck's two fragments
      !!
      !! The monomers are the deck's own `fragments`, so no new keyword is needed
      !! -- a SAPT deck is a geometry, a basis, and the partition that says which
      !! atoms are which monomer:
      !!
      !!     "model":     {"method": "sapt0", "basis": "6-31g"},
      !!     "molecules": [{"xyz": "dimer.xyz", "fragments": [[0,1,2],[3,4,5]], ...}]
      !!
      !! **Exactly two fragments.** SAPT partitions the Hamiltonian as
      !! `H_A + H_B + V`; there is no slot for a third monomer, and a cluster is
      !! a separate SAPT calculation per pair rather than one calculation. That
      !! is the theory rather than a limit of this code -- see
      !! `validation/check_sapt.f90`, which walks the pairs of a six-water prism.
      !!
      !! Rank zero only. The pairs of a cluster are the obvious thing to
      !! distribute, and a single pair is not.
      use mqc_libcint_bridge, only: run_libcint_sapt0
      use mqc_program_limits, only: N_SAPT_TERMS
      use mqc_elements, only: element_number_to_symbol
      use mqc_json_output_types, only: OUTPUT_MODE_UNFRAGMENTED
      type(driver_config_t), intent(in) :: config
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: rank
      logical, intent(in) :: write_output
      type(calculation_result_t), intent(out), optional :: result_out

      type(error_t) :: err
      type(json_output_data_t) :: json_data
      real(dp) :: terms(N_SAPT_TERMS)
      integer, allocatable :: z_a(:), z_b(:)
      real(dp), allocatable :: xyz_a(:, :), xyz_b(:, :)
      character(len=8), allocatable :: sym_a(:), sym_b(:)
      integer :: na, nb, i, a

      if (rank /= 0) return

      if (sys_geom%n_monomers /= 2) then
         call refuse(result_out, "SAPT: the system must have exactly two "// &
                     "fragments; this one has "//to_char(sys_geom%n_monomers)// &
                     ". SAPT is a two-body theory, so a cluster is one "// &
                     "calculation per pair.")
         return
      end if

      na = sys_geom%fragment_sizes(1)
      nb = sys_geom%fragment_sizes(2)
      allocate (z_a(na), xyz_a(3, na), sym_a(na))
      allocate (z_b(nb), xyz_b(3, nb), sym_b(nb))
      do i = 1, na
         a = sys_geom%fragment_atoms(i, 1) + 1        ! stored 0-based
         z_a(i) = sys_geom%element_numbers(a)
         xyz_a(:, i) = sys_geom%coordinates(:, a)
         sym_a(i) = element_number_to_symbol(z_a(i))
      end do
      do i = 1, nb
         a = sys_geom%fragment_atoms(i, 2) + 1
         z_b(i) = sys_geom%element_numbers(a)
         xyz_b(:, i) = sys_geom%coordinates(:, a)
         sym_b(i) = element_number_to_symbol(z_b(i))
      end do

      call run_libcint_sapt0(z_a, sym_a, xyz_a, z_b, sym_b, xyz_b, &
                             config%method_config%basis_set, terms, err)
      if (err%has_error()) then
         call refuse(result_out, "SAPT: "//err%get_message())
         return
      end if

      ! The response terms are what enters the total; the uncoupled ones and
      ! the S^2 exchange are printed beside them because they are what another
      ! code's output is quoted in, and a term can only be checked against the
      ! same term.
      call logger%info("============================================")
      call logger%info("  SAPT0 interaction energy, Hartree")
      call logger%info("    electrostatics        "//to_char(terms(1)))
      call logger%info("    exchange              "//to_char(terms(3)))
      call logger%info("      (S^2 approximation) "//to_char(terms(2)))
      call logger%info("    induction             "//to_char(terms(5)))
      call logger%info("      (uncoupled)         "//to_char(terms(4)))
      call logger%info("    exchange-induction    "//to_char(terms(7)))
      call logger%info("      (uncoupled)         "//to_char(terms(6)))
      call logger%info("    dispersion            "//to_char(terms(8)))
      call logger%info("    exchange-dispersion   "//to_char(terms(9)))
      call logger%info("    delta HF              "//to_char(terms(10)))
      call logger%info("    --")
      call logger%info("    supermolecular HF     "//to_char(terms(11)))
      call logger%info("    total                 "//to_char(terms(12)))
      call logger%info("============================================")

      if (present(result_out)) then
         ! The `scf` slot, being the reference `energy_t%total()` sums. This is an
         ! interaction energy and not an SCF energy; nothing here pretends otherwise.
         result_out%energy%scf = terms(N_SAPT_TERMS)
         result_out%has_energy = .true.
      end if

      if (write_output .and. .not. config%skip_json_output) then
         json_data%output_mode = OUTPUT_MODE_UNFRAGMENTED
         json_data%total_energy = terms(N_SAPT_TERMS)
         json_data%has_energy = .true.
         json_data%sapt_terms = terms
         json_data%has_sapt = .true.
         json_data%fragment_breakdown = config%fragment_breakdown
         call write_json_output(json_data)
         call json_data%destroy()
      end if

      deallocate (z_a, xyz_a, sym_a, z_b, xyz_b, sym_b)
   end subroutine run_sapt

   subroutine run_efp(config, sys_geom, rank, write_output, result_out)
      !! The interaction energy of a set of effective fragments
      !!
      !! Each fragment of the deck names a potential; the backend loads them, places
      !! each on the atoms the deck gave it, turns it into that orientation and
      !! evaluates the five EFP2 terms. There is no SCF here: the wavefunctions were
      !! solved when the potentials were made.
      !!
      !! The work itself is in `run_libcint_efp` rather than here, because all of it
      !! lives behind `MQC_ENABLE_LIBCINT` and this file compiles either way. The stub
      !! declines with the same signature and names the build option.
      !!
      !! Rank zero only, for now. The pair loop is the obvious thing to distribute and
      !! is already flat, but the terms are milliseconds on a dimer and shaping that
      !! before there is a cluster to shape it around would be guessing.
      use mqc_libcint_bridge, only: run_libcint_efp
      use mqc_program_limits, only: N_EFP_TERMS
      use mqc_json_output_types, only: OUTPUT_MODE_UNFRAGMENTED
      type(driver_config_t), intent(in) :: config
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: rank
      logical, intent(in) :: write_output
      type(calculation_result_t), intent(out), optional :: result_out

      type(error_t) :: err
      type(json_output_data_t) :: json_data
      real(dp) :: terms(N_EFP_TERMS)
      integer :: n

      if (rank /= 0) return

      n = sys_geom%n_monomers
      if (.not. allocated(config%fragment_potentials)) then
         call refuse(result_out, "EFP: no fragment carries a potential -- name "// &
                     "one per fragment in 'fragment_potentials' in the deck, or "// &
                     "through mqc_system_set_fragment_potentials from the C and "// &
                     "Python interfaces")
         return
      end if
      if (size(config%fragment_potentials) /= n) then
         call refuse(result_out, "EFP: the system has "//to_char(n)//" fragments "// &
                     "but "//to_char(size(config%fragment_potentials))//" potentials")
         return
      end if

      call run_libcint_efp(config%fragment_potentials, sys_geom%fragment_sizes, &
                           sys_geom%fragment_atoms, sys_geom%coordinates, terms, err)
      if (err%has_error()) then
         call refuse(result_out, "EFP: "//err%get_message())
         return
      end if

      call logger%info("============================================")
      call logger%info("  EFP2 interaction energy, Hartree")
      call logger%info("    electrostatics      "//to_char(terms(1)))
      call logger%info("    polarization        "//to_char(terms(2)))
      call logger%info("    exchange repulsion  "//to_char(terms(3)))
      call logger%info("    dispersion          "//to_char(terms(4)))
      call logger%info("    charge transfer     "//to_char(terms(5)))
      call logger%info("    total               "//to_char(terms(6)))
      call logger%info("============================================")

      if (present(result_out)) then
         ! In the `scf` slot because that is the one `energy_t%total()` sums as the
         ! reference, and an EFP interaction energy has no correlation correction
         ! sitting on top of it. It is not an SCF energy and nothing here pretends
         ! otherwise -- there is no wavefunction solved in this routine at all.
         result_out%energy%scf = terms(6)
         result_out%has_energy = .true.
      end if

      ! The same summary every other run leaves behind, so a driven run or the
      ! validation harness can read the number back rather than scrape the log.
      ! `UNFRAGMENTED` because that is the shape of what is written -- one energy for
      ! one system -- not a claim that the system has no fragments.
      if (write_output .and. .not. config%skip_json_output) then
         json_data%output_mode = OUTPUT_MODE_UNFRAGMENTED
         json_data%total_energy = terms(6)
         json_data%has_energy = .true.
         json_data%fragment_breakdown = config%fragment_breakdown
         call write_json_output(json_data)
         call json_data%destroy()
      end if
   end subroutine run_efp

   subroutine refuse(result_out, message)
      !! Log a refusal and, when there is a caller, put it on the result
      !!
      !! The two interaction-energy paths used to log and return. From `main`
      !! that is a red line on the terminal and an exit; from a session it was
      !! a clean zero handed back as a success, because `mqc_run` reads the
      !! energy off a result nothing had marked as failed. A zero interaction
      !! energy is a physically plausible number, which is what makes it the
      !! worst shape a failure can take here.
      type(calculation_result_t), intent(inout), optional :: result_out
      character(len=*), intent(in) :: message

      call logger%error(message)
      if (present(result_out)) then
         call result_out%error%set(ERROR_VALIDATION, message)
         result_out%has_error = .true.
      end if
   end subroutine refuse

   subroutine run_makefp(config, sys_geom, rank, result_out)
      !! Build an effective fragment potential for the whole system and write it
      !!
      !! Only rank zero does anything: the work is one SCF and a set of response
      !! solves, all of them threaded inside the integral backend, and the write is
      !! a single file. Distributing it would mean every rank recomputing the same
      !! potential to overwrite the same path.
      use mqc_calc_types, only: CALC_TYPE_MAKEFP
      use mqc_elements, only: element_number_to_symbol
      use mqc_libcint_bridge, only: run_libcint_makefp
      use mqc_io_helpers, only: get_basename
      type(driver_config_t), intent(in) :: config
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: rank
      type(calculation_result_t), intent(out), optional :: result_out
         !! Carries the failure, if there is one. There is no energy to put on
         !! it -- MAKEFP produces a file -- but a caller that gets neither an
         !! energy nor an error has no way to tell a written potential from a
         !! backend that declined.

      type(error_t) :: err
      character(len=8), allocatable :: symbols(:)
      character(len=:), allocatable :: path, name
      integer :: i
      real(dp), allocatable :: e_tol, d_tol
         !! Left unallocated when the deck did not name a tolerance, which makes
         !! the corresponding argument absent at the call below and leaves MAKEFP
         !! on its own tighter default. Allocating on assignment is the only way
         !! to pass "nothing" without branching on each flag separately.

      if (rank /= 0) return

      allocate (symbols(sys_geom%total_atoms))
      do i = 1, sys_geom%total_atoms
         symbols(i) = element_number_to_symbol(sys_geom%element_numbers(i))
      end do

      ! Named and placed off the deck, the way the JSON output is: a run on
      ! water.json leaves water.efp beside it.
      name = trim(get_basename())
      path = trim(name)//".efp"

      call logger%info("Building an effective fragment potential")
      ! `keywords.scf.density_fitting` and `keywords.scf.aux_basis_set` already
      ! existed for the SCF, and mean here what they mean there: fit the two-electron
      ! integrals against that auxiliary basis. What they reach in a MAKEFP run is the
      ! response Hessian, which is where the time goes.
      !
      ! `keywords.scf.tolerance` and `keywords.scf.density_tolerance` reach the SCF
      ! here only if the deck named them. MAKEFP's own 1e-10/1e-8 is deliberate --
      ! the multipoles and the response come off that density -- so it is not the
      ! shared 1e-6 default's to loosen. A deck that asks outright still wins: a
      ! user who set 1e-6 and watched this iterate past 37 steps was being ignored.
      !
      ! `keywords.efp` goes down by value instead, with no flags. Those four have a
      ! single default each -- `mqc_calculation_defaults` holds it and the solver
      ! reads the same constant -- so passing what a silent deck carries is passing
      ! the solver its own number, and there is nothing for a flag to distinguish.
      if (config%method_config%scf%energy_convergence_set) then
         e_tol = config%method_config%scf%energy_convergence
      end if
      if (config%method_config%scf%density_convergence_set) then
         d_tol = config%method_config%scf%density_convergence
      end if
      if (config%method_config%scf%density_fitting) then
         call run_libcint_makefp(sys_geom%element_numbers, symbols, sys_geom%coordinates, &
                                 config%method_config%basis_set, name, path, err, &
                                 charge=sys_geom%charge, verbose=.true., &
                                 aux_basis=trim(config%method_config%scf%aux_basis_set), &
                                 guess=trim(config%method_config%scf%guess), &
                                 energy_tol=e_tol, density_tol=d_tol, &
                                 vdwscl=config%method_config%efp%vdw_scale, &
                                 dynamic_tol=config%method_config%efp%dynamic_tolerance, &
                                 dynamic_maxiter=config%method_config%efp%dynamic_maxiter, &
                                 response=config%method_config%efp%response)
      else
         call run_libcint_makefp(sys_geom%element_numbers, symbols, sys_geom%coordinates, &
                                 config%method_config%basis_set, name, path, err, &
                                 charge=sys_geom%charge, verbose=.true., &
                                 guess=trim(config%method_config%scf%guess), &
                                 energy_tol=e_tol, density_tol=d_tol, &
                                 vdwscl=config%method_config%efp%vdw_scale, &
                                 dynamic_tol=config%method_config%efp%dynamic_tolerance, &
                                 dynamic_maxiter=config%method_config%efp%dynamic_maxiter, &
                                 response=config%method_config%efp%response)
      end if
      if (err%has_error()) then
         call refuse(result_out, "MAKEFP failed: "//err%get_message())
         return
      end if
      call logger%info("Wrote "//trim(path))
   end subroutine run_makefp

end module mqc_driver
