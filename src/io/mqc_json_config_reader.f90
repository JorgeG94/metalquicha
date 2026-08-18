!! Reader for MQC JSON input files
module mqc_json_config_reader
   !! The only input reader. Parses a JSON deck into `mqc_config_t`.
   !!
   !! Until 0.2.0 this program read a section-based `.mqc` format generated from
   !! the user's JSON by a Python helper, so ~2800 lines and a runtime Python
   !! dependency stood between an input file and this type. Both are gone.
   !! json-fortran was already required by the output writer and the basis
   !! reader, so reading JSON here costs no new dependency.
   !!
   !! Three behaviours are not obvious from the schema alone:
   !!
   !!   * **A molecule gives either `xyz` or `symbols` + `geometry`, never
   !!     both.** The `xyz` path is resolved relative to the directory holding
   !!     the JSON file, not the working directory, so an input deck and its
   !!     geometries move together.
   !!   * **`is_broken` is derived, not declared.** A bond is broken when the
   !!     set of fragments containing atom i differs from the set containing
   !!     atom j. Nothing in the JSON says so; it falls out of the fragment
   !!     definitions. Comparing *sets* rather than asking whether both atoms
   !!     share some fragment is what makes it right when fragments overlap.
   !!   * **Absent keys leave defaults alone.** The `optional_*` accessors
   !!     write only when the key is present, so `mqc_config_types` holds the
   !!     defaults and nothing here restates them. Route two keys through one
   !!     local and the first one's value survives into the second -- that bug
   !!     has been made once already.
   !!
   !! Atom indices are 0-based throughout, in the JSON and in the config, and
   !! are stored as read.
   !!
   !! Unknown keys, missing required ones, and the cross-checks between a
   !! molecule's parallel lists are `mqc_json_schema`'s job, run over the
   !! document before any of this. That separation is why the accessors below
   !! can stay as simple as they are: by the time they run, the deck is known
   !! to contain only keys this module knows about.
   use pic_io, only: to_char
   use pic_types, only: dp
   use mqc_program_limits, only: MAX_ELEMENT_SYMBOL_LEN, MAX_MBE_LEVEL
   use mqc_geometry, only: geometry_type
   use mqc_error, only: error_t, ERROR_IO, ERROR_PARSE, ERROR_VALIDATION
   use mqc_calc_types, only: calc_type_from_string
   use mqc_method_types, only: parse_method_string, method_spin_scaling, &
                               method_wants_density_fitting, method_type_to_string
   use mqc_cuest_iface, only: parse_backend_name, method_runs_on_cuest, BACKEND_CUEST
   use mqc_cuest_bridge, only: cuest_backend_available
   use mqc_config_types, only: mqc_config_t, input_fragment_t, bond_t
   use mqc_xyz_reader, only: read_xyz_file
   use mqc_json_schema, only: ensure_valid_json
   use json_module, only: json_file
   implicit none
   private

   public :: read_json_config_file  !! Parse a JSON input file into mqc_config_t
   public :: read_json_config_text  !! Parse settings from a JSON string, no molecules

   !> Element symbol width comes from `mqc_program_limits`, the same source
   !  `mqc_xyz_reader` uses, so a geometry given inline and one read from an
   !  .xyz file produce identical `elements` arrays rather than merely
   !  equivalent ones.

contains

   subroutine read_json_config_file(filename, config, error)
      !! Read and parse a JSON format input file
      character(len=*), intent(in) :: filename
      type(mqc_config_t), intent(out) :: config
      type(error_t), intent(out) :: error

      type(json_file) :: json
      logical :: file_exists

      inquire (file=filename, exist=file_exists)
      if (.not. file_exists) then
         call error%set(ERROR_IO, "Input file not found: "//trim(filename))
         return
      end if

      call json%initialize()
      call json%load(filename=filename)
      if (json%failed()) then
         call error%set(ERROR_PARSE, "Could not parse JSON input file: "//trim(filename))
         call json%destroy()
         return
      end if

      ! Validate before reading, so a misspelled key is reported as itself
      ! rather than silently taking the default of the key it was meant to be.
      ! One place to load and one to destroy; everything that can fail lives
      ! in the two routines below and simply returns.
      call ensure_valid_json(json, error)
      if (.not. error%has_error()) then
         call populate_config(json, directory_of(filename), config, error)
      end if
      call json%destroy()

      if (error%has_error()) call error%add_context("mqc_json_config_reader:"//trim(filename))
   end subroutine read_json_config_file

   subroutine read_json_config_text(text, config, error)
      !! Parse settings -- method, driver, keywords -- from a JSON string
      !!
      !! The same document a deck uses, less the molecules: this is how a
      !! caller that already holds a system says what to *do* with it. The
      !! string form matters because it is what crosses to the workers. A
      !! `method_config_t` is sixty-odd fields across seven nested sub-configs,
      !! and transcribing those into a broadcast means a field missed there is
      !! a worker quietly running a different method than the one asked for.
      !! Sending the document and parsing it on arrival has no such field list;
      !! it also puts the workers through the validator rather than around it.
      !!
      !! A molecules block is refused rather than ignored. It would otherwise
      !! be the geometry a caller believed they were running, silently losing
      !! to the one on the handle.
      character(len=*), intent(in) :: text
      type(mqc_config_t), intent(out) :: config
      type(error_t), intent(out) :: error

      type(json_file) :: json

      call json%initialize()
      call json%deserialize(text)
      if (json%failed()) then
         call error%set(ERROR_PARSE, "Could not parse JSON settings")
         call json%destroy()
         return
      end if

      call ensure_valid_json(json, error, settings_only=.true.)
      if (.not. error%has_error()) then
         ! No base directory: an xyz path is only ever reached through a
         ! molecules block, which this mode does not have.
         call populate_config(json, "", config, error, settings_only=.true.)
      end if
      call json%destroy()

      if (error%has_error()) call error%add_context("mqc_json_config_reader:settings")
   end subroutine read_json_config_text

   subroutine populate_config(json, base_dir, config, error, settings_only)
      !! Fill every section of the config from an already-loaded document
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: base_dir  !! Directory holding the JSON, for xyz paths
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error
      logical, intent(in), optional :: settings_only  !! Stop before molecules

      character(len=:), allocatable :: text
      integer :: n_mol, imol
      logical :: found, settings
      logical :: backend_named

      ! The two defaults that are not declared on the type itself.
      config%log_level = "info"
      config%fragment_breakdown = "csv"

      ! ---- schema ----------------------------------------------------------
      call require_string(json, "schema.name", config%schema_name, error)
      if (error%has_error()) return
      call require_string(json, "schema.version", config%schema_version, error)
      if (error%has_error()) return
      call optional_int(json, "schema.index_base", config%index_base)
      call optional_string(json, "schema.units", config%units)
      if (.not. allocated(config%units)) config%units = "angstrom"

      ! ---- model -----------------------------------------------------------
      call require_string(json, "model.method", text, error)
      if (error%has_error()) return
      config%method = parse_method_string(text)
      ! "scs-mp2" and "sos-mp2" parse to the same method type as "mp2", so the
      ! spin scaling has to be picked up here, while the spelling still exists.
      ! Left to the type alone they would run plain MP2 and report it as SCS.
      call method_spin_scaling(text, config%corr_scs, config%corr_scs_ss, config%corr_scs_os)
      ! Likewise "ri-mp2" and "df-mp2", which the method type cannot distinguish
      ! from "mp2". A later keyword can still turn it off.
      config%corr_density_fitting = method_wants_density_fitting(text)
      call optional_string(json, "model.basis", config%basis)
      call optional_string(json, "model.aux_basis", config%aux_basis)
      call optional_string(json, "model.functional", config%functional)

      ! ---- driver ----------------------------------------------------------
      call require_string(json, "driver", text, error)
      if (error%has_error()) return
      config%calc_type = calc_type_from_string(text)

      ! ---- system ----------------------------------------------------------
      ! Written straight into the config fields, which already hold the
      ! defaults set above. Routing them through a shared local instead would
      ! carry the previous key's value forward whenever one is absent, since
      ! `optional_string` leaves its target alone rather than clearing it.
      call optional_string(json, "system.logger.level", config%log_level)
      call optional_logical_seen(json, "system.gpu", config%gpu, config%gpu_set)
      call optional_logical(json, "system.skip_json_output", config%skip_json_output)
      call optional_logical(json, "system.unchecked_input", config%unchecked_input)
      call optional_string(json, "system.checkpoint", config%checkpoint_file)
      call optional_string(json, "system.fragment_breakdown", config%fragment_breakdown)

      ! ---- keywords --------------------------------------------------------
      call optional_int(json, "keywords.scf.maxiter", config%scf_maxiter)
      call optional_real(json, "keywords.scf.tolerance", config%scf_tolerance)
      call optional_logical(json, "keywords.scf.unrestricted", config%scf_unrestricted)
      call optional_string(json, "keywords.scf.guess", config%scf_guess)
      call optional_string(json, "keywords.guess.type", config%guess_type)
      call read_guess_steps(json, config, error)
      if (error%has_error()) return
      call optional_logical(json, "keywords.scf.allow_crap_scf", config%allow_crap_scf)
      call optional_logical(json, "keywords.scf.density_fitting", &
                            config%scf_density_fitting)
      call optional_logical(json, "keywords.correlation.freeze_core", &
                            config%corr_freeze_core)
      call optional_int(json, "keywords.correlation.n_frozen_core", &
                        config%corr_n_frozen_core)
      ! After the method name, so an explicit keyword wins over the spelling.
      call optional_logical(json, "keywords.correlation.density_fitting", &
                            config%corr_density_fitting)

      ! A fitted SCF needs a basis to fit against. There is a default, but it is a
      ! JK set that is wrong for some of what asks for fitting, so rather than pick
      ! one silently a fitted run must name its own. `aux_basis` is unallocated
      ! exactly when `model.aux_basis` was absent, which is the case to refuse.
      if (config%scf_density_fitting .and. .not. allocated(config%aux_basis)) then
         call error%set(ERROR_VALIDATION, "keywords.scf.density_fitting is on but no "// &
                        "auxiliary basis was given; set model.aux_basis")
         return
      end if

      call optional_logical(json, "keywords.correlation.scs", config%corr_scs)
      call optional_real(json, "keywords.correlation.scs_ss", config%corr_scs_ss)
      call optional_real(json, "keywords.correlation.scs_os", config%corr_scs_os)
      call optional_int(json, "keywords.cc.maxiter", config%cc_maxiter)
      call optional_real(json, "keywords.cc.tolerance", config%cc_tolerance)
      call optional_logical(json, "keywords.cc.diis", config%cc_diis)
      call optional_int(json, "keywords.cc.diis_size", config%cc_diis_size)
      ! Recorded as "was it named" as well as "what was it", because the method
      ! name is the usual source of this and only an explicit keyword may
      ! override it. `optional_logical` leaves its target alone when the key is
      ! absent, which cannot distinguish absent from false on its own.
      call optional_logical_seen(json, "keywords.cc.triples", config%cc_triples, &
                                 config%cc_triples_set)
      call optional_int(json, "keywords.dft.grid_level", config%dft_grid_level)
      ! The continuum is switched on by the presence of the block rather than by a
      ! flag inside it: a deck that names a dielectric wants solvent, and a
      ! separate "enabled" would let the two disagree.
      block
         ! Read into a deferred-length local, since the config field is fixed
         ! width and `optional_string` takes an allocatable.
         character(len=:), allocatable :: backend_text
         call json%get("backend", backend_text, backend_named)
         if (backend_named .and. allocated(backend_text)) config%backend = backend_text
      end block
      ! Both spellings are now in hand, so the GPU request can be settled and
      ! checked against what this build and this method can actually do.
      call resolve_gpu_request(config, backend_named, error)
      if (error%has_error()) return
      call json%info("keywords.pcm", found=found)
      config%pcm_enabled = found
      call optional_real(json, "keywords.pcm.dielectric", config%pcm_dielectric)
      call optional_int(json, "keywords.pcm.angular_points", config%pcm_angular_points)
      call optional_real(json, "keywords.pcm.radii_scale", config%pcm_radii_scale)
      call optional_real(json, "keywords.pcm.zeta", config%pcm_zeta)
      call optional_real(json, "keywords.pcm.tolerance", config%pcm_tolerance)
      call optional_int(json, "keywords.pcm.max_iter", config%pcm_max_iter)
      call optional_int(json, "keywords.dft.radial_points", config%dft_radial_points)
      call optional_int(json, "keywords.dft.angular_points", config%dft_angular_points)

      ! Both spellings of the displacement key are accepted, as the JSON
      ! generator used to allow.
      call optional_real(json, "keywords.hessian.finite_difference_displacement", &
                         config%hessian_displacement)
      call optional_real(json, "keywords.hessian.displacement", config%hessian_displacement)
      call optional_real(json, "keywords.hessian.temperature", config%hessian_temperature)
      call optional_real(json, "keywords.hessian.pressure", config%hessian_pressure)

      call optional_real(json, "keywords.aimd.dt", config%aimd_dt)
      call optional_real(json, "keywords.aimd.timestep", config%aimd_dt)
      call optional_int(json, "keywords.aimd.nsteps", config%aimd_nsteps)
      call optional_int(json, "keywords.aimd.steps", config%aimd_nsteps)
      call optional_real(json, "keywords.aimd.initial_temperature", config%aimd_initial_temperature)
      call optional_real(json, "keywords.aimd.temperature", config%aimd_initial_temperature)
      call optional_int(json, "keywords.aimd.output_frequency", config%aimd_output_frequency)
      call optional_int(json, "keywords.aimd.output_freq", config%aimd_output_frequency)

      ! Geometry optimization. Both spellings of the two that have an obvious
      ! synonym, the way the hessian and aimd blocks above accept theirs.
      call optional_int(json, "keywords.optimization.max_steps", config%opt_max_steps)
      call optional_int(json, "keywords.optimization.steps", config%opt_max_steps)
      call optional_real(json, "keywords.optimization.gradient_tolerance", &
                         config%opt_gradient_tolerance)
      call optional_real(json, "keywords.optimization.tolerance", config%opt_gradient_tolerance)
      call optional_real(json, "keywords.optimization.energy_tolerance", &
                         config%opt_energy_tolerance)
      call optional_real(json, "keywords.optimization.max_step", config%opt_max_step)
      call optional_string(json, "keywords.optimization.coordinates", config%opt_coordinates)
      call optional_string(json, "keywords.optimization.coordinate_system", config%opt_coordinates)
      call optional_string(json, "keywords.optimization.algorithm", config%opt_algorithm)
      call optional_string(json, "keywords.optimization.optimizer", config%opt_algorithm)
      call optional_int(json, "keywords.optimization.lbfgs_memory", config%opt_lbfgs_memory)
      call optional_int(json, "keywords.optimization.print_level", config%opt_print_level)
      call optional_logical(json, "keywords.optimization.trajectory", config%opt_trajectory)
      call optional_logical(json, "keywords.optimization.freeze_terms", config%opt_freeze_terms)

      call optional_string(json, "keywords.xtb.solvent", config%solvent)
      call optional_string(json, "keywords.xtb.solvation_model", config%solvation_model)
      ! Before 0.2.0 these were reachable only from a hand-written .mqc deck --
      ! the JSON generator had no field for either, though both are plumbed
      ! through to the xTB method. Retiring .mqc without adding them here
      ! would have quietly removed the only way to turn them off.
      call optional_logical(json, "keywords.xtb.use_cds", config%use_cds)
      call optional_logical(json, "keywords.xtb.use_shift", config%use_shift)
      call optional_real(json, "keywords.xtb.dielectric", config%dielectric)
      call optional_int(json, "keywords.xtb.cpcm_nang", config%cpcm_nang)
      call optional_real(json, "keywords.xtb.cpcm_rscale", config%cpcm_rscale)

      call read_fragmentation(json, config, error)
      if (error%has_error()) return

      ! ---- molecules -------------------------------------------------------
      settings = .false.
      if (present(settings_only)) settings = settings_only

      call json%info("molecules", found=found, n_children=n_mol)
      if (settings) then
         if (found) then
            call error%set(ERROR_VALIDATION, "settings must not define 'molecules': the "// &
                           "system comes from the handle, so a geometry given here would "// &
                           "be read, validated, and then quietly discarded")
         end if
         return
      end if
      if (.not. found .or. n_mol <= 0) then
         call error%set(ERROR_VALIDATION, "Missing or empty required key: molecules")
         return
      end if

      if (n_mol == 1) then
         ! Single-molecule mode: the fields live at the top of the config and
         ! nmol stays 0, which is what every consumer keys off.
         config%nmol = 0
         call read_molecule(json, "molecules(1)", base_dir, config%charge, &
                            config%multiplicity, config%geometry, config%nfrag, &
                            config%fragments, config%uniform_system, config%nbonds, &
                            config%nbroken, config%bonds, error)
         if (error%has_error()) return
      else
         config%nmol = n_mol
         allocate (config%molecules(n_mol))
         do imol = 1, n_mol
            call read_molecule(json, "molecules("//int_to_key(imol)//")", base_dir, &
                               config%molecules(imol)%charge, &
                               config%molecules(imol)%multiplicity, &
                               config%molecules(imol)%geometry, &
                               config%molecules(imol)%nfrag, &
                               config%molecules(imol)%fragments, &
                               config%molecules(imol)%uniform_system, &
                               config%molecules(imol)%nbonds, &
                               config%molecules(imol)%nbroken, &
                               config%molecules(imol)%bonds, error)
            if (error%has_error()) return
         end do
      end if
   end subroutine populate_config

   subroutine read_fragmentation(json, config, error)
      !! The keywords.fragmentation block, including per-level cutoffs
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      integer :: level, highest
      real(dp) :: cutoff
      logical :: found

      call json%info("keywords.fragmentation", found=found)
      if (.not. found) return

      call require_string(json, "keywords.fragmentation.method", config%frag_method, error)
      if (error%has_error()) return
      call optional_int(json, "keywords.fragmentation.level", config%frag_level)
      call optional_logical(json, "keywords.fragmentation.allow_overlapping_fragments", &
                            config%allow_overlapping_fragments)
      call optional_int(json, "keywords.fragmentation.max_intersection_level", &
                        config%max_intersection_level)
      call optional_string(json, "keywords.fragmentation.expansion", config%expansion_kind)
      call optional_string(json, "keywords.fragmentation.counterpoise", config%counterpoise)
      call optional_string(json, "keywords.fragmentation.far_field", config%fmo_far_field)
      call optional_real(json, "keywords.fragmentation.resppc", config%fmo_resppc)
      call optional_int(json, "keywords.fragmentation.max_outer", config%fmo_max_outer)
      call optional_real(json, "keywords.fragmentation.outer_tolerance", config%fmo_tolerance)
      call optional_int(json, "keywords.fragmentation.scf_max_iter", config%fmo_scf_max_iter)
      call optional_real(json, "keywords.fragmentation.scf_energy_tolerance", config%fmo_scf_energy_tol)
      call optional_real(json, "keywords.fragmentation.scf_density_tolerance", config%fmo_scf_density_tol)
      call optional_string(json, "keywords.fragmentation.embedding", config%embedding)
      call optional_string(json, "keywords.fragmentation.cutoff_method", config%cutoff_method)
      call optional_string(json, "keywords.fragmentation.distance_metric", config%distance_metric)
      call optional_int(json, "keywords.fragmentation.global_groups", config%global_groups)
      call optional_int(json, "keywords.fragmentation.nodes_per_group", config%nodes_per_group)

      call read_cutoffs(json, config, error)
   end subroutine read_fragmentation

   subroutine read_cutoffs(json, config, error)
      !! Per-level distance cutoffs from keywords.fragmentation.cutoffs
      !!
      !! Keys are either an n-mer name ("dimer", "trimer", ...) or the level as
      !! a decimal string ("2", "3", ...); both spellings mean the same level
      !! and the JSON may mix them. The result is an array indexed by level and
      !! sized to the highest one present, so `fragment_cutoffs(3)` is the
      !! trimer cutoff whichever spelling produced it.
      !!
      !! Both spellings are probed rather than enumerated because json-fortran
      !! offers no way to list an object's keys through the `json_file`
      !! interface.
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      character(len=*), parameter :: PREFIX = "keywords.fragmentation.cutoffs."
      integer :: level, highest
      real(dp) :: cutoff, previous
      logical :: found
      logical :: present_level(MAX_MBE_LEVEL)

      call json%info("keywords.fragmentation.cutoffs", found=found)
      if (.not. found) return

      highest = 0
      present_level = .false.
      do level = 2, MAX_MBE_LEVEL
         found = .false.
         if (len(nmer_name(level)) > 0) call json%info(PREFIX//nmer_name(level), found=found)
         if (.not. found) call json%info(PREFIX//int_to_key(level), found=found)
         if (found) then
            present_level(level) = .true.
            highest = level
         end if
      end do
      if (highest < 2) return

      ! Sized to every supported level and filled with the same -1 sentinel the
      ! .mqc parser uses, because downstream code reads "no cutoff at this
      ! level" off that sentinel and would take a 0 as "screen everything out".
      allocate (config%fragment_cutoffs(MAX_MBE_LEVEL))
      config%fragment_cutoffs = -1.0_dp

      previous = huge(1.0_dp)
      do level = 2, highest
         if (.not. present_level(level)) cycle
         cutoff = 0.0_dp
         if (len(nmer_name(level)) > 0) call optional_real(json, PREFIX//nmer_name(level), cutoff)
         call optional_real(json, PREFIX//int_to_key(level), cutoff)
         if (cutoff <= 0.0_dp) then
            call error%set(ERROR_VALIDATION, PREFIX//int_to_key(level)// &
                           " must be a positive number")
            return
         end if
         ! A higher-order cutoff above a lower-order one would screen in more
         ! trimers than dimers, which cannot be what anyone meant.
         if (cutoff > previous) then
            call error%set(ERROR_VALIDATION, "keywords.fragmentation.cutoffs: the level-"// &
                           int_to_key(level)//" cutoff exceeds a lower-order one; "// &
                           "cutoffs must decrease with n-mer level")
            return
         end if
         previous = cutoff
         config%fragment_cutoffs(level) = cutoff
      end do
   end subroutine read_cutoffs

   pure function nmer_name(level) result(name)
      !! The spelled-out key for an n-mer level, or "" past the named range
      integer, intent(in) :: level
      character(len=:), allocatable :: name

      select case (level)
      case (2)
         name = "dimer"
      case (3)
         name = "trimer"
      case (4)
         name = "tetramer"
      case (5)
         name = "pentamer"
      case (6)
         name = "hexamer"
      case (7)
         name = "heptamer"
      case (8)
         name = "octamer"
      case default
         name = ""
      end select
   end function nmer_name

   subroutine read_guess_steps(json, config, error)
      !! The basis-set-projection ladder, one entry per preliminary SCF
      !!
      !! Order is the order written: the first step is the cheapest basis and the
      !! last hands its density to the target one from `model.basis`, which is not
      !! listed here. Three SCFs -- STO-3G, 6-31G, cc-pVTZ -- is two steps.
      !!
      !! The schema has already refused a `subscf` that does not belong to a
      !! projection guess and a projection guess with no steps, so this only has
      !! to read what is there.
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: prefix
      integer :: n_steps, i
      logical :: found

      call json%info("keywords.guess.subscf.steps", found=found, n_children=n_steps)
      if (.not. found .or. n_steps < 1) return

      allocate (config%guess_steps(n_steps))
      do i = 1, n_steps
         prefix = "keywords.guess.subscf.steps("//trim(to_char(i))//")"
         call optional_string(json, prefix//".basis", config%guess_steps(i)%basis)
         if (.not. allocated(config%guess_steps(i)%basis)) then
            call error%set(ERROR_VALIDATION, "guess step "//trim(to_char(i))// &
                           " has no basis")
            return
         end if
         call optional_int(json, prefix//".maxiter", config%guess_steps(i)%maxiter)
         call optional_real(json, prefix//".tolerance", config%guess_steps(i)%tolerance)
      end do
   end subroutine read_guess_steps

   subroutine read_molecule(json, prefix, base_dir, charge, multiplicity, geom, &
                            nfrag, fragments, uniform, nbonds, nbroken, bonds, error)
      !! One entry of the molecules array
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: prefix    !! e.g. "molecules(1)"
      character(len=*), intent(in) :: base_dir  !! For resolving a relative xyz path
      integer, intent(out) :: charge, multiplicity
      type(geometry_type), intent(inout) :: geom
      integer, intent(out) :: nfrag
      type(input_fragment_t), allocatable, intent(out) :: fragments(:)
      logical, intent(out) :: uniform
      integer, intent(out) :: nbonds, nbroken
      type(bond_t), allocatable, intent(out) :: bonds(:)
      type(error_t), intent(out) :: error

      charge = 0
      multiplicity = 1
      nfrag = 0
      nbonds = 0
      nbroken = 0

      call read_molecule_geometry(json, prefix, base_dir, geom, error)
      if (error%has_error()) return

      call require_int(json, prefix//".molecular_charge", charge, error)
      if (error%has_error()) return
      call require_int(json, prefix//".molecular_multiplicity", multiplicity, error)
      if (error%has_error()) return
      if (multiplicity <= 0) then
         call error%set(ERROR_VALIDATION, prefix//".molecular_multiplicity must be >= 1")
         return
      end if

      call read_fragments(json, prefix, base_dir, geom%natoms, nfrag, fragments, uniform, error)
      if (error%has_error()) return

      call read_connectivity(json, prefix, geom%natoms, nfrag, fragments, &
                             nbonds, nbroken, bonds, error)
   end subroutine read_molecule

   subroutine read_molecule_geometry(json, prefix, base_dir, geom, error)
      !! Either an xyz file reference or inline symbols plus a flat coordinate list
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: prefix, base_dir
      type(geometry_type), intent(inout) :: geom
      type(error_t), intent(out) :: error

      character(len=:), allocatable :: xyz_rel, symbol
      real(dp), allocatable :: flat(:)
      integer :: natoms, iatom
      logical :: has_xyz, has_symbols, found

      call json%info(prefix//".xyz", found=has_xyz)
      call json%info(prefix//".symbols", found=has_symbols, n_children=natoms)

      if (has_xyz .and. has_symbols) then
         call error%set(ERROR_VALIDATION, prefix// &
                        ": give either xyz or symbols+geometry, not both")
         return
      end if

      if (has_xyz) then
         call require_string(json, prefix//".xyz", xyz_rel, error)
         if (error%has_error()) return
         ! Relative to the input deck, so a deck and its geometries travel
         ! together regardless of where the job is launched from.
         call read_xyz_file(join_path(base_dir, xyz_rel), geom, error)
         if (error%has_error()) call error%add_context(prefix//".xyz")
         return
      end if

      if (.not. has_symbols .or. natoms <= 0) then
         call error%set(ERROR_VALIDATION, prefix// &
                        ": must give either xyz or symbols+geometry")
         return
      end if

      call json%get(prefix//".geometry", flat, found)
      if (.not. found .or. .not. allocated(flat)) then
         call error%set(ERROR_VALIDATION, prefix//".geometry is missing")
         return
      end if
      if (size(flat) /= 3*natoms) then
         call error%set(ERROR_VALIDATION, prefix//".geometry: expected "// &
                        int_to_key(3*natoms)//" values, got "//int_to_key(size(flat)))
         return
      end if

      geom%natoms = natoms
      allocate (character(len=MAX_ELEMENT_SYMBOL_LEN) :: geom%elements(natoms))
      allocate (geom%coords(3, natoms))
      do iatom = 1, natoms
         call require_string(json, prefix//".symbols("//int_to_key(iatom)//")", symbol, error)
         if (error%has_error()) return
         geom%elements(iatom) = trim(adjustl(symbol))
         geom%coords(1:3, iatom) = flat(3*(iatom - 1) + 1:3*iatom)
      end do
   end subroutine read_molecule_geometry

   subroutine read_fragments(json, prefix, base_dir, natoms, nfrag, fragments, uniform, error)
      !! The fragments array, with its parallel charge, multiplicity and potential lists
      !!
      !! `fragment_potentials` is the one of those that may be given as a **single
      !! string rather than a list**, when `uniform_system` says every fragment is the
      !! same species. That is not sugar for its own sake: a ten-thousand-water cluster
      !! would otherwise repeat one filename ten thousand times, and a list that long is
      !! somewhere for a typo to hide rather than information.
      !!
      !! Potential paths are resolved against the deck, exactly as `xyz` is, so a
      !! deck refers to the potential beside it and works from any directory.
      !!
      !! `uniform_system` is checked rather than believed. Fragments claimed uniform
      !! must agree in atom count, so a deck that mixes species and says otherwise is
      !! refused here instead of being handed the wrong potential and run.
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: prefix
      character(len=*), intent(in) :: base_dir  !! For resolving relative potential paths
      integer, intent(in) :: natoms
      integer, intent(out) :: nfrag
      type(input_fragment_t), allocatable, intent(out) :: fragments(:)
      logical, intent(out) :: uniform
      type(error_t), intent(out) :: error

      integer, allocatable :: indices(:)
      character(len=:), allocatable :: frag_path, shared, one
      integer :: ifrag, n_here, n_pot
      logical :: found, has_list

      nfrag = 0
      uniform = .false.
      call json%info(prefix//".fragments", found=found, n_children=nfrag)
      if (.not. found .or. nfrag <= 0) then
         nfrag = 0
         return
      end if

      call optional_logical(json, prefix//".uniform_system", uniform)

      ! A single string broadcasts; a list is read per fragment. Which one was written
      ! is settled by asking for the string first -- `json%get` on a list into a scalar
      ! does not find one, so this distinguishes the two without inspecting types.
      call json%get(prefix//".fragment_potentials", shared, found)
      if (.not. (found .and. allocated(shared))) then
         if (allocated(shared)) deallocate (shared)
      end if
      n_pot = 0
      call json%info(prefix//".fragment_potentials", found=has_list, n_children=n_pot)
      has_list = has_list .and. .not. allocated(shared)
      if (has_list .and. n_pot /= nfrag) then
         call error%set(ERROR_VALIDATION, prefix// &
                        ": 'fragment_potentials' has "//int_to_key(n_pot)// &
                        " entries for "//int_to_key(nfrag)//" fragments")
         return
      end if
      if (allocated(shared) .and. .not. uniform .and. nfrag > 1) then
         call error%set(ERROR_VALIDATION, prefix// &
                        ": 'fragment_potentials' is a single potential for "// &
                        int_to_key(nfrag)//" fragments, which only makes sense with "// &
                        "'uniform_system': true")
         return
      end if

      allocate (fragments(nfrag))
      do ifrag = 1, nfrag
         frag_path = prefix//".fragments("//int_to_key(ifrag)//")"

         call json%get(frag_path, indices, found)
         if (.not. found .or. .not. allocated(indices)) then
            call error%set(ERROR_PARSE, frag_path//" is not a list of atom indices")
            return
         end if
         n_here = size(indices)
         if (n_here == 0) then
            call error%set(ERROR_VALIDATION, frag_path//" is empty")
            return
         end if
         if (any(indices < 0) .or. any(indices >= natoms)) then
            call error%set(ERROR_VALIDATION, frag_path// &
                           ": atom index out of range for natoms="//int_to_key(natoms)// &
                           " (indices are 0-based)")
            return
         end if
         fragments(ifrag)%indices = indices

         ! Charges and multiplicities default per fragment, so a missing list
         ! is neutral singlets rather than an error.
         fragments(ifrag)%charge = 0
         fragments(ifrag)%multiplicity = 1
         call optional_int(json, prefix//".fragment_charges("//int_to_key(ifrag)//")", &
                           fragments(ifrag)%charge)
         call optional_int(json, prefix//".fragment_multiplicities("//int_to_key(ifrag)//")", &
                           fragments(ifrag)%multiplicity)

         ! A fragment with no potential named for it stays quantum, so an absent
         ! entry is not an error -- it is the mixed QM/EFP case.
         if (allocated(shared)) then
            fragments(ifrag)%potential = join_path(base_dir, shared)
         else if (has_list) then
            if (allocated(one)) deallocate (one)
            call json%get(prefix//".fragment_potentials("//int_to_key(ifrag)//")", &
                          one, found)
            if (found .and. allocated(one)) then
               if (len_trim(one) > 0) then
                  fragments(ifrag)%potential = join_path(base_dir, trim(one))
               end if
            end if
         end if
      end do

      ! The uniformity claim, checked. Atom count is what a potential has to agree with
      ! for its own atoms to be superposable onto a fragment's, so it is the cheap
      ! necessary condition; the loader still has to check composition when it reads
      ! the file, which is where the elements become known.
      if (uniform) then
         do ifrag = 2, nfrag
            if (size(fragments(ifrag)%indices) /= size(fragments(1)%indices)) then
               call error%set(ERROR_VALIDATION, prefix// &
                              ": 'uniform_system' is true but fragment "// &
                              int_to_key(ifrag)//" has "// &
                              int_to_key(size(fragments(ifrag)%indices))// &
                              " atoms where the first has "// &
                              int_to_key(size(fragments(1)%indices)))
               return
            end if
         end do
      end if
   end subroutine read_fragments

   subroutine read_connectivity(json, prefix, natoms, nfrag, fragments, &
                                nbonds, nbroken, bonds, error)
      !! The connectivity list, deriving which bonds fragmentation breaks
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: prefix
      integer, intent(in) :: natoms, nfrag
      type(input_fragment_t), allocatable, intent(in) :: fragments(:)
      integer, intent(out) :: nbonds, nbroken
      type(bond_t), allocatable, intent(out) :: bonds(:)
      type(error_t), intent(out) :: error

      integer, allocatable :: triple(:)
      character(len=:), allocatable :: bond_path
      integer :: ibond
      logical :: found

      nbonds = 0
      nbroken = 0
      call json%info(prefix//".connectivity", found=found, n_children=nbonds)
      if (.not. found .or. nbonds <= 0) then
         nbonds = 0
         return
      end if

      allocate (bonds(nbonds))
      do ibond = 1, nbonds
         bond_path = prefix//".connectivity("//int_to_key(ibond)//")"
         call json%get(bond_path, triple, found)
         if (.not. found .or. .not. allocated(triple)) then
            call error%set(ERROR_PARSE, bond_path//" is not a list")
            return
         end if
         if (size(triple) /= 3) then
            call error%set(ERROR_VALIDATION, bond_path//" must be [i, j, order]")
            return
         end if
         if (any(triple(1:2) < 0) .or. any(triple(1:2) >= natoms)) then
            call error%set(ERROR_VALIDATION, bond_path// &
                           ": atom index out of range for natoms="//int_to_key(natoms))
            return
         end if
         if (triple(3) <= 0) then
            call error%set(ERROR_VALIDATION, bond_path//": bond order must be positive")
            return
         end if

         bonds(ibond)%atom_i = triple(1)
         bonds(ibond)%atom_j = triple(2)
         bonds(ibond)%order = triple(3)
         bonds(ibond)%is_broken = bond_is_broken(triple(1), triple(2), nfrag, fragments)
         if (bonds(ibond)%is_broken) nbroken = nbroken + 1
      end do
   end subroutine read_connectivity

   pure function bond_is_broken(atom_i, atom_j, nfrag, fragments) result(broken)
      !! Whether fragmentation severs the bond between two atoms
      !!
      !! Broken means the two atoms do not belong to the same set of fragments.
      !! Comparing *sets* rather than asking "are they both in some fragment"
      !! is what makes this correct for overlapping fragments, where an atom
      !! can legitimately be in more than one.
      integer, intent(in) :: atom_i, atom_j
      integer, intent(in) :: nfrag
      type(input_fragment_t), allocatable, intent(in) :: fragments(:)
      logical :: broken

      integer :: ifrag
      logical :: in_i, in_j

      broken = .false.
      if (nfrag <= 0 .or. .not. allocated(fragments)) return

      do ifrag = 1, nfrag
         in_i = any(fragments(ifrag)%indices == atom_i)
         in_j = any(fragments(ifrag)%indices == atom_j)
         if (in_i .neqv. in_j) then
            broken = .true.
            return
         end if
      end do
   end function bond_is_broken

   ! ==========================================================================
   !  Accessors
   !
   !  json-fortran reports a missing key through `found` rather than failing,
   !  so "required" and "optional" differ only in what happens when it is
   !  false. The optional forms leave the target at whatever default the config
   !  type declared, which is how the defaults stay in one place.
   ! ==========================================================================

   subroutine require_string(json, path, value, error)
      !! Fetch a string, or record a missing-key error
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      character(len=:), allocatable, intent(out) :: value
      type(error_t), intent(out) :: error

      logical :: found

      call json%get(path, value, found)
      if (.not. found .or. .not. allocated(value)) then
         call error%set(ERROR_VALIDATION, "Missing required key: "//path)
      end if
   end subroutine require_string

   subroutine require_int(json, path, value, error)
      !! Fetch an integer, or record a missing-key error
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      integer, intent(out) :: value
      type(error_t), intent(out) :: error

      logical :: found

      call json%get(path, value, found)
      if (.not. found) call error%set(ERROR_VALIDATION, "Missing required key: "//path)
   end subroutine require_int

   subroutine optional_string(json, path, value)
      !! Fetch a string if present, leaving `value` untouched otherwise
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      character(len=:), allocatable, intent(inout) :: value

      character(len=:), allocatable :: text
      logical :: found

      call json%get(path, text, found)
      if (found .and. allocated(text)) value = text
   end subroutine optional_string

   subroutine optional_int(json, path, value)
      !! Fetch an integer if present, leaving `value` at its default otherwise
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      integer, intent(inout) :: value

      integer :: found_value
      logical :: found

      call json%get(path, found_value, found)
      if (found) value = found_value
   end subroutine optional_int

   subroutine optional_real(json, path, value)
      !! Fetch a real if present, leaving `value` at its default otherwise
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      real(dp), intent(inout) :: value

      real(dp) :: found_value
      logical :: found

      call json%get(path, found_value, found)
      if (found) value = found_value
   end subroutine optional_real

   subroutine resolve_gpu_request(config, backend_named, error)
      !! Settle `system.gpu` against `backend`, then refuse a request nothing can honour
      !!
      !! Three things happen here, and all three happen before a calculation is
      !! set up rather than inside one. A deck that asks for a GPU it cannot have
      !! should be told so while it is still a deck: the alternative is finding
      !! out after fragmentation, once per fragment, from a bridge that computes
      !! nothing.
      !!
      !!   * `system.gpu` resolves into `backend`, so everything downstream reads
      !!     one field however the deck spelled the question. Naming both and
      !!     disagreeing -- `"gpu": true` beside `"backend": "cpu"` -- is refused
      !!     rather than resolved by precedence.
      !!   * asking for cuEST on a build without it is refused. This is asked of
      !!     the bridge that linked, not of a preprocessor symbol, so the answer
      !!     describes the binary rather than somebody's configure line.
      !!   * asking for cuEST for a method it has no implementation of is refused.
      !!     Substituting the CPU silently would report a provenance and a set of
      !!     timings that were not true.
      type(mqc_config_t), intent(inout) :: config
      logical, intent(in) :: backend_named  !! Whether the deck had a root `backend` key
      type(error_t), intent(out) :: error

      integer :: backend_kind

      if (config%gpu_set) then
         if (backend_named) then
            call parse_backend_name(config%backend, backend_kind, error)
            if (error%has_error()) return
            if ((backend_kind == BACKEND_CUEST) .neqv. config%gpu) then
               call error%set(ERROR_VALIDATION, "system.gpu and backend disagree: "// &
                              "system.gpu is "//trim(merge("true ", "false", config%gpu))// &
                              " but backend is '"//trim(config%backend)// &
                              "'. Name one of them, not both.")
               return
            end if
         end if
         if (config%gpu) then
            config%backend = "cuest"
         else
            config%backend = "libcint"
         end if
      end if

      call parse_backend_name(config%backend, backend_kind, error)
      if (error%has_error()) return
      if (backend_kind /= BACKEND_CUEST) return

      if (.not. cuest_backend_available()) then
         call error%set(ERROR_VALIDATION, gpu_asked_by(config)//" but this build has "// &
                        "no cuEST backend. Rebuild with CMake and "// &
                        "-DMQC_ENABLE_CUEST=ON -DCUEST_ROOT=<cuest archive>, or drop "// &
                        "the request to run on the CPU.")
         return
      end if

      if (.not. method_runs_on_cuest(config%method)) then
         call error%set(ERROR_VALIDATION, gpu_asked_by(config)//" but '"// &
                        trim(method_type_to_string(config%method))//"' has no cuEST "// &
                        "implementation; only Hartree-Fock and DFT are offloaded. "// &
                        "Run it on the CPU, or pick a method the GPU has.")
         return
      end if
   end subroutine resolve_gpu_request

   pure function gpu_asked_by(config) result(text)
      !! Which key asked for the GPU, so the refusal names what the deck wrote
      type(mqc_config_t), intent(in) :: config
      character(len=:), allocatable :: text

      if (config%gpu_set) then
         text = "system.gpu is true"
      else
         text = "backend '"//trim(config%backend)//"' was asked for"
      end if
   end function gpu_asked_by

   subroutine optional_logical_seen(json, path, value, seen)
      !! Fetch a boolean, and say whether the deck named it
      !!
      !! Needed wherever a default comes from somewhere other than the field's
      !! initialiser -- the method name, say -- because "absent" and "false" have
      !! to be told apart before the name can be allowed to lose to a keyword.
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      logical, intent(inout) :: value
      logical, intent(out) :: seen

      logical :: found_value

      call json%get(path, found_value, seen)
      if (seen) value = found_value
   end subroutine optional_logical_seen

   subroutine optional_logical(json, path, value)
      !! Fetch a boolean if present, leaving `value` at its default otherwise
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      logical, intent(inout) :: value

      logical :: found_value
      logical :: found

      call json%get(path, found_value, found)
      if (found) value = found_value
   end subroutine optional_logical

   ! ==========================================================================
   !  Small string utilities
   ! ==========================================================================

   pure function int_to_key(value) result(key)
      !! An integer as the string json-fortran's paths and keys want
      integer, intent(in) :: value
      character(len=:), allocatable :: key
      character(len=16) :: buffer

      write (buffer, "(I0)") value
      key = trim(buffer)
   end function int_to_key

   pure function directory_of(path) result(directory)
      !! The directory part of a path, or "." when there is none
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: directory

      integer :: slash

      slash = index(path, "/", back=.true.)
      if (slash > 0) then
         directory = path(1:slash - 1)
      else
         directory = "."
      end if
   end function directory_of

   pure function join_path(directory, relative) result(joined)
      !! Join a directory and a relative path
      !!
      !! An absolute `relative` wins, so a deck can point at a geometry
      !! somewhere else entirely.
      character(len=*), intent(in) :: directory, relative
      character(len=:), allocatable :: joined

      if (len_trim(relative) > 0) then
         if (relative(1:1) == "/") then
            joined = trim(relative)
            return
         end if
      end if
      joined = trim(directory)//"/"//trim(relative)
   end function join_path

end module mqc_json_config_reader
