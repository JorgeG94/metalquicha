!! Adapter module to convert mqc_config_t to internal driver structures
!! This module provides a bridge between the .mqc file format and the driver backend
module mqc_config_adapter
   !! Provides conversion utilities from mqc_config_t to driver-compatible structures
   use pic_types, only: dp, int32
   use mqc_config_parser, only: mqc_config_t
   use mqc_physical_fragment, only: system_geometry_t, to_bohr
   use mqc_elements, only: element_symbol_to_number
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_calculation_keywords, only: hessian_keywords_t, aimd_keywords_t, scf_keywords_t
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: driver_config_t  !! Minimal config for driver
   public :: config_to_driver, config_to_system_geometry
   public :: get_logger_level  !! Convert log level string to integer
   public :: check_fragment_overlap  !! Check for overlapping fragments (for testing)

   !! Runtime configuration for driver (internal use only)
   type :: driver_config_t
      ! Core calculation settings
      integer(int32) :: method      !! QC method constant
      integer(int32) :: calc_type   !! Calculation type constant

      ! Fragmentation settings
      integer :: nlevel = 0         !! Fragmentation level (0 = unfragmented)
      logical :: allow_overlapping_fragments = .false.  !! Enable GMBE for overlapping fragments
      integer :: max_intersection_level = 999  !! Maximum k-way intersection depth for GMBE (default: no limit)
      real(dp), allocatable :: fragment_cutoffs(:)  !! Distance cutoffs for n-mer screening (Angstrom)

      ! XTB solvation settings
      character(len=:), allocatable :: solvent  !! Solvent name or empty for gas phase
      character(len=:), allocatable :: solvation_model  !! "alpb" (default), "gbsa", or "cpcm"
      logical :: use_cds = .true.               !! Include CDS non-polar terms (not for CPCM)
      logical :: use_shift = .true.             !! Include solution state shift (not for CPCM)
      ! CPCM-specific settings
      real(dp) :: dielectric = -1.0_dp          !! Direct dielectric constant (-1 = use solvent lookup)
      integer :: cpcm_nang = 110                !! Number of angular grid points for CPCM
      real(dp) :: cpcm_rscale = 1.0_dp          !! Radii scaling factor for CPCM

      ! Calculation-specific keywords (structured)
      type(hessian_keywords_t) :: hessian  !! Hessian calculation keywords
      type(aimd_keywords_t) :: aimd        !! AIMD calculation keywords
      type(scf_keywords_t) :: scf          !! SCF calculation keywords
   end type driver_config_t

contains

   subroutine config_to_driver(mqc_config, driver_config, molecule_index)
      !! Convert mqc_config_t to minimal driver_config_t
      !! Extracts only the fields needed by the driver
      !! If molecule_index is provided, uses that molecule's fragment count
      type(mqc_config_t), intent(in) :: mqc_config
      type(driver_config_t), intent(out) :: driver_config
      integer, intent(in), optional :: molecule_index  !! Which molecule to use (for multi-molecule mode)

      integer :: nfrag_to_use

      ! Copy method and calc_type (already integers)
      driver_config%method = mqc_config%method
      driver_config%calc_type = mqc_config%calc_type

      ! Determine fragment count
      if (present(molecule_index)) then
         ! Multi-molecule mode: use specific molecule's fragment count
         if (molecule_index < 1 .or. molecule_index > mqc_config%nmol) then
            nfrag_to_use = 0
         else
            nfrag_to_use = mqc_config%molecules(molecule_index)%nfrag
         end if
      else
         ! Single molecule mode (backward compatible)
         nfrag_to_use = mqc_config%nfrag
      end if

      ! Set fragmentation level
      ! For unfragmented calculations (nfrag=0), nlevel must be 0
      if (nfrag_to_use == 0) then
         driver_config%nlevel = 0
      else
         driver_config%nlevel = mqc_config%frag_level
      end if

      ! Set GMBE overlapping fragments flag
      driver_config%allow_overlapping_fragments = mqc_config%allow_overlapping_fragments

      ! Set GMBE maximum intersection level
      driver_config%max_intersection_level = mqc_config%max_intersection_level

      ! Copy fragment distance cutoffs if present
      if (allocated(mqc_config%fragment_cutoffs)) then
         allocate (driver_config%fragment_cutoffs(size(mqc_config%fragment_cutoffs)))
         driver_config%fragment_cutoffs = mqc_config%fragment_cutoffs
      end if

      ! Copy XTB solvation settings
      if (allocated(mqc_config%solvent)) then
         driver_config%solvent = mqc_config%solvent
      end if
      if (allocated(mqc_config%solvation_model)) then
         driver_config%solvation_model = mqc_config%solvation_model
      end if
      driver_config%use_cds = mqc_config%use_cds
      driver_config%use_shift = mqc_config%use_shift
      ! Copy CPCM-specific settings
      driver_config%dielectric = mqc_config%dielectric
      driver_config%cpcm_nang = mqc_config%cpcm_nang
      driver_config%cpcm_rscale = mqc_config%cpcm_rscale

      ! Set calculation-specific keywords
      driver_config%hessian%displacement = mqc_config%hessian_displacement
      driver_config%hessian%temperature = mqc_config%hessian_temperature
      driver_config%hessian%pressure = mqc_config%hessian_pressure
      driver_config%aimd%dt = mqc_config%aimd_dt
      driver_config%aimd%nsteps = mqc_config%aimd_nsteps
      driver_config%aimd%initial_temperature = mqc_config%aimd_initial_temperature
      driver_config%aimd%output_frequency = mqc_config%aimd_output_frequency
      driver_config%scf%max_iterations = mqc_config%scf_maxiter
      driver_config%scf%convergence_threshold = mqc_config%scf_tolerance

   end subroutine config_to_driver

   subroutine config_to_system_geometry(mqc_config, sys_geom, error, molecule_index)
      !! Convert mqc_config_t geometry to system_geometry_t
      !! For unfragmented calculations (nfrag=0), treats entire system as single unit
      !! For fragmented calculations, currently assumes monomer-based fragmentation
      !! If molecule_index is provided, uses that specific molecule from multi-molecule mode
      type(mqc_config_t), intent(in) :: mqc_config
      type(system_geometry_t), intent(out) :: sys_geom
      type(error_t), intent(out) :: error
      integer, intent(in), optional :: molecule_index  !! Which molecule to use (for multi-molecule mode)

      integer :: i
      logical :: use_angstrom

      ! Determine units
      use_angstrom = .true.
      if (allocated(mqc_config%units)) then
         if (trim(mqc_config%units) == 'bohr') then
            use_angstrom = .false.
         end if
      end if

      ! Handle multi-molecule vs single molecule mode
      if (present(molecule_index)) then
         ! Multi-molecule mode: extract specific molecule
         if (molecule_index < 1 .or. molecule_index > mqc_config%nmol) then
            call error%set(ERROR_VALIDATION, "Invalid molecule_index in multi-molecule mode")
            return
         end if
         call molecule_to_system_geometry(mqc_config%molecules(molecule_index), &
                                          sys_geom, use_angstrom, mqc_config%allow_overlapping_fragments, error)
      else
         ! Single molecule mode (backward compatible)
         ! Check if geometry is loaded
         if (mqc_config%geometry%natoms == 0) then
            call error%set(ERROR_VALIDATION, "No geometry loaded in mqc_config")
            return
         end if

         if (mqc_config%nfrag == 0) then
            ! Unfragmented calculation: entire system is one "monomer"
            call geometry_to_system_unfragmented(mqc_config%geometry, sys_geom, use_angstrom)
            sys_geom%charge = mqc_config%charge
            sys_geom%multiplicity = mqc_config%multiplicity
         else
            ! Fragmented calculation with explicit fragments
            call geometry_to_system_fragmented(mqc_config, sys_geom, use_angstrom, error)
            if (error%has_error()) then
               call error%add_context("mqc_config_adapter:config_to_system_geometry")
               return
            end if
         end if
      end if

   end subroutine config_to_system_geometry

   subroutine geometry_to_system_unfragmented(geom, sys_geom, use_angstrom)
      !! Convert geometry to system_geometry_t for unfragmented calculation
      !! Treats entire system as a single monomer
      use mqc_geometry, only: geometry_type

      type(geometry_type), intent(in) :: geom
      type(system_geometry_t), intent(out) :: sys_geom
      logical, intent(in) :: use_angstrom

      integer :: i

      ! For unfragmented: n_monomers=1, atoms_per_monomer=natoms
      sys_geom%n_monomers = 1
      sys_geom%atoms_per_monomer = geom%natoms
      sys_geom%total_atoms = geom%natoms

      allocate (sys_geom%element_numbers(sys_geom%total_atoms))
      allocate (sys_geom%coordinates(3, sys_geom%total_atoms))

      ! Convert element symbols to atomic numbers
      do i = 1, sys_geom%total_atoms
         sys_geom%element_numbers(i) = element_symbol_to_number(geom%elements(i))
      end do

      ! Store coordinates (convert to Bohr if needed)
      if (use_angstrom) then
         sys_geom%coordinates = to_bohr(geom%coords)
      else
         sys_geom%coordinates = geom%coords
      end if

   end subroutine geometry_to_system_unfragmented

   subroutine initialize_fragmented_system(nfrag, geom, fragments, charge, multiplicity, &
                                           allow_overlapping, use_angstrom, sys_geom, error)
      !! Shared helper to initialize system_geometry_t for fragmented calculations
      !! Handles fragment allocation, size checking, and overlap validation
      use mqc_geometry, only: geometry_type
      use mqc_config_parser, only: input_fragment_t

      integer, intent(in) :: nfrag
      type(geometry_type), intent(in) :: geom
      type(input_fragment_t), intent(in) :: fragments(:)
      integer, intent(in) :: charge, multiplicity
      logical, intent(in) :: allow_overlapping
      logical, intent(in) :: use_angstrom
      type(system_geometry_t), intent(out) :: sys_geom
      type(error_t), intent(out) :: error

      integer :: i, j, atoms_in_first_frag, max_frag_size
      logical :: all_same_size

      ! Set up basic system geometry
      sys_geom%n_monomers = nfrag
      sys_geom%total_atoms = geom%natoms
      sys_geom%charge = charge
      sys_geom%multiplicity = multiplicity

      ! Allocate fragment info arrays
      allocate (sys_geom%fragment_sizes(nfrag))
      allocate (sys_geom%fragment_charges(nfrag))
      allocate (sys_geom%fragment_multiplicities(nfrag))

      ! Get fragment sizes
      max_frag_size = 0
      atoms_in_first_frag = size(fragments(1)%indices)
      all_same_size = .true.

      do i = 1, nfrag
         sys_geom%fragment_sizes(i) = size(fragments(i)%indices)
         sys_geom%fragment_charges(i) = fragments(i)%charge
         sys_geom%fragment_multiplicities(i) = fragments(i)%multiplicity
         max_frag_size = max(max_frag_size, sys_geom%fragment_sizes(i))
         if (sys_geom%fragment_sizes(i) /= atoms_in_first_frag) then
            all_same_size = .false.
         end if
      end do

      ! Allocate fragment_atoms array
      allocate (sys_geom%fragment_atoms(max_frag_size, nfrag))
      sys_geom%fragment_atoms = -1  ! Initialize with invalid index

      ! Store fragment atom indices (0-indexed from input file)
      do i = 1, nfrag
         do j = 1, sys_geom%fragment_sizes(i)
            sys_geom%fragment_atoms(j, i) = fragments(i)%indices(j)
         end do
      end do

      ! Check for overlapping fragments if not allowed
      if (.not. allow_overlapping) then
         call check_fragment_overlap(fragments, nfrag, error)
         if (error%has_error()) then
            call error%add_context("mqc_config_adapter:geometry_to_system_fragmented")
            return
         end if
      end if

      ! Set atoms_per_monomer: use common size if identical, else 0
      if (all_same_size) then
         sys_geom%atoms_per_monomer = atoms_in_first_frag
      else
         sys_geom%atoms_per_monomer = 0  ! Signal variable-sized fragments
      end if

      allocate (sys_geom%element_numbers(sys_geom%total_atoms))
      allocate (sys_geom%coordinates(3, sys_geom%total_atoms))

      ! Convert element symbols to atomic numbers
      do i = 1, sys_geom%total_atoms
         sys_geom%element_numbers(i) = element_symbol_to_number(geom%elements(i))
      end do

      ! Store coordinates (convert to Bohr if needed)
      if (use_angstrom) then
         sys_geom%coordinates = to_bohr(geom%coords)
      else
         sys_geom%coordinates = geom%coords
      end if

   end subroutine initialize_fragmented_system

   subroutine geometry_to_system_fragmented(mqc_config, sys_geom, use_angstrom, error)
      !! Convert geometry to system_geometry_t for fragmented calculation
      !! Supports both identical and variable-sized fragments
      type(mqc_config_t), intent(in) :: mqc_config
      type(system_geometry_t), intent(out) :: sys_geom
      logical, intent(in) :: use_angstrom
      type(error_t), intent(out) :: error

      call initialize_fragmented_system(mqc_config%nfrag, mqc_config%geometry, mqc_config%fragments, &
                                        mqc_config%charge, mqc_config%multiplicity, &
                                        mqc_config%allow_overlapping_fragments, use_angstrom, &
                                        sys_geom, error)

   end subroutine geometry_to_system_fragmented

   subroutine molecule_to_system_geometry(mol, sys_geom, use_angstrom, allow_overlapping, error)
      !! Convert a molecule_t to system_geometry_t
      !! Handles both unfragmented (nfrag=0) and fragmented molecules
      use mqc_config_parser, only: molecule_t

      type(molecule_t), intent(in) :: mol
      type(system_geometry_t), intent(out) :: sys_geom
      logical, intent(in) :: use_angstrom
      type(error_t), intent(out) :: error
      logical, intent(in) :: allow_overlapping

      ! Check if geometry is loaded
      if (mol%geometry%natoms == 0) then
         call error%set(ERROR_VALIDATION, "No geometry loaded in molecule")
         return
      end if

      if (mol%nfrag == 0) then
         ! Unfragmented molecule
         call geometry_to_system_unfragmented(mol%geometry, sys_geom, use_angstrom)
         sys_geom%charge = mol%charge
         sys_geom%multiplicity = mol%multiplicity
      else
         ! Fragmented molecule
         call initialize_fragmented_system(mol%nfrag, mol%geometry, mol%fragments, &
                                           mol%charge, mol%multiplicity, &
                                           allow_overlapping, use_angstrom, &
                                           sys_geom, error)
      end if

   end subroutine molecule_to_system_geometry

   function get_logger_level(level_string) result(level_int)
      !! Convert string log level to integer value
      !! This function uses the pic_logger constants
      use pic_logger, only: debug_level, verbose_level, info_level, performance_level, &
                            warning_level, error_level, knowledge_level
      character(len=*), intent(in) :: level_string
      integer :: level_int

      select case (trim(adjustl(level_string)))
      case ('debug', 'Debug', 'DEBUG')
         level_int = debug_level
      case ('verbose', 'Verbose', 'VERBOSE')
         level_int = verbose_level
      case ('info', 'Info', 'INFO')
         level_int = info_level
      case ('performance', 'Performance', 'PERFORMANCE')
         level_int = performance_level
      case ('warning', 'Warning', 'WARNING')
         level_int = warning_level
      case ('error', 'Error', 'ERROR')
         level_int = error_level
      case ('knowledge', 'Knowledge', 'KNOWLEDGE')
         level_int = knowledge_level
      case default
         ! Default to info level if unknown
         call logger%warning("Unknown log level string: "//level_string//". Defaulting to INFO level.")
         level_int = info_level
      end select
   end function get_logger_level

   subroutine check_fragment_overlap(fragments, nfrag, error)
      !! Check if any atoms appear in multiple fragments
      !! This is O(nfrag * natoms_per_frag^2) which is acceptable for typical fragment sizes
      use mqc_config_parser, only: input_fragment_t
      use pic_io, only: to_char

      type(input_fragment_t), intent(in) :: fragments(:)
      integer, intent(in) :: nfrag
      type(error_t), intent(out) :: error

      integer :: i, j, k, l
      integer :: atom_i, atom_j

      ! Compare each pair of fragments
      do i = 1, nfrag - 1
         do j = i + 1, nfrag
            ! Compare atoms in fragment i with atoms in fragment j
            do k = 1, size(fragments(i)%indices)
               atom_i = fragments(i)%indices(k)
               do l = 1, size(fragments(j)%indices)
                  atom_j = fragments(j)%indices(l)
                  if (atom_i == atom_j) then
                     ! Found overlapping atom
                    call error%set(ERROR_VALIDATION, "Overlapping fragments detected: fragments "//to_char(i)//" and "// &
                                    to_char(j)//" both contain atom "//to_char(atom_i)// &
                                    ". Set allow_overlapping_fragments = true to allow this.")
                     return
                  end if
               end do
            end do
         end do
      end do

   end subroutine check_fragment_overlap

end module mqc_config_adapter
