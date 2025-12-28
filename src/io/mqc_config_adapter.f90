!! Adapter module to convert mqc_config_t to internal driver structures
!! This module provides a bridge between the .mqc file format and the driver backend
module mqc_config_adapter
   !! Provides conversion utilities from mqc_config_t to driver-compatible structures
   use pic_types, only: dp, int32
   use mqc_config_parser, only: mqc_config_t
   use mqc_physical_fragment, only: system_geometry_t, to_bohr
   use mqc_elements, only: element_symbol_to_number
   implicit none
   private

   public :: driver_config_t  !! Minimal config for driver
   public :: config_to_driver, config_to_system_geometry
   public :: get_logger_level  !! Convert log level string to integer

   !! Minimal configuration for driver (internal use only)
   type :: driver_config_t
      integer(int32) :: method      !! QC method constant
      integer(int32) :: calc_type   !! Calculation type constant
      integer :: nlevel = 0         !! Fragmentation level (0 = unfragmented)
   end type driver_config_t

contains

   subroutine config_to_driver(mqc_config, driver_config)
      !! Convert mqc_config_t to minimal driver_config_t
      !! Extracts only the fields needed by the driver
      type(mqc_config_t), intent(in) :: mqc_config
      type(driver_config_t), intent(out) :: driver_config

      ! Copy method and calc_type (already integers)
      driver_config%method = mqc_config%method
      driver_config%calc_type = mqc_config%calc_type

      ! Set fragmentation level
      ! For unfragmented calculations (nfrag=0), nlevel must be 0
      if (mqc_config%nfrag == 0) then
         driver_config%nlevel = 0
      else
         driver_config%nlevel = mqc_config%frag_level
      end if

   end subroutine config_to_driver

   subroutine config_to_system_geometry(mqc_config, sys_geom, stat, errmsg)
      !! Convert mqc_config_t geometry to system_geometry_t
      !! For unfragmented calculations (nfrag=0), treats entire system as single unit
      !! For fragmented calculations, currently assumes monomer-based fragmentation
      type(mqc_config_t), intent(in) :: mqc_config
      type(system_geometry_t), intent(out) :: sys_geom
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      integer :: i
      logical :: use_angstrom

      stat = 0

      ! Check if geometry is loaded
      if (mqc_config%geometry%natoms == 0) then
         stat = 1
         errmsg = "No geometry loaded in mqc_config"
         return
      end if

      ! Determine units
      use_angstrom = .true.
      if (allocated(mqc_config%units)) then
         if (trim(mqc_config%units) == 'bohr') then
            use_angstrom = .false.
         end if
      end if

      if (mqc_config%nfrag == 0) then
         ! Unfragmented calculation: entire system is one "monomer"
         call geometry_to_system_unfragmented(mqc_config%geometry, sys_geom, use_angstrom)
      else
         ! Fragmented calculation with explicit fragments
         ! For now, we'll try to infer monomer size from first fragment
         ! TODO: This is a temporary solution - need better fragmentation support
         call geometry_to_system_fragmented(mqc_config, sys_geom, use_angstrom, stat, errmsg)
         if (stat /= 0) return
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

   subroutine geometry_to_system_fragmented(mqc_config, sys_geom, use_angstrom, stat, errmsg)
      !! Convert geometry to system_geometry_t for fragmented calculation
      !! Supports both identical and variable-sized fragments
      type(mqc_config_t), intent(in) :: mqc_config
      type(system_geometry_t), intent(out) :: sys_geom
      logical, intent(in) :: use_angstrom
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      integer :: i, j, atoms_in_first_frag, max_frag_size
      logical :: all_same_size

      stat = 0

      ! Set up basic system geometry
      sys_geom%n_monomers = mqc_config%nfrag
      sys_geom%total_atoms = mqc_config%geometry%natoms

      ! Allocate fragment info arrays
      allocate (sys_geom%fragment_sizes(mqc_config%nfrag))

      ! Get fragment sizes
      max_frag_size = 0
      atoms_in_first_frag = size(mqc_config%fragments(1)%indices)
      all_same_size = .true.

      do i = 1, mqc_config%nfrag
         sys_geom%fragment_sizes(i) = size(mqc_config%fragments(i)%indices)
         max_frag_size = max(max_frag_size, sys_geom%fragment_sizes(i))
         if (sys_geom%fragment_sizes(i) /= atoms_in_first_frag) then
            all_same_size = .false.
         end if
      end do

      ! Allocate fragment_atoms array
      allocate (sys_geom%fragment_atoms(max_frag_size, mqc_config%nfrag))
      sys_geom%fragment_atoms = -1  ! Initialize with invalid index

      ! Store fragment atom indices (0-indexed from input file)
      do i = 1, mqc_config%nfrag
         do j = 1, sys_geom%fragment_sizes(i)
            sys_geom%fragment_atoms(j, i) = mqc_config%fragments(i)%indices(j)
         end do
      end do

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
         sys_geom%element_numbers(i) = element_symbol_to_number(mqc_config%geometry%elements(i))
      end do

      ! Store coordinates (convert to Bohr if needed)
      if (use_angstrom) then
         sys_geom%coordinates = to_bohr(mqc_config%geometry%coords)
      else
         sys_geom%coordinates = mqc_config%geometry%coords
      end if

   end subroutine geometry_to_system_fragmented

   subroutine geometry_to_system(geom, sys_geom, use_angstrom)
      !! Simple conversion for backward compatibility
      !! Deprecated: use config_to_system_geometry instead
      use mqc_geometry, only: geometry_type

      type(geometry_type), intent(in) :: geom
      type(system_geometry_t), intent(out) :: sys_geom
      logical, intent(in) :: use_angstrom

      call geometry_to_system_unfragmented(geom, sys_geom, use_angstrom)

   end subroutine geometry_to_system

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
         level_int = info_level
      end select
   end function get_logger_level

end module mqc_config_adapter
