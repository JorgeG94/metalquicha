!! Configuration types an input file parses into
module mqc_config_types
   !! `mqc_config_t` and the pieces it is made of: one molecule, one input
   !! fragment, and the bonds between them.
   !!
   !! These live apart from any reader on purpose. Nearly everything that
   !! touches a config only wants the *type* -- the driver, the adapter, the
   !! MBE engine and the tests all import `mqc_config_t` and never parse
   !! anything -- so tying the type to one file format would make every
   !! consumer depend on that format's parser. It also means a second reader
   !! is a new module rather than a change to an existing one.
   !!
   !! Defaults are declared here, once, and a reader that finds no value for a
   !! key simply leaves the component alone. That is the whole mechanism: there
   !! is no separate table of defaults to drift out of step.
   use pic_types, only: dp, int32
   use mqc_method_types, only: METHOD_TYPE_GFN2
   use mqc_calc_types, only: CALC_TYPE_ENERGY
   use mqc_geometry, only: geometry_type
   use mqc_physical_fragment, only: bond_t
   use mqc_calculation_defaults, only: DEFAULT_DISPLACEMENT, DEFAULT_TEMPERATURE, &
                                       DEFAULT_PRESSURE, DEFAULT_AIMD_DT, &
                                       DEFAULT_AIMD_NSTEPS, DEFAULT_AIMD_TEMPERATURE, &
                                       DEFAULT_AIMD_OUTPUT_FREQ, DEFAULT_SCF_CONV, &
                                       DEFAULT_CPCM_NANG, DEFAULT_CPCM_RSCALE, &
                                       DEFAULT_FRAG_LEVEL, DEFAULT_MAX_INTERSECTION
   implicit none
   private

   public :: mqc_config_t       !! Everything one input file specifies
   public :: molecule_t         !! One molecule of a multi-molecule input
   public :: input_fragment_t   !! One fragment definition, as written in the input
   public :: bond_t             !! Re-exported so consumers need only this module

   type :: input_fragment_t
      !! Input fragment definition with charge, multiplicity, and atom indices
      !! This is the parsed representation from the input file, not the computational fragment
      integer :: charge = 0
      integer :: multiplicity = 1
      integer, allocatable :: indices(:)  !! Atom indices in this fragment
   contains
      procedure :: destroy => input_fragment_destroy
   end type input_fragment_t

   type :: molecule_t
      !! Single molecule definition with structure, geometry, fragments, and connectivity
      character(len=:), allocatable :: name  !! Optional molecule name

      ! Structure information
      integer :: charge = 0
      integer :: multiplicity = 1

      ! Geometry
      type(geometry_type) :: geometry

      ! Fragments
      integer :: nfrag = 0
      type(input_fragment_t), allocatable :: fragments(:)

      ! Connectivity
      integer :: nbonds = 0
      integer :: nbroken = 0
      type(bond_t), allocatable :: bonds(:)

   contains
      procedure :: destroy => molecule_destroy
   end type molecule_t

   type :: mqc_config_t
      !! Complete configuration from .mqc file

      ! Schema information
      character(len=:), allocatable :: schema_name
      character(len=:), allocatable :: schema_version
      integer :: index_base = 0  !! 0-based or 1-based indexing
      character(len=:), allocatable :: units  !! angstrom or bohr

      ! Model information
      integer(int32) :: method = METHOD_TYPE_GFN2
      character(len=:), allocatable :: basis
      character(len=:), allocatable :: aux_basis
      character(len=:), allocatable :: functional
      character(len=:), allocatable :: scf_guess
         !! Initial guess name from %scf
      logical :: scf_unrestricted = .false.
         !! Force UHF/UKS even when the shell is closed
         !! XC functional name, only meaningful when method = dft

      ! XTB solvation settings
      character(len=:), allocatable :: solvent  !! Solvent name (e.g., "water", "ethanol") or empty for gas phase
      character(len=:), allocatable :: solvation_model  !! Solvation model: "alpb" (default), "gbsa", or "cpcm"
      logical :: use_cds = .true.               !! Include non-polar CDS terms in solvation (not for CPCM)
      logical :: use_shift = .true.             !! Include solution state shift in solvation (not for CPCM)
      ! CPCM-specific settings
      real(dp) :: dielectric = -1.0_dp              !! Direct dielectric constant (-1 = use solvent lookup)
      integer :: cpcm_nang = DEFAULT_CPCM_NANG      !! Number of angular grid points for CPCM cavity
      real(dp) :: cpcm_rscale = DEFAULT_CPCM_RSCALE  !! Radii scaling factor for CPCM cavity

      ! Driver information
      integer(int32) :: calc_type = CALC_TYPE_ENERGY

      ! Multiple molecules support
      integer :: nmol = 0  !! Number of molecules (0 = single molecule mode for backward compatibility)
      type(molecule_t), allocatable :: molecules(:)  !! Array of molecules (if nmol > 0)

      ! Single molecule fields (backward compatibility - used if nmol == 0)
      ! Structure information
      integer :: charge = 0
      integer :: multiplicity = 1

      ! Geometry
      type(geometry_type) :: geometry

      ! Fragments
      integer :: nfrag = 0
      type(input_fragment_t), allocatable :: fragments(:)

      ! Connectivity
      integer :: nbonds = 0
      integer :: nbroken = 0
      type(bond_t), allocatable :: bonds(:)

      ! SCF settings
      integer :: scf_maxiter = 300              !! Using 300 (parser-specific, different from DEFAULT_SCF_MAXITER)
      real(dp) :: scf_tolerance = DEFAULT_SCF_CONV

      ! Hessian settings
      real(dp) :: hessian_displacement = DEFAULT_DISPLACEMENT  !! Finite difference displacement (Bohr)
      real(dp) :: hessian_temperature = DEFAULT_TEMPERATURE    !! Temperature for thermochemistry (K)
      real(dp) :: hessian_pressure = DEFAULT_PRESSURE          !! Pressure for thermochemistry (atm)

      ! AIMD settings
      real(dp) :: aimd_dt = DEFAULT_AIMD_DT                          !! Timestep (femtoseconds)
      integer :: aimd_nsteps = DEFAULT_AIMD_NSTEPS                   !! Number of MD steps (0 = no AIMD)
      real(dp) :: aimd_initial_temperature = DEFAULT_AIMD_TEMPERATURE  !! Initial temperature for velocity init (K)
      integer :: aimd_output_frequency = DEFAULT_AIMD_OUTPUT_FREQ    !! Write output every N steps

      ! Fragmentation settings
      character(len=:), allocatable :: frag_method  !! MBE, etc.
      integer :: frag_level = DEFAULT_FRAG_LEVEL
      logical :: allow_overlapping_fragments = .false.
      integer :: max_intersection_level = DEFAULT_MAX_INTERSECTION  !! Maximum k-way intersection depth for GMBE
      character(len=:), allocatable :: embedding
      character(len=:), allocatable :: cutoff_method
      character(len=:), allocatable :: distance_metric
      real(dp), allocatable :: fragment_cutoffs(:)  !! Distance cutoffs indexed by n-mer level (2=dimer, 3=trimer, etc.)
      integer :: global_groups = 0
      integer :: nodes_per_group = 0

      ! Logger settings (kept for compatibility)
      character(len=:), allocatable :: log_level

      ! Output control
      logical :: skip_json_output = .false.  !! Skip JSON output for large calculations
      character(len=:), allocatable :: fragment_breakdown
         !! Where the per-fragment MBE table goes: "csv", "json" or "none"

   contains
      procedure :: destroy => config_destroy
   end type mqc_config_t

contains

   subroutine molecule_destroy(this)
      !! Clean up allocated memory in molecule_t
      class(molecule_t), intent(inout) :: this
      integer :: i

      if (allocated(this%name)) deallocate (this%name)

      call this%geometry%destroy()

      if (allocated(this%fragments)) then
         do i = 1, size(this%fragments)
            call this%fragments(i)%destroy()
         end do
         deallocate (this%fragments)
      end if

      if (allocated(this%bonds)) deallocate (this%bonds)

   end subroutine molecule_destroy

   subroutine config_destroy(this)
      !! Clean up allocated memory in mqc_config_t
      class(mqc_config_t), intent(inout) :: this
      integer :: i

      if (allocated(this%schema_name)) deallocate (this%schema_name)
      if (allocated(this%schema_version)) deallocate (this%schema_version)
      if (allocated(this%units)) deallocate (this%units)
      if (allocated(this%basis)) deallocate (this%basis)
      if (allocated(this%aux_basis)) deallocate (this%aux_basis)
      if (allocated(this%functional)) deallocate (this%functional)
      if (allocated(this%log_level)) deallocate (this%log_level)
      if (allocated(this%fragment_breakdown)) deallocate (this%fragment_breakdown)
      if (allocated(this%frag_method)) deallocate (this%frag_method)
      if (allocated(this%embedding)) deallocate (this%embedding)
      if (allocated(this%cutoff_method)) deallocate (this%cutoff_method)
      if (allocated(this%distance_metric)) deallocate (this%distance_metric)
      if (allocated(this%fragment_cutoffs)) deallocate (this%fragment_cutoffs)

      call this%geometry%destroy()

      if (allocated(this%fragments)) then
         do i = 1, size(this%fragments)
            call this%fragments(i)%destroy()
         end do
         deallocate (this%fragments)
      end if

      if (allocated(this%bonds)) deallocate (this%bonds)

      ! Clean up molecules array (multi-molecule mode)
      if (allocated(this%molecules)) then
         do i = 1, size(this%molecules)
            call this%molecules(i)%destroy()
         end do
         deallocate (this%molecules)
      end if

   end subroutine config_destroy

   subroutine input_fragment_destroy(this)
      !! Clean up allocated memory in input_fragment_t
      class(input_fragment_t), intent(inout) :: this
      if (allocated(this%indices)) deallocate (this%indices)
   end subroutine input_fragment_destroy
end module mqc_config_types
