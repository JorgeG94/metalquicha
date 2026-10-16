!! What a deck asks for when it wants orbitals written as Gaussian cube files
module mqc_cube_config
   !! `cube_config_t` and the codes for which orbitals it selects. Backend
   !! neutral, so the configuration types can carry it without reaching into
   !! any backend.
   use pic_types, only: dp
   implicit none
   private

   public :: cube_config_t
   public :: CUBE_NONE, CUBE_OCCUPIED, CUBE_ALL, CUBE_LISTED
   public :: DEFAULT_CUBE_SPACING, DEFAULT_CUBE_MARGIN, DEFAULT_CUBE_PREFIX

   integer, parameter :: CUBE_NONE = 0
      !! No cube files. The default, and what an absent `cube` key means.
   integer, parameter :: CUBE_OCCUPIED = 1
      !! Every orbital the analysis does not call empty.
   integer, parameter :: CUBE_ALL = 2
      !! Every orbital, empty ones included.
   integer, parameter :: CUBE_LISTED = 3
      !! The orbitals named in `cube_config_t%indices`.

   real(dp), parameter :: DEFAULT_CUBE_SPACING = 0.2_dp
      !! Bohr between grid points.
   real(dp), parameter :: DEFAULT_CUBE_MARGIN = 4.0_dp
      !! Bohr of grid beyond the outermost nucleus on each side.
   character(len=*), parameter :: DEFAULT_CUBE_PREFIX = "quao"
      !! The tag in every file name, after the deck's output stem.

   type :: cube_config_t
      !! `properties.bonding_analysis.cube`
      integer :: selection = CUBE_NONE
         !! One of the `CUBE_*` codes.
      integer, allocatable :: indices(:)
         !! Orbital numbers as the printed report and the JSON rows number
         !! them: 1-based, counting the core orbitals the analysis leaves out,
         !! so the first quasi-atomic orbital of a first-row molecule is 2.
         !! Allocated when `selection` is `CUBE_LISTED`.
      logical :: single_file = .false.
         !! Write every selected orbital into one multi-orbital cube file
         !! rather than one file each.
      real(dp) :: spacing = DEFAULT_CUBE_SPACING
         !! Bohr between grid points along each axis. Positive.
      real(dp) :: margin = DEFAULT_CUBE_MARGIN
         !! Bohr beyond the outermost nucleus on each side. Not negative.
      character(len=64) :: prefix = DEFAULT_CUBE_PREFIX
         !! Tag put into every file name. Letters, digits, `_` and `-` only.
   end type cube_config_t

end module mqc_cube_config
