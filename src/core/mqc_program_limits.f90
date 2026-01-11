!! Program limits and default parameter, publics
module mqc_program_limits
   !! Contains compile-time limits and default values for the metalquicha program.
   !! These are tunable parameter that control memory allocation and algorithm behavior.
   use pic_types, only: dp
   implicit none
   private

   !---------------------------------------------------------------------------
   ! Many-Body Expansion Limits
   !---------------------------------------------------------------------------

   !> Maximum MBE truncation order (1-body, 2-body, ..., N-body)
   !> Higher orders require factorial growth in fragment combinations
   integer, parameter, public :: MAX_MBE_LEVEL = 10

   !---------------------------------------------------------------------------
   ! Numerical Differentiation Defaults
   !---------------------------------------------------------------------------

   !> Default step size for finite difference calculations (Bohr)
   !> ~0.005 Bohr = ~0.0026 Angstrom, suitable for Hessian/gradient FD
   real(dp), parameter, public :: DEFAULT_FD_DISPLACEMENT = 0.005_dp

   !---------------------------------------------------------------------------
   ! I/O Limits
   !---------------------------------------------------------------------------

   !> Maximum length for input file lines
   integer, parameter, public :: MAX_LINE_LENGTH = 1024

   !> Maximum length for element symbols (e.g., "He", "Uue")
   integer, parameter, public :: MAX_ELEMENT_SYMBOL_LEN = 4

   !> JSON output format for real numbers (scientific notation)
   !> Valid values: 'G', 'E', 'EN', 'ES' (json-fortran uses machine precision)
   character(len=*), parameter, public :: JSON_REAL_FORMAT = 'ES'

   !---------------------------------------------------------------------------
   ! Geometry/Structure Limits
   !---------------------------------------------------------------------------

   !> Minimum allowed distance between atoms (Bohr)
   !> Atoms closer than this are considered overlapping (error condition)
   real(dp), parameter, public :: MIN_ATOM_DISTANCE = 0.01_dp

end module mqc_program_limits
