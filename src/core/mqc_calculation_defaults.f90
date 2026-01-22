!! Centralized default values for calculation parameters
module mqc_calculation_defaults
   !! Provides compile-time constants for default calculation parameters.
   !! These defaults are used throughout the codebase when users don't specify values.
   !! This single source of truth prevents divergence between serial and parallel paths.
   use pic_types, only: dp
   implicit none
   private

   ! =========================================================================
   ! Hessian / Finite Difference
   ! =========================================================================
   real(dp), parameter, public :: DEFAULT_DISPLACEMENT = 0.005_dp  !! Bohr (~0.05 Angstrom)
   real(dp), parameter, public :: DEFAULT_TEMPERATURE = 298.15_dp  !! K (room temperature)
   real(dp), parameter, public :: DEFAULT_PRESSURE = 1.0_dp        !! atm (standard pressure)

   ! =========================================================================
   ! SCF
   ! =========================================================================
   integer, parameter, public :: DEFAULT_SCF_MAXITER = 100
   real(dp), parameter, public :: DEFAULT_SCF_CONV = 1.0e-6_dp
   logical, parameter, public :: DEFAULT_USE_DIIS = .true.

   ! =========================================================================
   ! AIMD
   ! =========================================================================
   real(dp), parameter, public :: DEFAULT_AIMD_DT = 1.0_dp         !! fs
   integer, parameter, public :: DEFAULT_AIMD_NSTEPS = 0
   real(dp), parameter, public :: DEFAULT_AIMD_TEMPERATURE = 300.0_dp  !! K
   integer, parameter, public :: DEFAULT_AIMD_OUTPUT_FREQ = 1

   ! =========================================================================
   ! XTB
   ! =========================================================================
   real(dp), parameter, public :: DEFAULT_XTB_ACCURACY = 0.01_dp
   integer, parameter, public :: DEFAULT_CPCM_NANG = 110
   real(dp), parameter, public :: DEFAULT_CPCM_RSCALE = 1.0_dp

   ! =========================================================================
   ! Fragmentation
   ! =========================================================================
   integer, parameter, public :: DEFAULT_FRAG_LEVEL = 1
   integer, parameter, public :: DEFAULT_MAX_INTERSECTION = 999

   ! Fragment type identifiers for MPI communication
   integer, parameter, public :: FRAGMENT_TYPE_MONOMERS = 0  !! Fragment specified by monomer indices (MBE)
   integer, parameter, public :: FRAGMENT_TYPE_ATOMS = 1     !! Fragment specified by atom list (GMBE/PIE)

end module mqc_calculation_defaults
