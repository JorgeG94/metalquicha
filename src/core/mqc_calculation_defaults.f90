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

   !> Angular points per atom on a continuum cavity, for the ab initio PCM.
   !>
   !> 110 is a Lebedev order and the usual starting point for a cavity, which
   !> needs far fewer points than an exchange-correlation grid does.
   integer, parameter, public :: DEFAULT_PCM_NANG = 110

   !> Cavity radius scaling. See `mqc_pcm_radii` for why 1.2 and not 1.0.
   real(dp), parameter, public :: DEFAULT_PCM_RSCALE = 1.2_dp

   !> Prefactor of the Gaussian switching exponent on a smooth cavity.
   !>
   !> **Unverified against cuEST's own convention.** A smooth continuum surface
   !> replaces each surface point by a normalised Gaussian, and 4.9 is the value
   !> Lange and Herbert [JCP 133, 244111 (2010)] fit for Lebedev cavities under
   !> the definition zeta_i = zeta * sqrt(n_angular) / R_i. cuEST takes the
   !> exponents themselves rather than this prefactor, so if its definition
   !> differs the cavity is smoothed by the wrong amount -- which changes the
   !> solvation energy without failing. Exposed as a keyword for exactly that
   !> reason: it can be checked on hardware without a rebuild.
   real(dp), parameter, public :: DEFAULT_PCM_ZETA = 4.9_dp

   ! =========================================================================
   ! Fragmentation
   ! =========================================================================
   integer, parameter, public :: DEFAULT_FRAG_LEVEL = 1
   integer, parameter, public :: DEFAULT_MAX_INTERSECTION = 999

   ! Fragment type identifiers for MPI communication
   integer, parameter, public :: FRAGMENT_TYPE_MONOMERS = 0  !! Fragment specified by monomer indices (MBE)
   integer, parameter, public :: FRAGMENT_TYPE_ATOMS = 1     !! Fragment specified by atom list (GMBE/PIE)

   !! Displacement code carried by every work task.
   !!
   !! DISP_WHOLE_FRAGMENT means "run the requested calc_type on the fragment as it
   !! stands" -- the original one-fragment-per-task behaviour, still used by every
   !! energy/gradient run and by the whole GMBE path. Any other value marks a task
   !! from a flattened Hessian queue: 0 is the undisplaced reference point and +/-k
   !! displaces coordinate k forward/backward (see build_hessian_task_table).
   integer, parameter, public :: DISP_WHOLE_FRAGMENT = huge(0)

end module mqc_calculation_defaults
