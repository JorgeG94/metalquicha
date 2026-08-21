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
   !> Density convergence a deck gets when it does not say. Matches what
   !> `method_config_t` already defaulted to, so naming the key changes nothing
   !> for a deck that leaves it alone.
   real(dp), parameter, public :: DEFAULT_SCF_DENSITY_CONV = 1.0e-6_dp
   logical, parameter, public :: DEFAULT_USE_DIIS = .true.

   ! =========================================================================
   ! AIMD
   ! =========================================================================
   real(dp), parameter, public :: DEFAULT_AIMD_DT = 1.0_dp         !! fs
   integer, parameter, public :: DEFAULT_AIMD_NSTEPS = 0
   real(dp), parameter, public :: DEFAULT_AIMD_TEMPERATURE = 300.0_dp  !! K
   integer, parameter, public :: DEFAULT_AIMD_OUTPUT_FREQ = 1

   ! =========================================================================
   ! Geometry optimization
   ! =========================================================================
   integer, parameter, public :: DEFAULT_OPT_MAX_STEPS = 100

   !> Convergence on the largest gradient component, Hartree/Bohr.
   !> DL-FIND's own default, restated here so a deck that never mentions it
   !> still gets a number this program can print and a user can look up.
   real(dp), parameter, public :: DEFAULT_OPT_GRADIENT_TOLERANCE = 4.5e-4_dp

   !> Negative means "leave the engine's own default alone". Used for the
   !> settings where DL-FIND's choice is better than a second one invented
   !> here: the step cap and the curvature history depend on the coordinate
   !> system, which this program does not want to reason about.
   real(dp), parameter, public :: DEFAULT_OPT_ENERGY_TOLERANCE = -1.0_dp
   real(dp), parameter, public :: DEFAULT_OPT_MAX_STEP = -1.0_dp
   integer, parameter, public :: DEFAULT_OPT_LBFGS_MEMORY = -1
   integer, parameter, public :: DEFAULT_OPT_PRINT_LEVEL = -1

   ! =========================================================================
   ! EFP / MAKEFP
   ! =========================================================================
   !> Innermost layer of the EFP screening grid, as a fraction of a van der
   !> Waals radius. GAMESS's `VDWSCL`.
   !>
   !> One name for what used to be `VDW_SCALE` in the screening grid and
   !> `vdwscl` on the potential writer, and now one number as well: the grid
   !> takes the scale the potential was built with, so the `VDWSCL` the SCREEN
   !> header reports is the one the fit was made against. They were only
   !> *equal* before, which was fine while nobody could change either -- and
   !> not fine the moment `keywords.efp.vdw_scale` existed, because it would
   !> have moved the header and left the grid where it was.
   real(dp), parameter, public :: DEFAULT_VDW_SCALE = 0.7_dp

   !> What the frequency-dependent response solve converges its residual to,
   !> relative to the right-hand side, when a deck does not name a tolerance.
   !>
   !> **The largest cheap lever on that solve.** Every iteration is four passes
   !> over the integrals, so the iteration count is the cost and this sets the
   !> count. The EFMO literature runs the equivalent CPHF and TDHF solves at
   !> `5e-5` and reports the resulting error in the total interaction energy at
   !> about `3e-6` Hartree -- Sattasathuchana et al., JCTC 20, 2445 (2024),
   !> Tables 2 and 3, where loosening from `1e-7` to `5e-5` cut their adenine
   !> wall time from 39.1 to 17.2 minutes.
   !>
   !> Left at `1e-7` all the same. That measurement is of their code and their
   !> accumulation, and what a potential built here does to an interaction
   !> energy still has not been measured. `keywords.efp.dynamic_tolerance` is
   !> how that measurement gets made without a rebuild, and it is the one EFP
   !> key that changes the numbers a potential reports rather than only how
   !> they are arrived at.
   !>
   !> Here rather than beside the solver so that the default a deck inherits
   !> and the default the solver falls back on are the same constant.
   real(dp), parameter, public :: DEFAULT_DYNAMIC_TOL = 1.0e-7_dp

   !> Iterations that solve gets before it gives up, per system.
   integer, parameter, public :: DEFAULT_DYNAMIC_MAXITER = 200

   !> How the frequency-dependent response is obtained: `keywords.efp.response`.
   !>
   !> `auto` is the rule the code followed before there was a keyword -- build
   !> `(A+B)` and `(A-B)` and factorise them, unless three `n_ov^2` matrices
   !> would pass `DENSE_OPERATOR_LIMIT`, and iterate matrix-free when they
   !> would. The two forcing values exist because that crossover is a statement
   !> about memory and not about wall clock, and the two routes could not be
   !> compared from a deck at all: `dense` builds the operator however large it
   !> is, `matrix_free` never builds it however small. They solve the same
   !> equations, so a small molecule run both ways prices the iterative route
   !> against an exact reference.
   integer, parameter, public :: EFP_RESPONSE_AUTO = 0
   integer, parameter, public :: EFP_RESPONSE_DENSE = 1
   integer, parameter, public :: EFP_RESPONSE_MATRIX_FREE = 2

   ! =========================================================================
   ! XTB
   ! =========================================================================
   real(dp), parameter, public :: DEFAULT_XTB_ACCURACY = 0.01_dp
   integer, parameter, public :: DEFAULT_CPCM_NANG = 110
   real(dp), parameter, public :: DEFAULT_CPCM_RSCALE = 1.0_dp

   !> Angular points per atom on a continuum cavity, for the ab initio PCM.
   !>
   !> 302, matching what PySCF and most codes use for a cavity. This was 110 on
   !> the assumption that a cavity needs far fewer points than an
   !> exchange-correlation grid, and that was measured to be wrong: on water in
   !> water, 110 points per atom gave a dielectric energy 21% short of the
   !> converged value (-8.07 mHartree against -10.21), while 302 agrees with
   !> PySCF's C-PCM to 1.5% and reproduces its solvated dipole to 0.0013 D. A
   !> continuum surface is a two-dimensional integral of a function with a cusp
   !> where spheres intersect, so it is not the easy quadrature it looks like.
   integer, parameter, public :: DEFAULT_PCM_NANG = 302

   !> Cavity radius scaling. See `mqc_pcm_radii` for why 1.2 and not 1.0.
   real(dp), parameter, public :: DEFAULT_PCM_RSCALE = 1.2_dp

   !> Prefactor of the Gaussian switching exponent on a smooth cavity.
   !>
   !> Used as zeta_i = zeta * sqrt(n_angular) / R_i, since the exponent has to
   !> scale with the spacing between surface points, which on a sphere of radius
   !> R carrying n points goes as R/sqrt(n).
   !>
   !> **Calibrated, not derived, and that is worth stating.** Lange and Herbert
   !> [JCP 133, 244111 (2010)] fit 4.9 for Lebedev cavities, but cuEST takes the
   !> exponents themselves rather than a prefactor and its convention is not
   !> documented in the bindings -- so 4.9 here is not necessarily 4.9 there. Two
   !> is the value that reproduces PySCF's C-PCM on water in water; at 4.9 the
   !> dielectric energy came out half. What that means is that this number
   !> absorbs whatever cuEST's definition actually is, so it should be treated as
   !> an empirical constant of this interface rather than a physical one, and it
   !> is a keyword so it can be revisited.
   real(dp), parameter, public :: DEFAULT_PCM_ZETA = 2.0_dp

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
