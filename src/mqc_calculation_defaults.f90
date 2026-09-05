!! Centralized default values for calculation parameters
module mqc_calculation_defaults
   !! Compile-time constants for the values a deck gets when it names none,
   !! kept in one place so the serial and parallel paths cannot diverge.
   use pic_types, only: dp
   implicit none
   private

   ! =========================================================================
   ! Hessian / Finite Difference
   ! =========================================================================
   real(dp), parameter, public :: DEFAULT_DISPLACEMENT = 0.005_dp  !! Bohr
   ! TODO(mqc): the Angstrom equivalent given here read ~0.05, where 0.005 Bohr
   ! is 0.0026. `DEFAULT_FD_DISPLACEMENT` in `mqc_program_limits` is the same
   ! 0.005 Bohr declared a second time, correctly documented, and used by
   ! nothing.
   real(dp), parameter, public :: DEFAULT_TEMPERATURE = 298.15_dp  !! K (room temperature)
   real(dp), parameter, public :: DEFAULT_PRESSURE = 1.0_dp        !! atm (standard pressure)
   real(dp), parameter, public :: DEFAULT_RESPONSE_TOL = 1.0e-9_dp
   !! Largest residual at which an analytic Hessian's coupled-perturbed solve
   !! stops, on the orbital response scaled by the energy denominators. The
   !! Hessian is linear in that response, so its error tracks this number.
   integer, parameter, public :: DEFAULT_RESPONSE_MAX_ITER = 50
   !! Krylov cycles the solve may take before it reports non-convergence.

   ! =========================================================================
   ! SCF
   ! =========================================================================
   integer, parameter, public :: DEFAULT_SCF_MAXITER = 100
   real(dp), parameter, public :: DEFAULT_SCF_CONV = 1.0e-9_dp
   !! Convergence threshold a deck gets when it does not name
   !! `keywords.scf.tolerance`.
   real(dp), parameter, public :: DEFAULT_SCF_DENSITY_CONV = 1.0e-6_dp
   !! Density convergence a deck gets when it does not say.
   integer, parameter, public :: DEFAULT_RESPONSE_BATCH = 12
   !! Densities contracted against one pass over the integrals in the
   !! frequency-dependent response.
   !!
   !! A machine number, not a physical one: a pass evaluates the same quartets
   !! whatever it contracts against them, so wide batches amortise the
   !! integrals, while each thread scatters into an accumulator that grows with
   !! the batch and with the square of the basis. Where the accumulator starts
   !! to win back depends on the cache, which is why
   !! `keywords.efp.response_batch` exists rather than a second guess.
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

   real(dp), parameter, public :: DEFAULT_OPT_GRADIENT_TOLERANCE = 4.5e-4_dp
   !! Convergence on the largest gradient component, Hartree/Bohr. DL-FIND's
   !! own default, restated here so a deck that never mentions it still gets a
   !! number this program can print.

   real(dp), parameter, public :: DEFAULT_OPT_ENERGY_TOLERANCE = -1.0_dp
   !! Negative means "leave the engine's own default alone", which is what the
   !! settings below it use: the step cap and the curvature history depend on
   !! the coordinate system, so DL-FIND chooses them.
   real(dp), parameter, public :: DEFAULT_OPT_MAX_STEP = -1.0_dp
   integer, parameter, public :: DEFAULT_OPT_LBFGS_MEMORY = -1
   integer, parameter, public :: DEFAULT_OPT_PRINT_LEVEL = -1

   ! =========================================================================
   ! EFP / MAKEFP
   ! =========================================================================
   real(dp), parameter, public :: DEFAULT_VDW_SCALE = 0.7_dp
   !! Innermost layer of the EFP screening grid, as a fraction of a van der
   !! Waals radius. GAMESS's `VDWSCL`.
   !!
   !! The grid takes the scale the potential was built with, so the `VDWSCL`
   !! the SCREEN header reports is the one the fit was made against.

   real(dp), parameter, public :: DEFAULT_DYNAMIC_TOL = 1.0e-7_dp
   !! What the frequency-dependent response solve converges its residual to,
   !! relative to the right-hand side, when a deck does not name a tolerance.
   !!
   !! **The largest cheap lever on that solve.** Every iteration is four passes
   !! over the integrals, so the iteration count is the cost and this sets the
   !! count. The EFMO literature runs the equivalent CPHF and TDHF solves at
   !! `5e-5` -- Sattasathuchana et al., JCTC 20, 2445 (2024) -- for an error of
   !! about `3e-6` Hartree in the total interaction energy.
   !!
   !! Left at `1e-7` because what a potential built here does to an interaction
   !! energy has not been measured. `keywords.efp.dynamic_tolerance` is the one
   !! EFP key that changes the numbers a potential reports rather than only how
   !! they are arrived at.

   integer, parameter, public :: DEFAULT_DYNAMIC_MAXITER = 200
   !! Iterations that solve gets before it gives up, per system.

   integer, parameter, public :: EFP_RESPONSE_AUTO = 0
   !! How the frequency-dependent response is obtained: `keywords.efp.response`.
   !!
   !! `auto` builds `(A+B)` and `(A-B)` and factorises them, unless three
   !! `n_ov^2` matrices would pass `DENSE_OPERATOR_LIMIT`, and iterates
   !! matrix-free when they would. `dense` builds the operator however large it
   !! is and `matrix_free` never builds it however small; the crossover is a
   !! statement about memory rather than wall clock, and the two routes solve
   !! the same equations.
   integer, parameter, public :: EFP_RESPONSE_DENSE = 1
   integer, parameter, public :: EFP_RESPONSE_MATRIX_FREE = 2

   ! =========================================================================
   ! XTB
   ! =========================================================================
   real(dp), parameter, public :: DEFAULT_XTB_ACCURACY = 0.01_dp
   integer, parameter, public :: DEFAULT_CPCM_NANG = 110
   real(dp), parameter, public :: DEFAULT_CPCM_RSCALE = 1.0_dp

   integer, parameter, public :: DEFAULT_PCM_NANG = 302
   !! Angular points per atom on a continuum cavity, for the ab initio PCM.
   !!
   !! 302, matching what PySCF and most codes use for a cavity. A continuum
   !! surface is a two-dimensional integral of a function with a cusp where
   !! spheres intersect, so it is not the easy quadrature it looks like: at 110
   !! points per atom the dielectric energy comes out 21% short.

   real(dp), parameter, public :: DEFAULT_PCM_RSCALE = 1.2_dp
   !! Cavity radius scaling. See `mqc_pcm_radii` for why 1.2 and not 1.0.

   real(dp), parameter, public :: DEFAULT_PCM_ZETA = 2.0_dp
   !! Prefactor of the Gaussian switching exponent on a smooth cavity.
   !!
   !! Used as zeta_i = zeta * sqrt(n_angular) / R_i, since the exponent has to
   !! scale with the spacing between surface points, which on a sphere of radius
   !! R carrying n points goes as R/sqrt(n).
   !!
   !! **Calibrated, not derived.** Lange and Herbert [JCP 133, 244111 (2010)]
   !! fit 4.9 for Lebedev cavities, but cuEST takes the exponents themselves
   !! rather than a prefactor and its convention is undocumented, so this number
   !! absorbs whatever that definition is. Treat it as an empirical constant of
   !! this interface rather than a physical one.

   ! =========================================================================
   ! Fragmentation
   ! =========================================================================
   integer, parameter, public :: DEFAULT_FRAG_LEVEL = 1
   integer, parameter, public :: DEFAULT_MAX_INTERSECTION = 999

   ! Fragment type identifiers for MPI communication
   integer, parameter, public :: FRAGMENT_TYPE_MONOMERS = 0  !! Fragment specified by monomer indices (MBE)
   integer, parameter, public :: FRAGMENT_TYPE_ATOMS = 1     !! Fragment specified by atom list (GMBE/PIE)

   integer, parameter, public :: DISP_WHOLE_FRAGMENT = huge(0)
   !! Displacement code carried by every work task.
   !!
   !! DISP_WHOLE_FRAGMENT means "run the requested calc_type on the fragment as
   !! it stands", which is every energy/gradient run and the whole GMBE path.
   !! Any other value marks a task from a flattened Hessian queue: 0 is the
   !! undisplaced reference point and +/-k displaces coordinate k
   !! forward/backward (see `build_hessian_task_table`).

end module mqc_calculation_defaults
