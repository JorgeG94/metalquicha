!! The types that say how a self-consistent field is set up and driven
module mqc_scf_types
   !! Every SCF in the program -- the supermolecular one, the deltaSCF ions
   !! behind a Fukui analysis, each rung of a basis-projection ladder, each
   !! fragment of an FMO run -- is configured by the types here.
   !!
   !! They live in `src/scf` rather than beside the deck types they are read
   !! from: `src/scf` depends on nothing but `mqc_error`, so `mqc_config_types`
   !! can name them without a cycle.
   use pic_types, only: dp
   implicit none
   private

   public :: guess_step_t       !! One rung of a basis-set-projection ladder
   public :: scf_numerics_t     !! How an SCF is driven, shared by every SCF in a run
   public :: deltascf_options_t  !! Convergence settings for a second SCF beside the first
   public :: print_scf_config    !! Echo a resolved SCF configuration

   type :: guess_step_t
      !! One rung of a basis-set-projection ladder
      !!
      !! Each rung converges an SCF in its own basis and hands the result to
      !! the next as a starting density. Convergence settings are per rung: an
      !! early one only has to land in the right basin.
      character(len=:), allocatable :: basis
      integer :: maxiter = 50
      real(dp) :: tolerance = 1.0e-5_dp
   end type guess_step_t

   type :: scf_numerics_t
      !! How a self-consistent field is driven to convergence, defined once
      !!
      !! The boundary is "how the iteration is driven". Anything selecting
      !! *which* wave function is being converged -- `unrestricted`,
      !! `density_fitting`, the basis, the frozen core -- stays in
      !! `scf_options_t`, because a second SCF over the same molecule has no
      !! business disagreeing about those.
      ! Split out of `scf_options_t`, and not extended from it, because
      ! `scf_options_t` carries a `properties_config_t` and that is where the
      ! ion settings live: a `deltascf_options_t` extending the full options
      ! type would contain itself.
      integer :: max_iter = 100
         !! Maximum SCF iterations
      real(dp) :: energy_tol = 1.0e-8_dp
         !! Energy convergence threshold
      real(dp) :: density_tol = 1.0e-6_dp
         !! Density matrix convergence threshold
      real(dp) :: grad_tol = 0.0_dp
         !! Commutator threshold; zero derives `sqrt(energy_tol)`.
      real(dp) :: linear_dependence = 0.0_dp
         !! Zero means the orthogonaliser's own cutoff.
      real(dp) :: level_shift = 0.0_dp
         !! Hartree added to the virtual block before each diagonalisation.
         !! Zero is off.
      logical :: use_diis = .true.
         !! Use DIIS acceleration
      integer :: diis_size = 8
         !! Number of Fock matrices for DIIS
      logical :: incremental_fock = .true.
         !! Build each Fock matrix from the density *change* since the last
         !! full build. Exact to the convergence threshold and several times
         !! cheaper late in an SCF; false forces a full build every iteration,
         !! which is the first thing to rule out when a run stalls.
      character(len=32) :: accelerator = "diis"
         !! 'diis' (the default), 'adiis' or 'ediis'. The energy-based pair
         !! runs only while the error is large and hands over to DIIS, so
         !! naming one asks for a different opening, not a different endgame.
      character(len=32) :: convergence_metric = "standard"
         !! Which measure decides this SCF has stopped. See
         !! `mqc_scf_convergence`; the default is the energy-and-commutator
         !! pair.
      logical :: allow_crap_scf = .false.
         !! Keep a non-converged SCF instead of failing
      character(len=32) :: guess = "auto"
         !! Initial guess: 'core', 'gwh', 'sac', 'sad', 'basis_set_projection'
         !! or 'auto', where the backend picks.
         !!
         !! **Not the same key as `mqc_config_t%fukui_guess`**, which takes
         !! 'neutral' or 'independent' and says whether the ion starts from the
         !! neutral's converged density at all. This one says how a density is
         !! built when it starts from nothing. The two are one level apart in
         !! the deck and spelled the same.
      type(guess_step_t), allocatable :: guess_steps(:)
         !! The basis ladder for 'basis_set_projection', one entry per
         !! preliminary SCF in order.
   end type scf_numerics_t

   type, extends(scf_numerics_t) :: deltascf_options_t
      !! Convergence settings for a second SCF run beside the first
      !!
      !! The ions behind a Fukui analysis are their own SCF problem: a charged
      !! species in a basis chosen for a neutral one, so whatever converged the
      !! neutral is a default rather than a rule.
      !!
      !! **Inheritance happens at read time, not here.** The reader seeds every
      !! inherited field from the resolved `keywords.scf` before it looks at
      !! the deck, so a key the deck does not name keeps the neutral's value
      !! and a key it does name wins. Nothing here carries an "unset" sentinel
      !! because unset is spelled by the field holding the inherited value.
      logical :: inherit_scf = .true.
         !! Seed from `keywords.scf` (the default) or from the defaults above.
         !! Only changes what the keys the deck does not name fall back to;
         !! keys it does name win either way.
   end type deltascf_options_t

contains

   subroutine print_scf_config(scf, label)
      !! Echo what an SCF was actually told to do
      !!
      !! Printed from the resolved configuration rather than from the deck, so
      !! it shows what the SCF will use and not what was asked for. Those
      !! differ wherever a default, an inheritance or a backend override sits
      !! between the two. A setting missing from this block did not arrive.
      use pic_logger, only: logger => global_logger
      class(scf_numerics_t), intent(in) :: scf
      character(len=*), intent(in) :: label
         !! Which SCF this is -- a Fukui or MakeFP run has more than one

      character(len=160) :: line

      call logger%info("  "//label//" configuration")
      write (line, "(a,i0,a,es9.2,a,es9.2)") &
         "    maxiter ", scf%max_iter, "   energy_tol ", scf%energy_tol, &
         "   density_tol ", scf%density_tol
      call logger%info(trim(line))
      ! The resolved threshold, not the field it came from: zero means the SCF
      ! derives it as `sqrt(energy_tol)`, so that is what is printed.
      if (scf%grad_tol > 0.0_dp) then
         write (line, "(a,es9.2,a)") "    commutator_tol ", scf%grad_tol, "  (stated)"
      else
         write (line, "(a,es9.2,a)") "    commutator_tol ", sqrt(scf%energy_tol), &
            "  (derived, sqrt of energy_tol)"
      end if
      call logger%info(trim(line))
      write (line, "(a,a,a,i0,a,l1)") &
         "    accelerator ", trim(scf%accelerator), "   diis_size ", scf%diis_size, &
         "   diis ", scf%use_diis
      call logger%info(trim(line))
      ! Same again: zero is a sentinel the orthogonaliser turns into its own
      ! cutoff, so printing 0.00E+00 would report a threshold nobody uses.
      if (scf%linear_dependence > 0.0_dp) then
         write (line, "(a,es9.2,a,es9.2)") &
            "    level_shift ", scf%level_shift, "   linear_dependence ", scf%linear_dependence
      else
         write (line, "(a,es9.2,a)") &
            "    level_shift ", scf%level_shift, "   linear_dependence  backend default"
      end if
      call logger%info(trim(line))
      write (line, "(a,a,a,l1,a,l1)") &
         "    guess ", trim(scf%guess), "   incremental_fock ", scf%incremental_fock, &
         "   allow_crap_scf ", scf%allow_crap_scf
      call logger%info(trim(line))
   end subroutine print_scf_config

end module mqc_scf_types
