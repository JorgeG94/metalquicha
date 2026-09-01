!! The types that say how a self-consistent field is set up and driven
module mqc_scf_types
   !! **Here rather than in `mqc_config_types`, where these started.**
   !!
   !! `guess_step_t` was put beside the deck types years ago, and when a shared
   !! SCF settings type was needed it followed that precedent -- so "how an SCF
   !! converges" ended up declared in the I/O layer, which is not what the I/O
   !! layer is for. A deck reads these; it does not own them.
   !!
   !! The direction works because `src/scf` depends on nothing but `mqc_error`,
   !! so `mqc_config_types` can name these without a cycle. Every SCF in the
   !! program -- the supermolecular one, the deltaSCF ions behind a Fukui
   !! analysis, each rung of a basis-projection ladder, each fragment of an FMO
   !! run -- is configured by the same type, and that is the point: a fifth
   !! parallel copy of the same fields is how settings went missing.
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
      !! Each rung converges an SCF in its own basis and hands the result to the
      !! next as a starting density. The convergence settings are per rung
      !! because the rungs are not doing the same job: an early one only has to
      !! land in the right basin and can stop loosely, where the last one before
      !! the target may as well be tight since it is cheap relative to what
      !! follows.
      character(len=:), allocatable :: basis
      integer :: maxiter = 50
      real(dp) :: tolerance = 1.0e-5_dp
   end type guess_step_t

   type :: scf_numerics_t
      !! How a self-consistent field is driven to convergence, defined once
      !!
      !! Split out of `scf_options_t` so that a *second* SCF in the same run
      !! can be configured without inheriting things that are not about
      !! convergence. The deltaSCF behind a Fukui analysis is the case that
      !! forced it: those ions need their own iteration limit and level shift,
      !! but they cannot have their own basis (the density difference is taken
      !! over the neutral's basis functions by construction) and they must not
      !! have their own `properties` -- a Fukui analysis nested inside a Fukui
      !! analysis is not a thing anyone is asking for.
      !!
      !! That last point is not only taste. `scf_options_t` carries a
      !! `properties_config_t`, and `properties_config_t` is where the ion
      !! settings live, so a `deltascf_options_t` extending the *full* options
      !! type would contain itself. Fortran takes that if the component is
      !! allocatable, and it would mean nothing.
      !!
      !! The boundary is "how the iteration is driven". Anything selecting
      !! *which* wave function is being converged -- `unrestricted`,
      !! `density_fitting`, the basis, the frozen core -- stays in
      !! `scf_options_t`, because a second SCF over the same molecule has no
      !! business disagreeing about those.
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
         !! Zero is off -- which is why this could never carry a sentinel for
         !! "unset", and why the ions used to spell theirs "any negative".
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
      character(len=32) :: convergence_metric = "standard"
         !! Which measure decides this SCF has stopped. See
         !! `mqc_scf_convergence`; the default is the energy-and-commutator
         !! pair every loop tested before that module existed.
         !! 'diis' (the default), 'adiis' or 'ediis'. The energy-based pair
         !! runs only while the error is large and hands over to DIIS, so
         !! naming one asks for a different opening, not a different endgame.
      logical :: allow_crap_scf = .false.
         !! Keep a non-converged SCF instead of failing
      character(len=32) :: guess = "auto"
         !! Initial guess: 'core', 'gwh', 'sac', 'sad', 'basis_set_projection'
         !! or 'auto', where the backend picks.
         !!
         !! **Not the same key as `mqc_config_t%fukui_guess`**, which takes
         !! 'neutral' or 'independent' and says whether the ion starts from the
         !! neutral's converged density at all. This one says how a density is
         !! built when it starts from nothing. They are one level apart in the
         !! deck and spelled the same, which is a trap: see `fukui_keys` in
         !! `mqc_json_schema`.
      type(guess_step_t), allocatable :: guess_steps(:)
         !! The basis ladder for 'basis_set_projection', one entry per
         !! preliminary SCF in order.
   end type scf_numerics_t

   type, extends(scf_numerics_t) :: deltascf_options_t
      !! Convergence settings for a second SCF run beside the first
      !!
      !! The ions behind a Fukui analysis are their own SCF problem and are
      !! allowed to be configured as one. They are harder than the neutral -- a
      !! charged species in a basis chosen to describe a neutral one -- so
      !! inheriting whatever converged the neutral is a default rather than a
      !! rule.
      !!
      !! **Inheritance happens at read time, not here.** The reader seeds every
      !! field above from the resolved `keywords.scf` before it looks at the
      !! deck, so a key the deck does not name keeps the neutral's value and a
      !! key it does name wins. That is why nothing here carries a sentinel:
      !! there is no "unset" state to encode, because "unset" is spelled by the
      !! field already holding the inherited value. The three fields this
      !! replaced each needed one, and `level_shift` could not have a clean one
      !! at all.
      logical :: inherit_scf = .true.
         !! Seed from `keywords.scf` (the default) or from the defaults above.
         !!
         !! False is for the case where the neutral was converged tightly, or
         !! with an accelerator chosen for it, and carrying that onto a charged
         !! species is the wrong starting point rather than a cautious one. It
         !! does not stop the deck naming keys explicitly -- those still win
         !! either way; it only changes what the unnamed ones fall back to.
   end type deltascf_options_t

contains

   subroutine print_scf_config(scf, label)
      !! Echo what an SCF was actually told to do
      !!
      !! Every silent-drop bug this code has had looks the same from outside: a
      !! key is set, the schema accepts it, and the run behaves as though it
      !! were never written. The Fukui ions lost six settings that way and the
      !! MakeFP SCF lost six more, and in both cases the only way to find out
      !! was to read the call site. A deck that can be echoed back cannot hide
      !! that -- if a setting is missing from this block, it did not arrive.
      !!
      !! Printed from the resolved configuration rather than from the deck, so
      !! it shows what the SCF will use and not what was asked for. Those differ
      !! wherever a default, an inheritance or a backend override sits between
      !! the two, which is exactly where the interesting bugs are.
      use pic_logger, only: logger => global_logger
      class(scf_numerics_t), intent(in) :: scf
      character(len=*), intent(in) :: label
         !! Which SCF this is -- there is more than one in a Fukui or MakeFP
         !! run, and two identical blocks with no names would be worse than
         !! none.

      character(len=160) :: line

      call logger%info("  "//label//" configuration")
      write (line, "(a,i0,a,es9.2,a,es9.2)") &
         "    maxiter ", scf%max_iter, "   energy_tol ", scf%energy_tol, &
         "   density_tol ", scf%density_tol
      call logger%info(trim(line))
      ! **The RESOLVED threshold, not the field it came from.** Zero means the
      ! SCF derives it as `sqrt(energy_tol)`, so that is what has to be
      ! printed: reporting the 1e-8 that a 1e-4 gate was derived from tells the
      ! reader the opposite of what the run did, and this block exists to be
      ! read against a run that stopped somewhere surprising.
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
