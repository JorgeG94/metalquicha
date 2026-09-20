!! Linear-response excitation energies for a closed-shell reference
module mqc_czt_tddft
   !! The Tamm-Dancoff and random-phase eigenproblems, singlet and triplet,
   !! over the response operator this backend already applies.
   !!
   !! Linear response asks for the lowest eigenvalues of
   !!
   !!     A_{ia,jb} = d_{ij} d_{ab} (e_a - e_i) + 2(ai|bj) - c_x (ab|ij) + 2 f_xc
   !!     B_{ia,jb} =                             2(ai|jb) - c_x (aj|ib) + 2 f_xc
   !!
   !! -- of `A` alone in the Tamm-Dancoff approximation, of the paired
   !! `[A B; B A]` problem in the full one -- and nothing here computes any of
   !! those integrals: `response_product` does, `mqc_davidson` and
   !! `mqc_czt_rpa_solver` find the eigenvalues. What this module is, is the
   !! arithmetic between the two.
   !!
   !! ## Getting `A` out of `(A+B)` and `(A-B)`
   !!
   !! `response_product` returns `(A+B)u` or `(A-B)u`, never `Au`:
   !!
   !!     (A+B) u = dEps u + 4(ai|bj) u - c_x[(ab|ij)+(aj|ib)] u + 4 f_xc u
   !!     (A-B) u = dEps u             - c_x[(ab|ij)-(aj|ib)] u
   !!
   !! whose half sum is `A` exactly, with no term surviving that should not.
   !! So `A u` is **two** calls, averaged -- and the full RPA is the same two
   !! calls kept apart, which is why the two methods here share one operator
   !! core and differ only in what they do with its two halves.
   !!
   !! The alternative is one call on the *unsymmetrised* transition density
   !! `C_vir u C_occ^T` -- PySCF's `_gen_tda_operation` -- which gives `A`
   !! directly but needs `build_fock_direct_nosym`, whose unfolded permutations
   !! cost about 2.7 times a symmetric build. Two symmetric passes cost 2.0 of
   !! the same unit, so the half sum is both cheaper and free of new integral
   !! code. The margin is wider than that in the two cases that matter most:
   !!
   !! * A **pure functional** has no exchange at all, so `(A-B)` reduces to the
   !!   orbital-energy diagonal and is not built -- one pass, against 2.7.
   !! * A **range-separated hybrid** makes two exchange passes per product
   !!   either way, so the nosym route's penalty applies to both of them.
   !!
   !! `(A-B)` never touches the quadrature, so the grid is walked once per
   !! product whichever route is taken, and with the kernel cache filled it is
   !! not re-evaluated at all.
   !!
   !! ## Triplets
   !!
   !! Two changes, both inside `response_product`, and neither a rescaling of
   !! the other operator:
   !!
   !!     A_T = dEps - c_x (ab|ij) + 2 f_xc^T,   f_xc^T = (f_aa - f_ab)/2
   !!
   !! -- **no Coulomb term at all**, and the spin-difference kernel in place
   !! of the spin-sum one. `(A - B)` does not move: it is exchange only, and
   !! exchange is same-spin. So the two halves become
   !!
   !!     (A+B)_T u = dEps u - c_x[(ab|ij)+(aj|ib)] u + 4 f_xc^T u
   !!     (A-B)_T u = (A-B) u
   !!
   !! which is one flag threaded from `response_core_t` through to the Fock
   !! build's `j_scale` and the kernel's spin channel, and every route above
   !! then works unchanged -- including `casida`, whose assumption is about
   !! `(A-B)` and so is as true for a triplet as for a singlet.
   !!
   !! A root at or below zero in the triplet manifold is not an excitation:
   !! it says the closed shell is a saddle point with respect to spin
   !! polarisation, and it is reported as an instability rather than as a
   !! small number.
   !!
   !! ## The three routes, and why the third exists
   !!
   !! `tda` is a Hermitian Davidson on `A`. `rpa` is the
   !! Stratmann-Scuseria-Frisch paired solver on the two halves. `casida` is
   !! the same spectrum as `rpa` reached a different way, and only for a
   !! **pure functional**: with no exact exchange `(A-B)` is exactly
   !! `diag(dEps)`, so the square root the paired solver has to take inside
   !! its subspace can be taken outright, and
   !!
   !!     dEps^{1/2} (A+B) dEps^{1/2} z = w^2 z
   !!
   !! is an ordinary symmetric eigenproblem for the same `w`. It costs one
   !! product per trial vector rather than two and it goes through a different
   !! solver, so it is an independent check on the paired one rather than a
   !! second spelling of it. It is not the production route: it says nothing
   !! about a hybrid, which is most of what anybody runs.
   !!
   !! ## Layout and normalisation
   !!
   !! Trial vectors are flat, `idx = (i-1)*n_vir + a` -- virtual fastest, which
   !! is `reshape` of the `(n_vir, n_occ)` rectangle and the layout
   !! `response_product`, `mqc_czt_ov_hessian` and the CI Davidson already
   !! share.
   !!
   !! **Every route returns `|X|^2 - |Y|^2 = 1/2`,** which in Tamm-Dancoff --
   !! where `Y` is zero -- reads `|X|^2 = 1/2`. It is the restricted
   !! closed-shell convention PySCF and Psi4 report, and the one the
   !! transition moments next door are written in. One convention rather than
   !! two, so nothing downstream has to ask which route produced what it was
   !! handed.
   !!
   !! ## What is not here
   !!
   !! A meta-GGA triplet kernel and an unrestricted reference. Transition
   !! moments, oscillator strengths and natural transition orbitals are in
   !! `mqc_czt_tddft_properties`, over the amplitudes this returns.
   !! `excited_decline_reason` in `mqc_czt_bridge` refuses what cannot be
   !! computed; the bridge refuses the rest by name rather than answering a
   !! different question.
   use pic_types, only: dp
   use pic_io, only: to_char
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   use mqc_program_limits, only: MAX_LINE_LENGTH
   use mqc_physical_constants, only: HARTREE_TO_EV
   use mqc_calculation_defaults, only: DEFAULT_RESPONSE_BATCH, DEFAULT_EXCITED_TOL, &
                                       STATE_SPIN_SINGLET, STATE_SPIN_TRIPLET, &
                                       STATE_SPIN_UNRESTRICTED
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_xc, only: xc_context_t, xc_kernel_cache_t, xc_kernel_cache_fill, &
                         xc_kernel_cache_uks_t, xc_kernel_cache_uks_fill
   use mqc_czt_response_product, only: response_product, response_product_uhf
   use mqc_davidson, only: davidson_flat, sigma_operator_t
   use mqc_czt_rpa_solver, only: paired_operator_t, rpa_solve, RPA_REASON_UNSTABLE_PLUS
   implicit none
   private

   public :: excitation_spectrum_t
   public :: response_core_t
   public :: tda_operator_t
   public :: rpa_operator_t
   public :: build_tda_operator
   public :: build_rpa_operator
   public :: tda_dense_matrix
   public :: rpa_dense_matrices
   public :: response_excitations
   public :: response_core_uhf_t
   public :: tda_operator_uhf_t
   public :: rpa_operator_uhf_t
   public :: build_tda_operator_uhf
   public :: build_rpa_operator_uhf
   public :: tda_dense_matrix_uhf
   public :: rpa_dense_matrices_uhf
   public :: response_excitations_uhf

   real(dp), parameter :: EXCITATION_FLOOR = 1.0e-3_dp
      !! Roots below this are not reported.
      !!
      !! A closed-shell singlet spectrum has nothing down here: the lowest
      !! excitation of even a small-gap molecule is tenths of a Hartree. What
      !! does appear is an artefact -- a rotation the reference is marginally
      !! unstable along, or a root the solver has not separated from zero --
      !! and neither is an excitation. The same threshold separates the
      !! degenerate block the guess is extended over.
      !!
      !! In the **triplet** manifold the same threshold means something
      !! stronger and is treated as an error rather than a filter: a root down
      !! here is the reference's own triplet instability.

   real(dp), parameter :: ROTATION_HINT = 5.0e-2_dp
      !! A reported root below this is flagged in the unrestricted state table.
      !!
      !! **A hint in the output, not a classification.** An open-shell
      !! reference carries rotations of its own singly-occupied orbitals, and
      !! the two routes disagree about them: the paired problem puts such a
      !! rotation at `omega^2` on the numerical zero and drops it below
      !! `EXCITATION_FLOOR`, while Tamm-Dancoff keeps it as a small positive
      !! root -- 6.7e-3 hartree on the OH radical -- because `A` alone is not
      !! the operator whose null space it lives in. So the same molecule can
      !! have a first Tamm-Dancoff state the paired spectrum does not have,
      !! and this marks the roots where that is worth checking. Nothing is
      !! dropped or relabelled on account of it; a genuine excitation this
      !! low, which a small-gap radical can have, is marked and reported.

   real(dp), parameter :: DEGENERACY_WINDOW = 1.0e-3_dp
      !! How close two orbital-energy gaps have to be for the guess to have to
      !! carry both.
      !!
      !! The starting vectors are unit vectors on the lowest gaps. Splitting a
      !! degenerate pair between the guess and the space outside it leaves the
      !! solver converging one partner against a subspace that cannot represent
      !! the other, so the count is extended over the whole degenerate block
      !! and the extra roots are found and then not reported.
      !!
      !! It is a window on the **gaps**, and so necessary rather than
      !! sufficient: the coupling reorders them, so a root of the spectrum can
      !! sit outside a guess sized on the gaps with nothing degenerate at the
      !! edge for this to catch. `GUESS_PER_ROOT` is what covers that.

   integer, parameter :: MAX_EXTRA_ROOTS = 8
      !! A cap on that extension. A highly symmetric molecule can put many gaps
      !! inside the window, and every one of them costs a converged root; eight
      !! covers a degeneracy no point group produces and stops a pathological
      !! case from turning a five-root request into a fifty-root solve.

   integer, parameter :: GUESS_PER_ROOT = 4
      !! Starting unit vectors per root converged.
      !!
      !! One per root is what a Hermitian Davidson starts from, and it is not
      !! enough here. The guess is picked on the gaps and what comes back is
      !! the spectrum; the coupling reorders the two, so the lowest `k` roots
      !! are not in general built on the `k` lowest gaps. A starting space
      !! that cannot reach one of them does not stall on it either: every
      !! root returned is converged to the tolerance asked for and is a true
      !! eigenpair, one of them simply a higher one, and nothing in the
      !! output says which.
      !!
      !! F2/cc-pVDZ asked for three singlets is the case in the suite. Its
      !! three lowest gaps are 0.759, 0.759 and 0.841 hartree; the 0.841 one
      !! carries the **fifth** root, and the pair carrying the third and
      !! fourth sits at 0.903, outside a three-vector guess and 0.06 hartree
      !! past anything `DEGENERACY_WINDOW` would have reached. The three
      !! roots come back as roots one, two and five. Asked for five the same
      !! solver finds all five, which is the signature of a starting space
      !! too narrow rather than of a solver defect.
      !!
      !! Four per root is what Psi4 uses; PySCF relies on the degeneracy
      !! window alone and is exposed to the same failure. The cost is
      !! matrix-vector products in the first iteration only -- the subspace
      !! that follows is driven by the residuals of the roots asked for, not
      !! by the width of the start.

   integer, parameter :: MAX_EXTRA_GUESS = 16
      !! A cap on the degeneracy extension of the guess, as `MAX_EXTRA_ROOTS`
      !! is of the root count. Wider, because a vector added here costs one
      !! product rather than a converged root.

   real(dp), parameter :: PYSCF_HARTREE_TO_EV = 27.21138602_dp
      !! What PySCF converts with, for the note the state table prints.
      !! This program uses `HARTREE_TO_EV`, which is the CODATA 2018 value;
      !! the two differ by 9.6e-10 relative, which is below the resolution of
      !! any excitation energy here but is the kind of thing a cross-code
      !! comparison in eV trips over.

   real(dp), parameter :: RHF_PAIRED_NORM = 0.5_dp
      !! `|X|^2 - |Y|^2` every route's amplitudes are handed back at.
      !!
      !! One half, not one. A closed-shell excitation is two spin-orbital
      !! excitations of equal weight and the spatial-orbital amplitude carries
      !! both, so the norm the spin-orbital problem sets to one comes to a half
      !! here. It is the convention PySCF and Psi4 report, and the one the
      !! transition moment `mu = 2 sum <i|r|a> (X+Y)` is written for; handing
      !! back unit-normalised amplitudes instead would put a factor of the
      !! square root of two into every oscillator strength downstream.
      !!
      !! Tamm-Dancoff has no `Y`, so the Davidson's unit eigenvector is scaled
      !! by the square root of this and the convention holds there too.

   type :: excitation_spectrum_t
      !! One converged spectrum, as `response_excitations` hands it back
      !!
      !! The four arrays are one object: they run over the same roots in the
      !! same order, they are allocated together, and three of them say
      !! nothing without the fourth. A consumer of the whole spectrum takes
      !! this rather than four separate arguments, and so cannot be handed a
      !! set of amplitudes belonging to a different solve than the energies.
      !!
      !! `response_excitations` still returns the four separately, because
      !! eleven call sites read them that way; packing is one assignment per
      !! array at the boundary.
      real(dp), allocatable :: excitations(:)
         !! (n_states) excitation energies, ascending, in Hartree
      integer, allocatable :: state_spin(:)
         !! (n_states) `STATE_SPIN_SINGLET` or `STATE_SPIN_TRIPLET` per root
      real(dp), allocatable :: x_amplitudes(:, :)
         !! (n_occ*n_vir, n_states) excitation amplitudes, virtual fastest
      real(dp), allocatable :: y_amplitudes(:, :)
         !! (n_occ*n_vir, n_states) de-excitation amplitudes; zero for `tda`
   end type excitation_spectrum_t

   type :: response_core_t
      !! Everything a response product needs, and nothing about which one
      !!
      !! `A`, `(A+B)`, `(A-B)` and the Casida reduction are four arithmetic
      !! rearrangements of the same two calls into `response_product`, so what
      !! they have in common -- the orbitals, the gaps, the exchange
      !! coefficients, the kernel cache, the batching and the flat-to-rectangle
      !! unpacking -- is here once, and each operator below is the few lines
      !! that are actually its own.
      !!
      !! The molecule and the exchange-correlation context are pointers because
      !! both outlive the solve and neither is cheap to copy; their targets
      !! have to outlive this object.
      type(czt_molecule_t), pointer :: mol => null()
      type(xc_context_t), pointer :: xc => null()
         !! Null for Hartree-Fock. Associated, `reference` is allocated too --
         !! the kernel is evaluated at a density, and a context without one
         !! reaches the Fock build with nothing to evaluate.
      real(dp), allocatable :: c_occ(:, :)     !! (n_ao, n_occ)
      real(dp), allocatable :: c_vir(:, :)     !! (n_ao, n_vir)
      real(dp), allocatable :: gaps(:, :)      !! (n_vir, n_occ), `e_a - e_i`
      real(dp), allocatable :: zero_h(:, :)    !! (n_ao, n_ao) of zeros
      real(dp), allocatable :: bounds(:, :)    !! Schwarz bounds
      real(dp), allocatable :: reference(:, :)  !! The converged density
      type(xc_kernel_cache_t) :: cache
         !! The kernel's coefficients over the whole grid, filled once by
         !! `build_response_core`. Unfilled for Hartree-Fock and ignored then.
      real(dp) :: k_scale = 1.0_dp
         !! Exact exchange the reference kept: one for Hartree-Fock, the
         !! mixing fraction for a hybrid, zero for a pure functional.
      real(dp) :: rs_k_lr = 0.0_dp
      real(dp) :: rs_omega = 0.0_dp
         !! A range-separated functional's attenuated second exchange pass.
      logical :: triplet = .false.
         !! Which manifold the next product belongs to. Everything else in
         !! this core -- the orbitals, the gaps, the exchange coefficients,
         !! the Schwarz screen, the filled cache -- is the same for both, so a
         !! solve over both manifolds flips this between them rather than
         !! building a second core and walking the quadrature again.
      integer :: n_occ = 0
      integer :: n_vir = 0
      integer :: batch = DEFAULT_RESPONSE_BATCH
         !! Trial vectors sharing one pass over the integrals.
      integer :: n_products = 0
         !! Applications of one half of the operator to one trial vector, for
         !! the cost line.
   contains
      procedure :: length => core_length
      procedure :: diagonal => core_diagonal
      procedure :: has_exchange => core_has_exchange
      procedure :: half_block => core_half_block
   end type response_core_t

   type, extends(sigma_operator_t) :: tda_operator_t
      !! `A` as something `mqc_davidson` will multiply a vector by
      type(response_core_t), allocatable :: core
         !! Allocated by the builder below, or moved in by `solve_manifold`.
         !!
         !! Allocatable and not a plain component so that lending a core to a
         !! route is a `move_alloc` and not an intrinsic assignment. The core
         !! carries `xc_kernel_cache_t` by value, whose thirteen grid-sized
         !! arrays an assignment deep-copies: on a solve over both manifolds
         !! that is a second copy of the whole quadrature alive for the length
         !! of the solve, three times over.
   contains
      procedure :: apply => tda_apply
      procedure :: apply_many => tda_apply_many
      procedure :: length => tda_length
      procedure :: diagonal => tda_diagonal
   end type tda_operator_t

   type, extends(paired_operator_t) :: rpa_operator_t
      !! `(A+B)` and `(A-B)` as something `mqc_czt_rpa_solver` will pair up
      type(response_core_t), allocatable :: core
         !! Allocated by the builder below, or moved in by `solve_manifold`.
         !!
         !! Allocatable and not a plain component so that lending a core to a
         !! route is a `move_alloc` and not an intrinsic assignment. The core
         !! carries `xc_kernel_cache_t` by value, whose thirteen grid-sized
         !! arrays an assignment deep-copies: on a solve over both manifolds
         !! that is a second copy of the whole quadrature alive for the length
         !! of the solve, three times over.
   contains
      procedure :: apply_plus => rpa_apply_plus
      procedure :: apply_minus => rpa_apply_minus
   end type rpa_operator_t

   type, extends(sigma_operator_t) :: casida_operator_t
      !! `dEps^{1/2}(A+B)dEps^{1/2}` for a pure functional, as a Davidson operator
      !!
      !! Not public: this is the cross-check route, reached through
      !! `response_excitations` with `method = "casida"` and refused for
      !! anything carrying exact exchange, where `(A-B)` is not the diagonal
      !! it assumes.
      type(response_core_t), allocatable :: core
         !! Allocated by the builder below, or moved in by `solve_manifold`.
         !!
         !! Allocatable and not a plain component so that lending a core to a
         !! route is a `move_alloc` and not an intrinsic assignment. The core
         !! carries `xc_kernel_cache_t` by value, whose thirteen grid-sized
         !! arrays an assignment deep-copies: on a solve over both manifolds
         !! that is a second copy of the whole quadrature alive for the length
         !! of the solve, three times over.
      real(dp), allocatable :: root_gaps(:)
         !! (n_ov) `sqrt(e_a - e_i)`, flat
   contains
      procedure :: apply => casida_apply
      procedure :: apply_many => casida_apply_many
   end type casida_operator_t

   type :: response_core_uhf_t
      !! Everything an unrestricted response product needs
      !!
      !! `response_core_t` with two of almost everything: two sets of
      !! orbitals, two gap rectangles, two reference densities and a kernel
      !! whose second derivative is spin resolved. A trial vector is the alpha
      !! rectangle flattened followed by the beta one, each virtual fastest.
      !!
      !! The molecule and the exchange-correlation context are pointers
      !! because both outlive the solve and neither is cheap to copy; their
      !! targets have to outlive this object.
      type(czt_molecule_t), pointer :: mol => null()
      type(xc_context_t), pointer :: xc => null()
         !! Null for Hartree-Fock. Associated, it is **spin-polarised** and
         !! both reference densities are allocated too.
      real(dp), allocatable :: c_occ_a(:, :), c_vir_a(:, :)   !! (n_ao, n_occ_a), (n_ao, n_vir_a)
      real(dp), allocatable :: c_occ_b(:, :), c_vir_b(:, :)
      real(dp), allocatable :: gaps_a(:, :), gaps_b(:, :)     !! (n_vir_s, n_occ_s), `e_a - e_i`
      real(dp), allocatable :: zero_h(:, :)                   !! (n_ao, n_ao) of zeros
      real(dp), allocatable :: bounds(:, :)                   !! Schwarz bounds
      real(dp), allocatable :: ref_a(:, :), ref_b(:, :)
         !! The converged spin densities, `C_sigma C_sigma^T` and not doubled
      type(xc_kernel_cache_uks_t) :: cache
         !! The polarised kernel's coefficients over the whole grid, filled
         !! once by `build_response_core_uhf`. Unfilled for Hartree-Fock.
      real(dp) :: k_scale = 1.0_dp
         !! Exact exchange the reference kept: one for Hartree-Fock, the
         !! mixing fraction for a hybrid, zero for a pure functional.
      real(dp) :: rs_k_lr = 0.0_dp
      real(dp) :: rs_omega = 0.0_dp
         !! A range-separated functional's attenuated second exchange pass.
      integer :: n_occ_a = 0
      integer :: n_vir_a = 0
      integer :: n_occ_b = 0
      integer :: n_vir_b = 0
      integer :: batch = DEFAULT_RESPONSE_BATCH
         !! Trial vectors sharing one pass over the integrals.
      integer :: n_products = 0
         !! Applications of one half of the operator to one trial vector.
   contains
      procedure :: length => core_uhf_length
      procedure :: diagonal => core_uhf_diagonal
      procedure :: has_exchange => core_uhf_has_exchange
      procedure :: half_block => core_uhf_half_block
      procedure :: destroy => core_uhf_destroy
   end type response_core_uhf_t

   type, extends(sigma_operator_t) :: tda_operator_uhf_t
      !! The unrestricted `A` as something `mqc_davidson` will multiply by
      type(response_core_uhf_t) :: core
   contains
      procedure :: apply => tda_uhf_apply
      procedure :: apply_many => tda_uhf_apply_many
      procedure :: length => tda_uhf_length
      procedure :: diagonal => tda_uhf_diagonal
   end type tda_operator_uhf_t

   type, extends(paired_operator_t) :: rpa_operator_uhf_t
      !! The unrestricted `(A+B)` and `(A-B)` for the paired solver
      type(response_core_uhf_t) :: core
   contains
      procedure :: apply_plus => rpa_uhf_apply_plus
      procedure :: apply_minus => rpa_uhf_apply_minus
   end type rpa_operator_uhf_t

contains

   pure function core_length(this) result(n)
      !! How long a trial vector is: `n_occ * n_vir`
      class(response_core_t), intent(in) :: this
      integer :: n

      n = this%n_occ*this%n_vir
   end function core_length

   function core_diagonal(this) result(diag)
      !! `e_a - e_i`, flat: what the solvers precondition on and start from
      !!
      !! Not `diag(A)`. The exact diagonal costs one Fock build per element,
      !! which is the whole matrix; the gaps are what PySCF preconditions
      !! with, and a preconditioner moves the iteration count rather than the
      !! eigenvalue.
      class(response_core_t), intent(in) :: this
      real(dp), allocatable :: diag(:)

      diag = reshape(this%gaps, [this%n_occ*this%n_vir])
   end function core_diagonal

   pure function core_has_exchange(this) result(yes)
      !! Whether `(A-B)` is anything but the orbital-energy diagonal
      !!
      !! `(A-B)` is exchange only. With no exchange of any range it is
      !! `diag(dEps)`, and building it would be an integral pass whose every
      !! quartet is multiplied by zero.
      class(response_core_t), intent(in) :: this
      logical :: yes

      yes = this%k_scale /= 0.0_dp .or. this%rs_k_lr /= 0.0_dp
   end function core_has_exchange

   subroutine core_half_block(this, vectors, images, minus, error)
      !! `(A+B)v` or `(A-B)v` for a block of flat trial vectors
      !!
      !! The block is what makes the batching worth having: one pass over the
      !! quartets and one walk of the quadrature serve every vector in it, so
      !! an iteration costs what a single product would if it were alone.
      !! `batch` caps the width because `build_fock_direct_many` is
      !! memory-bandwidth bound past a dozen or two densities.
      !!
      !! The two `response_product` calls differ only in whether `xc` and
      !! `reference` are passed, and they are passed together or not at all:
      !! neither a null pointer nor an unallocated array may reach a
      !! non-optional dummy.
      class(response_core_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)   !! (n_ov, n_vectors)
      real(dp), intent(out) :: images(:, :)   !! (n_ov, n_vectors)
      logical, intent(in) :: minus            !! `(A-B)` rather than `(A+B)`
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: u(:, :, :), au(:, :, :)
      integer, allocatable :: idx(:)
      integer :: n_ov, n_vec, width, first, last, w, m

      ! `images` is intent(out) and every exit below zeroes it; this one has
      ! to as well, or a caller that ignored an error it was already carrying
      ! reads whatever the array held.
      if (error%has_error()) then
         images = 0.0_dp
         return
      end if

      n_ov = this%length()
      n_vec = size(vectors, 2)
      if (n_ov < 1 .or. .not. associated(this%mol)) then
         call error%set(ERROR_VALIDATION, "the response operator was applied before "// &
                        "it was given a reference to apply itself over")
         images = 0.0_dp
         return
      end if
      if (size(vectors, 1) /= n_ov .or. size(images, 1) /= n_ov) then
         call error%set(ERROR_VALIDATION, "a trial vector handed to the response "// &
                        "operator is not the length of the occupied-virtual space")
         images = 0.0_dp
         return
      end if

      ! No integrals at all when the difference is the diagonal.
      if (minus .and. .not. this%has_exchange()) then
         do m = 1, n_vec
            images(:, m) = reshape(this%gaps, [n_ov])*vectors(:, m)
         end do
         this%n_products = this%n_products + n_vec
         return
      end if

      width = min(max(this%batch, 1), n_vec)
      allocate (u(this%n_vir, this%n_occ, width), au(this%n_vir, this%n_occ, width))
      allocate (idx(width))

      do first = 1, n_vec, width
         last = min(first + width - 1, n_vec)
         w = last - first + 1
         do m = 1, w
            u(:, :, m) = reshape(vectors(:, first + m - 1), [this%n_vir, this%n_occ])
            idx(m) = m
         end do

         au = 0.0_dp
         if (associated(this%xc)) then
            call response_product(this%mol, this%c_occ, this%c_vir, this%gaps, &
                                  this%zero_h, u, idx, w, minus, au, error, &
                                  direct=.true., bounds=this%bounds, &
                                  k_scale=this%k_scale, xc=this%xc, &
                                  reference=this%reference, rs_k_lr=this%rs_k_lr, &
                                  rs_omega=this%rs_omega, cache=this%cache, &
                                  triplet=this%triplet)
         else
            call response_product(this%mol, this%c_occ, this%c_vir, this%gaps, &
                                  this%zero_h, u, idx, w, minus, au, error, &
                                  direct=.true., bounds=this%bounds, &
                                  k_scale=this%k_scale, triplet=this%triplet)
         end if
         if (error%has_error()) exit

         do m = 1, w
            images(:, first + m - 1) = reshape(au(:, :, m), [n_ov])
         end do
         this%n_products = this%n_products + w
      end do

      if (error%has_error()) images = 0.0_dp
      deallocate (u, au, idx)
   end subroutine core_half_block

   pure function tda_length(this) result(n)
      !! How long a trial vector is: `n_occ * n_vir`
      class(tda_operator_t), intent(in) :: this
      integer :: n

      n = this%core%length()
   end function tda_length

   function tda_diagonal(this) result(diag)
      !! What the solver preconditions on and starts from
      class(tda_operator_t), intent(in) :: this
      real(dp), allocatable :: diag(:)

      diag = this%core%diagonal()
   end function tda_diagonal

   subroutine tda_apply_many(this, vectors, images, error)
      !! `A x = (1/2)[(A+B)x + (A-B)x]` for a block of trial vectors
      class(tda_operator_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)   !! (n_ov, n_vectors)
      real(dp), intent(out) :: images(:, :)   !! (n_ov, n_vectors)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: minus_image(:, :)

      if (error%has_error()) return

      call this%core%half_block(vectors, images, .false., error)
      if (error%has_error()) return

      allocate (minus_image(size(images, 1), size(images, 2)))
      call this%core%half_block(vectors, minus_image, .true., error)
      if (error%has_error()) then
         images = 0.0_dp
      else
         images = 0.5_dp*(images + minus_image)
      end if
      deallocate (minus_image)
   end subroutine tda_apply_many

   subroutine tda_apply(this, vector, image, error)
      !! `A x` for one vector, as a block of one
      class(tda_operator_t), intent(inout) :: this
      real(dp), intent(in) :: vector(:)
      real(dp), intent(out) :: image(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: block_in(:, :), block_out(:, :)

      allocate (block_in(size(vector), 1), block_out(size(image), 1))
      block_in(:, 1) = vector
      call this%apply_many(block_in, block_out, error)
      image = block_out(:, 1)
      deallocate (block_in, block_out)
   end subroutine tda_apply

   subroutine rpa_apply_plus(this, vectors, images, error)
      !! `(A+B)v` for a block of trial vectors
      class(rpa_operator_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)
      real(dp), intent(out) :: images(:, :)
      type(error_t), intent(inout) :: error

      call this%core%half_block(vectors, images, .false., error)
   end subroutine rpa_apply_plus

   subroutine rpa_apply_minus(this, vectors, images, error)
      !! `(A-B)v` for a block of trial vectors
      class(rpa_operator_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)
      real(dp), intent(out) :: images(:, :)
      type(error_t), intent(inout) :: error

      call this%core%half_block(vectors, images, .true., error)
   end subroutine rpa_apply_minus

   subroutine casida_apply_many(this, vectors, images, error)
      !! `dEps^{1/2}(A+B)dEps^{1/2} v` for a block of trial vectors
      class(casida_operator_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)
      real(dp), intent(out) :: images(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: scaled(:, :)
      integer :: m

      if (error%has_error()) return

      allocate (scaled(size(vectors, 1), size(vectors, 2)))
      do m = 1, size(vectors, 2)
         scaled(:, m) = this%root_gaps*vectors(:, m)
      end do
      call this%core%half_block(scaled, images, .false., error)
      if (error%has_error()) then
         images = 0.0_dp
      else
         do m = 1, size(images, 2)
            images(:, m) = this%root_gaps*images(:, m)
         end do
      end if
      deallocate (scaled)
   end subroutine casida_apply_many

   subroutine casida_apply(this, vector, image, error)
      !! The same for one vector, as a block of one
      class(casida_operator_t), intent(inout) :: this
      real(dp), intent(in) :: vector(:)
      real(dp), intent(out) :: image(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: block_in(:, :), block_out(:, :)

      allocate (block_in(size(vector), 1), block_out(size(image), 1))
      block_in(:, 1) = vector
      call this%apply_many(block_in, block_out, error)
      image = block_out(:, 1)
      deallocate (block_in, block_out)
   end subroutine casida_apply

   subroutine build_response_core(mol, orbitals, energies, n_occ, core, error, &
                                  xc, reference, bounds, batch, spin)
      !! Point a `response_core_t` at a converged closed-shell reference
      !!
      !! `mol` and `xc` are `target` and the core keeps pointers to them, so
      !! both have to outlive every product taken through it. `xc` and
      !! `reference` are one argument in two halves and are refused
      !! separately, for the reason `build_scf_ov_hessian` refuses them.
      !!
      !! The kernel cache is filled here rather than on first use: it is a
      !! property of the converged density and a response solve applies the
      !! kernel hundreds of times, on a quadrature that is most of a DFT run.
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals(:, :)     !! (n_ao, n_mo)
      real(dp), intent(in) :: energies(:)        !! (n_mo), ascending
      integer, intent(in) :: n_occ
      type(response_core_t), intent(out) :: core
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: reference(:, :)
         !! The converged SCF density the kernel is evaluated at.
      real(dp), intent(in), optional :: bounds(:, :)
         !! Schwarz bounds, computed here when the caller has none.
      integer, intent(in), optional :: batch
         !! Trial vectors per integral pass; `DEFAULT_RESPONSE_BATCH` absent.
      character(len=*), intent(in), optional :: spin
         !! Which manifold: `singlet` (the default), `triplet`, or `both`. It
         !! decides two things -- what the core applies now, and whether the
         !! kernel cache is filled with the triplet coefficients as well,
         !! which costs one more libxc pass over a grid whose basis functions
         !! are already in hand and is what makes `both` one quadrature rather
         !! than two.

      integer :: n_ao, n_mo, n_vir
      logical :: kernel_triplet
      character(len=16) :: manifold

      if (error%has_error()) return

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_vir = n_mo - n_occ

      if (n_occ < 1 .or. n_vir < 1) then
         call error%set(ERROR_VALIDATION, "an excitation needs at least one occupied "// &
                        "and one virtual orbital; there is nothing to excite between")
         return
      end if
      if (size(energies) < n_mo) then
         call error%set(ERROR_VALIDATION, "there are fewer orbital energies than "// &
                        "orbitals, so the response diagonal cannot be formed")
         return
      end if
      if (present(xc) .neqv. present(reference)) then
         call error%set(ERROR_VALIDATION, "the response operator was given an "// &
                        "exchange-correlation context without the reference density "// &
                        "its kernel is evaluated at, or the reverse; it needs both "// &
                        "or neither")
         return
      end if

      manifold = "singlet"
      if (present(spin)) manifold = spin
      select case (trim(manifold))
      case ("singlet")
         core%triplet = .false.
         kernel_triplet = .false.
      case ("triplet")
         core%triplet = .true.
         kernel_triplet = .true.
      case ("both")
         ! Starts on the singlets; the solve flips it for the second
         ! manifold, over the one cache that already holds both.
         core%triplet = .false.
         kernel_triplet = .true.
      case default
         call error%set(ERROR_VALIDATION, "the response operator was asked for the '"// &
                        trim(manifold)//"' manifold; it knows singlet, triplet and both")
         return
      end select
      if (kernel_triplet .and. present(xc)) then
         if (xc%any_mgga) then
            call error%set(ERROR_VALIDATION, "a triplet linear response over a "// &
                           "meta-GGA reference is not implemented: the triplet "// &
                           "kernel has no tau channels here, and the three it does "// &
                           "have would be a functional missing a term rather than "// &
                           "this one")
            return
         end if
      end if

      core%mol => mol
      core%n_occ = n_occ
      core%n_vir = n_vir
      core%c_occ = orbitals(:, 1:n_occ)
      core%c_vir = orbitals(:, n_occ + 1:n_mo)
      allocate (core%zero_h(n_ao, n_ao))
      core%zero_h = 0.0_dp
      call fill_gaps(energies, n_occ, n_mo, core%gaps)

      ! Here, and not after the Schwarz bounds and the kernel fill: it reads
      ! nothing but the gaps just formed, and the two it used to sit behind
      ! are the expensive part of building the operator.
      if (any(core%gaps <= 0.0_dp)) then
         call error%set(ERROR_VALIDATION, "an occupied orbital lies above a virtual "// &
                        "one, so these are not the aufbau orbitals and the gaps the "// &
                        "excitation solver preconditions on are not positive")
         return
      end if

      if (present(batch)) then
         if (batch > 0) core%batch = batch
      end if

      if (present(bounds)) then
         core%bounds = bounds
      else
         call schwarz_bounds(mol, core%bounds, error)
         if (error%has_error()) return
      end if

      if (present(xc)) then
         core%xc => xc
         core%reference = reference
         core%k_scale = xc%exx_fraction
         if (xc%range_separated) then
            core%rs_k_lr = xc%rs_k_lr
            core%rs_omega = xc%rs_omega
         end if
         call xc_kernel_cache_fill(xc, mol, reference, core%cache, error, &
                                   triplet=kernel_triplet)
         if (error%has_error()) return
      else
         core%k_scale = 1.0_dp
      end if
   end subroutine build_response_core

   subroutine build_tda_operator(mol, orbitals, energies, n_occ, operator, error, &
                                 xc, reference, bounds, batch, spin)
      !! Point a `tda_operator_t` at a converged closed-shell reference
      !!
      !! Everything the operator holds is in its `core`; this is the two lines
      !! that say which of the two halves the Davidson will be shown.
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals(:, :)     !! (n_ao, n_mo)
      real(dp), intent(in) :: energies(:)        !! (n_mo), ascending
      integer, intent(in) :: n_occ
      type(tda_operator_t), intent(out) :: operator
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: reference(:, :)
      real(dp), intent(in), optional :: bounds(:, :)
      integer, intent(in), optional :: batch
      character(len=*), intent(in), optional :: spin
         !! `singlet`, `triplet` or `both`; see `build_response_core`.

      allocate (operator%core)
      if (present(xc)) then
         call build_response_core(mol, orbitals, energies, n_occ, operator%core, &
                                  error, xc=xc, reference=reference, bounds=bounds, &
                                  batch=batch, spin=spin)
      else
         call build_response_core(mol, orbitals, energies, n_occ, operator%core, &
                                  error, bounds=bounds, batch=batch, spin=spin)
      end if
   end subroutine build_tda_operator

   subroutine build_rpa_operator(mol, orbitals, energies, n_occ, operator, error, &
                                 xc, reference, bounds, batch, spin)
      !! Point an `rpa_operator_t` at a converged closed-shell reference
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals(:, :)     !! (n_ao, n_mo)
      real(dp), intent(in) :: energies(:)        !! (n_mo), ascending
      integer, intent(in) :: n_occ
      type(rpa_operator_t), intent(out) :: operator
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: reference(:, :)
      real(dp), intent(in), optional :: bounds(:, :)
      integer, intent(in), optional :: batch
      character(len=*), intent(in), optional :: spin
         !! `singlet`, `triplet` or `both`; see `build_response_core`.

      allocate (operator%core)
      if (present(xc)) then
         call build_response_core(mol, orbitals, energies, n_occ, operator%core, &
                                  error, xc=xc, reference=reference, bounds=bounds, &
                                  batch=batch, spin=spin)
      else
         call build_response_core(mol, orbitals, energies, n_occ, operator%core, &
                                  error, bounds=bounds, batch=batch, spin=spin)
      end if
   end subroutine build_rpa_operator

   subroutine tda_dense_matrix(operator, a, error)
      !! The explicit `A`, from the operator applied to every unit vector
      !!
      !! `n_ov` matrix-vector products and an `n_ov` by `n_ov` array, so this
      !! is for a small case and for checking the solver against a dense
      !! diagonalisation -- which is how the operator is compared with another
      !! program element by element rather than through one summary of it.
      !! Nothing on the production path calls it.
      type(tda_operator_t), intent(inout) :: operator
      real(dp), allocatable, intent(out) :: a(:, :)
         !! (n_ov, n_ov), column `j` being `A e_j`
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: unit_vectors(:, :)
      integer :: n_ov

      if (error%has_error()) return

      n_ov = operator%length()
      call unit_matrix(n_ov, unit_vectors)
      allocate (a(n_ov, n_ov))
      call operator%apply_many(unit_vectors, a, error)
      deallocate (unit_vectors)
      if (error%has_error()) deallocate (a)
   end subroutine tda_dense_matrix

   subroutine rpa_dense_matrices(operator, aplus, aminus, error)
      !! The explicit `(A+B)` and `(A-B)`, from the operator the solver uses
      !!
      !! Two sets of `n_ov` products, for the same purpose as
      !! `tda_dense_matrix`: the dense reduction built from these is what the
      !! iterative solver is checked against, and building it out of the
      !! shipped operator rather than out of a second construction of the same
      !! physics is what makes the comparison say anything.
      type(rpa_operator_t), intent(inout) :: operator
      real(dp), allocatable, intent(out) :: aplus(:, :), aminus(:, :)
         !! (n_ov, n_ov) each, column `j` being the half applied to `e_j`
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: unit_vectors(:, :)
      integer :: n_ov

      if (error%has_error()) return

      n_ov = operator%core%length()
      call unit_matrix(n_ov, unit_vectors)
      allocate (aplus(n_ov, n_ov), aminus(n_ov, n_ov))
      call operator%apply_plus(unit_vectors, aplus, error)
      if (.not. error%has_error()) then
         call operator%apply_minus(unit_vectors, aminus, error)
      end if
      deallocate (unit_vectors)
      if (error%has_error()) deallocate (aplus, aminus)
   end subroutine rpa_dense_matrices

   subroutine unit_matrix(n, u)
      !! The `n` by `n` identity, as the block of every unit vector
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: u(:, :)

      integer :: j

      allocate (u(n, n))
      u = 0.0_dp
      do j = 1, n
         u(j, j) = 1.0_dp
      end do
   end subroutine unit_matrix

   function roots_to_solve(diagonal, n_states) result(n_solve)
      !! How many roots to converge so no degenerate block is cut in half
      !!
      !! `n_states`, extended over every gap within `DEGENERACY_WINDOW` of the
      !! `n_states`-th smallest and capped by `MAX_EXTRA_ROOTS` and by the
      !! space itself. The extras are converged and then not reported.
      real(dp), intent(in) :: diagonal(:)
      integer, intent(in) :: n_states
      integer :: n_solve

      real(dp) :: edge
      integer :: n, want

      n = size(diagonal)
      want = min(max(n_states, 1), n)
      edge = lowest_edge(diagonal, want)

      n_solve = count(diagonal <= edge + DEGENERACY_WINDOW)
      n_solve = min(n_solve, want + MAX_EXTRA_ROOTS, n)
      n_solve = max(n_solve, want)
   end function roots_to_solve

   function guess_count(diagonal, n_solve) result(n_guess)
      !! How many unit vectors the starting space needs
      !!
      !! `GUESS_PER_ROOT` per root converged, extended over every gap within
      !! `DEGENERACY_WINDOW` of the last one taken, and capped by
      !! `MAX_EXTRA_GUESS` and by the space itself. Counted before anything is
      !! allocated, because the subspace cap is sized from it.
      real(dp), intent(in) :: diagonal(:)
      integer, intent(in) :: n_solve
      integer :: n_guess

      real(dp) :: edge
      integer :: n, want

      n = size(diagonal)
      want = min(GUESS_PER_ROOT*max(n_solve, 1), n)
      edge = lowest_edge(diagonal, want)

      n_guess = min(count(diagonal <= edge + DEGENERACY_WINDOW), &
                    want + MAX_EXTRA_GUESS, n)
      n_guess = max(n_guess, want)
   end function guess_count

   function lowest_edge(diagonal, k) result(edge)
      !! The `k`-th smallest diagonal element
      !!
      !! A partial selection rather than a sort: `k` is a small multiple of
      !! the roots asked for and the space is not.
      real(dp), intent(in) :: diagonal(:)
      integer, intent(in) :: k
      real(dp) :: edge

      logical, allocatable :: taken(:)
      integer :: i, pick

      allocate (taken(size(diagonal)))
      taken = .false.
      edge = 0.0_dp
      do i = 1, k
         pick = lowest_free(diagonal, taken)
         taken(pick) = .true.
         edge = diagonal(pick)
      end do
      deallocate (taken)
   end function lowest_edge

   function lowest_free(diagonal, taken) result(pick)
      !! The smallest diagonal element not already spoken for
      real(dp), intent(in) :: diagonal(:)
      logical, intent(in) :: taken(:)
      integer :: pick

      integer :: i

      pick = 0
      do i = 1, size(diagonal)
         if (taken(i)) cycle
         if (pick == 0) then
            pick = i
         else if (diagonal(i) < diagonal(pick)) then
            pick = i
         end if
      end do
   end function lowest_free

   subroutine fill_guess(diagonal, n_guess, guess)
      !! Unit vectors on the `n_guess` lowest gaps
      real(dp), intent(in) :: diagonal(:)
      integer, intent(in) :: n_guess
      real(dp), intent(out) :: guess(:, :)
         !! (n_ov, n_guess), one starting vector a column

      logical, allocatable :: taken(:)
      integer :: k, pick

      allocate (taken(size(diagonal)))
      taken = .false.
      guess = 0.0_dp
      do k = 1, n_guess
         pick = lowest_free(diagonal, taken)
         taken(pick) = .true.
         guess(pick, k) = 1.0_dp
      end do
      deallocate (taken)
   end subroutine fill_guess

   subroutine solve_manifold(core, route, is_triplet, n_states, tol, subspace, &
                             values, x, y, error, max_iter, verbose)
      !! One manifold's excitation energies, with the artefacts already dropped
      !!
      !! The route's operator is built by copying `core` into it, rather than
      !! by constructing a second one: the orbitals, the Schwarz screen and
      !! the filled kernel cache are the same for both spins and for all three
      !! routes, and re-deriving them per manifold would be a second walk of
      !! the quadrature for nothing. The product count is carried back out so
      !! the caller's cost line covers every manifold it asked for.
      !!
      !! **A triplet root at the floor is an error, not a filter.** A singlet
      !! root that converges near zero is a rotation of the reference and is
      !! dropped with a warning, which is the right answer for a spectrum that
      !! simply has fewer states than were asked for. A triplet one says the
      !! closed shell is unstable against spin polarisation -- there is a
      !! lower-energy unrestricted solution, and every root of this operator
      !! is an expansion about a saddle point. Reporting the rest would be a
      !! spectrum of a reference nobody should be using.
      !!
      !! **The three routes reach that verdict in the same order.** Each
      !! returns its solver's own error first, then the shared convergence
      !! check, and only then the triplet floor -- so a solve that merely ran
      !! out of iterations is reported as that and not as an instability. The
      !! paired route's third mechanism, an imaginary frequency, is not an
      !! exception to the ordering: `rpa_solve` refuses a negative squared
      !! frequency itself, which is a statement about `(A+B)` on the subspace
      !! and true whether or not the solve had converged, and
      !! `name_the_instability` relabels only that one reason.
      type(response_core_t), allocatable, intent(inout) :: core
         !! Lent to the route's operator for the length of the solve and
         !! given back, rather than copied into it: see `tda_operator_t`.
      character(len=*), intent(in) :: route     !! `tda`, `rpa` or `casida`
      logical, intent(in) :: is_triplet
      integer, intent(in) :: n_states
      real(dp), intent(in) :: tol
      integer, intent(in) :: subspace
         !! Trial vectors kept before a collapse; zero takes each solver's
         !! own rule.
      real(dp), allocatable, intent(out) :: values(:)
         !! (n_found) excitation energies, ascending, in Hartree
      real(dp), allocatable, intent(out) :: x(:, :), y(:, :)
         !! (n_ov, n_found) the amplitudes, at the route's own normalisation
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      logical, intent(in), optional :: verbose

      type(tda_operator_t) :: tda
      type(rpa_operator_t) :: rpa
      type(casida_operator_t) :: casida
      real(dp), allocatable :: diagonal(:), raw(:), vectors(:, :), residuals(:)
      real(dp), allocatable :: xpy(:, :), xmy(:, :), all_x(:, :), all_y(:, :)
      real(dp), allocatable :: guess(:, :)
      character(len=:), allocatable :: manifold
      integer :: n_ov, n_solve, n_guess, iterations, products, n_found, k, keep, why
      integer :: cap
      logical :: converged

      if (error%has_error()) return

      core%triplet = is_triplet
      manifold = trim(manifold_word(is_triplet))
      diagonal = core%diagonal()
      n_ov = size(diagonal)
      n_solve = roots_to_solve(diagonal, n_states)
      n_guess = guess_count(diagonal, n_solve)
      allocate (all_x(n_ov, n_solve), all_y(n_ov, n_solve))

      ! The Davidson routes start on `GUESS_PER_ROOT` unit vectors a root
      ! rather than one, because the guess is picked on the gaps and what
      ! comes back is the spectrum, and the coupling reorders the two: see
      ! `GUESS_PER_ROOT`. The extra columns widen the starting subspace and
      ! do not ask for more roots. The paired route needs nothing here --
      ! `rpa_solve` builds a start of the same width out of its own copy of
      ! these rules.
      cap = davidson_subspace(subspace, n_solve, n_guess, n_ov)
      n_guess = min(n_guess, cap)

      ! Each route takes the core, works, and hands it back. `move_alloc`
      ! rather than assignment, and every branch below therefore runs to the
      ! end of its `case` instead of returning out of it -- a `return` from
      ! the middle would leave the caller's core unallocated for the second
      ! manifold.
      select case (trim(route))
      case ("rpa")
         call move_alloc(core, rpa%core)
         call rpa_solve(rpa, diagonal, n_solve, raw, xpy, xmy, residuals, &
                        iterations, products, converged, error, tolerance=tol, &
                        max_iterations=max_iter, max_subspace=subspace, &
                        verbose=verbose, label=manifold//" RPA iterations", &
                        reason=why)
         if (error%has_error()) then
            call name_the_instability(is_triplet, why, error)
         else
            ! `xpy . xmy = 1` out of the solver; the closed-shell convention
            ! is a half of that, and both vectors take the same factor so
            ! their half sum and half difference are `X` and `Y`.
            do k = 1, n_solve
               all_x(:, k) = sqrt(RHF_PAIRED_NORM)*0.5_dp*(xpy(:, k) + xmy(:, k))
               all_y(:, k) = sqrt(RHF_PAIRED_NORM)*0.5_dp*(xpy(:, k) - xmy(:, k))
            end do
         end if
         call move_alloc(rpa%core, core)
      case ("casida")
         call move_alloc(core, casida%core)
         if (casida%core%has_exchange()) then
            call error%set(ERROR_VALIDATION, "the Casida reduction assumes (A-B) is "// &
                           "the orbital-energy diagonal, which holds only for a "// &
                           "functional carrying no exact exchange; this reference "// &
                           "keeps a fraction of it, so ask for 'rpa'")
         else
            casida%root_gaps = sqrt(diagonal)
            ! The eigenvalue is `w^2` and so is the preconditioner: `dEps^2`
            ! is the diagonal of the reduced operator up to its two-electron
            ! part, the way `dEps` is the diagonal of `A`.
            allocate (guess(n_ov, n_guess))
            call fill_guess(diagonal, n_guess, guess)
            call davidson_flat(casida, diagonal*diagonal, n_solve, raw, vectors, &
                               residuals, iterations, products, converged, error, &
                               tolerance=tol, max_iterations=max_iter, &
                               max_subspace=cap, guess=guess, &
                               verbose=verbose, label=manifold//" Casida iterations", &
                               value_label="omega^2")
            deallocate (guess)
            if (.not. error%has_error()) then
               call casida_amplitudes(raw, vectors, diagonal, all_x, all_y)
               ! A negative `w^2` is the instability arriving as an imaginary
               ! frequency rather than as a failed factorisation, so it is
               ! clamped to zero here and caught by the floor below.
               raw = sqrt(max(raw, 0.0_dp))
            end if
         end if
         call move_alloc(casida%core, core)
      case default
         call move_alloc(core, tda%core)
         allocate (guess(n_ov, n_guess))
         call fill_guess(diagonal, n_guess, guess)
         call davidson_flat(tda, diagonal, n_solve, raw, vectors, residuals, &
                            iterations, products, converged, error, tolerance=tol, &
                            max_iterations=max_iter, max_subspace=cap, guess=guess, &
                            verbose=verbose, label=manifold//" TDA iterations", &
                            value_label="excitation")
         deallocate (guess)
         if (.not. error%has_error()) then
            ! The Davidson returns unit eigenvectors; the closed-shell
            ! convention is `|X|^2 = 1/2`, and scaling here rather than at
            ! every reader is what makes the three routes one convention.
            all_x = sqrt(RHF_PAIRED_NORM)*vectors
            all_y = 0.0_dp
         end if
         call move_alloc(tda%core, core)
      end select
      if (error%has_error()) return

      if (.not. converged) then
         call error%set(ERROR_GENERIC, "the "//manifold//" "//trim(route)//" solve "// &
                        "did not converge its roots in the iterations allowed; raise "// &
                        "keywords.excited_states.max_iter or loosen "// &
                        "keywords.excited_states.tolerance")
         return
      end if

      if (is_triplet .and. raw(1) <= EXCITATION_FLOOR) then
         call unstable_reference(raw(1), error)
         return
      end if

      ! Drop the artefacts, then the degeneracy padding, keeping the order.
      keep = 0
      do k = 1, n_solve
         if (raw(k) > EXCITATION_FLOOR) keep = keep + 1
      end do
      n_found = min(keep, n_states)
      allocate (values(n_found), x(n_ov, n_found), y(n_ov, n_found))
      keep = 0
      do k = 1, n_solve
         if (raw(k) <= EXCITATION_FLOOR) cycle
         keep = keep + 1
         if (keep > n_found) exit
         values(keep) = raw(k)
         x(:, keep) = all_x(:, k)
         y(:, keep) = all_y(:, k)
      end do

      ! Name every root the filter removed, with its value, and separate the
      ! two reasons one lands below the floor. A small positive root is a
      ! rotation the solver has not resolved from zero. A *negative* one is a
      ! different statement about the reference, not about the solver: the
      ! manifold's `A` is not positive definite, so the closed shell is a
      ! saddle point and every energy reported above it is an excitation of
      ! an unstable reference. Reported whether or not enough roots survived,
      ! because it is a fact about the ground state rather than about the
      ! request. A triplet never reaches here with a root at the floor --
      ! `unstable_reference` has already refused above.
      do k = 1, n_solve
         if (raw(k) > EXCITATION_FLOOR) cycle
         if (raw(k) < 0.0_dp) then
            call logger%warning("  "//manifold//" root "//to_char(k)// &
                                " converged to "//to_char(raw(k))//" hartree, which "// &
                                "is negative: the response matrix is not positive "// &
                                "definite, so this reference is a saddle point "// &
                                "rather than a minimum and the spectrum below is "// &
                                "taken from an unstable reference. "// &
                                "keywords.scf.stability examines the reference "// &
                                "itself.")
         else
            call logger%warning("  "//manifold//" root "//to_char(k)// &
                                " converged to "//to_char(raw(k))//" hartree, below "// &
                                "the "//to_char(EXCITATION_FLOOR)//" hartree floor: "// &
                                "a rotation of the reference the solver has not "// &
                                "separated from zero, not an excitation.")
         end if
      end do
      if (n_found < n_states) then
         call logger%warning("  only "//to_char(n_found)//" of the "// &
                             to_char(n_states)//" "//manifold//" roots asked for are "// &
                             "excitations; the rest converged below "// &
                             to_char(EXCITATION_FLOOR)//" hartree and are rotations "// &
                             "of the reference, not excited states")
      end if
   end subroutine solve_manifold

   pure function manifold_word(is_triplet) result(word)
      !! `singlet` or `triplet`, for a message
      logical, intent(in) :: is_triplet
      character(len=7) :: word

      word = "singlet"
      if (is_triplet) word = "triplet"
   end function manifold_word

   subroutine unstable_reference(lowest, error)
      !! Report a non-positive triplet root as what it is
      !!
      !! Concatenated rather than written into a buffer: the message is longer
      !! than an output record, and a format-directed write that overflows one
      !! is a run-time failure rather than a truncation.
      real(dp), intent(in) :: lowest   !! The offending root, in Hartree
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: what

      what = "its lowest triplet root is "//to_char(lowest)//" hartree, at or "// &
             "below zero"
      call error%set(ERROR_GENERIC, "the reference is triplet-unstable: "//what// &
                     ", so this closed shell is a saddle point against spin "// &
                     "polarisation and a lower unrestricted solution exists. What "// &
                     "the solver would report is an expansion about that saddle "// &
                     "point rather than an excitation spectrum; converge an "// &
                     "unrestricted reference.")
   end subroutine unstable_reference

   subroutine name_the_instability(is_triplet, reason, error)
      !! Say what a failed paired solve means when the manifold is a triplet
      !!
      !! An instability of `(A+B)` reaches the solver as a negative squared
      !! frequency, and the solver says so in its own terms. That is the right
      !! diagnosis of the arithmetic and the wrong one for a reader, who has
      !! asked for a spectrum and wants to know that the *reference* is what
      !! is wrong, and in which manifold. The solver's message is kept and
      !! prefixed rather than replaced.
      !!
      !! **Only that one reason is relabelled**, and it is read from `reason`
      !! rather than from the message text. `(A+B)` is where the Coulomb term
      !! and the exchange-correlation kernel live, so it is the half that
      !! differs between the manifolds and the only one a triplet run may
      !! claim. `(A-B)` is exchange alone -- the identical operator for both
      !! spins -- so calling its failure a triplet instability would tell the
      !! reader to converge an unrestricted reference on the strength of the
      !! manifold they happened to ask for, which is a diagnosis invented
      !! rather than read off. The same goes for every other way the solve
      !! fails: a LAPACK error, an exhausted subspace, a solve that ran out of
      !! iterations.
      logical, intent(in) :: is_triplet
      integer, intent(in) :: reason   !! One of the solver's `RPA_REASON_*`
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: was

      if (.not. is_triplet) return
      if (.not. error%has_error()) return
      if (reason /= RPA_REASON_UNSTABLE_PLUS) return
      was = error%get_message()
      call error%set(ERROR_GENERIC, "the reference is triplet-unstable: the paired "// &
                     "solve of the triplet manifold found an imaginary excitation "// &
                     "energy, which for a closed shell means it is a saddle point "// &
                     "against spin polarisation. The solver reported: "//was)
   end subroutine name_the_instability

   subroutine response_excitations(mol, orbitals, energies, n_occ, n_states, method, &
                                   spin, excitations, state_spin, x_amplitudes, &
                                   y_amplitudes, error, xc, reference, bounds, &
                                   tolerance, max_iter, max_subspace, batch, verbose)
      !! The lowest excitation energies of a closed shell
      !!
      !! What comes back is ascending, in Hartree above the reference, with
      !! every root below `EXCITATION_FLOOR` already dropped -- so `size` of
      !! it can be smaller than `n_states`, and a caller has to read the size
      !! rather than assume it. Always allocated when this returns without an
      !! error, empty included.
      !!
      !! **Normalisation does not depend on `method`.** All three routes
      !! return `sum(X^2) - sum(Y^2) = 1/2`, the restricted closed-shell
      !! convention. `tda` has no `Y` in its approximation, so its `Y` is
      !! zeros -- which says there is none rather than standing for something
      !! left unfilled -- and its `X` alone carries the half.
      !!
      !! With `spin = "both"` the two manifolds are solved over one core and
      !! one filled kernel cache, and the results are **interleaved by
      !! energy** rather than concatenated: `state_spin` is what says which
      !! root is which, and there are then up to `2 * n_states` of them.
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      integer, intent(in) :: n_states
         !! Roots per manifold, not in total.
      character(len=*), intent(in) :: method
         !! `tda`, `rpa`, or `casida` for the cross-check reduction, which
         !! needs a functional carrying no exact exchange.
      character(len=*), intent(in) :: spin
         !! `singlet`, `triplet` or `both`.
      real(dp), allocatable, intent(out) :: excitations(:)
         !! (n_found) excitation energies, ascending, in Hartree
      integer, allocatable, intent(out) :: state_spin(:)
         !! (n_found) `STATE_SPIN_SINGLET` or `STATE_SPIN_TRIPLET` per root
      real(dp), allocatable, intent(out) :: x_amplitudes(:, :)
         !! (n_occ*n_vir, n_found) excitation amplitudes, virtual fastest
      real(dp), allocatable, intent(out) :: y_amplitudes(:, :)
         !! (n_occ*n_vir, n_found) de-excitation amplitudes; zero for `tda`
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: reference(:, :)
      real(dp), intent(in), optional :: bounds(:, :)
      real(dp), intent(in), optional :: tolerance
         !! Residual norm a root is accepted at; `DEFAULT_EXCITED_TOL` absent.
      integer, intent(in), optional :: max_iter
      integer, intent(in), optional :: max_subspace
         !! Trial vectors kept before the subspace is collapsed. Zero or
         !! absent sizes it from the starting space, which is wider than the
         !! rule each solver would apply on its own; it is decided here rather
         !! than forwarded because a deck's unset `max_subspace` is a zero and
         !! a solver would read that as a subspace of no vectors. Set below
         !! the starting space, it truncates the guess.
      integer, intent(in), optional :: batch
      logical, intent(in), optional :: verbose
         !! A line per iteration. Each one is an integral pass.

      type(response_core_t), allocatable :: core
      real(dp), allocatable :: e_singlet(:), xs(:, :), ys(:, :)
      real(dp), allocatable :: e_triplet(:), xt(:, :), yt(:, :)
      real(dp) :: tol
      character(len=MAX_LINE_LENGTH) :: line
      character(len=16) :: route
      integer :: n_ov, subspace
      logical :: want_singlet, want_triplet

      if (error%has_error()) return
      ! Allocated empty rather than left unallocated: the contract above is
      ! that a caller reads the size, and a caller doing that on an
      ! unallocated array has no way to notice.
      if (n_states < 1) then
         allocate (excitations(0), state_spin(0), x_amplitudes(0, 0))
         allocate (y_amplitudes(0, 0))
         return
      end if

      route = trim(adjustl(method))
      if (route /= "tda" .and. route /= "rpa" .and. route /= "casida") then
         call error%set(ERROR_VALIDATION, "'"//trim(route)//"' is not a linear-"// &
                        "response method this backend knows; it has 'tda', 'rpa' "// &
                        "and the 'casida' cross-check")
         return
      end if
      select case (trim(spin))
      case ("singlet")
         want_singlet = .true.
         want_triplet = .false.
      case ("triplet")
         want_singlet = .false.
         want_triplet = .true.
      case ("both")
         want_singlet = .true.
         want_triplet = .true.
      case default
         call error%set(ERROR_VALIDATION, "keywords.excited_states.spin is '"// &
                        trim(spin)//"'; this solver knows singlet, triplet and both")
         return
      end select

      allocate (core)
      if (present(xc)) then
         call build_response_core(mol, orbitals, energies, n_occ, core, error, &
                                  xc=xc, reference=reference, bounds=bounds, &
                                  batch=batch, spin=trim(spin))
      else
         call build_response_core(mol, orbitals, energies, n_occ, core, error, &
                                  bounds=bounds, batch=batch, spin=trim(spin))
      end if
      if (error%has_error()) return

      n_ov = core%length()
      if (n_states > n_ov) then
         call error%set(ERROR_VALIDATION, "keywords.excited_states asked for "// &
                        to_char(n_states)//" roots, but this reference has only "// &
                        to_char(n_ov)//" single excitations to build them from")
         return
      end if

      tol = DEFAULT_EXCITED_TOL
      if (present(tolerance)) tol = tolerance
      subspace = 0
      if (present(max_subspace)) then
         if (max_subspace > 0) subspace = max_subspace
      end if

      if (want_singlet) then
         call solve_manifold(core, route, .false., n_states, tol, subspace, &
                             e_singlet, xs, ys, error, max_iter=max_iter, &
                             verbose=verbose)
         if (error%has_error()) return
      else
         allocate (e_singlet(0), xs(n_ov, 0), ys(n_ov, 0))
      end if
      if (want_triplet) then
         call solve_manifold(core, route, .true., n_states, tol, subspace, &
                             e_triplet, xt, yt, error, max_iter=max_iter, &
                             verbose=verbose)
         if (error%has_error()) return
      else
         allocate (e_triplet(0), xt(n_ov, 0), yt(n_ov, 0))
      end if

      call merge_manifolds(e_singlet, xs, ys, e_triplet, xt, yt, excitations, &
                           state_spin, x_amplitudes, y_amplitudes)

      write (line, "(a,a,a,i0,a,i0,a)") "  ", route_name(route), ": ", &
         size(excitations), " root(s) from ", core%n_products, &
         " matrix-vector products"
      call logger%info(trim(line))
      call log_state_table(route, excitations, state_spin, x_amplitudes, &
                           y_amplitudes, n_occ, n_ov/n_occ)
   end subroutine response_excitations

   subroutine merge_manifolds(e_singlet, xs, ys, e_triplet, xt, yt, excitations, &
                              state_spin, x_amplitudes, y_amplitudes)
      !! The two manifolds as one spectrum, ascending
      !!
      !! Both lists arrive sorted, so this is the merge step of a merge sort
      !! and nothing here is quadratic. A tie goes to the singlet, which only
      !! decides a print order: two roots that close cannot be told apart by
      !! their energies anyway, and `state_spin` says which is which.
      real(dp), intent(in) :: e_singlet(:), e_triplet(:)
      real(dp), intent(in) :: xs(:, :), ys(:, :), xt(:, :), yt(:, :)
      real(dp), allocatable, intent(out) :: excitations(:)
      integer, allocatable, intent(out) :: state_spin(:)
      real(dp), allocatable, intent(out) :: x_amplitudes(:, :), y_amplitudes(:, :)

      integer :: ns, nt, n_ov, is, it, k
      logical :: take_singlet

      ns = size(e_singlet)
      nt = size(e_triplet)
      n_ov = max(size(xs, 1), size(xt, 1))
      allocate (excitations(ns + nt), state_spin(ns + nt))
      allocate (x_amplitudes(n_ov, ns + nt), y_amplitudes(n_ov, ns + nt))

      is = 1
      it = 1
      do k = 1, ns + nt
         if (it > nt) then
            take_singlet = .true.
         else if (is > ns) then
            take_singlet = .false.
         else
            take_singlet = e_singlet(is) <= e_triplet(it)
         end if
         if (take_singlet) then
            excitations(k) = e_singlet(is)
            state_spin(k) = STATE_SPIN_SINGLET
            x_amplitudes(:, k) = xs(:, is)
            y_amplitudes(:, k) = ys(:, is)
            is = is + 1
         else
            excitations(k) = e_triplet(it)
            state_spin(k) = STATE_SPIN_TRIPLET
            x_amplitudes(:, k) = xt(:, it)
            y_amplitudes(:, k) = yt(:, it)
            it = it + 1
         end if
      end do
   end subroutine merge_manifolds

   function davidson_subspace(requested, n_solve, n_guess, n_ov) result(subspace)
      !! The Davidson's subspace cap, with a deck's unset zero resolved
      !!
      !! `davidson_flat` takes the cap as a plain integer and would read a
      !! zero as a subspace of no vectors, so its own default rule is restated
      !! here rather than forwarded. Shared by the Tamm-Dancoff and Casida
      !! routes, which both go through that solver; the paired solver takes
      !! the zero itself and has a rule of its own.
      !!
      !! Sized from the **starting space** plus the room the solver's own rule
      !! leaves for expansion, `n_solve` vectors an iteration and eight over.
      !! Taking `n_guess` rather than `n_solve` as the base is what keeps a
      !! start wider than the roots asked for from being collapsed away after
      !! one iteration. A cap the caller asked for is honoured as given, and
      !! set below the starting space it truncates the guess.
      integer, intent(in) :: requested   !! Zero for the default
      integer, intent(in) :: n_solve
      integer, intent(in) :: n_guess     !! Columns the starting space wants
      integer, intent(in) :: n_ov
      integer :: subspace

      subspace = max(n_guess + n_solve + 8, 16)
      if (requested > 0) subspace = requested
      subspace = max(min(subspace, n_ov), n_solve)
   end function davidson_subspace

   subroutine casida_amplitudes(w2, z, gaps, x, y)
      !! `X` and `Y` from the Casida eigenvectors
      !!
      !! `R = X+Y` is `dEps^{1/2} z` and, with `(A-B) = diag(dEps)`, the
      !! partner follows without another product: `(A-B)(X-Y) = w (X+Y)` reads
      !! `L_i = w R_i / dEps_i`. Scaled to `R . L = 1` and then to the
      !! closed-shell half, so these come back on the same footing as the
      !! paired solver's.
      real(dp), intent(in) :: w2(:)      !! (n_roots) the eigenvalues, `w^2`
      real(dp), intent(in) :: z(:, :)    !! (n_ov, n_roots) the eigenvectors
      real(dp), intent(in) :: gaps(:)    !! (n_ov) `dEps`
      real(dp), intent(out) :: x(:, :), y(:, :)

      real(dp), allocatable :: r(:), l(:)
      real(dp) :: w, inner
      integer :: k

      allocate (r(size(gaps)), l(size(gaps)))
      do k = 1, size(w2)
         w = sqrt(max(w2(k), 0.0_dp))
         r = sqrt(gaps)*z(:, k)
         if (w > 0.0_dp) then
            l = w*r/gaps
         else
            l = 0.0_dp
         end if
         inner = dot_product(r, l)
         if (inner > 0.0_dp) then
            r = r/sqrt(inner)
            l = l/sqrt(inner)
         end if
         x(:, k) = sqrt(RHF_PAIRED_NORM)*0.5_dp*(r + l)
         y(:, k) = sqrt(RHF_PAIRED_NORM)*0.5_dp*(r - l)
      end do
      deallocate (r, l)
   end subroutine casida_amplitudes

   pure function route_name(route) result(name)
      !! What a route is called in a log line
      character(len=*), intent(in) :: route
      character(len=:), allocatable :: name

      select case (trim(route))
      case ("tda")
         name = "Tamm-Dancoff"
      case ("rpa")
         name = "RPA"
      case default
         name = "RPA (Casida)"
      end select
   end function route_name

   subroutine log_state_table(route, excitations, state_spin, x, y, n_occ, n_vir)
      !! The spectrum, as a table with the orbitals each root is made of
      !!
      !! Amplitudes below `AMPLITUDE_FLOOR` are left out: a converged root of
      !! a molecule of any size has a handful of contributions above it and
      !! hundreds of numerical dust below, and printing the dust hides the
      !! assignment the table exists for. Orbitals are numbered from one over
      !! all of them, occupied and virtual together, so the indices match what
      !! every other table in this program prints.
      !!
      !! The `|X|^2-|Y|^2` column is the norm the paired problem conserves
      !! and the Tamm-Dancoff one is scaled to, and it is printed rather than
      !! asserted: a root whose value has drifted off the convention is one
      !! whose biorthonormalisation did not take, and every transition moment
      !! built from it afterwards would be wrong by that factor.
      !!
      !! The amplitudes printed are the ones stored, so a dominant single
      !! excitation reads about 0.707 rather than 1 -- `AMPLITUDE_FLOOR` is a
      !! floor on a vector of norm `sqrt(1/2)`, not on a unit one.
      character(len=*), intent(in) :: route
      real(dp), intent(in) :: excitations(:)
      integer, intent(in) :: state_spin(:)
         !! `STATE_SPIN_*` per root, printed as a column: with `spin = "both"`
         !! the two manifolds interleave, and a table of energies alone would
         !! not say which row is which.
      real(dp), intent(in) :: x(:, :), y(:, :)
      integer, intent(in) :: n_occ, n_vir

      real(dp), parameter :: AMPLITUDE_FLOOR = 0.1_dp
      character(len=MAX_LINE_LENGTH) :: line, piece
      real(dp) :: weight
      integer :: k, i, a, idx

      if (size(excitations) < 1) return

      ! Said once, beside the numbers it applies to: the two codes a
      ! cross-check is run against convert with a different Hartree.
      if (abs(HARTREE_TO_EV - PYSCF_HARTREE_TO_EV) > 0.0_dp) then
         write (line, "(a,f16.12,a,f14.8,a)") "  excitation energies in eV use ", &
            HARTREE_TO_EV, " eV/hartree; PySCF uses ", PYSCF_HARTREE_TO_EV, &
            ", which differs in the eighth decimal"
         call logger%info(trim(line))
      end if
      call logger%info("  "//route_name(route)//" amplitudes are normalised to "// &
                       "|X|^2 - |Y|^2 = 0.5")
      call logger%info("   state     spin       hartree           eV    "// &
                       "|X|^2-|Y|^2   dominant amplitudes")

      do k = 1, size(excitations)
         ! Unconditional, unlike the paired-only column this replaced: every
         ! route is scaled to the same convention now, so the number means
         ! the same thing on all three and is worth printing on all three.
         weight = dot_product(x(:, k), x(:, k)) - dot_product(y(:, k), y(:, k))
         write (line, "(a,i4,a9,f16.9,f13.4,f15.9,a)") "   ", k, &
            trim(spin_word(state_spin(k))), excitations(k), &
            excitations(k)*HARTREE_TO_EV, weight, "   "
         contributions: do i = 1, n_occ
            do a = 1, n_vir
               idx = (i - 1)*n_vir + a
               if (abs(x(idx, k)) < AMPLITUDE_FLOOR) cycle
               write (piece, "(i0,a,i0,a,f7.3,a)") i, " -> ", n_occ + a, " (", &
                  x(idx, k), ")  "
               ! Out of line: stop, rather than skip this one and keep going.
               ! A later, shorter contribution would fit and be printed where
               ! the omitted ones should have been, and the row would read as
               ! a complete list in orbital order when it is not one. The
               ! ellipsis says the row was cut, where there is room to say it.
               if (len_trim(line) + len_trim(piece) + 1 > len(line)) then
                  if (len_trim(line) + 4 <= len(line)) line = trim(line)//" ..."
                  exit contributions
               end if
               line = trim(line)//" "//trim(piece)
            end do
         end do contributions
         call logger%info(trim(line))
      end do
   end subroutine log_state_table

   pure function spin_word(code) result(word)
      !! The `STATE_SPIN_*` code as the word the table prints
      integer, intent(in) :: code
      character(len=8) :: word

      select case (code)
      case (STATE_SPIN_SINGLET)
         word = "singlet"
      case (STATE_SPIN_TRIPLET)
         word = "triplet"
      case default
         word = "unknown"
      end select
   end function spin_word

   ! ---------------------------------------------------------------------------
   ! The unrestricted response operator
   !
   ! Everything below is the same three routes -- Tamm-Dancoff, paired, and the
   ! dense matrices the tests compare -- over a reference that has two sets of
   ! orbitals instead of one. It is a second core rather than a flag on the
   ! first because almost every array doubles: two occupied blocks, two virtual
   ! blocks, two gap matrices, two reference densities and a kernel with three
   ! spin channels rather than one. A restricted solve would then carry the
   ! beta halves of all of it as empty arrays it never reads.
   !
   ! What is shared is what does not double: the solvers, `roots_to_solve`,
   ! `davidson_subspace`, `unit_matrix` and the reporting.
   !
   ! **No singlet or triplet.** An unrestricted reference is not a spin
   ! eigenfunction, so neither are the roots of its response operator, and the
   ! spin-resolved kernel has no combination that would separate them. The
   ! restricted manifolds are the two combinations `z_beta = +/- z_alpha` of
   ! this operator, which is what the cross-check test asserts, and there is no
   ! `spin` argument here to ask for one of them.
   ! ---------------------------------------------------------------------------

   subroutine core_uhf_destroy(this)
      !! Release everything the core holds and leave it as newly declared
      !!
      !! **The kernel cache is why this exists.** It is a dozen grid-sized
      !! arrays -- hundreds of megabytes on a large molecule -- and an
      !! operator takes a whole copy of the core, so between the copy and the
      !! solve there are two of them alive. Destroying the builder's copy as
      !! soon as the operator has taken it leaves one.
      !!
      !! The two pointers are nulled rather than deallocated: their targets
      !! are the caller's molecule and exchange-correlation context, which
      !! outlive this object by construction.
      class(response_core_uhf_t), intent(inout) :: this

      nullify (this%mol)
      nullify (this%xc)
      if (allocated(this%c_occ_a)) deallocate (this%c_occ_a)
      if (allocated(this%c_vir_a)) deallocate (this%c_vir_a)
      if (allocated(this%c_occ_b)) deallocate (this%c_occ_b)
      if (allocated(this%c_vir_b)) deallocate (this%c_vir_b)
      if (allocated(this%gaps_a)) deallocate (this%gaps_a)
      if (allocated(this%gaps_b)) deallocate (this%gaps_b)
      if (allocated(this%zero_h)) deallocate (this%zero_h)
      if (allocated(this%bounds)) deallocate (this%bounds)
      if (allocated(this%ref_a)) deallocate (this%ref_a)
      if (allocated(this%ref_b)) deallocate (this%ref_b)
      call this%cache%destroy()
      this%k_scale = 1.0_dp
      this%rs_k_lr = 0.0_dp
      this%rs_omega = 0.0_dp
      this%n_occ_a = 0
      this%n_vir_a = 0
      this%n_occ_b = 0
      this%n_vir_b = 0
      this%batch = DEFAULT_RESPONSE_BATCH
      this%n_products = 0
   end subroutine core_uhf_destroy

   pure function core_uhf_length(this) result(n)
      !! How long a trial vector is: both spin blocks end to end
      class(response_core_uhf_t), intent(in) :: this
      integer :: n

      n = this%n_occ_a*this%n_vir_a + this%n_occ_b*this%n_vir_b
   end function core_uhf_length

   function core_uhf_diagonal(this) result(diag)
      !! `e_a - e_i` of both spins, flat: what the solvers precondition on
      class(response_core_uhf_t), intent(in) :: this
      real(dp), allocatable :: diag(:)

      integer :: na

      na = this%n_occ_a*this%n_vir_a
      allocate (diag(this%length()))
      diag(1:na) = reshape(this%gaps_a, [na])
      diag(na + 1:) = reshape(this%gaps_b, [this%n_occ_b*this%n_vir_b])
   end function core_uhf_diagonal

   pure function core_uhf_has_exchange(this) result(yes)
      !! Whether `(A-B)` is anything but the orbital-energy diagonal
      class(response_core_uhf_t), intent(in) :: this
      logical :: yes

      yes = this%k_scale /= 0.0_dp .or. this%rs_k_lr /= 0.0_dp
   end function core_uhf_has_exchange

   subroutine core_uhf_half_block(this, vectors, images, minus, error)
      !! `(A+B)v` or `(A-B)v` for a block of flat spin-concatenated vectors
      !!
      !! The unpacking is the whole of it: a trial vector is the alpha
      !! rectangle followed by the beta one, each virtual fastest, and
      !! `response_product_uhf` wants them as two three-dimensional arrays.
      class(response_core_uhf_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)   !! (n_ov, n_vectors)
      real(dp), intent(out) :: images(:, :)   !! (n_ov, n_vectors)
      logical, intent(in) :: minus            !! `(A-B)` rather than `(A+B)`
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: ua(:, :, :), ub(:, :, :), aua(:, :, :), aub(:, :, :)
      real(dp), allocatable :: gaps(:)
      integer :: n_ov, n_vec, width, first, last, w, m, na, nb

      if (error%has_error()) return

      n_ov = this%length()
      na = this%n_occ_a*this%n_vir_a
      nb = this%n_occ_b*this%n_vir_b
      n_vec = size(vectors, 2)
      if (n_ov < 1 .or. .not. associated(this%mol)) then
         call error%set(ERROR_VALIDATION, "the unrestricted response operator was "// &
                        "applied before it was given a reference to apply itself over")
         images = 0.0_dp
         return
      end if
      if (size(vectors, 1) /= n_ov .or. size(images, 1) /= n_ov) then
         call error%set(ERROR_VALIDATION, "a trial vector handed to the unrestricted "// &
                        "response operator is not the length of the two "// &
                        "occupied-virtual blocks together")
         images = 0.0_dp
         return
      end if

      ! No integrals at all when the difference is the diagonal: `(A-B)` is
      ! exchange only and exchange is same-spin, so with none of it the two
      ! blocks are their own gaps and nothing couples them.
      if (minus .and. .not. this%has_exchange()) then
         gaps = this%diagonal()
         do m = 1, n_vec
            images(:, m) = gaps*vectors(:, m)
         end do
         this%n_products = this%n_products + n_vec
         return
      end if

      width = min(max(this%batch, 1), n_vec)
      allocate (ua(this%n_vir_a, this%n_occ_a, width), aua(this%n_vir_a, this%n_occ_a, width))
      allocate (ub(this%n_vir_b, this%n_occ_b, width), aub(this%n_vir_b, this%n_occ_b, width))

      do first = 1, n_vec, width
         last = min(first + width - 1, n_vec)
         w = last - first + 1
         do m = 1, w
            ua(:, :, m) = reshape(vectors(1:na, first + m - 1), [this%n_vir_a, this%n_occ_a])
            ub(:, :, m) = reshape(vectors(na + 1:, first + m - 1), [this%n_vir_b, this%n_occ_b])
         end do

         aua = 0.0_dp
         aub = 0.0_dp
         if (associated(this%xc)) then
            call response_product_uhf(this%mol, this%c_occ_a, this%c_vir_a, &
                                      this%c_occ_b, this%c_vir_b, this%gaps_a, &
                                      this%gaps_b, this%zero_h, this%bounds, &
                                      ua(:, :, 1:w), ub(:, :, 1:w), minus, &
                                      aua(:, :, 1:w), aub(:, :, 1:w), error, &
                                      k_scale=this%k_scale, xc=this%xc, &
                                      ref_a=this%ref_a, ref_b=this%ref_b, &
                                      rs_k_lr=this%rs_k_lr, rs_omega=this%rs_omega, &
                                      cache=this%cache)
         else
            call response_product_uhf(this%mol, this%c_occ_a, this%c_vir_a, &
                                      this%c_occ_b, this%c_vir_b, this%gaps_a, &
                                      this%gaps_b, this%zero_h, this%bounds, &
                                      ua(:, :, 1:w), ub(:, :, 1:w), minus, &
                                      aua(:, :, 1:w), aub(:, :, 1:w), error, &
                                      k_scale=this%k_scale)
         end if
         if (error%has_error()) exit

         do m = 1, w
            images(1:na, first + m - 1) = reshape(aua(:, :, m), [na])
            images(na + 1:, first + m - 1) = reshape(aub(:, :, m), [nb])
         end do
         this%n_products = this%n_products + w
      end do

      if (error%has_error()) images = 0.0_dp
      deallocate (ua, ub, aua, aub)
   end subroutine core_uhf_half_block

   pure function tda_uhf_length(this) result(n)
      !! How long a trial vector is
      class(tda_operator_uhf_t), intent(in) :: this
      integer :: n

      n = this%core%length()
   end function tda_uhf_length

   function tda_uhf_diagonal(this) result(diag)
      !! What the solver preconditions on and starts from
      class(tda_operator_uhf_t), intent(in) :: this
      real(dp), allocatable :: diag(:)

      diag = this%core%diagonal()
   end function tda_uhf_diagonal

   subroutine tda_uhf_apply_many(this, vectors, images, error)
      !! `A x = (1/2)[(A+B)x + (A-B)x]` for a block of trial vectors
      class(tda_operator_uhf_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)
      real(dp), intent(out) :: images(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: minus_image(:, :)

      if (error%has_error()) return

      call this%core%half_block(vectors, images, .false., error)
      if (error%has_error()) return

      allocate (minus_image(size(images, 1), size(images, 2)))
      call this%core%half_block(vectors, minus_image, .true., error)
      if (error%has_error()) then
         images = 0.0_dp
      else
         images = 0.5_dp*(images + minus_image)
      end if
      deallocate (minus_image)
   end subroutine tda_uhf_apply_many

   subroutine tda_uhf_apply(this, vector, image, error)
      !! `A x` for one vector, as a block of one
      class(tda_operator_uhf_t), intent(inout) :: this
      real(dp), intent(in) :: vector(:)
      real(dp), intent(out) :: image(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: block_in(:, :), block_out(:, :)

      allocate (block_in(size(vector), 1), block_out(size(image), 1))
      block_in(:, 1) = vector
      call this%apply_many(block_in, block_out, error)
      image = block_out(:, 1)
      deallocate (block_in, block_out)
   end subroutine tda_uhf_apply

   subroutine rpa_uhf_apply_plus(this, vectors, images, error)
      !! `(A+B)v` for a block of trial vectors
      class(rpa_operator_uhf_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)
      real(dp), intent(out) :: images(:, :)
      type(error_t), intent(inout) :: error

      call this%core%half_block(vectors, images, .false., error)
   end subroutine rpa_uhf_apply_plus

   subroutine rpa_uhf_apply_minus(this, vectors, images, error)
      !! `(A-B)v` for a block of trial vectors
      class(rpa_operator_uhf_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)
      real(dp), intent(out) :: images(:, :)
      type(error_t), intent(inout) :: error

      call this%core%half_block(vectors, images, .true., error)
   end subroutine rpa_uhf_apply_minus

   subroutine build_response_core_uhf(mol, orbitals_a, energies_a, n_occ_a, &
                                      orbitals_b, energies_b, n_occ_b, core, error, &
                                      xc, ref_a, ref_b, bounds, batch)
      !! Point a `response_core_uhf_t` at a converged unrestricted reference
      !!
      !! `mol` and `xc` are `target` and the core keeps pointers to them, so
      !! both have to outlive every product taken through it. `xc` and the two
      !! reference densities are one argument in three halves and are refused
      !! separately: a kernel is evaluated at a density, and a context arriving
      !! without one has nothing to evaluate.
      !!
      !! The kernel cache is filled here rather than on first use, for the
      !! reason the restricted core fills its own: it is a property of the
      !! converged densities, and a response solve applies the kernel hundreds
      !! of times on a quadrature that is most of a Kohn-Sham run.
      !!
      !! **A spin with nothing in it is allowed.** A high-spin reference can
      !! have no beta electrons at all -- triplet H2, quartet lithium -- and
      !! its alpha excitations are as well defined as any other doublet's.
      !! That spin simply contributes no rotations, so the trial vector is the
      !! alpha block alone and the beta half of every product is empty. What
      !! is refused is a reference with no single excitation in *either* spin,
      !! which is the case there is nothing to solve for.
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals_a(:, :)   !! (n_ao, n_mo)
      real(dp), intent(in) :: energies_a(:)      !! (n_mo), ascending
      integer, intent(in) :: n_occ_a
      real(dp), intent(in) :: orbitals_b(:, :), energies_b(:)
      integer, intent(in) :: n_occ_b
      type(response_core_uhf_t), intent(out) :: core
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional, target :: xc
         !! Spin-polarised, which the kernel checks.
      real(dp), intent(in), optional :: ref_a(:, :), ref_b(:, :)
         !! The converged spin densities, `C_sigma C_sigma^T` and not doubled.
      real(dp), intent(in), optional :: bounds(:, :)
         !! Schwarz bounds, computed here when the caller has none.
      integer, intent(in), optional :: batch
         !! Trial vectors per integral pass; `DEFAULT_RESPONSE_BATCH` absent.

      integer :: n_ao, n_mo

      if (error%has_error()) return

      n_ao = size(orbitals_a, 1)
      n_mo = size(orbitals_a, 2)

      if (n_occ_a < 0 .or. n_occ_b < 0 .or. n_occ_a > n_mo .or. n_occ_b > n_mo) then
         call error%set(ERROR_VALIDATION, "an unrestricted reference was handed to "// &
                        "the response operator with an occupation outside the "// &
                        "orbital space it spans")
         return
      end if
      if (n_occ_a*(n_mo - n_occ_a) + n_occ_b*(n_mo - n_occ_b) < 1) then
         call error%set(ERROR_VALIDATION, "an unrestricted excitation needs one "// &
                        "occupied and one virtual orbital in at least one spin; "// &
                        "this reference has nothing to excite between in either")
         return
      end if
      if (size(orbitals_b, 2) /= n_mo .or. size(orbitals_b, 1) /= n_ao) then
         call error%set(ERROR_VALIDATION, "the two spins of this reference span "// &
                        "different orbital spaces, so there is no common basis to "// &
                        "build a response operator over")
         return
      end if
      if (size(energies_a) < n_mo .or. size(energies_b) < n_mo) then
         call error%set(ERROR_VALIDATION, "there are fewer orbital energies than "// &
                        "orbitals, so the response diagonal cannot be formed")
         return
      end if
      if (present(xc) .neqv. (present(ref_a) .and. present(ref_b))) then
         call error%set(ERROR_VALIDATION, "the unrestricted response operator was "// &
                        "given an exchange-correlation context without both spin "// &
                        "densities its kernel is evaluated at, or the reverse; it "// &
                        "needs all three or none")
         return
      end if

      core%mol => mol
      core%n_occ_a = n_occ_a
      core%n_vir_a = n_mo - n_occ_a
      core%n_occ_b = n_occ_b
      core%n_vir_b = n_mo - n_occ_b
      core%c_occ_a = orbitals_a(:, 1:n_occ_a)
      core%c_vir_a = orbitals_a(:, n_occ_a + 1:n_mo)
      core%c_occ_b = orbitals_b(:, 1:n_occ_b)
      core%c_vir_b = orbitals_b(:, n_occ_b + 1:n_mo)
      allocate (core%zero_h(n_ao, n_ao))
      core%zero_h = 0.0_dp
      call fill_gaps(energies_a, n_occ_a, n_mo, core%gaps_a)
      call fill_gaps(energies_b, n_occ_b, n_mo, core%gaps_b)
      if (present(batch)) then
         if (batch > 0) core%batch = batch
      end if

      if (present(bounds)) then
         core%bounds = bounds
      else
         call schwarz_bounds(mol, core%bounds, error)
         if (error%has_error()) return
      end if

      if (present(xc)) then
         core%xc => xc
         core%ref_a = ref_a
         core%ref_b = ref_b
         core%k_scale = xc%exx_fraction
         if (xc%range_separated) then
            core%rs_k_lr = xc%rs_k_lr
            core%rs_omega = xc%rs_omega
         end if
         call xc_kernel_cache_uks_fill(xc, mol, ref_a, ref_b, core%cache, error)
         if (error%has_error()) return
      else
         core%k_scale = 1.0_dp
      end if

      if (any(core%gaps_a <= 0.0_dp) .or. any(core%gaps_b <= 0.0_dp)) then
         call error%set(ERROR_VALIDATION, "an occupied orbital of one spin lies "// &
                        "above a virtual one of the same spin, so these are not "// &
                        "the aufbau orbitals and the gaps the excitation solver "// &
                        "preconditions on are not positive")
         return
      end if
   end subroutine build_response_core_uhf

   subroutine fill_gaps(energies, n_occ, n_mo, gaps)
      !! `e_a - e_i` as an `(n_vir, n_occ)` rectangle
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ, n_mo
      real(dp), allocatable, intent(out) :: gaps(:, :)

      integer :: i, a

      allocate (gaps(n_mo - n_occ, n_occ))
      do i = 1, n_occ
         do a = 1, n_mo - n_occ
            gaps(a, i) = energies(n_occ + a) - energies(i)
         end do
      end do
   end subroutine fill_gaps

   subroutine build_tda_operator_uhf(mol, orbitals_a, energies_a, n_occ_a, &
                                     orbitals_b, energies_b, n_occ_b, operator, &
                                     error, xc, ref_a, ref_b, bounds, batch)
      !! Point a `tda_operator_uhf_t` at a converged unrestricted reference
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals_a(:, :), energies_a(:)
      integer, intent(in) :: n_occ_a
      real(dp), intent(in) :: orbitals_b(:, :), energies_b(:)
      integer, intent(in) :: n_occ_b
      type(tda_operator_uhf_t), intent(out) :: operator
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: ref_a(:, :), ref_b(:, :)
      real(dp), intent(in), optional :: bounds(:, :)
      integer, intent(in), optional :: batch

      if (present(xc)) then
         call build_response_core_uhf(mol, orbitals_a, energies_a, n_occ_a, &
                                      orbitals_b, energies_b, n_occ_b, operator%core, &
                                      error, xc=xc, ref_a=ref_a, ref_b=ref_b, &
                                      bounds=bounds, batch=batch)
      else
         call build_response_core_uhf(mol, orbitals_a, energies_a, n_occ_a, &
                                      orbitals_b, energies_b, n_occ_b, operator%core, &
                                      error, bounds=bounds, batch=batch)
      end if
   end subroutine build_tda_operator_uhf

   subroutine build_rpa_operator_uhf(mol, orbitals_a, energies_a, n_occ_a, &
                                     orbitals_b, energies_b, n_occ_b, operator, &
                                     error, xc, ref_a, ref_b, bounds, batch)
      !! Point an `rpa_operator_uhf_t` at a converged unrestricted reference
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals_a(:, :), energies_a(:)
      integer, intent(in) :: n_occ_a
      real(dp), intent(in) :: orbitals_b(:, :), energies_b(:)
      integer, intent(in) :: n_occ_b
      type(rpa_operator_uhf_t), intent(out) :: operator
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: ref_a(:, :), ref_b(:, :)
      real(dp), intent(in), optional :: bounds(:, :)
      integer, intent(in), optional :: batch

      if (present(xc)) then
         call build_response_core_uhf(mol, orbitals_a, energies_a, n_occ_a, &
                                      orbitals_b, energies_b, n_occ_b, operator%core, &
                                      error, xc=xc, ref_a=ref_a, ref_b=ref_b, &
                                      bounds=bounds, batch=batch)
      else
         call build_response_core_uhf(mol, orbitals_a, energies_a, n_occ_a, &
                                      orbitals_b, energies_b, n_occ_b, operator%core, &
                                      error, bounds=bounds, batch=batch)
      end if
   end subroutine build_rpa_operator_uhf

   subroutine tda_dense_matrix_uhf(operator, a, error)
      !! The explicit unrestricted `A`, from the operator applied to unit vectors
      !!
      !! Both spin blocks at once, so what comes back is the whole
      !! `(n_ov_a + n_ov_b)` square with the two diagonal blocks and the two
      !! coupling ones in it. For a small case and for checking the solver
      !! against a dense diagonalisation; nothing on the production path calls
      !! it.
      type(tda_operator_uhf_t), intent(inout) :: operator
      real(dp), allocatable, intent(out) :: a(:, :)
         !! (n_ov, n_ov), column `j` being `A e_j`
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: unit_vectors(:, :)
      integer :: n_ov

      if (error%has_error()) return

      n_ov = operator%core%length()
      call unit_matrix(n_ov, unit_vectors)
      allocate (a(n_ov, n_ov))
      call operator%apply_many(unit_vectors, a, error)
      deallocate (unit_vectors)
      if (error%has_error()) deallocate (a)
   end subroutine tda_dense_matrix_uhf

   subroutine rpa_dense_matrices_uhf(operator, aplus, aminus, error)
      !! The explicit unrestricted `(A+B)` and `(A-B)`, from the shipped operator
      type(rpa_operator_uhf_t), intent(inout) :: operator
      real(dp), allocatable, intent(out) :: aplus(:, :), aminus(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: unit_vectors(:, :)
      integer :: n_ov

      if (error%has_error()) return

      n_ov = operator%core%length()
      call unit_matrix(n_ov, unit_vectors)
      allocate (aplus(n_ov, n_ov), aminus(n_ov, n_ov))
      call operator%apply_plus(unit_vectors, aplus, error)
      if (.not. error%has_error()) then
         call operator%apply_minus(unit_vectors, aminus, error)
      end if
      deallocate (unit_vectors)
      if (error%has_error()) deallocate (aplus, aminus)
   end subroutine rpa_dense_matrices_uhf

   subroutine response_excitations_uhf(mol, orbitals_a, energies_a, n_occ_a, &
                                       orbitals_b, energies_b, n_occ_b, n_states, &
                                       method, excitations, state_spin, x_amplitudes, &
                                       y_amplitudes, error, xc, ref_a, ref_b, bounds, &
                                       tolerance, max_iter, max_subspace, batch, &
                                       verbose)
      !! The lowest excitation energies of an unrestricted reference
      !!
      !! What comes back is ascending, in Hartree above the reference, with
      !! every root below `EXCITATION_FLOOR` already dropped -- so `size` of it
      !! can be smaller than `n_states`, and a caller has to read the size.
      !!
      !! **One normalisation, both routes.** `sum_sigma (|X|^2 - |Y|^2) = 1`,
      !! which is what the paired solver returns and what a unit Tamm-Dancoff
      !! vector already satisfies. The restricted routes disagree with each
      !! other on this; here there is no closed-shell factor of two to put
      !! anywhere, so there is nothing for the two to disagree about.
      !!
      !! **Every root is labelled `STATE_SPIN_UNRESTRICTED`.** The reference is
      !! not a spin eigenfunction, so its excitations are not singlets or
      !! triplets, and there is no `spin` argument asking for one.
      !!
      !! **A root near zero is not always an artefact here.** A doublet's
      !! response operator has a rotation of the half-filled shell in it whose
      !! Tamm-Dancoff root is small but not zero -- 6.7e-3 hartree on the OH
      !! radical -- and the paired problem drops it as an `omega^2` at the
      !! numerical zero. Neither is a fault in the solver: Tamm-Dancoff keeps
      !! it because `A` alone is not the operator whose null space the rotation
      !! lives in. The floor below removes what is genuinely at zero and no
      !! more.
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals_a(:, :), energies_a(:)
      integer, intent(in) :: n_occ_a
      real(dp), intent(in) :: orbitals_b(:, :), energies_b(:)
      integer, intent(in) :: n_occ_b
      integer, intent(in) :: n_states
      character(len=*), intent(in) :: method   !! `tda` or `rpa`
      real(dp), allocatable, intent(out) :: excitations(:)
      integer, allocatable, intent(out) :: state_spin(:)
      real(dp), allocatable, intent(out) :: x_amplitudes(:, :)
         !! (n_ov_a + n_ov_b, n_found), alpha block then beta, virtual fastest
      real(dp), allocatable, intent(out) :: y_amplitudes(:, :)
         !! The same shape; zero for `tda`, which has no `Y`
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: ref_a(:, :), ref_b(:, :)
      real(dp), intent(in), optional :: bounds(:, :)
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iter
      integer, intent(in), optional :: max_subspace
      integer, intent(in), optional :: batch
      logical, intent(in), optional :: verbose

      type(response_core_uhf_t) :: core
      type(tda_operator_uhf_t) :: tda
      type(rpa_operator_uhf_t) :: rpa
      real(dp), allocatable :: diagonal(:), raw(:), vectors(:, :), residuals(:)
      real(dp), allocatable :: xpy(:, :), xmy(:, :), all_x(:, :), all_y(:, :)
      real(dp), allocatable :: guess(:, :)
      character(len=MAX_LINE_LENGTH) :: line
      character(len=16) :: route
      real(dp) :: tol
      integer :: n_ov, n_solve, n_guess, cap, subspace, iterations, products
      integer :: n_found, k, keep
      integer :: dropped, n_products, occ_a, vir_a, occ_b, vir_b
      logical :: converged

      if (error%has_error()) return
      ! Allocated empty rather than left unallocated: a caller reads the size,
      ! and one doing that on an unallocated array has no way to notice.
      if (n_states < 1) then
         allocate (excitations(0), state_spin(0), x_amplitudes(0, 0))
         allocate (y_amplitudes(0, 0))
         return
      end if

      route = trim(adjustl(method))
      if (route /= "tda" .and. route /= "rpa") then
         call error%set(ERROR_VALIDATION, "'"//trim(route)//"' is not an unrestricted "// &
                        "linear-response method this backend knows; it has 'tda' and "// &
                        "'rpa'. The 'casida' reduction is restricted-only, since it "// &
                        "assumes (A-B) is the orbital-energy diagonal")
         return
      end if

      if (present(xc)) then
         call build_response_core_uhf(mol, orbitals_a, energies_a, n_occ_a, &
                                      orbitals_b, energies_b, n_occ_b, core, error, &
                                      xc=xc, ref_a=ref_a, ref_b=ref_b, bounds=bounds, &
                                      batch=batch)
      else
         call build_response_core_uhf(mol, orbitals_a, energies_a, n_occ_a, &
                                      orbitals_b, energies_b, n_occ_b, core, error, &
                                      bounds=bounds, batch=batch)
      end if
      if (error%has_error()) return

      diagonal = core%diagonal()
      n_ov = size(diagonal)
      ! Read off before the core is released below; the state table needs them
      ! and nothing else does.
      occ_a = core%n_occ_a
      vir_a = core%n_vir_a
      occ_b = core%n_occ_b
      vir_b = core%n_vir_b
      if (n_states > n_ov) then
         call error%set(ERROR_VALIDATION, "keywords.excited_states asked for "// &
                        to_char(n_states)//" roots, but this reference has only "// &
                        to_char(n_ov)//" single excitations to build them from")
         return
      end if

      tol = DEFAULT_EXCITED_TOL
      if (present(tolerance)) tol = tolerance
      subspace = 0
      if (present(max_subspace)) then
         if (max_subspace > 0) subspace = max_subspace
      end if

      n_solve = roots_to_solve(diagonal, n_states)
      n_guess = guess_count(diagonal, n_solve)
      cap = davidson_subspace(subspace, n_solve, n_guess, n_ov)
      n_guess = min(n_guess, cap)
      allocate (all_x(n_ov, n_solve), all_y(n_ov, n_solve))

      ! The operator takes a whole copy of the core, kernel cache included, so
      ! the builder's copy is released here rather than at the end of the
      ! routine: one grid-sized cache is alive across the solve, not two.
      if (trim(route) == "rpa") then
         rpa%core = core
         call core%destroy()
         call rpa_solve(rpa, diagonal, n_solve, raw, xpy, xmy, residuals, iterations, &
                        products, converged, error, tolerance=tol, &
                        max_iterations=max_iter, max_subspace=subspace, &
                        verbose=verbose, label="unrestricted RPA iterations")
         n_products = rpa%core%n_products
         call rpa%core%destroy()
         if (error%has_error()) return
         ! `xpy . xmy = 1` out of the solver, which **is** the convention here:
         ! an unrestricted amplitude carries one spin orbital, not two, so
         ! there is no closed-shell half to apply.
         do k = 1, n_solve
            all_x(:, k) = 0.5_dp*(xpy(:, k) + xmy(:, k))
            all_y(:, k) = 0.5_dp*(xpy(:, k) - xmy(:, k))
         end do
      else
         tda%core = core
         call core%destroy()
         ! `GUESS_PER_ROOT` vectors a root rather than one, for the reason
         ! that constant gives: the guess is picked on the spin-orbital gaps
         ! and what comes back is the spectrum, and the coupling reorders the
         ! two. The paired route below needs nothing -- `rpa_solve` builds a
         ! start of the same width itself.
         allocate (guess(n_ov, n_guess))
         call fill_guess(diagonal, n_guess, guess)
         call davidson_flat(tda, diagonal, n_solve, raw, vectors, residuals, &
                            iterations, products, converged, error, tolerance=tol, &
                            max_iterations=max_iter, max_subspace=cap, guess=guess, &
                            verbose=verbose, label="unrestricted TDA iterations", &
                            value_label="excitation")
         deallocate (guess)
         n_products = tda%core%n_products
         call tda%core%destroy()
         if (error%has_error()) return
         all_x = vectors
         all_y = 0.0_dp
      end if

      if (.not. converged) then
         call error%set(ERROR_GENERIC, "the unrestricted "//trim(route)//" solve did "// &
                        "not converge its roots in the iterations allowed; raise "// &
                        "keywords.excited_states.max_iter or loosen "// &
                        "keywords.excited_states.tolerance")
         return
      end if

      ! Drop the artefacts, then the degeneracy padding, keeping the order.
      keep = 0
      do k = 1, n_solve
         if (raw(k) > EXCITATION_FLOOR) keep = keep + 1
      end do
      n_found = min(keep, n_states)
      allocate (excitations(n_found), state_spin(n_found))
      allocate (x_amplitudes(n_ov, n_found), y_amplitudes(n_ov, n_found))
      state_spin = STATE_SPIN_UNRESTRICTED
      keep = 0
      do k = 1, n_solve
         ! The negated `>`, not `<=`, so the two passes ask the same question
         ! of a root that is neither: the rotation of an open shell is an
         ! `omega^2` at the numerical zero of either sign, and a square root
         ! of it can come back as a NaN, which fails both comparisons.
         if (.not. (raw(k) > EXCITATION_FLOOR)) cycle
         keep = keep + 1
         if (keep > n_found) exit
         excitations(keep) = raw(k)
         x_amplitudes(:, keep) = all_x(:, k)
         y_amplitudes(:, keep) = all_y(:, k)
      end do

      if (n_found < n_states) then
         call logger%warning("  only "//to_char(n_found)//" of the "// &
                             to_char(n_states)//" roots asked for are excitations; "// &
                             "the rest converged below "//to_char(EXCITATION_FLOOR)// &
                             " hartree and are rotations of the reference, not "// &
                             "excited states")
      end if

      ! Say when a route dropped something, because the two routes drop
      ! different things and the spectrum alone does not show it. A doublet's
      ! rotation of its own open shell is an `omega^2` on the numerical zero,
      ! which the paired route loses here; Tamm-Dancoff keeps the same
      ! rotation as a small positive root and reports it as state 1. Without
      ! this line the two spectra look like a disagreement.
      dropped = 0
      do k = 1, n_solve
         if (.not. (raw(k) > EXCITATION_FLOOR)) dropped = dropped + 1
      end do
      if (dropped > 0) then
         call logger%info("  "//to_char(dropped)//" of the "//to_char(n_solve)// &
                          " converged root(s) are not above "// &
                          to_char(EXCITATION_FLOOR)//" hartree and are not "// &
                          "reported: on an open shell these are rotations of the "// &
                          "reference rather than excitations. Tamm-Dancoff keeps "// &
                          "such a rotation as a small root instead of dropping it, "// &
                          "so the two routes need not report the same first state.")
      end if

      write (line, "(a,a,a,i0,a,i0,a)") "  ", route_name(route), ": ", &
         size(excitations), " root(s) from ", n_products, &
         " matrix-vector products"
      call logger%info(trim(line))
      call log_state_table_uhf(route, excitations, x_amplitudes, y_amplitudes, &
                               occ_a, vir_a, occ_b, vir_b)
   end subroutine response_excitations_uhf

   subroutine log_state_table_uhf(route, excitations, x, y, n_occ_a, n_vir_a, &
                                  n_occ_b, n_vir_b)
      !! The unrestricted spectrum, with the spin orbitals each root is made of
      !!
      !! As the restricted table, with the dominant amplitudes tagged `a` or
      !! `b`: the two spin blocks have different orbital numbering, and a bare
      !! `4 -> 6` would name two different excitations depending on which half
      !! of the vector it came from.
      !!
      !! A root below `ROTATION_HINT` is marked with a star and the reason
      !! printed under the table: on an open shell that is where a rotation
      !! of the singly-occupied orbitals turns up, and it is the one place
      !! the Tamm-Dancoff and paired spectra of the same reference disagree
      !! about which state is first.
      character(len=*), intent(in) :: route
      real(dp), intent(in) :: excitations(:)
      real(dp), intent(in) :: x(:, :), y(:, :)
      integer, intent(in) :: n_occ_a, n_vir_a, n_occ_b, n_vir_b

      real(dp), parameter :: AMPLITUDE_FLOOR = 0.1_dp
      character(len=MAX_LINE_LENGTH) :: line, piece
      real(dp) :: weight
      integer :: k, i, a, idx, na
      logical :: paired, flagged

      if (size(excitations) < 1) return
      flagged = .false.
      paired = trim(route) /= "tda"
      na = n_occ_a*n_vir_a

      call logger%info("  unrestricted amplitudes are normalised to "// &
                       "sum_spin (|X|^2 - |Y|^2) = 1")
      if (paired) then
         call logger%info("   state       hartree           eV    "// &
                          "|X|^2-|Y|^2   dominant amplitudes")
      else
         call logger%info("   state       hartree           eV   dominant amplitudes")
      end if

      do k = 1, size(excitations)
         if (paired) then
            weight = dot_product(x(:, k), x(:, k)) - dot_product(y(:, k), y(:, k))
            write (line, "(a,i4,f16.9,f13.4,f15.9,a)") "   ", k, excitations(k), &
               excitations(k)*HARTREE_TO_EV, weight, "   "
         else
            write (line, "(a,i4,f16.9,f13.4,a)") "   ", k, excitations(k), &
               excitations(k)*HARTREE_TO_EV, "   "
         end if
         do i = 1, n_occ_a
            do a = 1, n_vir_a
               idx = (i - 1)*n_vir_a + a
               if (abs(x(idx, k)) < AMPLITUDE_FLOOR) cycle
               write (piece, "(i0,a,i0,a,f7.3,a)") i, "a -> ", n_occ_a + a, "a (", &
                  x(idx, k), ")  "
               if (len_trim(line) + len_trim(piece) + 1 > len(line)) cycle
               line = trim(line)//" "//trim(piece)
            end do
         end do
         do i = 1, n_occ_b
            do a = 1, n_vir_b
               idx = na + (i - 1)*n_vir_b + a
               if (abs(x(idx, k)) < AMPLITUDE_FLOOR) cycle
               write (piece, "(i0,a,i0,a,f7.3,a)") i, "b -> ", n_occ_b + a, "b (", &
                  x(idx, k), ")  "
               if (len_trim(line) + len_trim(piece) + 1 > len(line)) cycle
               line = trim(line)//" "//trim(piece)
            end do
         end do
         if (excitations(k) < ROTATION_HINT) then
            flagged = .true.
            if (len_trim(line) + 2 <= len(line)) line = trim(line)//" *"
         end if
         call logger%info(trim(line))
      end do

      if (flagged) then
         call logger%info("  * below "//to_char(ROTATION_HINT)//" hartree: on an "// &
                          "open shell a root this low is commonly a rotation of the "// &
                          "singly-occupied orbitals rather than an excitation. The "// &
                          "paired route puts such a rotation at the numerical zero "// &
                          "and drops it; Tamm-Dancoff reports it.")
      end if
   end subroutine log_state_table_uhf

end module mqc_czt_tddft
