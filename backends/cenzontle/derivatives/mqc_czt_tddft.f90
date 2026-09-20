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
   use mqc_calculation_defaults, only: DEFAULT_RESPONSE_BATCH, DEFAULT_EXCITED_TOL
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_xc, only: xc_context_t, xc_kernel_cache_t, xc_kernel_cache_fill
   use mqc_czt_response_product, only: response_product
   use mqc_davidson, only: davidson_flat, sigma_operator_t
   use mqc_czt_rpa_solver, only: paired_operator_t, rpa_solve
   use mqc_result_types, only: STATE_SPIN_SINGLET, STATE_SPIN_TRIPLET
   implicit none
   private

   public :: response_core_t
   public :: tda_operator_t
   public :: rpa_operator_t
   public :: build_tda_operator
   public :: build_rpa_operator
   public :: tda_dense_matrix
   public :: rpa_dense_matrices
   public :: response_excitations

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

   real(dp), parameter :: DEGENERACY_WINDOW = 1.0e-3_dp
      !! How close two orbital-energy gaps have to be for the guess to have to
      !! carry both.
      !!
      !! The starting vectors are unit vectors on the lowest gaps. Splitting a
      !! degenerate pair between the guess and the space outside it leaves the
      !! solver converging one partner against a subspace that cannot represent
      !! the other, which stalls rather than converging to the wrong answer --
      !! so the count is extended over the whole degenerate block and the extra
      !! roots are found and then not reported.

   integer, parameter :: MAX_EXTRA_ROOTS = 8
      !! A cap on that extension. A highly symmetric molecule can put many gaps
      !! inside the window, and every one of them costs a converged root; eight
      !! covers a degeneracy no point group produces and stops a pathological
      !! case from turning a five-root request into a fifty-root solve.

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
      type(response_core_t) :: core
   contains
      procedure :: apply => tda_apply
      procedure :: apply_many => tda_apply_many
      procedure :: length => tda_length
      procedure :: diagonal => tda_diagonal
   end type tda_operator_t

   type, extends(paired_operator_t) :: rpa_operator_t
      !! `(A+B)` and `(A-B)` as something `mqc_czt_rpa_solver` will pair up
      type(response_core_t) :: core
   contains
      procedure :: apply_plus => rpa_apply_plus
      procedure :: apply_minus => rpa_apply_minus
   end type rpa_operator_t

   type, extends(sigma_operator_t) :: casida_operator_t
      !! `dEps^{1/2}(A+B)dEps^{1/2}` for a pure functional, as a Davidson operator
      !!
      !! Not public: this is the cross-check route, reached through
      !! `singlet_excitations` with `method = "casida"` and refused for
      !! anything carrying exact exchange, where `(A-B)` is not the diagonal
      !! it assumes.
      type(response_core_t) :: core
      real(dp), allocatable :: root_gaps(:)
         !! (n_ov) `sqrt(e_a - e_i)`, flat
   contains
      procedure :: apply => casida_apply
      procedure :: apply_many => casida_apply_many
   end type casida_operator_t

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

      if (error%has_error()) return

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

      integer :: n_ao, n_mo, n_vir, i, a
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
      allocate (core%gaps(n_vir, n_occ))
      do i = 1, n_occ
         do a = 1, n_vir
            core%gaps(a, i) = energies(n_occ + a) - energies(i)
         end do
      end do
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

      if (any(core%gaps <= 0.0_dp)) then
         call error%set(ERROR_VALIDATION, "an occupied orbital lies above a virtual "// &
                        "one, so these are not the aufbau orbitals and the gaps the "// &
                        "excitation solver preconditions on are not positive")
         return
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

      real(dp), allocatable :: sorted(:)
      real(dp) :: best
      logical, allocatable :: taken(:)
      integer :: n, k, i, pick

      n = size(diagonal)
      n_solve = min(max(n_states, 1), n)
      allocate (taken(n), sorted(n_solve))
      taken = .false.
      ! A partial selection rather than a sort: `n_states` is small and the
      ! space is not.
      do k = 1, n_solve
         pick = 0
         best = 0.0_dp
         do i = 1, n
            if (taken(i)) cycle
            if (pick == 0 .or. diagonal(i) < best) then
               pick = i
               best = diagonal(i)
            end if
         end do
         taken(pick) = .true.
         sorted(k) = diagonal(pick)
      end do

      n_solve = count(diagonal <= sorted(n_solve) + DEGENERACY_WINDOW)
      n_solve = min(n_solve, min(n_states, n) + MAX_EXTRA_ROOTS, n)
      n_solve = max(n_solve, min(max(n_states, 1), n))
      deallocate (taken, sorted)
   end function roots_to_solve

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
      type(response_core_t), intent(inout) :: core
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
      character(len=:), allocatable :: manifold
      integer :: n_ov, n_solve, iterations, products, n_found, k, keep, imaginary
      logical :: converged

      if (error%has_error()) return

      core%triplet = is_triplet
      manifold = trim(manifold_word(is_triplet))
      diagonal = core%diagonal()
      n_ov = size(diagonal)
      n_solve = roots_to_solve(diagonal, n_states)
      allocate (all_x(n_ov, n_solve), all_y(n_ov, n_solve))

      select case (trim(route))
      case ("rpa")
         rpa%core = core
         imaginary = 0
         call rpa_solve(rpa, diagonal, n_solve, raw, xpy, xmy, residuals, &
                        iterations, products, converged, error, tolerance=tol, &
                        max_iterations=max_iter, max_subspace=subspace, &
                        verbose=verbose, label=manifold//" RPA iterations", &
                        imaginary_roots=imaginary)
         core%n_products = rpa%core%n_products
         if (error%has_error()) then
            call name_the_instability(is_triplet, error)
            return
         end if
         ! The paired solver skips an imaginary root rather than returning
         ! it, and only complains when it runs out of real ones. Asking for
         ! fewer roots than there are imaginary ones therefore succeeds and
         ! says nothing -- which for a triplet manifold would be the
         ! instability going unreported under a clean spectrum.
         if (is_triplet .and. imaginary > 0) then
            call unstable_reference(0.0_dp, error, imaginary)
            return
         end if
         ! `xpy . xmy = 1` out of the solver; the closed-shell convention is a
         ! half of that, and both vectors take the same factor so their half
         ! sum and half difference are `X` and `Y`.
         do k = 1, n_solve
            all_x(:, k) = sqrt(RHF_PAIRED_NORM)*0.5_dp*(xpy(:, k) + xmy(:, k))
            all_y(:, k) = sqrt(RHF_PAIRED_NORM)*0.5_dp*(xpy(:, k) - xmy(:, k))
         end do
      case ("casida")
         casida%core = core
         if (casida%core%has_exchange()) then
            call error%set(ERROR_VALIDATION, "the Casida reduction assumes (A-B) is "// &
                           "the orbital-energy diagonal, which holds only for a "// &
                           "functional carrying no exact exchange; this reference "// &
                           "keeps a fraction of it, so ask for 'rpa'")
            return
         end if
         casida%root_gaps = sqrt(diagonal)
         ! The eigenvalue is `w^2` and so is the preconditioner: `dEps^2` is
         ! the diagonal of the reduced operator up to its two-electron part,
         ! the way `dEps` is the diagonal of `A`.
         call davidson_flat(casida, diagonal*diagonal, n_solve, raw, vectors, &
                            residuals, iterations, products, converged, error, &
                            tolerance=tol, max_iterations=max_iter, &
                            max_subspace=davidson_subspace(subspace, n_solve, n_ov), &
                            verbose=verbose, label=manifold//" Casida iterations", &
                            value_label="omega^2")
         core%n_products = casida%core%n_products
         if (error%has_error()) return
         call casida_amplitudes(raw, vectors, diagonal, all_x, all_y)
         ! A negative `w^2` is the instability arriving as an imaginary
         ! frequency rather than as a failed factorisation, so it is clamped
         ! to zero here and caught by the floor below.
         raw = sqrt(max(raw, 0.0_dp))
      case default
         tda%core = core
         call davidson_flat(tda, diagonal, n_solve, raw, vectors, residuals, &
                            iterations, products, converged, error, tolerance=tol, &
                            max_iterations=max_iter, &
                            max_subspace=davidson_subspace(subspace, n_solve, n_ov), &
                            verbose=verbose, label=manifold//" TDA iterations", &
                            value_label="excitation")
         core%n_products = tda%core%n_products
         if (error%has_error()) return
         ! The Davidson returns unit eigenvectors; the closed-shell
         ! convention is `|X|^2 = 1/2`, and scaling here rather than at every
         ! reader is what makes the three routes one convention.
         all_x = sqrt(RHF_PAIRED_NORM)*vectors
         all_y = 0.0_dp
      end select

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

   subroutine unstable_reference(lowest, error, imaginary)
      !! Report a non-positive triplet root as what it is
      !!
      !! Concatenated rather than written into a buffer: the message is longer
      !! than an output record, and a format-directed write that overflows one
      !! is a run-time failure rather than a truncation.
      real(dp), intent(in) :: lowest   !! The offending root, in Hartree
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: imaginary
         !! Squared frequencies that came out negative, when that is how the
         !! instability showed itself rather than as a root at the floor.

      character(len=:), allocatable :: what

      if (present(imaginary)) then
         what = to_char(imaginary)//" of its squared triplet frequencies are negative"
      else
         what = "its lowest triplet root is "//to_char(lowest)//" hartree, at or "// &
                "below zero"
      end if
      call error%set(ERROR_GENERIC, "the reference is triplet-unstable: "//what// &
                     ", so this closed shell is a saddle point against spin "// &
                     "polarisation and a lower unrestricted solution exists. What "// &
                     "the solver would report is an expansion about that saddle "// &
                     "point rather than an excitation spectrum; converge an "// &
                     "unrestricted reference.")
   end subroutine unstable_reference

   subroutine name_the_instability(is_triplet, error)
      !! Say what a failed paired solve means when the manifold is a triplet
      !!
      !! The Stratmann-Scuseria-Frisch reduction factorises the projected
      !! `(A-B)`, and an instability reaches it as a factorisation that will
      !! not go through. That is the right diagnosis in the solver's own
      !! terms and the wrong one for a reader, who has asked for a spectrum
      !! and wants to know that the *reference* is what is wrong. The solver's
      !! message is kept and prefixed rather than replaced.
      !!
      !! Only a failure the solver itself called an instability is relabelled.
      !! The paired solve fails for other reasons -- a LAPACK error, a
      !! subspace too small for the roots asked for -- and calling one of
      !! those a triplet instability would be a diagnosis invented from the
      !! manifold rather than read off the arithmetic.
      logical, intent(in) :: is_triplet
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: was

      if (.not. is_triplet) return
      if (.not. error%has_error()) return
      was = error%get_message()
      if (index(was, "instability") == 0 .and. index(was, "unstable") == 0) return
      call error%set(ERROR_GENERIC, "the reference is triplet-unstable: the paired "// &
                     "solve of the triplet manifold could not be reduced, which for "// &
                     "a closed shell means it is a saddle point against spin "// &
                     "polarisation. The solver reported: "//was)
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
         !! absent reproduces each solver's own rule, which is why this is
         !! resolved here rather than forwarded: a deck's unset `max_subspace`
         !! is a zero, and a solver would read that as a subspace of no
         !! vectors.
      integer, intent(in), optional :: batch
      logical, intent(in), optional :: verbose
         !! A line per iteration. Each one is an integral pass.

      type(response_core_t) :: core
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

   function davidson_subspace(requested, n_solve, n_ov) result(subspace)
      !! The Davidson's subspace cap, with a deck's unset zero resolved
      !!
      !! `davidson_flat` takes the cap as a plain integer and would read a
      !! zero as a subspace of no vectors, so its own default rule is restated
      !! here rather than forwarded. Shared by the Tamm-Dancoff and Casida
      !! routes, which both go through that solver; the paired solver takes
      !! the zero itself and has a rule of its own.
      integer, intent(in) :: requested   !! Zero for the default
      integer, intent(in) :: n_solve, n_ov
      integer :: subspace

      subspace = max(2*n_solve + 8, 16)
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
         weight = dot_product(x(:, k), x(:, k)) - dot_product(y(:, k), y(:, k))
         write (line, "(a,i4,a9,f16.9,f13.4,f15.9,a)") "   ", k, &
            trim(spin_word(state_spin(k))), excitations(k), &
            excitations(k)*HARTREE_TO_EV, weight, "   "
         do i = 1, n_occ
            do a = 1, n_vir
               idx = (i - 1)*n_vir + a
               if (abs(x(idx, k)) < AMPLITUDE_FLOOR) cycle
               write (piece, "(i0,a,i0,a,f7.3,a)") i, " -> ", n_occ + a, " (", &
                  x(idx, k), ")  "
               if (len_trim(line) + len_trim(piece) + 1 > len(line)) cycle
               line = trim(line)//" "//trim(piece)
            end do
         end do
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

end module mqc_czt_tddft
