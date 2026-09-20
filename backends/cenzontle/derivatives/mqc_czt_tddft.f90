!! Linear-response excitation energies for a closed-shell reference
module mqc_czt_tddft
   !! The singlet Tamm-Dancoff and random-phase eigenproblems, over the
   !! response operator this backend already applies.
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
   !! **Amplitudes come back at different normalisations for the two
   !! approximations, and `singlet_excitations` says which.** Tamm-Dancoff has
   !! no `Y`, so its `X` is a unit vector; the paired problem conserves
   !! `|X|^2 - |Y|^2` and nothing else, and its amplitudes carry the
   !! restricted closed-shell convention `|X|^2 - |Y|^2 = 1/2` that PySCF,
   !! Psi4 and the transition moments of Layer 5 are written in.
   !!
   !! ## What is not here
   !!
   !! Triplets, oscillator strengths and an unrestricted reference.
   !! `excited_decline_reason` in `mqc_czt_bridge` refuses what cannot be
   !! computed; the bridge refuses the rest by name rather than answering a
   !! different question.

   ! TODO(mqc): the two amplitude normalisations above should be one.
   ! Tamm-Dancoff amplitudes are unit vectors because that is what Layer 2
   ! shipped, and halving them to match the paired route would silently move
   ! every printed amplitude of an existing calculation -- so until it is
   ! done, Layer 5 has to know which route produced what it is handed.
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
   implicit none
   private

   public :: response_core_t
   public :: tda_operator_t
   public :: rpa_operator_t
   public :: build_tda_operator
   public :: build_rpa_operator
   public :: tda_dense_matrix
   public :: rpa_dense_matrices
   public :: singlet_excitations

   real(dp), parameter :: EXCITATION_FLOOR = 1.0e-3_dp
      !! Roots below this are not reported.
      !!
      !! A closed-shell singlet spectrum has nothing down here: the lowest
      !! excitation of even a small-gap molecule is tenths of a Hartree. What
      !! does appear is an artefact -- a rotation the reference is marginally
      !! unstable along, or a root the solver has not separated from zero --
      !! and neither is an excitation. The same threshold separates the
      !! degenerate block the guess is extended over.

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
      !! `|X|^2 - |Y|^2` the paired amplitudes are handed back at.
      !!
      !! One half, not one. A closed-shell excitation is two spin-orbital
      !! excitations of equal weight and the spatial-orbital amplitude carries
      !! both, so the norm the spin-orbital problem sets to one comes to a half
      !! here. It is the convention PySCF and Psi4 report, and the one the
      !! transition moment `mu = 2 sum <i|r|a> (X+Y)` is written for; handing
      !! back unit-normalised amplitudes instead would put a factor of the
      !! square root of two into every oscillator strength downstream.

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
                                  rs_omega=this%rs_omega, cache=this%cache)
         else
            call response_product(this%mol, this%c_occ, this%c_vir, this%gaps, &
                                  this%zero_h, u, idx, w, minus, au, error, &
                                  direct=.true., bounds=this%bounds, &
                                  k_scale=this%k_scale)
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
                                  xc, reference, bounds, batch)
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

      integer :: n_ao, n_mo, n_vir, i, a

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
         call xc_kernel_cache_fill(xc, mol, reference, core%cache, error)
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
                                 xc, reference, bounds, batch)
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

      if (present(xc)) then
         call build_response_core(mol, orbitals, energies, n_occ, operator%core, &
                                  error, xc=xc, reference=reference, bounds=bounds, &
                                  batch=batch)
      else
         call build_response_core(mol, orbitals, energies, n_occ, operator%core, &
                                  error, bounds=bounds, batch=batch)
      end if
   end subroutine build_tda_operator

   subroutine build_rpa_operator(mol, orbitals, energies, n_occ, operator, error, &
                                 xc, reference, bounds, batch)
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

      if (present(xc)) then
         call build_response_core(mol, orbitals, energies, n_occ, operator%core, &
                                  error, xc=xc, reference=reference, bounds=bounds, &
                                  batch=batch)
      else
         call build_response_core(mol, orbitals, energies, n_occ, operator%core, &
                                  error, bounds=bounds, batch=batch)
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

   subroutine singlet_excitations(mol, orbitals, energies, n_occ, n_states, method, &
                                  excitations, x_amplitudes, y_amplitudes, error, xc, &
                                  reference, bounds, tolerance, max_iter, &
                                  max_subspace, batch, verbose)
      !! The lowest singlet excitation energies of a closed shell
      !!
      !! What comes back is ascending, in Hartree above the reference, with
      !! every root below `EXCITATION_FLOOR` already dropped -- so `size` of
      !! it can be smaller than `n_states`, and a caller has to read the size
      !! rather than assume it. Always allocated when this returns without an
      !! error, empty included.
      !!
      !! **Normalisation depends on `method`,** which is the one thing about
      !! this routine a caller cannot ignore. `tda` returns a unit `X` and a
      !! `Y` of zeros -- there is no `Y` in that approximation, and the zeros
      !! say so rather than standing for something left unfilled. `rpa` and
      !! `casida` return `X` and `Y` at `sum(X^2) - sum(Y^2) = 1/2`, the
      !! restricted closed-shell convention.
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      integer, intent(in) :: n_states
      character(len=*), intent(in) :: method
         !! `tda`, `rpa`, or `casida` for the cross-check reduction, which
         !! needs a functional carrying no exact exchange.
      real(dp), allocatable, intent(out) :: excitations(:)
         !! (n_found) excitation energies, ascending, in Hartree
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

      type(tda_operator_t) :: tda
      type(rpa_operator_t) :: rpa
      type(casida_operator_t) :: casida
      real(dp), allocatable :: diagonal(:), values(:), vectors(:, :), residuals(:)
      real(dp), allocatable :: xpy(:, :), xmy(:, :), all_x(:, :), all_y(:, :)
      real(dp) :: tol
      character(len=MAX_LINE_LENGTH) :: line
      character(len=16) :: route
      integer :: n_ov, n_solve, iterations, products, n_found, k, keep, subspace
      logical :: converged

      if (error%has_error()) return
      ! Allocated empty rather than left unallocated: the contract above is
      ! that a caller reads the size, and a caller doing that on an
      ! unallocated array has no way to notice.
      if (n_states < 1) then
         allocate (excitations(0), x_amplitudes(0, 0), y_amplitudes(0, 0))
         return
      end if

      route = trim(adjustl(method))
      if (route /= "tda" .and. route /= "rpa" .and. route /= "casida") then
         call error%set(ERROR_VALIDATION, "'"//trim(route)//"' is not a linear-"// &
                        "response method this backend knows; it has 'tda', 'rpa' "// &
                        "and the 'casida' cross-check")
         return
      end if

      select case (route)
      case ("tda")
         if (present(xc)) then
            call build_tda_operator(mol, orbitals, energies, n_occ, tda, error, &
                                    xc=xc, reference=reference, bounds=bounds, &
                                    batch=batch)
         else
            call build_tda_operator(mol, orbitals, energies, n_occ, tda, error, &
                                    bounds=bounds, batch=batch)
         end if
         if (error%has_error()) return
         diagonal = tda%core%diagonal()
      case ("rpa")
         if (present(xc)) then
            call build_rpa_operator(mol, orbitals, energies, n_occ, rpa, error, &
                                    xc=xc, reference=reference, bounds=bounds, &
                                    batch=batch)
         else
            call build_rpa_operator(mol, orbitals, energies, n_occ, rpa, error, &
                                    bounds=bounds, batch=batch)
         end if
         if (error%has_error()) return
         diagonal = rpa%core%diagonal()
      case default
         if (present(xc)) then
            call build_response_core(mol, orbitals, energies, n_occ, casida%core, &
                                     error, xc=xc, reference=reference, &
                                     bounds=bounds, batch=batch)
         else
            call build_response_core(mol, orbitals, energies, n_occ, casida%core, &
                                     error, bounds=bounds, batch=batch)
         end if
         if (error%has_error()) return
         if (casida%core%has_exchange()) then
            call error%set(ERROR_VALIDATION, "the Casida reduction assumes (A-B) is "// &
                           "the orbital-energy diagonal, which holds only for a "// &
                           "functional carrying no exact exchange; this reference "// &
                           "keeps a fraction of it, so ask for 'rpa'")
            return
         end if
         diagonal = casida%core%diagonal()
         casida%root_gaps = sqrt(diagonal)
      end select

      n_ov = size(diagonal)
      if (n_states > n_ov) then
         call error%set(ERROR_VALIDATION, "keywords.excited_states asked for "// &
                        to_char(n_states)//" roots, but this reference has only "// &
                        to_char(n_ov)//" single excitations to build them from")
         return
      end if

      n_solve = roots_to_solve(diagonal, n_states)
      tol = DEFAULT_EXCITED_TOL
      if (present(tolerance)) tol = tolerance
      subspace = 0
      if (present(max_subspace)) then
         if (max_subspace > 0) subspace = max_subspace
      end if

      allocate (all_x(n_ov, n_solve), all_y(n_ov, n_solve))

      select case (route)
      case ("rpa")
         call rpa_solve(rpa, diagonal, n_solve, values, xpy, xmy, residuals, &
                        iterations, products, converged, error, tolerance=tol, &
                        max_iterations=max_iter, max_subspace=subspace, &
                        verbose=verbose, label="RPA iterations")
         if (error%has_error()) return
         ! `xpy . xmy = 1` out of the solver; the closed-shell convention is a
         ! half of that, and both vectors take the same factor so their half
         ! sum and half difference are `X` and `Y`.
         do k = 1, n_solve
            all_x(:, k) = sqrt(RHF_PAIRED_NORM)*0.5_dp*(xpy(:, k) + xmy(:, k))
            all_y(:, k) = sqrt(RHF_PAIRED_NORM)*0.5_dp*(xpy(:, k) - xmy(:, k))
         end do
         products = rpa%core%n_products
      case ("casida")
         ! The eigenvalue is `w^2` and so is the preconditioner: `dEps^2` is
         ! the diagonal of the reduced operator up to its two-electron part,
         ! the way `dEps` is the diagonal of `A`.
         call davidson_flat(casida, diagonal*diagonal, n_solve, values, vectors, &
                            residuals, iterations, products, converged, error, &
                            tolerance=tol, max_iterations=max_iter, &
                            max_subspace=davidson_subspace(subspace, n_solve, n_ov), &
                            verbose=verbose, label="Casida iterations", &
                            value_label="omega^2")
         if (error%has_error()) return
         call casida_amplitudes(values, vectors, diagonal, all_x, all_y)
         values = sqrt(max(values, 0.0_dp))
         products = casida%core%n_products
      case default
         call davidson_flat(tda, diagonal, n_solve, values, vectors, residuals, &
                            iterations, products, converged, error, tolerance=tol, &
                            max_iterations=max_iter, &
                            max_subspace=davidson_subspace(subspace, n_solve, n_ov), &
                            verbose=verbose, label="TDA iterations", &
                            value_label="excitation")
         if (error%has_error()) return
         all_x = vectors
         all_y = 0.0_dp
         products = tda%core%n_products
      end select

      if (.not. converged) then
         call error%set(ERROR_GENERIC, "the "//trim(route)//" solve did not converge "// &
                        "its roots in the iterations allowed; raise "// &
                        "keywords.excited_states.max_iter or loosen "// &
                        "keywords.excited_states.tolerance")
         return
      end if

      ! Drop the artefacts, then the degeneracy padding, keeping the order.
      keep = 0
      do k = 1, n_solve
         if (values(k) > EXCITATION_FLOOR) keep = keep + 1
      end do
      n_found = min(keep, n_states)
      allocate (excitations(n_found), x_amplitudes(n_ov, n_found))
      allocate (y_amplitudes(n_ov, n_found))
      keep = 0
      do k = 1, n_solve
         if (values(k) <= EXCITATION_FLOOR) cycle
         keep = keep + 1
         if (keep > n_found) exit
         excitations(keep) = values(k)
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

      write (line, "(a,a,a,i0,a,i0,a)") "  ", route_name(route), ": ", n_found, &
         " singlet root(s) from ", products, " matrix-vector products"
      call logger%info(trim(line))
      call log_state_table(route, excitations, x_amplitudes, y_amplitudes, n_occ, &
                           n_ov/n_occ)
   end subroutine singlet_excitations

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

   subroutine log_state_table(route, excitations, x, y, n_occ, n_vir)
      !! The spectrum, as a table with the orbitals each root is made of
      !!
      !! Amplitudes below `AMPLITUDE_FLOOR` are left out: a converged root of
      !! a molecule of any size has a handful of contributions above it and
      !! hundreds of numerical dust below, and printing the dust hides the
      !! assignment the table exists for. Orbitals are numbered from one over
      !! all of them, occupied and virtual together, so the indices match what
      !! every other table in this program prints.
      !!
      !! The paired routes carry a `|X|^2-|Y|^2` column. It is the norm the
      !! problem conserves, and it is printed rather than asserted: a root
      !! whose value has drifted off the convention is one whose
      !! biorthonormalisation did not take, and every transition moment built
      !! from it afterwards would be wrong by that factor.
      character(len=*), intent(in) :: route
      real(dp), intent(in) :: excitations(:)
      real(dp), intent(in) :: x(:, :), y(:, :)
      integer, intent(in) :: n_occ, n_vir

      real(dp), parameter :: AMPLITUDE_FLOOR = 0.1_dp
      character(len=MAX_LINE_LENGTH) :: line, piece
      real(dp) :: weight
      integer :: k, i, a, idx
      logical :: paired

      if (size(excitations) < 1) return
      paired = trim(route) /= "tda"

      ! Said once, beside the numbers it applies to: the two codes a
      ! cross-check is run against convert with a different Hartree.
      if (abs(HARTREE_TO_EV - PYSCF_HARTREE_TO_EV) > 0.0_dp) then
         write (line, "(a,f16.12,a,f14.8,a)") "  excitation energies in eV use ", &
            HARTREE_TO_EV, " eV/hartree; PySCF uses ", PYSCF_HARTREE_TO_EV, &
            ", which differs in the eighth decimal"
         call logger%info(trim(line))
      end if
      if (paired) then
         call logger%info("  "//route_name(route)//" amplitudes are normalised to "// &
                          "|X|^2 - |Y|^2 = 0.5")
         call logger%info("   state       hartree           eV    |X|^2-|Y|^2   "// &
                          "dominant amplitudes")
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

end module mqc_czt_tddft
