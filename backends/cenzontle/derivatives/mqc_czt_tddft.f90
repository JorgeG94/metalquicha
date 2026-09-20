!! Tamm-Dancoff excitation energies for a closed-shell reference
module mqc_czt_tddft
   !! The singlet TDA eigenproblem, over the response operator this backend
   !! already applies.
   !!
   !! Linear response asks for the lowest eigenvalues of
   !!
   !!     A_{ia,jb} = d_{ij} d_{ab} (e_a - e_i) + 2(ai|bj) - c_x (ab|ij) + 2 f_xc
   !!
   !! and nothing here computes any of those integrals: `response_product` does,
   !! and `mqc_davidson` finds the eigenvalues. What this module is, is the
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
   !! So `A u` is **two** calls, averaged.
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
   !! ## Layout
   !!
   !! Trial vectors are flat, `idx = (i-1)*n_vir + a` -- virtual fastest, which
   !! is `reshape` of the `(n_vir, n_occ)` rectangle and the layout
   !! `response_product`, `mqc_czt_ov_hessian` and the CI Davidson already
   !! share. `Y` does not exist in the Tamm-Dancoff approximation; the full
   !! paired problem is Layer 4.
   !!
   !! ## What is not here
   !!
   !! Triplets, the full RPA, oscillator strengths and an unrestricted
   !! reference. `excited_decline_reason` in `mqc_czt_bridge` refuses what
   !! cannot be computed; the bridge refuses the rest by name rather than
   !! answering a different question.
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
   implicit none
   private

   public :: tda_operator_t
   public :: build_tda_operator
   public :: tda_dense_matrix
   public :: tda_singlet_excitations

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

   type, extends(sigma_operator_t) :: tda_operator_t
      !! `A` as something `mqc_davidson` will multiply a vector by
      !!
      !! Holds what a response product needs and nothing else. The molecule
      !! and the exchange-correlation context are pointers because both
      !! outlive the solve and neither is cheap to copy; their targets have to
      !! outlive this object.
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
         !! `build_tda_operator`. Unfilled for Hartree-Fock and ignored then.
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
         !! Matrix-vector products spent, for the cost line.
   contains
      procedure :: apply => tda_apply
      procedure :: apply_many => tda_apply_many
      procedure :: length => tda_length
      procedure :: diagonal => tda_diagonal
   end type tda_operator_t

contains

   pure function tda_length(this) result(n)
      !! How long a trial vector is: `n_occ * n_vir`
      class(tda_operator_t), intent(in) :: this
      integer :: n

      n = this%n_occ*this%n_vir
   end function tda_length

   function tda_diagonal(this) result(diag)
      !! What the solver preconditions on and starts from
      !!
      !! `e_a - e_i`, not `diag(A)`. The exact diagonal costs one Fock build
      !! per element, which is the whole matrix; the gaps are what PySCF
      !! preconditions a TDA with, and a preconditioner moves the iteration
      !! count rather than the eigenvalue.
      class(tda_operator_t), intent(in) :: this
      real(dp), allocatable :: diag(:)

      diag = reshape(this%gaps, [this%n_occ*this%n_vir])
   end function tda_diagonal

   subroutine tda_half(this, u, idx, nact, minus, au, error)
      !! One of `(A+B)u` and `(A-B)u`, with the optional arguments resolved
      !!
      !! The branch exists because `xc` and `reference` are present together
      !! or not at all, and neither a null pointer nor an unallocated array
      !! may be passed to a non-optional-shaped dummy. Everything else about
      !! the two calls is the same.
      class(tda_operator_t), intent(inout) :: this
      real(dp), intent(in) :: u(:, :, :)
      integer, intent(in) :: idx(:)
      integer, intent(in) :: nact
      logical, intent(in) :: minus
      real(dp), intent(inout) :: au(:, :, :)
      type(error_t), intent(inout) :: error

      if (error%has_error()) return

      if (associated(this%xc)) then
         call response_product(this%mol, this%c_occ, this%c_vir, this%gaps, &
                               this%zero_h, u, idx, nact, minus, au, error, &
                               direct=.true., bounds=this%bounds, &
                               k_scale=this%k_scale, xc=this%xc, &
                               reference=this%reference, rs_k_lr=this%rs_k_lr, &
                               rs_omega=this%rs_omega, cache=this%cache)
      else
         call response_product(this%mol, this%c_occ, this%c_vir, this%gaps, &
                               this%zero_h, u, idx, nact, minus, au, error, &
                               direct=.true., bounds=this%bounds, &
                               k_scale=this%k_scale)
      end if
   end subroutine tda_half

   subroutine tda_apply_many(this, vectors, images, error)
      !! `A x` for a block of trial vectors, in as few integral passes as it takes
      !!
      !! The block is what makes this worth overriding: one pass over the
      !! quartets and one walk of the quadrature serve every vector in it, so
      !! a Davidson iteration costs what a single product would if it were
      !! alone. `batch` caps the width because `build_fock_direct_many` is
      !! memory-bandwidth bound past a dozen or two densities.
      class(tda_operator_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)   !! (n_ov, n_vectors)
      real(dp), intent(out) :: images(:, :)   !! (n_ov, n_vectors)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: u(:, :, :), ap(:, :, :), am(:, :, :)
      integer, allocatable :: idx(:)
      integer :: n_ov, n_vec, width, first, last, w, m
      logical :: needs_minus

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
         call error%set(ERROR_VALIDATION, "the TDA operator was applied before it "// &
                        "was given a reference to apply itself over")
         images = 0.0_dp
         return
      end if
      if (size(vectors, 1) /= n_ov .or. size(images, 1) /= n_ov) then
         call error%set(ERROR_VALIDATION, "a trial vector handed to the TDA operator "// &
                        "is not the length of the occupied-virtual space")
         images = 0.0_dp
         return
      end if

      ! `(A-B)` is exchange only. With no exchange of any range it is the
      ! orbital-energy diagonal, and building it would be an integral pass
      ! whose every quartet is multiplied by zero.
      needs_minus = this%k_scale /= 0.0_dp .or. this%rs_k_lr /= 0.0_dp

      width = min(max(this%batch, 1), n_vec)
      allocate (u(this%n_vir, this%n_occ, width), ap(this%n_vir, this%n_occ, width))
      allocate (idx(width))
      if (needs_minus) allocate (am(this%n_vir, this%n_occ, width))

      do first = 1, n_vec, width
         last = min(first + width - 1, n_vec)
         w = last - first + 1
         do m = 1, w
            u(:, :, m) = reshape(vectors(:, first + m - 1), [this%n_vir, this%n_occ])
            idx(m) = m
         end do

         ap = 0.0_dp
         call tda_half(this, u, idx, w, .false., ap, error)
         if (error%has_error()) exit
         if (needs_minus) then
            am = 0.0_dp
            call tda_half(this, u, idx, w, .true., am, error)
            if (error%has_error()) exit
         end if

         do m = 1, w
            if (needs_minus) then
               images(:, first + m - 1) = 0.5_dp*(reshape(ap(:, :, m), [n_ov]) &
                                                  + reshape(am(:, :, m), [n_ov]))
            else
               images(:, first + m - 1) = 0.5_dp*(reshape(ap(:, :, m), [n_ov]) &
                                                  + reshape(this%gaps*u(:, :, m), [n_ov]))
            end if
         end do
         this%n_products = this%n_products + w
      end do

      if (error%has_error()) images = 0.0_dp
      deallocate (u, ap, idx)
      if (allocated(am)) deallocate (am)
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

   subroutine build_tda_operator(mol, orbitals, energies, n_occ, operator, error, &
                                 xc, reference, bounds, batch)
      !! Point a `tda_operator_t` at a converged closed-shell reference
      !!
      !! `mol` and `xc` are `target` and the operator keeps pointers to them,
      !! so both have to outlive every product taken through it. `xc` and
      !! `reference` are one argument in two halves and are refused
      !! separately, for the reason `build_scf_ov_hessian` refuses them.
      !!
      !! The kernel cache is filled here rather than on first use: it is a
      !! property of the converged density and a Davidson applies the kernel
      !! hundreds of times, on a quadrature that is most of a DFT run.
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals(:, :)     !! (n_ao, n_mo)
      real(dp), intent(in) :: energies(:)        !! (n_mo), ascending
      integer, intent(in) :: n_occ
      type(tda_operator_t), intent(out) :: operator
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
                        "orbitals, so the TDA diagonal cannot be formed")
         return
      end if
      if (present(xc) .neqv. present(reference)) then
         call error%set(ERROR_VALIDATION, "the TDA operator was given an "// &
                        "exchange-correlation context without the reference density "// &
                        "its kernel is evaluated at, or the reverse; it needs both "// &
                        "or neither")
         return
      end if

      operator%mol => mol
      operator%n_occ = n_occ
      operator%n_vir = n_vir
      operator%c_occ = orbitals(:, 1:n_occ)
      operator%c_vir = orbitals(:, n_occ + 1:n_mo)
      allocate (operator%zero_h(n_ao, n_ao))
      operator%zero_h = 0.0_dp
      allocate (operator%gaps(n_vir, n_occ))
      do i = 1, n_occ
         do a = 1, n_vir
            operator%gaps(a, i) = energies(n_occ + a) - energies(i)
         end do
      end do

      ! Here, and not after the Schwarz bounds and the kernel fill: it reads
      ! nothing but the gaps just formed, and the two it used to sit behind
      ! are the expensive part of building the operator.
      if (any(operator%gaps <= 0.0_dp)) then
         call error%set(ERROR_VALIDATION, "an occupied orbital lies above a virtual "// &
                        "one, so these are not the aufbau orbitals and the gaps the "// &
                        "excitation solver preconditions on are not positive")
         return
      end if

      if (present(batch)) then
         if (batch > 0) operator%batch = batch
      end if

      if (present(bounds)) then
         operator%bounds = bounds
      else
         call schwarz_bounds(mol, operator%bounds, error)
         if (error%has_error()) return
      end if

      if (present(xc)) then
         operator%xc => xc
         operator%reference = reference
         operator%k_scale = xc%exx_fraction
         if (xc%range_separated) then
            operator%rs_k_lr = xc%rs_k_lr
            operator%rs_omega = xc%rs_omega
         end if
         call xc_kernel_cache_fill(xc, mol, reference, operator%cache, error)
         if (error%has_error()) return
      else
         operator%k_scale = 1.0_dp
      end if
   end subroutine build_tda_operator

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
      integer :: n_ov, j

      if (error%has_error()) return

      n_ov = operator%length()
      allocate (unit_vectors(n_ov, n_ov), a(n_ov, n_ov))
      unit_vectors = 0.0_dp
      do j = 1, n_ov
         unit_vectors(j, j) = 1.0_dp
      end do

      call operator%apply_many(unit_vectors, a, error)
      deallocate (unit_vectors)
      if (error%has_error()) deallocate (a)
   end subroutine tda_dense_matrix

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

   subroutine tda_singlet_excitations(mol, orbitals, energies, n_occ, n_states, &
                                      excitations, amplitudes, error, xc, reference, &
                                      bounds, tolerance, max_iter, max_subspace, &
                                      batch, verbose)
      !! The lowest singlet Tamm-Dancoff excitation energies of a closed shell
      !!
      !! What comes back is ascending, in Hartree above the reference, with
      !! every root below `EXCITATION_FLOOR` already dropped -- so `size` of
      !! it can be smaller than `n_states`, and a caller has to read the size
      !! rather than assume it. Always allocated when this returns without an
      !! error, empty included.
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      integer, intent(in) :: n_states
      real(dp), allocatable, intent(out) :: excitations(:)
         !! (n_found) excitation energies, ascending, in Hartree
      real(dp), allocatable, intent(out) :: amplitudes(:, :)
         !! (n_occ*n_vir, n_found) normalised `X`, virtual fastest. `Y` does
         !! not exist in this approximation.
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: reference(:, :)
      real(dp), intent(in), optional :: bounds(:, :)
      real(dp), intent(in), optional :: tolerance
         !! Residual norm a root is accepted at; `DEFAULT_EXCITED_TOL` absent.
      integer, intent(in), optional :: max_iter
      integer, intent(in), optional :: max_subspace
         !! Trial vectors kept before the subspace is collapsed. Zero or
         !! absent reproduces the solver's own rule, which is stated here
         !! rather than forwarded because a deck's unset `max_subspace` is a
         !! zero and the solver would read that as a subspace of no vectors.
      integer, intent(in), optional :: batch
      logical, intent(in), optional :: verbose
         !! A line per Davidson iteration. Each one is an integral pass.

      type(tda_operator_t) :: operator
      real(dp), allocatable :: diagonal(:), values(:), vectors(:, :), residuals(:)
      real(dp) :: tol
      character(len=MAX_LINE_LENGTH) :: line
      integer :: n_ov, n_solve, iterations, products, n_found, k, keep, subspace
      logical :: converged

      if (error%has_error()) return
      ! Allocated empty rather than left unallocated: the contract above is
      ! that a caller reads the size, and a caller doing that on an
      ! unallocated array has no way to notice.
      if (n_states < 1) then
         allocate (excitations(0), amplitudes(0, 0))
         return
      end if

      if (present(xc)) then
         call build_tda_operator(mol, orbitals, energies, n_occ, operator, error, &
                                 xc=xc, reference=reference, bounds=bounds, batch=batch)
      else
         call build_tda_operator(mol, orbitals, energies, n_occ, operator, error, &
                                 bounds=bounds, batch=batch)
      end if
      if (error%has_error()) return

      n_ov = operator%length()
      if (n_states > n_ov) then
         call error%set(ERROR_VALIDATION, "keywords.excited_states asked for "// &
                        to_char(n_states)//" roots, but this reference has only "// &
                        to_char(n_ov)//" single excitations to build them from")
         return
      end if

      diagonal = operator%diagonal()
      n_solve = roots_to_solve(diagonal, n_states)
      tol = DEFAULT_EXCITED_TOL
      if (present(tolerance)) tol = tolerance
      subspace = max(2*n_solve + 8, 16)
      if (present(max_subspace)) then
         if (max_subspace > 0) subspace = max_subspace
      end if
      subspace = max(min(subspace, n_ov), n_solve)

      ! No `guess`: the solver with none starts on unit vectors at the lowest
      ! diagonal elements, and the diagonal here is the orbital-energy gaps,
      ! which is exactly the guess this wants. `n_solve` is what carries the
      ! degeneracy awareness.
      call davidson_flat(operator, diagonal, n_solve, values, vectors, residuals, &
                         iterations, products, converged, error, tolerance=tol, &
                         max_iterations=max_iter, max_subspace=subspace, &
                         verbose=verbose, label="TDA iterations", &
                         value_label="excitation")
      if (error%has_error()) return
      if (.not. converged) then
         call error%set(ERROR_GENERIC, "the Tamm-Dancoff solve did not converge its "// &
                        "roots in the iterations allowed; raise "// &
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
      allocate (excitations(n_found), amplitudes(n_ov, n_found))
      keep = 0
      do k = 1, n_solve
         if (values(k) <= EXCITATION_FLOOR) cycle
         keep = keep + 1
         if (keep > n_found) exit
         excitations(keep) = values(k)
         amplitudes(:, keep) = vectors(:, k)
      end do

      if (n_found < n_states) then
         call logger%warning("  only "//to_char(n_found)//" of the "// &
                             to_char(n_states)//" roots asked for are excitations; "// &
                             "the rest converged below "//to_char(EXCITATION_FLOOR)// &
                             " hartree and are rotations of the reference, not "// &
                             "excited states")
      end if

      write (line, "(a,i0,a,i0,a)") "  Tamm-Dancoff: ", n_found, &
         " singlet root(s) from ", operator%n_products, &
         " matrix-vector products"
      call logger%info(trim(line))
      call log_state_table(excitations, amplitudes, operator%n_occ, operator%n_vir)
   end subroutine tda_singlet_excitations

   subroutine log_state_table(excitations, amplitudes, n_occ, n_vir)
      !! The spectrum, as a table with the orbitals each root is made of
      !!
      !! Amplitudes below `AMPLITUDE_FLOOR` are left out: a converged root of
      !! a molecule of any size has a handful of contributions above it and
      !! hundreds of numerical dust below, and printing the dust hides the
      !! assignment the table exists for. Orbitals are numbered from one over
      !! all of them, occupied and virtual together, so the indices match what
      !! every other table in this program prints.
      real(dp), intent(in) :: excitations(:)
      real(dp), intent(in) :: amplitudes(:, :)
      integer, intent(in) :: n_occ, n_vir

      real(dp), parameter :: AMPLITUDE_FLOOR = 0.1_dp
      character(len=MAX_LINE_LENGTH) :: line, piece
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
      call logger%info("   state       hartree           eV   dominant amplitudes")

      do k = 1, size(excitations)
         write (line, "(a,i4,f16.9,f13.4,a)") "   ", k, excitations(k), &
            excitations(k)*HARTREE_TO_EV, "   "
         contributions: do i = 1, n_occ
            do a = 1, n_vir
               idx = (i - 1)*n_vir + a
               if (abs(amplitudes(idx, k)) < AMPLITUDE_FLOOR) cycle
               write (piece, "(i0,a,i0,a,f7.3,a)") i, " -> ", n_occ + a, " (", &
                  amplitudes(idx, k), ")  "
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

end module mqc_czt_tddft
