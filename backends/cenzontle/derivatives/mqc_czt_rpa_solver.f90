!! The paired eigenproblem of linear response, without either matrix
module mqc_czt_rpa_solver
   !! Stratmann-Scuseria-Frisch: the lowest `omega` of the random-phase
   !! approximation, from nothing but `(A+B)b` and `(A-B)b`.
   !!
   !! The problem is
   !!
   !!     [A  B][X]      [1   0][X]
   !!     [B  A][Y] = w  [0  -1][Y]
   !!
   !! which is `2n` and not Hermitian, and which nobody solves in that form.
   !! Written in the sum and difference vectors `R = X+Y` and `L = X-Y` it is
   !! the pair
   !!
   !!     (A+B) R = w L,      (A-B) L = w R
   !!
   !! and eliminating `L` gives `(A-B)^{1/2}(A+B)(A-B)^{1/2} T = w^2 T` with
   !! `R = (A-B)^{1/2} T`, which **is** Hermitian and `n` by `n`. The square
   !! root is the catch: forming it needs `(A-B)` as a matrix, and the whole
   !! point of an iterative solver is not having one.
   !!
   !! The resolution, which is the method, is to take the square root **inside
   !! the subspace**. `H1 = b^T (A+B) b` and `H2 = b^T (A-B) b` are small, so
   !! `H2^{1/2}` is a dense eigendecomposition of something the size of the
   !! trial space, and every root, every eigenvector and every residual is
   !! built from those two small matrices and the stored products. The large
   !! dimension appears only in the products themselves and in the three
   !! `n`-by-`subspace` arrays they are stored in.
   !!
   !! ## What this does that the reference implementations do not
   !!
   !! Psi4's `hamiltonian_solver` and PySCF's `real_eig` are the two written
   !! forms of this method. Against them:
   !!
   !! * **Products for the new vectors only.** Psi4 recomputes `(A+B)b` and
   !!   `(A-B)b` for the *whole* subspace on every iteration, which for an
   !!   operator whose product is a Fock build is the entire cost of the solve
   !!   multiplied by the iteration count. The products here are stored beside
   !!   the basis and carried through the collapse, where they survive exactly:
   !!   a collapsed vector is a linear combination of vectors whose images are
   !!   already in hand.
   !! * **Twice-repeated Gram-Schmidt.** One pass loses orthogonality when the
   !!   new direction is nearly in the subspace already, which near convergence
   !!   it always is. The CI Davidson in this program has done two passes since
   !!   it was written; Psi4's single pass is its known weak point.
   !! * **Stagnation is reported.** A solve that stops improving returns a named
   !!   error carrying the residual it reached, rather than spinning to
   !!   `max_iter` and rather than returning a number with nothing to say how
   !!   good it is.
   !! * **The paired residual is a maximum, not a sum.** A sum of two norms
   !!   converges a root at twice the tolerance it was asked for.
   !!
   !! ## The instability
   !!
   !! `(A-B)` is positive definite exactly when the reference is a minimum of
   !! the energy with respect to the rotations the operator spans. When it is
   !! not, `H2^{1/2}` does not exist, `w^2` is not what this solves for, and
   !! there is no excitation spectrum to report -- the reference itself is
   !! wrong. That is caught here on the projected `H2`, by name, rather than
   !! left to a square root of a negative number downstream.
   !!
   !! ## What comes back
   !!
   !! `R` and `L`, biorthonormalised so that `R . L = |X|^2 - |Y|^2 = 1` for
   !! each root. A caller wanting `X` and `Y` separately takes their half sum
   !! and half difference, and rescales to whatever normalisation convention it
   !! documents -- this module fixes one and says which.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_io, only: to_char
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   use mqc_convergence_report, only: convergence_header, convergence_footer
   use, intrinsic :: iso_fortran_env, only: output_unit, int64
   implicit none
   private

   public :: paired_operator_t
   public :: rpa_solve

   real(dp), parameter, public :: DEFAULT_RPA_TOLERANCE = 1.0e-6_dp
      !! On the larger of the two paired residual norms. Both reference codes
      !! floor near 1e-7 on a quadrature, so this is a default that can be
      !! reached rather than the CI Davidson's 1e-10.
   integer, parameter, public :: DEFAULT_RPA_MAX_ITERATIONS = 100

   integer, parameter :: GUESS_PER_ROOT = 4
      !! Starting vectors per root, on the lowest diagonal elements.
      !!
      !! One per root is what a Hermitian Davidson starts from and it is not
      !! enough here: the paired problem needs both `R` and `L` inside the
      !! subspace, the two differ once `B` is not zero, and a subspace that can
      !! hold only one of them converges the pair against its own projection.
      !! Four is what both reference implementations use.

   real(dp), parameter :: DEGENERACY_WINDOW = 1.0e-3_dp
      !! How close two diagonal elements have to be for the guess to carry
      !! both. Splitting a degenerate block across the edge of the starting
      !! space leaves one partner converging against a subspace that cannot
      !! represent the other.

   integer, parameter :: MAX_EXTRA_GUESS = 16
      !! A cap on that extension, so a highly symmetric system with many
      !! coincident gaps cannot turn a five-root request into a fifty-vector
      !! starting space.

   real(dp), parameter :: OMEGA2_FLOOR = 1.0e-10_dp
      !! Below this, `w^2` is not a root. Zero eigenvalues are the rotations
      !! the reference is flat along and negative ones an instability of
      !! `(A+B)`; neither is an excitation, and `sqrt` of the second is not a
      !! number.

   real(dp), parameter :: PRECONDITION_FLOOR = 1.0e-4_dp
      !! A diagonal element sitting on the current `w` would divide by
      !! nothing. Unlike the CI Davidson, which shifts the denominator to its
      !! own floor, the correction is taken **unscaled** there -- the residual
      !! itself, which is a direction the subspace does not yet hold and is
      !! what the step is for. Scaling it by 1e4 instead would add a vector
      !! whose entries are all one element.

   real(dp), parameter :: LINEAR_DEPENDENCE = 1.0e-8_dp
      !! How much of a normalised correction has to survive projection out of
      !! the subspace for it to be worth adding. Tested on the normalised
      !! vector, so the threshold does not scale with the residual and the
      !! subspace keeps growing as the solve converges.

   integer, parameter :: STAGNATION_WINDOW = 3
      !! Iterations without an improvement in the worst residual before the
      !! solve is declared stuck. Two is too few: a residual plateaus for an
      !! iteration when a new root enters the window and then drops again.

   integer, parameter :: STAGNATION_GRACE = 3
      !! Iterations before the check above is armed at all. The first few
      !! iterations of a paired solve can move the root ordering around, and a
      !! reordering looks like a residual that went up.

   real(dp), parameter :: SUBSPACE_MEMORY_LIMIT = 2.0e9_dp
      !! Bytes the default subspace cap allows the basis and its two product
      !! sets. Three arrays of `n` by `max_subspace`, so the default rule --
      !! which is a multiple of the root count and knows nothing about `n` --
      !! is cut back to fit rather than asking for terabytes on a large case.

   type, abstract :: paired_operator_t
      !! Whatever can apply both halves of the response operator to a block
      !!
      !! `(A+B)` and `(A-B)` and nothing else: no `A`, no `B`, no diagonal
      !! (which arrives as an argument, so an operator can be built before its
      !! preconditioner is known), and nothing that says what the vector
      !! indexes. The excitation solver, a frequency-dependent
      !! coupled-perturbed solve and a dense test matrix are all the same
      !! problem to this module.
      !!
      !! Both procedures take and return a block of columns and must pair them
      !! by position: column `i` of `images` is the operator applied to column
      !! `i` of `vectors`. An implementation is free to compute the block
      !! however it likes, and is expected to -- sharing one pass over the
      !! integrals across the block is the reason the interface is a block at
      !! all.
   contains
      procedure(paired_apply_i), deferred :: apply_plus
      procedure(paired_apply_i), deferred :: apply_minus
   end type paired_operator_t

   abstract interface
      subroutine paired_apply_i(this, vectors, images, error)
         import :: paired_operator_t, dp, error_t
         implicit none
         class(paired_operator_t), intent(inout) :: this
         real(dp), intent(in) :: vectors(:, :)
         real(dp), intent(out) :: images(:, :)
         type(error_t), intent(inout) :: error
      end subroutine paired_apply_i
   end interface

contains

   subroutine rpa_solve(operator, diagonal, n_roots, omega, xpy, xmy, residuals, &
                        iterations_taken, products, converged, error, tolerance, &
                        max_iterations, max_subspace, verbose, label)
      !! The `n_roots` lowest excitation energies of a paired response operator
      !!
      !! `diagonal` is what the corrections are preconditioned with and what
      !! the starting vectors are chosen from -- the orbital-energy differences
      !! for a linear response, which is the diagonal of both `A+B` and `A-B`
      !! up to the two-electron terms. A wrong preconditioner costs iterations
      !! and not accuracy.
      !!
      !! **Normalisation.** `xpy` and `xmy` come back biorthonormal,
      !! `sum(xpy(:,k)*xmy(:,k)) = 1`, which is `|X|^2 - |Y|^2 = 1` for root
      !! `k`. The restricted closed-shell convention that halves that is the
      !! caller's to apply, and to say it has.
      class(paired_operator_t), intent(inout) :: operator
      real(dp), intent(in) :: diagonal(:)
         !! (n) the operator's diagonal, used to precondition and to start
      integer, intent(in) :: n_roots
      real(dp), allocatable, intent(out) :: omega(:)
         !! (n_roots) excitation energies, ascending
      real(dp), allocatable, intent(out) :: xpy(:, :)
         !! (n, n_roots) the right-hand vectors `X+Y`
      real(dp), allocatable, intent(out) :: xmy(:, :)
         !! (n, n_roots) the left-hand vectors `X-Y`
      real(dp), allocatable, intent(out) :: residuals(:)
         !! (n_roots) the larger of the two paired residual norms, per root
      integer, intent(out) :: iterations_taken
      integer, intent(out) :: products
         !! Applications of one of the two halves to one trial vector. A new
         !! trial vector costs one of each, so this counts two per vector --
         !! which is what it costs, an operator that knows `(A-B)` is diagonal
         !! and answers it for free notwithstanding.
      logical, intent(out) :: converged
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      integer, intent(in), optional :: max_subspace
         !! Trial vectors held before the subspace is collapsed onto the
         !! current `R` and `L`. Absent or zero takes the rule below.
      logical, intent(in), optional :: verbose
         !! A line per iteration. Each one is an integral pass, or several.
      character(len=*), intent(in), optional :: label
         !! What the iteration table is called.

      real(dp), allocatable :: basis(:, :), sp(:, :), sm(:, :)
      real(dp), allocatable :: h1(:, :), h2(:, :), root(:, :), hss(:, :)
      real(dp), allocatable :: rss(:, :), lss(:, :), work(:, :)
      real(dp), allocatable :: w2(:), wl(:, :), wr(:, :), correction(:)
      real(dp), allocatable :: kept(:)
      integer(int64) :: tick, last, rate
      character(len=128) :: line
      real(dp) :: tol, best_worst, worst, norm, unstable_tol
      integer :: n, nmax, iterations, iteration, nsub, n_new, first_new
      integer :: k, i, added, stall, n_positive, n_flat
      logical :: loud
      logical, allocatable :: root_converged(:)

      if (error%has_error()) return

      n = size(diagonal)
      iterations_taken = 0
      products = 0
      converged = .false.

      tol = DEFAULT_RPA_TOLERANCE
      if (present(tolerance)) then
         if (tolerance > 0.0_dp) tol = tolerance
      end if
      iterations = DEFAULT_RPA_MAX_ITERATIONS
      if (present(max_iterations)) then
         if (max_iterations > 0) iterations = max_iterations
      end if

      if (n_roots < 1 .or. n_roots > n) then
         call error%set(ERROR_VALIDATION, "asked the paired response solver for "// &
                        to_char(n_roots)//" roots of a problem with "//to_char(n)// &
                        " dimensions")
         return
      end if

      loud = .false.
      if (present(verbose)) loud = verbose

      ! How wide the guess is, before anything is allocated: the subspace cap
      ! has to be able to hold it, and on a real case `n` is large enough that
      ! an array sized to find out would be the largest one in the program.
      n_new = guess_count(diagonal, n_roots)
      nmax = subspace_cap(n, n_roots, n_new, tol, max_subspace)
      allocate (basis(n, nmax), sp(n, nmax), sm(n, nmax))
      call fill_guess(diagonal, n_new, basis)
      allocate (omega(n_roots), residuals(n_roots), root_converged(n_roots))
      allocate (xpy(n, n_roots), xmy(n, n_roots), correction(n))
      omega = 0.0_dp
      residuals = huge(1.0_dp)
      xpy = 0.0_dp
      xmy = 0.0_dp

      if (present(label)) then
         call convergence_header(loud, label, "    iter          omega     "// &
                                 "residual   subspace   products       time", 74)
      else
         call convergence_header(loud, "RPA iterations", "    iter          omega     "// &
                                 "residual   subspace   products       time", 74)
      end if
      call system_clock(last, rate)

      nsub = 0
      first_new = 1
      best_worst = huge(1.0_dp)
      stall = 0

      do iteration = 1, iterations
         iterations_taken = iteration

         ! Both halves, for the vectors added since the last pass only. The
         ! images of the older ones are already in `sp` and `sm`, and stay
         ! correct across a collapse because a collapsed basis is a linear
         ! combination of the basis whose images those are.
         call operator%apply_plus(basis(:, first_new:first_new + n_new - 1), &
                                  sp(:, first_new:first_new + n_new - 1), error)
         if (error%has_error()) return
         call operator%apply_minus(basis(:, first_new:first_new + n_new - 1), &
                                   sm(:, first_new:first_new + n_new - 1), error)
         if (error%has_error()) return
         products = products + 2*n_new
         nsub = first_new + n_new - 1

         allocate (h1(nsub, nsub), h2(nsub, nsub))
         call pic_gemm(basis(:, 1:nsub), sp(:, 1:nsub), h1, transa="T")
         call pic_gemm(basis(:, 1:nsub), sm(:, 1:nsub), h2, transa="T")
         h1 = 0.5_dp*(h1 + transpose(h1))
         h2 = 0.5_dp*(h2 + transpose(h2))

         ! `(A-B)^{1/2}`, and the reason the whole method needs a stable
         ! reference. A projected eigenvalue at or below zero says the energy
         ! is not a minimum along some rotation the subspace holds, and no
         ! excitation energy follows from it.
         call square_root(h2, root, error)
         if (error%has_error()) then
            deallocate (h1, h2)
            return
         end if

         allocate (work(nsub, nsub), hss(nsub, nsub), w2(nsub))
         call pic_gemm(h1, root, work)
         call pic_gemm(root, work, hss)
         hss = 0.5_dp*(hss + transpose(hss))
         call symmetric_eigen(hss, w2, error)
         if (error%has_error()) then
            deallocate (h1, h2, root, work, hss, w2)
            return
         end if

         ! A negative `w^2` is not a small root to be stepped over. By
         ! Sylvester's law of inertia it is a negative eigenvalue of the
         ! projected `(A+B)`, which `square_root` cannot see because that only
         ! tests `(A-B)`: the two halves of one instability, and only one of
         ! them was being caught. Answering from the positive roots above it
         ! returns a converged, plausible spectrum with the lowest state
         ! missing and every index below it shifted up -- no error, no
         ! warning. It is refused here instead.
         !
         ! A root merely *at* zero is different and is stepped over, not
         ! refused: those are the rotations the reference is flat along, and
         ! an unrestricted reference carries one by construction. The
         ! threshold separates the two, scaled by the spectrum so it means the
         ! same thing for a valence problem and a core one.
         unstable_tol = -max(OMEGA2_FLOOR, epsilon(1.0_dp)*maxval(abs(w2)))
         if (any(w2 < unstable_tol)) then
            call error%set(ERROR_GENERIC, "the reference is unstable ((A+B) is not "// &
                           "positive definite): the projected matrix has a squared "// &
                           "frequency of "//to_char(minval(w2))//", so the excitation "// &
                           "energy there is imaginary. There is no excitation "// &
                           "spectrum of an unstable reference; reconverge it -- "// &
                           "keywords.scf.stability will say which rotation, and "// &
                           "keywords.scf.second_order can escape it -- before asking "// &
                           "for one")
            deallocate (h1, h2, root, work, hss, w2)
            return
         end if

         ! Ascending, so everything at or below the floor is at the front and
         ! the roots wanted are the next `n_roots`.
         n_positive = count(w2 > OMEGA2_FLOOR)
         n_flat = nsub - n_positive
         if (n_flat > 0 .and. loud) then
            write (line, '(a,i0,a)') "   skipped ", n_flat, &
               " rotation(s) of the reference at zero frequency"
            call logger%debug(trim(line))
         end if
         if (n_positive < n_roots) then
            call error%set(ERROR_GENERIC, "the paired response subspace holds only "// &
                           to_char(n_positive)//" positive squared frequencies, and "// &
                           to_char(n_roots)//" roots were asked for: the rest are "// &
                           "imaginary, which is an instability of the reference "// &
                           "rather than an excitation")
            deallocate (h1, h2, root, work, hss, w2)
            return
         end if

         allocate (rss(nsub, n_roots), lss(nsub, n_roots), kept(n_roots))
         do k = 1, n_roots
            kept(k) = sqrt(w2(nsub - n_positive + k))
         end do
         ! `R = H2^{1/2} T` in subspace coordinates, `L = (A+B) R / w`.
         call pic_gemm(root, hss(:, nsub - n_positive + 1:nsub - n_positive + n_roots), &
                       rss)
         call pic_gemm(h1, rss, lss)
         do k = 1, n_roots
            lss(:, k) = lss(:, k)/kept(k)
         end do
         call biorthonormalise(rss, lss, error)
         if (error%has_error()) then
            deallocate (h1, h2, root, work, hss, w2, rss, lss, kept)
            return
         end if
         omega = kept

         ! The vectors themselves and their paired residuals, all four of
         ! which are the stored products read through the subspace coefficients
         ! -- no new application of anything.
         call pic_gemm(basis(:, 1:nsub), rss, xpy)
         call pic_gemm(basis(:, 1:nsub), lss, xmy)
         allocate (wl(n, n_roots), wr(n, n_roots))
         call pic_gemm(sp(:, 1:nsub), rss, wl)
         call pic_gemm(sm(:, 1:nsub), lss, wr)
         do k = 1, n_roots
            wl(:, k) = wl(:, k) - omega(k)*xmy(:, k)
            wr(:, k) = wr(:, k) - omega(k)*xpy(:, k)
            residuals(k) = max(sqrt(dot_product(wr(:, k), wr(:, k))), &
                               sqrt(dot_product(wl(:, k), wl(:, k))))
            root_converged(k) = residuals(k) < tol
         end do
         worst = maxval(residuals)

         if (loud) then
            call system_clock(tick)
            write (line, "(i8,f15.9,es13.3,i11,i11,f11.2,a)") iteration, omega(1), &
               worst, nsub, products, real(tick - last, dp)/real(rate, dp), " s"
            call logger%info(trim(line))
            ! Redirected output is block buffered; without this the rows sit in
            ! a 4 kB buffer and arrive in lumps.
            flush (output_unit)
            last = tick
         end if

         if (all(root_converged)) then
            converged = .true.
            deallocate (h1, h2, root, work, hss, w2, rss, lss, kept, wl, wr)
            exit
         end if

         ! Collapse before growing, so the corrections below are added to the
         ! fresh subspace rather than thrown away by it.
         if (nsub + 2*n_roots > nmax) then
            call collapse(basis, sp, sm, rss, lss, nsub, error)
            if (error%has_error()) then
               deallocate (h1, h2, root, work, hss, w2, rss, lss, kept, wl, wr)
               return
            end if
         end if

         ! Two corrections per unconverged root, which is the cap: `2*n_roots`.
         added = 0
         do k = 1, n_roots
            if (root_converged(k)) cycle
            if (nsub + added >= nmax) exit
            do i = 1, 2
               if (nsub + added >= nmax) exit
               if (i == 1) then
                  correction = wr(:, k)
               else
                  correction = wl(:, k)
               end if
               call precondition(correction, omega(k), diagonal)
               call orthonormalise_against(basis, nsub + added, correction, norm)
               if (norm < LINEAR_DEPENDENCE) cycle
               added = added + 1
               basis(:, nsub + added) = correction
            end do
         end do

         deallocate (h1, h2, root, work, hss, w2, rss, lss, kept, wl, wr)

         if (worst < best_worst) then
            best_worst = worst
            stall = 0
         else
            stall = stall + 1
         end if

         if (added == 0) then
            call error%set(ERROR_GENERIC, "the paired response solve stopped short: "// &
                           "no correction survived projection out of a subspace of "// &
                           to_char(nsub)//" vectors, with the worst residual at "// &
                           to_char(worst)//" against a tolerance of "//to_char(tol)// &
                           ". The subspace already spans everything the residuals "// &
                           "point at, so more iterations will not help; loosen "// &
                           "keywords.excited_states.tolerance")
            exit
         end if
         if (stall >= STAGNATION_WINDOW .and. iteration > STAGNATION_GRACE) then
            call error%set(ERROR_GENERIC, "the paired response solve stagnated: the "// &
                           "worst residual has not improved on "//to_char(best_worst)// &
                           " in "//to_char(STAGNATION_WINDOW)//" iterations and the "// &
                           "tolerance asked for is "//to_char(tol)//". Something "// &
                           "below the tolerance is noise -- a quadrature grid, most "// &
                           "often -- so raise keywords.excited_states.tolerance to "// &
                           "above what the reference can resolve")
            exit
         end if

         ! A collapse above reset `nsub` and the corrections went in on top of
         ! it, so what is owed a product is exactly what was just added.
         first_new = nsub + 1
         n_new = added
      end do

      call convergence_footer(loud, converged, iterations_taken, "iterations", 74)

      ! The other two ways out of the loop -- stagnation and an exhausted
      ! subspace -- both raise. Falling out on `max_iterations` used to return
      ! normally with `converged` false, so a caller that forgot to test it
      ! got partly-converged roots that look like an answer. The shipped
      ! caller does test it; a reusable solver should not depend on that.
      if (.not. converged .and. .not. error%has_error()) then
         call error%set(ERROR_GENERIC, "the paired response solver used all "// &
                        to_char(iterations)//" iterations without converging; the "// &
                        "worst residual was "//to_char(best_worst)//" against a "// &
                        "tolerance of "//to_char(tol)//". Raise "// &
                        "keywords.excited_states.max_iter, or loosen "// &
                        "keywords.excited_states.tolerance")
      end if

      deallocate (basis, sp, sm, correction, root_converged)
   end subroutine rpa_solve

   function guess_count(diagonal, n_roots) result(n_guess)
      !! How many unit vectors the starting space needs
      !!
      !! `GUESS_PER_ROOT` per root, extended over every diagonal element
      !! within `DEGENERACY_WINDOW` of the last one taken so a degenerate
      !! block is not cut in half, and capped by `MAX_EXTRA_GUESS` and by the
      !! space itself. Counted before anything is allocated, because the
      !! subspace cap is sized from it.
      real(dp), intent(in) :: diagonal(:)
      integer, intent(in) :: n_roots
      integer :: n_guess

      logical, allocatable :: taken(:)
      real(dp) :: edge
      integer :: n, want, k, pick

      n = size(diagonal)
      want = min(GUESS_PER_ROOT*n_roots, n)

      allocate (taken(n))
      taken = .false.
      edge = 0.0_dp
      do k = 1, want
         pick = lowest_free(diagonal, taken)
         taken(pick) = .true.
         edge = diagonal(pick)
      end do
      deallocate (taken)

      n_guess = min(count(diagonal <= edge + DEGENERACY_WINDOW), &
                    want + MAX_EXTRA_GUESS, n)
      n_guess = max(n_guess, want)
   end function guess_count

   subroutine fill_guess(diagonal, n_guess, basis)
      !! Unit vectors on the `n_guess` lowest diagonal elements
      real(dp), intent(in) :: diagonal(:)
      integer, intent(in) :: n_guess
      real(dp), intent(inout) :: basis(:, :)
         !! (n, >= n_guess); columns `1:n_guess` are written and the rest left

      logical, allocatable :: taken(:)
      integer :: k, pick

      allocate (taken(size(diagonal)))
      taken = .false.
      basis(:, 1:n_guess) = 0.0_dp
      do k = 1, n_guess
         pick = lowest_free(diagonal, taken)
         taken(pick) = .true.
         basis(pick, k) = 1.0_dp
      end do
      deallocate (taken)
   end subroutine fill_guess

   function lowest_free(diagonal, taken) result(pick)
      !! The smallest diagonal element not already spoken for
      real(dp), intent(in) :: diagonal(:)
      logical, intent(in) :: taken(:)
      integer :: pick
         !! Zero when every element has been taken

      real(dp) :: best
      integer :: i

      pick = 0
      best = huge(1.0_dp)
      do i = 1, size(diagonal)
         if (taken(i)) cycle
         if (diagonal(i) < best) then
            best = diagonal(i)
            pick = i
         end if
      end do
   end function lowest_free

   function subspace_cap(n, n_roots, n_guess, tol, requested) result(nmax)
      !! How many trial vectors to hold before collapsing
      !!
      !! The rule, absent an explicit request, is Psi4's: a multiple of the
      !! root count that grows as the tolerance tightens, because a tighter
      !! solve takes more iterations and each one adds two vectors per root.
      !! It knows nothing about `n`, so it is then cut back to what
      !! `SUBSPACE_MEMORY_LIMIT` allows for the basis and its two product sets,
      !! and raised again to hold the guess plus one round of corrections --
      !! below that the solver collapses on its first iteration and every
      !! iteration after, which converges nothing.
      integer, intent(in) :: n, n_roots, n_guess
      real(dp), intent(in) :: tol
      integer, intent(in), optional :: requested
         !! Absent or zero takes the rule; anything else is honoured as given,
         !! subject only to the floor and to `n`.
      integer :: nmax

      integer :: by_memory, floor_vectors

      nmax = int(-log10(max(tol, tiny(1.0_dp)))*50.0_dp)*n_roots
      by_memory = int(SUBSPACE_MEMORY_LIMIT/(3.0_dp*8.0_dp*real(max(n, 1), dp)))
      nmax = min(nmax, max(by_memory, 1))
      if (present(requested)) then
         if (requested > 0) nmax = requested
      end if

      ! The floor: the guess, plus two corrections per root on top of it.
      floor_vectors = min(n_guess + 2*n_roots, n)
      nmax = max(nmax, floor_vectors)
      nmax = min(nmax, n)
   end function subspace_cap

   subroutine square_root(h, root, error)
      !! `H^{1/2}` of a small symmetric matrix, or the instability by name
      !!
      !! By eigendecomposition rather than a Cholesky, because the
      !! decomposition is what says *why* it failed: a Cholesky reports only
      !! that a leading minor was not positive, while the eigenvalue that went
      !! negative is the rotation the reference is unstable along.
      real(dp), intent(in) :: h(:, :)
      real(dp), allocatable, intent(out) :: root(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: vectors(:, :), values(:), scaled(:, :)
      integer :: m, k

      m = size(h, 1)
      allocate (vectors(m, m), values(m), scaled(m, m), root(m, m))
      vectors = h
      call symmetric_eigen(vectors, values, error)
      if (error%has_error()) then
         deallocate (vectors, values, scaled)
         return
      end if

      ! Relative to the spectrum, not against zero: an eigenvalue of 1e-18
      ! passes an absolute test, `sqrt` returns 1e-9, and the root is then
      ! numerically the square root of nothing -- `S H1 S` comes back as noise
      ! with no error raised. A near-instability is exactly what this is here
      ! to catch, so the threshold has to be able to see one.
      if (any(values <= max(1.0e-12_dp, epsilon(1.0_dp)*maxval(abs(values))))) then
         call error%set(ERROR_GENERIC, "the reference is unstable ((A-B) is not "// &
                        "positive definite): the projected matrix has an eigenvalue "// &
                        "of "//to_char(minval(values))//", so the energy falls along "// &
                        "one of the orbital rotations the response operator spans. "// &
                        "There is no excitation spectrum of an unstable reference; "// &
                        "reconverge the reference -- keywords.scf.stability will say "// &
                        "which rotation -- before asking for one")
         deallocate (vectors, values, scaled, root)
         return
      end if

      do k = 1, m
         scaled(:, k) = vectors(:, k)*sqrt(values(k))
      end do
      call pic_gemm(scaled, vectors, root, transb="T")
      deallocate (vectors, values, scaled)
   end subroutine square_root

   subroutine symmetric_eigen(a, values, error)
      !! Eigenvalues ascending, eigenvectors in place, of a small symmetric matrix
      real(dp), intent(inout) :: a(:, :)
      real(dp), intent(out) :: values(:)
      type(error_t), intent(inout) :: error

      integer :: info

      call pic_syev(a, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_GENERIC, "a subspace matrix of the paired response "// &
                        "solver could not be diagonalized (info = "//to_char(info)//")")
      end if
   end subroutine symmetric_eigen

   subroutine biorthonormalise(rss, lss, error)
      !! Scale each pair of columns so that `R . L` is one
      !!
      !! `R . L` is `|X|^2 - |Y|^2`, which is the norm the paired problem
      !! conserves and the only one it does -- neither `R` nor `L` has a length
      !! of its own worth fixing. Both columns take the same factor, so their
      !! ratio, which is what carries the physics, is untouched.
      real(dp), intent(inout) :: rss(:, :), lss(:, :)
      type(error_t), intent(inout) :: error

      real(dp) :: inner
      integer :: k

      do k = 1, size(rss, 2)
         inner = dot_product(rss(:, k), lss(:, k))
         if (inner <= 0.0_dp) then
            call error%set(ERROR_GENERIC, "a paired response root has |X|^2 - |Y|^2 "// &
                           "of "//to_char(inner)//", which is not positive: the "// &
                           "excitation carries more de-excitation than excitation "// &
                           "character, which a stable reference cannot produce")
            return
         end if
         rss(:, k) = rss(:, k)/sqrt(inner)
         lss(:, k) = lss(:, k)/sqrt(inner)
      end do
   end subroutine biorthonormalise

   subroutine precondition(vector, w, diagonal)
      !! Divide a residual by `w - diagonal`, in place
      !!
      !! The Davidson correction, which would be exact if the operator were
      !! its own diagonal. Where the denominator is too small to divide by, the
      !! residual is taken as it stands: it is still a direction the subspace
      !! does not hold, and amplifying one element of it by ten thousand is not
      !! a better one.
      real(dp), intent(inout) :: vector(:)
      real(dp), intent(in) :: w
      real(dp), intent(in) :: diagonal(:)

      real(dp) :: denominator
      integer :: i

      do i = 1, size(vector)
         denominator = w - diagonal(i)
         if (abs(denominator) >= PRECONDITION_FLOOR) vector(i) = vector(i)/denominator
      end do
   end subroutine precondition

   subroutine orthonormalise_against(basis, nsub, vector, norm)
      !! Project the first `nsub` basis columns out of `vector` and normalise it
      !!
      !! Normalised first, so the length returned measures how much of a *new
      !! direction* is there rather than how large the residual was -- a test
      !! that scaled with the residual would stop the subspace growing exactly
      !! when the solve needs it to. Twice, because one pass of Gram-Schmidt
      !! loses orthogonality when the vector is nearly in the subspace already.
      real(dp), intent(in) :: basis(:, :)
      integer, intent(in) :: nsub
      real(dp), intent(inout) :: vector(:)
      real(dp), intent(out) :: norm
         !! What survived, as a fraction of the incoming length. Below
         !! `LINEAR_DEPENDENCE` the vector is not worth adding and has not been
         !! normalised.

      real(dp) :: overlap
      integer :: j, pass

      norm = sqrt(dot_product(vector, vector))
      if (norm < tiny(1.0_dp)) then
         norm = 0.0_dp
         return
      end if
      vector = vector/norm

      do pass = 1, 2
         do j = 1, nsub
            overlap = dot_product(basis(:, j), vector)
            vector = vector - overlap*basis(:, j)
         end do
      end do

      norm = sqrt(dot_product(vector, vector))
      if (norm >= LINEAR_DEPENDENCE) vector = vector/norm
   end subroutine orthonormalise_against

   subroutine collapse(basis, sp, sm, rss, lss, nsub, error)
      !! Restart the subspace on the current `R` and `L`, products and all
      !!
      !! Both are kept, not just `R`: the paired problem needs the two, and a
      !! collapse onto `R` alone throws away the half of the current solution
      !! that carries the de-excitation amplitudes and has to rebuild it.
      !!
      !! The products survive exactly. Every new basis vector is a linear
      !! combination of the old ones with known images, so the same combination
      !! of `sp` and `sm` is its image -- no application of the operator is
      !! needed here, which is the whole reason the collapse is written in
      !! subspace coordinates.
      real(dp), intent(inout) :: basis(:, :), sp(:, :), sm(:, :)
      real(dp), intent(in) :: rss(:, :), lss(:, :)
         !! (nsub, n_roots) the subspace coordinates of `R` and `L`
      integer, intent(inout) :: nsub
         !! In, the subspace dimension; out, what is left after collapsing
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: candidate(:, :), keep(:, :), image(:, :)
      real(dp) :: norm
      integer :: n_roots, m, k, taken

      n_roots = size(rss, 2)
      m = 2*n_roots
      allocate (candidate(nsub, m), keep(nsub, m))
      candidate(:, 1:n_roots) = rss
      candidate(:, n_roots + 1:m) = lss

      ! Orthonormalised in subspace coordinates, which is the same thing as in
      ! the full space because the basis is orthonormal.
      taken = 0
      do k = 1, m
         call orthonormalise_against(keep, taken, candidate(:, k), norm)
         if (norm < LINEAR_DEPENDENCE) cycle
         taken = taken + 1
         keep(:, taken) = candidate(:, k)
      end do
      if (taken < 1) then
         call error%set(ERROR_GENERIC, "collapsing the paired response subspace left "// &
                        "no independent vectors, which cannot happen for a solve that "// &
                        "has roots")
         deallocate (candidate, keep)
         return
      end if

      allocate (image(size(basis, 1), taken))
      call pic_gemm(basis(:, 1:nsub), keep(:, 1:taken), image)
      basis(:, 1:taken) = image
      call pic_gemm(sp(:, 1:nsub), keep(:, 1:taken), image)
      sp(:, 1:taken) = image
      call pic_gemm(sm(:, 1:nsub), keep(:, 1:taken), image)
      sm(:, 1:taken) = image
      nsub = taken

      deallocate (candidate, keep, image)
   end subroutine collapse

end module mqc_czt_rpa_solver
