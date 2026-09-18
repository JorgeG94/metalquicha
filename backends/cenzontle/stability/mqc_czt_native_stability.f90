!! The SCF stability analysis, on this program's own Davidson
module mqc_czt_native_stability
   !! Ask whether a converged SCF is a minimum, and say which way is down if
   !! it is not -- with nothing optional in the build.
   !!
   !! The matrix is `mqc_czt_ov_hessian`'s. What diagonalises it is
   !! `mqc_davidson`, which is already in the program for the CI -- so a
   !! stability analysis and the second-order SCF that rests on the same
   !! Hessian are core functionality and work in a build with nothing optional
   !! in it. This is the path `keywords.scf.stability` takes and the one to add
   !! to. `native` is in the name because a second implementation over the same
   !! operator, borrowed from a library rather than written here, exists on a
   !! separate branch purely as a cross-check.
   !!
   !! ## Why `mqc_davidson` fits without being bent
   !!
   !! `davidson_flat` asks an operator for a matrix-vector product and a
   !! diagonal to precondition with, over a flat vector, and returns the lowest
   !! eigenpairs. That is the entire question a stability analysis asks. Three
   !! things about it are worth stating, because they are what a reader would
   !! otherwise have to check:
   !!
   !! * **Nothing in it is CI-specific.** `sigma_operator_t` is abstract and
   !!   the shaped `davidson_lowest` wrapper is where the determinants live;
   !!   the flat entry point below it knows only `diagonal` and `apply`.
   !! * **The preconditioner is the right one.** Davidson divides the residual
   !!   by `theta - diagonal`, and the diagonal here is the orbital-energy
   !!   differences -- the standard level-shifted diagonal preconditioner for
   !!   this operator, and the approximation the electronic Hessian is
   !!   diagonally dominant enough for.
   !! * **The starting vector lands where an instability is.** With no guess,
   !!   `initial_basis` takes the unit vector on the smallest diagonal element,
   !!   which over this space is the HOMO-LUMO rotation.
   !!
   !! Two properties fall out of using it rather than a library with a bare
   !! callback. `apply` carries an `error_t`, so a failed Fock build travels
   !! with its message and this module needs **no module-level state** -- no
   !! pointer to say which operator a callback is applying, no re-entry flag,
   !! and no reason two analyses cannot run at once. And Davidson's start is
   !! the diagonal rather than a random subspace, so an analysis is reproducible
   !! without a seed.
   !!
   !! ## What is reported
   !!
   !! The eigenvalue, **always** -- not only when the reference turns out to be
   !! a saddle. `mqc_czt_stability` can only recover it for an unstable
   !! reference, because the library it calls returns the eigenvector and not
   !! the eigenvalue; that is a gap in an interface and not a fact about the
   !! physics, so it is not reproduced here. The descent direction is still
   !! returned only when there is one: a minimum has no downhill rotation.
   !!
   !! ## What this does and does not decide
   !!
   !! `mqc_czt_ov_hessian`'s header, unchanged: this is the **real singlet**
   !! orbital-rotation Hessian, so a stable verdict is the RHF-to-RHF question
   !! and says nothing about a triplet or a complex instability.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   use mqc_program_limits, only: MAX_LINE_LENGTH
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_xc, only: xc_context_t
   use mqc_czt_hessian, only: nuclear_response_t
   use mqc_czt_ov_hessian, only: ov_hessian_t, stability_result_t, build_scf_ov_hessian
   use mqc_davidson, only: davidson_flat, sigma_operator_t
   implicit none
   private

   public :: native_scf_stability
   public :: native_stability_of_hessian
   public :: hessian_sigma_t

   type, extends(sigma_operator_t) :: hessian_sigma_t
      !! `ov_hessian_t` as something `mqc_davidson` will accept
      !!
      !! The whole adapter. `ov_hessian_t%apply` leaves its message on the
      !! operator because the eigensolver it was written for had nowhere to put
      !! one; `sigma_operator_t%apply` carries an `error_t`, so this moves the
      !! message onto it and the caller reads a failed Fock build the same way
      !! it reads every other failure.
      type(ov_hessian_t), pointer :: hessian => null()
   contains
      procedure :: apply => hessian_sigma_apply
   end type hessian_sigma_t

contains

   subroutine hessian_sigma_apply(this, vector, image, error)
      !! `(A+B) x`, with the message moved onto the error the solver carries
      class(hessian_sigma_t), intent(inout) :: this
      real(dp), intent(in) :: vector(:)
      real(dp), intent(out) :: image(:)
      type(error_t), intent(inout) :: error

      if (error%has_error()) return
      if (.not. associated(this%hessian)) then
         call error%set(ERROR_VALIDATION, "the stability analysis' eigensolver was "// &
                        "handed no electronic Hessian to apply")
         image = 0.0_dp
         return
      end if

      call this%hessian%apply(vector, image)
      if (this%hessian%error%has_error()) then
         error = this%hessian%error
         call error%add_context("applying the electronic Hessian for the "// &
                                "stability analysis")
      end if
   end subroutine hessian_sigma_apply

   subroutine native_scf_stability(mol, orbitals, energies, n_occ, result, error, xc, &
                                   reference, k_scale, rs_k_lr, rs_omega, bounds, &
                                   conv_tol, max_iter, verbose)
      !! Is this converged SCF a minimum?
      !!
      !! Builds the electronic Hessian over the non-redundant occupied-virtual
      !! rotations and hands it to `native_stability_of_hessian`. The optional
      !! exchange-correlation arguments are the ones the response operator
      !! needs for a Kohn-Sham reference; without them this is the
      !! Hartree-Fock response of whatever orbitals it was given, which for a
      !! Kohn-Sham reference is the wrong operator rather than an approximate
      !! one.
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      type(stability_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: reference(:, :)
      real(dp), intent(in), optional :: k_scale
      real(dp), intent(in), optional :: rs_k_lr, rs_omega
      real(dp), intent(in), optional :: bounds(:, :)
      real(dp), intent(in), optional :: conv_tol
      integer, intent(in), optional :: max_iter
      logical, intent(in), optional :: verbose

      ! `target`, and load-bearing: `hessian%response` points at this, so it
      ! has to be an entity that outlives the call it is used in.
      type(nuclear_response_t), target :: response
      type(ov_hessian_t), target :: hessian

      if (error%has_error()) return

      call build_scf_ov_hessian(mol, orbitals, energies, n_occ, response, hessian, &
                                error, xc=xc, reference=reference, k_scale=k_scale, &
                                rs_k_lr=rs_k_lr, rs_omega=rs_omega, bounds=bounds)
      if (error%has_error()) return

      call native_stability_of_hessian(hessian, result, error, conv_tol=conv_tol, &
                                       max_iter=max_iter, verbose=verbose)
   end subroutine native_scf_stability

   subroutine native_stability_of_hessian(hessian, result, error, conv_tol, max_iter, &
                                          max_subspace, verbose)
      !! The lowest eigenvalue of an electronic Hessian, and its sign
      !!
      !! Separate from `native_scf_stability` so that the eigensolver can be
      !! exercised against an operator whose matrix is known, which is how the
      !! convention this program works in is checked -- in
      !! `test/test_mqc_czt_stability.f90`, against a dense diagonalisation of
      !! the same operator.
      type(ov_hessian_t), intent(inout), target :: hessian
      type(stability_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: conv_tol
         !! Norm of the residual `(A+B) v - theta v` the lowest eigenpair is
         !! converged to. Also what decides the verdict: a curvature more
         !! negative than this is a real negative eigenvalue, one between
         !! `-conv_tol` and zero is the solver's own resolution and is not
         !! called an instability.
      integer, intent(in), optional :: max_iter
      integer, intent(in), optional :: max_subspace
         !! Trial vectors kept before the subspace is collapsed. Davidson's own
         !! default is sixteen, capped by the size of the rotation space.
      logical, intent(in), optional :: verbose
         !! A line per Davidson iteration. Off by default; every one of those
         !! lines is a Fock build.

      type(hessian_sigma_t) :: operator
      real(dp), allocatable :: diagonal(:), values(:), vectors(:, :), residuals(:)
      real(dp) :: tol
      integer :: n_param, iterations, products
      logical :: converged
      character(len=MAX_LINE_LENGTH) :: line

      if (error%has_error()) return

      n_param = hessian%length()
      if (n_param < 1) then
         call error%set(ERROR_VALIDATION, "there are no non-redundant orbital "// &
                        "rotations, so there is no curvature to examine")
         return
      end if
      diagonal = hessian%diagonal()
      if (size(diagonal) /= n_param) then
         call error%set(ERROR_VALIDATION, "the electronic Hessian's diagonal is not "// &
                        "the length of the rotation space it acts on")
         return
      end if

      tol = 1.0e-6_dp
      if (present(conv_tol)) tol = conv_tol

      call hessian%error%clear()
      hessian%n_apply = 0
      operator%hessian => hessian

      result%n_parameters = n_param

      call davidson_flat(operator, diagonal, 1, values, vectors, residuals, &
                         iterations, products, converged, error, &
                         tolerance=tol, max_iterations=max_iter, &
                         max_subspace=max_subspace, verbose=verbose, &
                         label="stability analysis", value_label="curvature")
      result%n_products = hessian%n_apply
      if (error%has_error()) return
      if (.not. converged) then
         call error%set(ERROR_GENERIC, "the stability analysis did not converge its "// &
                        "lowest eigenpair in the iterations allowed; raise "// &
                        "keywords.scf.stability_maxiter or loosen "// &
                        "keywords.scf.stability_tolerance")
         return
      end if

      result%ran = .true.
      ! Always, which is the point of this path existing. The eigenvalue is
      ! what a reader wants whether or not it turned out to be negative -- how
      ! stiff a minimum is, is as much of an answer as which way a saddle falls.
      result%lowest_curvature = values(1)
      result%has_curvature = .true.
      result%stable = values(1) > -tol
      ! Only for a saddle: there is no downhill direction from a minimum, and
      ! handing back the softest uphill one invites a caller to follow it.
      if (.not. result%stable) result%rotation = vectors(:, 1)

      if (result%stable) then
         write (line, "(a,i0,a,es12.4,a,i0,a)") "  the reference is a minimum with "// &
            "respect to the ", n_param, " real closed-shell orbital rotations: "// &
            "lowest curvature ", result%lowest_curvature, " hartree (", &
            result%n_products, " Hessian-vector products)"
         call logger%info(trim(line))
      else
         write (line, "(a,es12.4,a,i0,a)") "  the reference is a saddle point: lowest "// &
            "orbital-rotation curvature ", result%lowest_curvature, &
            " hartree (", result%n_products, " Hessian-vector products)"
         call logger%warning(trim(line))
      end if
   end subroutine native_stability_of_hessian

end module mqc_czt_native_stability
