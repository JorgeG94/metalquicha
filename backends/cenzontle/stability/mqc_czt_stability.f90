!! The SCF stability analysis, through OpenTrustRegion
module mqc_czt_stability
   !! Ask whether a converged SCF is a minimum, and say which way is down if
   !! it is not.
   !!
   !! The matrix is `mqc_czt_ov_hessian`'s; this module only diagonalises it,
   !! and it does not do that itself either. OpenTrustRegion's
   !! `stability_check` is a matrix-free Davidson -- Jacobi-Davidson once the
   !! subspace stops paying -- for the lowest eigenvalue of a symmetric
   !! operator, which is exactly the question, and it needs only a
   !! Hessian-vector product and a diagonal. No orbital update, no line search,
   !! no objective function: that is why it is the cheapest possible first
   !! consumer of the library.
   !!
   !! ## The module state below, and why it has to exist
   !!
   !! `stability_check` takes the Hessian-vector product as a bare procedure
   !! pointer:
   !!
   !!     subroutine hess_x_type(x, hess_x, error)
   !!
   !! There is no host-context argument. A Fortran procedure pointer cannot
   !! carry a closure, so the callback has no way to be told *which* Hessian it
   !! is applying except by looking at something outside its own arguments. The
   !! two private module variables below are that something, and they are as
   !! small as the constraint permits: a pointer to the operator in play, and a
   !! flag that refuses a second analysis while one is running rather than
   !! letting it quietly repoint the first one's callback.
   !!
   !! Everything else that would otherwise have joined them is on the operator
   !! instead -- `ov_hessian_t%error` carries the message a failed Fock build
   !! produced, and `ov_hessian_t%n_apply` the count -- so the state is a
   !! pointer and a boolean and not a pile of parallel arrays.
   !!
   !! **One argument upstream would remove both.** A `class(*)` or
   !! `type(c_ptr)` context passed through `stability_check` to `hess_x` is all
   !! it would take; the library is being consumed as its interface stands
   !! today, deliberately, and this note is here so that the state is read as
   !! the cost of that decision rather than as a design.
   !!
   !! The flag is not a substitute for thread safety and does not pretend to
   !! be: two analyses in two threads would race on it. Nothing here is called
   !! from inside a parallel region, and the Fock build underneath is itself
   !! threaded, so there is no reason for one to be.
   !!
   !! ## What the library does not report, and what is done about it
   !!
   !! `stability_check` returns a verdict and, when the verdict is negative,
   !! the descent direction. It does **not** return the eigenvalue: today's
   !! interface writes it into a log message and nothing else. So the lowest
   !! curvature is recovered here as a Rayleigh quotient of the direction the
   !! library did hand back -- one extra Hessian-vector product, exact because
   !! that direction *is* the eigenvector -- and is therefore available only
   !! for an unstable reference. A stable one reports `has_curvature` false
   !! rather than a number nobody computed. This is the other half of the
   !! upstream addition the note above describes.
   use pic_types, only: dp, default_int
   use pic_logger, only: logger => global_logger, error_level, warning_level, &
                         verbose_level, debug_level
   use mqc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   use mqc_program_limits, only: MAX_LINE_LENGTH
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_xc, only: xc_context_t
   use mqc_czt_hessian, only: nuclear_response_t
   use mqc_czt_ov_hessian, only: ov_hessian_t, stability_result_t, build_scf_ov_hessian
   use opentrustregion, only: otr_rp => rp, otr_ip => ip, &
                              stability_settings_type, stability_check, &
                              hess_x_type, logger_type, &
                              otr_error_stability_check => error_stability_check, &
                              otr_error_max_iter => error_stability_check_max_iter, &
                              otr_error_hess_x => error_hess_x, &
                              otr_silent => verbosity_silent, &
                              otr_error => verbosity_error, &
                              otr_warning => verbosity_warning, &
                              otr_info => verbosity_info, &
                              otr_debug => verbosity_debug
   implicit none
   private

   public :: scf_stability
   public :: stability_of_hessian
   public :: otr_available

   ! The whole of the module state the callback interface forces. See the
   ! header: one pointer to say which operator the callback is applying, one
   ! flag to refuse re-entering while it is set.
   type(ov_hessian_t), pointer :: active_hessian => null()
   logical :: analysis_running = .false.

   ! What the callback returns for the two things that can go wrong in it.
   ! OpenTrustRegion requires a positive code below 100 and adds its own origin
   ! to it, so these come back as `error_hess_x + code`.
   integer(otr_ip), parameter :: CALLBACK_NO_OPERATOR = 1_otr_ip
   integer(otr_ip), parameter :: CALLBACK_APPLY_FAILED = 2_otr_ip

contains

   pure function otr_available() result(available)
      !! Whether this build can run a stability analysis
      logical :: available

      available = .true.
   end function otr_available

   subroutine scf_stability(mol, orbitals, energies, n_occ, result, error, xc, &
                            reference, k_scale, rs_k_lr, rs_omega, bounds, &
                            conv_tol, max_iter, n_trial_vectors, seed)
      !! Is this converged SCF a minimum?
      !!
      !! Builds the electronic Hessian over the non-redundant occupied-virtual
      !! rotations and hands it to `stability_of_hessian`. The optional
      !! exchange-correlation arguments are the ones the response operator
      !! needs for a Kohn-Sham reference and mean the same thing here as they
      !! do in `mqc_czt_hessian`; without them this is the Hartree-Fock
      !! response of whatever orbitals it was given, which for a Kohn-Sham
      !! reference is the wrong operator rather than an approximate one.
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
      integer, intent(in), optional :: n_trial_vectors
      integer, intent(in), optional :: seed

      ! `target`, and load-bearing: `hessian%response` points at this, so it
      ! has to be an entity that outlives the call it is used in.
      type(nuclear_response_t), target :: response
      type(ov_hessian_t) :: hessian

      if (error%has_error()) return

      call build_scf_ov_hessian(mol, orbitals, energies, n_occ, response, hessian, &
                                error, xc=xc, reference=reference, k_scale=k_scale, &
                                rs_k_lr=rs_k_lr, rs_omega=rs_omega, bounds=bounds)
      if (error%has_error()) return

      call stability_of_hessian(hessian, result, error, conv_tol=conv_tol, &
                                max_iter=max_iter, n_trial_vectors=n_trial_vectors, &
                                seed=seed)
   end subroutine scf_stability

   subroutine stability_of_hessian(hessian, result, error, conv_tol, max_iter, &
                                   n_trial_vectors, seed)
      !! The lowest eigenvalue of an electronic Hessian, and its sign
      !!
      !! Separate from `scf_stability` so that the eigensolver can be exercised
      !! against an operator whose matrix is known -- which is how the
      !! convention this program hands OpenTrustRegion is checked, in
      !! `test/test_mqc_czt_stability.f90`, against a dense diagonalisation of
      !! the same operator.
      type(ov_hessian_t), intent(inout), target :: hessian
      type(stability_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: conv_tol
         !! Root-mean-square residual the lowest eigenpair is converged to.
      integer, intent(in), optional :: max_iter
      integer, intent(in), optional :: n_trial_vectors
         !! Random starting directions beyond the one the diagonal picks.
         !! OpenTrustRegion caps this at half the parameter count itself.
      integer, intent(in), optional :: seed
         !! Those directions are random, so a reproducible analysis needs a
         !! fixed seed. Defaulted by the library rather than left to the clock.

      type(stability_settings_type) :: settings
      procedure(hess_x_type), pointer :: hess_x_funptr
      procedure(logger_type), pointer :: log_funptr
      real(dp), allocatable :: h_diag(:), kappa(:), image(:)
      integer(otr_ip) :: otr_status
      integer :: n_param
      logical :: stable
      real(dp) :: norm2
      character(len=MAX_LINE_LENGTH) :: line

      if (error%has_error()) return

      if (analysis_running) then
         call error%set(ERROR_VALIDATION, "a stability analysis is already running; "// &
                        "the eigensolver's callback carries no context, so only one "// &
                        "can be in flight at a time")
         return
      end if

      n_param = hessian%length()
      if (n_param < 1) then
         call error%set(ERROR_VALIDATION, "there are no non-redundant orbital "// &
                        "rotations, so there is no curvature to examine")
         return
      end if
      h_diag = hessian%diagonal()
      if (size(h_diag) /= n_param) then
         call error%set(ERROR_VALIDATION, "the electronic Hessian's diagonal is not "// &
                        "the length of the rotation space it acts on")
         return
      end if

      call settings%init(otr_status)
      if (otr_status /= 0_otr_ip) then
         call error%set(ERROR_GENERIC, "OpenTrustRegion declined to initialise its "// &
                        "stability settings"//otr_code(otr_status))
         return
      end if

      ! Its own printing goes through this program's logger, so a stability
      ! analysis obeys `--verbose` and the log file like everything else rather
      ! than writing to standard output on its own account. The hook carries
      ! only the text -- the severity is consumed by the library's own
      ! `verbose` filter before the message is handed over -- so the level is
      ! chosen once, here, by what this program's logger is set to.
      log_funptr => otr_log
      settings%logger => log_funptr
      settings%verbose = otr_verbosity()

      ! Left null. The library's internal default preconditioner is already the
      ! level-shifted diagonal one, applied to exactly the `h_diag` below --
      ! the orbital-energy differences -- so supplying `precond` would be this
      ! program reimplementing what it is about to pass in anyway. `project` is
      ! null for a stronger reason: the parameterisation is the non-redundant
      ! rotations and nothing else, so there is no redundant component for a
      ! projection to remove. See `mqc_czt_ov_hessian`.
      if (present(conv_tol)) settings%conv_tol = real(conv_tol, otr_rp)
      if (present(max_iter)) settings%n_iter = int(max_iter, otr_ip)
      if (present(n_trial_vectors)) then
         settings%n_random_trial_vectors = int(n_trial_vectors, otr_ip)
      end if
      if (present(seed)) settings%seed = int(seed, otr_ip)

      hess_x_funptr => otr_hess_x
      call hessian%error%clear()
      hessian%n_apply = 0
      allocate (kappa(n_param))
      kappa = 0.0_dp

      active_hessian => hessian
      analysis_running = .true.
      call stability_check(real(h_diag, otr_rp), hess_x_funptr, stable, otr_status, &
                           settings, kappa)
      analysis_running = .false.
      active_hessian => null()

      result%n_parameters = n_param
      result%n_products = hessian%n_apply

      if (otr_status /= 0_otr_ip) then
         call translate_error(otr_status, hessian, error)
         return
      end if

      result%ran = .true.
      result%stable = stable

      if (.not. stable) then
         ! `kappa` is the eigenvector, so its Rayleigh quotient is the
         ! eigenvalue exactly and not an estimate of it. One more Fock build,
         ! spent because the library does not return the number itself.
         norm2 = dot_product(kappa, kappa)
         if (norm2 > 0.0_dp) then
            allocate (image(n_param))
            call hessian%apply(kappa, image)
            if (hessian%error%has_error()) then
               error = hessian%error
               return
            end if
            result%lowest_curvature = dot_product(kappa, image)/norm2
            result%has_curvature = .true.
            result%n_products = hessian%n_apply
            deallocate (image)
         end if
         result%rotation = kappa
      end if

      if (result%stable) then
         write (line, "(a,i0,a,i0,a)") "  the reference is a minimum with respect to "// &
            "the ", n_param, " real closed-shell orbital rotations (", &
            result%n_products, " Hessian-vector products)"
         call logger%info(trim(line))
      else if (result%has_curvature) then
         write (line, "(a,es12.4,a,i0,a)") "  the reference is a saddle point: lowest "// &
            "orbital-rotation curvature ", result%lowest_curvature, &
            " hartree (", result%n_products, " Hessian-vector products)"
         call logger%warning(trim(line))
      else
         call logger%warning("  the reference is a saddle point with respect to real "// &
                             "closed-shell orbital rotations")
      end if
   end subroutine stability_of_hessian

   function otr_verbosity() result(verbose)
      !! This program's log level, in OpenTrustRegion's five
      !!
      !! Read rather than configured: the library filters its own messages
      !! before the logger hook sees them, so leaving `verbose` at its default
      !! of silent would discard everything on the way out no matter what this
      !! program's logger was set to.
      integer(otr_ip) :: verbose

      integer(default_int) :: level

      call logger%configuration(level)
      if (level >= debug_level) then
         verbose = otr_debug
      else if (level >= verbose_level) then
         verbose = otr_info
      else if (level >= warning_level) then
         verbose = otr_warning
      else if (level >= error_level) then
         verbose = otr_error
      else
         verbose = otr_silent
      end if
   end function otr_verbosity

   subroutine translate_error(status, hessian, error)
      !! An OpenTrustRegion status code as something a user can read
      !!
      !! The library composes a code from an origin and a cause -- the origin
      !! being which of its routines failed, the cause being what the callback
      !! returned -- so the arithmetic below is on ranges rather than equality.
      !! A failure inside the Hessian-vector product is the interesting case:
      !! the cause is one of this program's own, and its message was left on
      !! the operator because the callback had nowhere to put it.
      integer(otr_ip), intent(in) :: status
      type(ov_hessian_t), intent(in) :: hessian
      type(error_t), intent(inout) :: error

      integer(otr_ip) :: cause

      if (status >= otr_error_hess_x .and. status < otr_error_hess_x + 100_otr_ip) then
         cause = status - otr_error_hess_x
         if (cause == CALLBACK_APPLY_FAILED .and. hessian%error%has_error()) then
            error = hessian%error
            call error%add_context("applying the electronic Hessian for the "// &
                                   "stability analysis")
            return
         end if
         if (cause == CALLBACK_NO_OPERATOR) then
            call error%set(ERROR_GENERIC, "the stability analysis' Hessian-vector "// &
                           "product was called with no operator registered, which "// &
                           "means two analyses ran at once")
            return
         end if
         call error%set(ERROR_GENERIC, "the electronic Hessian could not be applied "// &
                        "during the stability analysis"//otr_code(status))
         return
      end if

      if (status == otr_error_max_iter) then
         call error%set(ERROR_GENERIC, "the stability analysis did not converge its "// &
                        "lowest eigenpair in the iterations allowed; raise "// &
                        "keywords.scf.stability_maxiter or loosen "// &
                        "keywords.scf.stability_tolerance"//otr_code(status))
         return
      end if

      if (status >= otr_error_stability_check .and. &
          status < otr_error_stability_check + 100_otr_ip) then
         call error%set(ERROR_GENERIC, "OpenTrustRegion's stability check failed"// &
                        otr_code(status))
         return
      end if

      call error%set(ERROR_GENERIC, "OpenTrustRegion returned an error this program "// &
                     "does not recognise"//otr_code(status))
   end subroutine translate_error

   function otr_code(status) result(text)
      !! " (OpenTrustRegion error <n>)", so a report can be matched to the library
      integer(otr_ip), intent(in) :: status
      character(len=:), allocatable :: text

      character(len=32) :: buffer

      write (buffer, "(i0)") status
      text = " (OpenTrustRegion error "//trim(buffer)//")"
   end function otr_code

   subroutine otr_hess_x(x, hess_x, error)
      !! `(A+B) x`, as OpenTrustRegion's callback interface asks for it
      !!
      !! The interface is fixed by the library and carries no context, so the
      !! operator is reached through the module pointer; the header says why
      !! that pointer exists. The `target` attributes are the library's and are
      !! reproduced because a procedure pointer assignment will not compile
      !! without an exact match.
      real(otr_rp), intent(in), target :: x(:)
      real(otr_rp), intent(out), target :: hess_x(:)
      integer(otr_ip), intent(out) :: error

      error = 0_otr_ip

      if (.not. associated(active_hessian)) then
         hess_x = 0.0_otr_rp
         error = CALLBACK_NO_OPERATOR
         return
      end if

      call active_hessian%apply(real(x, dp), hess_x)
      if (active_hessian%error%has_error()) error = CALLBACK_APPLY_FAILED
   end subroutine otr_hess_x

   subroutine otr_log(message)
      !! OpenTrustRegion's own output, on this program's logger
      !!
      !! At info, not at the severity the library assigned: `print_message`
      !! filters on its own `verbose` setting and then hands the hook the text
      !! alone, so the severity is gone by the time it arrives. What survives
      !! is the filter, which `otr_verbosity` has already matched to this
      !! program's log level.
      character(*), intent(in) :: message

      call logger%info("  opentrustregion:"//trim(message))
   end subroutine otr_log

end module mqc_czt_stability
