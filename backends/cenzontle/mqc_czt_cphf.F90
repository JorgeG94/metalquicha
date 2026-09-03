!! Coupled-perturbed Hartree-Fock, and the static polarizability from it
module mqc_czt_cphf
   !! First-order orbital response to a one-electron perturbation.
   !!
   !! Writing the first-order orbitals as a rotation into the virtual space,
   !! `C_i^(1) = sum_a C_a U_ai`, a closed-shell reference gives
   !!
   !!     (eps_a - eps_i) U_ai + sum_bj [4(ai|bj) - (ab|ij) - (aj|ib)] U_bj = -h_ai
   !!
   !! with `h_ai` the perturbation in the occupied-virtual MO block. Only that
   !! block appears: a rotation among occupied orbitals alone cannot change the
   !! density.
   !!
   !! The two-electron term is one Fock build per iteration rather than an
   !! `n_occ^2 n_vir^2` matrix. Contracting the symmetrized response density
   !!
   !!     Dt = C_vir U C_occ^T + (C_vir U C_occ^T)^T
   !!
   !! as an SCF contracts a density, `G(D) = J - K/2`, gives half the bracket,
   !! so the two-electron contribution is `2 C_vir^T G(Dt) C_occ`.
   !!
   !! Solved by conjugate gradients preconditioned on the orbital-energy
   !! differences. The bracket is positive definite only where the SCF sits at a
   !! minimum, so a failure to converge indicts the reference -- a saddle point
   !! or broken symmetry -- and the guards below say so rather than iterating on.
   !!
   !! **Hartree-Fock unless an `xc` context is passed.** With a functional the
   !! response operator gains the exchange-correlation kernel; a Kohn-Sham
   !! reference solved without one silently gets the Hartree-Fock response of its
   !! orbitals.
   use, intrinsic :: iso_fortran_env, only: output_unit
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_getrf, pic_getrs
   use mqc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   use mqc_czt_integrals, only: czt_molecule_t, ket_transformed_pairs, build_df_mo_block
   use mqc_czt_multipole, only: multipole_matrices
   use mqc_czt_localize, only: boys_localize
   use mqc_czt_rhf, only: build_fock
   use mqc_czt_xc, only: xc_context_t, xc_kernel_apply
   use mqc_czt_direct, only: build_fock_direct, build_fock_direct_nosym, &
                             build_fock_direct_many, &
                             schwarz_bounds, direct_stats_t
   use pic_logger, only: logger => global_logger
   use mqc_program_limits, only: MAX_LINE_LENGTH
   use mqc_calculation_defaults, only: DEFAULT_DYNAMIC_TOL, DEFAULT_DYNAMIC_MAXITER, &
                                       DEFAULT_RESPONSE_BATCH, &
                                       EFP_RESPONSE_AUTO, EFP_RESPONSE_DENSE, &
                                       EFP_RESPONSE_MATRIX_FREE
   implicit none
   private

   public :: response_hessian_t
   public :: cphf_solve
   public :: static_polarizability
   public :: distributed_polarizability
   public :: dynamic_polarizability
   public :: distributed_dynamic_polarizability
   public :: distributed_dynamic_cross
   public :: casimir_polder_frequencies
   public :: build_hessian, build_hessian_df, build_hessian_mo
   public :: dynamic_response_iterative
   public :: static_response_dense
   public :: fitted_potential_general
   ! TODO(mqc): only `cphf_solve` and `static_polarizability` accept an `xc`
   ! context. Every frequency-dependent route -- `dynamic_polarizability`,
   ! `dynamic_response_iterative`, `response_operator_minus` and all three
   ! Hessian builds -- has no kernel at all, so a Kohn-Sham reference gets its
   ! Hartree-Fock response there with nothing in the output to say so.

   type :: response_hessian_t
      !! A built response Hessian, so several blocks can share one build
      !!
      !! `(A-B)` and the product `(A-B)(A+B)` are what a solve at finite
      !! frequency needs. `(A+B)` is kept separately for `nu = 0`, where the
      !! equations reduce to `(A+B) S = -2h` and going through the product would
      !! square the condition number.
      logical :: ready = .false.
      logical :: fitted = .false.
         !! Whether `(A+B)` came from fitted integrals rather than exact ones.
         !! `distributed_polarizability` declines a fitted Hessian and iterates
         !! instead, so the static block is not moved by the fitting error
         !! without the deck having asked for it.
      real(dp), allocatable :: aminus(:, :)
      real(dp), allocatable :: aplus(:, :)
      real(dp), allocatable :: product(:, :)
   contains
      procedure :: destroy => hessian_destroy
   end type response_hessian_t

   integer, parameter :: DEFAULT_MAX_ITER = 200  !! CG iterations before giving up

   integer, parameter :: HESSIAN_CHUNK = 64
      !! Columns of the Hessian built per pass over the integrals. At 115
      !! orbitals, 64 of them is 60 MB of densities.

   integer(int64), parameter :: SOLVE_BATCH_BYTES = 8_int64*1024_int64**3
      !! What the concurrent frequency solves may take, in bytes.

   real(dp), parameter :: MO_TRANSFORM_LIMIT = 16.0e9_dp
      !! What the exact MO transform may take at its peak, in bytes. Not
      !! comparable to `IN_CORE_LIMIT`: that one bounds a tensor kept for a
      !! whole solve, this one a few arrays held across one transform.

   real(dp), parameter :: DENSE_OPERATOR_LIMIT = 8.0e9_dp
      !! Above this many bytes the operator is never formed at all, whatever
      !! route would fill it, and `dynamic_response_iterative` takes over.
      !! `(A+B)`, `(A-B)` and their product are three `n_ov^2` matrices, and
      !! `n_ov` is the *product* of the occupied and virtual counts.

   integer, parameter :: IN_CORE_MAX_ORBITALS = 40
      !! Above this many orbitals, recompute the integrals rather than store
      !! them: storing wins below about forty and loses badly above it.

   real(dp), parameter, public :: IN_CORE_LIMIT = 4.0e9_dp
      !! Even below that, refuse to store a tensor larger than this many bytes.

   real(dp), parameter :: DEFAULT_TOL = 1.0e-11_dp
      !! Convergence on the residual norm, relative to the right-hand side.

   ! The matrix-free frequency-dependent solve takes `DEFAULT_DYNAMIC_TOL` and
   ! `DEFAULT_DYNAMIC_MAXITER` from `mqc_calculation_defaults` instead, because
   ! `keywords.efp.dynamic_tolerance` and `keywords.efp.dynamic_maxiter` name
   ! them from a deck.

   real(dp) :: prof_dens = 0.0_dp, prof_fock = 0.0_dp, prof_back = 0.0_dp
      !! Where a matrix-free pass spends itself -- forming the atomic-orbital
      !! densities, contracting them, transforming back -- accumulated across a
      !! whole solve and reported once at the end. Not threadsafe:
      !! `response_batch` is called from serial code and threads inside the
      !! build.
   integer :: prof_calls = 0  !! Passes accumulated into the three above

   integer, parameter, public :: N_CASIMIR_POLDER = 12
      !! Casimir-Polder quadrature points: how many imaginary frequencies a
      !! potential is tabulated at.

contains

   subroutine hessian_destroy(self)
      !! Release a built Hessian
      class(response_hessian_t), intent(inout) :: self
      if (allocated(self%aminus)) deallocate (self%aminus)
      if (allocated(self%aplus)) deallocate (self%aplus)
      if (allocated(self%product)) deallocate (self%product)
      self%ready = .false.
      ! `fitted` describes the matrices just released, so it has to go with
      ! them -- left set, it would claim the next Hessian built into this
      ! object came from fitted integrals whatever it was built from.
      self%fitted = .false.
   end subroutine hessian_destroy

   subroutine cphf_solve(mol, orbitals, orbital_energies, n_occ, perturbations, &
                         response, error, max_iter, tol, iterations, in_core, mo_rhs, &
                         eri_in, bmat, xc, density)
      !! Solve the coupled-perturbed equations for one or more perturbations
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)          !! MO coefficients, (n_ao, n_mo)
      real(dp), intent(in) :: orbital_energies(:)     !! (n_mo), Hartree
      integer, intent(in) :: n_occ
      real(dp), intent(in), optional :: perturbations(:, :, :)
         !! One-electron perturbation operators in the **AO** basis,
         !! `(n_ao, n_ao, n_perturbations)`. Exactly one of this and `mo_rhs`
         !! is given.
      real(dp), allocatable, intent(out) :: response(:, :, :)
         !! `U_ai` per perturbation, `(n_vir, n_occ, n_perturbations)`.
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: iterations
         !! Worst iteration count over the perturbations, for a caller that
         !! wants to notice the equations getting harder.
      logical, intent(in), optional :: in_core
         !! Store every integral instead of recomputing them. Default is direct.
      real(dp), intent(in), optional, target :: eri_in(:, :, :, :)
         !! An already-built two-electron tensor, used in place of building one.
         !! Implies `in_core`, and is pointed at rather than copied.
      type(xc_context_t), intent(inout), optional :: xc
         !! Solve the coupled-perturbed *Kohn-Sham* equations rather than the
         !! Hartree-Fock ones: the operator gains the exchange-correlation
         !! kernel. Give it whenever the reference was Kohn-Sham -- omitting it
         !! solves a different linear system and answers about a different
         !! energy.
      real(dp), intent(in), optional :: density(:, :)
         !! The converged reference density, where the kernel is evaluated.
      real(dp), intent(in), optional :: bmat(:, :)
         !! The fitted tensor `B(mu nu, P)`, in place of any four-index
         !! integrals. Not a storage choice like `in_core`: it makes the
         !! operator the fitted reference's own, so give it when and only when
         !! the SCF whose response this is was itself fitted. Mutually exclusive
         !! with `eri_in`.
      real(dp), intent(in), optional :: mo_rhs(:, :, :)
         !! The right-hand side already in the occupied-virtual MO block,
         !! `(n_vir, n_occ, n_perturbations)`. The Z-vector equation of an MP2
         !! gradient arrives this way, its right-hand side being a Lagrangian
         !! with no AO matrix it is the transform of.

      real(dp), allocatable, target :: eri_own(:, :, :, :)
      real(dp), pointer :: eri(:, :, :, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :), bounds(:, :)
      real(dp), allocatable :: gaps(:, :), rhs(:, :), x(:, :), r(:, :), z(:, :)
      real(dp), allocatable :: p(:, :), ap(:, :), work(:, :), zero_h(:, :)
      real(dp) :: rz, rz_new, pap, target_norm, step, use_tol
      integer :: n_ao, n_mo, n_vir, n_pert, ipert, a, i, iter, worst, limit
      logical :: direct
      character(len=32) :: text

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_vir = n_mo - n_occ

      if (n_occ < 1 .or. n_occ >= n_mo) then
         write (text, "(i0,a,i0)") n_occ, " of ", n_mo
         call error%set(ERROR_VALIDATION, "CPHF needs both occupied and virtual "// &
                        "orbitals; got "//trim(text))
         return
      end if
      if (size(orbital_energies) /= n_mo) then
         call error%set(ERROR_VALIDATION, "CPHF: one orbital energy per orbital")
         return
      end if
      if (present(perturbations) .eqv. present(mo_rhs)) then
         call error%set(ERROR_VALIDATION, "CPHF: give either AO perturbations or an "// &
                        "MO-basis right-hand side, not both and not neither")
         return
      end if
      if (present(perturbations)) then
         if (size(perturbations, 1) /= n_ao .or. size(perturbations, 2) /= n_ao) then
            call error%set(ERROR_VALIDATION, "CPHF: the perturbations are not n_ao square")
            return
         end if
         n_pert = size(perturbations, 3)
      else
         if (size(mo_rhs, 1) /= n_vir .or. size(mo_rhs, 2) /= n_occ) then
            call error%set(ERROR_VALIDATION, "CPHF: the right-hand side is not "// &
                           "(n_vir, n_occ)")
            return
         end if
         n_pert = size(mo_rhs, 3)
      end if
      limit = DEFAULT_MAX_ITER
      if (present(max_iter)) limit = max_iter
      use_tol = DEFAULT_TOL
      if (present(tol)) use_tol = tol

      allocate (c_occ(n_ao, n_occ), c_vir(n_ao, n_vir))
      c_occ = orbitals(:, 1:n_occ)
      c_vir = orbitals(:, n_occ + 1:n_mo)

      ! The orbital energy differences, which are both the diagonal of the
      ! response operator and the CG preconditioner.
      allocate (gaps(n_vir, n_occ))
      do i = 1, n_occ
         do a = 1, n_vir
            gaps(a, i) = orbital_energies(n_occ + a) - orbital_energies(i)
         end do
      end do
      if (minval(gaps) <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "CPHF: an occupied orbital lies above a "// &
                        "virtual one, so the reference is not a minimum and the "// &
                        "response operator is not positive definite")
         return
      end if

      ! Integral-direct by default: the stored tensor is n_ao^4, so it is the
      ! storage and not the solver that bounds how large a fragment can be done.
      direct = .true.
      if (present(in_core)) direct = .not. in_core
      if (present(eri_in)) direct = .false.
      if (present(bmat)) then
         ! Neither stored nor direct: fitting removes the four-index integrals
         ! rather than choosing where to keep them, so both of the others are
         ! left empty and `response_operator` takes a third path.
         if (present(eri_in)) then
            call error%set(ERROR_VALIDATION, "CPHF: given both a fitted tensor and an "// &
                           "exact one. They are different operators, not two routes to "// &
                           "the same one; pass whichever the reference was built with.")
            return
         end if
         direct = .false.
         allocate (eri_own(0, 0, 0, 0))
         eri => eri_own
         allocate (bounds(0, 0))
      else if (present(eri_in)) then
         eri => eri_in
         allocate (bounds(0, 0))
      else if (direct) then
         call schwarz_bounds(mol, bounds, error)
         if (error%has_error()) return
         allocate (eri_own(0, 0, 0, 0))
         eri => eri_own
      else
         call mol%eris(eri_own)
         eri => eri_own
         allocate (bounds(0, 0))
      end if
      allocate (zero_h(n_ao, n_ao))
      zero_h = 0.0_dp

      allocate (response(n_vir, n_occ, n_pert))
      allocate (rhs(n_vir, n_occ), x(n_vir, n_occ), r(n_vir, n_occ))
      allocate (z(n_vir, n_occ), p(n_vir, n_occ), ap(n_vir, n_occ))
      allocate (work(n_ao, n_occ))
      worst = 0

      do ipert = 1, n_pert
         ! h_ai, the perturbation in the occupied-virtual block -- or that
         ! block itself, when the caller assembled it there.
         if (present(perturbations)) then
            call pic_gemm(perturbations(:, :, ipert), c_occ, work)
            call pic_gemm(c_vir, work, rhs, transa="T")
         else
            rhs = mo_rhs(:, :, ipert)
         end if
         rhs = -rhs

         target_norm = sqrt(sum(rhs*rhs))
         if (target_norm <= 0.0_dp) then
            response(:, :, ipert) = 0.0_dp
            cycle
         end if
         target_norm = target_norm*use_tol

         ! Preconditioned CG from zero rather than from the uncoupled solution:
         ! with x = 0 the first iterate *is* the uncoupled answer, and the
         ! residual norm then starts at ||h|| so the threshold means what it says.
         x = 0.0_dp
         r = rhs
         z = r/gaps
         p = z
         rz = sum(r*z)

         iter = 0
         do iter = 1, limit
            if (present(bmat) .and. present(xc)) then
               call response_operator(mol, direct, eri, bounds, zero_h, c_occ, c_vir, &
                                      gaps, p, ap, error, bmat=bmat, xc=xc, &
                                      density=density)
            else if (present(bmat)) then
               call response_operator(mol, direct, eri, bounds, zero_h, c_occ, c_vir, &
                                      gaps, p, ap, error, bmat=bmat)
            else if (present(xc)) then
               call response_operator(mol, direct, eri, bounds, zero_h, c_occ, c_vir, &
                                      gaps, p, ap, error, xc=xc, density=density)
            else
               call response_operator(mol, direct, eri, bounds, zero_h, c_occ, c_vir, &
                                      gaps, p, ap, error)
            end if
            if (error%has_error()) return
            pap = sum(p*ap)
            if (pap <= 0.0_dp) then
               call error%set(ERROR_VALIDATION, "CPHF: the response operator is not "// &
                              "positive definite along a search direction, which means "// &
                              "the reference SCF is a saddle point")
               return
            end if
            step = rz/pap
            x = x + step*p
            r = r - step*ap
            if (sqrt(sum(r*r)) <= target_norm) exit
            z = r/gaps
            rz_new = sum(r*z)
            p = z + (rz_new/rz)*p
            rz = rz_new
         end do

         if (sqrt(sum(r*r)) > target_norm) then
            write (text, "(i0)") iter - 1
            call error%set(ERROR_GENERIC, "CPHF did not converge in "// &
                           trim(text)//" iterations")
            return
         end if

         response(:, :, ipert) = x
         worst = max(worst, iter)
      end do

      if (present(iterations)) iterations = worst
      nullify (eri)
      if (allocated(eri_own)) deallocate (eri_own)
      deallocate (bounds, c_occ, c_vir, gaps, rhs, x, r, z, p, ap, work, zero_h)
   end subroutine cphf_solve

   subroutine response_operator(mol, direct, eri, bounds, zero_h, c_occ, c_vir, &
                                gaps, u, au, error, bmat, xc, density)
      !! Apply the coupled-perturbed operator to a trial rotation
      type(czt_molecule_t), intent(in) :: mol
      logical, intent(in) :: direct          !! Recompute integrals rather than store them
      real(dp), intent(in) :: eri(:, :, :, :)   !! Zero-sized when `direct`
      real(dp), intent(in) :: bounds(:, :)      !! Zero-sized when not `direct`
      real(dp), intent(in) :: zero_h(:, :)   !! Zero, so the build returns only G
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :)
      real(dp), intent(in) :: gaps(:, :)
      real(dp), intent(in) :: u(:, :)
      real(dp), intent(out) :: au(:, :)
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional :: xc
         !! Present, the exchange-correlation kernel is added on top of the
         !! two-electron part, whichever integral source built it.
      real(dp), intent(in), optional :: density(:, :)
         !! The reference density the kernel is evaluated at. Required with `xc`.
      real(dp), intent(in), optional :: bmat(:, :)
         !! The fitted tensor `B(mu nu, P)`. Present, the operator is built from
         !! it and neither `eri` nor the direct build is touched.

      real(dp), allocatable :: dtilde(:, :), g(:, :), half(:, :), work(:, :)
      type(direct_stats_t) :: stats
      integer :: n_ao, n_occ
      real(dp) :: kf

      ! The exchange fraction the *reference* kept. Building full exchange for a
      ! pure functional makes the operator indefinite, and the solver then
      ! reports a saddle point that is not there.
      kf = 1.0_dp
      if (present(xc)) kf = xc%exx_fraction

      n_ao = size(c_occ, 1)
      n_occ = size(c_occ, 2)
      allocate (dtilde(n_ao, n_ao), g(n_ao, n_ao), half(n_ao, n_occ), work(n_ao, n_occ))

      ! Dt = C_vir U C_occ^T + transpose. Symmetrized because `build_fock`
      ! contracts a symmetric density, and because the physical first-order
      ! density is symmetric.
      call pic_gemm(c_vir, u, half)
      call pic_gemm(half, c_occ, dtilde, transb="T")
      dtilde = dtilde + transpose(dtilde)

      if (present(bmat)) then
         ! `half` is C_vir U, which is exactly the factor the fitted build
         ! wants: it never assembles Dt at all.
         call response_operator_df(bmat, half, c_occ, dtilde, g, k_scale=kf)
      else if (direct) then
         ! density_screen=.false. is not optional here. `dtilde` is the trial
         ! rotation the CG solver drives towards zero, so a density-weighted
         ! screen would tighten as the solve converges and the operator would
         ! stop being the fixed linear map the solver assumes. The plain Schwarz
         ! screen depends on the basis alone and keeps it constant.
         call build_fock_direct(mol, zero_h, dtilde, bounds, g, stats, error, &
                                density_screen=.false., k_scale=kf)
         if (error%has_error()) return
      else
         call build_fock(zero_h, eri, dtilde, g, k_scale=kf)
      end if

      ! The kernel, on top of whatever built `g` above.
      if (present(xc)) then
         if (.not. present(density)) then
            call error%set(ERROR_VALIDATION, "the response operator was given an "// &
                           "exchange-correlation context but no reference density to "// &
                           "evaluate its kernel at")
            return
         end if
         call xc_kernel_apply(xc, mol, density, dtilde, g, error)
         if (error%has_error()) return
      end if

      call pic_gemm(g, c_occ, work)
      call pic_gemm(c_vir, work, au, transa="T")
      au = gaps*u + 2.0_dp*au

      deallocate (dtilde, g, half, work)
   end subroutine response_operator

   subroutine response_operator_df(b, x, c_occ, dtilde, g, k_scale)
      !! `J - K/2` for a response density, from the fitted tensor
      !!
      !! Not `build_fock_df`: that one assumes an idempotent SCF density to get
      !! its occupied orbitals from, and a response density is symmetric but
      !! *indefinite*, so it has none. It arrives already factored instead,
      !!
      !!     Dt = X C_occ^T + C_occ X^T,     X = C_vir U
      !!
      !! which turns exchange into
      !!
      !!     K = sum_P [ (B_P X)(B_P C_occ)^T + (B_P C_occ)(B_P X)^T ]
      !!
      !! two n^2 n_occ products per auxiliary function. `Dt` is still passed,
      !! but only for the Coulomb term, where any density will do.
      real(dp), intent(in) :: b(:, :)        !! `B(mu nu, P)`, (n_ao^2, naux)
      real(dp), intent(in) :: x(:, :)        !! `C_vir U`, (n_ao, n_occ)
      real(dp), intent(in) :: c_occ(:, :)    !! (n_ao, n_occ)
      real(dp), intent(in) :: dtilde(:, :)   !! The assembled response density, for J
      real(dp), intent(out) :: g(:, :)
      real(dp), intent(in), optional :: k_scale
         !! The exchange fraction the reference kept. Absent is all of it.

      real(dp), allocatable :: coul(:, :), exch(:, :), bx(:, :), bc(:, :)
      real(dp) :: c_p, kf
      integer :: n, n_occ, naux, p

      kf = 1.0_dp
      if (present(k_scale)) kf = k_scale

      n = size(c_occ, 1)
      n_occ = size(c_occ, 2)
      naux = size(b, 2)
      allocate (coul(n, n), exch(n, n), bx(n, n_occ), bc(n, n_occ))

      coul = 0.0_dp
      exch = 0.0_dp
      do p = 1, naux
         associate (b_p => reshape(b(:, p), [n, n]))
            c_p = sum(b_p*dtilde)
            coul = coul + c_p*b_p
            call pic_gemm(b_p, x, bx)
            call pic_gemm(b_p, c_occ, bc)
            call pic_gemm(bx, bc, exch, transb="T", alpha=1.0_dp, beta=1.0_dp)
            call pic_gemm(bc, bx, exch, transb="T", alpha=1.0_dp, beta=1.0_dp)
         end associate
      end do

      g = coul - 0.5_dp*kf*exch
   end subroutine response_operator_df

   subroutine fitted_potential_general(b, dens, g, k_scale)
      !! `J - k K/2` from the fitted tensor, for a density with no structure
      !!
      !! The third of the three fitted builds, and the one that assumes nothing:
      !! `build_fock_df` needs an idempotent density and `response_operator_df`
      !! one that arrives already factored. An MP2 relaxed density is neither.
      !! Exchange goes through the general form, per auxiliary function,
      !!
      !!     J += (sum_uv B^P_uv D_uv) B^P,     K += B^P D B^P
      !!
      !! at `n^3 n_aux` rather than the `n^2 n_occ n_aux` the other two manage.
      real(dp), intent(in) :: b(:, :)      !! `B(mu nu, P)`, (n_ao^2, naux)
      real(dp), intent(in) :: dens(:, :)   !! Symmetric, otherwise arbitrary
      real(dp), intent(out) :: g(:, :)
      real(dp), intent(in), optional :: k_scale
         !! How much exact exchange the reference kept. Absent is all of it.

      real(dp), allocatable :: coul(:, :), exch(:, :), bd(:, :)
      real(dp) :: c_p, kf
      integer :: n, naux, p

      kf = 1.0_dp
      if (present(k_scale)) kf = k_scale

      n = size(dens, 1)
      naux = size(b, 2)
      allocate (coul(n, n), exch(n, n), bd(n, n))
      coul = 0.0_dp
      exch = 0.0_dp
      do p = 1, naux
         associate (b_p => reshape(b(:, p), [n, n]))
            c_p = sum(b_p*dens)
            coul = coul + c_p*b_p
            call pic_gemm(b_p, dens, bd)
            call pic_gemm(bd, b_p, exch, alpha=1.0_dp, beta=1.0_dp)
         end associate
      end do

      g = coul - 0.5_dp*kf*exch
   end subroutine fitted_potential_general

   subroutine response_operator_minus(mol, direct, eri, bounds, zero_h, c_occ, c_vir, &
                                      gaps, u, au, error)
      !! Apply `A - B`, the other half of a frequency-dependent response
      !!
      !! Same shape as `response_operator` and the same Fock build, on the
      !! *antisymmetrized* response density instead of the symmetrized one:
      !!
      !!     A + B   from   Dt = C_vir U C_occ^T + transpose
      !!     A - B   from   Dt = C_vir U C_occ^T - transpose
      !!
      !! The Coulomb term vanishes identically on the second, so `A - B` carries
      !! only exchange -- `delta (eps_a - eps_i) - (ab|ij) + (aj|ib)`.
      !!
      !! **Must go through `build_fock_direct_nosym`, not the fast build.** That
      !! one folds three of the eightfold permutations into a multiplicity
      !! factor, which doubles them where an antisymmetric density needs them to
      !! cancel.
      type(czt_molecule_t), intent(in) :: mol
      logical, intent(in) :: direct              !! Recompute integrals rather than store them
      real(dp), intent(in) :: eri(:, :, :, :)    !! Zero-sized when `direct`
      real(dp), intent(in) :: bounds(:, :)
      real(dp), intent(in) :: zero_h(:, :)
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :)
      real(dp), intent(in) :: gaps(:, :)
      real(dp), intent(in) :: u(:, :)
      real(dp), intent(out) :: au(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dtilde(:, :, :), g(:, :, :), half(:, :), work(:, :)
      type(direct_stats_t) :: stats
      integer :: n_ao, n_occ

      n_ao = size(c_occ, 1)
      n_occ = size(c_occ, 2)
      allocate (dtilde(n_ao, n_ao, 1), half(n_ao, n_occ), work(n_ao, n_occ))

      call pic_gemm(c_vir, u, half)
      call pic_gemm(half, c_occ, dtilde(:, :, 1), transb="T")
      dtilde(:, :, 1) = dtilde(:, :, 1) - transpose(dtilde(:, :, 1))

      ! In core the plain four-index contraction is already right for an
      ! antisymmetric density -- the Coulomb term vanishes on its own, since the
      ! integral is symmetric in the contracted pair -- so no `nosym` variant is
      ! needed on this path.
      if (direct) then
         call build_fock_direct_nosym(mol, zero_h, dtilde, bounds, g, stats, error)
         if (error%has_error()) return
      else
         allocate (g(size(dtilde, 1), size(dtilde, 2), 1))
         call build_fock(zero_h, eri, dtilde(:, :, 1), g(:, :, 1))
      end if

      call pic_gemm(g(:, :, 1), c_occ, work)
      call pic_gemm(c_vir, work, au, transa="T")
      au = gaps*u + 2.0_dp*au

      deallocate (dtilde, g, half, work)
   end subroutine response_operator_minus

   function casimir_polder_frequencies() result(nu)
      !! The twelve imaginary frequencies a potential is tabulated at, a.u.
      !!
      !! Gauss-Legendre on `t` in (-1,1) under `nu = w0 (1 + t)/(1 - t)` with
      !! `w0 = 0.3`, the Casimir-Polder quadrature GAMESS uses. Runs from
      !! 0.002792 to 32.239080.
      real(dp) :: nu(N_CASIMIR_POLDER)
      real(dp), parameter :: W0 = 0.3_dp
      !! 12-point Gauss-Legendre abscissae on (-1,1), positive half mirrored.
      real(dp), parameter :: T(N_CASIMIR_POLDER) = [ &
                             -0.9815606342467192_dp, -0.9041172563704749_dp, &
                             -0.7699026741943047_dp, -0.5873179542866175_dp, &
                             -0.3678314989981802_dp, -0.1252334085114689_dp, &
                             0.1252334085114689_dp, 0.3678314989981802_dp, &
                             0.5873179542866175_dp, 0.7699026741943047_dp, &
                             0.9041172563704749_dp, 0.9815606342467192_dp]
      integer :: i
      do i = 1, N_CASIMIR_POLDER
         nu(i) = W0*(1.0_dp + T(i))/(1.0_dp - T(i))
      end do
   end function casimir_polder_frequencies

   subroutine dynamic_polarizability(mol, orbitals, orbital_energies, n_occ, &
                                     frequencies, alpha, error, max_iter, tol, &
                                     response, perturbations, in_core, hessian, &
                                     progress, aux, route, allow_unconverged, batch)
      !! `alpha(i nu)` at each imaginary frequency
      !!
      !! On the imaginary axis the time-dependent equations are real and positive
      !! definite, with no poles. Writing `S = X + Y` and `D' = (X - Y)/i`,
      !!
      !!     (A+B) S + nu D' = -2h
      !!     (A-B) D' - nu S = 0
      !!
      !! and eliminating `D'` gives
      !!
      !!     [ (A+B) + nu^2 (A-B)^-1 ] S = -2h
      !!
      !! which is better conditioned as `nu` grows rather than worse.
      !!
      !! `alpha_kl(i nu) = -2 sum_ai h^k_ai S^l_ai`, which at `nu = 0` reduces to
      !! the static `-4 sum h U` because `S = 2U` there, and so must reproduce
      !! `static_polarizability` exactly.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: frequencies(:)      !! Imaginary frequencies, a.u.
      real(dp), allocatable, intent(out) :: alpha(:, :, :)
         !! `(n_pert, n_pert, n_frequencies)`, Bohr^3 for the dipole case.
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      real(dp), allocatable, intent(out), optional :: response(:, :, :, :)
         !! `S_ai` per perturbation per frequency, `(n_vir, n_occ, n_pert, n_freq)`,
         !! for a caller that wants to distribute it over localized orbitals rather
         !! than only take the total.
      type(response_hessian_t), intent(inout), optional :: hessian
         !! Reuse a built Hessian, or take one back to reuse. It depends only on
         !! the reference, not on the perturbations or the frequencies, so the
         !! three blocks of a potential can share one build.
      logical, intent(in), optional :: in_core
         !! Store the two-electron integrals instead of recomputing them each
         !! time the response operator is applied. Defaults to storing whenever
         !! the tensor fits `IN_CORE_LIMIT`, since this solver applies that
         !! operator frequencies times perturbations times outer times inner
         !! iterations.
      real(dp), intent(in), optional :: perturbations(:, :, :)
         !! One-electron operators in the AO basis. Defaults to the three dipole
         !! components; the six quadrupole components instead make `alpha` the
         !! quadrupole-quadrupole response, which is the `LMOQQPOL` block of a
         !! potential.
      logical, intent(in), optional :: progress
         !! Report where the two long phases are while they run. This routine is
         !! most of the wall clock on any real fragment and otherwise prints
         !! nothing for minutes at a time.
      type(czt_molecule_t), intent(in), optional :: aux
         !! An auxiliary basis. Supplied, the Hessian is built from fitted
         !! integrals rather than from `n_ov` Fock builds -- see
         !! `build_hessian_df`, and `validation/check_df_hessian` for what the
         !! approximation costs. Absent, the build is exact.
      integer, intent(in), optional :: route
         !! `EFP_RESPONSE_AUTO`, `_DENSE` or `_MATRIX_FREE`, which is what
         !! `keywords.efp.response` carries. Absent means `auto`, which is the
         !! size rule below and nothing besides.
      logical, intent(in), optional :: allow_unconverged
         !! Passed to `dynamic_response_iterative`: return whatever was reached
         !! instead of raising. The answer is wrong.
      integer, intent(in), optional :: batch
         !! Densities per integral pass, `keywords.efp.response_batch`.

      real(dp), allocatable :: dip(:, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :), gaps(:, :), h(:, :, :)
      real(dp), allocatable :: eri0(:, :, :, :), work(:, :)
      logical :: direct
      logical :: reuse, iterate
      integer :: take
      real(dp), allocatable :: aplus(:, :), aminus(:, :), product(:, :), lu(:, :)
      real(dp), allocatable :: rhs_flat(:, :), h_flat(:, :), solution(:, :)
      integer, allocatable :: ipiv(:), infos(:)
      integer :: n_ov, j, info, n_freq, concurrent
      real(dp), allocatable :: operators(:, :, :)
      real(dp) :: nu, rz, rz_new, pap, step, target_norm, use_tol
      integer :: n_ao, n_mo, n_vir, k, l, i, a, ifreq, iter, limit, n_pert
      character(len=32) :: text
      logical :: talk
      type(timer_type) :: clock

      character(len=MAX_LINE_LENGTH) :: line

      talk = .false.
      if (present(progress)) talk = progress
      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_vir = n_mo - n_occ
      limit = DEFAULT_MAX_ITER
      if (present(max_iter)) limit = max_iter
      use_tol = DEFAULT_TOL
      if (present(tol)) use_tol = tol

      if (n_occ < 1 .or. n_occ >= n_mo) then
         call error%set(ERROR_VALIDATION, "dynamic response needs both occupied "// &
                        "and virtual orbitals")
         return
      end if

      if (present(perturbations)) then
         operators = perturbations
      else
         call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, error)
         if (error%has_error()) return
         operators = dip
      end if
      n_pert = size(operators, 3)
      if (size(operators, 1) /= n_ao .or. size(operators, 2) /= n_ao) then
         call error%set(ERROR_VALIDATION, "dynamic response: the perturbations are "// &
                        "not n_ao square")
         return
      end if
      direct = .true.
      if (n_ao <= IN_CORE_MAX_ORBITALS .and. &
          real(n_ao, dp)**4*8.0_dp <= IN_CORE_LIMIT) direct = .false.
      if (present(in_core)) direct = .not. in_core
      if (direct) then
         call schwarz_bounds(mol, bounds, error)
         if (error%has_error()) return
         allocate (eri0(0, 0, 0, 0))
      else
         call mol%eris(eri0)
         allocate (bounds(0, 0))
      end if

      allocate (c_occ(n_ao, n_occ), c_vir(n_ao, n_vir), zero_h(n_ao, n_ao))
      c_occ = orbitals(:, 1:n_occ)
      c_vir = orbitals(:, n_occ + 1:n_mo)
      zero_h = 0.0_dp

      allocate (gaps(n_vir, n_occ))
      do i = 1, n_occ
         do a = 1, n_vir
            gaps(a, i) = orbital_energies(n_occ + a) - orbital_energies(i)
         end do
      end do
      if (minval(gaps) <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "dynamic response: the reference is not "// &
                        "a minimum")
         return
      end if

      allocate (h(n_vir, n_occ, n_pert), work(n_ao, n_occ))
      do k = 1, n_pert
         call pic_gemm(operators(:, :, k), c_occ, work)
         call pic_gemm(c_vir, work, h(:, :, k), transa="T")
      end do

      allocate (alpha(n_pert, n_pert, size(frequencies)))
      if (present(response)) then
         allocate (response(n_vir, n_occ, n_pert, size(frequencies)))
      end if
      ! Built once, then solved densely. Materialising `(A+B)` and `(A-B)` costs
      ! `n_ov` Fock builds; after that every frequency and every perturbation is
      ! dense linear algebra on matrices already in hand. The product
      ! `(A-B)(A+B)` is the same at every frequency -- only the `nu^2` on the
      ! diagonal changes -- so it is formed once too, and each frequency is one
      ! LU factorisation and one solve per perturbation.
      n_ov = n_vir*n_occ
      reuse = .false.
      if (present(hessian)) reuse = hessian%ready
      take = EFP_RESPONSE_AUTO
      if (present(route)) take = route

      ! Too big to build, so never built. A fitted response is excluded because
      ! `build_hessian_df` is bounded by the auxiliary basis rather than by
      ! `n_ov`, and an already-built Hessian because the expensive part is
      ! behind it. `dense` builds the operator whatever the size rule says;
      ! `matrix_free` declines to build it even where it would fit, and declines
      ! the auxiliary basis and any Hessian on hand with it.
      iterate = .false.
      if (.not. reuse) then
         iterate = .not. present(aux) .and. &
                   3.0_dp*real(n_ov, dp)**2*8.0_dp > DENSE_OPERATOR_LIMIT
      end if
      select case (take)
      case (EFP_RESPONSE_AUTO)
      case (EFP_RESPONSE_DENSE)
         iterate = .false.
      case (EFP_RESPONSE_MATRIX_FREE)
         iterate = .true.
      case default
         call error%set(ERROR_VALIDATION, "dynamic response: unknown response route. "// &
                        "Accepted: auto, dense, matrix_free")
         return
      end select

      ! Said out loud only where a deck overruled the size rule. Which route ran
      ! is visible either way; that a keyword is the reason is not.
      if (talk .and. take /= EFP_RESPONSE_AUTO) then
         if (iterate) then
            call logger%info("        response route: matrix free, by keywords.efp.response")
            if (present(aux)) call logger%info("          the auxiliary basis fits a Hessian "// &
                                               "this route never builds, so it goes unused")
            if (reuse) call logger%info("          a built Hessian was on hand and goes unused too")
         else
            call logger%info("        response route: dense, by keywords.efp.response")
         end if
      end if

      if (iterate) then
         call dynamic_response_iterative(mol, direct, eri0, bounds, zero_h, c_occ, &
                                         c_vir, gaps, h, frequencies, alpha, error, &
                                         max_iter=max_iter, tol=tol, &
                                         response=response, progress=talk, &
                                         allow_unconverged=allow_unconverged, &
                                         batch=batch)
         return
      end if
      if (reuse) then
         if (size(hessian%product, 1) /= n_ov) then
            call error%set(ERROR_VALIDATION, "dynamic response: the supplied Hessian "// &
                           "is the wrong size for this reference")
            return
         end if
      end if

      ! Skipped entirely on reuse. The `if (reuse)` block above only checks the
      ! supplied Hessian's size; these builds then ran anyway and their results
      ! were discarded below in favour of `hessian%aplus`/`%aminus`, so passing
      ! a built Hessian across a potential's three blocks saved one GEMM and
      ! none of the `n_ov` Fock builds it exists to save.
      if (reuse) then
         continue
      else if (present(aux)) then
         call build_hessian_df(mol, aux, c_occ, c_vir, gaps, aplus, aminus, error, &
                               progress=talk)
         if (present(hessian)) hessian%fitted = .true.
      else if (mo_transform_fits(mol%nao, n_occ, n_vir, direct)) then
         ! Exact integrals, assembled the way the fitted build assembles them.
         ! The column build is kept for the systems whose AO tensor will not fit.
         call build_hessian_mo(mol, eri0, c_occ, c_vir, gaps, aplus, aminus, error, &
                               progress=talk)
      else
         call build_hessian(mol, direct, eri0, bounds, zero_h, c_occ, c_vir, gaps, &
                            aplus, aminus, HESSIAN_CHUNK, error, progress=talk)
      end if
      if (error%has_error()) return

      allocate (lu(n_ov, n_ov), ipiv(n_ov))
      allocate (rhs_flat(n_ov, n_pert), h_flat(n_ov, n_pert))
      if (.not. reuse) then
         allocate (product(n_ov, n_ov))
         call pic_gemm(aminus, aplus, product)
      end if

      do l = 1, n_pert
         h_flat(:, l) = reshape(h(:, :, l), [n_ov])
      end do
      ! The right-hand side is `-2 (A-B) h`, which the multiplication through by
      ! `(A-B)` leaves behind.
      if (reuse) then
         call pic_gemm(hessian%aminus, h_flat, rhs_flat)
      else
         call pic_gemm(aminus, h_flat, rhs_flat)
      end if
      rhs_flat = -2.0_dp*rhs_flat

      allocate (solution(n_ov, n_pert))

      ! The frequencies are independent -- each its own shifted factorization of
      ! the same matrix -- so they are threaded over rather than inside the
      ! factorization, which keeps the sequential BLAS sequential and needs no
      ! threaded LAPACK. How many at once is a memory question and nothing else:
      ! `getrf` factorizes in place, so each carries its own copy of the operator.
      n_freq = size(frequencies)
      concurrent = concurrent_solves(n_ov, n_freq)
      allocate (infos(n_freq))
      infos = 0
      if (talk) then
         write (line, "(A,I0,A,I0,A,I0,A,I0,A)") "        solving ", n_freq, &
            " frequencies x ", n_pert, " perturbations over ", n_ov, " pairs, ", &
            concurrent, " at a time"
         call logger%info(trim(line))
         call clock%start()
      end if
      deallocate (lu, ipiv, solution)

      !$omp parallel do default(none) num_threads(concurrent) schedule(dynamic) &
      !$omp    private(ifreq, lu, ipiv, solution, j, k, l, info) &
      !$omp    shared(n_freq, n_ov, n_pert, n_vir, n_occ, reuse, hessian, product, &
      !$omp           aplus, frequencies, rhs_flat, h_flat, alpha, response, infos)
      do ifreq = 1, n_freq
         allocate (lu(n_ov, n_ov), ipiv(n_ov), solution(n_ov, n_pert))
         if (frequencies(ifreq) == 0.0_dp) then
            ! The static limit, solved as the static equation `(A+B) S = -2h`
            ! rather than as the zero-frequency member of the family: going
            ! through the product would square the condition number.
            if (reuse) then
               lu = hessian%aplus
            else
               lu = aplus
            end if
            solution = -2.0_dp*h_flat
         else
            if (reuse) then
               lu = hessian%product
            else
               lu = product
            end if
            do j = 1, n_ov
               lu(j, j) = lu(j, j) + frequencies(ifreq)**2
            end do
            solution = rhs_flat
         end if
         call pic_getrf(lu, ipiv, info)
         if (info /= 0) then
            infos(ifreq) = 1
         else
            call pic_getrs(lu, ipiv, solution, info=info)
            if (info /= 0) infos(ifreq) = 2
         end if

         if (infos(ifreq) == 0) then
            do l = 1, n_pert
               do k = 1, n_pert
                  alpha(k, l, ifreq) = -2.0_dp*sum(h_flat(:, k)*solution(:, l))
               end do
               if (present(response)) then
                  response(:, :, l, ifreq) = reshape(solution(:, l), [n_vir, n_occ])
               end if
            end do
         end if
         deallocate (lu, ipiv, solution)
      end do
      !$omp end parallel do

      ! Read after the region rather than raised inside it: `error_t` carries a
      ! message, and threads setting one would write over each other.
      if (any(infos == 1)) then
         call error%set(ERROR_GENERIC, "dynamic response: the operator is singular")
         return
      end if
      if (any(infos == 2)) then
         call error%set(ERROR_GENERIC, "dynamic response: the solve failed")
         return
      end if
      deallocate (infos)
      if (talk) call tick(clock, n_freq, n_freq, "frequency")

      ! Hand the build back if the caller wants it, rather than freeing it.
      if (present(hessian) .and. .not. reuse) then
         call move_alloc(aminus, hessian%aminus)
         call move_alloc(aplus, hessian%aplus)
         call move_alloc(product, hessian%product)
         hessian%ready = .true.
      else if (.not. reuse) then
         deallocate (aplus, aminus, product)
      end if
      deallocate (rhs_flat, h_flat)

      deallocate (bounds, zero_h, c_occ, c_vir, gaps, h, work, eri0, operators)
      if (allocated(dip)) deallocate (dip)
   end subroutine dynamic_polarizability

   subroutine build_hessian(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                            aplus, aminus, chunk, error, progress)
      !! `(A+B)` and `(A-B)` as explicit matrices over the occupied-virtual space
      !!
      !! Column `j` is the operator applied to the unit vector `e_j`, so the
      !! whole thing is `n_ov` applications, batched so one pass over the
      !! integrals serves `chunk` columns. Every frequency and perturbation
      !! afterwards is dense linear algebra on matrices already in hand.
      !!
      !! `n_ov` is `n_vir * n_occ`, and it is its square that limits this rather
      !! than the basis: a thousand occupied-virtual pairs needs 16 MB for the
      !! two matrices, ten thousand needs 1.6 GB.
      type(czt_molecule_t), intent(in) :: mol
      logical, intent(in) :: direct
      real(dp), intent(in) :: eri(:, :, :, :)
      real(dp), intent(in) :: bounds(:, :), zero_h(:, :)
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), allocatable, intent(out) :: aplus(:, :), aminus(:, :)
      integer, intent(in) :: chunk
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: progress

      real(dp), allocatable :: u(:, :, :), au(:, :, :)
      integer, allocatable :: idx(:)
      integer :: n_vir, n_occ, n_ov, first, last, nthis, m, col, a, i
      logical :: talk
      type(timer_type) :: clock

      character(len=MAX_LINE_LENGTH) :: line

      talk = .false.
      if (present(progress)) talk = progress
      n_vir = size(gaps, 1)
      n_occ = size(gaps, 2)
      n_ov = n_vir*n_occ
      allocate (aplus(n_ov, n_ov), aminus(n_ov, n_ov))
      allocate (u(n_vir, n_occ, chunk), au(n_vir, n_occ, chunk), idx(chunk))
      if (talk) then
         write (line, "(A,I0,A,I0,A,F0.1,A)") "        Hessian: ", n_ov, " columns in chunks of ", chunk, &
            ", ", 2.0_dp*real(n_ov, dp)**2*8.0_dp/1.0e6_dp, " MB"
         call logger%info(trim(line))
         call clock%start()
      end if

      first = 1
      do while (first <= n_ov)
         last = min(first + chunk - 1, n_ov)
         nthis = last - first + 1
         u = 0.0_dp
         do m = 1, nthis
            col = first + m - 1
            a = mod(col - 1, n_vir) + 1
            i = (col - 1)/n_vir + 1
            u(a, i, m) = 1.0_dp
            idx(m) = m
         end do

         call response_batch(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                             u, idx, nthis, .false., au, error)
         if (error%has_error()) return
         do m = 1, nthis
            aplus(:, first + m - 1) = reshape(au(:, :, m), [n_ov])
         end do

         call response_batch(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                             u, idx, nthis, .true., au, error)
         if (error%has_error()) return
         do m = 1, nthis
            aminus(:, first + m - 1) = reshape(au(:, :, m), [n_ov])
         end do

         if (talk) call tick(clock, last, n_ov, "column")
         first = last + 1
      end do
      deallocate (u, au, idx)
   end subroutine build_hessian

   subroutine build_hessian_df(orb, aux, c_occ, c_vir, gaps, aplus, aminus, error, &
                               progress)
      !! `(A+B)` and `(A-B)` from fitted integrals instead of from Fock builds
      !!
      !! The operators are, written out,
      !!
      !!     (A+B)_{ai,bj} = d Deps + 4(ai|bj) - (ab|ij) - (aj|ib)
      !!     (A-B)_{ai,bj} = d Deps          - (ab|ij) + (aj|ib)
      !!
      !! so with the integrals fitted every term is a matrix product over the
      !! auxiliary index and `build_hessian`'s column loop disappears. `(ai|bj)`
      !! is `B_ov B_ov^T`, already in the `n_ov` layout the solver uses because
      !! `build_df_mo_block` runs its compound index left-fastest; `(aj|ib)` is
      !! that same matrix read at permuted indices; `(ab|ij)` is `B_vv B_oo^T`.
      !! Two gemms of `n_ov^2 naux` replace `n_ov` Fock builds.
      !!
      !! **The result is not the same matrix** -- fitting is an approximation,
      !! and `validation/check_df_hessian` measures it against the exact build.
      !!
      !! Memory is what limits this: `n_ov^2` for each operator, and `n_ov^2`
      !! again for each of the two integral matrices while they are assembled.
      type(czt_molecule_t), intent(in) :: orb, aux
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), allocatable, intent(out) :: aplus(:, :), aminus(:, :)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: progress

      ! TODO(mqc): `gap_flat`, `term_abij`, `term_ajib`, `row`, `col`, `a`, `i`,
      ! `b` and `j` are never read here -- leftovers from before
      ! `assemble_hessian` was factored out of this routine and `build_hessian_mo`.
      real(dp), allocatable :: bov(:, :), bvv(:, :), boo(:, :)
      real(dp), allocatable :: coul(:, :), exch(:, :), gap_flat(:)
      real(dp) :: term_abij, term_ajib
      integer :: n_vir, n_occ, n_ov, naux, row, col, a, i, b, j
      logical :: talk
      type(timer_type) :: clock

      character(len=MAX_LINE_LENGTH) :: line

      talk = .false.
      if (present(progress)) talk = progress
      n_vir = size(gaps, 1)
      n_occ = size(gaps, 2)
      n_ov = n_vir*n_occ
      if (talk) call clock%start()

      call build_df_mo_block(orb, aux, c_vir, c_occ, bov, error)
      if (error%has_error()) return
      call build_df_mo_block(orb, aux, c_vir, c_vir, bvv, error)
      if (error%has_error()) return
      call build_df_mo_block(orb, aux, c_occ, c_occ, boo, error)
      if (error%has_error()) return
      naux = size(bov, 2)
      if (talk) then
         write (line, "(A,I0,A,I0,A,F0.1,A)") "        fitted Hessian: ", n_ov, " pairs, ", naux, &
            " auxiliary functions, ", 4.0_dp*real(n_ov, dp)**2*8.0_dp/1.0e6_dp, " MB"
         call logger%info(trim(line))
         call tick(clock, 1, 3, "step")
      end if

      allocate (coul(n_ov, n_ov))
      call pic_gemm(bov, bov, coul, transb="T")
      allocate (exch(n_vir*n_vir, n_occ*n_occ))
      call pic_gemm(bvv, boo, exch, transb="T")
      deallocate (bov, bvv, boo)
      if (talk) call tick(clock, 2, 3, "step")

      call assemble_hessian(coul, exch, gaps, aplus, aminus)
      deallocate (coul, exch)
      if (talk) call tick(clock, 3, 3, "step")
   end subroutine build_hessian_df

   function concurrent_solves(n_ov, n_freq) result(concurrent)
      !! How many frequency solves may run at once
      !!
      !! One, when the BLAS threads itself -- see `MQC_SEQUENTIAL_BLAS` in the
      !! top-level CMakeLists. Otherwise a memory question: `getrf` factorizes in
      !! place, so every concurrent solve needs its own `n_ov^2` copy of the
      !! operator.
      use omp_lib, only: omp_get_max_threads
      integer, intent(in) :: n_ov
      integer, intent(in) :: n_freq
      integer :: concurrent

      integer(int64) :: per_solve

#ifndef MQC_SEQUENTIAL_BLAS
      ! A threaded BLAS already fills the machine inside each factorization, and
      ! nesting a parallel loop around it takes concurrency away rather than
      ! adding it -- MKL serialises a call made from a parallel region.
      concurrent = 1
      return
#endif
      per_solve = int(n_ov, int64)**2*8_int64
      concurrent = int(max(1_int64, SOLVE_BATCH_BYTES/max(per_solve, 1_int64)))
      concurrent = min(concurrent, n_freq, omp_get_max_threads())
   end function concurrent_solves

   function mo_transform_fits(n_ao, n_occ, n_vir, direct) result(fits)
      !! Can the whole Hessian be had by transformation rather than column by column
      !!
      !! There is no `n_ao^4` term: `ket_transformed_pairs` contracts the ket as
      !! the quartets come out, so the largest thing held is the half-transformed
      !! pair block. Weighed at the peak, the stage where the pair blocks and the
      !! first bra transform are both live. `(A+B)` and `(A-B)` are left out
      !! because the column build allocates those too.
      integer, intent(in) :: n_ao
      integer, intent(in) :: n_occ
      integer, intent(in) :: n_vir
      logical, intent(in) :: direct
      logical :: fits

      real(dp) :: right, peak

      if (.not. direct) then
         fits = .true.
         return
      end if
      ! Both ket blocks, `(b j)` and `(i j)`, share the pass and so are both
      ! live, as are the two operators and the product the solve goes on to form.
      right = real(n_vir, dp)*real(n_occ, dp) + real(n_occ, dp)**2
      peak = real(n_ao, dp)**2*right + real(n_vir, dp)*real(n_ao, dp)*right &
             + 4.0_dp*(real(n_vir, dp)*real(n_occ, dp))**2 &
             + (real(n_vir, dp)*real(n_occ, dp))**2
      fits = peak*8.0_dp <= MO_TRANSFORM_LIMIT
   end function mo_transform_fits

   subroutine build_hessian_mo(mol, eri_in, c_occ, c_vir, gaps, aplus, aminus, error, &
                               progress)
      !! `(A+B)` and `(A-B)` from exact integrals, transformed rather than probed
      !!
      !! What is wanted is two blocks of MO integrals, and those come from one
      !! pass over the AO integrals followed by dense linear algebra --
      !! `build_hessian_df`'s argument applied to exact integrals instead of
      !! fitted ones, sharing `assemble_hessian` so the only difference between
      !! the two is where the integrals came from. `(aj|ib)` needs no block of
      !! its own: it is `coul` read at permuted indices.
      !!
      !! Four quarter-transforms, none forming anything bigger than
      !! `n_ao^3 * n_occ`, where the column build costs `n_quartets * n_ov`
      !! however its columns are chunked.
      !!
      !! Memory is what bounds it, and it is the AO tensor rather than anything
      !! here: `n_ao^4`, checked against `IN_CORE_LIMIT` by the caller.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in), target :: eri_in(:, :, :, :)
         !! Zero-sized to have this build its own with one integral pass
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), allocatable, intent(out) :: aplus(:, :), aminus(:, :)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: progress

      real(dp), allocatable :: half(:, :), pair_ov(:, :), pair_oo(:, :)
      real(dp), allocatable :: step_ov(:, :), step_oo(:, :)
      real(dp), allocatable :: coul(:, :), exch(:, :)
      integer :: n_ao, n_vir, n_occ, n_ov, j
      logical :: talk, direct_ket
      type(timer_type) :: clock

      character(len=MAX_LINE_LENGTH) :: line

      talk = .false.
      if (present(progress)) talk = progress
      n_ao = size(c_occ, 1)
      n_vir = size(gaps, 1)
      n_occ = size(gaps, 2)
      n_ov = n_vir*n_occ
      if (talk) call clock%start()

      allocate (pair_ov(n_ao**2, n_ov), pair_oo(n_ao**2, n_occ*n_occ))
      direct_ket = size(eri_in) == 0
      if (talk) then
         write (line, "(A,I0,A,F0.1,A)") "        exact Hessian: ", n_ov, &
            " pairs by transformation, ", &
            (real(n_ao, dp)**2*real(n_ov, dp) + 4.0_dp*real(n_ov, dp)**2)*8.0_dp/1.0e6_dp, " MB"
         call logger%info(trim(line))
         call tick(clock, 1, 3, "step")
      end if

      ! One and two: the ket pair onto the MO blocks. From a stored tensor this
      ! is two gemms over its `n_ao^3` leading extent; without one it is a pass
      ! over the quartets that contracts as they come out, so the `n_ao^4` never
      ! exists.
      if (direct_ket) then
         call ket_transformed_pairs(mol, c_vir, c_occ, pair_ov, pair_oo)
      else
         allocate (half(n_ao**3, n_occ))
         call pic_gemm(reshape(eri_in, [n_ao**3, n_ao]), c_occ, half)
         do j = 1, n_occ
            call pic_gemm(reshape(half(:, j), [n_ao**2, n_ao]), c_vir, &
                          pair_ov(:, 1 + (j - 1)*n_vir:j*n_vir))
            call pic_gemm(reshape(half(:, j), [n_ao**2, n_ao]), c_occ, &
                          pair_oo(:, 1 + (j - 1)*n_occ:j*n_occ))
         end do
         deallocate (half)
      end if
      if (talk) call tick(clock, 2, 3, "step")

      ! Three and four: the bra pair, the same way round. `step_*` holds the first
      ! AO index transformed, with the second still in the AO basis.
      allocate (step_ov(n_vir, n_ao*n_ov), step_oo(n_vir, n_ao*n_occ*n_occ))
      call pic_gemm(c_vir, reshape(pair_ov, [n_ao, n_ao*n_ov]), step_ov, transa="T")
      call pic_gemm(c_vir, reshape(pair_oo, [n_ao, n_ao*n_occ*n_occ]), step_oo, transa="T")
      deallocate (pair_ov, pair_oo)

      allocate (coul(n_ov, n_ov), exch(n_vir*n_vir, n_occ*n_occ))
      call transform_second(step_ov, c_occ, n_vir, n_ao, n_ov, coul)
      call transform_second(step_oo, c_vir, n_vir, n_ao, n_occ*n_occ, exch)
      deallocate (step_ov, step_oo)

      call assemble_hessian(coul, exch, gaps, aplus, aminus)
      deallocate (coul, exch)
      if (talk) call tick(clock, 3, 3, "step")
   end subroutine build_hessian_mo

   subroutine transform_second(step, c, n_vir, n_ao, n_right, out)
      !! Transform the second AO index of `(a nu | right)` with `c`
      !!
      !! `step` is `(n_vir, n_ao * n_right)` with `nu` running fastest inside the
      !! second extent, so for each column of the right-hand pair this is a small
      !! gemm of `(n_vir, n_ao)` against `c`. The result's compound index runs
      !! left-fastest -- `a` before `i` -- which is the layout `assemble_hessian`
      !! and the solver both read.
      real(dp), intent(in) :: step(:, :)
      real(dp), intent(in) :: c(:, :)
      integer, intent(in) :: n_vir, n_ao, n_right
      real(dp), intent(out) :: out(:, :)

      real(dp), allocatable :: slab(:, :)
      integer :: r, n_left2

      n_left2 = size(c, 2)
      !$omp parallel default(none) private(r, slab) &
      !$omp    shared(step, c, out, n_vir, n_ao, n_right, n_left2)
      allocate (slab(n_vir, n_left2))
      !$omp do schedule(static)
      do r = 1, n_right
         call pic_gemm(step(:, 1 + (r - 1)*n_ao:r*n_ao), c, slab)
         out(:, r) = reshape(slab, [n_vir*n_left2])
      end do
      !$omp end do
      deallocate (slab)
      !$omp end parallel
   end subroutine transform_second

   subroutine assemble_hessian(coul, exch, gaps, aplus, aminus)
      !! `(A+B)` and `(A-B)` from the two MO integral blocks they are made of
      !!
      !!     (A+B)_{ai,bj} = d Deps + 4(ai|bj) - (ab|ij) - (aj|ib)
      !!     (A-B)_{ai,bj} = d Deps          - (ab|ij) + (aj|ib)
      !!
      !! Shared by the fitted build and the exact one, which differ only in where
      !! `coul` and `exch` came from, so whatever `validation/check_df_hessian`
      !! measures between them is the fitting error and nothing else.
      real(dp), intent(in) :: coul(:, :)   !! `(ai|bj)`, indexed `[a+(i-1)nv, b+(j-1)nv]`
      real(dp), intent(in) :: exch(:, :)   !! `(ab|ij)`, indexed `[a+(b-1)nv, i+(j-1)no]`
      real(dp), intent(in) :: gaps(:, :)
      real(dp), allocatable, intent(out) :: aplus(:, :), aminus(:, :)

      real(dp), allocatable :: gap_flat(:)
      real(dp) :: term_abij, term_ajib
      integer :: n_vir, n_occ, n_ov, row, col, a, i, b, j

      n_vir = size(gaps, 1)
      n_occ = size(gaps, 2)
      n_ov = n_vir*n_occ

      allocate (aplus(n_ov, n_ov), aminus(n_ov, n_ov))
      !$omp parallel do default(none) private(col, row, a, i, b, j, term_abij, term_ajib) &
      !$omp shared(n_ov, n_vir, n_occ, coul, exch, aplus, aminus) schedule(static)
      do col = 1, n_ov
         b = mod(col - 1, n_vir) + 1
         j = (col - 1)/n_vir + 1
         do row = 1, n_ov
            a = mod(row - 1, n_vir) + 1
            i = (row - 1)/n_vir + 1
            term_abij = exch(a + (b - 1)*n_vir, i + (j - 1)*n_occ)
            ! `(aj|ib)`: the pair `(a,j)` against the pair `(b,i)`, both of which are
            ! virtual-occupied, so both are indices into the same matrix.
            term_ajib = coul(a + (j - 1)*n_vir, b + (i - 1)*n_vir)
            aplus(row, col) = 4.0_dp*coul(row, col) - term_abij - term_ajib
            aminus(row, col) = term_ajib - term_abij
         end do
      end do
      !$omp end parallel do

      gap_flat = reshape(gaps, [n_ov])
      do col = 1, n_ov
         aplus(col, col) = aplus(col, col) + gap_flat(col)
         aminus(col, col) = aminus(col, col) + gap_flat(col)
      end do
      deallocate (gap_flat)
   end subroutine assemble_hessian

   subroutine tick(clock, done, total, unit_name)
      !! One progress line: how far in, how long it has taken, what is left
      !!
      !! The time left assumes the remaining work costs what the finished work
      !! did, which screening makes only approximately true.
      type(timer_type), intent(inout) :: clock
      integer, intent(in) :: done, total
      character(len=*), intent(in) :: unit_name

      real(dp) :: so_far, left

      ! Read the clock without stopping it: `start` resets, so a stop/start pair
      ! here would leave `so_far` measuring only the interval since the previous
      ! line.
      character(len=MAX_LINE_LENGTH) :: line

      so_far = clock%get_elapsed_time()
      left = so_far*real(total - done, dp)/real(max(done, 1), dp)
      write (line, "(A,I6,A,I0,1X,A,A,F8.1,A,F8.1,A)") "        ", done, " of ", total, trim(unit_name), &
         "s   ", so_far, " s in, ~", left, " s left"
      call logger%info(trim(line))
      flush (output_unit)
   end subroutine tick

   subroutine dynamic_response_iterative(mol, direct, eri, bounds, zero_h, c_occ, &
                                         c_vir, gaps, h, frequencies, alpha, error, &
                                         max_iter, tol, response, progress, &
                                         allow_unconverged, batch)
      !! The frequency-dependent response without ever forming its operator
      !!
      !! For the sizes where `2 n_ov^2` of storage and an `n_ov^3` factorisation
      !! per frequency are out of reach. Each iteration is two passes over the
      !! integrals however many pairs there are, and the storage is a handful of
      !! vectors per system.
      !!
      !! **Every system iterates together, which is the whole economy.** The
      !! integrals a direct build recomputes do not depend on which density is
      !! being contracted, so one pass serves every frequency and perturbation at
      !! once, down to the index mask that drops a system out of the pass when it
      !! converges.
      !!
      !! **Bi-conjugate gradient stabilised, because the operator is not
      !! symmetric.** Eliminating `D'` from the coupled pair leaves
      !!
      !!     [ (A-B)(A+B) + nu^2 ] S = -2 (A-B) h
      !!
      !! and a product of two symmetric matrices is not symmetric. The
      !! alternative form `[(A+B) + nu^2 (A-B)^-1] S = -2h` is symmetric and
      !! would take CG, but each application needs an inner solve whose iteration
      !! count differs per system, which is what stops the batch sharing a pass.
      !! BiCGSTAB keeps the batch in lockstep at two applications per iteration.
      !!
      !! **Zero frequency is a different equation and is solved as one**,
      !! `(A+B) S = -2h`; going through the product would square the condition
      !! number. Those systems skip the `(A-B)` application, which the index mask
      !! expresses without a second code path.
      !!
      !! Density screening was tried and loses: a trial vector here is
      !! `C_vir U C_occ^T`, delocalised even when the density that generated it
      !! is not, so the test finds nothing negligible and only costs.
      use pic_blas_interfaces, only: pic_gemm
      type(czt_molecule_t), intent(in) :: mol
      logical, intent(in) :: direct
      real(dp), intent(in) :: eri(:, :, :, :)
      real(dp), intent(in) :: bounds(:, :), zero_h(:, :)
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), intent(in) :: h(:, :, :)              !! (n_vir, n_occ, n_pert)
      real(dp), intent(in) :: frequencies(:)
      real(dp), allocatable, intent(out) :: alpha(:, :, :)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      real(dp), allocatable, intent(out), optional :: response(:, :, :, :)
      logical, intent(in), optional :: progress
      logical, intent(in), optional :: allow_unconverged
         !! Return whatever was reached instead of raising. See
         !! `efp_config_t%allow_crap_response` -- the answer is wrong.
      integer, intent(in), optional :: batch
         !! Densities per integral pass. The right value is a property of the
         !! machine and the basis, not of the physics -- see
         !! `efp_config_t%response_batch`.

      real(dp), allocatable :: x(:, :, :), r(:, :, :), r0(:, :, :), p(:, :, :)
      real(dp), allocatable :: v(:, :, :), s(:, :, :), t(:, :, :), rhs(:, :, :)
      real(dp), allocatable :: ph(:, :, :), sh(:, :, :), precon(:, :, :)
      real(dp), allocatable :: rho(:), rho_old(:), omega(:), alpha_bi(:), beta(:)
      real(dp), allocatable :: nu2(:), rnorm(:), bnorm(:)
      integer, allocatable :: pert_of(:), freq_of(:), live(:), nonzero(:)
      logical, allocatable :: done(:)
      integer :: n_vir, n_occ, n_pert, n_freq, n_sys, m, k, l, it, cycles, eff_batch
      integer :: nlive, nnz
      real(dp) :: threshold, worst
      logical :: talk
      type(timer_type) :: clock
      character(len=MAX_LINE_LENGTH) :: line

      talk = .false.
      if (present(progress)) talk = progress
      n_vir = size(gaps, 1)
      n_occ = size(gaps, 2)
      n_pert = size(h, 3)
      n_freq = size(frequencies)
      n_sys = n_freq*n_pert
      cycles = DEFAULT_DYNAMIC_MAXITER
      if (present(max_iter)) cycles = max_iter
      threshold = DEFAULT_DYNAMIC_TOL
      if (present(tol)) threshold = tol

      allocate (x(n_vir, n_occ, n_sys), r(n_vir, n_occ, n_sys))
      allocate (r0(n_vir, n_occ, n_sys), p(n_vir, n_occ, n_sys))
      allocate (v(n_vir, n_occ, n_sys), s(n_vir, n_occ, n_sys))
      allocate (t(n_vir, n_occ, n_sys), rhs(n_vir, n_occ, n_sys))
      allocate (ph(n_vir, n_occ, n_sys), sh(n_vir, n_occ, n_sys))
      allocate (precon(n_vir, n_occ, n_sys))
      allocate (rho(n_sys), rho_old(n_sys), omega(n_sys), alpha_bi(n_sys), beta(n_sys))
      allocate (nu2(n_sys), rnorm(n_sys), bnorm(n_sys), done(n_sys))
      allocate (pert_of(n_sys), freq_of(n_sys), live(n_sys), nonzero(n_sys))

      m = 0
      do k = 1, n_freq
         do l = 1, n_pert
            m = m + 1
            freq_of(m) = k
            pert_of(m) = l
            nu2(m) = frequencies(k)**2
         end do
      end do

      ! Resolved once so the banner reports the same width the passes use.
      eff_batch = DEFAULT_RESPONSE_BATCH
      if (present(batch)) then
         if (batch > 0) eff_batch = batch
      end if

      ! The right-hand side. `-2 (A-B) h` at every nonzero frequency, `-2 h` at
      ! zero -- one batched application of `(A-B)` covers all of them, because
      ! the zero-frequency systems simply are not in the mask.
      do m = 1, n_sys
         rhs(:, :, m) = h(:, :, pert_of(m))
      end do
      nnz = 0
      do m = 1, n_sys
         if (nu2(m) /= 0.0_dp) then
            nnz = nnz + 1
            nonzero(nnz) = m
         end if
      end do
      prof_dens = 0.0_dp
      prof_fock = 0.0_dp
      prof_back = 0.0_dp
      prof_calls = 0
      if (talk) then
         write (line, "(A,I0,A,I0,A,I0,A)") "        solving ", n_freq, &
            " frequencies x ", n_pert, " perturbations over ", n_vir*n_occ, &
            " pairs, matrix free"
         call logger%info(trim(line))
         write (line, "(A,I0,A,I0,A,F0.2,A)") "          ", n_sys, &
            " systems in flight, ", 4*((n_sys + eff_batch - 1)/eff_batch), &
            " integral passes per iteration, ", &
            11.0_dp*real(n_vir*n_occ, dp)*real(n_sys, dp)*8.0_dp/1.0e9_dp, " GB of vectors"
         call logger%info(trim(line))
         write (line, "(A,I0,A)") "          building the right-hand side, ", &
            (nnz + 11)/12, " integral passes"
         call logger%info(trim(line))
         flush (output_unit)
         call clock%start()
      end if

      ! Chunked like every other pass: handing all the systems to one call is a
      ! density block and a per-thread `n_ao^2` accumulator of tens of gigabytes.
      if (nnz > 0) then
         call batched_in_chunks(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                                rhs, nonzero, nnz, .true., t, error, batch)
         if (error%has_error()) return
         do m = 1, nnz
            rhs(:, :, nonzero(m)) = t(:, :, nonzero(m))
         end do
      end if
      rhs = -2.0_dp*rhs

      ! Jacobi on the operator's diagonal: `(A-B)(A+B)` is dominated by the
      ! square of the orbital energy gap, and the shift adds to it. The static
      ! solve's preconditioner, squared.
      do m = 1, n_sys
         if (nu2(m) == 0.0_dp) then
            precon(:, :, m) = 1.0_dp/gaps
         else
            precon(:, :, m) = 1.0_dp/(gaps**2 + nu2(m))
         end if
      end do

      x = 0.0_dp
      r = rhs
      r0 = r
      p = 0.0_dp
      v = 0.0_dp
      rho_old = 1.0_dp
      alpha_bi = 1.0_dp
      omega = 1.0_dp
      done = .false.
      do m = 1, n_sys
         bnorm(m) = sqrt(sum(rhs(:, :, m)**2))
         if (bnorm(m) <= 0.0_dp) then
            bnorm(m) = 1.0_dp
            done(m) = .true.
         end if
      end do

      if (talk) then
         write (line, "(A,F0.1,A,ES8.1,A,I0,A)") "          right-hand side done in ", &
            clock%get_elapsed_time(), " s; iterating to ", threshold, &
            ", at most ", cycles, " iterations"
         call logger%info(trim(line))
         flush (output_unit)
         call clock%start()
      end if

      do it = 1, cycles
         nlive = 0
         do m = 1, n_sys
            if (.not. done(m)) then
               nlive = nlive + 1
               live(nlive) = m
            end if
         end do
         if (nlive == 0) exit

         do m = 1, nlive
            k = live(m)
            rho(k) = sum(r0(:, :, k)*r(:, :, k))
            if (it == 1) then
               p(:, :, k) = r(:, :, k)
            else
               beta(k) = (rho(k)/rho_old(k))*(alpha_bi(k)/omega(k))
               p(:, :, k) = r(:, :, k) + beta(k)*(p(:, :, k) - omega(k)*v(:, :, k))
            end if
            ph(:, :, k) = precon(:, :, k)*p(:, :, k)
         end do

         call apply_dynamic(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                            ph, live, nlive, nu2, v, error, batch)
         if (error%has_error()) return

         do m = 1, nlive
            k = live(m)
            worst = sum(r0(:, :, k)*v(:, :, k))
            if (worst == 0.0_dp) then
               done(k) = .true.
               cycle
            end if
            alpha_bi(k) = rho(k)/worst
            s(:, :, k) = r(:, :, k) - alpha_bi(k)*v(:, :, k)
            sh(:, :, k) = precon(:, :, k)*s(:, :, k)
         end do

         call apply_dynamic(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                            sh, live, nlive, nu2, t, error, batch)
         if (error%has_error()) return

         do m = 1, nlive
            k = live(m)
            if (done(k)) cycle
            worst = sum(t(:, :, k)**2)
            if (worst == 0.0_dp) then
               omega(k) = 0.0_dp
            else
               omega(k) = sum(t(:, :, k)*s(:, :, k))/worst
            end if
            x(:, :, k) = x(:, :, k) + alpha_bi(k)*ph(:, :, k) + omega(k)*sh(:, :, k)
            r(:, :, k) = s(:, :, k) - omega(k)*t(:, :, k)
            rnorm(k) = sqrt(sum(r(:, :, k)**2))/bnorm(k)
            if (rnorm(k) < threshold) done(k) = .true.
            if (omega(k) == 0.0_dp) done(k) = .true.
            rho_old(k) = rho(k)
         end do

         ! Every iteration, not every fifth: one of these is minutes at the
         ! sizes this route exists for, and what is worth watching is whether the
         ! residual is still falling.
         if (talk) then
            worst = 0.0_dp
            if (any(.not. done)) worst = maxval(rnorm, mask=.not. done)
            write (line, "(A,I4,A,F9.1,A,I4,A,ES9.2)") "          iteration ", it, &
               "   ", clock%get_elapsed_time(), " s in, ", count(.not. done), &
               " systems live, worst residual ", worst
            call logger%info(trim(line))
            flush (output_unit)
         end if
      end do

      if (any(.not. done)) then
         ! Loud even when it is allowed: nothing downstream of here can tell an
         ! unconverged polarizability from a converged one by looking.
         ! TODO(mqc): `goto 100` jumps out of the warning branch, which the house
         ! style forbids outright.
         if (present(allow_unconverged)) then
            if (allow_unconverged) then
               write (line, "(A,I0,A,ES9.2,A)") "          WARNING: ", &
                  count(.not. done), " systems did not converge, worst residual ", &
                  maxval(rnorm, mask=.not. done), " -- the potential is wrong"
               call logger%warning(trim(line))
               flush (output_unit)
               goto 100
            end if
         end if
         call error%set(ERROR_VALIDATION, "the frequency-dependent response did not "// &
                        "converge. The operator is positive definite when the "// &
                        "reference is a minimum, so a reference that is not one is "// &
                        "the first thing to check.")
         return
      end if
100   continue

      if (talk) then
         write (line, "(A,I0,A)") "          where the passes went (", prof_calls, &
            " of them, CPU seconds):"
         call logger%info(trim(line))
         write (line, "(A,F10.1,A,F10.1,A,F10.1)") &
            "            densities ", prof_dens, "   Fock ", prof_fock, &
            "   back-transform ", prof_back
         call logger%info(trim(line))
         flush (output_unit)
      end if

      allocate (alpha(n_pert, n_pert, n_freq))
      if (present(response)) allocate (response(n_vir, n_occ, n_pert, n_freq))
      do m = 1, n_sys
         k = freq_of(m)
         l = pert_of(m)
         do it = 1, n_pert
            alpha(it, l, k) = -2.0_dp*sum(h(:, :, it)*x(:, :, m))
         end do
         if (present(response)) response(:, :, l, k) = x(:, :, m)
      end do
   end subroutine dynamic_response_iterative

   subroutine apply_dynamic(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                            u, live, nlive, nu2, au, error, width)
      !! `[(A-B)(A+B) + nu^2] u`, or `(A+B) u` where the frequency is zero
      !!
      !! Two passes over the integrals for the whole batch, which is the point.
      !! The second serves only the systems at nonzero frequency; the rest have
      !! their answer after the first, because at `nu = 0` the equation never had
      !! the `(A-B)` factor in it.
      type(czt_molecule_t), intent(in) :: mol
      logical, intent(in) :: direct
      real(dp), intent(in) :: eri(:, :, :, :)
      real(dp), intent(in) :: bounds(:, :), zero_h(:, :)
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), intent(in) :: u(:, :, :)
      integer, intent(in) :: live(:)
      integer, intent(in) :: nlive
      real(dp), intent(in) :: nu2(:)
      real(dp), intent(inout) :: au(:, :, :)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: width

      real(dp), allocatable :: q1(:, :, :)
      integer, allocatable :: shifted(:)
      integer :: m, k, nshift

      if (nlive <= 0) return
      allocate (q1(size(u, 1), size(u, 2), size(u, 3)))
      q1 = 0.0_dp

      call batched_in_chunks(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                             u, live, nlive, .false., q1, error, width)
      if (error%has_error()) return

      allocate (shifted(nlive))
      nshift = 0
      do m = 1, nlive
         k = live(m)
         if (nu2(k) /= 0.0_dp) then
            nshift = nshift + 1
            shifted(nshift) = k
         else
            au(:, :, k) = q1(:, :, k)
         end if
      end do

      if (nshift > 0) then
         call batched_in_chunks(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                                q1, shifted, nshift, .true., au, error, width)
         if (error%has_error()) return
         do m = 1, nshift
            k = shifted(m)
            au(:, :, k) = au(:, :, k) + nu2(k)*u(:, :, k)
         end do
      end if
   end subroutine apply_dynamic

   subroutine batched_in_chunks(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                                u, idx, nact, minus, au, error, width)
      !! `response_batch`, a dozen densities at a time rather than all of them
      !!
      !! **More is not better past about twelve, and it is worse.** The batched
      !! Fock build shares one pass over the quartets across every density it is
      !! given, but each thread scatters into an accumulator sized `n_ao^2` by the
      !! batch width, which grows with the batch while the integral pass it is
      !! amortising does not. The win saturates near fourfold at a dozen
      !! densities and reverses by twenty-four.
      type(czt_molecule_t), intent(in) :: mol
      logical, intent(in) :: direct
      real(dp), intent(in) :: eri(:, :, :, :)
      real(dp), intent(in) :: bounds(:, :), zero_h(:, :)
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), intent(in) :: u(:, :, :)
      integer, intent(in) :: idx(:)
      integer, intent(in) :: nact
      logical, intent(in) :: minus
      real(dp), intent(inout) :: au(:, :, :)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: width
         !! Densities per integral pass. Absent or non-positive takes the
         !! built-in default; `keywords.efp.response_batch` overrides it.

      integer :: max_batch

      integer :: first, last

      max_batch = DEFAULT_RESPONSE_BATCH
      if (present(width)) then
         if (width > 0) max_batch = width
      end if
      first = 1
      do while (first <= nact)
         last = min(first + max_batch - 1, nact)
         call response_batch(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                             u, idx(first:last), last - first + 1, minus, au, error)
         if (error%has_error()) return
         first = last + 1
      end do
   end subroutine batched_in_chunks

   subroutine response_batch(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                             u, idx, nact, minus, au, error)
      !! `(A+B)u` or `(A-B)u` for many vectors in **one** integral pass
      !!
      !! Every right-hand side of the dynamic response shares one operator and
      !! differs only in the frequency shift, so a direct build recomputes the
      !! same integrals for all of them and one pass over the quartets serves the
      !! whole batch.
      !!
      !! `idx(1:nact)` selects which vectors are still wanted, so a converged
      !! system stops costing anything rather than riding along.
      !!
      !! `minus` picks the antisymmetric combination and the build that does not
      !! fold permutations, which is what `A - B` needs; otherwise the symmetric
      !! one.
      type(czt_molecule_t), intent(in) :: mol
      logical, intent(in) :: direct
      real(dp), intent(in) :: eri(:, :, :, :)
      real(dp), intent(in) :: bounds(:, :), zero_h(:, :)
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), intent(in) :: u(:, :, :)
      integer, intent(in) :: idx(:)
      integer, intent(in) :: nact
      logical, intent(in) :: minus
      real(dp), intent(inout) :: au(:, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dens(:, :, :), g(:, :, :), half(:, :), work(:, :)
      real(dp) :: t0, t1
      type(direct_stats_t) :: stats
      integer :: n_ao, n_occ, m, j

      if (nact <= 0) return
      n_ao = size(c_occ, 1)
      n_occ = size(c_occ, 2)
      allocate (dens(n_ao, n_ao, nact), half(n_ao, n_occ), work(n_ao, n_occ))
      call cpu_time(t0)

      do m = 1, nact
         j = idx(m)
         call pic_gemm(c_vir, u(:, :, j), half)
         call pic_gemm(half, c_occ, dens(:, :, m), transb="T")
         if (minus) then
            dens(:, :, m) = dens(:, :, m) - transpose(dens(:, :, m))
         else
            dens(:, :, m) = dens(:, :, m) + transpose(dens(:, :, m))
         end if
      end do

      call cpu_time(t1)
      prof_dens = prof_dens + (t1 - t0)
      t0 = t1

      if (direct) then
         if (minus) then
            call build_fock_direct_nosym(mol, zero_h, dens, bounds, g, stats, error)
         else
            call build_fock_direct_many(mol, zero_h, dens, bounds, g, stats, error)
         end if
         if (error%has_error()) return
      else
         allocate (g(n_ao, n_ao, nact))
         do m = 1, nact
            call build_fock(zero_h, eri, dens(:, :, m), g(:, :, m))
         end do
      end if

      call cpu_time(t1)
      prof_fock = prof_fock + (t1 - t0)
      t0 = t1
      do m = 1, nact
         j = idx(m)
         call pic_gemm(g(:, :, m), c_occ, work)
         call pic_gemm(c_vir, work, au(:, :, j), transa="T")
         au(:, :, j) = gaps*u(:, :, j) + 2.0_dp*au(:, :, j)
      end do

      call cpu_time(t1)
      prof_back = prof_back + (t1 - t0)
      prof_calls = prof_calls + 1
      deallocate (dens, g, half, work)
   end subroutine response_batch

   subroutine distributed_dynamic_cross(mol, orbitals, orbital_energies, n_occ, &
                                        frequencies, measure, respond, tensors, &
                                        centroids, error, n_core, max_iter, tol, hessian, &
                                        progress, aux, route, allow_unconverged, batch)
      !! Mixed-multipole dynamic response, per localized orbital and frequency
      !!
      !!     alpha^(i)_{km}(i nu) = -2 sum_a h^{measure,k}_{ai} S^{respond,m}_{ai}
      !!
      !! `measure` is the operator the response is read off with, `respond` the
      !! one that drives it. All three dynamic blocks of a potential are cases of
      !! this: dipole against dipole gives `DYNAMIC POLARIZABLE POINTS`, dipole
      !! against quadrupole the dipole-quadrupole block, quadrupole against
      !! quadrupole `LMOQQPOL`.
      !!
      !! **The quadrupole to pass is the traceless Buckingham form**, not the raw
      !! second moment, and the quadrupole-quadrupole case carries an extra factor
      !! of one third that the mixed case does not. Both are applied by the
      !! caller; see `mqc_efp_potential`.
      !!
      !! **Per orbital, which operator measures is not a free choice.** One order
      !! gives `h P M^-1 h` and the other `h M^-1 P h`, and the projector onto the
      !! localized set does not commute with the response operator, so the two
      !! differ. Summed over every orbital `P` is the identity and they agree,
      !! which is also why each per-orbital tensor is asymmetric and only the
      !! total is symmetric.
      !!
      !! Phase conventions do not enter: each component is quadratic in its
      !! localized orbital, so flipping an orbital's sign leaves it unchanged.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: frequencies(:)
      real(dp), intent(in) :: measure(:, :, :)   !! AO operators, (n_ao, n_ao, n_measure)
      real(dp), intent(in) :: respond(:, :, :)   !! AO operators, (n_ao, n_ao, n_respond)
      real(dp), allocatable, intent(out) :: tensors(:, :, :, :)
         !! `(n_measure, n_respond, n_localized, n_frequencies)`
      real(dp), allocatable, intent(out) :: centroids(:, :)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_core, max_iter
      real(dp), intent(in), optional :: tol
      type(response_hessian_t), intent(inout), optional :: hessian
         !! Passed straight through, so the blocks of one potential share a build.
      logical, intent(in), optional :: progress
         !! Passed straight through as well: the solver is where the time goes.
      type(czt_molecule_t), intent(in), optional :: aux
         !! Passed straight through: fit the Hessian rather than build it exactly.
      integer, intent(in), optional :: route
         !! Passed straight through: whether the response operator is built at
         !! all. `keywords.efp.response`.
      logical, intent(in), optional :: allow_unconverged
         !! Take an unconverged response rather than refusing it.
      integer, intent(in), optional :: batch
         !! Passed straight through: densities per integral pass.

      real(dp), allocatable :: alpha(:, :, :), s_all(:, :, :, :), localized(:, :)
      real(dp), allocatable :: s(:, :), sc(:, :), w(:, :), h_loc(:, :, :)
      real(dp), allocatable :: s_loc(:, :), c_occ(:, :), c_vir(:, :), work(:, :)
      integer :: n_ao, n_mo, n_vir, n_lmo, core, k, m, i, a, ifreq

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_vir = n_mo - n_occ
      core = 0
      if (present(n_core)) core = n_core
      n_lmo = n_occ - core
      if (core < 0 .or. n_lmo < 1) then
         call error%set(ERROR_VALIDATION, "no orbitals left to distribute over")
         return
      end if

      call boys_localize(mol, orbitals(:, core + 1:n_occ), n_lmo, localized, &
                         centroids, error)
      if (error%has_error()) return

      call dynamic_polarizability(mol, orbitals, orbital_energies, n_occ, frequencies, &
                                  alpha, error, max_iter=max_iter, tol=tol, &
                                  response=s_all, perturbations=respond, &
                                  hessian=hessian, progress=progress, aux=aux, &
                                  route=route, allow_unconverged=allow_unconverged, batch=batch)
      if (error%has_error()) return

      allocate (c_occ(n_ao, n_occ), c_vir(n_ao, n_vir))
      c_occ = orbitals(:, 1:n_occ)
      c_vir = orbitals(:, n_occ + 1:n_mo)
      call mol%overlap(s)
      allocate (sc(n_ao, n_lmo), w(n_occ, n_lmo))
      call pic_gemm(s, localized, sc)
      call pic_gemm(c_occ, sc, w, transa="T")

      allocate (h_loc(n_vir, n_lmo, size(measure, 3)), work(n_ao, n_lmo))
      do k = 1, size(measure, 3)
         call pic_gemm(measure(:, :, k), localized, work)
         call pic_gemm(c_vir, work, h_loc(:, :, k), transa="T")
      end do

      allocate (s_loc(n_vir, n_lmo))
      allocate (tensors(size(measure, 3), size(respond, 3), n_lmo, size(frequencies)))
      do ifreq = 1, size(frequencies)
         do m = 1, size(respond, 3)
            call pic_gemm(s_all(:, :, m, ifreq), w, s_loc)
            do i = 1, n_lmo
               do k = 1, size(measure, 3)
                  tensors(k, m, i, ifreq) = 0.0_dp
                  do a = 1, n_vir
                     tensors(k, m, i, ifreq) = tensors(k, m, i, ifreq) &
                                               - 2.0_dp*h_loc(a, i, k)*s_loc(a, i)
                  end do
               end do
            end do
         end do
      end do

      deallocate (alpha, s_all, localized, s, sc, w, h_loc, s_loc, c_occ, c_vir, work)
   end subroutine distributed_dynamic_cross

   subroutine distributed_dynamic_polarizability(mol, orbitals, orbital_energies, &
                                                 n_occ, frequencies, tensors, &
                                                 centroids, error, n_core, max_iter, &
                                                 tol, hessian, progress, aux)
      !! One tensor per localized orbital at each imaginary frequency
      !!
      !! What `DYNAMIC POLARIZABLE POINTS` carries: the Casimir-Polder integral
      !! over these is the `C6` between two fragments. The projection is
      !! `distributed_polarizability`'s -- insert `W W^T` into the occupied sum --
      !! applied to `S` at each frequency instead of to the static `U`, so the
      !! `nu = 0` member must reproduce the static tensors orbital by orbital.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: frequencies(:)
      real(dp), allocatable, intent(out) :: tensors(:, :, :, :)
         !! `(3, 3, n_localized, n_frequencies)`, Bohr^3.
      real(dp), allocatable, intent(out) :: centroids(:, :)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_core
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      type(response_hessian_t), intent(inout), optional :: hessian
         !! Passed straight through, so the blocks of one potential share a build.
      logical, intent(in), optional :: progress
         !! Passed straight through as well: the solver is where the time goes.
      type(czt_molecule_t), intent(in), optional :: aux
         !! Passed straight through: fit the Hessian rather than build it exactly.

      real(dp), allocatable :: alpha(:, :, :), s_all(:, :, :, :), localized(:, :)
      real(dp), allocatable :: dip(:, :, :), s(:, :), sc(:, :), w(:, :)
      real(dp), allocatable :: h_loc(:, :, :), s_loc(:, :), c_occ(:, :), c_vir(:, :)
      real(dp), allocatable :: work(:, :)
      integer :: n_ao, n_mo, n_vir, n_lmo, core, k, l, i, a, ifreq

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_vir = n_mo - n_occ
      core = 0
      if (present(n_core)) core = n_core
      n_lmo = n_occ - core
      if (core < 0 .or. n_lmo < 1) then
         call error%set(ERROR_VALIDATION, "no orbitals left to distribute over")
         return
      end if

      call boys_localize(mol, orbitals(:, core + 1:n_occ), n_lmo, localized, &
                         centroids, error)
      if (error%has_error()) return

      call dynamic_polarizability(mol, orbitals, orbital_energies, n_occ, frequencies, &
                                  alpha, error, max_iter=max_iter, tol=tol, &
                                  response=s_all, hessian=hessian, &
                                  progress=progress, aux=aux)
      if (error%has_error()) return

      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, error)
      if (error%has_error()) return

      allocate (c_occ(n_ao, n_occ), c_vir(n_ao, n_vir))
      c_occ = orbitals(:, 1:n_occ)
      c_vir = orbitals(:, n_occ + 1:n_mo)
      call mol%overlap(s)
      allocate (sc(n_ao, n_lmo), w(n_occ, n_lmo))
      call pic_gemm(s, localized, sc)
      call pic_gemm(c_occ, sc, w, transa="T")

      allocate (h_loc(n_vir, n_lmo, 3), work(n_ao, n_lmo), s_loc(n_vir, n_lmo))
      do k = 1, 3
         call pic_gemm(dip(:, :, k), localized, work)
         call pic_gemm(c_vir, work, h_loc(:, :, k), transa="T")
      end do

      allocate (tensors(3, 3, n_lmo, size(frequencies)))
      do ifreq = 1, size(frequencies)
         do l = 1, 3
            call pic_gemm(s_all(:, :, l, ifreq), w, s_loc)
            do i = 1, n_lmo
               do k = 1, 3
                  tensors(k, l, i, ifreq) = 0.0_dp
                  do a = 1, n_vir
                     tensors(k, l, i, ifreq) = tensors(k, l, i, ifreq) &
                                               - 2.0_dp*h_loc(a, i, k)*s_loc(a, i)
                  end do
               end do
            end do
         end do
      end do

      deallocate (alpha, s_all, localized, dip, s, sc, w, h_loc, s_loc)
      deallocate (c_occ, c_vir, work)
   end subroutine distributed_dynamic_polarizability

   subroutine static_response_dense(aplus, dip, c_occ, c_vir, u, error)
      !! `U` from `(A+B) U = -h`, factorized rather than iterated
      !!
      !! The same equation `cphf_solve` runs conjugate gradients on, solved
      !! directly because `(A+B)` is already in hand; `test_mqc_czt_cphf`
      !! pins the two together.
      !!
      !! Not the zero-frequency member of the dynamic family. That one is reached
      !! by multiplying through by `(A-B)`, and its `nu = 0` limit does *not*
      !! reproduce this to better than 4e-5.
      real(dp), intent(in) :: aplus(:, :)
      real(dp), intent(in) :: dip(:, :, :)      !! AO perturbations, `(n_ao, n_ao, n_pert)`
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :)
      real(dp), allocatable, intent(out) :: u(:, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: lu(:, :), rhs(:, :), work(:, :), block_h(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n_ao, n_occ, n_vir, n_ov, n_pert, k, info

      n_ao = size(c_occ, 1)
      n_occ = size(c_occ, 2)
      n_vir = size(c_vir, 2)
      n_ov = n_vir*n_occ
      n_pert = size(dip, 3)

      allocate (u(n_vir, n_occ, n_pert))
      allocate (lu(n_ov, n_ov), ipiv(n_ov), rhs(n_ov, n_pert))
      allocate (work(n_ao, n_occ), block_h(n_vir, n_occ))

      do k = 1, n_pert
         call pic_gemm(dip(:, :, k), c_occ, work)
         call pic_gemm(c_vir, work, block_h, transa="T")
         rhs(:, k) = -reshape(block_h, [n_ov])
      end do

      lu = aplus
      call pic_getrf(lu, ipiv, info)
      if (info /= 0) then
         call error%set(ERROR_GENERIC, "static response: (A+B) is singular")
         return
      end if
      call pic_getrs(lu, ipiv, rhs, info=info)
      if (info /= 0) then
         call error%set(ERROR_GENERIC, "static response: the solve failed")
         return
      end if

      do k = 1, n_pert
         u(:, :, k) = reshape(rhs(:, k), [n_vir, n_occ])
      end do
      deallocate (lu, ipiv, rhs, work, block_h)
   end subroutine static_response_dense

   subroutine distributed_polarizability(mol, orbitals, orbital_energies, n_occ, &
                                         tensors, centroids, error, n_core, &
                                         max_iter, tol, iterations, in_core, hessian)
      !! One polarizability tensor per localized orbital, at its centroid
      !!
      !! Where in the molecule the response happens, so an induced dipole can be
      !! placed on each bond and lone pair and polarized by the local field.
      !!
      !! **The decomposition is exact.** The occupied index of
      !! `alpha_kl = -4 sum_ai h^k_ai U^l_ai` is summed over, and a rotation `W`
      !! among the occupied orbitals is orthogonal, so inserting `W W^T` changes
      !! nothing:
      !!
      !!     alpha_kl = -4 sum_a sum_i' (h^k W)_ai' (U^l W)_ai'
      !!
      !! and each `i'` is one localized orbital's share. Nothing is fitted or
      !! partitioned by distance.
      !!
      !! **Each tensor is asymmetric, and that is not a bug.** The contraction
      !! pairs perturbation `k` with response `l`, and for a single orbital those
      !! are different objects; only the sum over orbitals restores symmetry. All
      !! nine components are returned and none is averaged away.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: tensors(:, :, :)
         !! `(3, 3, n_localized)`, Bohr^3, in the order the localized orbitals
         !! come back -- which is the order `centroids` uses too.
      real(dp), allocatable, intent(out) :: centroids(:, :)
         !! `(3, n_localized)`, Bohr. Where each tensor belongs.
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_core
         !! Occupied orbitals to leave out, lowest first. Default zero, where
         !! GAMESS excludes them, so a reference file wants its count.
         !!
         !! The response itself is always solved over *all* occupied orbitals,
         !! because the core does polarize; only the distribution drops it. What
         !! comes back therefore sums to the response projected onto the orbitals
         !! kept, exactly, for any `n_core`.
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: iterations
      logical, intent(in), optional :: in_core
      type(response_hessian_t), intent(in), optional :: hessian
         !! A built Hessian to solve against instead of iterating. The static
         !! response is `(A+B) U = -h`, so once the dynamic blocks have built
         !! `(A+B)` the same equation is one factorization away. A *fitted*
         !! Hessian is declined -- see `fitted` on the type.

      real(dp), allocatable :: dip(:, :, :), u(:, :, :), localized(:, :)
      real(dp), allocatable :: s(:, :), w(:, :), u_loc(:, :, :), h_loc(:, :, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :), work(:, :), sc(:, :)
      integer :: n_ao, n_mo, n_vir, n_lmo, core, k, l, i, a
      logical :: dense
      character(len=32) :: text

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_vir = n_mo - n_occ
      core = 0
      if (present(n_core)) core = n_core
      n_lmo = n_occ - core

      if (core < 0 .or. n_lmo < 1) then
         write (text, "(i0,a,i0)") core, " of ", n_occ
         call error%set(ERROR_VALIDATION, "cannot distribute a polarizability over "// &
                        "no orbitals: n_core is "//trim(text))
         return
      end if

      ! The localized orbitals to distribute over, and where they sit.
      call boys_localize(mol, orbitals(:, core + 1:n_occ), n_lmo, localized, &
                         centroids, error)
      if (error%has_error()) return

      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, error)
      if (error%has_error()) return

      allocate (c_occ(n_ao, n_occ), c_vir(n_ao, n_vir))
      c_occ = orbitals(:, 1:n_occ)
      c_vir = orbitals(:, n_occ + 1:n_mo)

      ! The response over the whole occupied space -- see `n_core` above.
      dense = .false.
      if (present(hessian)) dense = hessian%ready .and. .not. hessian%fitted
      if (dense) then
         call static_response_dense(hessian%aplus, dip, c_occ, c_vir, u, error)
         if (present(iterations)) iterations = 0
      else
         call cphf_solve(mol, orbitals, orbital_energies, n_occ, dip, u, error, &
                         max_iter=max_iter, tol=tol, iterations=iterations, &
                         in_core=in_core)
      end if
      if (error%has_error()) return

      ! W maps the canonical occupied orbitals onto the localized ones. With a
      ! core excluded it is an isometry rather than a rotation, and `W W^T` is
      ! then the projector onto the orbitals kept, which is why the sum rule
      ! stays exact.
      call mol%overlap(s)
      allocate (sc(n_ao, n_lmo), w(n_occ, n_lmo))
      call pic_gemm(s, localized, sc)
      call pic_gemm(c_occ, sc, w, transa="T")

      allocate (h_loc(n_vir, n_lmo, 3), u_loc(n_vir, n_lmo, 3))
      allocate (work(n_ao, n_lmo))
      do k = 1, 3
         ! h in the localized basis, straight from the localized orbitals rather
         ! than by transforming the canonical block.
         call pic_gemm(dip(:, :, k), localized, work)
         call pic_gemm(c_vir, work, h_loc(:, :, k), transa="T")
         call pic_gemm(u(:, :, k), w, u_loc(:, :, k))
      end do

      allocate (tensors(3, 3, n_lmo))
      do i = 1, n_lmo
         do l = 1, 3
            do k = 1, 3
               tensors(k, l, i) = 0.0_dp
               do a = 1, n_vir
                  tensors(k, l, i) = tensors(k, l, i) &
                                     - 4.0_dp*h_loc(a, i, k)*u_loc(a, i, l)
               end do
            end do
         end do
      end do

      deallocate (dip, u, localized, s, w, u_loc, h_loc, c_occ, c_vir, work, sc)
   end subroutine distributed_polarizability

   subroutine static_polarizability(mol, orbitals, orbital_energies, n_occ, alpha, &
                                    error, max_iter, tol, iterations, in_core, xc, density)
      !! The dipole polarizability tensor, in Bohr^3
      !!
      !! `alpha_kl = d mu_k / d F_l`, from the response to the three components
      !! of a uniform field. The perturbation for a field `F` is `+ r . F`, since
      !! the electronic dipole operator is `-r`, and contracting the first-order
      !! density against `-r_k` gives
      !!
      !!     alpha_kl = -4 sum_ai (r_k)_ai U^l_ai
      !!
      !! which is `4 h^k A^-1 h^l` once `U = -A^-1 h`, and therefore positive
      !! definite.
      !!
      !! **The origin does not matter.** Shifting it adds a multiple of the
      !! overlap to each dipole matrix, and the occupied-virtual block of the
      !! overlap is exactly zero because the orbitals are orthogonal, so there is
      !! no origin argument to choose.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), intent(out) :: alpha(3, 3)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: iterations
      logical, intent(in), optional :: in_core
      type(xc_context_t), intent(inout), optional :: xc
         !! Makes this a coupled-perturbed *Kohn-Sham* polarizability: the
         !! response operator gains the exchange-correlation kernel. Omitting it
         !! for a Kohn-Sham reference answers a different question.
      real(dp), intent(in), optional :: density(:, :)
         !! The reference density the kernel is evaluated at. Required with `xc`.

      real(dp), allocatable :: dip(:, :, :), u(:, :, :), h(:, :), work(:, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :)
      integer :: n_ao, n_mo, n_vir, k, l

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_vir = n_mo - n_occ

      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, error)
      if (error%has_error()) return

      ! With the kernel when this is a Kohn-Sham reference, without it
      ! otherwise: `density` has to travel with `xc` and means nothing without it.
      if (present(xc)) then
         call cphf_solve(mol, orbitals, orbital_energies, n_occ, dip, u, error, &
                         max_iter=max_iter, tol=tol, iterations=iterations, &
                         in_core=in_core, xc=xc, density=density)
      else
         call cphf_solve(mol, orbitals, orbital_energies, n_occ, dip, u, error, &
                         max_iter=max_iter, tol=tol, iterations=iterations, &
                         in_core=in_core)
      end if
      if (error%has_error()) return

      allocate (c_occ(n_ao, n_occ), c_vir(n_ao, n_vir))
      c_occ = orbitals(:, 1:n_occ)
      c_vir = orbitals(:, n_occ + 1:n_mo)
      allocate (h(n_vir, n_occ), work(n_ao, n_occ))

      do k = 1, 3
         call pic_gemm(dip(:, :, k), c_occ, work)
         call pic_gemm(c_vir, work, h, transa="T")
         do l = 1, 3
            alpha(k, l) = -4.0_dp*sum(h*u(:, :, l))
         end do
      end do

      deallocate (dip, u, h, work, c_occ, c_vir)
   end subroutine static_polarizability

end module mqc_czt_cphf
