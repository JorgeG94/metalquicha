!! The analytic RI-MP2 correlation Hessian, with nothing four-index anywhere
module mqc_czt_ri_mp2_hessian
   !! The second derivative of a fitted MP2 correlation energy over an exact
   !! restricted Hartree-Fock reference, transcribed from the lean path of
   !! `tools/cpu_validation/ri_mp2_hessian.py`. Read that header first: it
   !! carries the formulation, every convention that cost time, and the checks.
   !!
   !! It is the ladder of `mqc_czt_mp2_hessian` with one change of principle:
   !! the correlation functional contracts the fitted tensor
   !! `(pr|P) J^-1_PQ (Q|qs)` and everything belonging to the reference keeps the
   !! exact one. Each use of a four-index tensor there reduces here to one of
   !! two things.
   !!
   !! * The reference operator `T1(D) = C^T (2J[D] - K[D]) C` of a symmetric MO
   !!   matrix, and its first derivative at fixed orbitals. Both go through the
   !!   direct Fock builders (`response_mean_field`, `h1_contract`) in batches
   !!   over the perturbations, so every stage below is a handful of integral
   !!   passes and never a `nao^4` array.
   !! * The fitted tensor, only ever against the cumulant, as the three-index
   !!   objects `Gamma^P_ia`, `gamma_PQ` and `Theta^Q_ia`.
   !!
   !! **Memory.** Resident once: the three-index MO integrals and fitting
   !! coefficients (`n^2 P` each, `P` auxiliary functions), the amplitudes
   !! (`o^2 v^2`) and a few `n^2` matrices. Resident for all `3N` perturbations:
   !! `n^2` matrices only. Transient per batch of atoms: the three-index
   !! first-derivative integrals of those atoms (`3 nao^2 P` each), their MO
   !! images and the three-index densities. Nothing of size `3N * o v P` exists;
   !! the contraction of those densities with every other nucleus's derivative
   !! integrals is one sweep per batch (`contract_first_3c`).
   !!
   !! The fitted energy differentiated is the one `czt_ri_mp2_gradient` reports
   !! with `fit_response = .false.`: an exact reference, a fitted correlation.
   !! The inverse metric is the pseudo-inverse that gradient uses (eigenvalues
   !! below `FIT_METRIC_THRESHOLD` dropped); the derivative `-J^-1 dJ J^-1` is exact
   !! only when none is, so a metric that drops a mode is refused.
   use, intrinsic :: iso_fortran_env, only: int64
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_memory, only: memory_budget
   use mqc_czt_integrals, only: czt_molecule_t, atom_ao_blocks, three_centre, &
                                two_centre, build_df_shell_table, max_block, shell_dim, &
                                metric_inverse_sqrt
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_ecp, only: ecp_refuses_derivatives
   use mqc_czt_cphf, only: cphf_solve
   use mqc_czt_response_product, only: response_mean_field
   use mqc_czt_gradient, only: two_centre_deriv
   use mqc_czt_xc, only: xc_context_t
   use mqc_czt_xc_hessian, only: xc_potential_deriv, xc_potential_hessian
   use mqc_czt_mp2_hessian, only: ks_pair_kernel, ks_kernel_deriv_stack, ks_dl_pair_kernel
   use mqc_czt_multipole, only: multipole_matrices, dipole_integral_derivatives
   use mqc_czt_hess_ints, only: hess_1e_block, hess_rinv_block, hess_2e_contract, &
                                hess_3c2c_contract, h1_contract, &
                                HESS_KIN_II, HESS_KIN_IJ, HESS_NUC_II, HESS_NUC_IJ, &
                                HESS_OVLP_II, HESS_OVLP_IJ, HESS_RINV_II, HESS_RINV_IJ
   use mqc_czt_hessian, only: hcore_deriv_atom, overlap_deriv_atom, solve_mo1_batch
   use libcint_fortran, only: LIBCINT_ATOM_OF, libcint_3c2e_ip1_sph, libcint_3c2e_ip1_cart, &
                              libcint_3c2e_ip2_sph, libcint_3c2e_ip2_cart
   implicit none
   private

   public :: ri_mp2_correlation_hessian
   public :: ri_hessian_parts_t
   public :: contract_first_3c
   public :: inverse_metric

   integer, parameter :: DEFAULT_ATOM_BATCH_CAP = 1000000

   type :: ri_hessian_parts_t
      !! The four additive pieces of the correlation Hessian, `(3N, 3N)` with
      !! `X = 3*(atom-1) + component` the row and `Y` the column, for a test that
      !! wants to locate a disagreement
      real(dp), allocatable :: pass1_reference(:, :)
         !! `D_rel . h^(XY)`, `I' . S^(XY)` and the Fock fold `D_rel (x) P_occ`
      real(dp), allocatable :: pass1_fitted(:, :)
         !! Second derivatives of `(mn|P)`, `(P|Q)` against `Gamma`, `gamma`. The
         !! pass-1 cross terms from the first-order change of the fitting
         !! coefficients at fixed cumulant ride in `response`, in the same
         !! sweep as the pair sum they are a different density of
      real(dp), allocatable :: orbital(:, :)
         !! `2 U^Y_ai X~^(X)_ai + S^(Y)_pq I''~^(X)_pq [+ P^(X)_pq f^(Y)_pq]`
      real(dp), allocatable :: response(:, :)
         !! `d_Y D_rel . f^(X) + d_Y Gamma . R^(X) + d_Y I' . S^(X)`
   end type ri_hessian_parts_t

   type :: ri_state_t
      !! Everything the stages share
      type(czt_molecule_t), pointer :: mol => null()
      type(czt_molecule_t), pointer :: aux => null()
      integer :: nao = 0, n = 0, no = 0, nf = 0, na = 0, nv = 0, np = 0
      integer :: natm = 0, n3 = 0, batch = 1
      logical :: dh = .false.
         !! A double hybrid's perturbative term over a Kohn-Sham reference
      real(dp) :: kf = 1.0_dp
         !! Exact-exchange fraction of the reference's operator
      type(xc_context_t), pointer :: xc => null()
      real(dp), allocatable :: scf_density(:, :)
      real(dp), allocatable :: c(:, :), eps(:)
      real(dp), allocatable :: bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: jinv(:, :)
      real(dp), allocatable :: imo(:, :, :), cmo(:, :, :)
         !! `(pq|P)` and `c^P_pq = sum_Q J^-1_PQ (pq|Q)`, all MO pairs, `(n, n, np)`
      real(dp), allocatable :: t2(:, :, :, :), u2(:, :)
         !! Amplitudes `(i,j,a,b)` and `u_ijab = 2 t - t^swap` as `(ia, jb)`
      real(dp), allocatable :: iov2(:, :), cov2(:, :)
         !! The active-occupied virtual block of `imo` and `cmo`, `(na*nv, np)`
      real(dp), allocatable :: tt_o(:, :), ut_o(:, :), tp_v(:, :), up_v(:, :)
         !! `t` and `u` rearranged so the one-particle densities are matrix products
      real(dp), allocatable :: gam3(:, :, :), gam2(:, :)
         !! `Gamma^P_ia = u c`, `(na, nv, np)`, and `gamma_PQ`
      real(dp), allocatable :: d0(:, :), drel(:, :), ip(:, :), pco(:, :), z(:, :)
      real(dp), allocatable :: zsym(:, :), psym(:, :)
      real(dp), allocatable :: t1_d0(:, :), t1_drel(:, :), t1_z(:, :), t1_p(:, :)
   end type ri_state_t

contains

   subroutine ri_mp2_correlation_hessian(mol, aux, coeff, orbital_energies, n_occ, &
                                         hess_corr, error, n_frozen, tol, batch, &
                                         atom_batch, parts, xc, scf_density, pt2_scale, &
                                         dipole_derivatives)
      !! The correlation block of the RI-MP2 Hessian, `(3, 3, natm, natm)`
      !!
      !! Add it to the RHF Hessian (`rhf_hessian`) for the total. `hess_corr(a, b,
      !! A, B)` is the derivative along component `a` of atom `A` and then
      !! component `b` of atom `B`.
      type(czt_molecule_t), intent(in), target :: mol
      type(czt_molecule_t), intent(in), target :: aux   !! Same atoms
      real(dp), intent(in) :: coeff(:, :)               !! C, `(n_ao, n_mo)`
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ                      !! Doubly occupied, core included
      real(dp), allocatable, intent(out) :: hess_corr(:, :, :, :)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_frozen
      real(dp), intent(in), optional :: tol
         !! Tolerance of every Z-vector and coupled-perturbed solve; `1e-12`
      integer, intent(in), optional :: batch
         !! Perturbations per integral pass in the Fock-type builds and the
         !! coupled-perturbed solve; `hessian_response_batch`'s default
      integer, intent(in), optional :: atom_batch
         !! Atoms whose three-index first-derivative integrals are held at once.
         !! Absent, sized from `memory_budget`; present, taken as given, which is
         !! how a test runs several batches on a molecule that fits in one.
      type(ri_hessian_parts_t), intent(out), optional :: parts
      type(xc_context_t), intent(inout), optional, target :: xc
         !! The functional of a Kohn-Sham reference. Present, this is a *double
         !! hybrid's* perturbative term: the correlation block alone, over that
         !! reference, with exact exchange at the functional's fraction and the
         !! exchange-correlation kernel in every response, scaled once by
         !! `pt2_scale`. The reference's own Hessian is `ks_hessian`'s. LDA and GGA
         !! only; no frozen core.
      real(dp), intent(in), optional :: scf_density(:, :)
         !! The converged reference density, where the kernel is evaluated.
         !! Required with `xc`.
      real(dp), intent(in), optional :: pt2_scale
         !! The functional's PT2 coefficient, applied once to the finished block.
         !! The Z-vector solves are linear, so scaling inputs as well would count
         !! it twice.
      real(dp), allocatable, intent(out), optional :: dipole_derivatives(:, :)
         !! The perturbative term's share of `d mu_a / dR_(X,b)`, `(3, 3N)`,
         !! scaled by `pt2_scale`. Requires `xc`: over Hartree-Fock an MP2 dipole
         !! derivative is a different quantity.

      type(ri_state_t) :: st
      real(dp) :: cscale
      real(dp), allocatable :: ddip(:, :)
      real(dp), allocatable :: h1(:, :, :, :), s1(:, :, :, :), mo1(:, :, :, :)
      real(dp), allocatable :: fx(:, :, :), sx(:, :, :)
      real(dp), allocatable :: pass1_ref(:, :, :, :), pass1_fit(:, :, :, :)
      real(dp), allocatable :: p1_ref(:, :), p1_fit(:, :), pair(:, :), orb(:, :), resp(:, :)
      real(dp) :: use_tol
      integer :: use_batch, nab, x, y, ax, cx, by, cy

      if (error%has_error()) return
      if (ecp_refuses_derivatives(mol%core_electrons, "RI-MP2 Hessian", error)) return

      use_tol = 1.0e-12_dp
      if (present(tol)) use_tol = tol
      use_batch = 12
      if (present(batch)) use_batch = max(1, batch)

      cscale = 1.0_dp
      if (present(pt2_scale)) cscale = pt2_scale
      if (present(dipole_derivatives) .and. .not. present(xc)) then
         call error%set(ERROR_VALIDATION, "the RI-MP2 Hessian's dipole derivative was "// &
                        "asked for without a functional: over Hartree-Fock an MP2 dipole "// &
                        "derivative is a different quantity, and this returns the double "// &
                        "hybrid's or nothing.")
         return
      end if
      call setup_state(st, mol, aux, coeff, orbital_energies, n_occ, n_frozen, use_batch, &
                       error, xc=xc, scf_density=scf_density)
      if (error%has_error()) return
      call build_fitted_objects(st, error)
      if (error%has_error()) return
      call solve_unperturbed(st, use_tol, error)
      if (error%has_error()) return

      if (present(atom_batch)) then
         nab = max(1, min(atom_batch, st%natm))
      else
         nab = default_atom_batch(st)
      end if

      call skeleton_matrices(st, h1, s1, fx, sx, error)
      if (error%has_error()) return
      if (st%dh) then
         call solve_mo1_batch(mol, st%c, st%eps, st%no, h1, s1, mo1, error, xc=xc, &
                              reference=scf_density, k_scale=st%kf, max_iter=200, &
                              tol=use_tol, batch=use_batch, bounds=st%bounds)
      else
         call solve_mo1_batch(mol, st%c, st%eps, st%no, h1, s1, mo1, error, max_iter=200, &
                              tol=use_tol, batch=use_batch, bounds=st%bounds)
      end if
      if (error%has_error()) return
      deallocate (h1, s1)

      if (present(dipole_derivatives)) then
         call perturbed_ladder(st, fx, sx, mo1, nab, use_tol, pair, orb, resp, error, ddip)
      else
         call perturbed_ladder(st, fx, sx, mo1, nab, use_tol, pair, orb, resp, error)
      end if
      if (error%has_error()) return

      ! ---- pass 1 -----------------------------------------------------------
      allocate (pass1_ref(3, 3, st%natm, st%natm), pass1_fit(3, 3, st%natm, st%natm))
      pass1_ref = 0.0_dp
      pass1_fit = 0.0_dp
      call pass1_reference_share(st, pass1_ref, error)
      if (error%has_error()) return
      call pass1_fitted_share(st, pass1_fit, error)
      if (error%has_error()) return

      allocate (hess_corr(3, 3, st%natm, st%natm))
      allocate (p1_ref(st%n3, st%n3), p1_fit(st%n3, st%n3))
      ! Once, here and nowhere else: the gradient's own rule, one order up.
      do y = 1, st%n3
         by = (y - 1)/3 + 1
         cy = y - 3*(by - 1)
         do x = 1, st%n3
            ax = (x - 1)/3 + 1
            cx = x - 3*(ax - 1)
            p1_ref(x, y) = pass1_ref(cx, cy, ax, by)
            p1_fit(x, y) = pass1_fit(cx, cy, ax, by)
            hess_corr(cx, cy, ax, by) = cscale*(pass1_ref(cx, cy, ax, by) &
                                                + pass1_fit(cx, cy, ax, by) &
                                                + orb(x, y) + pair(x, y) + resp(x, y))
         end do
      end do
      if (present(dipole_derivatives)) then
         ddip = cscale*ddip
         call move_alloc(ddip, dipole_derivatives)
      end if
      if (present(parts)) then
         call move_alloc(p1_ref, parts%pass1_reference)
         call move_alloc(p1_fit, parts%pass1_fitted)
         call move_alloc(orb, parts%orbital)
         resp = resp + pair
         call move_alloc(resp, parts%response)
         parts%pass1_reference = cscale*parts%pass1_reference
         parts%pass1_fitted = cscale*parts%pass1_fitted
         parts%orbital = cscale*parts%orbital
         parts%response = cscale*parts%response
      end if
   end subroutine ri_mp2_correlation_hessian

   ! ==========================================================================
   ! setup
   ! ==========================================================================

   subroutine setup_state(st, mol, aux, coeff, eps, n_occ, n_frozen, batch, error, xc, &
                          scf_density)
      type(ri_state_t), intent(out) :: st
      type(czt_molecule_t), intent(in), target :: mol, aux
      real(dp), intent(in) :: coeff(:, :), eps(:)
      integer, intent(in) :: n_occ, batch
      integer, intent(in), optional :: n_frozen
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: scf_density(:, :)

      st%mol => mol
      st%aux => aux
      st%nao = mol%nao
      st%n = size(coeff, 2)
      st%no = n_occ
      st%nf = 0
      if (present(n_frozen)) st%nf = n_frozen
      st%na = st%no - st%nf
      st%nv = st%n - st%no
      st%np = aux%nao
      st%natm = mol%natm
      st%n3 = 3*mol%natm
      st%batch = batch
      if (size(coeff, 1) /= st%nao .or. size(eps) /= st%n) then
         call error%set(ERROR_VALIDATION, "the RI-MP2 Hessian was handed orbitals that "// &
                        "do not fit the basis.")
         return
      end if
      if (st%nf < 0 .or. st%na < 1 .or. st%nv < 1) then
         call error%set(ERROR_VALIDATION, "the RI-MP2 Hessian needs at least one "// &
                        "correlated occupied and one virtual orbital.")
         return
      end if
      if (aux%natm /= mol%natm) then
         call error%set(ERROR_VALIDATION, "the RI-MP2 Hessian's auxiliary basis is on a "// &
                        "different number of atoms than the orbital basis.")
         return
      end if
      if (mol%cartesian .neqv. aux%cartesian) then
         call error%set(ERROR_VALIDATION, "the RI-MP2 Hessian's orbital and auxiliary "// &
                        "bases must be in the same angular form.")
         return
      end if
      if (present(xc)) then
         if (.not. present(scf_density)) then
            call error%set(ERROR_VALIDATION, "the double-hybrid RI-MP2 Hessian needs the "// &
                           "converged reference density alongside the functional: the "// &
                           "kernel is evaluated at it.")
            return
         end if
         if (xc%any_mgga .or. xc%range_separated .or. xc%nlc_b /= 0.0_dp &
             .or. xc%nlc_c /= 0.0_dp) then
            call error%set(ERROR_VALIDATION, "the double-hybrid RI-MP2 Hessian covers LDA "// &
                           "and GGA functionals: the kernel's third derivative, the "// &
                           "screened exchange and the non-local term of the others are "// &
                           "not in the response.")
            return
         end if
         if (st%nf > 0) then
            call error%set(ERROR_VALIDATION, "the double-hybrid RI-MP2 Hessian has not "// &
                           "been validated with a frozen core: the core-active rotation "// &
                           "is not checked against a Kohn-Sham operator.")
            return
         end if
         st%dh = .true.
         st%xc => xc
         st%kf = xc%exx_fraction
         st%scf_density = scf_density
      end if
      st%c = coeff
      st%eps = eps
      call schwarz_bounds(mol, st%bounds, error)
      if (error%has_error()) return
      allocate (st%zero_h(st%nao, st%nao))
      st%zero_h = 0.0_dp
   end subroutine setup_state

   function default_atom_batch(st) result(nab)
      !! Atoms per three-index derivative batch, from what the machine has
      !!
      !! A batch holds the AO first-derivative integrals of its atoms
      !! (`3 nao^2 P` doubles each) and, per component, an MO image and a
      !! density (`n^2 P + o v P`), so three times the AO figure is a fair price.
      type(ri_state_t), intent(in) :: st
      integer :: nab

      real(dp) :: per_atom, budget

      per_atom = 8.0_dp*3.0_dp*real(st%nao, dp)**2*real(st%np, dp)*3.0_dp
      budget = memory_budget(4.0e9_dp, 0.25_dp)
      nab = int(max(1.0_dp, min(real(DEFAULT_ATOM_BATCH_CAP, dp), budget/per_atom)))
      nab = min(nab, st%natm)
   end function default_atom_batch

   subroutine inverse_metric(metric, jinv, n_dropped, error)
      !! `J^-1` as the square of `metric_inverse_sqrt`, and the modes it dropped
      !!
      !! The fit `czt_ri_mp2_gradient` differentiates, so the two cannot disagree
      !! about which modes survive.
      real(dp), intent(in) :: metric(:, :)
      real(dp), allocatable, intent(out) :: jinv(:, :)
      integer, intent(out) :: n_dropped
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: half(:, :)

      call metric_inverse_sqrt(metric, half, error, n_dropped=n_dropped)
      if (error%has_error()) return
      allocate (jinv(size(metric, 1), size(metric, 1)))
      call pic_gemm(half, half, jinv)
   end subroutine inverse_metric

   subroutine build_fitted_objects(st, error)
      !! `(pq|P)`, `c^P_pq`, the amplitudes, `Gamma`, `gamma`
      type(ri_state_t), intent(inout), target :: st
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: three(:, :), metric(:, :), block(:, :), half(:, :)
      real(dp), allocatable :: kia(:, :)
      real(dp), pointer :: i2(:, :), c2(:, :)
      integer :: n_dropped, p, i, j, a, b, n, no, nf, na, nv, np, nao

      n = st%n
      no = st%no
      nf = st%nf
      na = st%na
      nv = st%nv
      np = st%np
      nao = st%nao

      ! TODO(mqc): `imo`, `cmo` and the per-perturbation `d(pq|P)` are held for all
      ! MO pairs, `n^2 P` each, where only the (n, v) and (o, n) blocks are read
      ! by `gp_term`; keeping those two would cut them by up to a factor of two
      ! to four at the sizes where they matter.
      call three_centre(st%mol, st%aux, three)
      allocate (st%imo(n, n, np), block(nao, nao), half(nao, n))
      do p = 1, np
         block = reshape(three(:, p), [nao, nao])
         call pic_gemm(block, st%c, half)
         call pic_gemm(st%c, half, st%imo(:, :, p), transa="T")
      end do
      deallocate (three, block, half)

      call two_centre(st%aux, metric)
      call inverse_metric(metric, st%jinv, n_dropped, error)
      if (error%has_error()) return
      if (n_dropped > 0) then
         call error%set(ERROR_VALIDATION, "the RI-MP2 Hessian differentiates the "// &
                        "inverse of the auxiliary metric as -J^-1 dJ J^-1, which is only "// &
                        "exact when the fit drops no mode, and this metric has "// &
                        trim(int_text(n_dropped))//" below the threshold; use a "// &
                        "better-conditioned auxiliary basis or finite differences of the "// &
                        "gradient.")
         return
      end if
      deallocate (metric)

      allocate (st%cmo(n, n, np))
      i2(1:n*n, 1:np) => st%imo
      c2(1:n*n, 1:np) => st%cmo
      call pic_gemm(i2, st%jinv, c2)

      ! The amplitudes, from the occupied-virtual block.
      allocate (st%iov2(na*nv, np), st%cov2(na*nv, np), kia(na*nv, na*nv))
      do p = 1, np
         do a = 1, nv
            do i = 1, na
               st%iov2(i + na*(a - 1), p) = st%imo(nf + i, no + a, p)
               st%cov2(i + na*(a - 1), p) = st%cmo(nf + i, no + a, p)
            end do
         end do
      end do
      call pic_gemm(st%iov2, st%cov2, kia, transb="T")
      allocate (st%t2(na, na, nv, nv), st%u2(na*nv, na*nv))
      do b = 1, nv
         do a = 1, nv
            do j = 1, na
               do i = 1, na
                  st%t2(i, j, a, b) = kia(i + na*(a - 1), j + na*(b - 1)) &
                                      /(st%eps(nf + i) + st%eps(nf + j) &
                                        - st%eps(no + a) - st%eps(no + b))
               end do
            end do
         end do
      end do
      do b = 1, nv
         do a = 1, nv
            do j = 1, na
               do i = 1, na
                  st%u2(i + na*(a - 1), j + na*(b - 1)) = 2.0_dp*st%t2(i, j, a, b) &
                                                          - st%t2(i, j, b, a)
               end do
            end do
         end do
      end do
      allocate (st%gam3(na, nv, np), st%gam2(np, np))
      block
         real(dp), allocatable :: gtmp(:, :)
         allocate (gtmp(na*nv, np))
         call pic_gemm(st%u2, st%cov2, gtmp)
         st%gam3 = reshape(gtmp, [na, nv, np])
         call pic_gemm(gtmp, st%cov2, st%gam2, transa="T")
      end block
      deallocate (kia)
      call arrange_amplitudes(st, st%t2, st%tt_o, st%tp_v)
      block
         real(dp), allocatable :: uu(:, :, :, :)
         allocate (uu(na, na, nv, nv))
         do b = 1, nv
            do a = 1, nv
               uu(:, :, a, b) = 2.0_dp*st%t2(:, :, a, b) - st%t2(:, :, b, a)
            end do
         end do
         call arrange_amplitudes(st, uu, st%ut_o, st%up_v)
      end block
   end subroutine build_fitted_objects

   function int_text(i) result(text)
      integer, intent(in) :: i
      character(len=:), allocatable :: text
      character(len=16) :: buf

      write (buf, "(i0)") i
      text = trim(buf)
   end function int_text

   ! ==========================================================================
   ! the reference operator and the fitted Lagrangian terms
   ! ==========================================================================

   subroutine t1_batch(st, dmo, t1, error, kernel)
      !! `T1(D) = C^T (2J[D] - k K[D]) C + vk(D)/2` for a stack of symmetric MO
      !! matrices
      !!
      !! `response_mean_field` returns `J - kK/2` of an AO density, hence the
      !! factor of two; `k` is the reference's exchange fraction. For a
      !! Kohn-Sham reference `vk` is `ks_pair_kernel`, the kernel's share of the
      !! fold `sum_rs D_rs (L_rpsq + L_rqsp)`, of which `T1` carries half because
      !! the fold is `T1 + T1^T`. The stack goes through in chunks of `st%batch`.
      type(ri_state_t), intent(in) :: st
      real(dp), intent(in) :: dmo(:, :, :)
      real(dp), intent(out) :: t1(:, :, :)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: kernel
         !! Leave the kernel's share out. Default is to carry it, for a
         !! Kohn-Sham reference; the total derivative asks for the exchange part
         !! alone because it differentiates the kernel itself.

      real(dp), allocatable :: dens(:, :, :), g(:, :, :), half(:, :), vk(:, :)
      integer :: first, last, wide, q, nset
      logical :: with_kernel

      if (error%has_error()) return
      with_kernel = st%dh
      if (present(kernel)) with_kernel = kernel .and. st%dh
      nset = size(dmo, 3)
      allocate (half(st%nao, st%n))
      first = 1
      do while (first <= nset)
         last = min(first + st%batch - 1, nset)
         wide = last - first + 1
         allocate (dens(st%nao, st%nao, wide))
         do q = 1, wide
            call pic_gemm(st%c, dmo(:, :, first + q - 1), half)
            call pic_gemm(half, st%c, dens(:, :, q), transb="T")
         end do
         call response_mean_field(st%mol, dens, st%zero_h, g, error, direct=.true., &
                                  bounds=st%bounds, k_scale=st%kf)
         if (error%has_error()) return
         do q = 1, wide
            call pic_gemm(g(:, :, q), st%c, half)
            call pic_gemm(st%c, half, t1(:, :, first + q - 1), transa="T", alpha=2.0_dp)
         end do
         deallocate (dens, g)
         first = last + 1
      end do
      if (with_kernel) then
         allocate (vk(st%n, st%n))
         do q = 1, nset
            call ks_pair_kernel(st%xc, st%mol, st%c, st%scf_density, dmo(:, :, q), vk, error)
            if (error%has_error()) return
            t1(:, :, q) = t1(:, :, q) + 0.5_dp*vk
         end do
      end if
   end subroutine t1_batch

   subroutine t1x_stack(st, dmo, stack, error, kd)
      !! `T1x(D)[X]` for every nuclear coordinate: the first derivative at fixed
      !! orbitals, `(n, n, 3N)`
      !!
      !! For a Kohn-Sham reference the kernel's skeleton derivative is in it, at
      !! half weight, and `kd` returns the whole of that derivative
      !! (`ks_kernel_deriv_stack`) for the total derivative to reuse.
      type(ri_state_t), intent(in) :: st
      real(dp), intent(in) :: dmo(:, :)
      real(dp), intent(out) :: stack(:, :, :)
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: kd(:, :, :)

      real(dp), allocatable :: dao(:, :), half(:, :), h1(:, :, :, :), kds(:, :, :)
      integer :: a, k, x

      allocate (dao(st%nao, st%nao), half(st%nao, st%n))
      call pic_gemm(st%c, dmo, half)
      call pic_gemm(half, st%c, dao, transb="T")
      call h1_contract(st%mol, dao, h1, error, k_scale=st%kf, bounds=st%bounds)
      if (error%has_error()) return
      do a = 1, st%natm
         do k = 1, 3
            x = 3*(a - 1) + k
            call pic_gemm(h1(:, :, k, a), st%c, half)
            call pic_gemm(st%c, half, stack(:, :, x), transa="T", alpha=2.0_dp)
         end do
      end do
      if (st%dh) then
         call ks_kernel_deriv_stack(st%xc, st%mol, st%c, st%scf_density, dmo, kds, error)
         if (error%has_error()) return
         stack = stack + 0.5_dp*kds
         if (present(kd)) call move_alloc(kds, kd)
      end if
   end subroutine t1x_stack

   subroutine gp_term(st, a, g, out)
      !! `out[p,i] = sum_{P,a} A^P_pa G^P_ia` (i active) and
      !! `out[p,a] = sum_{P,i} A^P_ip G^P_ia`
      !!
      !! The whole of every two-particle Lagrangian term. `a` is an integral-like
      !! `(n, n, np)` tensor and `g` a density-like `(na, nv, np)` one.
      type(ri_state_t), intent(in) :: st
      real(dp), intent(in) :: a(:, :, :), g(:, :, :)
      real(dp), intent(out) :: out(:, :)

      real(dp), allocatable :: acc(:, :)
      integer :: p, n, no, nf, na, nv

      n = st%n
      no = st%no
      nf = st%nf
      na = st%na
      nv = st%nv
      out = 0.0_dp
      allocate (acc(n, n))
      do p = 1, st%np
         call pic_gemm(a(:, no + 1:n, p), g(:, :, p), acc(:, 1:na), transb="T")
         out(:, nf + 1:no) = out(:, nf + 1:no) + acc(:, 1:na)
         call pic_gemm(a(nf + 1:no, :, p), g(:, :, p), acc(:, 1:nv), transa="T")
         out(:, no + 1:n) = out(:, no + 1:n) + acc(:, 1:nv)
      end do
   end subroutine gp_term

   subroutine theta_term(st, a_ov, j1, th)
      !! `Theta^Q_ia = sum_jb u_ijab A^Q_jb - sum_R J1_QR Gamma^R_ia`
      type(ri_state_t), intent(in) :: st
      real(dp), intent(in) :: a_ov(:, :, :)   !! `(na, nv, np)`
      real(dp), intent(in) :: j1(:, :)
      real(dp), intent(out) :: th(:, :, :)

      real(dp), allocatable :: a2(:, :), g2(:, :), t2(:, :)
      integer :: na, nv, np

      na = st%na
      nv = st%nv
      np = st%np
      allocate (a2(na*nv, np), g2(na*nv, np), t2(na*nv, np))
      a2 = reshape(a_ov, [na*nv, np])
      g2 = reshape(st%gam3, [na*nv, np])
      call pic_gemm(st%u2, a2, t2)
      call pic_gemm(g2, j1, t2, transb="T", alpha=-1.0_dp, beta=1.0_dp)
      th = reshape(t2, [na, nv, np])
   end subroutine theta_term

   subroutine corr_lagrangian(st, a_full, j1, out, dgam3)
      !! `4 [gp(A, Gamma) + gp(c, Theta(A, J1)) (+ gp(I, dGamma))]`
      !!
      !! The correlation functional's whole contribution to a Lagrangian-shaped
      !! matrix: a first-derivative (or total-derivative) integral tensor `A` with
      !! its metric derivative `J1`, and optionally the response of the
      !! three-index density.
      type(ri_state_t), intent(in) :: st
      real(dp), intent(in) :: a_full(:, :, :), j1(:, :)
      real(dp), intent(out) :: out(:, :)
      real(dp), intent(in), optional :: dgam3(:, :, :)

      real(dp), allocatable :: th(:, :, :), part(:, :)
      integer :: no, nf, n

      n = st%n
      no = st%no
      nf = st%nf
      allocate (th(st%na, st%nv, st%np), part(n, n))
      call theta_term(st, a_full(nf + 1:no, no + 1:n, :), j1, th)
      call gp_term(st, a_full, st%gam3, out)
      call gp_term(st, st%cmo, th, part)
      out = out + part
      if (present(dgam3)) then
         call gp_term(st, st%imo, dgam3, part)
         out = out + part
      end if
      out = 4.0_dp*out
   end subroutine corr_lagrangian

   ! ==========================================================================
   ! the unperturbed ladder
   ! ==========================================================================

   subroutine solve_response(st, rhs, response, tol, error)
      !! `(n_vir, n_occ, k)` coupled-perturbed solves against the reference's own
      !! orbital Hessian: Hartree-Fock, or Kohn-Sham with the kernel
      type(ri_state_t), intent(in) :: st
      real(dp), intent(in) :: rhs(:, :, :)
      real(dp), allocatable, intent(out) :: response(:, :, :)
      real(dp), intent(in) :: tol
      type(error_t), intent(inout) :: error

      if (st%dh) then
         call cphf_solve(st%mol, st%c, st%eps, st%no, response=response, error=error, &
                         mo_rhs=rhs, tol=tol, max_iter=200, xc=st%xc, density=st%scf_density)
      else
         call cphf_solve(st%mol, st%c, st%eps, st%no, response=response, error=error, &
                         mo_rhs=rhs, tol=tol, max_iter=200)
      end if
   end subroutine solve_response

   subroutine solve_unperturbed(st, tol, error)
      !! Amplitudes to relaxed density and Lagrangian, from the three-index
      !! objects and the reference operator only
      type(ri_state_t), intent(inout) :: st
      real(dp), intent(in) :: tol
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: ip0(:, :), t1p(:, :, :), dd(:, :, :), x(:, :)
      real(dp), allocatable :: rhs(:, :, :), zsol(:, :, :), tmp(:, :, :), psym3(:, :, :)
      real(dp), allocatable :: t3(:, :)
      integer :: n, no, nf, na, nv, i, a, c

      n = st%n
      no = st%no
      nf = st%nf
      na = st%na
      nv = st%nv

      ! ---- the unrelaxed density, d0 ---------------------------------------
      allocate (st%d0(n, n))
      st%d0 = 0.0_dp
      allocate (t3(na, na))
      call pic_gemm(st%tt_o, st%ut_o, t3, transb="T", alpha=-2.0_dp)
      st%d0(nf + 1:no, nf + 1:no) = t3
      deallocate (t3)
      allocate (t3(nv, nv))
      call pic_gemm(st%up_v, st%tp_v, t3, transb="T", alpha=2.0_dp)
      st%d0(no + 1:n, no + 1:n) = t3
      deallocate (t3)

      ! ---- I'(d0), core<->active rotation, Z-vector -------------------------
      allocate (tmp(n, n, 1), t1p(n, n, 1), ip0(n, n))
      tmp(:, :, 1) = st%d0
      call t1_batch(st, tmp, t1p, error)
      if (error%has_error()) return
      allocate (st%t1_d0(n, n))
      st%t1_d0 = t1p(:, :, 1)
      call lagrangian_unrelaxed(st, st%d0, st%t1_d0, ip0)

      if (nf > 0) then
         allocate (st%pco(nf, na))
         do i = 1, na
            do c = 1, nf
               st%pco(c, i) = (ip0(c, nf + i) - ip0(nf + i, c))/(st%eps(c) - st%eps(nf + i))
            end do
         end do
         allocate (st%psym(n, n))
         st%psym = 0.0_dp
         do i = 1, na
            do c = 1, nf
               st%psym(nf + i, c) = st%pco(c, i)
               st%psym(c, nf + i) = st%pco(c, i)
            end do
         end do
         tmp(:, :, 1) = st%psym
         call t1_batch(st, tmp, t1p, error)
         if (error%has_error()) return
         allocate (st%t1_p(n, n))
         st%t1_p = t1p(:, :, 1)
      end if

      allocate (x(no, nv), rhs(nv, no, 1))
      do a = 1, nv
         do i = 1, no
            x(i, a) = ip0(i, no + a) - ip0(no + a, i)
            if (nf > 0) x(i, a) = x(i, a) + st%t1_p(no + a, i)
            rhs(a, i, 1) = -x(i, a)
         end do
      end do
      call solve_response(st, rhs, zsol, tol, error)
      if (error%has_error()) return
      allocate (st%z(no, nv))
      st%z = transpose(zsol(:, :, 1))

      allocate (st%drel(n, n))
      st%drel = st%d0
      if (nf > 0) then
         do i = 1, na
            do c = 1, nf
               st%drel(nf + i, c) = st%drel(nf + i, c) + st%pco(c, i)
               st%drel(c, nf + i) = st%drel(c, nf + i) + st%pco(c, i)
            end do
         end do
      end if
      do a = 1, nv
         do i = 1, no
            st%drel(no + a, i) = st%drel(no + a, i) - st%z(i, a)
            st%drel(i, no + a) = st%drel(i, no + a) - st%z(i, a)
         end do
      end do
      allocate (st%zsym(n, n))
      st%zsym = 0.0_dp
      do a = 1, nv
         do i = 1, no
            st%zsym(i, no + a) = st%z(i, a)
            st%zsym(no + a, i) = st%z(i, a)
         end do
      end do

      allocate (dd(n, n, 2), psym3(n, n, 2))
      dd(:, :, 1) = st%drel
      dd(:, :, 2) = st%zsym
      call t1_batch(st, dd, psym3, error)
      if (error%has_error()) return
      allocate (st%t1_drel(n, n), st%t1_z(n, n), st%ip(n, n))
      st%t1_drel = psym3(:, :, 1)
      st%t1_z = psym3(:, :, 2)
      call lagrangian_unrelaxed(st, st%drel, st%t1_drel, st%ip)
   end subroutine solve_unperturbed

   subroutine lagrangian_unrelaxed(st, d, t1d, ip)
      !! `I'_pq = -1/2 [ eps_p (D + D^T)_pq + 4 gp + (T1(D) + T1(D)^T)_pq for q occ ]`
      type(ri_state_t), intent(in) :: st
      real(dp), intent(in) :: d(:, :), t1d(:, :)
      real(dp), intent(out) :: ip(:, :)

      real(dp), allocatable :: g(:, :)
      integer :: p, q

      allocate (g(st%n, st%n))
      call gp_term(st, st%imo, st%gam3, g)
      do q = 1, st%n
         do p = 1, st%n
            ip(p, q) = st%eps(p)*(d(p, q) + d(q, p)) + 4.0_dp*g(p, q)
         end do
      end do
      do q = 1, st%no
         do p = 1, st%n
            ip(p, q) = ip(p, q) + t1d(p, q) + t1d(q, p)
         end do
      end do
      ip = -0.5_dp*ip
   end subroutine lagrangian_unrelaxed

   subroutine arrange_amplitudes(st, t, mat_o, mat_v)
      !! `t(i,j,a,b)` as the two matrices the one-particle densities multiply
      !!
      !! `mat_o(i, m+na*(e-1+nv*(f-1))) = t(i,m,e,f)` and
      !! `mat_v(a, k+na*(m-1+na*(e-1))) = t(k,m,a,e)`, so that
      !! `D_ij = -2 sum_mef t_imef u_jmef` and `D_ab = 2 sum_kme u_kmae t_kmbe`
      !! are each one product.
      type(ri_state_t), intent(in) :: st
      real(dp), intent(in) :: t(:, :, :, :)
      real(dp), allocatable, intent(out) :: mat_o(:, :), mat_v(:, :)

      integer :: na, nv, i, m, e, f, k, a

      na = st%na
      nv = st%nv
      allocate (mat_o(na, na*nv*nv), mat_v(nv, na*na*nv))
      do f = 1, nv
         do e = 1, nv
            do m = 1, na
               do i = 1, na
                  mat_o(i, m + na*(e - 1 + nv*(f - 1))) = t(i, m, e, f)
               end do
            end do
         end do
      end do
      do e = 1, nv
         do m = 1, na
            do k = 1, na
               do a = 1, nv
                  mat_v(a, k + na*(m - 1 + na*(e - 1))) = t(k, m, a, e)
               end do
            end do
         end do
      end do
   end subroutine arrange_amplitudes

   ! ==========================================================================
   ! first-order skeletons
   ! ==========================================================================

   subroutine skeleton_matrices(st, h1, s1, fx, sx, error)
      !! `h1`, `s1` in the AO basis for the coupled-perturbed solve, and the MO
      !! images `f^(X) = C^T h1 C`, `S^(X) = C^T s1 C`, `(n, n, 3N)`
      type(ri_state_t), intent(in) :: st
      real(dp), allocatable, intent(out) :: h1(:, :, :, :), s1(:, :, :, :)
      real(dp), allocatable, intent(out) :: fx(:, :, :), sx(:, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dhf(:, :), hc(:, :, :), sa(:, :, :), half(:, :)
      integer :: a, k, x

      allocate (dhf(st%nao, st%nao), half(st%nao, st%no))
      call pic_gemm(st%c(:, 1:st%no), st%c(:, 1:st%no), dhf, transb="T", alpha=2.0_dp)
      call h1_contract(st%mol, dhf, h1, error, k_scale=st%kf, bounds=st%bounds)
      if (error%has_error()) return
      ! The Kohn-Sham mean field also keeps part of itself on the quadrature grid;
      ! `h1_contract` differentiates integrals only.
      if (st%dh) then
         call xc_potential_deriv(st%xc, st%mol, st%scf_density, h1, error)
         if (error%has_error()) return
      end if
      allocate (s1(st%nao, st%nao, 3, st%natm))
      do a = 1, st%natm
         call hcore_deriv_atom(st%mol, a, hc, error)
         call overlap_deriv_atom(st%mol, a, sa, error)
         if (error%has_error()) return
         h1(:, :, :, a) = h1(:, :, :, a) + hc
         s1(:, :, :, a) = sa
         deallocate (hc, sa)
      end do
      deallocate (half)
      allocate (fx(st%n, st%n, st%n3), sx(st%n, st%n, st%n3), half(st%nao, st%n))
      do a = 1, st%natm
         do k = 1, 3
            x = 3*(a - 1) + k
            call pic_gemm(h1(:, :, k, a), st%c, half)
            call pic_gemm(st%c, half, fx(:, :, x), transa="T")
            call pic_gemm(s1(:, :, k, a), st%c, half)
            call pic_gemm(st%c, half, sx(:, :, x), transa="T")
         end do
      end do
   end subroutine skeleton_matrices

   ! ==========================================================================
   ! three-index first-derivative integrals
   ! ==========================================================================

   subroutine deriv_3c_range(orb, aux, ia0, ia1, d3, error)
      !! `d(mu nu|P)/dR_{A,c}` for the atoms `ia0..ia1`, `(nao, nao, P, 3, nat)`
      !!
      !! The full nuclear derivative: minus the electron-coordinate integrals
      !! `ip1` for a bra or ket function on the atom and `ip2` for an auxiliary
      !! one. The `nu` term is the transpose of the `mu` one, so the `ip1` sweep
      !! writes only rows its own shell owns and the transpose is added once at
      !! the end; the `ip2` sweep follows and shares no element with it.
      type(czt_molecule_t), intent(in) :: orb, aux
      integer, intent(in) :: ia0, ia1
      real(dp), allocatable, intent(out) :: d3(:, :, :, :, :)
      type(error_t), intent(inout) :: error

      integer, allocatable :: bas(:, :), orb_atom(:), aux_atom(:)
      real(dp), allocatable :: env(:), buf(:), slab(:, :)
      integer :: shls(4)
      integer :: dummy, nbas_orb, nbas_aux, ish, jsh, ksh, di, dj, dk, io, jo, ko
      integer :: i, j, k, comp, ret, idx, mx, at, nat, kk

      if (error%has_error()) return
      nbas_orb = orb%nbas
      nbas_aux = aux%nbas
      nat = ia1 - ia0 + 1
      call build_df_shell_table(orb, aux, bas, env, dummy)
      allocate (orb_atom(nbas_orb), aux_atom(nbas_aux))
      do ish = 1, nbas_orb
         orb_atom(ish) = orb%bas(LIBCINT_ATOM_OF, ish) + 1
      end do
      do ksh = 1, nbas_aux
         aux_atom(ksh) = aux%bas(LIBCINT_ATOM_OF, ksh) + 1
      end do
      mx = max(max_block(orb), max_block(aux))
      allocate (d3(orb%nao, orb%nao, aux%nao, 3, nat))
      d3 = 0.0_dp

      !$omp parallel default(none) &
      !$omp    shared(d3, orb, aux, bas, env, dummy, nbas_orb, nbas_aux, mx, orb_atom, &
      !$omp           aux_atom, ia0, ia1) &
      !$omp    private(buf, shls, ish, jsh, ksh, di, dj, dk, io, jo, ko, i, j, k, comp, ret, &
      !$omp            idx, at)
      allocate (buf(mx**3*3))
      !$omp do collapse(2) schedule(dynamic)
      do ish = 1, nbas_orb
         do jsh = 1, nbas_orb
            if (orb_atom(ish) < ia0 .or. orb_atom(ish) > ia1) cycle
            at = orb_atom(ish) - ia0 + 1
            di = shell_dim(orb%cartesian, ish - 1, bas)
            io = orb%shell_offset(ish)
            dj = shell_dim(orb%cartesian, jsh - 1, bas)
            jo = orb%shell_offset(jsh)
            do ksh = 1, nbas_aux
               dk = shell_dim(orb%cartesian, nbas_orb + ksh - 1, bas)
               ko = aux%shell_offset(ksh)
               shls = [ish - 1, jsh - 1, nbas_orb + ksh - 1, dummy - 1]
               if (orb%cartesian) then
                  ret = libcint_3c2e_ip1_cart(buf, shls, orb%atm, orb%natm, bas, dummy, env)
               else
                  ret = libcint_3c2e_ip1_sph(buf, shls, orb%atm, orb%natm, bas, dummy, env)
               end if
               if (ret == 0) cycle
               do comp = 1, 3
                  do k = 1, dk
                     do j = 1, dj
                        do i = 1, di
                           idx = i + di*(j - 1 + dj*(k - 1 + dk*(comp - 1)))
                           d3(io + i, jo + j, ko + k, comp, at) = -buf(idx)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end do
      deallocate (buf)
      !$omp end parallel

      allocate (slab(orb%nao, orb%nao))
      do at = 1, nat
         do comp = 1, 3
            do kk = 1, aux%nao
               slab = d3(:, :, kk, comp, at)
               d3(:, :, kk, comp, at) = slab + transpose(slab)
            end do
         end do
      end do
      deallocate (slab)

      !$omp parallel default(none) &
      !$omp    shared(d3, orb, aux, bas, env, dummy, nbas_orb, nbas_aux, mx, aux_atom, ia0, ia1) &
      !$omp    private(buf, shls, ish, jsh, ksh, di, dj, dk, io, jo, ko, i, j, k, comp, ret, &
      !$omp            idx, at)
      allocate (buf(mx**3*3))
      !$omp do collapse(2) schedule(dynamic)
      do ish = 1, nbas_orb
         do jsh = 1, nbas_orb
            di = shell_dim(orb%cartesian, ish - 1, bas)
            io = orb%shell_offset(ish)
            dj = shell_dim(orb%cartesian, jsh - 1, bas)
            jo = orb%shell_offset(jsh)
            do ksh = 1, nbas_aux
               if (aux_atom(ksh) < ia0 .or. aux_atom(ksh) > ia1) cycle
               at = aux_atom(ksh) - ia0 + 1
               dk = shell_dim(orb%cartesian, nbas_orb + ksh - 1, bas)
               ko = aux%shell_offset(ksh)
               shls = [ish - 1, jsh - 1, nbas_orb + ksh - 1, dummy - 1]
               if (orb%cartesian) then
                  ret = libcint_3c2e_ip2_cart(buf, shls, orb%atm, orb%natm, bas, dummy, env)
               else
                  ret = libcint_3c2e_ip2_sph(buf, shls, orb%atm, orb%natm, bas, dummy, env)
               end if
               if (ret == 0) cycle
               do comp = 1, 3
                  do k = 1, dk
                     do j = 1, dj
                        do i = 1, di
                           idx = i + di*(j - 1 + dj*(k - 1 + dk*(comp - 1)))
                           d3(io + i, jo + j, ko + k, comp, at) = &
                              d3(io + i, jo + j, ko + k, comp, at) - buf(idx)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end do
      deallocate (buf)
      !$omp end parallel
   end subroutine deriv_3c_range

   subroutine contract_first_3c(orb, aux, c_act, c_vir, dens, pair, error)
      !! `pair(X, y) = sum_{mn P} D^y_{mn P} d(mn|P)/dX` for every nuclear
      !! coordinate `X`, with `D^y_{mn P} = sum_{ia} C_mi C_na dens(i,a,P,y)`
      !!
      !! One sweep over the first-derivative shell triples serves the whole batch
      !! of densities. For each auxiliary shell the AO slab of every density is
      !! built once, and each triple's `ip1` and `ip2` blocks are contracted with
      !! all of them as they are generated, so no density meets an integral it
      !! has to regenerate. The `nu` derivative reads the density transposed.
      type(czt_molecule_t), intent(in) :: orb, aux
      real(dp), intent(in) :: c_act(:, :), c_vir(:, :)   !! `(nao, na)`, `(nao, nv)`
      real(dp), intent(in) :: dens(:, :, :, :)           !! `(na, nv, P, ny)`
      real(dp), intent(out) :: pair(:, :)                !! `(3 natm, ny)`
      type(error_t), intent(inout) :: error

      integer, allocatable :: bas(:, :), orb_atom(:), aux_atom(:)
      real(dp), allocatable :: env(:), buf1(:), buf2(:), ploc(:, :), slab(:, :, :, :)
      real(dp), allocatable :: half(:, :), tmp(:, :)
      real(dp) :: s1, s2
      integer :: shls(4)
      integer :: dummy, nbas_orb, nbas_aux, ish, jsh, ksh, di, dj, dk, io, jo, ko
      integer :: i, j, k, comp, ret1, ret2, idx, mx, ny, y, nao, na, nv, natm

      if (error%has_error()) return
      nao = orb%nao
      na = size(c_act, 2)
      nv = size(c_vir, 2)
      ny = size(dens, 4)
      natm = orb%natm
      nbas_orb = orb%nbas
      nbas_aux = aux%nbas
      pair = 0.0_dp
      call build_df_shell_table(orb, aux, bas, env, dummy)
      allocate (orb_atom(nbas_orb), aux_atom(nbas_aux))
      do ish = 1, nbas_orb
         orb_atom(ish) = orb%bas(LIBCINT_ATOM_OF, ish) + 1
      end do
      do ksh = 1, nbas_aux
         aux_atom(ksh) = aux%bas(LIBCINT_ATOM_OF, ksh) + 1
      end do
      mx = max(max_block(orb), max_block(aux))
      allocate (half(nao, nv), tmp(nao, nao))

      do ksh = 1, nbas_aux
         dk = aux%shell_offset(ksh + 1) - aux%shell_offset(ksh)
         ko = aux%shell_offset(ksh)
         allocate (slab(nao, nao, dk, ny))
         do y = 1, ny
            do k = 1, dk
               call pic_gemm(c_act, dens(:, :, ko + k, y), half)
               call pic_gemm(half, c_vir, slab(:, :, k, y), transb="T")
            end do
         end do

         !$omp parallel default(none) &
         !$omp    shared(slab, orb, bas, env, dummy, nbas_orb, ksh, dk, ko, mx, ny, natm, &
         !$omp           orb_atom, aux_atom, pair) &
         !$omp    private(buf1, buf2, ploc, shls, ish, jsh, di, dj, io, jo, i, j, k, comp, &
         !$omp            ret1, ret2, idx, y, s1, s2)
         allocate (buf1(mx**3*3), buf2(mx**3*3), ploc(3*natm, ny))
         ploc = 0.0_dp
         !$omp do collapse(2) schedule(dynamic)
         do ish = 1, nbas_orb
            do jsh = 1, nbas_orb
               di = shell_dim(orb%cartesian, ish - 1, bas)
               io = orb%shell_offset(ish)
               dj = shell_dim(orb%cartesian, jsh - 1, bas)
               jo = orb%shell_offset(jsh)
               shls = [ish - 1, jsh - 1, nbas_orb + ksh - 1, dummy - 1]
               if (orb%cartesian) then
                  ret1 = libcint_3c2e_ip1_cart(buf1, shls, orb%atm, orb%natm, bas, dummy, env)
                  ret2 = libcint_3c2e_ip2_cart(buf2, shls, orb%atm, orb%natm, bas, dummy, env)
               else
                  ret1 = libcint_3c2e_ip1_sph(buf1, shls, orb%atm, orb%natm, bas, dummy, env)
                  ret2 = libcint_3c2e_ip2_sph(buf2, shls, orb%atm, orb%natm, bas, dummy, env)
               end if
               do y = 1, ny
                  do comp = 1, 3
                     s1 = 0.0_dp
                     s2 = 0.0_dp
                     do k = 1, dk
                        do j = 1, dj
                           do i = 1, di
                              idx = i + di*(j - 1 + dj*(k - 1 + dk*(comp - 1)))
                              if (ret1 /= 0) s1 = s1 + buf1(idx)* &
                                                  (slab(io + i, jo + j, k, y) &
                                                   + slab(jo + j, io + i, k, y))
                              if (ret2 /= 0) s2 = s2 + buf2(idx)*slab(io + i, jo + j, k, y)
                           end do
                        end do
                     end do
                     ploc(3*(orb_atom(ish) - 1) + comp, y) = &
                        ploc(3*(orb_atom(ish) - 1) + comp, y) - s1
                     ploc(3*(aux_atom(ksh) - 1) + comp, y) = &
                        ploc(3*(aux_atom(ksh) - 1) + comp, y) - s2
                  end do
               end do
            end do
         end do
         !$omp end do
         !$omp critical
         pair = pair + ploc
         !$omp end critical
         deallocate (buf1, buf2, ploc)
         !$omp end parallel
         deallocate (slab)
      end do
   end subroutine contract_first_3c

   ! ==========================================================================
   ! the perturbed ladder
   ! ==========================================================================

   subroutine perturbed_ladder(st, fx, sx, mo1, nab, tol, pair, orb, resp, error, ddip)
      !! Orbital response and density response, in stages that batch every
      !! Fock-type build over all `3N` perturbations
      !!
      !!     A  rotations `U^Y`, `d_Y f`, and the total derivatives of `T1(D)`
      !!     B  per batch of atoms: three-index derivatives, `d_Y t2`, `d_Y D`,
      !!        the fitted Lagrangian terms and the pair sum
      !!     C  `d_Y I'` at the unrelaxed density, the Z-vector right-hand side,
      !!        the core<->active rotation
      !!     D  one batched Z-vector solve, then the solved blocks
      !!     E  the skeleton Lagrangian's carriers and the orbital response
      type(ri_state_t), intent(inout) :: st
      real(dp), intent(in) :: fx(:, :, :), sx(:, :, :), mo1(:, :, :, :)
      integer, intent(in) :: nab
      real(dp), intent(in) :: tol
      real(dp), allocatable, intent(out) :: pair(:, :), orb(:, :), resp(:, :)
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: ddip(:, :)
         !! The perturbative term's dipole derivative, unscaled; a double hybrid's

      real(dp), allocatable :: u_all(:, :, :), df_all(:, :, :), t1f(:, :, :), work(:, :, :)
      real(dp), allocatable :: ddmat(:, :, :), kd_d0(:, :, :), kd_z(:, :, :), kd_drel(:, :, :)
      real(dp), allocatable :: mo1_half(:, :)
      real(dp), allocatable :: dt1_d0(:, :, :), dt1_z(:, :, :), dt1_drel(:, :, :), dt1_p(:, :, :)
      real(dp), allocatable :: t1x_d0(:, :, :), t1x_z(:, :, :), t1x_drel(:, :, :), t1x_p(:, :, :)
      real(dp), allocatable :: ddg_all(:, :, :), lagc_pert(:, :, :), lagc_carr(:, :, :)
      real(dp), allocatable :: di_all(:, :, :), ddrel_all(:, :, :), rhs1(:, :, :)
      real(dp), allocatable :: t1dd(:, :, :), dps_all(:, :, :), t1dps(:, :, :)
      real(dp), allocatable :: rhs(:, :, :), zx(:, :, :), zm_all(:, :, :), t1zm(:, :, :)
      real(dp), allocatable :: ip2c(:, :, :), d3(:, :, :, :, :), d3b(:, :, :, :), d2b(:, :, :)
      real(dp), allocatable :: p3(:, :), dip(:, :), tm(:, :), dpco(:, :), pf_all(:, :, :)
      real(dp), allocatable :: ip_all(:, :, :), t1pf(:, :, :), xt(:, :), it(:, :), pof(:, :)
      real(dp), allocatable :: pfm(:, :), tmp(:, :)
      integer, allocatable :: offsets(:), counts(:), own_aux(:)
      real(dp) :: acc
      integer :: n, no, nf, na, nv, np, n3, natm
      integer :: y, ax, k, i, a, c, j, p, q, ia0, ia1, nb, ny, iy, at, comp, y0
      logical :: frozen

      n = st%n
      no = st%no
      nf = st%nf
      na = st%na
      nv = st%nv
      np = st%np
      n3 = st%n3
      natm = st%natm
      frozen = nf > 0

      ! ---- A: rotations and Fock responses ----------------------------------
      allocate (u_all(n, n, n3), df_all(n, n, n3), t1f(n, n, n3), work(n, n, n3))
      work = 0.0_dp
      do y = 1, n3
         ax = (y - 1)/3 + 1
         k = y - 3*(ax - 1)
         work(1:no, 1:no, y) = -sx(1:no, 1:no, y)
         work(no + 1:n, 1:no, y) = mo1(no + 1:n, :, k, ax)
         work(1:no, no + 1:n, y) = transpose(mo1(no + 1:n, :, k, ax))
      end do
      call t1_batch(st, work, t1f, error)
      if (error%has_error()) return
      if (st%dh) then
         ! The reference density's matrix response, `2 C (mx + mx^T) C^T`, which is
         ! `2 C work C^T`: what the kernel's coefficients move with.
         allocate (ddmat(st%nao, st%nao, n3), mo1_half(st%nao, n))
         do y = 1, n3
            call pic_gemm(st%c, work(:, :, y), mo1_half)
            call pic_gemm(mo1_half, st%c, ddmat(:, :, y), transb="T", alpha=2.0_dp)
         end do
         deallocate (mo1_half)
      end if
      do y = 1, n3
         ax = (y - 1)/3 + 1
         k = y - 3*(ax - 1)
         u_all(:, :, y) = 0.0_dp
         u_all(:, 1:no, y) = mo1(:, :, k, ax)
         u_all(no + 1:n, no + 1:n, y) = -0.5_dp*sx(no + 1:n, no + 1:n, y)
         u_all(1:no, no + 1:n, y) = -sx(1:no, no + 1:n, y) - transpose(mo1(no + 1:n, :, k, ax))
         if (frozen) then
            do i = nf + 1, no
               do c = 1, nf
                  acc = fx(c, i, y) - sx(c, i, y)*st%eps(i) + t1f(c, i, y)
                  u_all(c, i, y) = -acc/(st%eps(c) - st%eps(i))
                  u_all(i, c, y) = -sx(c, i, y) - u_all(c, i, y)
               end do
            end do
         end if
         do q = 1, n
            do p = 1, n
               df_all(p, q, y) = fx(p, q, y) + u_all(p, q, y)*st%eps(p) &
                                 + u_all(q, p, y)*st%eps(q) + t1f(p, q, y)
            end do
         end do
      end do
      deallocate (t1f, work)

      allocate (t1x_d0(n, n, n3), t1x_z(n, n, n3), t1x_drel(n, n, n3))
      call t1x_stack(st, st%d0, t1x_d0, error, kd=kd_d0)
      if (error%has_error()) return
      call t1x_stack(st, st%zsym, t1x_z, error, kd=kd_z)
      if (error%has_error()) return
      call t1x_stack(st, st%drel, t1x_drel, error, kd=kd_drel)
      if (error%has_error()) return
      allocate (dt1_d0(n, n, n3), dt1_z(n, n, n3), dt1_drel(n, n, n3))
      if (st%dh) then
         call total_derivative_t1(st, u_all, st%d0, st%t1_d0, t1x_d0, dt1_d0, error, &
                                  kd=kd_d0, ddmat=ddmat)
         if (error%has_error()) return
         call total_derivative_t1(st, u_all, st%zsym, st%t1_z, t1x_z, dt1_z, error, &
                                  kd=kd_z, ddmat=ddmat)
         if (error%has_error()) return
         call total_derivative_t1(st, u_all, st%drel, st%t1_drel, t1x_drel, dt1_drel, error, &
                                  kd=kd_drel, ddmat=ddmat)
         if (error%has_error()) return
         deallocate (kd_d0, kd_z, kd_drel, ddmat)
      else
         call total_derivative_t1(st, u_all, st%d0, st%t1_d0, t1x_d0, dt1_d0, error)
         call total_derivative_t1(st, u_all, st%zsym, st%t1_z, t1x_z, dt1_z, error)
         call total_derivative_t1(st, u_all, st%drel, st%t1_drel, t1x_drel, dt1_drel, error)
         if (error%has_error()) return
      end if
      if (frozen) then
         allocate (t1x_p(n, n, n3), dt1_p(n, n, n3))
         call t1x_stack(st, st%psym, t1x_p, error)
         if (error%has_error()) return
         call total_derivative_t1(st, u_all, st%psym, st%t1_p, t1x_p, dt1_p, error)
         if (error%has_error()) return
      end if

      ! ---- B: per batch of atoms ---------------------------------------------
      call two_centre_deriv(st%aux, ip2c)
      allocate (offsets(natm), counts(natm), own_aux(np))
      call atom_ao_blocks(st%aux, offsets, counts)
      do at = 1, natm
         own_aux(offsets(at) + 1:offsets(at) + counts(at)) = at
      end do
      allocate (ddg_all(n, n, n3), lagc_pert(n, n, n3), lagc_carr(n, n, n3), pair(n3, n3))
      pair = 0.0_dp
      ia0 = 1
      do while (ia0 <= natm)
         ia1 = min(ia0 + nab - 1, natm)
         nb = ia1 - ia0 + 1
         ny = 3*nb
         call deriv_3c_range(st%mol, st%aux, ia0, ia1, d3, error)
         if (error%has_error()) return
         allocate (d3b(na, nv, np, ny), d2b(np, np, ny), p3(n3, ny))
         do iy = 1, ny
            at = ia0 + (iy - 1)/3
            comp = iy - 3*((iy - 1)/3)
            y = 3*(at - 1) + comp
            call member_quantities(st, d3(:, :, :, comp, at - ia0 + 1), u_all(:, :, y), &
                                   df_all(:, :, y), ip2c(:, :, comp), own_aux, at, &
                                   lagc_carr(:, :, y), lagc_pert(:, :, y), &
                                   ddg_all(:, :, y), d3b(:, :, :, iy), d2b(:, :, iy))
         end do
         deallocate (d3)
         d3b = 4.0_dp*d3b
         call contract_first_3c(st%mol, st%aux, st%c(:, st%nf + 1:st%no), st%c(:, st%no + 1:n), &
                                d3b, p3, error)
         if (error%has_error()) return
         y0 = 3*(ia0 - 1)
         do iy = 1, ny
            pair(:, y0 + iy) = p3(:, iy)
            ! The two-index half, -2 D2 . J^(X), with J^(X)(P,Q) the minus of
            ! the electron-coordinate derivative on whichever function X moves.
            do comp = 1, 3
               do at = 1, natm
                  acc = 0.0_dp
                  do q = 1, np
                     do p = offsets(at) + 1, offsets(at) + counts(at)
                        acc = acc + (d2b(p, q, iy) + d2b(q, p, iy))*ip2c(p, q, comp)
                     end do
                  end do
                  pair(3*(at - 1) + comp, y0 + iy) = pair(3*(at - 1) + comp, y0 + iy) &
                                                     + 2.0_dp*acc
               end do
            end do
         end do
         deallocate (d3b, d2b, p3)
         ia0 = ia1 + 1
      end do
      deallocate (ip2c, offsets, counts, own_aux)

      ! ---- C: the Lagrangian at the unrelaxed density -------------------------
      allocate (t1dd(n, n, n3), rhs1(no, nv, n3), di_all(n, n, n3), ddrel_all(n, n, n3))
      call t1_batch(st, ddg_all, t1dd, error)
      if (error%has_error()) return
      allocate (dip(n, n), tm(n, n), tmp(n, n))
      if (frozen) allocate (dps_all(n, n, n3), dpco(nf, na))
      do y = 1, n3
         call lagrangian_perturbed(st, df_all(:, :, y), st%d0, ddg_all(:, :, y), &
                                   lagc_pert(:, :, y), t1dd(:, :, y) + dt1_d0(:, :, y), dip)
         do a = 1, nv
            do i = 1, no
               rhs1(i, a, y) = dip(i, no + a) - dip(no + a, i) - dt1_z(no + a, i, y)
            end do
         end do
         ! - z df_vv^T + df_oo z
         call pic_gemm(st%z, df_all(no + 1:n, no + 1:n, y), tm(1:no, 1:nv), transb="T", &
                       alpha=-1.0_dp)
         rhs1(:, :, y) = rhs1(:, :, y) + tm(1:no, 1:nv)
         call pic_gemm(df_all(1:no, 1:no, y), st%z, tm(1:no, 1:nv))
         rhs1(:, :, y) = rhs1(:, :, y) + tm(1:no, 1:nv)
         if (frozen) then
            do i = 1, na
               do c = 1, nf
                  acc = dip(c, nf + i) - dip(nf + i, c)
                  acc = acc - dot_product(df_all(c, 1:nf, y), st%pco(:, i))
                  acc = acc + dot_product(st%pco(c, :), df_all(nf + 1:no, nf + i, y))
                  dpco(c, i) = acc/(st%eps(c) - st%eps(nf + i))
               end do
            end do
            dps_all(:, :, y) = 0.0_dp
            do i = 1, na
               do c = 1, nf
                  dps_all(nf + i, c, y) = dpco(c, i)
                  dps_all(c, nf + i, y) = dpco(c, i)
               end do
            end do
         end if
      end do
      if (frozen) then
         allocate (t1dps(n, n, n3))
         call t1_batch(st, dps_all, t1dps, error)
         if (error%has_error()) return
         do y = 1, n3
            do a = 1, nv
               do i = 1, no
                  rhs1(i, a, y) = rhs1(i, a, y) + t1dps(no + a, i, y) + dt1_p(no + a, i, y)
               end do
            end do
            ddg_all(:, :, y) = ddg_all(:, :, y) + dps_all(:, :, y)
            t1dd(:, :, y) = t1dd(:, :, y) + t1dps(:, :, y)
         end do
         deallocate (t1dps, dps_all, dpco)
      end if
      do y = 1, n3
         call lagrangian_perturbed(st, df_all(:, :, y), st%drel, ddg_all(:, :, y), &
                                   lagc_pert(:, :, y), t1dd(:, :, y) + dt1_drel(:, :, y), dip)
         di_all(:, :, y) = dip
         ddrel_all(:, :, y) = ddg_all(:, :, y)
      end do
      deallocate (t1dd, dip, tm, tmp, lagc_pert, ddg_all, dt1_d0, dt1_z, dt1_drel)
      if (frozen) deallocate (dt1_p, t1x_p)

      ! ---- D: one batched Z-vector solve --------------------------------------
      allocate (rhs(nv, no, n3))
      do y = 1, n3
         rhs(:, :, y) = transpose(rhs1(:, :, y))
      end do
      call solve_response(st, rhs, zx, tol, error)
      if (error%has_error()) return
      deallocate (rhs, rhs1)
      allocate (zm_all(n, n, n3), t1zm(n, n, n3))
      zm_all = 0.0_dp
      do y = 1, n3
         do a = 1, nv
            do i = 1, no
               zm_all(i, no + a, y) = zx(a, i, y)
               zm_all(no + a, i, y) = zx(a, i, y)
            end do
         end do
      end do
      deallocate (zx)
      call t1_batch(st, zm_all, t1zm, error)
      if (error%has_error()) return
      do y = 1, n3
         ddrel_all(:, :, y) = ddrel_all(:, :, y) + zm_all(:, :, y)
         do q = 1, n
            do p = 1, n
               di_all(p, q, y) = di_all(p, q, y) - 0.5_dp*st%eps(p)*(zm_all(p, q, y) + zm_all(q, p, y))
            end do
         end do
         do q = 1, no
            do p = 1, n
               di_all(p, q, y) = di_all(p, q, y) - 0.5_dp*(t1zm(p, q, y) + t1zm(q, p, y))
            end do
         end do
      end do
      deallocate (zm_all, t1zm)

      ! ---- E: the skeleton Lagrangian's carriers and the orbital response ------
      allocate (ip_all(n, n, n3))
      do y = 1, n3
         call pic_gemm(fx(:, :, y), st%drel, ip_all(:, :, y), transb="T")
         call pic_gemm(fx(:, :, y), st%drel, ip_all(:, :, y), transa="T", beta=1.0_dp)
         call pic_gemm(sx(:, :, y), st%ip, ip_all(:, :, y), transb="T", beta=1.0_dp)
         call pic_gemm(sx(:, :, y), st%ip, ip_all(:, :, y), transa="T", beta=1.0_dp)
         ip_all(:, :, y) = ip_all(:, :, y) + lagc_carr(:, :, y)
         do q = 1, no
            do p = 1, n
               ip_all(p, q, y) = ip_all(p, q, y) + t1x_drel(p, q, y) + t1x_drel(q, p, y)
            end do
         end do
         ip_all(:, :, y) = -0.5_dp*ip_all(:, :, y)
      end do
      deallocate (lagc_carr, t1x_drel, t1x_d0, t1x_z)
      if (frozen) then
         allocate (pf_all(n, n, n3), t1pf(n, n, n3), pof(no, no))
         pf_all = 0.0_dp
         do y = 1, n3
            pof = 0.0_dp
            do q = nf + 1, no
               do p = 1, nf
                  pof(p, q) = (ip_all(p, q, y) - ip_all(q, p, y))/(st%eps(p) - st%eps(q))
                  pof(q, p) = pof(p, q)
               end do
            end do
            pf_all(1:no, 1:no, y) = pof
         end do
         call t1_batch(st, pf_all, t1pf, error)
         if (error%has_error()) return
      end if
      allocate (orb(n3, n3), xt(n, n), it(n, n))
      do p = 1, n3
         ! carriers for X = p
         do q = 1, n
            do j = 1, n
               xt(j, q) = ip_all(q, j, p) - ip_all(j, q, p)
            end do
         end do
         it = ip_all(:, :, p)
         it(no + 1:n, 1:no) = transpose(ip_all(1:no, no + 1:n, p))
         if (frozen) then
            xt(no + 1:n, 1:no) = xt(no + 1:n, 1:no) + t1pf(no + 1:n, 1:no, p)
            do q = 1, no
               do j = 1, no
                  it(j, q) = it(j, q) - pf_all(j, q, p)*st%eps(q) - t1pf(j, q, p)
               end do
            end do
         end if
         do y = 1, n3
            ax = (y - 1)/3 + 1
            k = y - 3*(ax - 1)
            acc = 2.0_dp*sum(mo1(no + 1:n, :, k, ax)*xt(no + 1:n, 1:no)) &
                  + sum(sx(:, :, y)*it)
            if (frozen) acc = acc + sum(pf_all(:, :, p)*fx(:, :, y))
            orb(p, y) = acc
         end do
      end do
      deallocate (ip_all, xt, it)
      if (frozen) deallocate (pf_all, t1pf, pof)

      ! ---- the density response, contracted with the first-order skeletons ----
      allocate (resp(n3, n3))
      do y = 1, n3
         do p = 1, n3
            resp(p, y) = pair(p, y) + sum(ddrel_all(:, :, y)*fx(:, :, p)) &
                         + sum(di_all(:, :, y)*sx(:, :, p))
         end do
      end do
      resp = resp - pair
      if (present(ddip)) then
         call dipole_derivatives_of_density(st, u_all, ddrel_all, ddip, error)
         if (error%has_error()) return
      end if
      deallocate (ddrel_all, di_all, u_all, df_all)
   end subroutine perturbed_ladder

   subroutine dipole_derivatives_of_density(st, u_all, ddrel, ddip, error)
      !! `d mu_a / dR_(X,b) = - Tr(dD_AO/dR r_a) - Tr(D_AO d r_a/dR)` for the
      !! relaxed correlation density, `(3, 3N)`, unscaled
      !!
      !! The correlation term reaches the dipole through the relaxed density
      !! alone; there is no nuclear term, which belongs to the reference. The AO
      !! density's derivative moves the orbitals too:
      !! `dD_AO = C (dD_MO + U D + D U^T) C^T`.
      type(ri_state_t), intent(in) :: st
      real(dp), intent(in) :: u_all(:, :, :), ddrel(:, :, :)
      real(dp), allocatable, intent(out) :: ddip(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dip(:, :, :), dipd(:, :, :, :, :), drel_ao(:, :), dd_mo(:, :)
      real(dp), allocatable :: dd_ao(:, :), half(:, :)
      real(dp) :: origin(3)
      integer :: ia, a, b, y

      if (error%has_error()) return
      origin = 0.0_dp
      do ia = 1, st%natm
         origin = origin + st%mol%coords(:, ia)
      end do
      origin = origin/real(st%natm, dp)
      call multipole_matrices(st%mol, origin, 1, dip, error)
      if (error%has_error()) return
      call dipole_integral_derivatives(st%mol, origin, dipd, error)
      if (error%has_error()) return
      allocate (drel_ao(st%nao, st%nao), dd_ao(st%nao, st%nao), half(st%nao, st%n))
      allocate (dd_mo(st%n, st%n), ddip(3, st%n3))
      call pic_gemm(st%c, st%drel, half)
      call pic_gemm(half, st%c, drel_ao, transb="T")
      do ia = 1, st%natm
         do b = 1, 3
            y = 3*(ia - 1) + b
            dd_mo = ddrel(:, :, y)
            call pic_gemm(u_all(:, :, y), st%drel, dd_mo, beta=1.0_dp)
            call pic_gemm(st%drel, u_all(:, :, y), dd_mo, transb="T", beta=1.0_dp)
            call pic_gemm(st%c, dd_mo, half)
            call pic_gemm(half, st%c, dd_ao, transb="T")
            do a = 1, 3
               ddip(a, y) = -(sum(dd_ao*dip(:, :, a)) + sum(drel_ao*dipd(:, :, a, b, ia)))
            end do
         end do
      end do
   end subroutine dipole_derivatives_of_density

   subroutine total_derivative_t1(st, u_all, d, t1d, t1x, dt1, error, kd, ddmat)
      !! `d_Y T1(D)` at a fixed MO matrix `D`, for every `Y`
      !!
      !!     T1x(D)[Y] + T1(U D + D U^T) + U^T T1(D) + T1(D) U
      !!
      !! over a Kohn-Sham reference the exchange part of that is built from the
      !! exchange-only pieces, and the kernel's share is half of
      !! `ks_dl_pair_kernel`, which carries the same four rotations and
      !! skeleton derivative plus the term for the kernel's coefficients moving
      !! with the reference density (`ddmat`, `(nao, nao, 3N)`).
      type(ri_state_t), intent(in) :: st
      real(dp), intent(in) :: u_all(:, :, :), d(:, :), t1d(:, :), t1x(:, :, :)
      real(dp), intent(out) :: dt1(:, :, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: kd(:, :, :), ddmat(:, :, :)

      real(dp), allocatable :: m(:, :, :), tm(:, :, :), w(:, :), vk(:, :), t1hf(:, :)
      real(dp), allocatable :: wk(:, :)
      integer :: y, n

      n = st%n
      allocate (m(n, n, st%n3), tm(n, n, st%n3), w(n, n), t1hf(n, n))
      do y = 1, st%n3
         call pic_gemm(u_all(:, :, y), d, m(:, :, y))
         call pic_gemm(d, u_all(:, :, y), m(:, :, y), transb="T", beta=1.0_dp)
      end do
      call t1_batch(st, m, tm, error, kernel=.false.)
      if (error%has_error()) return
      t1hf = t1d
      if (st%dh) then
         allocate (vk(n, n), wk(n, n))
         call ks_pair_kernel(st%xc, st%mol, st%c, st%scf_density, d, vk, error)
         if (error%has_error()) return
         t1hf = t1d - 0.5_dp*vk
      end if
      do y = 1, st%n3
         call pic_gemm(u_all(:, :, y), t1hf, w, transa="T")
         dt1(:, :, y) = t1x(:, :, y) + tm(:, :, y) + w
         call pic_gemm(t1hf, u_all(:, :, y), dt1(:, :, y), beta=1.0_dp)
         if (st%dh) then
            dt1(:, :, y) = dt1(:, :, y) - 0.5_dp*kd(:, :, y)
            call ks_dl_pair_kernel(st%xc, st%mol, st%c, st%scf_density, d, vk, kd(:, :, y), &
                                   u_all(:, :, y), ddmat(:, :, y), wk, error)
            if (error%has_error()) return
            dt1(:, :, y) = dt1(:, :, y) + 0.5_dp*wk
         end if
      end do
   end subroutine total_derivative_t1

   subroutine lagrangian_perturbed(st, df, d, dd, lagc, t1sum, dip)
      !! `d_Y I'` at the density `d` whose response is `dd`
      !!
      !!     -1/2 [ d_Y f (D + D^T) + eps (dD + dD^T) + corr
      !!            + (T + T^T) for q occupied ],   T = T1(dD) + d_Y T1(D)
      type(ri_state_t), intent(in) :: st
      real(dp), intent(in) :: df(:, :), d(:, :), dd(:, :), lagc(:, :), t1sum(:, :)
      real(dp), intent(out) :: dip(:, :)

      real(dp), allocatable :: dsym(:, :)
      integer :: p, q

      allocate (dsym(st%n, st%n))
      dsym = d + transpose(d)
      call pic_gemm(df, dsym, dip)
      do q = 1, st%n
         do p = 1, st%n
            dip(p, q) = dip(p, q) + st%eps(p)*(dd(p, q) + dd(q, p)) + lagc(p, q)
         end do
      end do
      do q = 1, st%no
         do p = 1, st%n
            dip(p, q) = dip(p, q) + t1sum(p, q) + t1sum(q, p)
         end do
      end do
      dip = -0.5_dp*dip
   end subroutine lagrangian_perturbed

   subroutine member_quantities(st, ao, u, df, ip2c_c, own_aux, atom, lagc_carr, &
                                lagc_pert, ddg, d3y, d2y)
      !! Everything one perturbation needs from its three-index first-derivative
      !! integrals
      !!
      !! `ao` is `d(mn|P)/dR` for this coordinate. From it: the skeleton and total
      !! MO integrals, the fitted Lagrangian terms of the orbital carriers and of
      !! the density response, `d_Y t2` and `d_Y D`, and the two three-index
      !! densities the pair sum contracts -- `d_Y Gamma~ + J^-1 Theta` and the
      !! matching two-index one. The latter is what pass 1's cross terms and the
      !! response's `d_Y Gamma . R^(X)` both are, so they share one sweep.
      type(ri_state_t), intent(in) :: st
      real(dp), intent(in) :: ao(:, :, :), u(:, :), df(:, :), ip2c_c(:, :)
      integer, intent(in) :: own_aux(:)
      integer, intent(in) :: atom
      real(dp), intent(out) :: lagc_carr(:, :), lagc_pert(:, :), ddg(:, :)
      real(dp), intent(out) :: d3y(:, :, :), d2y(:, :)

      real(dp), allocatable :: jx(:, :), isk(:, :, :), di3(:, :, :), half(:, :), w(:, :)
      real(dp), allocatable :: diov2(:, :), cj(:, :), dk2(:, :), dt(:, :, :, :), du(:, :, :, :)
      real(dp), allocatable :: dg2(:, :), dg3(:, :, :), th(:, :, :), th2(:, :), skov2(:, :)
      real(dp), allocatable :: dgsk(:, :), dcy(:, :), gam2d(:, :), t3(:, :)
      real(dp), allocatable :: dtt_o(:, :), dut_o(:, :), dtp_v(:, :), dup_v(:, :)
      real(dp) :: acc
      integer :: n, no, nf, na, nv, np, p, q, i, j, a, b, c, k

      n = st%n
      no = st%no
      nf = st%nf
      na = st%na
      nv = st%nv
      np = st%np

      allocate (jx(np, np))
      do q = 1, np
         do p = 1, np
            jx(p, q) = 0.0_dp
            if (own_aux(p) == atom) jx(p, q) = jx(p, q) - ip2c_c(p, q)
            if (own_aux(q) == atom) jx(p, q) = jx(p, q) - ip2c_c(q, p)
         end do
      end do

      allocate (isk(n, n, np), di3(n, n, np), half(st%nao, n), w(n, n))
      do p = 1, np
         call pic_gemm(ao(:, :, p), st%c, half)
         call pic_gemm(st%c, half, isk(:, :, p), transa="T")
         call pic_gemm(u, st%imo(:, :, p), w, transa="T")
         call pic_gemm(st%imo(:, :, p), u, w, beta=1.0_dp)
         di3(:, :, p) = isk(:, :, p) + w
      end do
      call corr_lagrangian(st, isk, jx, lagc_carr)

      ! d_Y (ia|jb)_RI = dI c^T + c dI^T - c J^(Y) c^T
      allocate (diov2(na*nv, np), skov2(na*nv, np), cj(na*nv, np), dk2(na*nv, na*nv))
      do p = 1, np
         do a = 1, nv
            do i = 1, na
               diov2(i + na*(a - 1), p) = di3(nf + i, no + a, p)
               skov2(i + na*(a - 1), p) = isk(nf + i, no + a, p)
            end do
         end do
      end do
      call pic_gemm(st%cov2, jx, cj)
      call pic_gemm(diov2, st%cov2, dk2, transb="T")
      call pic_gemm(st%cov2, diov2, dk2, transb="T", beta=1.0_dp)
      call pic_gemm(cj, st%cov2, dk2, transb="T", alpha=-1.0_dp, beta=1.0_dp)

      allocate (dt(na, na, nv, nv))
      do b = 1, nv
         do a = 1, nv
            do j = 1, na
               do i = 1, na
                  acc = dk2(i + na*(a - 1), j + na*(b - 1))
                  do c = 1, nv
                     acc = acc + df(no + a, no + c)*st%t2(i, j, c, b) &
                           + df(no + b, no + c)*st%t2(i, j, a, c)
                  end do
                  do k = 1, na
                     acc = acc - df(nf + i, nf + k)*st%t2(k, j, a, b) &
                           - df(nf + j, nf + k)*st%t2(i, k, a, b)
                  end do
                  dt(i, j, a, b) = acc/(st%eps(nf + i) + st%eps(nf + j) &
                                        - st%eps(no + a) - st%eps(no + b))
               end do
            end do
         end do
      end do

      allocate (du(na, na, nv, nv))
      do b = 1, nv
         do a = 1, nv
            du(:, :, a, b) = 2.0_dp*dt(:, :, a, b) - dt(:, :, b, a)
         end do
      end do
      call arrange_amplitudes(st, dt, dtt_o, dtp_v)
      call arrange_amplitudes(st, du, dut_o, dup_v)
      ddg = 0.0_dp
      allocate (t3(na, na))
      call pic_gemm(dtt_o, st%ut_o, t3, transb="T", alpha=-2.0_dp)
      call pic_gemm(st%tt_o, dut_o, t3, transb="T", alpha=-2.0_dp, beta=1.0_dp)
      ddg(nf + 1:no, nf + 1:no) = t3
      deallocate (t3)
      allocate (t3(nv, nv))
      call pic_gemm(st%up_v, dtp_v, t3, transb="T", alpha=2.0_dp)
      call pic_gemm(dup_v, st%tp_v, t3, transb="T", alpha=2.0_dp, beta=1.0_dp)
      ddg(no + 1:n, no + 1:n) = t3
      deallocate (t3)

      ! d_Y Gamma~ = du c, and the Lagrangian term with the density's response in
      allocate (dg2(na*nv, np), dg3(na, nv, np))
      block
         real(dp), allocatable :: du2(:, :)
         allocate (du2(na*nv, na*nv))
         do b = 1, nv
            do a = 1, nv
               do j = 1, na
                  do i = 1, na
                     du2(i + na*(a - 1), j + na*(b - 1)) = du(i, j, a, b)
                  end do
               end do
            end do
         end do
         call pic_gemm(du2, st%cov2, dg2)
      end block
      dg3 = reshape(dg2, [na, nv, np])
      call corr_lagrangian(st, di3, jx, lagc_pert, dgam3=dg3)

      ! The skeleton change of the fitting coefficients at fixed cumulant
      allocate (th(na, nv, np), th2(na*nv, np), dgsk(na*nv, np), dcy(na*nv, np))
      call theta_term(st, reshape(skov2, [na, nv, np]), jx, th)
      th2 = reshape(th, [na*nv, np])
      call pic_gemm(th2, st%jinv, dgsk)
      cj = skov2 - cj
      call pic_gemm(cj, st%jinv, dcy)

      d3y = reshape(dg2 + dgsk, [na, nv, np])
      allocate (gam2d(na*nv, np))
      gam2d = reshape(st%gam3, [na*nv, np])
      call pic_gemm(dg2, st%cov2, d2y, transa="T")
      call pic_gemm(dgsk, st%cov2, d2y, transa="T", beta=1.0_dp)
      call pic_gemm(gam2d, dcy, d2y, transa="T", beta=1.0_dp)
   end subroutine member_quantities

   ! ==========================================================================
   ! pass 1
   ! ==========================================================================

   subroutine pass1_reference_share(st, hess, error)
      !! The reference operator's share of pass 1, against the relaxed density
      !! and the energy-weighted matrix
      !!
      !!     D_rel . h^(XY) + I' . S^(XY) + D_rel . (J - K/2)^(XY)[D_HF]
      !!
      !! Statement for statement `mp2_skeleton_hessian`'s one-electron deposits,
      !! without its reference density and without `gamma_eff_ao`: the Fock fold
      !! is `hess_2e_contract` with the relaxed density on the left and the
      !! reference's on the right, so no four-index array is formed.
      type(ri_state_t), intent(in) :: st
      real(dp), intent(inout) :: hess(:, :, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: drel_ao(:, :), w_ao(:, :), dhf(:, :), half(:, :), sym(:, :)
      real(dp), allocatable :: h2(:, :, :), hab(:, :, :), s2(:, :, :), sab(:, :, :)
      real(dp), allocatable :: tmp(:, :, :), r2(:, :, :), rab(:, :, :)
      real(dp), allocatable :: cross_c(:, :, :)
      integer, allocatable :: owner(:), offsets(:), counts(:)
      real(dp) :: total_c(3, 3)
      real(dp) :: wc, zc
      integer :: nao, natm, comp, a, b, iao, jao, ia, ja, c

      if (error%has_error()) return
      nao = st%nao
      natm = st%natm
      allocate (drel_ao(nao, nao), w_ao(nao, nao), dhf(nao, nao), half(nao, st%n))
      allocate (sym(st%n, st%n))
      call pic_gemm(st%c, st%drel, half)
      call pic_gemm(half, st%c, drel_ao, transb="T")
      sym = 0.5_dp*(st%ip + transpose(st%ip))
      call pic_gemm(st%c, sym, half)
      call pic_gemm(half, st%c, w_ao, transb="T")
      call pic_gemm(st%c(:, 1:st%no), st%c(:, 1:st%no), dhf, transb="T", alpha=2.0_dp)

      allocate (offsets(natm), counts(natm), owner(nao))
      call atom_ao_blocks(st%mol, offsets, counts)
      owner = 0
      do c = 1, natm
         owner(offsets(c) + 1:offsets(c) + counts(c)) = c
      end do

      call hess_1e_block(st%mol, HESS_KIN_II, h2, error)
      call hess_1e_block(st%mol, HESS_NUC_II, tmp, error)
      if (error%has_error()) return
      h2 = h2 + tmp
      deallocate (tmp)
      call hess_1e_block(st%mol, HESS_KIN_IJ, hab, error)
      call hess_1e_block(st%mol, HESS_NUC_IJ, tmp, error)
      if (error%has_error()) return
      hab = hab + tmp
      deallocate (tmp)
      call hess_1e_block(st%mol, HESS_OVLP_II, s2, error)
      call hess_1e_block(st%mol, HESS_OVLP_IJ, sab, error)
      if (error%has_error()) return

      do comp = 1, 9
         a = (comp - 1)/3 + 1
         b = comp - 3*(a - 1)
         do jao = 1, nao
            ja = owner(jao)
            do iao = 1, nao
               ia = owner(iao)
               wc = drel_ao(iao, jao)
               hess(a, b, ia, ia) = hess(a, b, ia, ia) + 2.0_dp*wc*h2(iao, jao, comp)
               hess(a, b, ia, ja) = hess(a, b, ia, ja) + 2.0_dp*wc*hab(iao, jao, comp)
               wc = w_ao(iao, jao)
               hess(a, b, ia, ia) = hess(a, b, ia, ia) + 2.0_dp*wc*s2(iao, jao, comp)
               hess(a, b, ia, ja) = hess(a, b, ia, ja) + 2.0_dp*wc*sab(iao, jao, comp)
            end do
         end do
      end do
      deallocate (h2, hab, s2, sab)

      allocate (cross_c(3, 3, natm))
      do c = 1, natm
         call hess_rinv_block(st%mol, c, HESS_RINV_II, r2, error)
         call hess_rinv_block(st%mol, c, HESS_RINV_IJ, rab, error)
         if (error%has_error()) return
         zc = st%mol%charges(c)
         cross_c = 0.0_dp
         do comp = 1, 9
            a = (comp - 1)/3 + 1
            b = comp - 3*(a - 1)
            do jao = 1, nao
               do iao = 1, nao
                  ia = owner(iao)
                  wc = r2(iao, jao, comp) + rab(iao, jao, comp)
                  cross_c(a, b, ia) = cross_c(a, b, ia) + drel_ao(iao, jao)*wc
               end do
            end do
         end do
         deallocate (r2, rab)
         total_c = 0.0_dp
         do ia = 1, natm
            total_c = total_c + cross_c(:, :, ia)
            do b = 1, 3
               do a = 1, 3
                  hess(a, b, ia, c) = hess(a, b, ia, c) + 2.0_dp*zc*cross_c(a, b, ia)
                  hess(a, b, c, ia) = hess(a, b, c, ia) + 2.0_dp*zc*cross_c(b, a, ia)
               end do
            end do
         end do
         hess(:, :, c, c) = hess(:, :, c, c) - 2.0_dp*zc*total_c
      end do
      deallocate (cross_c)

      call hess_2e_contract(st%mol, drel_ao, hess, error, k_scale=st%kf, density2=dhf, &
                            bounds=st%bounds)
      if (error%has_error()) return
      ! The reference operator's second derivative contains the exchange-
      ! correlation potential's, contracted against the relaxed correlation
      ! density: the one pass-1 term with no Hartree-Fock counterpart. Grid held
      ! fixed, as everywhere on this ladder.
      if (st%dh) call xc_potential_hessian(st%xc, st%mol, st%scf_density, drel_ao, hess, error)
   end subroutine pass1_reference_share

   subroutine pass1_fitted_share(st, hess, error)
      !! `4 Gamma^P_mn (mn|P)^(XY) - 2 gamma_PQ (P|Q)^(XY)`
      type(ri_state_t), intent(in) :: st
      real(dp), intent(inout) :: hess(:, :, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: g3(:, :, :), half(:, :), g2(:, :)
      integer :: p

      if (error%has_error()) return
      allocate (g3(st%nao, st%nao, st%np), half(st%nao, st%nv), g2(st%np, st%np))
      do p = 1, st%np
         call pic_gemm(st%c(:, st%nf + 1:st%no), st%gam3(:, :, p), half)
         call pic_gemm(half, st%c(:, st%no + 1:st%n), g3(:, :, p), transb="T", alpha=4.0_dp)
      end do
      g2 = -2.0_dp*st%gam2
      call hess_3c2c_contract(st%mol, st%aux, hess, error, gamma3=g3, gamma2=g2)
   end subroutine pass1_fitted_share

end module mqc_czt_ri_mp2_hessian
