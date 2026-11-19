module test_mqc_czt_df_hessian
   !! The analytic Hessian of a density-fitted Hartree-Fock energy.
   !!
   !! Held to the one thing it has to be: the derivative of the fitted
   !! gradient. A five-point stencil of `czt_scf_gradient` with the same
   !! auxiliary basis is a reference that shares no line with the Hessian --
   !! the gradient contracts first-derivative integrals against densities of
   !! its own, and knows nothing of the coupled-perturbed solve, the
   !! second-derivative integrals or the product term this adds.
   !!
   !! That matters most for the term a fitted Hessian can lose without looking
   !! wrong: the response of the fit itself to a moved auxiliary basis. Carried
   !! only to first order -- PySCF's `auxbasis_response = 1` -- the Hessian is
   !! still symmetric and still translationally invariant, and is 1e-4 of an
   !! element from the differenced gradient on water. The validation suite holds the
   !! same Hessian to PySCF's; this holds it to our own gradient, which is what
   !! a geometry optimisation and a differenced Hessian actually use.
   !!
   !! The geometry is off every axis and plane, so no block is zero for a
   !! reason that has nothing to do with the integrals.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use pic_blas_interfaces, only: pic_gemm
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, atom_ao_blocks, &
                                three_centre, two_centre, metric_inverse_sqrt
   use mqc_czt_hess_ints, only: hess_3c2c_contract
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_gradient, only: czt_scf_gradient, three_centre_deriv, two_centre_deriv
   use mqc_czt_hessian, only: df_rhf_hessian, rhf_hessian, response_controls_t
   use mqc_czt_df_hessian, only: df_metric_dropped, df_hessian_terms_t, &
                                 df_two_electron_derivatives, df_fit_to_mo, &
                                 df_response_mo
   use mqc_czt_cphf, only: fitted_potential_general
   implicit none
   private

   public :: collect_mqc_czt_df_hessian_tests

   integer, parameter :: N_ATOMS = 3
   integer, parameter :: N_COORD = 3*N_ATOMS
   integer, parameter :: N_OCC = 5
   integer, parameter :: N_ELECTRONS = 10
   real(dp), parameter :: WATER(3, N_ATOMS) = reshape( &
                          [0.05_dp, 0.02_dp, -0.03_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.04_dp, 1.4308_dp, 1.1078_dp], [3, N_ATOMS])
   character(len=2), parameter :: SYM(N_ATOMS) = ["O ", "H ", "H "]
   integer, parameter :: NUMBERS(N_ATOMS) = [8, 1, 1]

   !! Bohr. Small enough for the stencil's own error to sit near 1e-9, large
   !! enough that the SCF's convergence does not show through the division.
   real(dp), parameter :: STEP = 2.0e-3_dp
   !! Hartree/Bohr^2. The stencil, not the Hessian, sets this.
   real(dp), parameter :: TOL_FD = 3.0e-7_dp
   !! What separates a fitted Hessian from an exact one on this molecule, at
   !! the least: the two must differ, or the fit was never differentiated.
   real(dp), parameter :: MIN_FIT_EFFECT = 1.0e-5_dp
   !! On the SCF's commutator, so the density and not only the energy is
   !! converged.
   real(dp), parameter :: SCF_GRAD_TOL = 1.0e-10_dp
   !! Exact exchange kept by Hartree-Fock, by a typical hybrid, and by a pure
   !! functional.
   real(dp), parameter :: EXCHANGE_FRACTIONS(3) = [1.0_dp, 0.25_dp, 0.0_dp]
   !! `three_centre_deriv`'s selector for `(nabla mu nu|P)`
   integer, parameter :: DERIV_ORBITAL = 1
   !! and for `(mu nu|nabla P)`
   integer, parameter :: DERIV_AUXILIARY = 2
   !! Two orderings of the same sums over the same integrals: rounding only,
   !! through an inverse square root of the metric. Measured at 6e-11.
   real(dp), parameter :: TOL_SAME = 1.0e-9_dp

contains

   subroutine collect_mqc_czt_df_hessian_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("df_hessian_is_the_fitted_gradients_derivative", &
                               against_gradient_differences), &
                  new_unittest("df_hessian_with_d_functions", &
                               against_gradient_differences_polarised), &
                  new_unittest("df_hessian_is_not_the_exact_one", differs_from_exact), &
                  new_unittest("df_hessian_metric_keeps_every_mode", metric_is_whole), &
                  new_unittest("df_hessian_terms_match_the_dense_form", matches_dense), &
                  new_unittest("df_response_in_the_mo_basis", response_matches_general) &
                  ]
   end subroutine collect_mqc_czt_df_hessian_tests

   subroutine converge(geo, basis, aux_basis, mol, aux, scf, err)
      !! A fitted closed-shell SCF on water at `geo`, converged hard
      real(dp), intent(in) :: geo(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      type(czt_molecule_t), intent(out) :: mol, aux
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      call build_czt_molecule(NUMBERS, SYM, geo, basis, mol, err)
      if (err%has_error()) return
      call build_czt_molecule(NUMBERS, SYM, geo, aux_basis, aux, err)
      if (err%has_error()) return
      ! The commutator is named because the density is what gets used: left to
      ! derive itself from the energy tolerance it stops near 1e-7, and the
      ! stencil below then divides the difference between two such stops by a
      ! step of 2e-3. Without it this test passed on one thread and failed on
      ! four about one run in two, by up to 1e-6.
      call run_czt_rhf(mol, N_ELECTRONS, 300, 1.0e-14_dp, 1.0e-12_dp, .false., scf, err, &
                       aux=aux, grad_tol=SCF_GRAD_TOL)
   end subroutine converge

   subroutine fitted_gradient(geo, basis, aux_basis, grad, err)
      !! The fitted SCF gradient at one geometry, flattened
      real(dp), intent(in) :: geo(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      real(dp), intent(out) :: grad(:)
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      real(dp), allocatable :: g(:, :)

      call converge(geo, basis, aux_basis, mol, aux, scf, err)
      if (err%has_error()) return
      call czt_scf_gradient(mol, scf%density, orbitals=scf%orbitals, &
                            orbital_energies=scf%orbital_energies, n_occupied=N_OCC, &
                            gradient=g, error=err, aux=aux)
      if (err%has_error()) return
      grad = reshape(g, [size(grad)])
      call mol%destroy()
      call aux%destroy()
   end subroutine fitted_gradient

   subroutine to_matrix(h4, matrix)
      !! `(3, 3, natm, natm)` to `(3N, 3N)`, row `3*(A-1)+a`, column `3*(B-1)+b`
      real(dp), intent(in) :: h4(:, :, :, :)
      real(dp), intent(out) :: matrix(:, :)

      integer :: ia, ib, a, b

      do ib = 1, size(h4, 4)
         do ia = 1, size(h4, 3)
            do b = 1, 3
               do a = 1, 3
                  matrix(3*(ia - 1) + a, 3*(ib - 1) + b) = h4(a, b, ia, ib)
               end do
            end do
         end do
      end do
   end subroutine to_matrix

   subroutine differences_case(error, basis, aux_basis)
      !! The analytic Hessian against a five-point stencil of the fitted gradient
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: basis, aux_basis

      real(dp), parameter :: SHIFT(4) = [2.0_dp, 1.0_dp, -1.0_dp, -2.0_dp]
      real(dp), parameter :: WEIGHT(4) = [-1.0_dp, 8.0_dp, -8.0_dp, 1.0_dp]
      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: hess4(:, :, :, :)
      real(dp) :: geo(3, N_ATOMS), g(N_COORD, 4), fd(N_COORD, N_COORD)
      real(dp) :: analytic(N_COORD, N_COORD)
      real(dp) :: dev
      integer :: y, k, atom, comp

      call converge(WATER, basis, aux_basis, mol, aux, scf, err)
      call check(error,.not. err%has_error(), "the SCF must converge: "//err%get_full_trace())
      if (allocated(error)) return
      call df_rhf_hessian(mol, aux, scf%density, scf%orbitals, scf%orbital_energies, &
                          N_OCC, hess4, err, &
                          controls=response_controls_t(tol=1.0e-13_dp))
      call check(error,.not. err%has_error(), "the Hessian must run: "//err%get_full_trace())
      if (allocated(error)) return
      call to_matrix(hess4, analytic)

      call check(error, maxval(abs(analytic - transpose(analytic))) < 1.0e-9_dp, &
                 "the fitted Hessian is not symmetric")
      if (allocated(error)) return

      do y = 1, N_COORD
         atom = (y - 1)/3 + 1
         comp = y - 3*(atom - 1)
         do k = 1, 4
            geo = WATER
            geo(comp, atom) = geo(comp, atom) + SHIFT(k)*STEP
            call fitted_gradient(geo, basis, aux_basis, g(:, k), err)
            call check(error,.not. err%has_error(), "the displaced gradient failed")
            if (allocated(error)) return
         end do
         fd(:, y) = (WEIGHT(1)*g(:, 1) + WEIGHT(2)*g(:, 2) + WEIGHT(3)*g(:, 3) &
                     + WEIGHT(4)*g(:, 4))/(12.0_dp*STEP)
      end do
      fd = 0.5_dp*(fd + transpose(fd))
      dev = maxval(abs(analytic - fd))
      call check(error, dev < TOL_FD, "the fitted Hessian is not the derivative of the "// &
                 "fitted gradient")
      if (allocated(error)) write (*, "(a,es10.3)") "largest deviation ", dev
      call mol%destroy()
      call aux%destroy()
   end subroutine differences_case

   subroutine against_gradient_differences(error)
      !! s and p functions only, fitted in a set made for Coulomb and exchange
      type(error_type), allocatable, intent(out) :: error

      call differences_case(error, "6-31g", "cc-pvtz-jkfit")
   end subroutine against_gradient_differences

   subroutine against_gradient_differences_polarised(error)
      !! d functions on oxygen, and a different auxiliary family
      type(error_type), allocatable, intent(out) :: error

      call differences_case(error, "cc-pvdz", "cc-pvdz-rifit")
   end subroutine against_gradient_differences_polarised

   subroutine differs_from_exact(error)
      !! A fitted Hessian that equalled the exact one would be the exact one
      !! with fitted orbitals, which is the mistake of differentiating the
      !! wrong energy. The two are a few 1e-4 apart here.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: fitted(:, :, :, :), exact(:, :, :, :)

      call converge(WATER, "6-31g", "cc-pvdz-rifit", mol, aux, scf, err)
      call check(error,.not. err%has_error(), "the SCF must converge")
      if (allocated(error)) return
      call df_rhf_hessian(mol, aux, scf%density, scf%orbitals, scf%orbital_energies, &
                          N_OCC, fitted, err)
      call check(error,.not. err%has_error(), "the fitted Hessian must run")
      if (allocated(error)) return
      call rhf_hessian(mol, scf%density, scf%orbitals, scf%orbital_energies, N_OCC, &
                       exact, err)
      call check(error,.not. err%has_error(), "the exact Hessian must run")
      if (allocated(error)) return

      call check(error, maxval(abs(fitted - exact)) > MIN_FIT_EFFECT, &
                 "the fitted Hessian is indistinguishable from the exact one")
      call mol%destroy()
      call aux%destroy()
   end subroutine differs_from_exact

   subroutine metric_is_whole(error)
      !! The auxiliary sets used above invert without dropping a mode, which is
      !! the condition under which the analytic route is taken at all.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: aux
      type(error_t) :: err
      integer :: n_dropped

      call build_czt_molecule(NUMBERS, SYM, WATER, "cc-pvtz-jkfit", aux, err)
      call check(error,.not. err%has_error(), "the auxiliary basis must build")
      if (allocated(error)) return
      call df_metric_dropped(aux, n_dropped, err)
      call check(error,.not. err%has_error(), "the metric must factor")
      if (allocated(error)) return
      call check(error, n_dropped == 0, "cc-pvtz-jkfit on water should drop no mode")
      call aux%destroy()
   end subroutine metric_is_whole

   subroutine dense_two_electron_derivatives(orb, aux, density, c_occ, terms, error, k_scale)
      !! `df_two_electron_derivatives` as the formulas read
      !!
      !! With `c = M^-1 V` and, per perturbation, `U^x = V^x - M^x c`:
      !!
      !!     (mn|ls)^x  = V^x_mn . c_ls + c_mn . U^x_ls
      !!     E^xy       = sum Gamma [ V^xy . c - 1/2 c . M^xy . c ]
      !!                + sum Gamma   U^x . M^-1 . U^y
      !!
      !! `V^x`, `M^x` and `U^x` are built whole for each perturbation and
      !! contracted exactly as written, `naux nao^3` each and unthreaded. This
      !! is the reference: the module's routine never forms any of the three,
      !! so the two share the integrals and nothing else. Same arguments, and
      !! the same `b`, `h1` and `hess` filled.
      type(czt_molecule_t), intent(in) :: orb, aux
      real(dp), intent(in) :: density(:, :)   !! Closed shell, carrying its factor of two
      real(dp), intent(in) :: c_occ(:, :)     !! (nao, n_occ)
      type(df_hessian_terms_t), intent(out) :: terms
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: k_scale
         !! Exact-exchange fraction; one (Hartree-Fock) when absent.

      real(dp), allocatable :: v(:, :), metric(:, :), mih(:, :), minv(:, :), c(:, :)
      real(dp), allocatable :: d(:), dvec(:), dc(:, :, :), dcd(:, :), t(:, :)
      real(dp), allocatable :: gamma3(:, :, :), gamma2(:, :)
      real(dp), allocatable :: ip1(:, :, :, :), ip2(:, :, :, :), ip2c(:, :, :)
      real(dp), allocatable :: vx(:, :, :), mx(:, :), u(:, :), ubar(:), jx(:)
      real(dp), allocatable :: kx(:, :), work(:, :), half(:, :), y(:, :)
      real(dp), allocatable :: ut(:, :), yt(:, :, :)
      integer, allocatable :: orb_first(:), orb_count(:), aux_first(:), aux_count(:)
      real(dp) :: kf
      integer :: nao, naux, natm, n_occ, n_dropped, ia, ja, x, yy, p, q, mu, nu, ip, jp
      integer :: m0, m1, a0, a1

      if (error%has_error()) return
      kf = 1.0_dp
      if (present(k_scale)) kf = k_scale
      nao = orb%nao
      naux = aux%nao
      natm = orb%natm
      n_occ = size(c_occ, 2)
      if (aux%natm /= natm) then
         call error%set(ERROR_VALIDATION, "the fitted Hessian needs the auxiliary "// &
                        "basis on the same atoms as the orbital one.")
         return
      end if

      ! ---- the fit ---------------------------------------------------------
      call three_centre(orb, aux, v)
      call two_centre(aux, metric)
      call metric_inverse_sqrt(metric, mih, error, n_dropped=n_dropped)
      if (error%has_error()) return
      if (n_dropped > 0) then
         call error%set(ERROR_VALIDATION, "the auxiliary metric drops a mode, and "// &
                        "-M^-1 dM M^-1 is the derivative of its inverse only when it "// &
                        "drops none.")
         return
      end if
      allocate (minv(naux, naux), c(nao*nao, naux), terms%b(nao*nao, naux))
      call pic_gemm(mih, mih, minv)
      call pic_gemm(v, minv, c)
      call pic_gemm(v, mih, terms%b)
      deallocate (metric, minv)

      ! ---- the two fixed densities of the second-derivative integrals ------
      ! gamma3 = D d_P - k/2 D c_P D, gamma2 = -1/2 [d d^T - k/2 Tr(D c_P D c_Q)].
      allocate (dvec(nao*nao), d(naux))
      dvec = reshape(density, [nao*nao])
      d = matmul(dvec, c)

      allocate (dc(nao, nao, naux), dcd(nao*nao, naux), work(nao, nao))
      do p = 1, naux
         call pic_gemm(density, reshape(c(:, p), [nao, nao]), dc(:, :, p))
         call pic_gemm(dc(:, :, p), density, work)
         dcd(:, p) = reshape(work, [nao*nao])
      end do
      allocate (t(naux, naux))
      call pic_gemm(dcd, c, t, transa="T")

      allocate (gamma3(nao, nao, naux), gamma2(naux, naux))
      do p = 1, naux
         gamma3(:, :, p) = d(p)*density - 0.5_dp*kf*reshape(dcd(:, p), [nao, nao])
      end do
      do q = 1, naux
         do p = 1, naux
            gamma2(p, q) = -0.5_dp*(d(p)*d(q) - 0.5_dp*kf*t(p, q))
         end do
      end do
      deallocate (dcd, t)

      allocate (terms%hess(3, 3, natm, natm))
      terms%hess = 0.0_dp
      call hess_3c2c_contract(orb, aux, terms%hess, error, gamma3=gamma3, gamma2=gamma2)
      if (error%has_error()) return
      deallocate (gamma3, gamma2)

      ! ---- first derivatives, one perturbation at a time -------------------
      call three_centre_deriv(orb, aux, DERIV_ORBITAL, ip1)
      call three_centre_deriv(orb, aux, DERIV_AUXILIARY, ip2)
      call two_centre_deriv(aux, ip2c)

      allocate (orb_first(natm), orb_count(natm), aux_first(natm), aux_count(natm))
      call atom_ao_blocks(orb, orb_first, orb_count)
      call atom_ao_blocks(aux, aux_first, aux_count)

      allocate (terms%h1(nao, nao, 3, natm))
      allocate (vx(nao, nao, naux), mx(naux, naux), u(nao*nao, naux), ubar(naux))
      allocate (jx(nao*nao), kx(nao, nao), half(n_occ, nao), y(n_occ*n_occ, naux))
      allocate (ut(naux, 3*natm), yt(n_occ*n_occ, naux, 3*natm))

      do ia = 1, natm
         m0 = orb_first(ia) + 1
         m1 = orb_first(ia) + orb_count(ia)
         a0 = aux_first(ia) + 1
         a1 = aux_first(ia) + aux_count(ia)
         do x = 1, 3
            ip = 3*(ia - 1) + x

            ! `V^x`. The library's nabla acts on the electron, so the derivative
            ! with respect to a centre is minus it, for whichever of the three
            ! functions sits on this atom.
            vx = 0.0_dp
            do p = 1, naux
               do nu = 1, nao
                  do mu = m0, m1
                     vx(mu, nu, p) = vx(mu, nu, p) - ip1(mu, nu, p, x)
                     vx(nu, mu, p) = vx(nu, mu, p) - ip1(mu, nu, p, x)
                  end do
               end do
            end do
            if (a1 >= a0) vx(:, :, a0:a1) = vx(:, :, a0:a1) - ip2(:, :, a0:a1, x)

            ! `M^x`, the same rule for the metric's two functions.
            mx = 0.0_dp
            do q = 1, naux
               do p = a0, a1
                  mx(p, q) = mx(p, q) - ip2c(p, q, x)
                  mx(q, p) = mx(q, p) - ip2c(p, q, x)
               end do
            end do

            ! `U^x = V^x - M^x c`.
            u = reshape(vx, [nao*nao, naux])
            call pic_gemm(c, mx, u, alpha=-1.0_dp, beta=1.0_dp)
            ubar = matmul(dvec, u)

            ! Coulomb: `J^x = V^x d + c (U^x : D)`.
            jx = matmul(reshape(vx, [nao*nao, naux]), d) + matmul(c, ubar)

            ! Exchange: `K^x = sum_P [ V^x_P (D c_P) + (D c_P)^T U^x_P ]`.
            kx = 0.0_dp
            do p = 1, naux
               call pic_gemm(vx(:, :, p), dc(:, :, p), kx, alpha=1.0_dp, beta=1.0_dp)
               call pic_gemm(dc(:, :, p), reshape(u(:, p), [nao, nao]), kx, transa="T", &
                             alpha=1.0_dp, beta=1.0_dp)
            end do
            terms%h1(:, :, x, ia) = reshape(jx, [nao, nao]) - 0.5_dp*kf*kx

            ! What the product term of `E^xy` needs of this perturbation, with
            ! `M^-1/2` already on it: `U^x : D`, and `U^x` over occupied pairs.
            ut(:, ip) = matmul(mih, ubar)
            do p = 1, naux
               call pic_gemm(c_occ, reshape(u(:, p), [nao, nao]), half, transa="T")
               y(:, p) = reshape(matmul(half, c_occ), [n_occ*n_occ])
            end do
            call pic_gemm(y, mih, yt(:, :, ip))
         end do
      end do

      ! ---- sum Gamma U^x M^-1 U^y ------------------------------------------
      ! `D = 2 C C^T`, so the exchange trace is four times the sum over occupied
      ! pairs, and with its -k/2 that is -2k.
      do ja = 1, natm
         do ia = 1, natm
            do yy = 1, 3
               do x = 1, 3
                  ip = 3*(ia - 1) + x
                  jp = 3*(ja - 1) + yy
                  terms%hess(x, yy, ia, ja) = terms%hess(x, yy, ia, ja) &
                                              + dot_product(ut(:, ip), ut(:, jp)) &
                                              - 2.0_dp*kf*sum(yt(:, :, ip)*yt(:, :, jp))
               end do
            end do
         end do
      end do
   end subroutine dense_two_electron_derivatives

   subroutine matches_dense(error)
      !! The routine the Hessian calls against the dense one it restates
      !!
      !! The dense form builds `V^x`, `M^x` and `U^x` whole for each
      !! perturbation and contracts them as the formulas read. The other never
      !! forms any of the three. They share the integrals and nothing else, so
      !! agreement here is every index convention and every sign of the
      !! rearrangement at once: the Fock derivative for each of nine
      !! perturbations, and all eighty-one second derivatives.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(df_hessian_terms_t) :: lean, dense
      integer :: k

      call converge(WATER, "cc-pvdz", "cc-pvdz-rifit", mol, aux, scf, err)
      call check(error,.not. err%has_error(), "the SCF must converge: "//err%get_full_trace())
      if (allocated(error)) return
      ! Hartree-Fock, a hybrid's fraction, and a pure functional's none: the
      ! exchange terms enter the two forms through different intermediates, so
      ! a fraction applied to one and not the other shows at 0.25 and a term
      ! that ignores it shows at zero.
      do k = 1, size(EXCHANGE_FRACTIONS)
         call df_two_electron_derivatives(mol, aux, scf%density, scf%orbitals(:, 1:N_OCC), &
                                          lean, err, k_scale=EXCHANGE_FRACTIONS(k))
         call dense_two_electron_derivatives(mol, aux, scf%density, &
                                             scf%orbitals(:, 1:N_OCC), dense, err, &
                                             k_scale=EXCHANGE_FRACTIONS(k))
         call check(error,.not. err%has_error(), "both must run: "//err%get_full_trace())
         if (allocated(error)) return

         write (*, "(a,f5.2,a,3es10.2)") "    k = ", EXCHANGE_FRACTIONS(k), &
            ", lean - dense, max |b|, |h1|, |hess|: ", &
            maxval(abs(lean%b - dense%b)), maxval(abs(lean%h1 - dense%h1)), &
            maxval(abs(lean%hess - dense%hess))
         call check(error, maxval(abs(lean%b - dense%b)) < TOL_SAME, &
                    "the fitted tensor differs between the two forms")
         if (allocated(error)) return
         call check(error, maxval(abs(lean%h1 - dense%h1)) < TOL_SAME, &
                    "the Fock derivative differs between the two forms")
         if (allocated(error)) return
         call check(error, maxval(abs(lean%hess - dense%hess)) < TOL_SAME, &
                    "the fixed-density second derivative differs between the two forms")
         if (allocated(error)) return
         call lean%destroy()
         call dense%destroy()
      end do
      call mol%destroy()
      call aux%destroy()
   end subroutine matches_dense

   subroutine response_matches_general(error)
      !! The response operator's fitted build against the general-density one
      !!
      !! `df_response_mo` never forms a density. The reference does: it builds
      !! `D1 = 2 C (mo1 E^T + E mo1^T) C^T` in the atomic basis, hands it to
      !! `fitted_potential_general`, which assumes nothing about it, and
      !! transforms the result back. The trial vectors fill occupied and
      !! virtual rows alike, because the solver's do.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N_VEC = 2
      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(df_hessian_terms_t) :: terms
      real(dp), allocatable :: mo1(:, :, :), gmo(:, :, :), dens(:, :), g(:, :)
      real(dp), allocatable :: half(:, :), ref(:, :)
      real(dp) :: dev
      integer :: n_ao, n_mo, m, k

      call converge(WATER, "6-31g", "cc-pvtz-jkfit", mol, aux, scf, err)
      call check(error,.not. err%has_error(), "the SCF must converge: "//err%get_full_trace())
      if (allocated(error)) return
      call df_two_electron_derivatives(mol, aux, scf%density, scf%orbitals(:, 1:N_OCC), &
                                       terms, err)
      call check(error,.not. err%has_error(), "the terms must build: "//err%get_full_trace())
      if (allocated(error)) return
      call df_fit_to_mo(terms, scf%orbitals)

      n_ao = size(scf%orbitals, 1)
      n_mo = size(scf%orbitals, 2)
      allocate (mo1(n_mo, N_OCC, N_VEC), gmo(n_mo, N_OCC, N_VEC), dens(n_ao, n_ao), &
                g(n_ao, n_ao), half(n_ao, N_OCC), ref(n_mo, N_OCC))
      call random_number(mo1)
      mo1 = mo1 - 0.5_dp
      dev = 0.0_dp
      do k = 1, size(EXCHANGE_FRACTIONS)
         call df_response_mo(terms%b_mo, mo1, gmo, k_scale=EXCHANGE_FRACTIONS(k))
         do m = 1, N_VEC
            half = 2.0_dp*matmul(scf%orbitals, mo1(:, :, m))
            dens = matmul(half, transpose(scf%orbitals(:, 1:N_OCC)))
            dens = dens + transpose(dens)
            call fitted_potential_general(terms%b, dens, g, k_scale=EXCHANGE_FRACTIONS(k))
            ref = matmul(transpose(scf%orbitals), matmul(g, scf%orbitals(:, 1:N_OCC)))
            dev = max(dev, maxval(abs(gmo(:, :, m) - ref)))
         end do
      end do
      write (*, "(a,es10.2)") "    MO response - general density: ", dev
      call check(error, dev < TOL_SAME, "the MO-basis response differs from the general build")
      call terms%destroy()
      call mol%destroy()
      call aux%destroy()
   end subroutine response_matches_general

end module test_mqc_czt_df_hessian

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_df_hessian, only: collect_mqc_czt_df_hessian_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_df_hessian", collect_mqc_czt_df_hessian_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
