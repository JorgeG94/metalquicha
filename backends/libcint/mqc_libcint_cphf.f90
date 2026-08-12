!! Coupled-perturbed Hartree-Fock, and the static polarizability from it
module mqc_libcint_cphf
   !! First-order orbital response to a one-electron perturbation.
   !!
   !! An SCF is variational, so a perturbation does not simply shift the
   !! orbitals by first-order perturbation theory: changing them changes the
   !! Fock operator, which changes them again. The self-consistency is what
   !! makes the equations *coupled*, and solving them is the difference between
   !! an uncoupled estimate that is wrong by tens of percent and a derivative
   !! that is exact.
   !!
   !! **The equations.** Writing the first-order orbitals as a rotation into the
   !! virtual space, `C_i^(1) = sum_a C_a U_ai`, the stationarity of the energy
   !! through first order gives, for a closed-shell reference,
   !!
   !!     (eps_a - eps_i) U_ai + sum_bj [4(ai|bj) - (ab|ij) - (aj|ib)] U_bj = -h_ai
   !!
   !! with `h_ai` the perturbation in the occupied-virtual MO block. Only that
   !! block appears: a rotation among occupied orbitals alone cannot change the
   !! density, which is the same fact `..._localize` relies on.
   !!
   !! **The two-electron term is a Fock build, and that is the whole trick.**
   !! Forming that bracket as a matrix is `n_occ^2 n_vir^2` and unnecessary.
   !! Assemble instead the symmetrized response density
   !!
   !!     Dt = C_vir U C_occ^T + (C_vir U C_occ^T)^T
   !!
   !! and contract it exactly as an SCF contracts a density, `G(D) = J - K/2`.
   !! Expanding `C_vir^T G(Dt) C_occ` term by term gives
   !!
   !!     2(ai|bj) - (ab|ij)/2 - (aj|ib)/2
   !!
   !! which is half the bracket, so the two-electron contribution is
   !! `2 C_vir^T G(Dt) C_occ` and one existing Fock build per iteration is all
   !! the machinery needed. The same substitution is what lets analytical
   !! Hessians and MP2 gradients reuse this later.
   !!
   !! **Why conjugate gradients.** The bracket above is symmetric -- swapping
   !! `(ai)` with `(bj)` maps each of the three terms to itself -- and it is
   !! positive definite whenever the SCF sits at a minimum rather than a saddle.
   !! So the natural solver is preconditioned CG with the orbital-energy
   !! differences as the preconditioner -- the same denominators the uncoupled
   !! sum-over-states expression uses, so the preconditioner is exactly "ignore
   !! the coupling" and CG spends its iterations on nothing else. It converges in
   !! a handful. If it ever fails to here, the reference is the suspect and not
   !! the solver: a broken-symmetry or saddle-point SCF is the usual reason, and
   !! the positive-definiteness guards below say so rather than iterating on.
   !!
   !! **Restricted Hartree-Fock only.** With a functional the response operator
   !! acquires the exchange-correlation kernel, the second derivative of the
   !! energy density, which is a different object from anything the SCF already
   !! builds -- so a hybrid handed to this module would silently return the
   !! Hartree-Fock response of the Kohn-Sham orbitals. The caller has to keep
   !! that straight; there is no functional argument to get wrong. EFP2, which
   !! is what wants this, is defined over a Hartree-Fock reference anyway.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_libcint_multipole, only: multipole_matrices
   use mqc_libcint_localize, only: boys_localize
   use mqc_libcint_rhf, only: build_fock
   use mqc_libcint_direct, only: build_fock_direct, schwarz_bounds, direct_stats_t
   implicit none
   private

   public :: cphf_solve
   public :: static_polarizability
   public :: distributed_polarizability

   !> CG iterations before giving up.
   integer, parameter :: DEFAULT_MAX_ITER = 200

   !> Convergence on the residual norm, relative to the right-hand side.
   real(dp), parameter :: DEFAULT_TOL = 1.0e-11_dp

contains

   subroutine cphf_solve(mol, orbitals, orbital_energies, n_occ, perturbations, &
                         response, error, max_iter, tol, iterations, in_core)
      !! Solve the coupled-perturbed equations for one or more perturbations
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)          !! MO coefficients, (n_ao, n_mo)
      real(dp), intent(in) :: orbital_energies(:)     !! (n_mo), Hartree
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: perturbations(:, :, :)
         !! One-electron perturbation operators in the **AO** basis,
         !! `(n_ao, n_ao, n_perturbations)`. AO rather than MO because every
         !! caller has them that way and the occupied-virtual transform is this
         !! routine's business.
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
         !! Present and true is the reference path, not the production one -- see
         !! the note where it is used.

      real(dp), allocatable :: eri(:, :, :, :), c_occ(:, :), c_vir(:, :), bounds(:, :)
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
      if (size(perturbations, 1) /= n_ao .or. size(perturbations, 2) /= n_ao) then
         call error%set(ERROR_VALIDATION, "CPHF: the perturbations are not n_ao square")
         return
      end if

      n_pert = size(perturbations, 3)
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

      ! Integral-direct by default. The stored tensor is n_ao^4 -- 2 GB for uracil
      ! in 6-31G*, where the response equations themselves need tens of MB -- so it
      ! is the storage that decides how large a fragment can be done, not anything
      ! about the solver. It is kept because it is the path validated against a
      ! dense exact solve, and is therefore what the direct build gets checked
      ! against, exactly as `run_libcint_rhf` keeps its own `in_core`.
      direct = .true.
      if (present(in_core)) direct = .not. in_core
      if (direct) then
         call schwarz_bounds(mol, bounds, error)
         if (error%has_error()) return
         allocate (eri(0, 0, 0, 0))
      else
         call mol%eris(eri)
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
         ! h_ai, the perturbation in the occupied-virtual block.
         call pic_gemm(perturbations(:, :, ipert), c_occ, work)
         call pic_gemm(c_vir, work, rhs, transa="T")
         rhs = -rhs

         target_norm = sqrt(sum(rhs*rhs))
         if (target_norm <= 0.0_dp) then
            response(:, :, ipert) = 0.0_dp
            cycle
         end if
         target_norm = target_norm*use_tol

         ! Preconditioned CG. The starting guess is zero rather than the
         ! uncoupled solution on purpose: with x = 0 the first iterate *is* the
         ! uncoupled answer, so nothing is lost, and the residual norm then
         ! starts at ||h|| and the convergence threshold means what it says.
         x = 0.0_dp
         r = rhs
         z = r/gaps
         p = z
         rz = sum(r*z)

         iter = 0
         do iter = 1, limit
            call response_operator(mol, direct, eri, bounds, zero_h, c_occ, c_vir, &
                                   gaps, p, ap, error)
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
      deallocate (eri, bounds, c_occ, c_vir, gaps, rhs, x, r, z, p, ap, work, zero_h)
   end subroutine cphf_solve

   subroutine response_operator(mol, direct, eri, bounds, zero_h, c_occ, c_vir, &
                                gaps, u, au, error)
      !! Apply the coupled-perturbed operator to a trial rotation
      type(libcint_molecule_t), intent(in) :: mol
      logical, intent(in) :: direct          !! Recompute integrals rather than store them
      real(dp), intent(in) :: eri(:, :, :, :)   !! Zero-sized when `direct`
      real(dp), intent(in) :: bounds(:, :)      !! Zero-sized when not `direct`
      real(dp), intent(in) :: zero_h(:, :)   !! Zero, so the build returns only G
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :)
      real(dp), intent(in) :: gaps(:, :)
      real(dp), intent(in) :: u(:, :)
      real(dp), intent(out) :: au(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dtilde(:, :), g(:, :), half(:, :), work(:, :)
      type(direct_stats_t) :: stats
      integer :: n_ao, n_occ

      n_ao = size(c_occ, 1)
      n_occ = size(c_occ, 2)
      allocate (dtilde(n_ao, n_ao), g(n_ao, n_ao), half(n_ao, n_occ), work(n_ao, n_occ))

      ! Dt = C_vir U C_occ^T + transpose. Symmetrized because `build_fock`
      ! contracts a symmetric density, and because the physical first-order
      ! density is symmetric -- the two halves are the orbital and its conjugate.
      call pic_gemm(c_vir, u, half)
      call pic_gemm(half, c_occ, dtilde, transb="T")
      dtilde = dtilde + transpose(dtilde)

      if (direct) then
         call build_fock_direct(mol, zero_h, dtilde, bounds, g, stats, error)
         if (error%has_error()) return
      else
         call build_fock(zero_h, eri, dtilde, g)
      end if

      call pic_gemm(g, c_occ, work)
      call pic_gemm(c_vir, work, au, transa="T")
      au = gaps*u + 2.0_dp*au

      deallocate (dtilde, g, half, work)
   end subroutine response_operator

   subroutine distributed_polarizability(mol, orbitals, orbital_energies, n_occ, &
                                         tensors, centroids, error, n_core, &
                                         max_iter, tol, iterations, in_core)
      !! One polarizability tensor per localized orbital, at its centroid
      !!
      !! What an effective fragment potential wants: not the molecule's total
      !! response but where in the molecule it happens, so that an induced dipole
      !! can be placed on each bond and lone pair and polarized by the local
      !! field. This is M5 of the EFP plan.
      !!
      !! **The decomposition is exact and needs no new physics.** The occupied
      !! index of `alpha_kl = -4 sum_ai h^k_ai U^l_ai` is summed over, and a
      !! rotation `W` among the occupied orbitals is orthogonal, so inserting
      !! `W W^T` changes nothing:
      !!
      !!     alpha_kl = -4 sum_a sum_i' (h^k W)_ai' (U^l W)_ai'
      !!
      !! and each `i'` is one localized orbital's share. Nothing is fitted or
      !! partitioned by distance -- the localization does all the work, which is
      !! why M2 had to come first.
      !!
      !! **Each tensor is asymmetric, and that is not a bug.** The contraction
      !! pairs perturbation `k` with response `l`, and for a single orbital those
      !! are different objects; only the sum over orbitals restores symmetry.
      !! GAMESS's reference potential shows exactly this -- 0.17 asymmetry per
      !! O-H bond tensor, 4.6e-6 in their sum -- so all nine components are
      !! returned and none is averaged away. Symmetrizing here would agree with
      !! the total and disagree with every reference file.
      type(libcint_molecule_t), intent(in) :: mol
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
         !! Occupied orbitals to leave out, lowest first. Default zero.
         !!
         !! GAMESS excludes them: its water potential has five occupied orbitals
         !! and four polarizable points, and the four tensors sum to 4.865822
         !! against 4.866587 for the whole occupied space. The missing 7.7e-4 is
         !! diagonal-only, which is what a spherical 1s should contribute, so
         !! nothing is wrong with either number -- they are answers to different
         !! questions. Matching a reference file means asking GAMESS's.
         !!
         !! The response itself is always solved over *all* occupied orbitals,
         !! because the core does polarize; only the distribution drops it. So the
         !! sum of what comes back is the response projected onto the orbitals
         !! kept, which is what makes the sum rule exact for any `n_core`.
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: iterations
      logical, intent(in), optional :: in_core

      real(dp), allocatable :: dip(:, :, :), u(:, :, :), localized(:, :)
      real(dp), allocatable :: s(:, :), w(:, :), u_loc(:, :, :), h_loc(:, :, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :), work(:, :), sc(:, :)
      integer :: n_ao, n_mo, n_vir, n_lmo, core, k, l, i, a
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

      ! The response over the whole occupied space -- see `n_core` above.
      call cphf_solve(mol, orbitals, orbital_energies, n_occ, dip, u, error, &
                      max_iter=max_iter, tol=tol, iterations=iterations, &
                      in_core=in_core)
      if (error%has_error()) return

      allocate (c_occ(n_ao, n_occ), c_vir(n_ao, n_vir))
      c_occ = orbitals(:, 1:n_occ)
      c_vir = orbitals(:, n_occ + 1:n_mo)

      ! W maps the canonical occupied orbitals onto the localized ones. It is an
      ! isometry rather than a rotation when a core is excluded, and `W W^T` is
      ! then the projector onto the orbitals kept -- which is exactly why the sum
      ! rule stays exact instead of becoming approximate.
      call mol%overlap(s)
      allocate (sc(n_ao, n_lmo), w(n_occ, n_lmo))
      call pic_gemm(s, localized, sc)
      call pic_gemm(c_occ, sc, w, transa="T")

      allocate (h_loc(n_vir, n_lmo, 3), u_loc(n_vir, n_lmo, 3))
      allocate (work(n_ao, n_lmo))
      do k = 1, 3
         ! h in the localized basis, straight from the localized orbitals rather
         ! than by transforming the canonical block -- fewer places to transpose.
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
                                    error, max_iter, tol, iterations, in_core)
      !! The dipole polarizability tensor, in Bohr^3
      !!
      !! `alpha_kl = d mu_k / d F_l`, from the response to the three components
      !! of a uniform field. The perturbation to the electronic Hamiltonian for
      !! a field `F` is `+ r . F`, since the electronic dipole operator is `-r`;
      !! the first-order density is `2 (C_vir U C_occ^T + transpose)`, and
      !! contracting it against `-r_k` gives
      !!
      !!     alpha_kl = -4 sum_ai (r_k)_ai U^l_ai
      !!
      !! which is `4 h^k A^-1 h^l` once `U = -A^-1 h`, and therefore positive
      !! definite. That sign is derived rather than fitted, and the finite-field
      !! check confirms it independently -- worth stating because an overall sign
      !! or a factor of two is the most common defect in a response property and
      !! the least visible, a wrong factor still looking like a polarizability.
      !!
      !! **The origin does not matter, and that is provable rather than
      !! conventional.** Shifting it adds a multiple of the overlap to each
      !! dipole matrix, and the occupied-virtual block of the overlap is exactly
      !! zero because the orbitals are orthogonal. So no origin argument is
      !! offered: there is nothing for the caller to choose.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), intent(out) :: alpha(3, 3)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: iterations
      logical, intent(in), optional :: in_core

      real(dp), allocatable :: dip(:, :, :), u(:, :, :), h(:, :), work(:, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :)
      integer :: n_ao, n_mo, n_vir, k, l

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_vir = n_mo - n_occ

      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, error)
      if (error%has_error()) return

      call cphf_solve(mol, orbitals, orbital_energies, n_occ, dip, u, error, &
                      max_iter=max_iter, tol=tol, iterations=iterations, &
                      in_core=in_core)
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

end module mqc_libcint_cphf
