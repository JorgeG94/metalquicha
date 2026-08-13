!! Coupled-perturbed Hartree-Fock, and the static polarizability from it
module mqc_libcint_cphf
   !! First-order orbital response to a one-electron perturbation.
   !!
   !! An SCF is variational, so a perturbation does not simply shift the orbitals:
   !! changing them changes the Fock operator, which changes them again. That
   !! self-consistency is what makes the equations *coupled*, and it is worth tens of
   !! percent against an uncoupled estimate.
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
   !! **Why conjugate gradients.** That bracket is symmetric and positive definite
   !! whenever the SCF sits at a minimum, so preconditioned CG fits, with the
   !! orbital-energy differences as the preconditioner -- the preconditioner is then
   !! exactly "ignore the coupling" and CG spends its handful of iterations on
   !! nothing else. If it fails to converge the reference is the suspect, usually a
   !! saddle point or broken symmetry, and the guards below say so rather than
   !! iterating on.
   !!
   !! **Restricted Hartree-Fock only.** With a functional the response operator gains
   !! the exchange-correlation kernel, which nothing the SCF already builds supplies,
   !! so a hybrid handed to this module would silently get the Hartree-Fock response
   !! of Kohn-Sham orbitals. There is no functional argument to get wrong; the caller
   !! has to keep it straight. EFP2 wants Hartree-Fock anyway.
   use iso_fortran_env, only: output_unit
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_getrf, pic_getrs
   use mqc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   use mqc_libcint_integrals, only: libcint_molecule_t, build_df_mo_block
   use mqc_libcint_multipole, only: multipole_matrices
   use mqc_libcint_localize, only: boys_localize
   use mqc_libcint_rhf, only: build_fock
   use mqc_libcint_direct, only: build_fock_direct, build_fock_direct_nosym, &
                                 build_fock_direct_many, &
                                 schwarz_bounds, direct_stats_t
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
   !> Exposed side by side so a check can hold both builds of the same matrices.
   public :: build_hessian, build_hessian_df

   !> CG iterations before giving up.
   type :: response_hessian_t
      !! A built response Hessian, so several blocks can share one build
      !!
      !! Only `(A-B)` and the product `(A-B)(A+B)` are kept: those are all the solve
      !! needs, and `(A+B)` on its own is not wanted again.
      logical :: ready = .false.
      real(dp), allocatable :: aminus(:, :)
      real(dp), allocatable :: product(:, :)
   contains
      procedure :: destroy => hessian_destroy
   end type response_hessian_t

   integer, parameter :: DEFAULT_MAX_ITER = 200

   !> Columns of the Hessian built per pass over the integrals. Large enough that
   !> integral generation is amortised, small enough that the batch of densities
   !> stays a sensible size -- at 115 orbitals, 64 of them is 60 MB.
   integer, parameter :: HESSIAN_CHUNK = 64

   !> Above this many orbitals, recompute the integrals rather than store them.
   !>
   !> Set from measurement, not from memory. `validation/bench_fock` times both paths
   !> per density, batched the way the Hessian build batches them, on water:
   !>
   !>      orbitals   batched direct   stored contraction
   !>            19        1.10e-04             2.07e-05
   !>            24        1.76e-04             7.89e-05
   !>            58        2.84e-03             8.09e-03
   !>           115        3.24e-02             2.46e-01
   !>
   !> so storing wins below about forty orbitals and loses badly above it -- by
   !> nearly eight times at 115. The direct build is screened, threaded and exploits
   !> permutational symmetry; the stored contraction is a quadruple loop running at a
   !> few GFlop/s. An earlier version of this threshold was a two gigabyte memory
   !> budget, about 126 orbitals, which sent everything between forty and that down
   !> the slower path.
   integer, parameter :: IN_CORE_MAX_ORBITALS = 40

   !> Even below that, refuse to store what will not fit.
   real(dp), parameter :: IN_CORE_LIMIT = 4.0e9_dp

   !> Convergence on the residual norm, relative to the right-hand side.
   real(dp), parameter :: DEFAULT_TOL = 1.0e-11_dp

   !> Casimir-Polder quadrature points, which is how many imaginary frequencies a
   !> potential is tabulated at.
   integer, parameter, public :: N_CASIMIR_POLDER = 12

contains

   subroutine hessian_destroy(self)
      !! Release a built Hessian
      class(response_hessian_t), intent(inout) :: self
      if (allocated(self%aminus)) deallocate (self%aminus)
      if (allocated(self%product)) deallocate (self%product)
      self%ready = .false.
   end subroutine hessian_destroy

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

   subroutine response_operator_minus(mol, direct, eri, bounds, zero_h, c_occ, c_vir, &
                                      gaps, u, au, error)
      !! Apply `A - B`, the other half of a frequency-dependent response
      !!
      !! Same shape as `response_operator` above and the same Fock build, on the
      !! *antisymmetrized* response density instead of the symmetrized one:
      !!
      !!     A + B   from   Dt = C_vir U C_occ^T + transpose
      !!     A - B   from   Dt = C_vir U C_occ^T - transpose
      !!
      !! The Coulomb term vanishes identically on the second, because the integral
      !! is symmetric in its ket pair while the density is not, so `A - B` comes out
      !! carrying only exchange -- `delta (eps_a - eps_i) - (ab|ij) + (aj|ib)`.
      !!
      !! It must go through `build_fock_direct_nosym`, not the fast build: that one
      !! folds three of the eightfold permutations into a multiplicity factor, which
      !! doubles them where an antisymmetric density needs them to cancel. See the
      !! note on that routine, and the test that asserts the fast one really is
      !! wrong here.
      type(libcint_molecule_t), intent(in) :: mol
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
      !! The twelve imaginary frequencies a potential is tabulated at
      !!
      !! Gauss-Legendre on `t` in (-1,1) under the substitution
      !! `nu = w0 (1 + t)/(1 - t)` with `w0 = 0.3`, which is the Casimir-Polder
      !! quadrature GAMESS uses. Confirmed against the reference potential: the
      !! twelve values it labels its dynamic blocks with are exactly these, from
      !! 0.002792 to 32.239080.
      real(dp) :: nu(N_CASIMIR_POLDER)
      real(dp), parameter :: W0 = 0.3_dp
      !> 12-point Gauss-Legendre abscissae on (-1,1), positive half mirrored.
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
                                     progress)
      !! `alpha(i nu)` at each imaginary frequency
      !!
      !! **Imaginary frequency is the friendly case.** The time-dependent equations
      !! couple the `+omega` and `-omega` blocks, and at real frequency the result
      !! is non-symmetric with poles at every excitation energy -- which is why
      !! GAMESS solves it with GMRES against an explicit Hessian on disk. On the
      !! imaginary axis the same equations fold into something real and positive
      !! definite. Writing `S = X + Y` and `D' = (X - Y)/i`,
      !!
      !!     (A+B) S + nu D' = -2h
      !!     (A-B) D' - nu S = 0
      !!
      !! and eliminating `D'` gives
      !!
      !!     [ (A+B) + nu^2 (A-B)^-1 ] S = -2h
      !!
      !! whose operator is a sum of positive definite pieces -- no poles anywhere on
      !! the axis, and *better* conditioned as `nu` grows rather than worse. So
      !! conjugate gradients still applies, with an inner solve against `A - B` for
      !! each outer product.
      !!
      !! `alpha_kl(i nu) = -2 sum_ai h^k_ai S^l_ai`, which at `nu = 0` reduces to
      !! the static `-4 sum h U` because `S = 2U` there. That limit is the check
      !! worth having: it must reproduce `static_polarizability` exactly, and it
      !! costs nothing to run.
      type(libcint_molecule_t), intent(in) :: mol
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
         !! Reuse a built Hessian, or take one back to reuse. It depends only on the
         !! reference -- not on the perturbations or the frequencies -- so the three
         !! blocks of a potential can share one build. Without this each of them
         !! rebuilds the same matrices: 55 seconds apiece at 115 orbitals.
      logical, intent(in), optional :: in_core
         !! Store the two-electron integrals instead of recomputing them each time
         !! the response operator is applied. Defaults to storing them whenever the
         !! tensor fits `IN_CORE_LIMIT`, because this solver applies that operator a
         !! great many times: every outer iteration runs an inner solve of its own,
         !! so the count is frequencies times perturbations times outer times inner.
         !! For a fragment-sized molecule recomputing all of it each time dominates
         !! everything else the potential costs.
      real(dp), intent(in), optional :: perturbations(:, :, :)
         !! One-electron operators in the AO basis. Defaults to the three dipole
         !! components, which is what `alpha` means; supply the six quadrupole
         !! components instead and `alpha` becomes the quadrupole-quadrupole
         !! response, which is what the `LMOQQPOL` block of a potential carries.
         !! The solver does not care -- only the right-hand sides change.
      logical, intent(in), optional :: progress
         !! Report where the two long phases are while they run. On a fragment of any
         !! size this routine is most of the wall clock and prints nothing for minutes
         !! at a time, which is indistinguishable from a hang.

      real(dp), allocatable :: dip(:, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :), gaps(:, :), h(:, :, :)
      real(dp), allocatable :: eri0(:, :, :, :), work(:, :)
      logical :: direct
      ! Every system's vectors carry a batch index: all frequencies and all
      ! perturbations are in flight together.
      logical :: reuse
      real(dp), allocatable :: aplus(:, :), aminus(:, :), product(:, :), lu(:, :)
      real(dp), allocatable :: rhs_flat(:, :), h_flat(:, :), solution(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n_ov, j, info
      real(dp), allocatable :: operators(:, :, :)
      real(dp) :: nu, rz, rz_new, pap, step, target_norm, use_tol
      integer :: n_ao, n_mo, n_vir, k, l, i, a, ifreq, iter, limit, n_pert
      character(len=32) :: text
      logical :: talk
      type(timer_type) :: clock

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
      ! **Built once, then solved densely.** Materialising `(A+B)` and `(A-B)` costs
      ! `n_ov` Fock builds. After that every frequency and every perturbation is
      ! dense linear algebra on matrices already in hand: no further integrals, no
      ! iterations, and nothing that can fail to converge.
      !
      ! The product `(A-B)(A+B)` is the same at every frequency -- only the `nu^2` on
      ! the diagonal changes -- so it is formed once too. Each frequency is then one
      ! LU factorisation and one solve per perturbation.
      !
      ! This is what makes the frequency and perturbation counts free. Matrix-free,
      ! twelve frequencies by nine operators by however many iterations each needed
      ! was thousands of builds; here it is `n_ov` of them and some BLAS.
      n_ov = n_vir*n_occ
      reuse = .false.
      if (present(hessian)) reuse = hessian%ready
      if (reuse) then
         if (size(hessian%product, 1) /= n_ov) then
            call error%set(ERROR_VALIDATION, "dynamic response: the supplied Hessian "// &
                           "is the wrong size for this reference")
            return
         end if
      else
         call build_hessian(mol, direct, eri0, bounds, zero_h, c_occ, c_vir, gaps, &
                            aplus, aminus, HESSIAN_CHUNK, error, progress=talk)
         if (error%has_error()) return
      end if

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

      if (talk) then
         write (*, "(A,I0,A,I0,A,I0,A)") "        solving ", size(frequencies), &
            " frequencies x ", n_pert, " perturbations over ", n_ov, " pairs"
         call clock%start()
      end if
      do ifreq = 1, size(frequencies)
         if (reuse) then
            lu = hessian%product
         else
            lu = product
         end if
         do j = 1, n_ov
            lu(j, j) = lu(j, j) + frequencies(ifreq)**2
         end do
         solution = rhs_flat
         call pic_getrf(lu, ipiv, info)
         if (info /= 0) then
            call error%set(ERROR_GENERIC, "dynamic response: the operator is singular")
            return
         end if
         call pic_getrs(lu, ipiv, solution, info=info)
         if (info /= 0) then
            call error%set(ERROR_GENERIC, "dynamic response: the solve failed")
            return
         end if

         do l = 1, n_pert
            do k = 1, n_pert
               alpha(k, l, ifreq) = -2.0_dp*sum(h_flat(:, k)*solution(:, l))
            end do
            if (present(response)) then
               response(:, :, l, ifreq) = reshape(solution(:, l), [n_vir, n_occ])
            end if
         end do
         if (talk) call tick(clock, ifreq, size(frequencies), "frequency")
      end do

      ! Hand the build back if the caller wants it, rather than freeing it.
      if (present(hessian) .and. .not. reuse) then
         call move_alloc(aminus, hessian%aminus)
         call move_alloc(product, hessian%product)
         hessian%ready = .true.
         if (allocated(aplus)) deallocate (aplus)
      else if (.not. reuse) then
         deallocate (aplus, aminus, product)
      end if
      deallocate (lu, ipiv, rhs_flat, h_flat, solution)

      deallocate (bounds, zero_h, c_occ, c_vir, gaps, h, work, eri0, operators)
      if (allocated(dip)) deallocate (dip)
   end subroutine dynamic_polarizability

   subroutine build_hessian(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                            aplus, aminus, chunk, error, progress)
      !! `(A+B)` and `(A-B)` as explicit matrices over the occupied-virtual space
      !!
      !! **Why materialise an operator a matrix-free solver was built to avoid.** The
      !! response is wanted at twelve frequencies for up to nine perturbations, and
      !! every one of those solves applies the same two operators. Matrix-free, the
      !! cost is the number of iterations times the number of systems times a Fock
      !! build. Built once, it is `n_ov` builds full stop -- independent of how many
      !! frequencies and perturbations follow, because each of those becomes dense
      !! linear algebra on a matrix that is already in hand.
      !!
      !! Column `j` is the operator applied to the unit vector `e_j`, so the whole
      !! thing is `n_ov` applications, done in batches so one pass over the integrals
      !! serves `chunk` columns.
      !!
      !! `n_ov` is `n_vir * n_occ` -- 550 for water in cc-pVQZ, where the two matrices
      !! come to 4.8 MB. It is the square that limits this, not the basis: a fragment
      !! with a thousand occupied-virtual pairs needs 16 MB and one with ten thousand
      !! needs 1.6 GB.
      type(libcint_molecule_t), intent(in) :: mol
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

      talk = .false.
      if (present(progress)) talk = progress
      n_vir = size(gaps, 1)
      n_occ = size(gaps, 2)
      n_ov = n_vir*n_occ
      allocate (aplus(n_ov, n_ov), aminus(n_ov, n_ov))
      allocate (u(n_vir, n_occ, chunk), au(n_vir, n_occ, chunk), idx(chunk))
      if (talk) then
         write (*, "(A,I0,A,I0,A,F0.1,A)") "        Hessian: ", n_ov, &
            " columns in chunks of ", chunk, ", ", &
            2.0_dp*real(n_ov, dp)**2*8.0_dp/1.0e6_dp, " MB"
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
      !! `build_hessian` gets each column by applying the operator to a unit vector,
      !! which costs an integral pass per chunk of columns -- `n_ov` Fock builds for
      !! the whole matrix, and that is the whole cost of a potential. The operators
      !! are, written out,
      !!
      !!     (A+B)_{ai,bj} = d Deps + 4(ai|bj) - (ab|ij) - (aj|ib)
      !!     (A-B)_{ai,bj} = d Deps          - (ab|ij) + (aj|ib)
      !!
      !! so with the integrals fitted, `B^P_pq`, every term is a matrix product over
      !! the auxiliary index and the column loop disappears:
      !!
      !!   * `(ai|bj)` is `B_ov B_ov^T` -- one gemm, and because
      !!     `build_df_mo_block` runs its compound index left-fastest that product
      !!     already sits in the `n_ov` layout the solver uses, with no repacking.
      !!   * `(aj|ib)` needs no product of its own. It is the *same* matrix read at
      !!     permuted indices: `(aj|ib)` is element `[(a,j),(b,i)]` of it. What looks
      !!     like a third integral class is a transpose of the first.
      !!   * `(ab|ij)` is `B_vv B_oo^T`, a second gemm of the same flop count, whose
      !!     result is indexed by `(ab)` and `(ij)` and so is scattered rather than
      !!     read straight off.
      !!
      !! Two gemms of `n_ov^2 naux` replace `n_ov` Fock builds. The result is not the
      !! same matrix: fitting is an approximation, and how good it is is what
      !! `validation/check_df_hessian` measures against the exact build rather than
      !! against GAMESS, so the fitting error is not confused with anything else.
      !!
      !! The memory is unchanged and is what limits this -- `n_ov^2` for each
      !! operator, plus `n_ov^2` again for each of the two integral matrices while
      !! they are assembled.
      type(libcint_molecule_t), intent(in) :: orb, aux
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), allocatable, intent(out) :: aplus(:, :), aminus(:, :)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: progress

      real(dp), allocatable :: bov(:, :), bvv(:, :), boo(:, :)
      real(dp), allocatable :: coul(:, :), exch(:, :), gap_flat(:)
      real(dp) :: term_abij, term_ajib
      integer :: n_vir, n_occ, n_ov, naux, row, col, a, i, b, j
      logical :: talk
      type(timer_type) :: clock

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
         write (*, "(A,I0,A,I0,A,F0.1,A)") "        fitted Hessian: ", n_ov, &
            " pairs, ", naux, " auxiliary functions, ", &
            4.0_dp*real(n_ov, dp)**2*8.0_dp/1.0e6_dp, " MB"
         call tick(clock, 1, 3, "step")
      end if

      allocate (coul(n_ov, n_ov))
      call pic_gemm(bov, bov, coul, transb="T")
      allocate (exch(n_vir*n_vir, n_occ*n_occ))
      call pic_gemm(bvv, boo, exch, transb="T")
      deallocate (bov, bvv, boo)
      if (talk) call tick(clock, 2, 3, "step")

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
            ! virtual-occupied, so both are indices into the same fitted matrix.
            term_ajib = coul(a + (j - 1)*n_vir, b + (i - 1)*n_vir)
            aplus(row, col) = 4.0_dp*coul(row, col) - term_abij - term_ajib
            aminus(row, col) = term_ajib - term_abij
         end do
      end do
      !$omp end parallel do
      deallocate (coul, exch)

      gap_flat = reshape(gaps, [n_ov])
      do col = 1, n_ov
         aplus(col, col) = aplus(col, col) + gap_flat(col)
         aminus(col, col) = aminus(col, col) + gap_flat(col)
      end do
      deallocate (gap_flat)
      if (talk) call tick(clock, 3, 3, "step")
   end subroutine build_hessian_df

   subroutine tick(clock, done, total, unit_name)
      !! One progress line: how far in, how long it has taken, what is left
      !!
      !! The estimate assumes the remaining work costs what the finished work did,
      !! which holds for the Hessian columns -- every chunk is the same integral pass
      !! -- and for the frequency solves, where each is the same factorization. It is
      !! printed as an estimate rather than a promise because screening makes some
      !! chunks cheaper than others.
      type(timer_type), intent(inout) :: clock
      integer, intent(in) :: done, total
      character(len=*), intent(in) :: unit_name

      real(dp) :: so_far, left

      ! Read the clock without stopping it: `start` resets, so a stop/start pair here
      ! would leave `so_far` measuring only the interval since the previous line and
      ! the estimate would be built from one chunk rather than all of them.
      so_far = clock%get_elapsed_time()
      left = so_far*real(total - done, dp)/real(max(done, 1), dp)
      write (*, "(A,I6,A,I0,1X,A,A,F8.1,A,F8.1,A)") "        ", done, " of ", total, &
         trim(unit_name), "s   ", so_far, " s in, ~", left, " s left"
      flush (output_unit)
   end subroutine tick

   subroutine response_batch(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                             u, idx, nact, minus, au, error)
      !! `(A+B)u` or `(A-B)u` for many vectors in **one** integral pass
      !!
      !! The reason this exists. Every right-hand side of the dynamic response shares
      !! one operator and differs only in the frequency shift, so the integrals a
      !! direct build recomputes are the same integrals for all of them. Contracting
      !! a batch of densities in a single pass over the quartets amortizes that cost
      !! across the batch -- with twelve frequencies and nine perturbations the pass
      !! count drops by two orders of magnitude, without storing anything.
      !!
      !! `idx(1:nact)` selects which of the vectors are still wanted, so a converged
      !! system stops costing anything rather than riding along.
      !!
      !! `minus` picks the antisymmetric combination and the build that does not fold
      !! permutations, which is what `A - B` needs; otherwise the symmetric one.
      type(libcint_molecule_t), intent(in) :: mol
      logical, intent(in) :: direct
      real(dp), intent(in) :: eri(:, :, :, :)
      real(dp), intent(in) :: bounds(:, :), zero_h(:, :)
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), intent(in) :: u(:, :, :)
      integer, intent(in) :: idx(:), nact
      logical, intent(in) :: minus
      real(dp), intent(inout) :: au(:, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dens(:, :, :), g(:, :, :), half(:, :), work(:, :)
      type(direct_stats_t) :: stats
      integer :: n_ao, n_occ, m, j

      if (nact <= 0) return
      n_ao = size(c_occ, 1)
      n_occ = size(c_occ, 2)
      allocate (dens(n_ao, n_ao, nact), half(n_ao, n_occ), work(n_ao, n_occ))

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

      do m = 1, nact
         j = idx(m)
         call pic_gemm(g(:, :, m), c_occ, work)
         call pic_gemm(c_vir, work, au(:, :, j), transa="T")
         au(:, :, j) = gaps*u(:, :, j) + 2.0_dp*au(:, :, j)
      end do

      deallocate (dens, g, half, work)
   end subroutine response_batch

   subroutine apply_dynamic_batch(mol, direct, bounds, zero_h, eri, c_occ, c_vir, &
                                  gaps, nu, u, idx, nact, au, limit, tol, error)
      !! `[(A+B) + nu^2 (A-B)^-1] u` for many vectors, each with its own frequency
      !!
      !! The inner solves run in lockstep for the same reason the outer one does:
      !! their operator is shared, so one pass serves every system still iterating.
      !! They converge at different rates, and a system that has converged drops out
      !! of the batch instead of being carried.
      type(libcint_molecule_t), intent(in) :: mol
      logical, intent(in) :: direct
      real(dp), intent(in) :: bounds(:, :), zero_h(:, :)
      real(dp), intent(in) :: eri(:, :, :, :)
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), intent(in) :: nu(:)
      real(dp), intent(in) :: u(:, :, :)
      integer, intent(in) :: idx(:), nact
      real(dp), intent(inout) :: au(:, :, :)
      integer, intent(in) :: limit
      real(dp), intent(in) :: tol
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: x(:, :, :), r(:, :, :), z(:, :, :), pv(:, :, :), ap(:, :, :)
      real(dp), allocatable :: rz(:), target_norm(:)
      integer, allocatable :: inner_idx(:)
      integer :: n_vir, n_occ, nb, m, j, iter, ninner
      real(dp) :: pap, step, rz_new

      call response_batch(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                          u, idx, nact, .false., au, error)
      if (error%has_error()) return

      nb = size(u, 3)
      if (all(nu == 0.0_dp)) return
      n_vir = size(u, 1)
      n_occ = size(u, 2)
      allocate (x(n_vir, n_occ, nb), r(n_vir, n_occ, nb), z(n_vir, n_occ, nb))
      allocate (pv(n_vir, n_occ, nb), ap(n_vir, n_occ, nb))
      allocate (rz(nb), target_norm(nb), inner_idx(nb))

      ninner = 0
      do m = 1, nact
         j = idx(m)
         if (nu(j) == 0.0_dp) cycle
         ninner = ninner + 1
         inner_idx(ninner) = j
         x(:, :, j) = 0.0_dp
         r(:, :, j) = u(:, :, j)
         z(:, :, j) = r(:, :, j)/gaps
         pv(:, :, j) = z(:, :, j)
         rz(j) = sum(r(:, :, j)*z(:, :, j))
         target_norm(j) = sqrt(sum(u(:, :, j)**2))*min(tol, 1.0e-13_dp)
      end do

      do iter = 1, limit
         if (ninner == 0) exit
         call response_batch(mol, direct, eri, bounds, zero_h, c_occ, c_vir, gaps, &
                             pv, inner_idx, ninner, .true., ap, error)
         if (error%has_error()) return
         m = 0
         do j = 1, ninner
            associate (jj => inner_idx(j))
               pap = sum(pv(:, :, jj)*ap(:, :, jj))
               if (pap <= 0.0_dp) then
                  call error%set(ERROR_VALIDATION, "A - B is not positive definite, "// &
                                 "so the reference is not a minimum")
                  return
               end if
               step = rz(jj)/pap
               x(:, :, jj) = x(:, :, jj) + step*pv(:, :, jj)
               r(:, :, jj) = r(:, :, jj) - step*ap(:, :, jj)
               if (sqrt(sum(r(:, :, jj)**2)) > target_norm(jj)) then
                  z(:, :, jj) = r(:, :, jj)/gaps
                  rz_new = sum(r(:, :, jj)*z(:, :, jj))
                  pv(:, :, jj) = z(:, :, jj) + (rz_new/rz(jj))*pv(:, :, jj)
                  rz(jj) = rz_new
                  m = m + 1
                  inner_idx(m) = jj
               end if
            end associate
         end do
         ninner = m
      end do

      do m = 1, nact
         j = idx(m)
         if (nu(j) == 0.0_dp) cycle
         au(:, :, j) = au(:, :, j) + nu(j)*nu(j)*x(:, :, j)
      end do

      deallocate (x, r, z, pv, ap, rz, target_norm, inner_idx)
   end subroutine apply_dynamic_batch

   subroutine apply_dynamic(mol, direct, bounds, zero_h, eri0, c_occ, c_vir, gaps, nu, u, au, &
                            scratch, limit, tol, error)
      !! `[(A+B) + nu^2 (A-B)^-1] u`, the inner solve included
      !!
      !! The `(A-B)` inverse is applied by its own conjugate gradients. Nested, and
      !! the inner tolerance is deliberately far tighter than the outer one: an
      !! inexact inner solve makes the outer operator non-linear, and conjugate
      !! gradients on a non-linear operator loses its conjugacy quietly rather than
      !! failing. `A - B` is well conditioned -- for water it needs a handful of
      !! iterations -- so tightening it is nearly free.
      type(libcint_molecule_t), intent(in) :: mol
      logical, intent(in) :: direct
      real(dp), intent(in) :: bounds(:, :), zero_h(:, :)
      real(dp), intent(in) :: eri0(:, :, :, :)
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), intent(in) :: nu
      real(dp), intent(in) :: u(:, :)
      real(dp), intent(out) :: au(:, :)
      real(dp), intent(inout) :: scratch(:, :)
      integer, intent(in) :: limit
      real(dp), intent(in) :: tol
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: r(:, :), z(:, :), p(:, :), ap(:, :)
      real(dp) :: rz, rz_new, pap, step, target_norm
      integer :: iter

      call response_operator(mol, direct, eri0, bounds, zero_h, c_occ, c_vir, gaps, &
                             u, au, error)
      if (error%has_error()) return
      if (nu == 0.0_dp) return

      ! (A-B) x = u, by conjugate gradients on the same preconditioner.
      allocate (r, source=u)
      allocate (z, mold=u)
      allocate (p, mold=u)
      allocate (ap, mold=u)
      scratch = 0.0_dp
      target_norm = sqrt(sum(u*u))*min(tol, 1.0e-13_dp)
      z = r/gaps
      p = z
      rz = sum(r*z)
      do iter = 1, limit
         call response_operator_minus(mol, direct, eri0, bounds, zero_h, c_occ, c_vir, gaps, p, ap, &
                                      error)
         if (error%has_error()) return
         pap = sum(p*ap)
         if (pap <= 0.0_dp) then
            call error%set(ERROR_VALIDATION, "A - B is not positive definite, so the "// &
                           "reference SCF is unstable")
            return
         end if
         step = rz/pap
         scratch = scratch + step*p
         r = r - step*ap
         if (sqrt(sum(r*r)) <= target_norm) exit
         z = r/gaps
         rz_new = sum(r*z)
         p = z + (rz_new/rz)*p
         rz = rz_new
      end do

      au = au + nu*nu*scratch
      deallocate (r, z, p, ap)
   end subroutine apply_dynamic

   subroutine distributed_dynamic_cross(mol, orbitals, orbital_energies, n_occ, &
                                        frequencies, measure, respond, tensors, &
                                        centroids, error, n_core, max_iter, tol, hessian, &
                                        progress)
      !! Mixed-multipole dynamic response, per localized orbital and frequency
      !!
      !!     alpha^(i)_{km}(i nu) = -2 sum_a h^{measure,k}_{ai} S^{respond,m}_{ai}
      !!
      !! `measure` is the operator the response is read off with, `respond` the one
      !! that drives it. All three dynamic blocks of a potential are cases of this
      !! and the solver does not distinguish them: dipole against dipole gives
      !! `DYNAMIC POLARIZABLE POINTS`, dipole against quadrupole the
      !! dipole-quadrupole block, quadrupole against quadrupole `LMOQQPOL`.
      !!
      !! **The quadrupole to pass is the traceless Buckingham form**, not the raw
      !! second moment, and for the quadrupole-quadrupole case GAMESS carries an
      !! extra factor of one third that the mixed case does not. Both are applied by
      !! the caller; see `mqc_efp_potential`, which also holds the write-time
      !! translation from the centre of mass to each centroid.
      !!
      !! **Per orbital, which operator measures is not a free choice.** One order
      !! gives `h P M^-1 h` and the other `h M^-1 P h`, and the projector onto the
      !! localized set does not commute with the response operator, so the two differ
      !! -- by 1.4e-01 against 4.5e-01 for the dipole-quadrupole block. Summed over
      !! every orbital `P` is the identity and they agree, which is what makes a
      !! summed comparison the way to test the quantity independently of the
      !! decomposition. It is also why each per-orbital tensor is asymmetric and only
      !! the total is symmetric.
      !!
      !! Validated against GAMESS: dipole-dipole to 8e-05 elementwise, and the mixed
      !! and quadrupole-quadrupole blocks through the dispersion energies GAMESS
      !! computes from an emitted potential -- `E7` to 4e-10, `E8` to 2e-07.
      !!
      !! Sign conventions do not enter: each component is quadratic in its localized
      !! orbital, so flipping an orbital's phase leaves it unchanged. That is why
      !! these compare elementwise against a reference where the Fock matrix cannot.
      type(libcint_molecule_t), intent(in) :: mol
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
                                  hessian=hessian, progress=progress)
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
                                                 tol, hessian, progress)
      !! One tensor per localized orbital at each imaginary frequency
      !!
      !! What `DYNAMIC POLARIZABLE POINTS` carries, and what dispersion is built
      !! from: the Casimir-Polder integral over these is the `C6` between two
      !! fragments. The projection is the same one `distributed_polarizability`
      !! uses -- insert `W W^T` into the occupied sum -- applied to `S` at each
      !! frequency instead of to the static `U`.
      !!
      !! The `nu = 0` member of the set must reproduce the static distributed
      !! tensors exactly, orbital by orbital, which is the check that this
      !! projection and that one agree.
      type(libcint_molecule_t), intent(in) :: mol
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
                                  progress=progress)
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
