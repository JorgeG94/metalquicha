!! The analytic Hessian, for a restricted Hartree-Fock reference
module mqc_libcint_hessian
   !! Second derivatives of the energy with respect to nuclear coordinates,
   !! without displacing anything.
   !!
   !! The finite-difference Hessian this replaces costs `6N+1` gradient
   !! evaluations and inherits each one's convergence noise amplified by `1/h`.
   !! That amplification lands hardest on the low-frequency modes, which is
   !! exactly where the rigid-rotor harmonic-oscillator partition function is
   !! most sensitive -- so the noise ends up in the thermochemistry numbers
   !! people quote. It also makes a transition-state search unreliable, since
   !! that needs one negative eigenvalue whose sign a noisy near-zero mode can
   !! flip.
   !!
   !! **Built in the pieces the standard decomposition uses**, which is also
   !! how PySCF's `hessian.rhf` is arranged, so the two can be compared stage
   !! by stage rather than only at the end: the nuclear repulsion term, the
   !! per-atom perturbation that drives the coupled-perturbed equations, the
   !! response solve itself, and the explicit second-derivative assembly.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, atom_ao_blocks
   use mqc_libcint_gradient, only: one_electron_deriv, iprinv_deriv_at, &
                                   DERIV_KIN, DERIV_NUC, DERIV_OVLP
   use mqc_libcint_hess_ints, only: eri_ip1_block, hess_1e_block, hess_2e_block, &
                                    hess_rinv_block, hess_2e_contract, h1_contract, &
                                    HESS_OVLP_II, HESS_OVLP_IJ, HESS_KIN_II, HESS_KIN_IJ, &
                                    HESS_NUC_II, HESS_NUC_IJ, HESS_RINV_II, HESS_RINV_IJ, &
                                    HESS_ERI_II, HESS_ERI_IJ, HESS_ERI_IK
   use mqc_libcint_xc, only: xc_context_t, xc_kernel_apply, xc_add_potential
   use mqc_libcint_xc_hessian, only: xc_potential_deriv
   use mqc_libcint_direct, only: build_fock_direct, build_fock_direct_many, &
                                 schwarz_bounds, direct_stats_t
   use mqc_libcint_response, only: response_operator_t, solve_response
   use pic_blas_interfaces, only: pic_gemm
   implicit none
   private

   public :: hcore_deriv_atom
   public :: make_h1_atom
   public :: overlap_deriv_atom
   public :: nuclear_response_t
   public :: solve_mo1_atom
   public :: solve_mo1_batch
   public :: nuclear_repulsion_hessian
   public :: partial_hessian
   public :: response_hessian
   public :: rhf_hessian
   public :: ks_hessian
   public :: hessian_to_matrix

   type, extends(response_operator_t) :: nuclear_response_t
      !! The electronic Hessian, applied to a response that has occupied rows
      !!
      !! The trial vector is `n_mo` by `n_occ` flattened, **not** virtual by
      !! occupied. That is the whole difference from the field case: the
      !! occupied-occupied block is not zero here, it is fixed by orthonormality
      !! at minus half the overlap derivative, and it contributes a density of
      !! its own that a virtual-by-occupied layout has nowhere to put.
      !!
      !! Everything the operator needs is held rather than rebuilt: the Schwarz
      !! bounds cost a pass over shell pairs and would otherwise be recomputed
      !! once per iteration.
      type(libcint_molecule_t), pointer :: mol => null()
      real(dp), allocatable :: orbitals(:, :)     !! (n_ao, n_mo)
      real(dp), allocatable :: c_occ(:, :)        !! (n_ao, n_occ)
      real(dp), allocatable :: energies(:)        !! (n_mo)
      real(dp), allocatable :: bounds(:, :)
      real(dp), allocatable :: zero_h(:, :)
      integer :: n_occ = 0
      integer :: n_mo = 0
      integer :: n_pert = 1
         !! How many perturbations travel together. The trial vector is
         !! `n_mo` by `n_occ` by this, and the whole point of it being more
         !! than one is that a Fock build costs an integral pass whether it
         !! contracts one density or a dozen.
      type(xc_context_t), pointer :: xc => null()
         !! The exchange-correlation context when the reference was Kohn-Sham.
         !! Null for Hartree-Fock, which is what makes the kernel term opt-in
         !! rather than a branch on the functional name.
      real(dp), allocatable :: reference(:, :)
         !! The converged density the kernel is evaluated at. The kernel is a
         !! second derivative of the energy, so it has a point to be taken at,
         !! and it is not the trial density being contracted.
      real(dp) :: k_scale = 1.0_dp
      real(dp) :: rs_k_lr = 0.0_dp
      real(dp) :: rs_omega = 0.0_dp
         !! A range-separated functional's second exchange term:
         !! `rs_k_lr` of the exchange built against `erf(omega r)/r`. Zero
         !! omega is the ordinary case and costs nothing.
         !! Exact exchange in the response operator. One for Hartree-Fock, zero
         !! for a pure functional, the mixing fraction for a hybrid.
   contains
      procedure :: apply => nuclear_apply
      procedure :: length => nuclear_length
   end type nuclear_response_t

contains

   subroutine hcore_deriv_atom(mol, iatom, hcore_a, error)
      !! `dH_core/dR_A`, the one-electron Hamiltonian moved by one atom
      !!
      !! Two contributions that look alike and are not. Moving atom `A` moves
      !! the basis functions centred on it, which is a derivative of the
      !! integral's bra and ket and only touches that atom's block. It also
      !! moves the **nucleus**, and every electron feels that wherever its
      !! orbital sits -- the Hellmann-Feynman term, which involves no
      !! basis-function derivative at all and is what differentiating the
      !! origin of `1/|r-R|` gives.
      !!
      !! Symmetrised at the end because the library puts the nabla on the bra,
      !! so what comes back is one half of a derivative that is symmetric in
      !! its two indices.
      !!
      !! The sign is the library's: its `ip` integrals carry a nabla on the bra,
      !! and the derivative with respect to the atom the bra sits on is minus
      !! that.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: iatom
      real(dp), allocatable, intent(out) :: hcore_a(:, :, :)   !! (n_ao, n_ao, 3)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: kin(:, :, :), nuc(:, :, :), vrinv(:, :, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: nao, comp, p0, p1

      if (error%has_error()) return
      if (iatom < 1 .or. iatom > mol%natm) then
         call error%set(ERROR_VALIDATION, "the core Hamiltonian derivative was asked "// &
                        "for with respect to an atom outside the molecule.")
         return
      end if

      nao = mol%nao
      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)

      call one_electron_deriv(mol, kin, DERIV_KIN)
      call one_electron_deriv(mol, nuc, DERIV_NUC)
      allocate (vrinv(nao, nao, 3), hcore_a(nao, nao, 3))

      call iprinv_deriv_at(mol, iatom, vrinv)
      vrinv = -mol%charges(iatom)*vrinv

      p0 = offsets(iatom) + 1
      p1 = offsets(iatom) + counts(iatom)
      if (counts(iatom) > 0) then
         vrinv(p0:p1, :, :) = vrinv(p0:p1, :, :) - (kin(p0:p1, :, :) + nuc(p0:p1, :, :))
      end if

      do comp = 1, 3
         hcore_a(:, :, comp) = vrinv(:, :, comp) + transpose(vrinv(:, :, comp))
      end do

      deallocate (kin, nuc, vrinv, offsets, counts)
   end subroutine hcore_deriv_atom

   subroutine make_h1_atom(mol, density, eri_ip1, iatom, h1, error)
      !! The perturbation that drives the coupled-perturbed equations
      !!
      !!     h1_A = dH_core/dR_A + dV_HF/dR_A
      !!
      !! with the mean field differentiated with respect to the same atom. It
      !! needs **first** derivatives only -- the second derivatives belong to
      !! the explicit part of the Hessian, not here, which is easy to assume the
      !! other way round.
      !!
      !! **A quartet contributes through every index that sits on this atom.**
      !! `int2e_ip1` puts the nabla on the first index alone, so the derivative
      !! with respect to atom `A` is assembled by permuting each of the four
      !! positions into first place in turn and keeping the ones whose orbital
      !! belongs to `A`. The permutations used are the ordinary symmetries of
      !! an undifferentiated integral -- `(ij|kl) = (ji|kl) = (kl|ij)` -- which
      !! hold because the nabla is applied *after* the permutation, not before.
      !!
      !! The sign is the library's, as everywhere else here: `ip` carries a
      !! nabla on the bra and the derivative with respect to that centre is
      !! minus it.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)          !! Total AO density
      real(dp), intent(in) :: eri_ip1(:, :, :, :, :)  !! From `eri_ip1_block`
      integer, intent(in) :: iatom
      real(dp), allocatable, intent(out) :: h1(:, :, :)   !! (n_ao, n_ao, 3)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: hcore_a(:, :, :), vhf(:, :, :)
      integer, allocatable :: offsets(:), counts(:), owner(:)
      integer :: nao, comp, mu, nu, la, si, a, p0, p1
      real(dp) :: d

      if (error%has_error()) return

      nao = mol%nao
      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)

      ! Which atom each basis function belongs to, so the inner loop can ask
      ! rather than search.
      allocate (owner(nao))
      owner = 0
      do a = 1, mol%natm
         p0 = offsets(a) + 1
         p1 = offsets(a) + counts(a)
         if (counts(a) > 0) owner(p0:p1) = a
      end do

      allocate (vhf(nao, nao, 3))
      vhf = 0.0_dp

      ! J - K/2, differentiated. Written as a plain quadruple loop over the
      ! whole basis: this is the readable form, and the shell-driven one that
      ! replaces it will be checked against it.
      do comp = 1, 3
         do si = 1, nao
            do la = 1, nao
               do nu = 1, nao
                  do mu = 1, nao
                     d = density(la, si)
                     if (abs(d) < 1.0e-14_dp) cycle
                     ! Coulomb: `(mu nu | la si)`, with the nabla on whichever
                     ! index sits on this atom, each permuted into first place.
                     if (owner(mu) == iatom) then
                        vhf(mu, nu, comp) = vhf(mu, nu, comp) &
                                            - d*eri_ip1(mu, nu, la, si, comp)
                     end if
                     if (owner(nu) == iatom) then
                        vhf(mu, nu, comp) = vhf(mu, nu, comp) &
                                            - d*eri_ip1(nu, mu, la, si, comp)
                     end if
                     if (owner(la) == iatom) then
                        vhf(mu, nu, comp) = vhf(mu, nu, comp) &
                                            - d*eri_ip1(la, si, mu, nu, comp)
                     end if
                     if (owner(si) == iatom) then
                        vhf(mu, nu, comp) = vhf(mu, nu, comp) &
                                            - d*eri_ip1(si, la, mu, nu, comp)
                     end if
                     ! Exchange: `(mu la | si nu)`, half weight for a closed
                     ! shell, where the total density already carries its two.
                     if (owner(mu) == iatom) then
                        vhf(mu, nu, comp) = vhf(mu, nu, comp) &
                                            + 0.5_dp*d*eri_ip1(mu, la, si, nu, comp)
                     end if
                     if (owner(la) == iatom) then
                        vhf(mu, nu, comp) = vhf(mu, nu, comp) &
                                            + 0.5_dp*d*eri_ip1(la, mu, si, nu, comp)
                     end if
                     if (owner(si) == iatom) then
                        vhf(mu, nu, comp) = vhf(mu, nu, comp) &
                                            + 0.5_dp*d*eri_ip1(si, nu, mu, la, comp)
                     end if
                     if (owner(nu) == iatom) then
                        vhf(mu, nu, comp) = vhf(mu, nu, comp) &
                                            + 0.5_dp*d*eri_ip1(nu, si, mu, la, comp)
                     end if
                  end do
               end do
            end do
         end do
      end do

      call hcore_deriv_atom(mol, iatom, hcore_a, error)
      if (error%has_error()) return

      allocate (h1(nao, nao, 3))
      h1 = vhf + hcore_a

      deallocate (vhf, hcore_a, offsets, counts, owner)
   end subroutine make_h1_atom

   subroutine overlap_deriv_atom(mol, iatom, s1, error)
      !! `dS/dR_A`, the overlap matrix moved by one atom
      !!
      !! **The reason a nuclear Hessian needs a different coupled-perturbed
      !! solve from every other perturbation in this code.** An electric field
      !! does not move the basis functions, so the overlap is unchanged and its
      !! derivative is zero; displacing a nucleus drags its functions with it,
      !! and the orbitals must stay orthonormal while that happens. That
      !! constraint is what puts a non-zero occupied-occupied block into the
      !! orbital response, and this matrix is where it comes from.
      !!
      !! Only functions centred on `A` move, so the derivative is the
      !! differentiated overlap restricted to that atom's rows, plus its
      !! transpose in the same columns -- one term for the bra moving and one
      !! for the ket. The sign is the library's, as everywhere else here.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: iatom
      real(dp), allocatable, intent(out) :: s1(:, :, :)   !! (n_ao, n_ao, 3)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: ds(:, :, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: comp, p0, p1

      if (error%has_error()) return
      if (iatom < 1 .or. iatom > mol%natm) then
         call error%set(ERROR_VALIDATION, "the overlap derivative was asked for with "// &
                        "respect to an atom outside the molecule.")
         return
      end if

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)

      call one_electron_deriv(mol, ds, DERIV_OVLP)
      ds = -ds

      allocate (s1(mol%nao, mol%nao, 3))
      s1 = 0.0_dp
      p0 = offsets(iatom) + 1
      p1 = offsets(iatom) + counts(iatom)
      if (counts(iatom) > 0) then
         do comp = 1, 3
            s1(p0:p1, :, comp) = s1(p0:p1, :, comp) + ds(p0:p1, :, comp)
            s1(:, p0:p1, comp) = s1(:, p0:p1, comp) + transpose(ds(p0:p1, :, comp))
         end do
      end if

      deallocate (ds, offsets, counts)
   end subroutine overlap_deriv_atom

   pure function nuclear_length(this) result(n)
      !! `n_mo` by `n_occ`, occupied rows included
      class(nuclear_response_t), intent(in) :: this
      integer :: n

      n = this%n_mo*this%n_occ*this%n_pert
   end function nuclear_length

   subroutine nuclear_apply(this, vector, image, error)
      !! The two-electron response a trial vector induces, scaled and flattened
      !!
      !! Four steps and each is somewhere the conventions can go wrong. The
      !! trial vector becomes an atomic-orbital density; that density is
      !! **symmetrised**, which is what makes the fast Fock build correct here
      !! -- the antisymmetrised variant belongs to the frequency-dependent
      !! `A - B` block and needs a different build entirely. The resulting `G`
      !! goes back to the molecular basis, and finally the virtual rows are
      !! divided by their orbital energy gaps while the occupied rows are set to
      !! zero.
      !!
      !! **The occupied rows are zeroed rather than solved for.** They are
      !! already known -- orthonormality fixes them -- so the iteration must not
      !! move them, and the solver puts them back from the right-hand side on
      !! every pass.
      !!
      !! Allocates its scratch on every call, and this runs once per iteration.
      !! At the sizes this is correct for that is invisible against the Fock
      !! build it wraps; on anything large the arrays want to live in the type
      !! alongside the bounds. Left as it is because moving them there before
      !! the assembly exists would be guessing at what the assembly needs.
      class(nuclear_response_t), intent(inout) :: this
      real(dp), intent(in) :: vector(:)
      real(dp), intent(out) :: image(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: mo1(:, :, :), half(:, :), dens(:, :, :), g(:, :, :)
      real(dp), allocatable :: vxc(:, :)
      real(dp), allocatable :: g_lr(:, :, :)
      real(dp), allocatable :: work(:, :), gmo(:, :)
      type(direct_stats_t) :: stats
      integer :: n_ao, i, a, p

      if (error%has_error()) return

      n_ao = size(this%orbitals, 1)
      allocate (mo1(this%n_mo, this%n_occ, this%n_pert))
      mo1 = reshape(vector, [this%n_mo, this%n_occ, this%n_pert])

      ! The densities these responses imply, symmetrised, one per perturbation.
      allocate (half(n_ao, this%n_occ), dens(n_ao, n_ao, this%n_pert))
      do p = 1, this%n_pert
         ! **Twice, for double occupancy.** The trial vector describes how one
         ! spatial orbital moves; the density it perturbs holds two electrons in
         ! that orbital. Leaving the factor out halves the response, which does not
         ! stop the iteration converging -- it converges to the wrong answer, and
         ! the first-order density comes out roughly a third too large rather than
         ! obviously doubled, because the error feeds back through the coupling.
         call pic_gemm(this%orbitals, mo1(:, :, p), half, alpha=2.0_dp, beta=0.0_dp)
         call pic_gemm(half, this%c_occ, dens(:, :, p), transb="T")
         dens(:, :, p) = dens(:, :, p) + transpose(dens(:, :, p))
      end do

      ! **One integral pass for the whole batch.** This is the entire reason the
      ! perturbations travel together: the quartets do not depend on which
      ! density is being contracted, and evaluating them is what the build
      ! costs. `build_fock_direct_many` screens on the Schwarz bound alone, so
      ! a batch is bit-for-bit what the same densities give one at a time --
      ! which matters here, because a response density is not a density and
      ! screening on its magnitude would be wrong.
      call build_fock_direct_many(this%mol, this%zero_h, dens, this%bounds, g, stats, error, &
                                  k_scale=this%k_scale)
      if (error%has_error()) return

      ! The long-range half of a range-separated functional's exchange. No
      ! Coulomb: the full-range pass above already supplied it.
      if (this%rs_omega > 0.0_dp) then
         if (.not. allocated(g_lr)) allocate (g_lr(size(g, 1), size(g, 2), size(g, 3)))
         call build_fock_direct_many(this%mol, this%zero_h, dens, this%bounds, g_lr, &
                                     stats, error, k_scale=this%rs_k_lr, j_scale=0.0_dp, &
                                     omega=this%rs_omega)
         if (error%has_error()) return
         g = g + g_lr
      end if

      ! The exchange-correlation kernel, for a Kohn-Sham reference. `J - K/2`
      ! above is the whole two-electron response of a Hartree-Fock operator; a
      ! functional's response also contains the second derivative of E_xc with
      ! respect to the density, which is what `xc_kernel_apply` contracts
      ! against the trial density. Leaving it out does not fail -- the solve
      ! still converges -- it converges to the wrong orbital response, and the
      ! Hessian is then wrong by a few percent with nothing to object.
      if (associated(this%xc)) then
         if (.not. allocated(vxc)) allocate (vxc(size(dens, 1), size(dens, 2)))
         do p = 1, size(dens, 3)
            ! Zeroed per perturbation: `xc_kernel_apply` accumulates into its
            ! output rather than assigning, so a shared buffer carries the
            ! previous perturbation's kernel into this one.
            vxc = 0.0_dp
            call xc_kernel_apply(this%xc, this%mol, this%reference, dens(:, :, p), vxc, error)
            if (error%has_error()) return
            g(:, :, p) = g(:, :, p) + vxc
         end do
      end if
      if (error%has_error()) return

      ! Back to the molecular basis, occupied columns only.
      allocate (work(n_ao, this%n_occ), gmo(this%n_mo, this%n_occ))
      do p = 1, this%n_pert
         call pic_gemm(g(:, :, p), this%c_occ, work)
         call pic_gemm(this%orbitals, work, gmo, transa="T")
         do i = 1, this%n_occ
            do a = 1, this%n_mo
               if (a <= this%n_occ) then
                  gmo(a, i) = 0.0_dp
               else
                  gmo(a, i) = gmo(a, i)/(this%energies(i) - this%energies(a))
               end if
            end do
         end do
         mo1(:, :, p) = gmo
      end do

      image = reshape(mo1, [this%n_mo*this%n_occ*this%n_pert])
      deallocate (mo1, half, dens, g, work, gmo)
   end subroutine nuclear_apply

   subroutine solve_mo1_atom(mol, orbitals, energies, n_occ, h1, s1, mo1, error, &
                             max_iter, tol)
      !! The first-order orbitals for one atom's displacement
      !!
      !! Assembles the right-hand side and hands it to the shared solver.
      !!
      !! **The occupied-occupied block is not solved for.** Orthonormality fixes
      !! it at minus half the overlap derivative, so it goes into the
      !! right-hand side and the operator zeroes it on every pass, which leaves
      !! it exactly where it was put. The virtual-occupied block is the only
      !! unknown, and it is coupled to the fixed block through the density the
      !! whole vector induces -- which is the reason the occupied rows have to
      !! be carried at all rather than dropped.
      !!
      !! The right-hand side is `h1 - s1 e_i` divided by the orbital energy
      !! gap, with the sign that makes the response enter the fixed point with
      !! the opposite one. That asymmetry is not a typo: the base is the
      !! uncoupled solution of an equation whose response term sits on the other
      !! side of the equals sign.
      type(libcint_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: h1(:, :, :)    !! (n_ao, n_ao, 3), from `make_h1_atom`
      real(dp), intent(in) :: s1(:, :, :)    !! (n_ao, n_ao, 3), from `overlap_deriv_atom`
      real(dp), allocatable, intent(out) :: mo1(:, :, :)
         !! (n_mo, n_occ, 3), the first-order orbitals in the molecular basis
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol

      ! Rebuilt per atom, including the Schwarz bounds and a copy of the
      ! orbitals. One pass over shell pairs per atom where one per molecule
      ! would do, which is worth hoisting when the assembly above this decides
      ! how it wants to loop -- and not before.
      type(nuclear_response_t) :: operator
      real(dp), allocatable :: rhs(:), answer(:), h1mo(:, :), s1mo(:, :)
      real(dp), allocatable :: work(:, :), base(:, :)
      integer :: n_ao, n_mo, comp, a, i

      if (error%has_error()) return

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)

      operator%mol => mol
      operator%orbitals = orbitals
      operator%c_occ = orbitals(:, 1:n_occ)
      operator%energies = energies
      operator%n_occ = n_occ
      operator%n_mo = n_mo
      allocate (operator%zero_h(n_ao, n_ao))
      operator%zero_h = 0.0_dp
      call schwarz_bounds(mol, operator%bounds, error)
      if (error%has_error()) return

      allocate (mo1(n_mo, n_occ, 3))
      allocate (h1mo(n_mo, n_occ), s1mo(n_mo, n_occ), base(n_mo, n_occ))
      allocate (work(n_ao, n_occ), rhs(n_mo*n_occ))

      do comp = 1, 3
         call pic_gemm(h1(:, :, comp), operator%c_occ, work)
         call pic_gemm(orbitals, work, h1mo, transa="T")
         call pic_gemm(s1(:, :, comp), operator%c_occ, work)
         call pic_gemm(orbitals, work, s1mo, transa="T")

         do i = 1, n_occ
            do a = 1, n_mo
               if (a <= n_occ) then
                  ! Fixed by orthonormality, not solved for.
                  base(a, i) = -0.5_dp*s1mo(a, i)
               else
                  base(a, i) = -(h1mo(a, i) - s1mo(a, i)*energies(i)) &
                               /(energies(a) - energies(i))
               end if
            end do
         end do

         rhs = reshape(base, [n_mo*n_occ])
         call solve_response(operator, rhs, rhs, answer, error, &
                             max_iter=max_iter, tol=tol)
         if (error%has_error()) return
         mo1(:, :, comp) = reshape(answer, [n_mo, n_occ])
         deallocate (answer)
      end do

      deallocate (h1mo, s1mo, base, work, rhs)
   end subroutine solve_mo1_atom

   subroutine nuclear_repulsion_hessian(atomic_numbers, coordinates, hess, error)
      !! Second derivatives of the nuclear repulsion
      !!
      !!     E_nn = sum_{A<B} Z_A Z_B / R_AB
      !!
      !! The only part of the Hessian with no electrons in it, and the only part
      !! that can be written down. For `A /= B` the block is
      !!
      !!     Z_A Z_B [ delta_ij / R^3 - 3 r_i r_j / R^5 ]
      !!
      !! and the diagonal block is minus the sum of the others in its row --
      !! which is not an economy but the statement that translating the whole
      !! molecule costs nothing, imposed rather than hoped for.
      !!
      !! `(3, 3, n_atoms, n_atoms)` with the Cartesian pair innermost, which is
      !! the layout the electronic part will want to add into.
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)      !! (3, n_atoms), Bohr
      real(dp), allocatable, intent(out) :: hess(:, :, :, :)
      type(error_t), intent(inout) :: error

      real(dp) :: r(3), block(3, 3)
      real(dp) :: dist, zz
      integer :: natm, a, b, i, j

      if (error%has_error()) return

      natm = size(atomic_numbers)
      if (size(coordinates, 2) /= natm) then
         call error%set(ERROR_VALIDATION, "the coordinates and the atomic numbers "// &
                        "describe different numbers of atoms.")
         return
      end if

      allocate (hess(3, 3, natm, natm))
      hess = 0.0_dp

      do a = 1, natm
         do b = 1, natm
            if (a == b) cycle
            r = coordinates(:, a) - coordinates(:, b)
            dist = norm2(r)
            if (dist < 1.0e-10_dp) then
               call error%set(ERROR_VALIDATION, "two atoms are on top of each other, "// &
                              "so the nuclear repulsion has no second derivative.")
               deallocate (hess)
               return
            end if
            zz = real(atomic_numbers(a), dp)*real(atomic_numbers(b), dp)
            do j = 1, 3
               do i = 1, 3
                  block(i, j) = -3.0_dp*zz*r(i)*r(j)/dist**5
                  if (i == j) block(i, j) = block(i, j) + zz/dist**3
               end do
            end do
            hess(:, :, a, b) = block
            ! The diagonal takes the negative of every off-diagonal partner, so
            ! each row sums to zero and a rigid translation is free.
            hess(:, :, a, a) = hess(:, :, a, a) - block
         end do
      end do
   end subroutine nuclear_repulsion_hessian

   subroutine partial_hessian(mol, density, weighted, hess, error, k_scale, &
                              rs_k_lr, rs_omega)
      !! The Hessian of the energy expression with the orbitals held fixed
      !!
      !! Everything in `d2E/dRdR` that survives when the density is not allowed
      !! to relax: the second derivatives of the integrals, contracted with the
      !! density that the unperturbed SCF converged to. The rest -- the response
      !! -- is what `solve_mo1_atom` is for, and is added on top of this.
      !!
      !! **Assembled by depositing rather than by slicing.** Each derivative
      !! belongs to whichever atom the differentiated basis function sits on, so
      !! the loop runs over all indices once and adds each term into the atom
      !! pair its own indices name. The usual arrangement instead fixes a pair
      !! of atoms and slices the integrals to their rows, which is faster and
      !! puts the ownership rules in the slice bounds, where an off-by-one is a
      !! wrong Hessian rather than a crash. Depositing keeps the rule next to
      !! the term it belongs to.
      !!
      !! **The permutational folding is real and is worth stating.** The
      !! two-electron sum has sixteen terms, one per ordered pair of the four
      !! indices. Contracted against a weight symmetric under all eight
      !! permutations of `(mu nu|lambda sigma)`, they collapse to three: four
      !! copies of the both-on-one-index term, four of the same-pair term and
      !! eight of the cross-pair term. That is where the 4, 4 and 8 below come
      !! from -- not from spin, and not from counting a term twice.
      !!
      !! The symmetric weight is
      !!
      !!     G = 1/2 D_mn D_ls - 1/8 (D_ms D_ln + D_ml D_ns)
      !!
      !! which is the Coulomb-minus-half-exchange contraction with its exchange
      !! half written in both of the ways the eightfold symmetry allows, so that
      !! no permutation of the four indices can tell it apart from itself.
      !!
      !! Cost is `nao^4` times nine, three times over, because that is what
      !! `hess_2e_block` returns. It is a correct implementation and not yet a
      !! usable one; the fix is to drive the contraction from shells and never
      !! form the array, which changes nothing above this line.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)    !! Closed shell, carrying its factor of two
      real(dp), intent(in) :: weighted(:, :)   !! `2 sum_i eps_i C_i C_i^T`
      real(dp), allocatable, intent(out) :: hess(:, :, :, :)   !! (3, 3, natm, natm)
      real(dp), intent(in), optional :: k_scale
         !! Exact-exchange fraction; one (Hartree-Fock) by default.
      real(dp), intent(in), optional :: rs_k_lr, rs_omega
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: s2(:, :, :), sab(:, :, :)
      real(dp), allocatable :: h2(:, :, :), hab(:, :, :)
      real(dp), allocatable :: tmp(:, :, :)
      real(dp), allocatable :: r2(:, :, :), rab(:, :, :)
      integer, allocatable :: owner(:), offsets(:), counts(:)
      real(dp), allocatable :: cross(:, :, :)
      real(dp) :: total(3, 3)
      real(dp) :: w, zc
      integer :: nao, natm, iao, jao, ia, ja, a, b, comp, c

      if (error%has_error()) return

      nao = mol%nao
      natm = mol%natm
      if (size(density, 1) /= nao .or. size(weighted, 1) /= nao) then
         call error%set(ERROR_VALIDATION, "the density and the basis disagree about "// &
                        "how many orbitals there are.")
         return
      end if

      allocate (offsets(natm), counts(natm))
      call atom_ao_blocks(mol, offsets, counts)
      allocate (owner(nao))
      owner = 0
      do c = 1, natm
         do iao = offsets(c) + 1, offsets(c) + counts(c)
            owner(iao) = c
         end do
      end do

      allocate (hess(3, 3, natm, natm))
      hess = 0.0_dp

      ! ---- one electron, both derivatives on basis centres -----------------
      !
      ! `ipipnuc` and `ipnucip` carry the sum over every nucleus and its charge
      ! already, so these are the terms where the nuclei sat still and only the
      ! functions moved. The terms where a nucleus moved come further down and
      ! need the operator pinned one atom at a time.
      call hess_1e_block(mol, HESS_KIN_II, h2, error)
      call hess_1e_block(mol, HESS_NUC_II, tmp, error)
      if (error%has_error()) return
      h2 = h2 + tmp
      deallocate (tmp)

      call hess_1e_block(mol, HESS_KIN_IJ, hab, error)
      call hess_1e_block(mol, HESS_NUC_IJ, tmp, error)
      if (error%has_error()) return
      hab = hab + tmp
      deallocate (tmp)

      call hess_1e_block(mol, HESS_OVLP_II, s2, error)
      call hess_1e_block(mol, HESS_OVLP_IJ, sab, error)
      if (error%has_error()) return

      ! Twice, because the bra and the ket each contribute and their two terms
      ! are the same number: relabelling the two indices maps one onto the
      ! other, and both the density and the operator are symmetric.
      do comp = 1, 9
         a = (comp - 1)/3 + 1
         b = comp - 3*(a - 1)
         do jao = 1, nao
            ja = owner(jao)
            do iao = 1, nao
               ia = owner(iao)
               w = density(iao, jao)
               hess(a, b, ia, ia) = hess(a, b, ia, ia) + 2.0_dp*w*h2(iao, jao, comp)
               hess(a, b, ia, ja) = hess(a, b, ia, ja) + 2.0_dp*w*hab(iao, jao, comp)
               w = weighted(iao, jao)
               hess(a, b, ia, ia) = hess(a, b, ia, ia) - 2.0_dp*w*s2(iao, jao, comp)
               hess(a, b, ia, ja) = hess(a, b, ia, ja) - 2.0_dp*w*sab(iao, jao, comp)
            end do
         end do
      end do
      deallocate (h2, hab, s2, sab)

      ! ---- one electron, at least one derivative on a nucleus --------------
      !
      ! `-Z_C / |r - R_C|` depends on three positions, not two: the bra centre,
      ! the ket centre and the nucleus. Differentiating the nucleus is the
      ! Hellmann-Feynman direction and produces no basis derivative at all, so
      ! it reaches atom pairs that the block above cannot see -- a nucleus with
      ! no basis functions on it would still have a Hessian row.
      !
      ! `cross(a, b, A)` collects, for the atom `c` whose nucleus is moving, the
      ! part belonging to basis functions centred on `A`. Both the bra and the
      ! ket contribute and fold onto each other, hence the factor of two; the
      ! sign is the derivative of a basis function with respect to its own
      ! centre being minus the derivative with respect to the electron.
      allocate (cross(3, 3, natm))
      do c = 1, natm
         call hess_rinv_block(mol, c, HESS_RINV_II, r2, error)
         call hess_rinv_block(mol, c, HESS_RINV_IJ, rab, error)
         if (error%has_error()) return
         zc = mol%charges(c)

         cross = 0.0_dp
         do comp = 1, 9
            a = (comp - 1)/3 + 1
            b = comp - 3*(a - 1)
            do jao = 1, nao
               do iao = 1, nao
                  ia = owner(iao)
                  cross(a, b, ia) = cross(a, b, ia) + density(iao, jao) &
                                    *(r2(iao, jao, comp) + rab(iao, jao, comp))
               end do
            end do
         end do
         deallocate (r2, rab)

         ! `w = -Z_C` is the charge and sign of the electron-nucleus
         ! attraction, which `hess_rinv_block` deliberately leaves off.
         w = -zc
         total = 0.0_dp
         do ia = 1, natm
            total = total + cross(:, :, ia)
            do b = 1, 3
               do a = 1, 3
                  hess(a, b, ia, c) = hess(a, b, ia, c) - 2.0_dp*w*cross(a, b, ia)
                  hess(a, b, c, ia) = hess(a, b, c, ia) - 2.0_dp*w*cross(b, a, ia)
               end do
            end do
         end do
         ! Both derivatives on the same nucleus. Translational invariance of
         ! `1/|r - R_C|` in its three arguments turns that into the sum of the
         ! four basis-derivative terms, which is what `total` already holds.
         hess(:, :, c, c) = hess(:, :, c, c) + 2.0_dp*w*total
      end do
      deallocate (cross)

      ! ---- two electron ----------------------------------------------------
      !
      ! Driven from shells and contracted on the spot rather than assembled and
      ! then read back: the three arrays this would otherwise form are `nao^4`
      ! times nine each. The deposit rules live there with the loop that uses
      ! them.
      call hess_2e_contract(mol, density, hess, error, k_scale=k_scale)
      if (error%has_error()) return
      if (present(rs_omega)) then
         if (rs_omega > 0.0_dp) then
            ! The attenuated exchange, with no Coulomb of its own.
            call hess_2e_contract(mol, density, hess, error, k_scale=rs_k_lr, &
                                  j_scale=0.0_dp, omega=rs_omega)
         end if
      end if
      if (error%has_error()) return

      deallocate (owner, offsets, counts)
   end subroutine partial_hessian

   subroutine response_hessian(mol, density, orbitals, energies, n_occ, hess, error, &
                               xc, reference, k_scale, rs_k_lr, rs_omega, &
                               max_iter, tol)
      !! What the Hessian gains from letting the density relax
      !!
      !! The gradient needs no density derivative -- that is what makes the
      !! Hellmann-Feynman-plus-Pulay form exact for a variational wavefunction.
      !! The Hessian does, because differentiating that gradient a second time
      !! reaches the density this time:
      !!
      !!     dD/dR_B contracted with dF/dR_A, minus dW/dR_B with dS/dR_A
      !!
      !! where the derivatives on the right are the skeleton ones -- integrals
      !! only -- and the ones on the left are total.
      !!
      !! **`W` is taken as `D F D / 2` rather than as orbital energies.** The
      !! two are the same object: `D F D / 2 = 2 C_occ eps C_occ^T`, because
      !! `C_occ^T F C_occ` is diagonal at convergence. Writing it the first way
      !! means its derivative needs only `dD/dR` and `dF/dR`, both of which are
      !! already here, instead of the derivative of the Lagrange multiplier
      !! matrix -- an occupied-occupied object that stops being diagonal the
      !! moment the molecule is perturbed, and whose relation to the overlap
      !! derivative is exactly the kind of factor of a half that this code has
      !! no independent way to check.
      !!
      !! The coupled-perturbed equations are solved once per atom, three
      !! components at a time, which is where essentially all of the time goes.
      type(libcint_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: density(:, :)     !! Closed shell, carrying its factor of two
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: hess(:, :, :, :)   !! (3, 3, natm, natm)
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: reference(:, :)
      real(dp), intent(in), optional :: k_scale
      real(dp), intent(in), optional :: rs_k_lr, rs_omega
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol

      real(dp), allocatable :: h1(:, :, :, :), s1(:, :, :, :)
      real(dp), allocatable :: d1(:, :, :, :), w1(:, :, :, :)
      real(dp), allocatable :: mo1(:, :, :, :), s1a(:, :, :), hcore_a(:, :, :)
      real(dp), allocatable :: hcore(:, :), fock(:, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: vxc_ref(:, :), fock_lr(:, :)
      real(dp), allocatable :: h1_lr(:, :, :, :)
      real(dp) :: e_xc_ref, n_xc_ref
      real(dp), allocatable :: g1all(:, :, :), f1(:, :), half(:, :), c_occ(:, :)
      real(dp), allocatable :: tmp(:, :), work(:, :)
      real(dp), allocatable :: outer(:, :), middle(:, :)
      type(direct_stats_t) :: stats
      integer :: nao, natm, ia, ja, a, b

      if (error%has_error()) return

      nao = mol%nao
      natm = mol%natm
      allocate (c_occ(nao, n_occ))
      c_occ = orbitals(:, 1:n_occ)

      call mol%core_hamiltonian(hcore)
      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) return
      allocate (fock(nao, nao), zero_h(nao, nao))
      zero_h = 0.0_dp
      ! The converged Fock, and for a Kohn-Sham reference that means *its* Fock:
      ! exact exchange at the functional's fraction and the exchange-correlation
      ! potential added on. `W = D F D / 2` is built from this, so a
      ! Hartree-Fock matrix here gives a wrong energy-weighted density and an
      ! error spread across every atom pair -- large, smooth, and symmetric
      ! enough to look like physics.
      call build_fock_direct(mol, hcore, density, bounds, fock, stats, error, &
                             density_screen=.false., k_scale=k_scale)
      if (error%has_error()) return
      if (present(rs_omega)) then
         if (rs_omega > 0.0_dp) then
            allocate (fock_lr(nao, nao))
            call build_fock_direct(mol, zero_h, density, bounds, fock_lr, stats, error, &
                                   density_screen=.false., k_scale=rs_k_lr, &
                                   j_scale=0.0_dp, omega=rs_omega)
            if (error%has_error()) return
            fock = fock + fock_lr
            deallocate (fock_lr)
         end if
      end if
      if (present(xc)) then
         allocate (vxc_ref(nao, nao))
         vxc_ref = 0.0_dp
         call xc_add_potential(xc, mol, density, vxc_ref, e_xc_ref, n_xc_ref, error)
         if (error%has_error()) return
         fock = fock + vxc_ref
         deallocate (vxc_ref)
      end if

      ! Every atom's two-electron perturbation in one pass over the shell
      ! quartets, rather than `make_h1_atom` per atom over a stored `nao^4`
      ! array. The core Hamiltonian derivative is per-atom either way and is
      ! added below.
      call h1_contract(mol, density, h1, error, k_scale=k_scale)
      if (error%has_error()) return
      if (present(rs_omega)) then
         if (rs_omega > 0.0_dp) then
            ! Into a temporary and added: `h1_contract` allocates its output and
            ! zeroes it, so a second call in place would discard the full-range
            ! pass rather than adding to it. `hess_2e_contract` next door
            ! accumulates, which is exactly the asymmetry that makes this easy
            ! to get wrong.
            call h1_contract(mol, density, h1_lr, error, k_scale=rs_k_lr, &
                             j_scale=0.0_dp, omega=rs_omega)
            if (error%has_error()) return
            h1 = h1 + h1_lr
            deallocate (h1_lr)
         end if
      end if
      if (error%has_error()) return
      ! For a Kohn-Sham reference the Fock derivative also carries the
      ! exchange-correlation potential's own, which `h1_contract` does not
      ! produce -- it knows only about integrals.
      if (present(xc)) call xc_potential_deriv(xc, mol, density, h1, error)
      if (error%has_error()) return

      allocate (s1(nao, nao, 3, natm))
      allocate (d1(nao, nao, 3, natm), w1(nao, nao, 3, natm))
      allocate (half(nao, n_occ), tmp(nao, nao), work(nao, nao))
      allocate (f1(nao, nao), outer(nao, nao), middle(nao, nao))

      do ia = 1, natm
         call hcore_deriv_atom(mol, ia, hcore_a, error)
         call overlap_deriv_atom(mol, ia, s1a, error)
         if (error%has_error()) return
         h1(:, :, :, ia) = h1(:, :, :, ia) + hcore_a
         s1(:, :, :, ia) = s1a
         deallocate (hcore_a, s1a)
      end do

      ! Every perturbation at once, in chunks, so the Fock builds inside the
      ! iteration are shared rather than repeated `3 natm` times.
      call solve_mo1_batch(mol, orbitals, energies, n_occ, h1, s1, mo1, error, &
                           xc=xc, reference=reference, k_scale=k_scale, &
                           rs_k_lr=rs_k_lr, rs_omega=rs_omega, &
                           max_iter=max_iter, tol=tol)
      if (error%has_error()) return

      do ia = 1, natm
         do a = 1, 3
            ! The first-order density, the same expression the finite-difference
            ! test in `test_mqc_hess_ints` pins directly.
            call pic_gemm(orbitals, mo1(:, :, a, ia), half, alpha=2.0_dp, beta=0.0_dp)
            call pic_gemm(half, c_occ, tmp, transb="T")
            tmp = tmp + transpose(tmp)
            d1(:, :, a, ia) = tmp
         end do
      end do

      ! Their mean fields in one batch, chunked the same way and for the same
      ! reason as the solve above: `3 natm` separate builds pay for the same
      ! integrals `3 natm` times.
      call mean_field_batch(mol, d1, bounds, zero_h, g1all, error, &
                            xc=xc, reference=reference, k_scale=k_scale, &
                            rs_k_lr=rs_k_lr, rs_omega=rs_omega)
      if (error%has_error()) return

      do ia = 1, natm
         do a = 1, 3
            tmp = d1(:, :, a, ia)
            ! `dF/dR` in full: the skeleton derivative plus the mean field of
            ! the relaxed density. Leaving the second term out is the same
            ! mistake as leaving the response out of the Hessian, made one
            ! level down.
            f1 = h1(:, :, a, ia) + g1all(:, :, 3*(ia - 1) + a)

            ! W = D F D / 2, differentiated by the product rule: three terms,
            ! all of which get the half. The first and last are transposes of
            ! each other and the middle one is already symmetric, which is what
            ! makes the sum symmetric -- but *only* if the middle one is not
            ! symmetrised a second time. Adding it before a `(X + X^T)/2` gives
            ! it twice its weight, which is a wrong `dW/dR` that is still
            ! symmetric, still translationally invariant, and wrong by tens of
            ! percent in the Hessian.
            call pic_gemm(tmp, fock, work)
            call pic_gemm(work, density, outer)
            call pic_gemm(density, f1, work)
            call pic_gemm(work, density, middle)
            w1(:, :, a, ia) = 0.5_dp*(outer + transpose(outer) + middle)
         end do
      end do

      allocate (hess(3, 3, natm, natm))
      hess = 0.0_dp
      do ja = 1, natm
         do ia = 1, natm
            do b = 1, 3
               do a = 1, 3
                  hess(a, b, ia, ja) = sum(d1(:, :, b, ja)*h1(:, :, a, ia)) &
                                       - sum(w1(:, :, b, ja)*s1(:, :, a, ia))
               end do
            end do
         end do
      end do

      deallocate (h1, s1, d1, w1, hcore, fock, bounds, zero_h)
      deallocate (g1all, f1, half, c_occ, tmp, work, outer, middle)
   end subroutine response_hessian

   subroutine rhf_hessian(mol, atomic_numbers, density, orbitals, energies, n_occ, &
                          hess, error, max_iter, tol)
      !! The analytic Hessian of a converged restricted Hartree-Fock energy
      !!
      !! Three pieces that are checked three different ways and only mean
      !! anything added together: the nuclear repulsion, which is arithmetic;
      !! the explicit integral second derivatives against the fixed density;
      !! and the response.
      !!
      !! Returned as `(3, 3, n_atoms, n_atoms)`. A vibrational analysis wants
      !! `(3*n_atoms, 3*n_atoms)`, which is a reshape the caller can do -- and
      !! which is deliberately not done here, because the four-index form is
      !! the one where a wrong atom pair is visible.
      type(libcint_molecule_t), intent(in), target :: mol
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: hess(:, :, :, :)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol

      real(dp), allocatable :: weighted(:, :), part(:, :, :, :), resp(:, :, :, :)
      integer :: i, nao

      if (error%has_error()) return

      nao = mol%nao
      allocate (weighted(nao, nao))
      weighted = 0.0_dp
      do i = 1, n_occ
         weighted = weighted + 2.0_dp*energies(i) &
                    *matmul(reshape(orbitals(:, i), [nao, 1]), &
                            reshape(orbitals(:, i), [1, nao]))
      end do

      call nuclear_repulsion_hessian(atomic_numbers, mol%coords, hess, error)
      call partial_hessian(mol, density, weighted, part, error)
      call response_hessian(mol, density, orbitals, energies, n_occ, resp, error, &
                            max_iter=max_iter, tol=tol)
      if (error%has_error()) return

      hess = hess + part + resp
      deallocate (weighted, part, resp)
   end subroutine rhf_hessian

   subroutine ks_hessian(mol, atomic_numbers, density, orbitals, energies, n_occ, &
                         xc, k_scale, hess, error, max_iter, tol)
      !! The analytic Hessian of a converged Kohn-Sham energy
      !!
      !! `rhf_hessian` with the exchange-correlation parts put back. Three
      !! things change and each is a place a Kohn-Sham Hessian is commonly got
      !! wrong:
      !!
      !! **Exact exchange is scaled.** A pure functional has none, so the
      !! exchange half of the explicit two-electron term and of the response
      !! operator both vanish; a hybrid keeps its mixing fraction. Left at one,
      !! the Hessian is a Hartree-Fock one with a functional's density in it.
      !!
      !! **The exchange-correlation second derivatives are added** -- the
      !! explicit term, `xc_hessian`.
      !!
      !! **The response operator gains the kernel.** Without it the
      !! coupled-perturbed solve still converges, to the wrong orbital
      !! response, and the result is wrong by a few percent with nothing to
      !! object to it.
      !!
      !! The grid is held fixed throughout, as `xc_hessian` sets out.
      !!
      !! Validated against `pyscf.hessian.rks` with `grid_response = False`:
      !! water at STO-3G and LDA exchange agrees to **1.8e-8** on every entry.
      !!
      !! Three bugs were found getting there, and all three produce a Hessian
      !! that is symmetric, translationally invariant and wrong -- so they are
      !! worth naming rather than leaving in the history:
      !!
      !!   * The relaxed mean field feeding `dW/dR` needed the same kernel and
      !!     exchange fraction as the response operator. Without it the two
      !!     halves of the response disagree about what the Fock operator is.
      !!   * `xc_kernel_apply` accumulates into its output rather than
      !!     assigning. A buffer shared across perturbations carried each one's
      !!     kernel into the next: 0.49 of translational-invariance violation on
      !!     its own.
      !!   * **The converged Fock that `W = D F D / 2` is built from was a
      !!     Hartree-Fock matrix** -- full exact exchange, no `V_xc`. That was
      !!     the last 0.12, spread smoothly across every atom pair.
      !!
      !! LDA, GGA, meta-GGA and hybrids of each. Range-separated hybrids and
      !! VV10 are refused; the refusal below records what is known.
      use mqc_libcint_xc_hessian, only: xc_hessian
      type(libcint_molecule_t), intent(in), target :: mol
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      type(xc_context_t), intent(inout) :: xc
      real(dp), intent(in) :: k_scale   !! Exact-exchange fraction of the functional
      real(dp), allocatable, intent(out) :: hess(:, :, :, :)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol

      real(dp), allocatable :: weighted(:, :), part(:, :, :, :), resp(:, :, :, :)
      integer :: i, nao

      if (error%has_error()) return

      ! **VV10 is refused rather than approximated.** The second derivative of
      ! the non-local term is not implemented here, and it is a real piece of
      ! work rather than a gap in the plumbing: PySCF implements all three parts
      ! of it -- `_get_enlc_deriv2` for the energy, `_get_vnlc_deriv1` for the Fock
      ! derivative and `get_vnlc_resp` for the orbital-Hessian response --
      ! behind dedicated C kernels (`VXC_vv10nlc_hessian_eval_*`), because
      ! VV10 is a double integral over the density and its second derivative is
      ! a double-grid object.
      !
      ! On water/STO-3G the term is small: our `b97m-v` Hessian sits 1.4e-3
      ! from PySCF's, which is the size of the missing contribution and also
      ! the size of meta-GGA grid noise. That it is small here is not a reason
      ! to omit it silently -- it is worth 43 mHartree in the energy, and
      ! nothing says it stays small on a dispersion-bound system, which is what
      ! VV10 is for.
      if (xc%nlc_b > 0.0_dp) then
         call error%set(ERROR_VALIDATION, "an analytic Hessian for a functional with "// &
                        "VV10 non-local correlation is not available: the second "// &
                        "derivative of the non-local term is not implemented, and "// &
                        "omitting it would give a plausible wrong answer. Use the "// &
                        "semi-numerical path.")
         return
      end if

      nao = mol%nao
      allocate (weighted(nao, nao))
      weighted = 0.0_dp
      do i = 1, n_occ
         weighted = weighted + 2.0_dp*energies(i) &
                    *matmul(reshape(orbitals(:, i), [nao, 1]), &
                            reshape(orbitals(:, i), [1, nao]))
      end do

      call nuclear_repulsion_hessian(atomic_numbers, mol%coords, hess, error)
      call partial_hessian(mol, density, weighted, part, error, k_scale=k_scale, &
                           rs_k_lr=xc%rs_k_lr, rs_omega=xc%rs_omega)
      if (error%has_error()) return
      call xc_hessian(xc, mol, density, part, error)
      if (error%has_error()) return
      call response_hessian(mol, density, orbitals, energies, n_occ, resp, error, &
                            xc=xc, reference=density, k_scale=k_scale, &
                            rs_k_lr=xc%rs_k_lr, rs_omega=xc%rs_omega, &
                            max_iter=max_iter, tol=tol)
      if (error%has_error()) return

      hess = hess + part + resp
      deallocate (weighted, part, resp)
   end subroutine ks_hessian

   subroutine hessian_to_matrix(hess, matrix)
      !! `(3, 3, natm, natm)` to the `(3N, 3N)` a vibrational analysis reads
      !!
      !! Atom slowest, Cartesian fastest -- the order
      !! `finite_diff_hessian_from_gradients` already builds, so the two
      !! Hessian paths hand the same layout to the same analysis.
      !!
      !! Worth being a routine rather than six lines at the call site, because
      !! only one of the ways to get it wrong is harmless and it is not the
      !! obvious one. Swapping **both** index pairs -- `hess(b, a, ja, ia)` --
      !! is the genuine transpose, and a Hessian is symmetric, so it changes
      !! nothing; measured, and it passes. Swapping only the atoms is not a
      !! transpose of anything, and neither is interleaving the two indices
      !! Cartesian-slowest. Both of those move every frequency and both are
      !! caught. The rule to keep is the layout, not a symmetry argument about
      !! it.
      real(dp), intent(in) :: hess(:, :, :, :)
      real(dp), allocatable, intent(out) :: matrix(:, :)

      integer :: natm, ia, ja, a, b

      natm = size(hess, 3)
      allocate (matrix(3*natm, 3*natm))
      do ja = 1, natm
         do ia = 1, natm
            do b = 1, 3
               do a = 1, 3
                  matrix(3*(ia - 1) + a, 3*(ja - 1) + b) = hess(a, b, ia, ja)
               end do
            end do
         end do
      end do
   end subroutine hessian_to_matrix

   subroutine solve_mo1_batch(mol, orbitals, energies, n_occ, h1, s1, mo1, error, &
                              xc, reference, k_scale, rs_k_lr, rs_omega, &
                              max_iter, tol)
      !! The first-order orbitals for every atom, a dozen perturbations at a time
      !!
      !! Same equations `solve_mo1_atom` solves and the same right-hand side;
      !! what changes is how many of them are in flight when the Fock build
      !! runs. A direct build costs an integral pass whichever density it
      !! contracts, so solving `3 natm` perturbations one at a time pays for
      !! the integrals `3 natm` times over.
      !!
      !! **Chunked at a dozen rather than all at once**, which is
      !! `build_fock_direct_many`'s own measured advice and worth repeating
      !! because the shape of it is not obvious: the win saturates near
      !! fourfold and then *reverses*, since the accumulator that makes the
      !! reuse possible grows with the batch while the integral it reuses does
      !! not. Past about a dozen densities it stops fitting in cache and the
      !! build becomes memory-bound. So this is not "batch everything"; it is
      !! "batch enough".
      !!
      !! The Schwarz bounds are built once here rather than once per atom,
      !! which is the wart `solve_mo1_atom` documents and could not fix while
      !! it was the only caller.
      type(libcint_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: h1(:, :, :, :)   !! (n_ao, n_ao, 3, natm)
      real(dp), intent(in) :: s1(:, :, :, :)
      real(dp), allocatable, intent(out) :: mo1(:, :, :, :)  !! (n_mo, n_occ, 3, natm)
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: reference(:, :)
      real(dp), intent(in), optional :: k_scale
      real(dp), intent(in), optional :: rs_k_lr, rs_omega
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol

      integer, parameter :: MAX_BATCH = 12

      type(nuclear_response_t) :: operator
      real(dp), allocatable :: rhs(:), answer(:), h1mo(:, :), s1mo(:, :)
      real(dp), allocatable :: work(:, :), base(:, :, :)
      integer :: n_ao, n_mo, natm, n_pert, first, last, wide, p, q, ia, comp, a, i
      integer :: n_chunks, per_chunk

      if (error%has_error()) return

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      natm = size(h1, 4)
      n_pert = 3*natm

      ! `xc` and `reference` are one argument in two halves: the kernel apply
      ! contracts the trial density against the functional's second derivative
      ! *at* the reference density, so a caller supplying the context without
      ! the density it was built at would reach `xc_kernel_apply` through
      ! `nuclear_apply` with `operator%reference` unallocated. Refuse the pair
      ! rather than index it.
      if (present(xc) .neqv. present(reference)) then
         call error%set(ERROR_VALIDATION, "the coupled-perturbed solver was given an "// &
                        "exchange-correlation context without the reference density its "// &
                        "kernel is evaluated at, or the reverse; it needs both or neither")
         return
      end if

      operator%mol => mol
      operator%orbitals = orbitals
      operator%c_occ = orbitals(:, 1:n_occ)
      operator%energies = energies
      operator%n_occ = n_occ
      operator%n_mo = n_mo
      allocate (operator%zero_h(n_ao, n_ao))
      operator%zero_h = 0.0_dp
      if (present(xc)) operator%xc => xc
      if (present(reference)) operator%reference = reference
      if (present(k_scale)) operator%k_scale = k_scale
      if (present(rs_k_lr)) operator%rs_k_lr = rs_k_lr
      if (present(rs_omega)) operator%rs_omega = rs_omega
      call schwarz_bounds(mol, operator%bounds, error)
      if (error%has_error()) return

      allocate (mo1(n_mo, n_occ, 3, natm))
      allocate (h1mo(n_mo, n_occ), s1mo(n_mo, n_occ), work(n_ao, n_occ))

      ! Balanced rather than greedy. Eighteen perturbations taken twelve at a
      ! time leaves a chunk of six, and the small chunk pays a full integral
      ! pass for half the reuse; nine and nine costs the same two passes and
      ! amortises both of them.
      n_chunks = (n_pert + MAX_BATCH - 1)/MAX_BATCH
      per_chunk = (n_pert + n_chunks - 1)/n_chunks

      first = 1
      do while (first <= n_pert)
         last = min(first + per_chunk - 1, n_pert)
         wide = last - first + 1
         allocate (base(n_mo, n_occ, wide))

         do q = 1, wide
            p = first + q - 1
            ia = (p - 1)/3 + 1
            comp = p - 3*(ia - 1)

            call pic_gemm(h1(:, :, comp, ia), operator%c_occ, work)
            call pic_gemm(orbitals, work, h1mo, transa="T")
            call pic_gemm(s1(:, :, comp, ia), operator%c_occ, work)
            call pic_gemm(orbitals, work, s1mo, transa="T")

            do i = 1, n_occ
               do a = 1, n_mo
                  if (a <= n_occ) then
                     ! Fixed by orthonormality, not solved for.
                     base(a, i, q) = -0.5_dp*s1mo(a, i)
                  else
                     base(a, i, q) = -(h1mo(a, i) - s1mo(a, i)*energies(i)) &
                                     /(energies(a) - energies(i))
                  end if
               end do
            end do
         end do

         operator%n_pert = wide
         rhs = reshape(base, [n_mo*n_occ*wide])
         call solve_response(operator, rhs, rhs, answer, error, &
                             max_iter=max_iter, tol=tol)
         if (error%has_error()) return

         base = reshape(answer, [n_mo, n_occ, wide])
         do q = 1, wide
            p = first + q - 1
            ia = (p - 1)/3 + 1
            comp = p - 3*(ia - 1)
            mo1(:, :, comp, ia) = base(:, :, q)
         end do
         deallocate (answer, base)
         first = last + 1
      end do

      deallocate (h1mo, s1mo, work, rhs)
   end subroutine solve_mo1_batch

   subroutine mean_field_batch(mol, d1, bounds, zero_h, g1, error, xc, reference, &
                               k_scale, rs_k_lr, rs_omega)
      !! `G(D')` for every first-order density, a chunk at a time
      !!
      !! Flattens `(n_ao, n_ao, 3, natm)` to `(n_ao, n_ao, 3*natm)` and hands
      !! it to the many-density build in the same chunks the coupled-perturbed
      !! solve uses, and for the same reason: past about a dozen the
      !! accumulator stops fitting in cache and the reuse turns into a loss.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: d1(:, :, :, :)
      real(dp), intent(in) :: bounds(:, :)
      real(dp), intent(in) :: zero_h(:, :)
      real(dp), allocatable, intent(out) :: g1(:, :, :)
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional :: xc
      real(dp), intent(in), optional :: reference(:, :)
      real(dp), intent(in), optional :: k_scale
         !! As the response operator's. The relaxed mean field entering
         !! `dW/dR` is the same object the coupled-perturbed operator
         !! applies, so it needs the same exchange fraction and the same
         !! kernel -- getting one right and not the other is a Hessian
         !! that is wrong while looking symmetric and plausible.
      real(dp), intent(in), optional :: rs_k_lr, rs_omega

      integer, parameter :: MAX_BATCH = 12

      real(dp), allocatable :: chunk(:, :, :), out(:, :, :), vxc_chunk(:, :)
      real(dp), allocatable :: out_lr(:, :, :)
      type(direct_stats_t) :: stats
      integer :: nao, natm, n_pert, first, last, wide, p, q, ia, a, ic
      integer :: n_chunks, per_chunk

      if (error%has_error()) return

      ! Both halves of the kernel or neither, as the coupled-perturbed solver
      ! demands: with only one supplied the loop below would quietly build a
      ! mean field without the exchange-correlation response, and `dW/dR` would
      ! then disagree with the operator that produced the orbitals it uses.
      if (present(xc) .neqv. present(reference)) then
         call error%set(ERROR_VALIDATION, "the relaxed mean field was given an "// &
                        "exchange-correlation context without the reference density its "// &
                        "kernel is evaluated at, or the reverse; it needs both or neither")
         return
      end if

      nao = size(d1, 1)
      natm = size(d1, 4)
      n_pert = 3*natm
      allocate (g1(nao, nao, n_pert))

      n_chunks = (n_pert + MAX_BATCH - 1)/MAX_BATCH
      per_chunk = (n_pert + n_chunks - 1)/n_chunks

      first = 1
      do while (first <= n_pert)
         last = min(first + per_chunk - 1, n_pert)
         wide = last - first + 1
         allocate (chunk(nao, nao, wide))
         do q = 1, wide
            p = first + q - 1
            ia = (p - 1)/3 + 1
            a = p - 3*(ia - 1)
            chunk(:, :, q) = d1(:, :, a, ia)
         end do
         call build_fock_direct_many(mol, zero_h, chunk, bounds, out, stats, error, &
                                     k_scale=k_scale)
         if (error%has_error()) return
         if (present(rs_omega)) then
            if (rs_omega > 0.0_dp) then
               if (.not. allocated(out_lr)) allocate (out_lr, mold=out)
               call build_fock_direct_many(mol, zero_h, chunk, bounds, out_lr, stats, &
                                           error, k_scale=rs_k_lr, j_scale=0.0_dp, &
                                           omega=rs_omega)
               if (error%has_error()) return
               out = out + out_lr
            end if
         end if

         ! The same kernel the response operator applies. `dW/dR` is built from
         ! this mean field, so leaving it out here while including it there
         ! makes the two halves of the response disagree about what the Fock
         ! operator is.
         if (present(xc) .and. present(reference)) then
            if (.not. allocated(vxc_chunk)) allocate (vxc_chunk(size(chunk, 1), size(chunk, 2)))
            do ic = 1, size(chunk, 3)
               vxc_chunk = 0.0_dp
               call xc_kernel_apply(xc, mol, reference, chunk(:, :, ic), vxc_chunk, error)
               if (error%has_error()) return
               out(:, :, ic) = out(:, :, ic) + vxc_chunk
            end do
         end if
         g1(:, :, first:last) = out
         deallocate (chunk, out)
         first = last + 1
      end do
   end subroutine mean_field_batch

end module mqc_libcint_hessian
