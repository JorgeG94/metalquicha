!! The analytic Hessian, for a restricted Hartree-Fock reference
module mqc_czt_hessian
   !! Second derivatives of the energy with respect to nuclear coordinates,
   !! without displacing anything. Replaces a `6N+1` finite-difference Hessian
   !! and the `1/h` amplification of gradient noise that goes with it, which
   !! lands hardest on the low-frequency modes the thermochemistry depends on.
   !!
   !! **Built in the pieces the standard decomposition uses**, which is also
   !! how PySCF's `hessian.rhf` is arranged, so the two can be compared stage
   !! by stage: the nuclear repulsion term, the per-atom perturbation that
   !! drives the coupled-perturbed equations, the response solve itself, and
   !! the explicit second-derivative assembly.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, atom_ao_blocks
   use mqc_czt_gradient, only: one_electron_deriv, iprinv_deriv_at, &
                               DERIV_KIN, DERIV_NUC, DERIV_OVLP
   use mqc_czt_hess_ints, only: eri_ip1_block, hess_1e_block, hess_2e_block, &
                                hess_rinv_block, hess_2e_contract, h1_contract, &
                                HESS_OVLP_II, HESS_OVLP_IJ, HESS_KIN_II, HESS_KIN_IJ, &
                                HESS_NUC_II, HESS_NUC_IJ, HESS_RINV_II, HESS_RINV_IJ, &
                                HESS_ERI_II, HESS_ERI_IJ, HESS_ERI_IK
   use mqc_czt_xc, only: xc_context_t, xc_kernel_apply, xc_add_potential, &
                         vv10_kernel_apply
   use mqc_czt_xc_hessian, only: xc_potential_deriv, vv10_potential_deriv
   use mqc_czt_direct, only: build_fock_direct, build_fock_direct_many, &
                             schwarz_bounds, direct_stats_t
   use mqc_czt_response, only: response_operator_t, solve_response
   use pic_blas_interfaces, only: pic_gemm
   use mqc_czt_ecp, only: ecp_refuses_derivatives
   ! TODO(mqc): `eri_ip1_block`, `hess_2e_block`, `HESS_ERI_II`, `HESS_ERI_IJ`
   ! and `HESS_ERI_IK` above are imported and used nowhere; the two-electron
   ! work moved to `hess_2e_contract` and `h1_contract` and the imports stayed.
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
      type(czt_molecule_t), pointer :: mol => null()
      real(dp), allocatable :: orbitals(:, :)     !! (n_ao, n_mo)
      real(dp), allocatable :: c_occ(:, :)        !! (n_ao, n_occ)
      real(dp), allocatable :: energies(:)        !! (n_mo)
      real(dp), allocatable :: bounds(:, :)
      real(dp), allocatable :: zero_h(:, :)
      integer :: n_occ = 0
      integer :: n_mo = 0
      integer :: n_pert = 1
         !! How many perturbations travel together; the trial vector is `n_mo`
         !! by `n_occ` by this. A Fock build costs an integral pass whether it
         !! contracts one density or a dozen.
      type(xc_context_t), pointer :: xc => null()
         !! The exchange-correlation context when the reference was Kohn-Sham,
         !! null for Hartree-Fock, which is what makes the kernel term opt-in.
      real(dp), allocatable :: reference(:, :)
         !! The converged density the kernel is evaluated at, which is not the
         !! trial density being contracted.
      real(dp) :: k_scale = 1.0_dp
         !! Exact exchange in the response operator. One for Hartree-Fock, zero
         !! for a pure functional, the mixing fraction for a hybrid.
      real(dp) :: rs_k_lr = 0.0_dp
      real(dp) :: rs_omega = 0.0_dp
         !! A range-separated functional's second exchange term:
         !! `rs_k_lr` of the exchange built against `erf(omega r)/r`. Zero
         !! omega is the ordinary case and costs nothing.
   contains
      procedure :: apply => nuclear_apply
      procedure :: length => nuclear_length
   end type nuclear_response_t

contains

   subroutine hcore_deriv_atom(mol, iatom, hcore_a, error)
      !! `dH_core/dR_A`, the one-electron Hamiltonian moved by one atom
      !!
      !! Two contributions that look alike and are not. Moving atom `A` moves
      !! the basis functions centred on it, which touches that atom's block
      !! alone. It also moves the **nucleus**, and every electron feels that
      !! wherever its orbital sits -- the Hellmann-Feynman term, which involves
      !! no basis-function derivative at all.
      !!
      !! Symmetrised at the end because the library puts the nabla on the bra,
      !! so what comes back is one half of a derivative symmetric in its two
      !! indices. The sign is the library's: the derivative with respect to the
      !! atom the bra sits on is minus its `ip` integral.
      type(czt_molecule_t), intent(in) :: mol
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

   subroutine make_h1_atom(mol, density, eri_ip1, iatom, h1, error, k_scale)
      !! The perturbation that drives the coupled-perturbed equations
      !!
      !!     h1_A = dH_core/dR_A + dV_HF/dR_A
      !!
      !! with the mean field differentiated with respect to the same atom. It
      !! needs **first** derivatives only -- the second derivatives belong to
      !! the explicit part of the Hessian, not here.
      !!
      !! **A quartet contributes through every index that sits on this atom.**
      !! `int2e_ip1` puts the nabla on the first index alone, so the derivative
      !! with respect to atom `A` is assembled by permuting each of the four
      !! positions into first place in turn and keeping the ones whose orbital
      !! belongs to `A`. The permutations used are the ordinary symmetries of
      !! an undifferentiated integral -- `(ij|kl) = (ji|kl) = (kl|ij)` -- which
      !! hold because the nabla is applied *after* the permutation.
      !!
      !! The sign is the library's, as everywhere else here.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)          !! Total AO density
      real(dp), intent(in) :: eri_ip1(:, :, :, :, :)  !! From `eri_ip1_block`
      integer, intent(in) :: iatom
      real(dp), allocatable, intent(out) :: h1(:, :, :)   !! (n_ao, n_ao, 3)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: k_scale
         !! Exact-exchange fraction of the mean field being differentiated.
         !! Absent is one, Hartree-Fock; a Kohn-Sham reference passes its
         !! functional's fraction. The exchange-correlation potential's own
         !! derivative is a grid quantity and no business of this routine --
         !! `xc_potential_deriv` accumulates it on top, at the caller.

      real(dp), allocatable :: hcore_a(:, :, :), vhf(:, :, :)
      integer, allocatable :: offsets(:), counts(:), owner(:)
      integer :: nao, comp, mu, nu, la, si, a, p0, p1
      real(dp) :: d, khalf

      if (error%has_error()) return
      khalf = 0.5_dp
      if (present(k_scale)) khalf = 0.5_dp*k_scale

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

      ! J - K/2, differentiated, as a plain quadruple loop over the whole basis.
      ! `h1_contract` is the shell-driven form the assembly actually uses.
      !
      ! **Two copies of the loop, split on the exchange fraction.** A variable
      ! multiplier on the exchange terms shifts the compiler's contraction and
      ! scheduling by an ulp even at a fraction of exactly one -- measured on
      ! this Hessian's bit-for-bit regression -- so full exchange keeps the
      ! literal statements it always had, and only a genuinely scaled
      ! reference pays for the variable.
      if (khalf == 0.5_dp) then
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
      else
         do comp = 1, 3
            do si = 1, nao
               do la = 1, nao
                  do nu = 1, nao
                     do mu = 1, nao
                        d = density(la, si)
                        if (abs(d) < 1.0e-14_dp) cycle
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
                        ! Exchange, at the fraction the reference kept.
                        if (owner(mu) == iatom) then
                           vhf(mu, nu, comp) = vhf(mu, nu, comp) &
                                               + khalf*d*eri_ip1(mu, la, si, nu, comp)
                        end if
                        if (owner(la) == iatom) then
                           vhf(mu, nu, comp) = vhf(mu, nu, comp) &
                                               + khalf*d*eri_ip1(la, mu, si, nu, comp)
                        end if
                        if (owner(si) == iatom) then
                           vhf(mu, nu, comp) = vhf(mu, nu, comp) &
                                               + khalf*d*eri_ip1(si, nu, mu, la, comp)
                        end if
                        if (owner(nu) == iatom) then
                           vhf(mu, nu, comp) = vhf(mu, nu, comp) &
                                               + khalf*d*eri_ip1(nu, si, mu, la, comp)
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end if

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
      !! leaves the overlap unchanged; displacing a nucleus drags its functions
      !! with it, and the orbitals must stay orthonormal while that happens.
      !! That constraint is what puts a non-zero occupied-occupied block into
      !! the orbital response.
      !!
      !! Only functions centred on `A` move, so the derivative is the
      !! differentiated overlap restricted to that atom's rows, plus its
      !! transpose in the same columns. The sign is the library's.
      type(czt_molecule_t), intent(in) :: mol
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
      !! The trial density is **symmetrised**, which is what makes the fast Fock
      !! build correct here -- the antisymmetrised variant belongs to the
      !! frequency-dependent `A - B` block and needs a different build entirely.
      !!
      !! **The occupied rows are zeroed rather than solved for.** They are
      !! already known -- orthonormality fixes them -- so the iteration must not
      !! move them, and the solver puts them back from the right-hand side on
      !! every pass. The virtual rows are divided by their orbital energy gaps.
      class(nuclear_response_t), intent(inout) :: this
      real(dp), intent(in) :: vector(:)
      real(dp), intent(out) :: image(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: mo1(:, :, :), half(:, :), dens(:, :, :), g(:, :, :)
      real(dp), allocatable :: vxc(:, :), vnl(:, :, :)
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
         ! Twice, for double occupancy: the trial vector describes how one
         ! spatial orbital moves, and the density it perturbs holds two
         ! electrons in it. Leaving the factor out still converges, to a wrong
         ! answer that is not obviously doubled because the error feeds back
         ! through the coupling.
         call pic_gemm(this%orbitals, mo1(:, :, p), half, alpha=2.0_dp, beta=0.0_dp)
         call pic_gemm(half, this%c_occ, dens(:, :, p), transb="T")
         dens(:, :, p) = dens(:, :, p) + transpose(dens(:, :, p))
      end do

      ! One integral pass for the whole batch. `build_fock_direct_many` screens
      ! on the Schwarz bound alone, so a batch is bit-for-bit what the same
      ! densities give one at a time -- which matters here, because a response
      ! density is not a density and screening on its magnitude would be wrong.
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

      ! The exchange-correlation kernel, for a Kohn-Sham reference: the second
      ! derivative of E_xc with respect to the density, contracted against the
      ! trial density. Leaving it out does not fail, it converges to the wrong
      ! orbital response.
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
         ! The non-local kernel, once for the whole batch rather than inside
         ! the loop above: `vv10_kernel_apply`'s pair sweep is O(npts^2)
         ! whether it carries one trial or a dozen. It accumulates, like the
         ! semilocal kernel, hence the zeroed buffer of its own.
         if (this%xc%nlc_b /= 0.0_dp .or. this%xc%nlc_c /= 0.0_dp) then
            allocate (vnl(size(dens, 1), size(dens, 2), size(dens, 3)))
            vnl = 0.0_dp
            call vv10_kernel_apply(this%xc, this%mol, this%reference, dens, vnl, error)
            if (error%has_error()) return
            g = g + vnl
            deallocate (vnl)
         end if
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
      !! `solve_mo1_batch` solves the same equations for every atom at once and
      !! is what the assembly uses.
      !!
      !! **The occupied-occupied block is not solved for.** Orthonormality fixes
      !! it at minus half the overlap derivative, so it goes into the
      !! right-hand side and the operator zeroes it on every pass. It is carried
      !! rather than dropped because the virtual-occupied block couples to it
      !! through the density the whole vector induces.
      !!
      !! The right-hand side is `h1 - s1 e_i` divided by the orbital energy
      !! gap, with the sign that makes the response enter the fixed point with
      !! the opposite one.
      type(czt_molecule_t), intent(in), target :: mol
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
      ! orbitals: one pass over shell pairs per atom where one per molecule
      ! would do. `solve_mo1_batch` builds them once.
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
      !! and the diagonal block is minus the sum of the others in its row, which
      !! imposes rather than hopes for translational invariance.
      !!
      !! `(3, 3, n_atoms, n_atoms)` with the Cartesian pair innermost, the
      !! layout the electronic part adds into.
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
      !! **The one-electron terms are assembled by depositing rather than by
      !! slicing.** Each derivative belongs to whichever atom the differentiated
      !! basis function sits on, so the loop runs over all indices once and adds
      !! each term into the atom pair its own indices name, which keeps the
      !! ownership rule next to the term rather than in a slice bound.
      !!
      !! The two-electron terms go through `hess_2e_contract`, which drives the
      !! same deposit from shells and never forms an `nao^4` array.
      type(czt_molecule_t), intent(in) :: mol
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
      ! functions moved.
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
      ! `-Z_C / |r - R_C|` depends on three positions: the bra centre, the ket
      ! centre and the nucleus. Differentiating the nucleus produces no basis
      ! derivative at all, so it reaches atom pairs the block above cannot see
      ! -- a nucleus with no basis functions on it still has a Hessian row.
      !
      ! `cross(a, b, A)` collects, for the atom `c` whose nucleus is moving, the
      ! part belonging to basis functions centred on `A`. Bra and ket fold onto
      ! each other, hence the factor of two.
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
      ! Driven from shells and contracted on the spot: the three arrays this
      ! would otherwise form are `nao^4` times nine each.
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
                               max_iter, tol, dipole_derivatives)
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
      !! two are the same object -- `D F D / 2 = 2 C_occ eps C_occ^T`, because
      !! `C_occ^T F C_occ` is diagonal at convergence -- but the first form's
      !! derivative needs only `dD/dR` and `dF/dR`, both already here, instead
      !! of the derivative of the Lagrange multiplier matrix.
      !!
      !! The coupled-perturbed solve is where essentially all of the time goes.
      type(czt_molecule_t), intent(in), target :: mol
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
      real(dp), allocatable, intent(out), optional :: dipole_derivatives(:, :)
         !! `d mu_a / dR_(X,b)`, `(3, 3*natm)` in atomic units, the Cartesian
         !! index fastest -- the layout `finite_diff_dipole_derivatives`
         !! produces and `compute_vibrational_analysis` reads for infrared
         !! intensities.
         !!
         !! Here rather than in a routine of its own because the expensive part
         !! is already done. What a dipole derivative needs beyond one-electron
         !! integrals is the perturbed density, and that is exactly what the
         !! coupled-perturbed solve below produces for the Hessian. Asking for
         !! it separately would pay for the whole response a second time, which
         !! is most of what an analytic Hessian costs.

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
      ! Hartree-Fock matrix here spreads a smooth, symmetric error across every
      ! atom pair.
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
      ! quartets. The core Hamiltonian derivative is per-atom either way and is
      ! added below.
      call h1_contract(mol, density, h1, error, k_scale=k_scale)
      if (error%has_error()) return
      if (present(rs_omega)) then
         if (rs_omega > 0.0_dp) then
            ! Into a temporary and added: `h1_contract` zeroes its output where
            ! `hess_2e_contract` next door accumulates, so a second call in
            ! place would discard the full-range pass.
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
      ! produce. Both of these accumulate into `h1`, and `vv10_potential_deriv`
      ! returns untouched when the functional carries no VV10.
      if (present(xc)) then
         call xc_potential_deriv(xc, mol, density, h1, error)
         if (error%has_error()) return
         call vv10_potential_deriv(xc, mol, density, h1, error)
      end if
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
            ! The first-order density.
            call pic_gemm(orbitals, mo1(:, :, a, ia), half, alpha=2.0_dp, beta=0.0_dp)
            call pic_gemm(half, c_occ, tmp, transb="T")
            tmp = tmp + transpose(tmp)
            d1(:, :, a, ia) = tmp
         end do
      end do

      ! The dipole derivative, while `d1` is in hand and before the mean-field
      ! batch reuses the buffers. Nothing below needs it and it needs nothing
      ! below.
      if (present(dipole_derivatives)) then
         call assemble_dipole_derivatives(mol, density, d1, dipole_derivatives, error)
         if (error%has_error()) return
      end if

      ! Their mean fields in one batch, chunked as the solve above is.
      call mean_field_batch(mol, d1, bounds, zero_h, g1all, error, &
                            xc=xc, reference=reference, k_scale=k_scale, &
                            rs_k_lr=rs_k_lr, rs_omega=rs_omega)
      if (error%has_error()) return

      do ia = 1, natm
         do a = 1, 3
            tmp = d1(:, :, a, ia)
            ! `dF/dR` in full: the skeleton derivative plus the mean field of
            ! the relaxed density.
            f1 = h1(:, :, a, ia) + g1all(:, :, 3*(ia - 1) + a)

            ! W = D F D / 2, differentiated by the product rule: three terms,
            ! all of which get the half. The first and last are transposes of
            ! each other and the middle one is already symmetric, so it must
            ! **not** be symmetrised a second time -- doing so doubles its
            ! weight and gives a `dW/dR` that is still symmetric, still
            ! translationally invariant, and wrong.
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
                          hess, error, max_iter, tol, dipole_derivatives)
      !! The analytic Hessian of a converged restricted Hartree-Fock energy
      !!
      !! Three pieces, only meaningful added together: the nuclear repulsion,
      !! the explicit integral second derivatives against the fixed density,
      !! and the response.
      !!
      !! Returned as `(3, 3, n_atoms, n_atoms)`. A vibrational analysis wants
      !! `(3*n_atoms, 3*n_atoms)`, which `hessian_to_matrix` produces.
      type(czt_molecule_t), intent(in), target :: mol
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: hess(:, :, :, :)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      real(dp), allocatable, intent(out), optional :: dipole_derivatives(:, :)
         !! `d mu_a / dR_(X,b)`, `(3, 3*natm)`, forwarded from
         !! `response_hessian` -- see there for why it comes out beside the
         !! Hessian rather than from a routine of its own. Present, an infrared
         !! intensity costs two one-electron integrals on top of a Hessian that
         !! was being computed anyway; absent, nothing extra is built.

      real(dp), allocatable :: weighted(:, :), part(:, :, :, :), resp(:, :, :, :)
      integer :: i, nao

      if (error%has_error()) return
      if (ecp_refuses_derivatives(mol%core_electrons, "nuclear Hessian", error)) return

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
                            max_iter=max_iter, tol=tol, &
                            dipole_derivatives=dipole_derivatives)
      if (error%has_error()) return

      hess = hess + part + resp
      deallocate (weighted, part, resp)
   end subroutine rhf_hessian

   subroutine ks_hessian(mol, atomic_numbers, density, orbitals, energies, n_occ, &
                         xc, k_scale, hess, error, max_iter, tol, dipole_derivatives)
      !! The analytic Hessian of a converged Kohn-Sham energy
      !!
      !! `rhf_hessian` with the exchange-correlation parts put back. Three
      !! things change:
      !!
      !! **Exact exchange is scaled.** A pure functional has none, so the
      !! exchange half of the explicit two-electron term and of the response
      !! operator both vanish; a hybrid keeps its mixing fraction.
      !!
      !! **The exchange-correlation second derivatives are added** -- the
      !! explicit term, `xc_hessian`.
      !!
      !! **The response operator gains the kernel.** Without it the
      !! coupled-perturbed solve still converges, to the wrong orbital response.
      !!
      !! The grid is held fixed throughout, as `xc_hessian` sets out, which is
      !! `grid_response = False` in PySCF's terms.
      !!
      !! Covers LDA, GGA, meta-GGA, hybrids of each, range-separated hybrids,
      !! and functionals carrying VV10 non-local correlation.
      ! TODO(mqc): no `ecp_refuses_derivatives` guard, where `rhf_hessian` has
      ! one, so a Kohn-Sham Hessian on a molecule with an effective core
      ! potential returns wrong derivatives instead of refusing.
      use mqc_czt_xc_hessian, only: xc_hessian, vv10_hessian
      type(czt_molecule_t), intent(in), target :: mol
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
      real(dp), allocatable, intent(out), optional :: dipole_derivatives(:, :)
         !! `d mu_a / dR_(X,b)`, `(3, 3*natm)`, forwarded from
         !! `response_hessian` -- see there for why it comes out beside the
         !! Hessian rather than from a routine of its own. Present, an infrared
         !! intensity costs two one-electron integrals on top of a Hessian that
         !! was being computed anyway; absent, nothing extra is built.

      real(dp), allocatable :: weighted(:, :), part(:, :, :, :), resp(:, :, :, :)
      integer :: i, nao

      if (error%has_error()) return
      if (ecp_refuses_derivatives(mol%core_electrons, "nuclear Hessian", error)) return

      ! VV10 contributes in the same three places as the semilocal term:
      ! `vv10_hessian` in the explicit part alongside `xc_hessian`,
      ! `vv10_potential_deriv` in the perturbed Fock alongside
      ! `xc_potential_deriv`, and `vv10_kernel_apply` in the response operator
      ! and the relaxed mean field alongside `xc_kernel_apply`.
      !
      ! The NLC grid is held fixed like the semilocal one, which departs from
      ! PySCF: its `_get_vnlc_deriv1` hard-codes the NLC grid response even
      ! when the rest of its Hessian omits grid response, so the two Hessians
      ! differ by exactly that term, shrinking as the NLC grid is refined.
      !
      ! TODO(mqc): `vv10_kernel_apply` rebuilds kernel intermediates that
      ! depend only on the SCF density on every call, so the solve pays one of
      ! its two pair sweeps per iteration for nothing.

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
      ! The non-local term's explicit second derivative, on the NLC grid.
      ! Accumulates into `part` as `xc_hessian` does, and returns untouched
      ! for a functional without VV10.
      call vv10_hessian(xc, mol, density, part, error)
      if (error%has_error()) return
      call response_hessian(mol, density, orbitals, energies, n_occ, resp, error, &
                            xc=xc, reference=density, k_scale=k_scale, &
                            rs_k_lr=xc%rs_k_lr, rs_omega=xc%rs_omega, &
                            max_iter=max_iter, tol=tol, &
                            dipole_derivatives=dipole_derivatives)
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
      !! **The rule to keep is the layout, not a symmetry argument about it.**
      !! Swapping both index pairs is the genuine transpose and changes nothing,
      !! but swapping only the atoms, or interleaving the two indices
      !! Cartesian-slowest, moves every frequency.
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
      !! contracts.
      !!
      !! **Chunked at `MAX_BATCH` rather than all at once.** The win saturates
      !! near fourfold and then *reverses*: the accumulator that makes the reuse
      !! possible grows with the batch while the integral it reuses does not,
      !! so past about a dozen densities the build becomes memory-bound.
      !!
      !! The Schwarz bounds are built once here rather than once per atom.
      ! TODO(mqc): `MAX_BATCH` is spelled twice, here and in
      ! `mean_field_batch`, for a figure both routines document as coming from
      ! the same cache limit. Changing one leaves the other behind.
      type(czt_molecule_t), intent(in), target :: mol
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

      ! `xc` and `reference` are one argument in two halves: the context without
      ! the density its kernel is evaluated at would reach `xc_kernel_apply`
      ! through `nuclear_apply` with `operator%reference` unallocated.
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

      ! Balanced rather than greedy: eighteen perturbations taken twelve at a
      ! time leaves a chunk of six that pays a full integral pass for half the
      ! reuse, where nine and nine costs the same two passes and amortises both.
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
      !! Flattens `(n_ao, n_ao, 3, natm)` to `(n_ao, n_ao, 3*natm)` and hands it
      !! to the many-density build in the same chunks, and for the same reason,
      !! as `solve_mo1_batch`.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: d1(:, :, :, :)
      real(dp), intent(in) :: bounds(:, :)
      real(dp), intent(in) :: zero_h(:, :)
      real(dp), allocatable, intent(out) :: g1(:, :, :)
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional :: xc
      real(dp), intent(in), optional :: reference(:, :)
      real(dp), intent(in), optional :: k_scale
         !! **The response operator's, and it has to match.** The relaxed mean
         !! field entering `dW/dR` is the same object the coupled-perturbed
         !! operator applies; getting one right and not the other gives a
         !! Hessian that is wrong while looking symmetric and plausible.
      real(dp), intent(in), optional :: rs_k_lr, rs_omega

      integer, parameter :: MAX_BATCH = 12

      real(dp), allocatable :: chunk(:, :, :), out(:, :, :), vxc_chunk(:, :)
      real(dp), allocatable :: vnl_chunk(:, :, :), out_lr(:, :, :)
      type(direct_stats_t) :: stats
      integer :: nao, natm, n_pert, first, last, wide, p, q, ia, a, ic
      integer :: n_chunks, per_chunk

      if (error%has_error()) return

      ! Both halves of the kernel or neither: with only one supplied the loop
      ! below builds a mean field without the exchange-correlation response, and
      ! `dW/dR` then disagrees with the operator that produced its orbitals.
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

         ! The same kernel the response operator applies, for the reason
         ! `k_scale` above gives.
         if (present(xc) .and. present(reference)) then
            if (.not. allocated(vxc_chunk)) allocate (vxc_chunk(size(chunk, 1), size(chunk, 2)))
            do ic = 1, size(chunk, 3)
               vxc_chunk = 0.0_dp
               call xc_kernel_apply(xc, mol, reference, chunk(:, :, ic), vxc_chunk, error)
               if (error%has_error()) return
               out(:, :, ic) = out(:, :, ic) + vxc_chunk
            end do
            ! The non-local kernel, once per chunk for the reason
            ! `nuclear_apply` gives.
            if (xc%nlc_b /= 0.0_dp .or. xc%nlc_c /= 0.0_dp) then
               allocate (vnl_chunk(size(chunk, 1), size(chunk, 2), size(chunk, 3)))
               vnl_chunk = 0.0_dp
               call vv10_kernel_apply(xc, mol, reference, chunk, vnl_chunk, error)
               if (error%has_error()) return
               out = out + vnl_chunk
               deallocate (vnl_chunk)
            end if
         end if
         g1(:, :, first:last) = out
         deallocate (chunk, out)
         first = last + 1
      end do
   end subroutine mean_field_batch

   subroutine assemble_dipole_derivatives(mol, density, d1, ddip_dr, error)
      !! `d mu_a / dR_(X,b)` from the perturbed densities the response produced
      !!
      !! The dipole is nuclear less electronic,
      !!
      !!     mu_a = sum_A Z_A R_(A,a) - Tr(D r_a)
      !!
      !! so differentiating it gives three terms and every one of them matters:
      !!
      !!     d mu_a / dR_(X,b) = Z_X delta_ab
      !!                       - Tr(D  d r_a / dR_(X,b))
      !!                       - Tr(D^(X,b) r_a)
      !!
      !! the first because the nucleus itself moved, the second because the
      !! basis functions ride on it, and the third because the electrons relax.
      !! Dropping the third is the same mistake as dropping the response from a
      !! Hessian: what is left is smooth, plausible, and wrong by tens of
      !! percent on exactly the modes an infrared spectrum is read for.
      !!
      !! **`mol%charges` is the right charge here**, ECP-reduced and all. It is
      !! not the atomic number, and elsewhere in this backend that distinction
      !! is a bug -- the exchange-correlation grid wants the element and had to
      !! have `core_electrons` added back. Here the reduced charge is correct
      !! and adding the core back would be the error: the electrons those cores
      !! hold are not in `density` either, so the two halves of the dipole have
      !! to be missing the same thing for the neutral molecule to come out
      !! neutral.
      use mqc_czt_multipole, only: multipole_matrices, dipole_integral_derivatives
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: density(:, :)      !! Closed shell, carrying its two
      real(dp), intent(in) :: d1(:, :, :, :)     !! (n_ao, n_ao, 3, natm)
      real(dp), allocatable, intent(out) :: ddip_dr(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dip(:, :, :), ddip(:, :, :, :, :)
      real(dp) :: origin(3)
      integer :: natm, ia, a, b, col

      natm = mol%natm

      ! The origin cancels out of a derivative -- shifting it moves every
      ! column by the same constant times the total charge, and that constant
      ! is what the translational sum rule removes -- so this is a choice about
      ! conditioning rather than about physics, and the nuclear centroid keeps
      ! `r` small over the basis.
      origin = 0.0_dp
      do ia = 1, natm
         origin = origin + mol%coords(:, ia)
      end do
      if (natm > 0) origin = origin/real(natm, dp)

      call multipole_matrices(mol, origin, 1, dip, error)
      if (error%has_error()) return
      call dipole_integral_derivatives(mol, origin, ddip, error)
      if (error%has_error()) return

      allocate (ddip_dr(3, 3*natm))
      ddip_dr = 0.0_dp

      do ia = 1, natm
         do b = 1, 3
            col = 3*(ia - 1) + b
            do a = 1, 3
               ddip_dr(a, col) = -sum(density*ddip(:, :, a, b, ia)) &
                                 - sum(d1(:, :, b, ia)*dip(:, :, a))
            end do
            ddip_dr(b, col) = ddip_dr(b, col) + mol%charges(ia)
         end do
      end do

      deallocate (dip, ddip)
   end subroutine assemble_dipole_derivatives

end module mqc_czt_hessian
