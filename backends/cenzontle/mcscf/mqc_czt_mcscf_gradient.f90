!! Analytic nuclear gradients for an optimised MCSCF wave function
module mqc_czt_mcscf_gradient
   !! dE/dR for a converged CASSCF or ORMAS-SCF, in Hartree/Bohr
   !!
   !! **No Z-vector.** A converged MCSCF is stationary with respect to both the
   !! orbital rotations and the CI coefficients, so the response terms are
   !! identically zero and the gradient is a contraction of differentiated
   !! integrals against densities that are already in hand.
   !!
   !! That stationarity is a precondition rather than a property: at a point
   !! where the orbital gradient is not zero this returns a confident number that
   !! is not the derivative of anything.
   !!
   !! **The two-electron term is split by how many active indices it carries.**
   !! Contracting the whole two-particle density in the AO basis is `n_ao^4` for
   !! every term, including the mean-field ones. The inactive-inactive and
   !! inactive-active blocks are separable, being products of one-particle
   !! densities, so
   !!
   !! ```
   !! E_2e = 1/2 D.(J - K/2)(D)  +  1/2 sum_tuvw  ddm2_tuvw (tu|vw)
   !! ```
   !!
   !! with `D` the total one-particle density and `ddm2` the active
   !! two-particle density with its own mean-field part removed,
   !!
   !! ```
   !! ddm2_tuvw = dm2_tuvw - dm1_tu dm1_vw + 1/2 dm1_tv dm1_uw
   !! ```
   !!
   !! The first term is what a closed-shell SCF gradient contracts and goes
   !! through the same routine. Only the second needs the general four-index
   !! machinery, and it is built from `n_active^4` rather than from the basis.
   !!
   !! Subtracting that mean-field part is not optional bookkeeping: leave it in
   !! and the active-active energy is counted twice. The failure is invisible at
   !! `n_active = 0`, where `ddm2` is empty and this reduces exactly to the
   !! closed-shell gradient.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, shell_dim, atom_ao_blocks
   use mqc_czt_gradient, only: nuclear_repulsion_gradient, one_electron_deriv, &
                               iprinv_deriv_at, two_electron_deriv, &
                               DERIV_OVLP, DERIV_KIN, DERIV_NUC
   use mqc_czt_mp2_gradient, only: two_electron_mp2_terms
   use mqc_czt_mcscf, only: generalized_fock, mcscf_fock_t
   use pic_blas_interfaces, only: pic_gemm
   implicit none
   private

   public :: czt_mcscf_gradient
   public :: cumulant_two_particle_density   !! Exposed for the tests

   real(dp), parameter :: BLOCK_TARGET = 2.0e8_dp
      !! Bytes the first-index block of the AO two-particle density may reach, so
      !! that the ceiling is set by choice rather than by the basis.

contains

   subroutine czt_mcscf_gradient(mol, orbitals, n_inactive, n_active, dm1, dm2, &
                                 gradient, error)
      !! The gradient of a stationary MCSCF energy
      !!
      !! `dm1` and `dm2` are the spin-traced active densities as
      !! `active_space_rdms` builds them, and they must be the ones belonging to
      !! `orbitals`. On a run that stopped early the two are one orbital step
      !! apart -- consistent with each other but not with the orbitals -- so the
      !! caller checks convergence before arriving here.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)       !! (n_ao, n_mo)
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: dm1(:, :)            !! (n_active, n_active)
      real(dp), intent(in) :: dm2(:, :, :, :)      !! (n_active^4), chemist order
      real(dp), allocatable, intent(out) :: gradient(:, :)   !! (3, n_atoms)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: d_inactive(:, :), d_active(:, :), d_total(:, :)
      real(dp), allocatable :: weighted(:, :), work(:, :)
      real(dp), allocatable :: c_active(:, :), ddm2(:, :, :, :)
      real(dp), allocatable :: s1(:, :, :), kin(:, :, :), h1(:, :, :)
      real(dp), allocatable :: vhf(:, :, :), vrinv(:, :, :), hcore_a(:, :, :)
      type(mcscf_fock_t) :: fock
      integer, allocatable :: offsets(:), counts(:)
      integer :: n_ao, n_mo, iatom, comp, p0, p1

      if (error%has_error()) return

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)

      if (n_ao /= mol%nao) then
         call error%set(ERROR_VALIDATION, "the MCSCF gradient was given orbitals "// &
                        "that do not span this basis.")
         return
      end if
      if (n_inactive + n_active > n_mo) then
         call error%set(ERROR_VALIDATION, "more inactive plus active orbitals than "// &
                        "the basis provides.")
         return
      end if
      if (n_active > 0) then
         if (size(dm1, 1) /= n_active .or. size(dm2, 1) /= n_active) then
            call error%set(ERROR_VALIDATION, "the active densities do not match the "// &
                           "active space they are said to belong to.")
            return
         end if
      end if

      allocate (gradient(3, mol%natm))
      gradient = 0.0_dp

      call nuclear_repulsion_gradient(mol, gradient)

      ! The one-particle density, in two pieces because the gradient needs them
      ! separately later and the sum immediately.
      allocate (d_inactive(n_ao, n_ao), d_active(n_ao, n_ao), d_total(n_ao, n_ao))
      d_inactive = 0.0_dp
      d_active = 0.0_dp
      if (n_inactive > 0) then
         call pic_gemm(orbitals(:, 1:n_inactive), orbitals(:, 1:n_inactive), d_inactive, &
                       transb="T", alpha=2.0_dp, beta=0.0_dp)
      end if
      if (n_active > 0) then
         allocate (c_active(n_ao, n_active))
         c_active = orbitals(:, n_inactive + 1:n_inactive + n_active)
         allocate (work(n_ao, n_active))
         call pic_gemm(c_active, dm1, work, beta=0.0_dp)
         call pic_gemm(work, c_active, d_active, transb="T", beta=0.0_dp)
         deallocate (work)
      end if
      d_total = d_inactive + d_active

      ! The Pulay term contracts against the generalised Fock matrix, which is
      ! the MCSCF energy-weighted density: for a closed shell with no active
      ! space it is `2 eps_i` on the diagonal, which is what the SCF gradient
      ! builds by hand from the orbital energies. Symmetrised because it is
      ! symmetric only at a stationary point, and rounding should not decide
      ! whether the Pulay term is.
      call generalized_fock(mol, orbitals, n_inactive, n_active, dm1, dm2, fock, error)
      if (error%has_error()) return

      allocate (weighted(n_ao, n_ao), work(n_ao, n_mo))
      call pic_gemm(orbitals, 0.5_dp*(fock%general + transpose(fock%general)), work, &
                    beta=0.0_dp)
      call pic_gemm(work, orbitals, weighted, transb="T", beta=0.0_dp)
      deallocate (work)

      ! The sign convention is libcint's, and the same one the SCF gradient
      ! states: `ip` integrals carry a nabla on the bra, and the derivative with
      ! respect to the atom the bra sits on is minus that.
      call one_electron_deriv(mol, s1, DERIV_OVLP)
      s1 = -s1
      call one_electron_deriv(mol, kin, DERIV_KIN)
      call one_electron_deriv(mol, h1, DERIV_NUC)
      h1 = -(kin + h1)
      deallocate (kin)

      ! The separable two-electron term: the inactive-inactive and
      ! inactive-active blocks together, which is the closed-shell contraction
      ! over the total density.
      call two_electron_deriv(mol, d_total, vhf, error)
      if (error%has_error()) return

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)

      allocate (vrinv(n_ao, n_ao, 3), hcore_a(n_ao, n_ao, 3))
      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)

         ! Moving an atom moves the basis functions centred on it and moves the
         ! nucleus every electron is attracted to. The first is the `h1` block,
         ! the second is Hellmann-Feynman and involves no basis derivative.
         call iprinv_deriv_at(mol, iatom, vrinv)
         vrinv = -mol%charges(iatom)*vrinv
         if (counts(iatom) > 0) then
            vrinv(p0:p1, :, :) = vrinv(p0:p1, :, :) + h1(p0:p1, :, :)
         end if
         do comp = 1, 3
            hcore_a(:, :, comp) = vrinv(:, :, comp) + transpose(vrinv(:, :, comp))
         end do

         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + sum(hcore_a(:, :, comp)*d_total)
         end do

         if (counts(iatom) == 0) cycle

         ! Twice, because the nabla sat on the bra and the ket's contribution is
         ! the same number by the symmetry of the undifferentiated integrals.
         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + 2.0_dp*sum(vhf(p0:p1, :, comp)*d_total(p0:p1, :)) &
                                    - 2.0_dp*sum(s1(p0:p1, :, comp)*weighted(p0:p1, :))
         end do
      end do
      deallocate (vrinv, hcore_a, vhf, s1, h1)

      ! What is left is the part of the active two-particle density that is not
      ! a product of one-particle densities. Nothing to do when there is no
      ! active space, and this module then reduces to the closed-shell gradient.
      if (n_active > 0) then
         call cumulant_two_particle_density(dm1, dm2, ddm2)
         call active_two_electron_gradient(mol, c_active, ddm2, gradient, error)
         if (error%has_error()) return
      end if
   end subroutine czt_mcscf_gradient

   pure subroutine cumulant_two_particle_density(dm1, dm2, ddm2)
      !! `dm2` with the part already counted as a mean-field product removed
      !!
      !! The separable term contracted above is `1/2 D.(J - K/2)(D)` over the
      !! *total* density, and expanding that leaves an active-active piece
      !! `1/2 sum [dm1_tu dm1_vw - 1/2 dm1_tv dm1_uw] (tu|vw)`. Subtracting it
      !! here keeps the two contractions a partition of the energy.
      !!
      !! Chemist ordering throughout, matching `active_space_rdms` and the
      !! `(tu|vw)` the generalised Fock contracts.
      real(dp), intent(in) :: dm1(:, :)
      real(dp), intent(in) :: dm2(:, :, :, :)
      real(dp), allocatable, intent(out) :: ddm2(:, :, :, :)

      integer :: n, t, u, v, w

      n = size(dm1, 1)
      allocate (ddm2(n, n, n, n))
      do w = 1, n
         do v = 1, n
            do u = 1, n
               do t = 1, n
                  ddm2(t, u, v, w) = dm2(t, u, v, w) - dm1(t, u)*dm1(v, w) &
                                     + 0.5_dp*dm1(t, v)*dm1(u, w)
               end do
            end do
         end do
      end do
   end subroutine cumulant_two_particle_density

   subroutine active_two_electron_gradient(mol, c_active, ddm2, gradient, error)
      !! The genuinely non-separable term, a block of the first index at a time
      !!
      !! The density is transformed out of the active space and into the AO
      !! basis one block at a time and handed to the same contraction the MP2
      !! gradient uses, which expects exactly this: a general four-index array
      !! that is *not* symmetric in the ket pair, and a block of the first index
      !! identified by its shell range.
      !!
      !! Blocks are cut on shell boundaries because a shell's functions share a
      !! quartet and cannot be split across two passes.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: c_active(:, :)       !! (n_ao, n_active)
      real(dp), intent(in) :: ddm2(:, :, :, :)
      real(dp), intent(inout) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: gamma_blk(:, :, :, :), eri_blk(:, :, :, :)
      real(dp), allocatable :: dummy_vhf1(:, :, :, :), hf_density(:, :)
      integer :: n_ao, ish, ish_lo, ish_hi, p_lo, p_hi, np, per_block

      if (error%has_error()) return
      n_ao = size(c_active, 1)

      ! `with_reference` false leaves these alone; they exist because an absent
      ! optional cannot be named where these are named.
      allocate (dummy_vhf1(1, 1, 1, 1), hf_density(n_ao, n_ao))
      dummy_vhf1 = 0.0_dp
      hf_density = 0.0_dp

      per_block = max(1, int(BLOCK_TARGET/(2.0_dp*real(n_ao, dp)**3*8.0_dp)))

      ish_lo = 1
      do while (ish_lo <= mol%nbas)
         p_lo = mol%shell_offset(ish_lo) + 1
         ish_hi = ish_lo
         do ish = ish_lo, mol%nbas
            p_hi = mol%shell_offset(ish) + shell_dim(mol%cartesian, ish - 1, mol%bas)
            if (ish > ish_lo .and. p_hi - p_lo + 1 > per_block) exit
            ish_hi = ish
         end do
         p_hi = mol%shell_offset(ish_hi) + shell_dim(mol%cartesian, ish_hi - 1, mol%bas)
         np = p_hi - p_lo + 1

         call gamma_block(c_active, ddm2, p_lo, p_hi, gamma_blk)

         allocate (eri_blk(np, n_ao, n_ao, n_ao))
         eri_blk = 0.0_dp
         call two_electron_mp2_terms(mol, gamma_blk, hf_density, gradient, dummy_vhf1, &
                                     ish_lo=ish_lo, ish_hi=ish_hi, &
                                     p_offset=p_lo - 1, eri_blk=eri_blk, &
                                     with_gamma=.true., with_reference=.false.)
         deallocate (eri_blk, gamma_blk)
         ish_lo = ish_hi + 1
      end do

      deallocate (dummy_vhf1, hf_density)
   end subroutine active_two_electron_gradient

   subroutine gamma_block(c_active, ddm2, p_lo, p_hi, gamma_blk)
      !! `ddm2` transformed to the AO basis over one range of the first index
      !!
      !! **Four gemms and no copies.** Each step contracts the *first* index and
      !! leaves the new one last, so after four the result is already `(p q r s)`
      !! and never has to be put in order. The copying is what costs here, not the
      !! arithmetic: the last intermediate alone is `np n_ao^3`.
      !!
      !! The buffers are therefore flat and viewed through pointers with their
      !! bounds remapped, as `build_two_particle_density` does next door. The
      !! memory order is identical either way, so a `reshape` between steps would
      !! only be telling the compiler what it already had.
      real(dp), intent(in), target :: c_active(:, :)
      real(dp), intent(in), target, contiguous :: ddm2(:, :, :, :)
         !! `contiguous` because its bounds are remapped below, and a rank
         !! remapping needs a target the compiler knows is not a stride.
      integer, intent(in) :: p_lo, p_hi
      real(dp), allocatable, target, intent(out) :: gamma_blk(:, :, :, :)

      real(dp), allocatable, target :: buf1(:), buf2(:)
      real(dp), pointer :: src(:, :), dst(:, :)
      integer :: n_ao, n_act, np, need

      n_ao = size(c_active, 1)
      n_act = size(c_active, 2)
      np = p_hi - p_lo + 1

      allocate (gamma_blk(np, n_ao, n_ao, n_ao))
      need = max(n_act**3*np, n_act*np*n_ao*n_ao, n_act**2*np*n_ao)
      allocate (buf1(need), buf2(need))

      ! (t u v w) -> (u v w p)
      src(1:n_act, 1:n_act**3) => ddm2
      dst(1:n_act**3, 1:np) => buf1
      call pic_gemm(src, c_active(p_lo:p_hi, :), dst, transa="T", transb="T", &
                    beta=0.0_dp)

      ! (u v w p) -> (v w p q)
      src(1:n_act, 1:n_act**2*np) => buf1
      dst(1:n_act**2*np, 1:n_ao) => buf2
      call pic_gemm(src, c_active, dst, transa="T", transb="T", beta=0.0_dp)

      ! (v w p q) -> (w p q r)
      src(1:n_act, 1:n_act*np*n_ao) => buf2
      dst(1:n_act*np*n_ao, 1:n_ao) => buf1
      call pic_gemm(src, c_active, dst, transa="T", transb="T", beta=0.0_dp)

      ! (w p q r) -> (p q r s), straight into the result
      src(1:n_act, 1:np*n_ao*n_ao) => buf1
      dst(1:np*n_ao*n_ao, 1:n_ao) => gamma_blk
      call pic_gemm(src, c_active, dst, transa="T", transb="T", beta=0.0_dp)

      nullify (src, dst)
      deallocate (buf1, buf2)
   end subroutine gamma_block

end module mqc_czt_mcscf_gradient
