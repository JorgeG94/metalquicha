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
                                   DERIV_KIN, DERIV_NUC
   implicit none
   private

   public :: hcore_deriv_atom

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

end module mqc_libcint_hessian
