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
   use mqc_libcint_hess_ints, only: eri_ip1_block
   implicit none
   private

   public :: hcore_deriv_atom
   public :: make_h1_atom
   public :: overlap_deriv_atom

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

end module mqc_libcint_hessian
