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
   use mqc_libcint_direct, only: build_fock_direct, schwarz_bounds, direct_stats_t
   use mqc_libcint_response, only: response_operator_t, solve_response
   use pic_blas_interfaces, only: pic_gemm
   implicit none
   private

   public :: hcore_deriv_atom
   public :: make_h1_atom
   public :: overlap_deriv_atom
   public :: nuclear_response_t
   public :: solve_mo1_atom
   public :: nuclear_repulsion_hessian

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

      n = this%n_mo*this%n_occ
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

      real(dp), allocatable :: mo1(:, :), half(:, :), dens(:, :), g(:, :), work(:, :)
      real(dp), allocatable :: gmo(:, :)
      type(direct_stats_t) :: stats
      integer :: n_ao, i, a

      if (error%has_error()) return

      n_ao = size(this%orbitals, 1)
      allocate (mo1(this%n_mo, this%n_occ))
      mo1 = reshape(vector, [this%n_mo, this%n_occ])

      ! The density this response implies, symmetrised.
      allocate (half(n_ao, this%n_occ), dens(n_ao, n_ao))
      ! **Twice, for double occupancy.** The trial vector describes how one
      ! spatial orbital moves; the density it perturbs holds two electrons in
      ! that orbital. Leaving the factor out halves the response, which does not
      ! stop the iteration converging -- it converges to the wrong answer, and
      ! the first-order density comes out roughly a third too large rather than
      ! obviously doubled, because the error feeds back through the coupling.
      call pic_gemm(this%orbitals, mo1, half, alpha=2.0_dp, beta=0.0_dp)
      call pic_gemm(half, this%c_occ, dens, transb="T")
      dens = dens + transpose(dens)

      allocate (g(n_ao, n_ao))
      call build_fock_direct(this%mol, this%zero_h, dens, this%bounds, g, stats, &
                             error, density_screen=.false.)
      if (error%has_error()) return

      ! Back to the molecular basis, occupied columns only.
      allocate (work(n_ao, this%n_occ), gmo(this%n_mo, this%n_occ))
      call pic_gemm(g, this%c_occ, work)
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

      image = reshape(gmo, [this%n_mo*this%n_occ])
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

end module mqc_libcint_hessian
