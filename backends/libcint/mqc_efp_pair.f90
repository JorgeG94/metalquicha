!! A molecule spanning two fragments, and the overlaps between them
module mqc_efp_pair
   !! Exchange repulsion, charge transfer and the damping on dispersion all come from
   !! the same thing: overlaps between *one fragment's* localized orbitals and
   !! *another's*. Those are integrals over both fragments' basis functions at once,
   !! so they need a molecule that spans the pair -- which is what this builds.
   !!
   !! **Nothing new is needed in the integral layer.** `libcint_molecule_t` already
   !! has `build`, which takes a `molecular_basis_type` directly;
   !! `build_libcint_molecule` is only the wrapper that looks a basis up by name. So
   !! the work here is assembling the basis object, not plumbing integrals.
   !!
   !! The basis comes from the potential's own `PROJECTION BASIS SET`, recovered by
   !! `mqc_efp_read` with GAMESS's primitive normalization divided back out. Using the
   !! *file's* basis rather than looking up the name it was computed with matters: a
   !! potential is a self-contained object, and the shipped GAMESS library potentials
   !! do not all name a basis this program has.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cgto, only: molecular_basis_type, ANGULAR_FORM_CARTESIAN
   use mqc_elements, only: element_number_to_symbol
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_efp_read, only: efp_fragment_t
   use mqc_efp_potential, only: from_gamess_ao_order, frozen_core
   implicit none
   private

   public :: fragment_basis
   public :: two_fragment_molecule
   public :: fragment_molecule
   public :: fragment_lmo
   public :: exchange_repulsion

contains

   subroutine fragment_basis(frag, basis, error)
      !! A basis object for one fragment, from the projection basis it carries
      !!
      !! Cartesian, because that is what a potential this program writes is: the
      !! ordering maps in `mqc_efp_potential` cover Cartesian s, p, d and f, and a
      !! spherical potential is refused at write time. A GAMESS library potential
      !! read here would need that assumption revisited.
      type(efp_fragment_t), intent(in) :: frag
      type(molecular_basis_type), intent(out) :: basis
      type(error_t), intent(inout) :: error

      integer :: natm, at, sh, k, count, first

      if (.not. frag%has_basis) then
         call error%set(ERROR_VALIDATION, "efp: this fragment carries no projection "// &
                        "basis, so no molecule can be built from it")
         return
      end if

      natm = frag%n_atoms
      call basis%allocate_elements(natm)
      do at = 1, natm
         ! How many shells this atom owns. The reader tags each shell with the atom
         ! it was read under, in file order, so a count is enough.
         count = 0
         do sh = 1, frag%n_shells
            if (frag%shell_atom(sh) == at) count = count + 1
         end do
         if (count == 0) then
            call error%set(ERROR_VALIDATION, "efp: an atom of the projection basis "// &
                           "carries no shells")
            return
         end if
         basis%elements(at)%element = trim(element_number_to_symbol(nint(frag%charge(at))))
         basis%elements(at)%angular_form = ANGULAR_FORM_CARTESIAN
         call basis%elements(at)%allocate_shells(count)
         count = 0
         do sh = 1, frag%n_shells
            if (frag%shell_atom(sh) /= at) cycle
            count = count + 1
            first = frag%shell_first(sh)
            call basis%elements(at)%shells(count)%allocate_arrays(frag%shell_nprim(sh))
            basis%elements(at)%shells(count)%ang_mom = frag%shell_l(sh)
            basis%elements(at)%shells(count)%nfunc = frag%shell_nprim(sh)
            do k = 1, frag%shell_nprim(sh)
               basis%elements(at)%shells(count)%exponents(k) = frag%prim_expo(first + k - 1)
               basis%elements(at)%shells(count)%coefficients(k) = frag%prim_coef(first + k - 1)
            end do
         end do
      end do
   end subroutine fragment_basis

   function exchange_repulsion(frag_a, frag_b, offset_a, offset_b, error) result(energy)
      !! Pauli exchange repulsion between two fragments
      !!
      !! From `EXREP` in GAMESS's `efpaul.src`, which is the energy routine -- `ODM` in
      !! the same file is its gradient and carries factors that only mean something in a
      !! derivative. Three terms over pairs of localized orbitals on different
      !! fragments, with `S` and `T` in the localized-orbital basis:
      !!
      !!     XR1 = -2 sqrt(2/pi) sum_ij sqrt(-ln|S_ij|) S_ij^2 / R_ij
      !!     XR2 = -sum_ij S_ij ( sum_k F^A_ik S_kj + sum_l F^B_jl S_il - 2 T_ij )
      !!     XR3 = +sum_ij S_ij^2 ( V_i + V_j - 1/R_ij )
      !!     E   = 2 (XR1 + XR2 + XR3)
      !!
      !! `V_i` is the electrostatic potential at orbital `i`'s centroid from the *other*
      !! fragment: minus its nuclear charges over their distances, plus two over the
      !! distance to each of its orbital centroids. The factor of two outside, and the
      !! two inside `V`, are `ICOEFF` and `JCOEFF` at `MLSWTCH = 1`, which is the
      !! closed-shell case.
      !!
      !! The first term is skipped where `|S_ij| <= 1e-7` -- its logarithm diverges as
      !! the overlap vanishes, and the `S^2` in front does not save it in floating
      !! point. GAMESS skips a whole fragment pair when every `|S_ij| <= 1e-6`; that is
      !! a screening decision rather than a numerical one and is not copied here.
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      type(error_t), intent(inout) :: error
      real(dp) :: energy

      real(dp), parameter :: RT2PI = 0.7978845608028654_dp   !! sqrt(2/pi)
      real(dp), parameter :: S_FLOOR = 1.0e-7_dp
      type(libcint_molecule_t) :: pair
      real(dp), allocatable :: s_ao(:, :), t_ao(:, :), lmo_a(:, :), lmo_b(:, :)
      real(dp), allocatable :: s_lmo(:, :), t_lmo(:, :), work(:, :)
      real(dp), allocatable :: cen_a(:, :), cen_b(:, :), v_a(:), v_b(:)
      integer :: n_ao_a, na, nb, i, j, n_lmo_a, n_lmo_b
      real(dp) :: xr1, xr2, xr3, sij, sij2, tij, rij, fij, fji
      integer :: k, l

      energy = 0.0_dp
      if (.not. (frag_a%has_fock .and. frag_b%has_fock)) then
         call error%set(ERROR_VALIDATION, "efp: exchange repulsion needs the LMO Fock "// &
                        "matrix of both fragments")
         return
      end if

      call two_fragment_molecule(frag_a, frag_b, offset_a, offset_b, pair, n_ao_a, error)
      if (error%has_error()) return
      call pair%overlap(s_ao)
      call pair%kinetic(t_ao)

      ! Each fragment's orbitals live on its own functions, so they are padded into the
      ! pair's space rather than transformed with a rectangular block: that way the
      ! transform is one gemm and the block boundary appears once, here.
      n_lmo_a = frag_a%n_lmo_proj
      n_lmo_b = frag_b%n_lmo_proj
      call padded_lmo(frag_a, pair, 0, n_ao_a, lmo_a, error)
      if (error%has_error()) return
      call padded_lmo(frag_b, pair, n_ao_a, pair%nao - n_ao_a, lmo_b, error)
      if (error%has_error()) return

      allocate (work(pair%nao, n_lmo_b), s_lmo(n_lmo_a, n_lmo_b), t_lmo(n_lmo_a, n_lmo_b))
      call pic_gemm(s_ao, lmo_b, work)
      call pic_gemm(lmo_a, work, s_lmo, transa="T")
      call pic_gemm(t_ao, lmo_b, work)
      call pic_gemm(lmo_a, work, t_lmo, transa="T")

      na = frag_a%n_atoms
      nb = frag_b%n_atoms
      allocate (cen_a(3, n_lmo_a), cen_b(3, n_lmo_b))
      do i = 1, n_lmo_a
         cen_a(:, i) = frag_a%pol_points(:, i) + offset_a
      end do
      do j = 1, n_lmo_b
         cen_b(:, j) = frag_b%pol_points(:, j) + offset_b
      end do

      ! The potential at each centroid from the other fragment, built once rather than
      ! inside the pair loop.
      allocate (v_a(n_lmo_a), v_b(n_lmo_b))
      do i = 1, n_lmo_a
         v_a(i) = other_potential(cen_a(:, i), frag_b, offset_b, cen_b)
      end do
      do j = 1, n_lmo_b
         v_b(j) = other_potential(cen_b(:, j), frag_a, offset_a, cen_a)
      end do

      xr1 = 0.0_dp
      xr2 = 0.0_dp
      xr3 = 0.0_dp
      do j = 1, n_lmo_b
         do i = 1, n_lmo_a
            sij = s_lmo(i, j)
            sij2 = sij*sij
            tij = t_lmo(i, j)
            rij = sqrt(sum((cen_a(:, i) - cen_b(:, j))**2))

            if (abs(sij) > S_FLOOR) then
               xr1 = xr1 - 2.0_dp*RT2PI*sqrt(-log(abs(sij)))*sij2/rij
            end if

            fij = 0.0_dp
            do k = 1, n_lmo_a
               fij = fij + frag_a%fock_lmo(i, k)*s_lmo(k, j)
            end do
            fji = 0.0_dp
            do l = 1, n_lmo_b
               fji = fji + frag_b%fock_lmo(j, l)*s_lmo(i, l)
            end do
            xr2 = xr2 - sij*(fij + fji - 2.0_dp*tij)

            xr3 = xr3 + sij2*(v_a(i) + v_b(j) - 1.0_dp/rij)
         end do
      end do
      energy = 2.0_dp*(xr1 + xr2 + xr3)

      call pair%destroy()
      deallocate (s_ao, t_ao, lmo_a, lmo_b, s_lmo, t_lmo, work, cen_a, cen_b, v_a, v_b)
   end function exchange_repulsion

   function other_potential(point, frag, offset, centroids) result(v)
      !! The potential at `point` from `frag`: its nuclei, then its orbital centroids
      real(dp), intent(in) :: point(3)
      type(efp_fragment_t), intent(in) :: frag
      real(dp), intent(in) :: offset(3)
      real(dp), intent(in) :: centroids(:, :)
      real(dp) :: v

      integer :: k

      v = 0.0_dp
      do k = 1, frag%n_atoms
         ! Valence charge, not the full nuclear one: the localized orbitals a potential
         ! carries are valence only, so the core electrons have to be counted as
         ! screening their own nucleus. This is the same number the PROJECTION BASIS SET
         ! atom header carries.
         v = v - (frag%charge(k) - 2.0_dp*real(frozen_core([nint(frag%charge(k))]), dp)) &
             /sqrt(sum((point - frag%points(:, k) - offset)**2))
      end do
      do k = 1, size(centroids, 2)
         v = v + 2.0_dp/sqrt(sum((point - centroids(:, k))**2))
      end do
   end function other_potential

   subroutine padded_lmo(frag, pair, offset_ao, n_own, padded, error)
      !! One fragment's orbitals in the pair's full AO space, zero elsewhere
      type(efp_fragment_t), intent(in) :: frag
      type(libcint_molecule_t), intent(in) :: pair
      integer, intent(in) :: offset_ao, n_own
      real(dp), allocatable, intent(out) :: padded(:, :)
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: own
      real(dp), allocatable :: lmo(:, :)

      call fragment_molecule(frag, [0.0_dp, 0.0_dp, 0.0_dp], own, error)
      if (error%has_error()) return
      if (own%nao /= n_own) then
         call error%set(ERROR_VALIDATION, "efp: a fragment's own basis does not match "// &
                        "its block of the pair")
         call own%destroy()
         return
      end if
      call fragment_lmo(frag, own, lmo, error)
      call own%destroy()
      if (error%has_error()) return
      allocate (padded(pair%nao, frag%n_lmo_proj))
      padded = 0.0_dp
      padded(offset_ao + 1:offset_ao + n_own, :) = lmo
      deallocate (lmo)
   end subroutine padded_lmo

   subroutine fragment_molecule(frag, offset, mol, error)
      !! One fragment as a molecule of its own, from the basis it carries
      !!
      !! Useful in its own right and as the reference the pair is checked against: the
      !! pair's diagonal block has to be this molecule's overlap.
      type(efp_fragment_t), intent(in) :: frag
      real(dp), intent(in) :: offset(3)
      type(libcint_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: error

      type(molecular_basis_type) :: basis
      real(dp), allocatable :: coords(:, :)
      integer, allocatable :: z(:)
      integer :: at

      call fragment_basis(frag, basis, error)
      if (error%has_error()) return
      allocate (coords(3, frag%n_atoms), z(frag%n_atoms))
      do at = 1, frag%n_atoms
         coords(:, at) = frag%points(:, at) + offset
         z(at) = nint(frag%charge(at))
      end do
      call mol%build(z, coords, basis, error, force_cartesian=.true.)
      deallocate (coords, z)
   end subroutine fragment_molecule

   subroutine fragment_lmo(frag, mol, lmo, error)
      !! The fragment's localized orbitals in *our* AO order
      !!
      !! The file stores them in GAMESS's, so the d and f permutation and its
      !! normalizations have to be undone -- `from_gamess_ao_order`, which lives beside
      !! the forward map so the two cannot drift apart.
      !!
      !! The check that this is right is that the orbitals come back orthonormal
      !! against the fragment's own overlap: `C^T S C = I`. Nothing weaker would
      !! notice a permutation applied in the wrong direction, because a wrongly
      !! permuted set is still a perfectly good-looking matrix.
      type(efp_fragment_t), intent(in) :: frag
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), allocatable, intent(out) :: lmo(:, :)
      type(error_t), intent(inout) :: error

      if (.not. frag%has_lmo) then
         call error%set(ERROR_VALIDATION, "efp: this fragment carries no projection "// &
                        "wavefunction")
         return
      end if
      if (frag%nao_proj /= mol%nao) then
         call error%set(ERROR_VALIDATION, "efp: the projection wavefunction and the "// &
                        "molecule disagree on the number of basis functions")
         return
      end if
      call from_gamess_ao_order(mol, frag%lmo_gamess, lmo, error)
   end subroutine fragment_lmo

   subroutine two_fragment_molecule(frag_a, frag_b, offset_a, offset_b, mol, &
                                    n_ao_a, error)
      !! One molecule covering both fragments, placed
      !!
      !! The atoms of `a` come first and then those of `b`, so the overlap matrix
      !! this yields is block structured: the leading `n_ao_a` rows and columns are
      !! `a`'s own overlap, the trailing block is `b`'s, and the off-diagonal block is
      !! what the inter-fragment terms want. `n_ao_a` is returned for exactly that
      !! reason -- without it the caller cannot find the block it came for.
      !!
      !! Only the atoms are included. A fragment's expansion points also sit at bond
      !! midpoints, and those carry multipoles but no basis functions, which is why
      !! the count here is `n_atoms` and not `n_points`.
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      type(libcint_molecule_t), intent(out) :: mol
      integer, intent(out) :: n_ao_a
      type(error_t), intent(inout) :: error

      type(molecular_basis_type) :: basis_a, basis_b, both
      real(dp), allocatable :: coords(:, :)
      integer, allocatable :: z(:)
      integer :: na, nb, at

      n_ao_a = 0
      call fragment_basis(frag_a, basis_a, error)
      if (error%has_error()) return
      call fragment_basis(frag_b, basis_b, error)
      if (error%has_error()) return

      na = frag_a%n_atoms
      nb = frag_b%n_atoms
      call both%allocate_elements(na + nb)
      do at = 1, na
         both%elements(at) = basis_a%elements(at)
      end do
      do at = 1, nb
         both%elements(na + at) = basis_b%elements(at)
      end do

      allocate (coords(3, na + nb), z(na + nb))
      do at = 1, na
         coords(:, at) = frag_a%points(:, at) + offset_a
         z(at) = nint(frag_a%charge(at))
      end do
      do at = 1, nb
         coords(:, na + at) = frag_b%points(:, at) + offset_b
         z(na + at) = nint(frag_b%charge(at))
      end do

      ! Cartesian explicitly, not left to the basis object to declare. A potential
      ! this program writes is Cartesian by construction -- the ordering maps cover
      ! Cartesian s through f and a spherical one is refused at write time -- and
      ! setting the angular form per element is not enough: the molecule asks the
      ! basis as a whole, which came back spherical and gave five d functions where
      ! the potential has six.
      call mol%build(z, coords, both, error, force_cartesian=.true.)
      if (error%has_error()) then
         deallocate (coords, z)
         return
      end if

      ! Where a's functions stop, counted from the basis rather than assumed to be
      ! half: the two fragments need not be the same species.
      n_ao_a = 0
      do at = 1, na
         n_ao_a = n_ao_a + basis_a%elements(at)%num_basis_functions()
      end do

      deallocate (coords, z)
   end subroutine two_fragment_molecule

end module mqc_efp_pair
