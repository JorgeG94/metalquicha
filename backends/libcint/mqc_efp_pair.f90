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
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cgto, only: molecular_basis_type, ANGULAR_FORM_CARTESIAN
   use mqc_elements, only: element_number_to_symbol
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_efp_read, only: efp_fragment_t
   implicit none
   private

   public :: fragment_basis
   public :: two_fragment_molecule

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
