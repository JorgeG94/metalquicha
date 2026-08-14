!! Placing a fragment in an arbitrary orientation
module mqc_efp_rotate
   !! A potential carries the geometry it was made at. A deck puts its atoms
   !! somewhere else, in general rotated as well as translated, and *everything* the
   !! potential stores is expressed in its own frame: the multipoles, the
   !! polarizability tensors at every rank, and the localized orbitals. Placing the
   !! fragment means rotating all of it.
   !!
   !! **The orbitals are the part that is easy to forget.** Cartesian `d` and `f`
   !! basis functions rotate among themselves, so a coefficient vector expressed in a
   !! basis at one orientation is wrong in a basis at another. GAMESS does this too --
   !! it is why `FRGROT` takes `PROVEC` and `CTVEC` as arguments and works through
   !! them shell by shell (`efinp.src:2090` onward). Skipping it leaves exchange
   !! repulsion, charge transfer and every damped dispersion term quietly wrong, the
   !! damping being an overlap between localized orbitals on different fragments.
   !!
   !! **The rotation itself is three atoms, not Euler angles.** `ROTMAT`
   !! (`efinp.src:5681`) builds an orthonormal frame from two edges of the triangle
   !! three points make, and the rotation is the one carrying the potential's own
   !! frame onto the deck's. There are no angle conventions to pin and no chirality
   !! ambiguity: the third axis is a cross product.
   !!
   !! Our route through the orbitals is shorter than GAMESS's. It rotates them in its
   !! own AO convention, where each Cartesian component carries its own
   !! normalization, hence the `2/sqrt(3)` factors dividing in and out of `FRGROT`.
   !! We convert to libcint's convention first, where every component of a shell
   !! shares one normalization constant, so the rotation there is the plain monomial
   !! expansion with no normalization ratios in it at all.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_efp_read, only: efp_fragment_t
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim
   use libcint_fortran, only: LIBCINT_ANG_OF
   use mqc_efp_potential, only: from_gamess_ao_order, to_gamess_ao_order
   use mqc_efp_pair, only: fragment_molecule
   implicit none
   private

   public :: superpose
   public :: rotate_fragment
   public :: cartesian_rotation

   !> Highest angular momentum a projection basis is handled at. `mqc_efp_potential`
   !> refuses to write anything higher, so a shell beyond this is a corrupt file
   !> rather than a case to support.
   integer, parameter :: MAX_L = 4

contains

   subroutine superpose(own, deck, rot, translation, rmsd, error)
      !! The rigid transform carrying a potential's own atoms onto a deck's
      !!
      !! `deck(:,k) ~ rot . own(:,k) + translation`. Built the way `ROTMAT` builds it,
      !! from the frame three atoms define rather than from a least-squares fit: an
      !! EFP fragment is rigid, so the transform is exact and there is nothing to fit.
      !! `rmsd` is then a *test* of that rigidity rather than a residual to minimise,
      !! and the caller is expected to refuse a placement where it is not tiny --
      !! which is what catches a potential paired with the wrong species, atoms listed
      !! in a different order, or a geometry that has been relaxed since.
      real(dp), intent(in) :: own(:, :)    !! (3, n) the potential's own atoms
      real(dp), intent(in) :: deck(:, :)   !! (3, n) where the deck puts them
      real(dp), intent(out) :: rot(3, 3)
      real(dp), intent(out) :: translation(3)
      real(dp), intent(out) :: rmsd
      type(error_t), intent(inout) :: error

      real(dp) :: t_own(3, 3), t_deck(3, 3), own_centre(3), deck_centre(3)
      real(dp) :: gap(3)
      integer :: n, k
      logical :: ok

      rot = 0.0_dp
      translation = 0.0_dp
      rmsd = 0.0_dp
      do k = 1, 3
         rot(k, k) = 1.0_dp
      end do

      n = size(own, 2)
      if (size(deck, 2) /= n) then
         call error%set(ERROR_VALIDATION, "efp: the deck and the potential disagree "// &
                        "on how many atoms this fragment has")
         return
      end if

      ! One atom is a translation and nothing else -- there is no orientation to fix,
      ! and any rotation about it is as good as the identity. Two atoms leave a free
      ! rotation about their axis, which the fragment's own tensors are *not*
      ! invariant under, so it is refused rather than guessed.
      if (n == 1) then
         translation = deck(:, 1) - own(:, 1)
         return
      end if
      if (n == 2) then
         call error%set(ERROR_VALIDATION, "efp: a two-atom fragment leaves a free "// &
                        "rotation about its own axis, which its multipoles and "// &
                        "polarizabilities are not invariant under. Placing one needs "// &
                        "an orientation the deck does not carry.")
         return
      end if

      call orientation_frame(own(:, 1), own(:, 2), own(:, 3), t_own, ok)
      if (.not. ok) then
         call error%set(ERROR_VALIDATION, "efp: this fragment's first three atoms "// &
                        "are collinear, so they fix no orientation")
         return
      end if
      call orientation_frame(deck(:, 1), deck(:, 2), deck(:, 3), t_deck, ok)
      if (.not. ok) then
         call error%set(ERROR_VALIDATION, "efp: this fragment's first three atoms "// &
                        "are collinear in the deck, so they fix no orientation")
         return
      end if

      ! Both frames are orthonormal, so carrying one onto the other is one matrix
      ! product and the result is orthogonal by construction rather than by
      ! orthogonalising afterwards.
      rot = matmul(t_deck, transpose(t_own))

      own_centre = 0.0_dp
      deck_centre = 0.0_dp
      do k = 1, n
         own_centre = own_centre + own(:, k)
         deck_centre = deck_centre + deck(:, k)
      end do
      own_centre = own_centre/real(n, dp)
      deck_centre = deck_centre/real(n, dp)
      translation = deck_centre - matmul(rot, own_centre)

      do k = 1, n
         gap = deck(:, k) - (matmul(rot, own(:, k)) + translation)
         rmsd = rmsd + dot_product(gap, gap)
      end do
      rmsd = sqrt(rmsd/real(n, dp))
   end subroutine superpose

   pure subroutine orientation_frame(a, b, c, t, ok)
      !! An orthonormal frame from three points, as `ROTMAT` builds it
      !!
      !! First axis along `b - a`; second is `c - a` with that component projected
      !! out; third is their cross product, which fixes the handedness so the result
      !! is a rotation and never a reflection.
      real(dp), intent(in) :: a(3), b(3), c(3)
      real(dp), intent(out) :: t(3, 3)
      logical, intent(out) :: ok

      real(dp) :: e1(3), e2(3), e3(3), n1, n2, n3, overlap
      real(dp), parameter :: TINY_AXIS = 1.0e-10_dp

      ok = .false.
      t = 0.0_dp

      e1 = b - a
      n1 = sqrt(dot_product(e1, e1))
      if (n1 < TINY_AXIS) return
      e1 = e1/n1

      e2 = c - a
      overlap = dot_product(e2, e1)
      e2 = e2 - overlap*e1
      n2 = sqrt(dot_product(e2, e2))
      if (n2 < TINY_AXIS) return       ! c lies on the a-b line
      e2 = e2/n2

      e3(1) = e1(2)*e2(3) - e1(3)*e2(2)
      e3(2) = e1(3)*e2(1) - e1(1)*e2(3)
      e3(3) = e1(1)*e2(2) - e1(2)*e2(1)
      n3 = sqrt(dot_product(e3, e3))
      if (n3 < TINY_AXIS) return
      e3 = e3/n3

      t(:, 1) = e1
      t(:, 2) = e2
      t(:, 3) = e3
      ok = .true.
   end subroutine orientation_frame

   subroutine rotate_fragment(frag, rot, error)
      !! Turn a fragment, and everything it carries, into a new orientation
      !!
      !! Positions are rotated about the origin and *not* translated: the placement's
      !! shift stays an offset the energy routines add, which is the interface they
      !! already have. So a placed fragment is `rot . stored + offset`.
      !!
      !! What is invariant and deliberately untouched: the Fock matrix over localized
      !! orbitals, which is expressed in the orbitals themselves rather than in space;
      !! the basis exponents and contraction coefficients, which are radial; and the
      !! monopoles.
      type(efp_fragment_t), intent(inout) :: frag
      real(dp), intent(in) :: rot(3, 3)
      type(error_t), intent(inout) :: error

      integer :: k, f

      ! Positions and the one vector-valued moment.
      if (allocated(frag%points)) frag%points = matmul(rot, frag%points)
      if (allocated(frag%centroids)) frag%centroids = matmul(rot, frag%centroids)
      if (allocated(frag%pol_points)) frag%pol_points = matmul(rot, frag%pol_points)
      if (allocated(frag%dipole)) frag%dipole = matmul(rot, frag%dipole)

      if (allocated(frag%quadrupole)) then
         do k = 1, size(frag%quadrupole, 2)
            frag%quadrupole(:, k) = rotate_packed_quadrupole(frag%quadrupole(:, k), rot)
         end do
      end if
      if (allocated(frag%octopole)) then
         do k = 1, size(frag%octopole, 2)
            frag%octopole(:, k) = rotate_packed_octopole(frag%octopole(:, k), rot)
         end do
      end if

      if (allocated(frag%static_pol)) then
         do k = 1, size(frag%static_pol, 3)
            frag%static_pol(:, :, k) = rank2(frag%static_pol(:, :, k), rot)
         end do
      end if
      if (allocated(frag%dyn_pol)) then
         do f = 1, size(frag%dyn_pol, 4)
            do k = 1, size(frag%dyn_pol, 3)
               frag%dyn_pol(:, :, k, f) = rank2(frag%dyn_pol(:, :, k, f), rot)
            end do
         end do
      end if

      ! The two higher dispersion blocks are stored flat in the file's own slot
      ! order, so they are unpacked, rotated and repacked rather than rotated in
      ! place. Their index conventions are GAMESS's and are documented where they are
      ! consumed, in `mqc_efp_pair`.
      if (allocated(frag%dipquad)) then
         do f = 1, size(frag%dipquad, 3)
            do k = 1, size(frag%dipquad, 2)
               frag%dipquad(:, k, f) = rotate_flat_rank3(frag%dipquad(:, k, f), rot)
            end do
         end do
      end if
      if (allocated(frag%quadquad)) then
         do f = 1, size(frag%quadquad, 3)
            do k = 1, size(frag%quadquad, 2)
               frag%quadquad(:, k, f) = rotate_flat_rank4(frag%quadquad(:, k, f), rot)
            end do
         end do
      end if

      call rotate_orbitals(frag, rot, error)
   end subroutine rotate_fragment

   subroutine rotate_orbitals(frag, rot, error)
      !! The localized and charge-transfer orbitals, through their AO coefficients
      !!
      !! Cartesian `d` and `f` functions mix under rotation, so a coefficient vector
      !! is only meaningful together with the orientation of the basis it is expressed
      !! in. Both stored sets are in GAMESS's AO order and normalization, so each is
      !! converted into ours, rotated shell by shell, and converted back -- leaving
      !! the fragment storing what it always stored, in a new orientation.
      type(efp_fragment_t), intent(inout) :: frag
      real(dp), intent(in) :: rot(3, 3)
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      real(dp), allocatable :: ours(:, :), spun(:, :), back(:, :)

      if (.not. (frag%has_lmo .or. frag%has_ctvec)) return
      if (.not. frag%has_basis) then
         call error%set(ERROR_VALIDATION, "efp: rotating this fragment's orbitals "// &
                        "needs its projection basis, which it does not carry")
         return
      end if

      call fragment_molecule(frag, [0.0_dp, 0.0_dp, 0.0_dp], mol, error)
      if (error%has_error()) return

      if (frag%has_lmo) then
         call from_gamess_ao_order(mol, frag%lmo_gamess, ours, error)
         if (error%has_error()) return
         call spin_coefficients(mol, ours, rot, spun, error)
         if (error%has_error()) return
         call to_gamess_ao_order(mol, spun, back, error)
         if (error%has_error()) return
         frag%lmo_gamess = back
         deallocate (ours, spun, back)
      end if

      if (frag%has_ctvec) then
         call from_gamess_ao_order(mol, frag%ctvec_gamess, ours, error)
         if (error%has_error()) return
         call spin_coefficients(mol, ours, rot, spun, error)
         if (error%has_error()) return
         call to_gamess_ao_order(mol, spun, back, error)
         if (error%has_error()) return
         frag%ctvec_gamess = back
         deallocate (ours, spun, back)
      end if

      call mol%destroy()
   end subroutine rotate_orbitals

   subroutine spin_coefficients(mol, coefficients, rot, spun, error)
      !! Apply the AO rotation to a set of coefficient columns
      !!
      !! Block diagonal over shells: a rotation preserves the radial part, so a shell
      !! maps onto itself and never onto another with the same angular momentum.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coefficients(:, :)
      real(dp), intent(in) :: rot(3, 3)
      real(dp), allocatable, intent(out) :: spun(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: d(:, :)
      integer :: ish, l, off, dim

      allocate (spun(size(coefficients, 1), size(coefficients, 2)))
      spun = coefficients
      do ish = 1, mol%nbas
         l = mol%bas(LIBCINT_ANG_OF, ish)
         if (l == 0) cycle                    ! an s function is invariant
         if (l > MAX_L) then
            call error%set(ERROR_VALIDATION, "efp: a projection basis shell above "// &
                           "g cannot be rotated by this code")
            return
         end if
         off = mol%shell_offset(ish)
         dim = shell_dim(mol%cartesian, ish - 1, mol%bas)
         call cartesian_rotation(l, rot, d)
         if (size(d, 1) /= dim) then
            call error%set(ERROR_VALIDATION, "efp: a shell's size does not match "// &
                           "its angular momentum, so the basis is spherical where "// &
                           "this code assumes Cartesian")
            return
         end if
         spun(off + 1:off + dim, :) = matmul(d, coefficients(off + 1:off + dim, :))
         deallocate (d)
      end do
   end subroutine spin_coefficients

   pure function rank2(a, rot) result(b)
      !! `B = R A R^T`
      real(dp), intent(in) :: a(3, 3), rot(3, 3)
      real(dp) :: b(3, 3)

      b = matmul(rot, matmul(a, transpose(rot)))
   end function rank2

   pure function rotate_packed_quadrupole(packed, rot) result(out)
      !! Six stored components, ordered `xx yy zz xy xz yz`
      real(dp), intent(in) :: packed(6), rot(3, 3)
      real(dp) :: out(6)

      integer, parameter :: QI(6) = [1, 2, 3, 1, 1, 2]
      integer, parameter :: QJ(6) = [1, 2, 3, 2, 3, 3]
      real(dp) :: full(3, 3)
      integer :: s

      full = 0.0_dp
      do s = 1, 6
         full(QI(s), QJ(s)) = packed(s)
         full(QJ(s), QI(s)) = packed(s)
      end do
      full = rank2(full, rot)
      do s = 1, 6
         out(s) = full(QI(s), QJ(s))
      end do
   end function rotate_packed_quadrupole

   pure function rotate_packed_octopole(packed, rot) result(out)
      !! Ten stored components, `xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz`
      !!
      !! GAMESS's order, read off how `efelec.src` unpacks `EFOCT` -- the same order
      !! `mqc_efp_interaction` spreads into a full tensor, and repeated here rather
      !! than shared because a rotation needs the full tensor anyway.
      real(dp), intent(in) :: packed(10), rot(3, 3)
      real(dp) :: out(10)

      integer, parameter :: OI(10) = [1, 2, 3, 1, 1, 1, 2, 1, 2, 1]
      integer, parameter :: OJ(10) = [1, 2, 3, 1, 1, 2, 2, 3, 3, 2]
      integer, parameter :: OK(10) = [1, 2, 3, 2, 3, 2, 3, 3, 3, 3]
      real(dp) :: full(3, 3, 3), spun(3, 3, 3)
      integer :: s, i, j, k

      full = 0.0_dp
      do s = 1, 10
         call spread_symmetric(full, OI(s), OJ(s), OK(s), packed(s))
      end do
      spun = 0.0_dp
      do k = 1, 3
         do j = 1, 3
            do i = 1, 3
               spun(i, j, k) = contract3(full, rot, i, j, k)
            end do
         end do
      end do
      do s = 1, 10
         out(s) = spun(OI(s), OJ(s), OK(s))
      end do
   end function rotate_packed_octopole

   pure subroutine spread_symmetric(full, i, j, k, value)
      !! One stored component into every index permutation it stands for
      real(dp), intent(inout) :: full(3, 3, 3)
      integer, intent(in) :: i, j, k
      real(dp), intent(in) :: value

      full(i, j, k) = value
      full(i, k, j) = value
      full(j, i, k) = value
      full(j, k, i) = value
      full(k, i, j) = value
      full(k, j, i) = value
   end subroutine spread_symmetric

   pure function contract3(t, rot, a, b, c) result(v)
      real(dp), intent(in) :: t(3, 3, 3), rot(3, 3)
      integer, intent(in) :: a, b, c
      real(dp) :: v
      integer :: i, j, k

      v = 0.0_dp
      do k = 1, 3
         do j = 1, 3
            do i = 1, 3
               v = v + rot(a, i)*rot(b, j)*rot(c, k)*t(i, j, k)
            end do
         end do
      end do
   end function contract3

   pure function rotate_flat_rank3(flat, rot) result(out)
      !! 27 slots, `slot = (i-1)*9 + (j-1)*3 + k`
      real(dp), intent(in) :: flat(27), rot(3, 3)
      real(dp) :: out(27)

      real(dp) :: full(3, 3, 3)
      integer :: i, j, k

      do k = 1, 3
         do j = 1, 3
            do i = 1, 3
               full(i, j, k) = flat((i - 1)*9 + (j - 1)*3 + k)
            end do
         end do
      end do
      do k = 1, 3
         do j = 1, 3
            do i = 1, 3
               out((i - 1)*9 + (j - 1)*3 + k) = contract3(full, rot, i, j, k)
            end do
         end do
      end do
   end function rotate_flat_rank3

   pure function rotate_flat_rank4(flat, rot) result(out)
      !! 81 slots, `slot = (i-1)*27 + (j-1)*9 + (k-1)*3 + l`
      real(dp), intent(in) :: flat(81), rot(3, 3)
      real(dp) :: out(81)

      real(dp) :: full(3, 3, 3, 3), step1(3, 3, 3, 3)
      integer :: a, b, c, e, i, j, k, l

      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               do i = 1, 3
                  full(i, j, k, l) = flat((i - 1)*27 + (j - 1)*9 + (k - 1)*3 + l)
               end do
            end do
         end do
      end do

      ! One index at a time: four passes of 3^5 rather than one of 3^8.
      step1 = 0.0_dp
      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               do a = 1, 3
                  do i = 1, 3
                     step1(a, j, k, l) = step1(a, j, k, l) + rot(a, i)*full(i, j, k, l)
                  end do
               end do
            end do
         end do
      end do
      full = 0.0_dp
      do l = 1, 3
         do k = 1, 3
            do b = 1, 3
               do j = 1, 3
                  do a = 1, 3
                     full(a, b, k, l) = full(a, b, k, l) + rot(b, j)*step1(a, j, k, l)
                  end do
               end do
            end do
         end do
      end do
      step1 = 0.0_dp
      do l = 1, 3
         do c = 1, 3
            do k = 1, 3
               do b = 1, 3
                  do a = 1, 3
                     step1(a, b, c, l) = step1(a, b, c, l) + rot(c, k)*full(a, b, k, l)
                  end do
               end do
            end do
         end do
      end do
      full = 0.0_dp
      do e = 1, 3
         do l = 1, 3
            do c = 1, 3
               do b = 1, 3
                  do a = 1, 3
                     full(a, b, c, e) = full(a, b, c, e) + rot(e, l)*step1(a, b, c, l)
                  end do
               end do
            end do
         end do
      end do

      do l = 1, 3
         do k = 1, 3
            do j = 1, 3
               do i = 1, 3
                  out((i - 1)*27 + (j - 1)*9 + (k - 1)*3 + l) = full(i, j, k, l)
               end do
            end do
         end do
      end do
   end function rotate_flat_rank4

   subroutine cartesian_rotation(l, rot, d)
      !! How a Cartesian shell's components mix under a rotation
      !!
      !! For a component `x^a y^b z^c`, rotating the function means evaluating it at
      !! `R^T r`, and each factor becomes a linear form in `x, y, z` built from a
      !! column of `R`. Multiplying the three expansions out and collecting monomials
      !! gives the matrix.
      !!
      !! **There are no normalization ratios here, and that is a property of the
      !! convention rather than an omission.** libcint gives every Cartesian component
      !! of a shell the same normalization constant -- the one belonging to
      !! `(l,0,0)` -- so the ratio between any two is 1 and the monomial expansion *is*
      !! the transformation. In a convention that normalizes each component
      !! separately, as GAMESS's does, factors of `sqrt(3)` and friends appear; that
      !! is what `from_gamess_ao_order` has already divided out by the time this runs.
      integer, intent(in) :: l
      real(dp), intent(in) :: rot(3, 3)
      real(dp), allocatable, intent(out) :: d(:, :)

      integer :: dim, p, q
      integer :: ea(3), sub(3, 3), tot(3)
      integer :: i1, j1, i2, j2, i3, j3, t
      real(dp) :: coef
      integer, allocatable :: expo(:, :)

      dim = (l + 1)*(l + 2)/2
      allocate (d(dim, dim))
      d = 0.0_dp
      call cartesian_exponents(l, expo)

      do p = 1, dim
         ea = expo(:, p)
         ! Distribute each factor's exponent over x, y and z. Three nested
         ! multinomials, one per Cartesian axis of the *rotated* argument.
         do i1 = 0, ea(1)
            do j1 = 0, ea(1) - i1
               sub(:, 1) = [i1, j1, ea(1) - i1 - j1]
               do i2 = 0, ea(2)
                  do j2 = 0, ea(2) - i2
                     sub(:, 2) = [i2, j2, ea(2) - i2 - j2]
                     do i3 = 0, ea(3)
                        do j3 = 0, ea(3) - i3
                           sub(:, 3) = [i3, j3, ea(3) - i3 - j3]

                           coef = 1.0_dp
                           do t = 1, 3
                              coef = coef*multinomial(sub(:, t)) &
                                     *rot(1, t)**sub(1, t) &
                                     *rot(2, t)**sub(2, t) &
                                     *rot(3, t)**sub(3, t)
                           end do
                           tot = sub(:, 1) + sub(:, 2) + sub(:, 3)
                           q = exponent_slot(l, tot, expo)
                           d(q, p) = d(q, p) + coef
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      deallocate (expo)
   end subroutine cartesian_rotation

   pure subroutine cartesian_exponents(l, expo)
      !! libcint's Cartesian component order: `x` descending, then `y`
      integer, intent(in) :: l
      integer, allocatable, intent(out) :: expo(:, :)

      integer :: dim, i, j, n

      dim = (l + 1)*(l + 2)/2
      allocate (expo(3, dim))
      n = 0
      do i = l, 0, -1
         do j = l - i, 0, -1
            n = n + 1
            expo(:, n) = [i, j, l - i - j]
         end do
      end do
   end subroutine cartesian_exponents

   pure function exponent_slot(l, e, expo) result(slot)
      !! Which component carries these exponents
      integer, intent(in) :: l, e(3)
      integer, intent(in) :: expo(:, :)
      integer :: slot

      integer :: n

      slot = 0
      do n = 1, size(expo, 2)
         if (all(expo(:, n) == e)) then
            slot = n
            return
         end if
      end do
   end function exponent_slot

   pure function multinomial(e) result(m)
      !! `(e1+e2+e3)! / (e1! e2! e3!)`
      integer, intent(in) :: e(3)
      real(dp) :: m

      m = factorial(e(1) + e(2) + e(3))/(factorial(e(1))*factorial(e(2))*factorial(e(3)))
   end function multinomial

   pure function factorial(n) result(f)
      integer, intent(in) :: n
      real(dp) :: f
      integer :: k

      f = 1.0_dp
      do k = 2, n
         f = f*real(k, dp)
      end do
   end function factorial

end module mqc_efp_rotate
