!! What the quasi-atomic orbitals turn out to be, and the table that says so
module mqc_libcint_quao_report
   !! The construction in `mqc_libcint_quao` produces a set of orbitals, one per
   !! atom, and a matrix of populations and bond orders. It does not say which
   !! orbital is a sigma bond and which is a lone pair, and without that the
   !! output is a numbered matrix that has to be read by hand against a picture
   !! of the molecule. This module assigns those labels and prints the report
   !! GAMESS prints -- every orbital pair with a kinetic bond order above a
   !! threshold, sorted by magnitude, each end annotated with its atom, its
   !! occupation, the atom it is bonded to, and what kind of orbital it is.
   !!
   !! **Every label here is empirical.** Nothing in Papers I or II defines a
   !! sigma orbital; the classification is a set of thresholds on occupations,
   !! bond orders and directions, chosen by the reference implementation to
   !! reproduce what a chemist would say. They are reproduced rather than
   !! reinvented -- `LOCAL_ORIENBOS_ASSIGN_ORBTYP` in GAMESS's `locsvd.src`, with
   !! its `BOTOL` array -- so that a disagreement with GAMESS is a disagreement
   !! about the wave function rather than about where a threshold was put. The
   !! numbers themselves carry no theory and should not be treated as if they
   !! do: an orbital half a percent the wrong side of `TOL_LONE_PAIR_OCC` is not
   !! a different orbital, it is a different label on the same orbital.
   !!
   !! The quantities being labelled -- populations, bond orders, kinetic bond
   !! orders -- are not empirical at all, and are validated against GAMESS to
   !! 1e-6 in `test/test_mqc_quao.f90`. It is worth keeping the two apart.
   use pic_types, only: dp
   use pic_io, only: to_char
   use pic_logger, only: logger => global_logger
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_libcint_multipole, only: multipole_matrices
   use mqc_libcint_quao, only: quao_result_t, mo_aambs_overlap
   use libcint_fortran, only: LIBCINT_ANG_OF, LIBCINT_ATOM_OF
   implicit none
   private

   public :: quao_labels_t
   public :: label_quasi_atomic_orbitals
   public :: print_quao_report
   public :: quao_type_name
   public :: QUAO_TYPE_LPMOD, QUAO_TYPE_NONE, QUAO_TYPE_LONE_PAIR, QUAO_TYPE_RADICAL
   public :: QUAO_TYPE_SIGMA, QUAO_TYPE_PI, QUAO_TYPE_RDLP, QUAO_TYPE_RDNV
   public :: QUAO_TYPE_NVMOD, QUAO_TYPE_NV
   public :: QUAO_MAX_L

   ! GAMESS's `BOTOL`, in its own order. The comments say what each one decides,
   ! which the names do not: `TOLNV` is not a tolerance on anything, it is the
   ! occupation below which an orbital counts as empty.
   real(dp), parameter :: TOL_STRONG_BOND = 0.40_dp
      !! `TOLSBO`. A bond order between two atoms at or above this makes the
      !! pair a "strong bond", and only strong bonds get orbitals labelled from
      !! them. Everything else in the table is reported without being classified.
   real(dp), parameter :: TOL_LONE_PAIR_OCC = 1.80_dp
      !! `TOLLPOCC`. An orbital this close to doubly occupied is a lone pair,
      !! whatever else it is doing.
   real(dp), parameter :: TOL_RADICAL_OCC = 0.40_dp
      !! `TOLRDOCC`. How far from singly occupied an orbital may be and still be
      !! called a radical.
   real(dp), parameter :: TOL_OVERLAP = 0.20_dp
      !! `TOLOVRLP`. The cosine between an orbital's own axis and the line
      !! joining the two atoms, above which the orbital points along the bond and
      !! is therefore sigma rather than pi. This is a geometric test, not a
      !! symmetry one -- see `assign_orbital_types`.
   real(dp), parameter :: TOL_EMPTY = 0.40_dp
      !! `TOLNV`. Occupation below which an orbital is empty for labelling
      !! purposes.
   real(dp), parameter :: TOL_PRINT = 1.0e-3_dp
      !! `TOLDEN`. Rows below this in magnitude are not printed. A formyl
      !! chloride table stops after 27 of its 78 orbital pairs.

   ! Column widths for the report, so the header and the rows cannot drift
   ! apart. Every field is written into a substring of exactly this length.
   integer, parameter :: WIDTH_NUMBER = 12     !! A bond order or an occupation
   integer, parameter :: WIDTH_ORB = 6         !! An orbital index
   integer, parameter :: WIDTH_ATOM = 7        !! A symbol and an atom index
   integer, parameter :: WIDTH_PARTNERS = 12   !! The bonded-to list, plus a lead space
   integer, parameter :: WIDTH_TYPE = 8        !! The orbital type, plus a lead space

   integer, parameter :: QUAO_MAX_L = 3
      !! s, p, d, f -- GAMESS's `INUMTYP`. The accurate atomic minimal basis has
      !! nothing higher, since it holds free-atom occupied shells only.

   ! Orbital types, numbered as GAMESS numbers them so the two can be compared
   ! without a translation table.
   integer, parameter :: QUAO_TYPE_LPMOD = -1
      !! Lone-pair occupancy, but engaged in a strong bond. A lone pair that is
      !! not quite one.
   integer, parameter :: QUAO_TYPE_NONE = 0
      !! Unclassified. Reached only on a transition metal, where the sigma/pi
      !! test is skipped -- see `skip_direction_test`.
   integer, parameter :: QUAO_TYPE_LONE_PAIR = 1
   integer, parameter :: QUAO_TYPE_RADICAL = 2
   integer, parameter :: QUAO_TYPE_SIGMA = 3
   integer, parameter :: QUAO_TYPE_PI = 4
   integer, parameter :: QUAO_TYPE_RDLP = 6      !! Reduced lone pair
   integer, parameter :: QUAO_TYPE_RDNV = 7      !! Reduced, and nearly empty
   integer, parameter :: QUAO_TYPE_NVMOD = 8     !! Empty, but in a strong bond
   integer, parameter :: QUAO_TYPE_NV = 9        !! Empty

   type :: quao_labels_t
      !! What each quasi-atomic orbital is, and what it is bonded to
      integer, allocatable :: orbital_type(:)
         !! (n_quao), one of the `QUAO_TYPE_*` codes
      integer, allocatable :: dominant_l(:)
         !! (n_quao), 0 to `QUAO_MAX_L`. Distinguishes an s lone pair from a p
         !! one, which is the difference between SLP and PLP in the report and
         !! is chemically the whole content of the label.
      real(dp), allocatable :: angular_character(:, :)
         !! (0:QUAO_MAX_L, n_quao), fractions summing to one. How much s, p, d
         !! and f each orbital is, measured against its own atom's free-atom
         !! orbitals.
      real(dp), allocatable :: direction(:, :)
         !! (3, n_quao), unit vectors. The axis the orbital is elongated along.
      integer, allocatable :: partner(:, :)
         !! (n_quao, n_quao). `partner(i, 1:partner_count(i))` are the orbitals
         !! `i` is strongly bonded to. More than one is normal: a sigma orbital
         !! in a strained ring has two.
      integer, allocatable :: partner_count(:)
      integer :: n_strong_bonds = 0
   end type quao_labels_t

contains

   subroutine label_quasi_atomic_orbitals(quao, mol, aambs, mixed, aambs_overlap, &
                                          atomic_numbers, coordinates, labels, error)
      !! Classify every quasi-atomic orbital
      !!
      !! Three things have to be known before a label can be put on an orbital:
      !! its occupation, which the construction already produced; which way it
      !! points, which needs second-moment integrals; and how much s and p
      !! character it has, which needs the projection back onto the free-atom
      !! orbitals. This assembles all three and then applies GAMESS's rules.
      type(quao_result_t), intent(in) :: quao
      type(libcint_molecule_t), intent(in) :: mol
         !! The orbital basis the SCF ran in, for the multipole integrals
      type(libcint_molecule_t), intent(in) :: aambs
         !! The free-atom minimal basis, for the angular momentum of each of its
         !! functions
      real(dp), intent(in) :: mixed(:, :)           !! < AO | AAMBS >
      real(dp), intent(in) :: aambs_overlap(:, :)   !! < AAMBS | AAMBS >
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)     !! (3, natm), Bohr
      type(quao_labels_t), intent(out) :: labels
      type(error_t), intent(inout) :: error

      if (error%has_error()) return
      if (.not. quao%oriented) then
         call error%set(ERROR_VALIDATION, "the orbital labels are defined for the "// &
                        "oriented quasi-atomic orbitals. Before orientation an atom's "// &
                        "bonding is spread arbitrarily over its orbitals, so there is "// &
                        "no single orbital to call the sigma bond.")
         return
      end if

      call quao_directions(mol, quao, coordinates, labels%direction, error)
      if (error%has_error()) return
      call quao_angular_character(quao, aambs, mixed, aambs_overlap, &
                                  labels%angular_character, labels%dominant_l, error)
      if (error%has_error()) return
      call assign_orbital_types(quao, atomic_numbers, coordinates, labels, error)
   end subroutine label_quasi_atomic_orbitals

   subroutine quao_directions(mol, quao, coordinates, direction, error)
      !! The axis each quasi-atomic orbital is elongated along
      !!
      !! The second moment of the orbital's own density about its own nucleus,
      !!
      !!     M_ab = < i | (r - R)_a (r - R)_b | i >
      !!
      !! diagonalized; the eigenvector of the largest eigenvalue is the direction
      !! the orbital sticks out in. For a p-like orbital that is unambiguous and
      !! for an s-like one it is noise -- which does not matter, because the only
      !! consumer is the sigma/pi test and s-like orbitals are classified as lone
      !! pairs before it runs.
      !!
      !! GAMESS gets the same axis from the electronic quadrupole moment, whose
      !! eigenvectors are the same as this tensor's -- a traceless quadrupole is
      !! `3M - tr(M)` and subtracting a multiple of the identity does not move an
      !! eigenvector. It arrives at the *smallest* eigenvalue rather than the
      !! largest because the electron charge is negative, which in its source
      !! costs three separate sign flips. This is the same statement without them.
      !!
      !! One dipole and one quadrupole build, not one per atom: the moment about
      !! `R` follows from the moments about a common origin by
      !! `M_ab(R) = <r_a r_b> - R_a <r_b> - R_b <r_a> + R_a R_b`, the last term
      !! using `<i|i> = 1`, which holds because the quasi-atomic orbitals are
      !! orthonormal.
      type(libcint_molecule_t), intent(in) :: mol
      type(quao_result_t), intent(in) :: quao
      real(dp), intent(in) :: coordinates(:, :)
      real(dp), allocatable, intent(out) :: direction(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dipole(:, :, :), quadrupole(:, :, :), work(:, :)
      real(dp), allocatable :: first(:, :), second(:, :)
      real(dp) :: origin(3), shift(3), moment(3, 3), values(3)
      integer :: i, a, b, comp, natm, info

      if (error%has_error()) return
      natm = size(coordinates, 2)

      ! The centroid, so the shifted terms stay the size of the molecule rather
      ! than the size of whatever coordinate frame the input happened to use.
      do a = 1, 3
         origin(a) = sum(coordinates(a, :))/real(natm, dp)
      end do

      call multipole_matrices(mol, origin, 1, dipole, error)
      if (error%has_error()) return
      call multipole_matrices(mol, origin, 2, quadrupole, error)
      if (error%has_error()) return

      allocate (first(3, quao%n_quao), second(9, quao%n_quao))
      allocate (work(mol%nao, quao%n_quao))
      do comp = 1, 3
         call pic_gemm(dipole(:, :, comp), quao%orbitals, work)
         do i = 1, quao%n_quao
            first(comp, i) = sum(quao%orbitals(:, i)*work(:, i))
         end do
      end do
      do comp = 1, 9
         call pic_gemm(quadrupole(:, :, comp), quao%orbitals, work)
         do i = 1, quao%n_quao
            second(comp, i) = sum(quao%orbitals(:, i)*work(:, i))
         end do
      end do

      allocate (direction(3, quao%n_quao))
      do i = 1, quao%n_quao
         shift = coordinates(:, quao%atom_of(i)) - origin
         do a = 1, 3
            do b = 1, 3
               moment(a, b) = second(3*(a - 1) + b, i) &
                              - shift(a)*first(b, i) - shift(b)*first(a, i) &
                              + shift(a)*shift(b)
            end do
         end do
         ! Symmetric to rounding, and `pic_syev` reads one triangle anyway.
         call pic_syev(moment, values, jobz="V", uplo="U", info=info)
         if (info /= 0) then
            call error%set(ERROR_VALIDATION, "the second-moment tensor of "// &
                           "quasi-atomic orbital "//to_char(i)//" could not be "// &
                           "diagonalized (info = "//to_char(info)//")")
            return
         end if
         direction(:, i) = moment(:, 3)
      end do

      deallocate (dipole, quadrupole, work, first, second)
   end subroutine quao_directions

   subroutine quao_angular_character(quao, aambs, mixed, aambs_overlap, &
                                     fraction, dominant, error)
      !! How much s, p, d and f character each quasi-atomic orbital has
      !!
      !! Measured against its own atom's free-atom orbitals: project the
      !! quasi-atomic orbital onto the orthogonalized minimal basis, keep the
      !! components on the atom it belongs to, and bin the squares by angular
      !! momentum. The result is normalized within the atom, so it reports the
      !! *shape* of the orbital and not how atomic it is -- that second number is
      !! `quao%atomic_character` and is a different question.
      !!
      !! The dominant momentum is what separates an s lone pair from a p one in
      !! the report, and that distinction is the whole chemical content of the
      !! label: a p lone pair points somewhere and can conjugate, an s one does
      !! not and cannot.
      type(quao_result_t), intent(in) :: quao
      type(libcint_molecule_t), intent(in) :: aambs
      real(dp), intent(in) :: mixed(:, :)
      real(dp), intent(in) :: aambs_overlap(:, :)
      real(dp), allocatable, intent(out) :: fraction(:, :)   !! (0:QUAO_MAX_L, n_quao)
      integer, allocatable, intent(out) :: dominant(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: projection(:, :)
      integer, allocatable :: function_l(:), function_atom(:)
      real(dp) :: total
      integer :: i, j, l, ish, iatom, ao, n_mbs

      if (error%has_error()) return
      call mo_aambs_overlap(quao%orbitals, mixed, aambs_overlap, projection, error)
      if (error%has_error()) return

      n_mbs = aambs%nao
      allocate (function_l(n_mbs), function_atom(n_mbs))
      do ish = 1, aambs%nbas
         l = aambs%bas(LIBCINT_ANG_OF, ish)
         iatom = aambs%bas(LIBCINT_ATOM_OF, ish) + 1
         do ao = aambs%shell_offset(ish) + 1, aambs%shell_offset(ish + 1)
            function_l(ao) = l
            function_atom(ao) = iatom
         end do
      end do

      if (maxval(function_l) > QUAO_MAX_L) then
         call error%set(ERROR_VALIDATION, "the atomic minimal basis has an l = "// &
                        to_char(maxval(function_l))//" shell, past the f functions "// &
                        "this labelling knows about. A free-atom occupied shell "// &
                        "should not go higher, so the basis file is wrong.")
         return
      end if

      allocate (fraction(0:QUAO_MAX_L, quao%n_quao), dominant(quao%n_quao))
      fraction = 0.0_dp
      do i = 1, quao%n_quao
         do j = 1, n_mbs
            if (function_atom(j) /= quao%atom_of(i)) cycle
            fraction(function_l(j), i) = fraction(function_l(j), i) + projection(i, j)**2
         end do
         total = sum(fraction(:, i))
         if (total <= 0.0_dp) then
            call error%set(ERROR_VALIDATION, "quasi-atomic orbital "//to_char(i)// &
                           " has no projection at all onto the free-atom orbitals of "// &
                           "the atom it was assigned to. The orbital-to-atom map and "// &
                           "the minimal-basis layout disagree.")
            return
         end if
         fraction(:, i) = fraction(:, i)/total
         dominant(i) = maxloc(fraction(:, i), dim=1) - 1
      end do

      deallocate (projection, function_l, function_atom)
   end subroutine quao_angular_character

   subroutine assign_orbital_types(quao, atomic_numbers, coordinates, labels, error)
      !! GAMESS's `LOCAL_ORIENBOS_ASSIGN_ORBTYP`, in order
      !!
      !! The rules run in a fixed sequence and each one only touches orbitals no
      !! earlier rule claimed, so the order is the algorithm:
      !!
      !! 1. Find the strong bonds -- interatomic pairs with a bond order at or
      !!    above `TOL_STRONG_BOND` -- in descending order of magnitude.
      !! 2. An orbital in a strong bond that is nearly doubly occupied is a lone
      !!    pair doing bonding work (LPMOD); one that is nearly empty is the
      !!    antibonding partner (NVMOD).
      !! 3. Any remaining nearly-doubly-occupied orbital is an ordinary lone pair.
      !! 4. The rest of each strong bond is sigma or pi, decided geometrically.
      !! 5. Anything left near single occupancy is a radical.
      !! 6. Everything else is sorted by occupation alone.
      !!
      !! **Sigma versus pi is decided by geometry, not by symmetry.** The test is
      !! the angle between the orbital's own axis and the line joining the two
      !! atoms, averaged over the two ends of the bond. That is why it needs the
      !! second moments, and why it says something even in a molecule with no
      !! symmetry at all -- where "sigma" and "pi" are not group-theoretic labels
      !! and cannot be assigned by any exact criterion.
      type(quao_result_t), intent(in) :: quao
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)
      type(quao_labels_t), intent(inout) :: labels
      type(error_t), intent(inout) :: error

      integer, allocatable :: bond_i(:), bond_j(:), order(:)
      real(dp), allocatable :: strength(:)
      real(dp) :: axis(3)
      real(dp) :: distance, along, occupation
      integer :: n, i, j, k, ib, n_pairs, n_strong

      if (error%has_error()) return
      n = quao%n_quao
      allocate (labels%orbital_type(n), labels%partner_count(n), labels%partner(n, n))
      labels%orbital_type = QUAO_TYPE_NONE
      labels%partner_count = 0
      labels%partner = 0

      ! ---- 1. the strong bonds, strongest first -----------------------------
      n_pairs = n*(n - 1)/2
      allocate (bond_i(n_pairs), bond_j(n_pairs), strength(n_pairs))
      n_strong = 0
      do i = 1, n
         do j = i + 1, n
            if (quao%atom_of(i) == quao%atom_of(j)) cycle
            if (abs(quao%population_bond_order(i, j)) < TOL_STRONG_BOND) cycle
            n_strong = n_strong + 1
            bond_i(n_strong) = i
            bond_j(n_strong) = j
            strength(n_strong) = abs(quao%population_bond_order(i, j))
         end do
      end do
      call sort_by_magnitude(strength(1:n_strong), order)
      labels%n_strong_bonds = n_strong

      ! ---- 2. lone pairs and empty orbitals that are nonetheless bonding ----
      do k = 1, n_strong
         ib = order(k)
         do j = 1, 2
            i = merge(bond_i(ib), bond_j(ib), j == 1)
            if (labels%orbital_type(i) /= QUAO_TYPE_NONE) cycle
            occupation = abs(quao%population_bond_order(i, i))
            if (occupation >= TOL_LONE_PAIR_OCC) labels%orbital_type(i) = QUAO_TYPE_LPMOD
            if (occupation < TOL_EMPTY) labels%orbital_type(i) = QUAO_TYPE_NVMOD
         end do
      end do

      ! ---- 3. ordinary lone pairs -------------------------------------------
      do i = 1, n
         if (labels%orbital_type(i) /= QUAO_TYPE_NONE) cycle
         if (abs(quao%population_bond_order(i, i)) < TOL_LONE_PAIR_OCC) cycle
         labels%orbital_type(i) = QUAO_TYPE_LONE_PAIR
      end do

      ! ---- 4. sigma and pi ---------------------------------------------------
      do k = 1, n_strong
         ib = order(k)
         i = bond_i(ib)
         j = bond_j(ib)
         if (labels%orbital_type(i) /= QUAO_TYPE_NONE .and. &
             labels%orbital_type(j) /= QUAO_TYPE_NONE) cycle

         axis = coordinates(:, quao%atom_of(i)) - coordinates(:, quao%atom_of(j))
         distance = sqrt(sum(axis**2))
         if (distance <= 0.0_dp) then
            call error%set(ERROR_VALIDATION, "atoms "//to_char(quao%atom_of(i))// &
                           " and "//to_char(quao%atom_of(j))//" are at the same point, "// &
                           "so the bond has no direction to compare an orbital with.")
            return
         end if
         axis = axis/distance
         ! Averaged over the two ends, so a bond does not come out sigma from one
         ! atom's point of view and pi from the other's.
         along = 0.5_dp*(abs(dot_product(axis, labels%direction(:, i))) &
                         + abs(dot_product(axis, labels%direction(:, j))))

         call set_bond_type(labels, quao, atomic_numbers, i, along)
         call set_bond_type(labels, quao, atomic_numbers, j, along)
      end do

      ! ---- 5. radicals -------------------------------------------------------
      do i = 1, n
         if (labels%orbital_type(i) /= QUAO_TYPE_NONE) cycle
         if (abs(1.0_dp - quao%population_bond_order(i, i)) > TOL_RADICAL_OCC) cycle
         labels%orbital_type(i) = QUAO_TYPE_RADICAL
      end do

      ! ---- 6. whatever is left, by occupation alone -------------------------
      do i = 1, n
         if (labels%orbital_type(i) /= QUAO_TYPE_NONE) cycle
         occupation = abs(quao%population_bond_order(i, i))
         if (occupation >= 1.0_dp) then
            labels%orbital_type(i) = QUAO_TYPE_RDLP
         else if (occupation >= TOL_EMPTY) then
            labels%orbital_type(i) = QUAO_TYPE_RDNV
         else
            labels%orbital_type(i) = QUAO_TYPE_NV
         end if
      end do

      ! ---- and who each bonding orbital is bonded to ------------------------
      ! Only the types that describe an interaction get a partner recorded. A
      ! lone pair keeps its parentheses empty even when it appears in a strong
      ! bond, which is the point: the report then shows a large bond order
      ! between an orbital that is bonded and one that is not, which is what
      ! delocalization looks like from here.
      do k = 1, n_strong
         ib = order(k)
         call add_partner(labels, bond_i(ib), bond_j(ib))
         call add_partner(labels, bond_j(ib), bond_i(ib))
      end do

      deallocate (bond_i, bond_j, strength, order)
   end subroutine assign_orbital_types

   subroutine set_bond_type(labels, quao, atomic_numbers, iorb, along)
      !! Sigma or pi for one end of one bond, if it is still unassigned
      type(quao_labels_t), intent(inout) :: labels
      type(quao_result_t), intent(in) :: quao
      integer, intent(in) :: atomic_numbers(:)
      integer, intent(in) :: iorb
      real(dp), intent(in) :: along     !! |cos| between orbital axis and bond axis

      if (labels%orbital_type(iorb) /= QUAO_TYPE_NONE) return
      if (skip_direction_test(atomic_numbers(quao%atom_of(iorb)))) return
      if (along >= TOL_OVERLAP) then
         labels%orbital_type(iorb) = QUAO_TYPE_SIGMA
      else
         labels%orbital_type(iorb) = QUAO_TYPE_PI
      end if
   end subroutine set_bond_type

   pure function skip_direction_test(atomic_number) result(skip)
      !! Whether an atom's orbitals are left unclassified
      !!
      !! The d block, in GAMESS's `LOCAL_ORIENBOS_SKIP`. A transition metal's
      !! bonding orbitals are d-heavy and a d orbital's second moment does not
      !! point along a bond in the way a p orbital's does -- a d(z^2) lobe and a
      !! d(xy) cloverleaf give quite different axes for orbitals doing the same
      !! job -- so sigma and pi are not assigned there at all rather than
      !! assigned wrongly. Those orbitals appear in the table with no type.
      integer, intent(in) :: atomic_number
      logical :: skip

      skip = (atomic_number >= 21 .and. atomic_number <= 30) .or. &
             (atomic_number >= 39 .and. atomic_number <= 48) .or. &
             (atomic_number >= 57 .and. atomic_number <= 80) .or. &
             (atomic_number >= 89)
   end function skip_direction_test

   subroutine add_partner(labels, iorb, jorb)
      !! Record `jorb` as something `iorb` is bonded to, if `iorb` is a bond
      type(quao_labels_t), intent(inout) :: labels
      integer, intent(in) :: iorb, jorb

      select case (labels%orbital_type(iorb))
      case (QUAO_TYPE_LPMOD, QUAO_TYPE_SIGMA, QUAO_TYPE_PI, QUAO_TYPE_NVMOD)
         labels%partner_count(iorb) = labels%partner_count(iorb) + 1
         labels%partner(iorb, labels%partner_count(iorb)) = jorb
      case default
         ! A lone pair, a radical or an empty orbital keeps its list empty even
         ! when it appears in a strong bond. That is deliberate rather than an
         ! omission: the report then shows a large bond order between an orbital
         ! that is bonded and one that is not, which is what delocalization looks
         ! like from here.
      end select
   end subroutine add_partner

   subroutine sort_by_magnitude(values, order)
      !! Indices of `values`, largest absolute value first
      real(dp), intent(in) :: values(:)
      integer, allocatable, intent(out) :: order(:)

      integer :: n, i, j, pick, swap

      n = size(values)
      allocate (order(n))
      do i = 1, n
         order(i) = i
      end do
      do i = 1, n - 1
         pick = i
         do j = i + 1, n
            if (abs(values(order(j))) > abs(values(order(pick)))) pick = j
         end do
         swap = order(i)
         order(i) = order(pick)
         order(pick) = swap
      end do
   end subroutine sort_by_magnitude

   pure function quao_type_name(orbital_type, dominant_l) result(name)
      !! The label the report prints
      !!
      !! `dominant_l` only matters for the lone-pair types, where it is the
      !! difference between SLP and PLP.
      !!
      !! GAMESS prints nothing at all for the types past PI, and prints its
      !! no-partner marker `NWB   0` for the unclassified ones. Both are given
      !! names here: a blank column in a table means the code forgot, and there
      !! is no way to tell that from an orbital that genuinely has no type.
      integer, intent(in) :: orbital_type
      integer, intent(in) :: dominant_l
      character(len=7) :: name

      select case (orbital_type)
      case (QUAO_TYPE_LPMOD)
         name = shell_letter(dominant_l)//"LPMOD"
      case (QUAO_TYPE_LONE_PAIR)
         name = shell_letter(dominant_l)//"LP"
      case (QUAO_TYPE_RADICAL)
         name = "RADICAL"
      case (QUAO_TYPE_SIGMA)
         name = "SIGMA"
      case (QUAO_TYPE_PI)
         name = "PI"
      case (QUAO_TYPE_RDLP)
         name = "RDLP"
      case (QUAO_TYPE_RDNV)
         name = "RDNV"
      case (QUAO_TYPE_NVMOD)
         name = "NVMOD"
      case (QUAO_TYPE_NV)
         name = "NV"
      case default
         name = "NOTYPE"
      end select
   end function quao_type_name

   pure function shell_letter(l) result(letter)
      integer, intent(in) :: l
      character(len=1) :: letter

      select case (l)
      case (0)
         letter = "S"
      case (1)
         letter = "P"
      case (2)
         letter = "D"
      case (3)
         letter = "F"
      case default
         letter = "?"
      end select
   end function shell_letter

   ! ---------------------------------------------------------------------------
   ! The report
   ! ---------------------------------------------------------------------------

   subroutine print_quao_report(verbose, quao, labels, interference, &
                                element_symbols, n_core)
      !! Every orbital pair worth looking at, strongest kinetic bond first
      !!
      !! Three blocks, the same three GAMESS prints and in the same order: pairs
      !! sorted by kinetic bond order, orbitals sorted by occupation, and the
      !! s/p/d/f composition of each orbital. The columns carry the same
      !! quantities under the same headings; the spacing is this code's, so the
      !! two tables are read side by side rather than diffed.
      !!
      !! The sort is on the kinetic bond order and not on the bond order,
      !! deliberately. A bond order is a population and says how much density
      !! two orbitals share; the kinetic bond order is an energy and says what
      !! that sharing is worth. They disagree about ordering more often than one
      !! would expect -- in formyl chloride the C-Cl bond order beats the C-H one
      !! and its kinetic bond order beats it by more than twice as much, because
      !! the C-H kinetic integral is smaller.
      logical, intent(in) :: verbose
      type(quao_result_t), intent(in) :: quao
      type(quao_labels_t), intent(in) :: labels
      real(dp), intent(in) :: interference(:, :)
         !! `p * T` in hartree, Paper II eq (1) -- the second output of
         !! `kinetic_bond_orders`, not the kcal/mol one. This is the column
         !! GAMESS heads KEI-BO, so the report quotes the same quantity.
      character(len=*), intent(in) :: element_symbols(:)
      integer, intent(in) :: n_core
         !! Chemical-core orbitals, which the construction never sees. Added to
         !! every printed orbital index so the numbering matches the molecular
         !! orbitals it came from -- and matches GAMESS, which does the same.

      real(dp), allocatable :: magnitude(:), occupation(:)
      integer, allocatable :: pair_i(:), pair_j(:), order(:)
      character(len=512) :: line
      ! Field by field, so a width here and a width in `append_orbital` cannot
      ! be changed independently without the mismatch being obvious.
      character(len=*), parameter :: HEADER_I = &
                                     " ORB I"//"       OCC I"//"  ATM I"// &
                                     "   BONDED TO "//" ORBTYP "
      character(len=*), parameter :: HEADER_J = &
                                     " ORB J"//"       OCC J"//"  ATM J"// &
                                     "   BONDED TO "//" ORBTYP "
      real(dp) :: total
      integer :: n, i, j, k, np, l, pos

      if (.not. verbose) return
      n = quao%n_quao

      ! ---- by kinetic bond order --------------------------------------------
      call logger%info("")
      call logger%info("  quasi-atomic orbital pairs, by kinetic bond order")
      call logger%info("")
      call logger%info("    BOND ORDER      KEI-BO"//HEADER_I//HEADER_J)

      allocate (pair_i(n*(n - 1)/2), pair_j(n*(n - 1)/2), magnitude(n*(n - 1)/2))
      np = 0
      do i = 1, n
         do j = i + 1, n
            if (abs(interference(i, j)) < TOL_PRINT) cycle
            np = np + 1
            pair_i(np) = i
            pair_j(np) = j
            magnitude(np) = interference(i, j)
         end do
      end do
      call sort_by_magnitude(magnitude(1:np), order)

      do k = 1, np
         i = pair_i(order(k))
         j = pair_j(order(k))
         call start_row(line, pos, quao%population_bond_order(i, j), interference(i, j))
         call append_orbital(line, pos, quao, labels, element_symbols, n_core, i)
         call append_orbital(line, pos, quao, labels, element_symbols, n_core, j)
         call logger%info(trim(line))
      end do
      if (np == 0) call logger%info("    (no pair above the print threshold)")

      ! ---- by occupation -----------------------------------------------------
      call logger%info("")
      call logger%info("  quasi-atomic orbitals, by occupation")
      call logger%info("")
      call logger%info("    OCCUPATION   KEI-OCC I"//" ORB I"//"  ATM I"// &
                       "   BONDED TO "//" ORBTYP ")

      allocate (occupation(n))
      do i = 1, n
         occupation(i) = quao%population_bond_order(i, i)
      end do
      deallocate (order)
      call sort_by_magnitude(occupation, order)
      do k = 1, n
         i = order(k)
         if (abs(occupation(i)) < TOL_PRINT) cycle
         call start_row(line, pos, quao%population_bond_order(i, i), interference(i, i))
         call append_orbital(line, pos, quao, labels, element_symbols, n_core, i, &
                             with_occupation=.false.)
         call logger%info(trim(line))
      end do

      ! The kinetic interference between orbitals on *different* atoms, which is
      ! the part of the kinetic energy that exists only because the atoms are
      ! bonded. Intra-atomic terms are excluded because a free atom has them too.
      total = 0.0_dp
      do i = 1, n
         do j = 1, n
            if (quao%atom_of(i) == quao%atom_of(j)) cycle
            total = total + interference(i, j)
         end do
      end do
      write (line, "(a,f20.10)") "  interatomic KEI-BO sum ", total
      call logger%info("")
      call logger%info(trim(line))

      ! ---- composition -------------------------------------------------------
      call logger%info("")
      call logger%info("  quasi-atomic orbital composition")
      call logger%info("")
      call logger%info("     ORB I    PERCENT S    PERCENT P    PERCENT D    PERCENT F")
      do i = 1, n
         line = ""
         write (line(1:10), "(i10)") n_core + i
         pos = 11
         do l = 0, QUAO_MAX_L
            write (line(pos:pos + 12), "(f13.7)") labels%angular_character(l, i)
            pos = pos + 13
         end do
         call logger%info(trim(line))
      end do

      deallocate (pair_i, pair_j, magnitude, occupation, order)
   end subroutine print_quao_report

   subroutine start_row(line, pos, first, second)
      !! Blank the line and lay down the two leading numbers
      character(len=*), intent(out) :: line
      integer, intent(out) :: pos
      real(dp), intent(in) :: first, second

      line = ""
      write (line(3:2 + 2*WIDTH_NUMBER), "(2f12.7)") first, second
      pos = 3 + 2*WIDTH_NUMBER
   end subroutine start_row

   subroutine append_orbital(line, pos, quao, labels, element_symbols, n_core, iorb, &
                             with_occupation)
      !! One orbital's columns: index, occupation, atom, partners, type
      !!
      !! The cursor is carried rather than recovered with `len_trim`, because a
      !! fixed-width field that ends in blanks -- which the orbital type usually
      !! does -- cannot be found again that way, and the columns after it drift
      !! by however long the last label happened to be.
      character(len=*), intent(inout) :: line
      integer, intent(inout) :: pos
      type(quao_result_t), intent(in) :: quao
      type(quao_labels_t), intent(in) :: labels
      character(len=*), intent(in) :: element_symbols(:)
      integer, intent(in) :: n_core
      integer, intent(in) :: iorb
      logical, intent(in), optional :: with_occupation
         !! Default true. The occupation table has printed it already.

      character(len=WIDTH_PARTNERS) :: partners
      character(len=8) :: field
      integer :: k, iatom
      logical :: occupation

      occupation = .true.
      if (present(with_occupation)) occupation = with_occupation

      write (line(pos:pos + WIDTH_ORB - 1), "(i6)") n_core + iorb
      pos = pos + WIDTH_ORB
      if (occupation) then
         write (line(pos:pos + WIDTH_NUMBER - 1), "(f12.7)") &
            quao%population_bond_order(iorb, iorb)
         pos = pos + WIDTH_NUMBER
      end if

      iatom = quao%atom_of(iorb)
      write (line(pos:pos + WIDTH_ATOM - 1), "(1x,a2,i4)") &
         adjustl(element_symbols(iatom)), iatom
      pos = pos + WIDTH_ATOM

      ! The atoms this orbital is bonded to. Overflows its column when an
      ! orbital has three or more partners, which is rare and is better than
      ! truncating the list.
      if (labels%partner_count(iorb) == 0) then
         partners = "(  --  )"
      else
         partners = "("
         do k = 1, labels%partner_count(iorb)
            if (k > 1) partners = trim(partners)//","
            iatom = quao%atom_of(labels%partner(iorb, k))
            write (field, "(1x,a2,i4)") adjustl(element_symbols(iatom)), iatom
            partners = trim(partners)//field
         end do
         partners = trim(partners)//" )"
      end if
      write (line(pos:pos + WIDTH_PARTNERS), "(1x,a)") partners
      pos = pos + WIDTH_PARTNERS + 1

      write (line(pos:pos + WIDTH_TYPE - 1), "(1x,a7)") &
         quao_type_name(labels%orbital_type(iorb), labels%dominant_l(iorb))
      pos = pos + WIDTH_TYPE
   end subroutine append_orbital

end module mqc_libcint_quao_report
