!! What the quasi-atomic orbitals turn out to be, and the table that says so
module mqc_czt_quao_report
   !! The construction in `mqc_czt_quao` produces a set of orbitals, one per
   !! atom, and a matrix of populations and bond orders, but does not say which
   !! orbital is a sigma bond and which is a lone pair. This module assigns
   !! those labels and prints the report: each orbital pair above a threshold,
   !! sorted by magnitude, annotated with its atom, occupation, bonding partner
   !! and kind.
   !!
   !! **Every label here is empirical.** Nothing defines a sigma orbital; the
   !! classification is a set of thresholds on occupations, bond orders and
   !! directions, whose values were chosen so that the labels reproduce the
   !! ones GAMESS prints. An orbital half a percent the wrong side of
   !! `TOL_LONE_PAIR_OCC` is not a different orbital, it is a different label
   !! on the same one.
   !!
   !! The quantities being labelled -- populations, bond orders, kinetic bond
   !! orders -- are not empirical, and are validated against GAMESS to 1e-6 in
   !! `test/test_mqc_quao.f90`.
   use pic_types, only: dp
   use pic_io, only: to_char
   use pic_logger, only: logger => global_logger
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_multipole, only: multipole_matrices
   use mqc_czt_quao, only: quao_result_t, mo_aambs_overlap
   use mqc_quao_rows, only: quao_rows_t, quao_type_name, QUAO_ROW_BOND, &
                            QUAO_ROW_DELOCALIZATION, QUAO_TYPE_LPMOD, QUAO_TYPE_NONE, &
                            QUAO_TYPE_LONE_PAIR, QUAO_TYPE_RADICAL, QUAO_TYPE_SIGMA, &
                            QUAO_TYPE_PI, QUAO_TYPE_RDLP, QUAO_TYPE_RDNV, QUAO_TYPE_NVMOD, &
                            QUAO_TYPE_NV
   use libcint_fortran, only: LIBCINT_ANG_OF, LIBCINT_ATOM_OF
   implicit none
   private

   public :: quao_labels_t
   public :: label_quasi_atomic_orbitals
   public :: print_quao_report
   public :: build_quao_rows
   public :: quao_type_name
   public :: QUAO_TYPE_LPMOD, QUAO_TYPE_NONE, QUAO_TYPE_LONE_PAIR, QUAO_TYPE_RADICAL
   public :: QUAO_TYPE_SIGMA, QUAO_TYPE_PI, QUAO_TYPE_RDLP, QUAO_TYPE_RDNV
   public :: QUAO_TYPE_NVMOD, QUAO_TYPE_NV
   public :: QUAO_MAX_L

   ! The labelling thresholds. Their values reproduce GAMESS's printed labels
   ! and are kept for that comparison; what each one decides is below.
   real(dp), parameter :: TOL_STRONG_BOND = 0.40_dp
      !! A bond order between two atoms at or above this makes the pair a
      !! "strong bond", and only strong bonds get orbitals labelled from them.
      !! Everything else in the table is reported without being classified.
   real(dp), parameter :: TOL_LONE_PAIR_OCC = 1.80_dp
      !! An orbital this close to doubly occupied is a lone pair, whatever else
      !! it is doing.
   real(dp), parameter :: TOL_RADICAL_OCC = 0.40_dp
      !! How far from singly occupied an orbital may be and still be called a
      !! radical.
   real(dp), parameter :: TOL_OVERLAP = 0.20_dp
      !! The cosine between an orbital's own axis and the line joining the two
      !! atoms, above which the orbital points along the bond and is therefore
      !! sigma rather than pi. This is a geometric test, not a symmetry one --
      !! see `assign_orbital_types`.
   real(dp), parameter :: TOL_EMPTY = 0.40_dp
      !! Occupation below which an orbital is empty for labelling purposes.
   real(dp), parameter :: TOL_PRINT = 1.0e-3_dp
      !! Rows below this in magnitude are not printed.

   ! Column widths for the report, so the header and the rows cannot drift
   ! apart. Every field is written into a substring of exactly this length.
   integer, parameter :: WIDTH_NUMBER = 12     !! A bond order or an occupation
   integer, parameter :: WIDTH_ORB = 6         !! An orbital index
   integer, parameter :: WIDTH_ATOM = 7        !! A symbol and an atom index
   integer, parameter :: WIDTH_PARTNERS = 12   !! The bonded-to list, plus a lead space
   integer, parameter :: WIDTH_TYPE = 8        !! The orbital type, plus a lead space

   integer, parameter :: QUAO_MAX_L = 3
      !! s, p, d, f. The accurate atomic minimal basis has nothing higher, since it holds free-atom occupied shells only.

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
      !! Three things go into a label: the occupation, which the construction
      !! already produced; which way the orbital points, which needs
      !! second-moment integrals; and how much s and p character it has, which
      !! needs the projection back onto the free-atom orbitals. Refused unless
      !! `quao` is oriented.
      type(quao_result_t), intent(in) :: quao
      type(czt_molecule_t), intent(in) :: mol
         !! The orbital basis the SCF ran in, for the multipole integrals
      type(czt_molecule_t), intent(in) :: aambs
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
      !! diagonalized; the eigenvector of the largest eigenvalue is the
      !! direction the orbital sticks out in. For an s-like orbital it is
      !! noise, which does not matter: the only consumer is the sigma/pi test,
      !! and an s-like orbital is normally a lone pair, which is labelled from
      !! its occupation without consulting the axis.
      !!
      !! The axis taken from the electronic quadrupole moment is the same one:
      !! a traceless quadrupole is `3M - tr(M)`, and subtracting a multiple of
      !! the identity does not move an eigenvector. There it belongs to the
      !! *smallest* eigenvalue, because the electron charge is negative.
      !!
      !! One dipole and one quadrupole build, not one per atom: the moment about
      !! `R` follows from the moments about a common origin by
      !! `M_ab(R) = <r_a r_b> - R_a <r_b> - R_b <r_a> + R_a R_b`, the last term
      !! using `<i|i> = 1`, which holds because the quasi-atomic orbitals are
      !! orthonormal.
      type(czt_molecule_t), intent(in) :: mol
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
      ! than of whatever coordinate frame the input used.
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
      !! *shape* of the orbital and not how atomic it is -- that second number
      !! is `quao%atomic_character`.
      !!
      !! The dominant momentum separates an s lone pair from a p one in the
      !! report: a p lone pair points somewhere and can conjugate, an s one
      !! does not.
      type(quao_result_t), intent(in) :: quao
      type(czt_molecule_t), intent(in) :: aambs
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
      !! Label every orbital from its own occupation and its strong bonds
      !!
      !! A *strong bond* is an interatomic orbital pair whose bond order reaches
      !! `TOL_STRONG_BOND`. They are ranked once, largest first, and each
      !! orbital's list of strong partners follows that ranking, so its first
      !! partner is its strongest. With that list in hand every orbital is
      !! labelled on its own by `classify_orbital`; no label depends on another
      !! orbital's label.
      !!
      !! **Sigma versus pi is decided by geometry, not by symmetry.** The test
      !! is the angle between the orbitals' axes and the line joining the two
      !! atoms of the orbital's strongest bond, averaged over the two ends, so
      !! it says something even in a molecule with no symmetry at all.
      type(quao_result_t), intent(in) :: quao
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)
      type(quao_labels_t), intent(inout) :: labels
      type(error_t), intent(inout) :: error

      integer, allocatable :: ends(:, :), ranking(:), n_bonded(:), bonded(:, :)
      real(dp), allocatable :: bond_order(:)
      real(dp) :: alignment
      integer :: n_orb, orb, other, rank, n_bonds
      logical :: axis_known

      if (error%has_error()) return
      n_orb = quao%n_quao

      call find_strong_bonds(quao, ends, bond_order)
      n_bonds = size(bond_order)
      call sort_by_magnitude(bond_order, ranking)

      allocate (n_bonded(n_orb), bonded(n_orb, n_orb))
      n_bonded = 0
      bonded = 0
      do rank = 1, n_bonds
         associate (one => ends(1, ranking(rank)), two => ends(2, ranking(rank)))
            n_bonded(one) = n_bonded(one) + 1
            bonded(one, n_bonded(one)) = two
            n_bonded(two) = n_bonded(two) + 1
            bonded(two, n_bonded(two)) = one
         end associate
      end do

      allocate (labels%orbital_type(n_orb), labels%partner_count(n_orb))
      allocate (labels%partner(n_orb, n_orb))
      labels%n_strong_bonds = n_bonds
      labels%partner_count = 0
      labels%partner = 0
      do orb = 1, n_orb
         axis_known = .not. axis_is_ambiguous(atomic_numbers(quao%atom_of(orb)))
         alignment = 0.0_dp
         if (n_bonded(orb) > 0 .and. axis_known) then
            other = bonded(orb, 1)
            call bond_alignment(quao, coordinates, labels%direction, orb, other, &
                                alignment, error)
            if (error%has_error()) return
         end if
         labels%orbital_type(orb) = classify_orbital(quao%population_bond_order(orb, orb), &
                                                     n_bonded(orb) > 0, axis_known, &
                                                     alignment)
         ! Only an orbital labelled as taking part in a bond lists what it is
         ! bonded to. A lone pair, radical or empty orbital that nonetheless
         ! has a large bond order keeps an empty list, and the report shows a
         ! bonded orbital against one that is not: delocalization.
         if (records_partners(labels%orbital_type(orb))) then
            labels%partner_count(orb) = n_bonded(orb)
            labels%partner(orb, :) = bonded(orb, :)
         end if
      end do
   end subroutine assign_orbital_types

   subroutine find_strong_bonds(quao, ends, bond_order)
      !! Every interatomic orbital pair whose bond order reaches `TOL_STRONG_BOND`
      type(quao_result_t), intent(in) :: quao
      integer, allocatable, intent(out) :: ends(:, :)      !! (2, n_bonds), lower index first
      real(dp), allocatable, intent(out) :: bond_order(:)  !! (n_bonds), magnitudes

      integer, allocatable :: found(:, :)
      real(dp), allocatable :: size_of(:)
      integer :: a, b, n_found

      allocate (found(2, max(1, quao%n_quao*(quao%n_quao - 1)/2)))
      allocate (size_of(size(found, 2)))
      n_found = 0
      do a = 1, quao%n_quao
         do b = a + 1, quao%n_quao
            if (quao%atom_of(a) == quao%atom_of(b)) cycle
            if (abs(quao%population_bond_order(a, b)) < TOL_STRONG_BOND) cycle
            n_found = n_found + 1
            found(:, n_found) = [a, b]
            size_of(n_found) = abs(quao%population_bond_order(a, b))
         end do
      end do
      ends = found(:, 1:n_found)
      bond_order = size_of(1:n_found)
   end subroutine find_strong_bonds

   subroutine bond_alignment(quao, coordinates, direction, orb, other, alignment, error)
      !! How nearly two bonded orbitals point along the line between their atoms
      !!
      !! The mean of `|cos|` between each orbital's axis and the internuclear
      !! axis. Taken over both ends, so a bond cannot be sigma from one atom's
      !! side and pi from the other's.
      type(quao_result_t), intent(in) :: quao
      real(dp), intent(in) :: coordinates(:, :)
      real(dp), intent(in) :: direction(:, :)   !! (3, n_quao), unit vectors
      integer, intent(in) :: orb, other
      real(dp), intent(out) :: alignment
      type(error_t), intent(inout) :: error

      real(dp) :: axis(3)
      real(dp) :: length

      alignment = 0.0_dp
      axis = coordinates(:, quao%atom_of(orb)) - coordinates(:, quao%atom_of(other))
      length = norm2(axis)
      if (length <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "atoms "//to_char(quao%atom_of(orb))// &
                        " and "//to_char(quao%atom_of(other))//" are at the same point, "// &
                        "so the bond has no direction to compare an orbital with.")
         return
      end if
      axis = axis/length
      alignment = 0.5_dp*(abs(dot_product(axis, direction(:, orb))) + &
                          abs(dot_product(axis, direction(:, other))))
   end subroutine bond_alignment

   pure function classify_orbital(population, strongly_bonded, axis_known, alignment) &
      result(code)
      !! The `QUAO_TYPE_*` code of one orbital
      !!
      !! Occupation is read on the 0-to-2 scale of a spatial orbital.
      !!
      !! | in a strong bond | occupation             | label              |
      !! |------------------|------------------------|--------------------|
      !! | yes              | `>= TOL_LONE_PAIR_OCC` | LPMOD (bonding lone pair) |
      !! | yes              | `<  TOL_EMPTY`         | NVMOD (bonding empty)     |
      !! | yes              | between, axis known    | sigma or pi by `alignment` |
      !! | no               | `>= TOL_LONE_PAIR_OCC` | lone pair          |
      !!
      !! Anything else -- including a bonded orbital on an atom whose axis is
      !! not trusted -- is labelled by occupation alone: radical within
      !! `TOL_RADICAL_OCC` of one electron, otherwise RDLP from one electron up,
      !! RDNV down to `TOL_EMPTY`, and nonvalence below it.
      real(dp), intent(in) :: population    !! The orbital's diagonal density element
      logical, intent(in) :: strongly_bonded
      logical, intent(in) :: axis_known
      real(dp), intent(in) :: alignment     !! From `bond_alignment`, its strongest bond
      integer :: code

      real(dp) :: occupation

      occupation = abs(population)
      if (strongly_bonded) then
         if (occupation >= TOL_LONE_PAIR_OCC) then
            code = QUAO_TYPE_LPMOD
            return
         end if
         if (occupation < TOL_EMPTY) then
            code = QUAO_TYPE_NVMOD
            return
         end if
         if (axis_known) then
            code = merge(QUAO_TYPE_SIGMA, QUAO_TYPE_PI, alignment >= TOL_OVERLAP)
            return
         end if
      else if (occupation >= TOL_LONE_PAIR_OCC) then
         code = QUAO_TYPE_LONE_PAIR
         return
      end if

      if (abs(population - 1.0_dp) <= TOL_RADICAL_OCC) then
         code = QUAO_TYPE_RADICAL
      else if (occupation >= 1.0_dp) then
         code = QUAO_TYPE_RDLP
      else if (occupation >= TOL_EMPTY) then
         code = QUAO_TYPE_RDNV
      else
         code = QUAO_TYPE_NV
      end if
   end function classify_orbital

   pure function axis_is_ambiguous(atomic_number) result(ambiguous)
      !! Whether an atom's orbitals get no sigma/pi label from their axis
      !!
      !! Atoms with d or f valence shells: transition metals, lanthanides and
      !! actinides. The second-moment axis of a d orbital is not a reliable
      !! guide to its bonding role -- a d(z^2) lobe and a d(xy) cloverleaf doing
      !! the same job point in unrelated directions -- so such an orbital is
      !! labelled by its occupation instead.
      integer, intent(in) :: atomic_number
      logical :: ambiguous

      select case (atomic_number)
      case (21:30, 39:48, 57:80, 89:)
         ambiguous = .true.
      case default
         ambiguous = .false.
      end select
   end function axis_is_ambiguous

   pure function records_partners(code) result(records)
      !! Whether an orbital with this label lists the orbitals it is bonded to
      integer, intent(in) :: code
      logical :: records

      select case (code)
      case (QUAO_TYPE_LPMOD, QUAO_TYPE_SIGMA, QUAO_TYPE_PI, QUAO_TYPE_NVMOD)
         records = .true.
      case default
         records = .false.
      end select
   end function records_partners

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

   ! ---------------------------------------------------------------------------
   ! The report
   ! ---------------------------------------------------------------------------

   subroutine build_quao_rows(quao, labels, interference, kinetic_bond_order, n_core, &
                              n_atoms, threshold, rows)
      !! The bonds and delocalization tables of `print_quao_report`, as data
      !!
      !! The same pairs in the same order: every pair whose kinetic
      !! interference reaches `TOL_PRINT` hartree, strongest first, bonds before
      !! delocalization, and delocalization weaker than `threshold` dropped. A
      !! delocalization row has its donor, the fuller orbital, at end 1. Bond
      !! orders are stored as magnitudes, their sign being a phase convention.
      !!
      !! The atom-pair sums cover every orbital pair, above threshold or not.
      type(quao_result_t), intent(in) :: quao
      type(quao_labels_t), intent(in) :: labels
      real(dp), intent(in) :: interference(:, :)
         !! (n_quao, n_quao), hartree
      real(dp), intent(in) :: kinetic_bond_order(:, :)
         !! (n_quao, n_quao), kcal/mol
      integer, intent(in) :: n_core
         !! Chemical-core orbitals, added to every orbital index
      integer, intent(in) :: n_atoms
      real(dp), intent(in) :: threshold
         !! kcal/mol
      type(quao_rows_t), intent(out) :: rows

      integer, allocatable :: pair_i(:), pair_j(:), order(:), kept(:)
      logical, allocatable :: bond(:)
      integer :: n, np, i, j, k, row, pass, side, orb, iatom, jatom

      n = quao%n_quao
      call ranked_pairs(interference, n, pair_i, pair_j, order, np)
      allocate (bond(np), kept(np))
      do k = 1, np
         bond(k) = is_bond(quao, labels, pair_i(order(k)), pair_j(order(k)))
      end do

      ! Bonds first, then delocalization, each in the ranked order.
      rows%n = 0
      do pass = 1, 2
         do k = 1, np
            i = pair_i(order(k))
            j = pair_j(order(k))
            if (pass == 1 .neqv. bond(k)) cycle
            if (pass == 2 .and. abs(kinetic_bond_order(i, j)) < threshold) cycle
            rows%n = rows%n + 1
            kept(rows%n) = k
         end do
      end do

      allocate (rows%kind(rows%n), rows%orbital(2, rows%n), rows%atom(2, rows%n), &
                rows%orbital_type(2, rows%n), rows%dominant_l(2, rows%n), &
                rows%partner_atom(2, rows%n), rows%occupation(2, rows%n), &
                rows%bond_order(rows%n), rows%kinetic_bond_order(rows%n))
      do row = 1, rows%n
         k = kept(row)
         i = pair_i(order(k))
         j = pair_j(order(k))
         if (bond(k)) then
            rows%kind(row) = QUAO_ROW_BOND
         else
            rows%kind(row) = QUAO_ROW_DELOCALIZATION
            if (quao%population_bond_order(j, j) > quao%population_bond_order(i, i)) then
               i = pair_j(order(k))
               j = pair_i(order(k))
            end if
         end if
         do side = 1, 2
            orb = i
            if (side == 2) orb = j
            rows%orbital(side, row) = n_core + orb
            rows%atom(side, row) = quao%atom_of(orb)
            rows%orbital_type(side, row) = labels%orbital_type(orb)
            rows%dominant_l(side, row) = labels%dominant_l(orb)
            rows%partner_atom(side, row) = 0
            if (labels%partner_count(orb) > 0) then
               rows%partner_atom(side, row) = quao%atom_of(labels%partner(orb, 1))
            end if
            rows%occupation(side, row) = quao%population_bond_order(orb, orb)
         end do
         rows%bond_order(row) = abs(quao%population_bond_order(i, j))
         rows%kinetic_bond_order(row) = kinetic_bond_order(i, j)
      end do

      allocate (rows%atom_bond_index(n_atoms, n_atoms), &
                rows%atom_kinetic_bond_order(n_atoms, n_atoms))
      rows%atom_bond_index = 0.0_dp
      rows%atom_kinetic_bond_order = 0.0_dp
      do j = 1, n
         jatom = quao%atom_of(j)
         do i = 1, n
            iatom = quao%atom_of(i)
            if (iatom == jatom) cycle
            rows%atom_bond_index(iatom, jatom) = rows%atom_bond_index(iatom, jatom) + &
                                                 quao%population_bond_order(i, j)**2
            rows%atom_kinetic_bond_order(iatom, jatom) = &
               rows%atom_kinetic_bond_order(iatom, jatom) + kinetic_bond_order(i, j)
         end do
      end do
      rows%threshold = threshold
      rows%orientation_stalled = quao%orientation_stalled
      deallocate (pair_i, pair_j, order, bond, kept)
   end subroutine build_quao_rows

   subroutine ranked_pairs(interference, n, pair_i, pair_j, order, np)
      !! Orbital pairs above `TOL_PRINT`, and their order strongest first
      !!
      !! `pair_i(order(k))`, `pair_j(order(k))` is the k-th strongest pair,
      !! with `pair_i < pair_j`. Shared by the report and the rows so the two
      !! cannot select differently.
      real(dp), intent(in) :: interference(:, :)
      integer, intent(in) :: n
      integer, allocatable, intent(out) :: pair_i(:), pair_j(:), order(:)
      integer, intent(out) :: np

      real(dp), allocatable :: magnitude(:)
      integer :: i, j

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
      deallocate (magnitude)
   end subroutine ranked_pairs

   subroutine print_quao_report(verbose, quao, labels, interference, &
                                element_symbols, n_core, kinetic_bond_order, threshold)
      !! The bonding picture, grouped by what each number means
      !!
      !! **Not GAMESS's layout**, which is one row per orbital pair with every
      !! column duplicated under an `I` and a `J` heading. The same information
      !! is regrouped into three blocks. **Bonds** are pairs where both ends are
      !! bonding orbitals pointing at each other -- the covalent skeleton, one
      !! row per bond. **Delocalization** is everything else above threshold,
      !! where a lone pair donates into a bond. **Orbitals** is the per-orbital
      !! table, with the hybridisation folded in.
      !!
      !! Both interaction blocks are sorted on the kinetic bond order rather
      !! than the bond order: a bond order is a population and says how much
      !! density two orbitals share, the kinetic bond order is an energy and
      !! says what that sharing is worth. The two disagree about ordering.
      logical, intent(in) :: verbose
      type(quao_result_t), intent(in) :: quao
      type(quao_labels_t), intent(in) :: labels
      real(dp), intent(in) :: interference(:, :)
         !! `p * T` in hartree, Paper II eq (1). What GAMESS prints as KEI-BO,
         !! kept so its table and this one can still be compared.
      character(len=*), intent(in) :: element_symbols(:)
      integer, intent(in) :: n_core
         !! Chemical-core orbitals, which the construction never sees. Added to
         !! every printed orbital index so the numbering matches the molecular
         !! orbitals it came from, as GAMESS's does.
      real(dp), intent(in), optional :: kinetic_bond_order(:, :)
         !! Paper II eq (2), in kcal/mol. Shown as the headline energy when
         !! present, that being the scale comparable with tabulated bond
         !! energies.
      real(dp), intent(in), optional :: threshold
         !! kcal/mol. Pairs weaker than this are counted and not printed.
         !! Ignored without `kinetic_bond_order`, that being the quantity it is
         !! a threshold on; the `TOL_PRINT` hartree floor is then all there is.

      real(dp), allocatable :: occupation(:)
      integer, allocatable :: pair_i(:), pair_j(:), order(:)
      character(len=160) :: line
      character(len=24) :: left, right
      character(len=18) :: field_a, field_b
      real(dp) :: total, energy, cutoff
      integer :: n, i, j, k, np, printed, printed_i, printed_j, suppressed
      logical :: bonded

      if (.not. verbose) return
      n = quao%n_quao
      ! A threshold on the kcal/mol quantity needs that quantity; without it
      ! the hartree floor is what filters.
      cutoff = 0.0_dp
      if (present(threshold) .and. present(kinetic_bond_order)) cutoff = threshold

      call ranked_pairs(interference, n, pair_i, pair_j, order, np)

      ! ---- the covalent skeleton --------------------------------------------
      call logger%info("")
      call logger%info("  bonds")
      call logger%info("     bond              type        order    kcal/mol   orbitals")
      printed = 0
      do k = 1, np
         i = pair_i(order(k))
         j = pair_j(order(k))
         if (.not. is_bond(quao, labels, i, j)) cycle
         printed = printed + 1
         call atom_text(quao, element_symbols, i, left)
         call atom_text(quao, element_symbols, j, right)
         energy = interference(i, j)
         if (present(kinetic_bond_order)) energy = kinetic_bond_order(i, j)
         ! Assigned to fixed-length buffers rather than written with an `a16`
         ! edit descriptor: that pads a short string on the *left*, so the
         ! column comes out ragged. Assignment pads on the right.
         field_a = trim(left)//" - "//trim(right)
         field_b = quao_type_name(labels%orbital_type(i), labels%dominant_l(i))
         write (line, "(4x,a18,a11,f10.4,f11.2,i7,a,i0)") field_a, field_b, &
            abs(quao%population_bond_order(i, j)), energy, n_core + i, " / ", n_core + j
         call logger%info(trim(line))
      end do
      if (printed == 0) call logger%info("    (none above threshold)")

      ! ---- everything else that couples ------------------------------------
      call logger%info("")
      call logger%info("  delocalization")
      call logger%info("     donor             into              "// &
                       "  order    kcal/mol   orbitals")
      printed = 0
      suppressed = 0
      do k = 1, np
         i = pair_i(order(k))
         j = pair_j(order(k))
         if (is_bond(quao, labels, i, j)) cycle
         if (present(kinetic_bond_order)) then
            if (abs(kinetic_bond_order(i, j)) < cutoff) then
               suppressed = suppressed + 1
               cycle
            end if
         end if
         printed = printed + 1
         ! Donor first, decided by occupation rather than by index: the arrow
         ! only means something if the fuller orbital is the one giving.
         if (quao%population_bond_order(i, i) >= quao%population_bond_order(j, j)) then
            call interaction_text(quao, labels, element_symbols, i, left)
            call interaction_text(quao, labels, element_symbols, j, right)
            printed_i = i
            printed_j = j
         else
            call interaction_text(quao, labels, element_symbols, j, left)
            call interaction_text(quao, labels, element_symbols, i, right)
            printed_i = j
            printed_j = i
         end if
         energy = interference(i, j)
         if (present(kinetic_bond_order)) energy = kinetic_bond_order(i, j)
         field_a = left
         field_b = right
         write (line, "(4x,a18,a18,f10.4,f11.2,i7,a,i0)") field_a, field_b, &
            quao%population_bond_order(i, j), energy, &
            n_core + printed_i, " / ", n_core + printed_j
         call logger%info(trim(line))
      end do
      if (printed == 0) call logger%info("    (none above threshold)")
      if (suppressed > 0) then
         write (line, "(a,i0,a,f0.2,a)") "    (", suppressed, " weaker than ", &
            cutoff, " kcal/mol not shown)"
         call logger%info(trim(line))
      end if

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
      write (line, "(a,f12.6,a)") "     interatomic total ", total, " hartree"
      call logger%info(trim(line))

      ! ---- one row per orbital ----------------------------------------------
      call logger%info("")
      call logger%info("  orbitals")
      call logger%info("      orb  atom   type        bonded to  "// &
                       " occupation     %s     %p     %d")
      allocate (occupation(n))
      do i = 1, n
         occupation(i) = quao%population_bond_order(i, i)
      end do
      deallocate (order)
      call sort_by_magnitude(occupation, order)
      do k = 1, n
         i = order(k)
         call atom_text(quao, element_symbols, i, left)
         call partner_text(quao, labels, element_symbols, i, right)
         field_a = left
         field_b = quao_type_name(labels%orbital_type(i), labels%dominant_l(i))
         write (line, "(4x,i5,2x,a7,a12,a10,f11.4,3f7.1)") &
            n_core + i, field_a(1:7), field_b(1:12), right(1:10), occupation(i), &
            100.0_dp*labels%angular_character(0, i), &
            100.0_dp*labels%angular_character(1, i), &
            100.0_dp*labels%angular_character(2, i)
         call logger%info(trim(line))
      end do

      deallocate (pair_i, pair_j, occupation, order)
   end subroutine print_quao_report

   pure function is_bond(quao, labels, i, j) result(bond)
      !! Whether a pair is a covalent bond rather than a delocalization
      !!
      !! Both ends must be bonding orbitals, and each must name the other's
      !! atom as what it is bonded to. A lone pair overlapping a bond fails the
      !! first test; two bonding orbitals on atoms that are not bonded to each
      !! other fail the second.
      type(quao_result_t), intent(in) :: quao
      type(quao_labels_t), intent(in) :: labels
      integer, intent(in) :: i, j
      logical :: bond

      integer :: k

      bond = .false.
      if (quao%atom_of(i) == quao%atom_of(j)) return
      select case (labels%orbital_type(i))
      case (QUAO_TYPE_SIGMA, QUAO_TYPE_PI)
      case default
         return
      end select
      select case (labels%orbital_type(j))
      case (QUAO_TYPE_SIGMA, QUAO_TYPE_PI)
      case default
         return
      end select

      do k = 1, labels%partner_count(i)
         if (quao%atom_of(labels%partner(i, k)) == quao%atom_of(j)) bond = .true.
      end do
      if (.not. bond) return
      bond = .false.
      do k = 1, labels%partner_count(j)
         if (quao%atom_of(labels%partner(j, k)) == quao%atom_of(i)) bond = .true.
      end do
   end function is_bond

   subroutine atom_text(quao, element_symbols, iorb, text)
      !! "C 1" -- the atom an orbital belongs to
      type(quao_result_t), intent(in) :: quao
      character(len=*), intent(in) :: element_symbols(:)
      integer, intent(in) :: iorb
      character(len=*), intent(out) :: text

      write (text, "(a,i0)") trim(adjustl(element_symbols(quao%atom_of(iorb))))//" ", &
         quao%atom_of(iorb)
   end subroutine atom_text

   subroutine partner_text(quao, labels, element_symbols, iorb, text)
      !! "(O 2)" -- what an orbital is bonded to, or blank
      type(quao_result_t), intent(in) :: quao
      type(quao_labels_t), intent(in) :: labels
      character(len=*), intent(in) :: element_symbols(:)
      integer, intent(in) :: iorb
      character(len=*), intent(out) :: text

      character(len=24) :: atom
      integer :: k

      text = ""
      if (labels%partner_count(iorb) == 0) return
      do k = 1, labels%partner_count(iorb)
         call atom_text(quao, element_symbols, labels%partner(iorb, k), atom)
         if (k == 1) then
            text = "("//trim(atom)
         else
            text = trim(text)//","//trim(atom)
         end if
      end do
      text = trim(text)//")"
   end subroutine partner_text

   subroutine interaction_text(quao, labels, element_symbols, iorb, text)
      !! "O 2 PLP" or "C 1-Cl 3" -- one end of a delocalization, said in words
      !!
      !! A lone pair is named by its atom and kind; a bonding orbital by the
      !! bond it belongs to, since naming only the atom would lose which of its
      !! bonds is accepting the density.
      type(quao_result_t), intent(in) :: quao
      type(quao_labels_t), intent(in) :: labels
      character(len=*), intent(in) :: element_symbols(:)
      integer, intent(in) :: iorb
      character(len=*), intent(out) :: text

      character(len=24) :: atom, partner
      integer :: kind

      call atom_text(quao, element_symbols, iorb, atom)
      kind = labels%orbital_type(iorb)
      if ((kind == QUAO_TYPE_SIGMA .or. kind == QUAO_TYPE_PI) .and. &
          labels%partner_count(iorb) > 0) then
         call atom_text(quao, element_symbols, labels%partner(iorb, 1), partner)
         text = trim(atom)//"-"//trim(partner)//" "// &
                trim(quao_type_name(kind, labels%dominant_l(iorb)))
      else
         text = trim(atom)//" "//trim(quao_type_name(kind, labels%dominant_l(iorb)))
      end if
   end subroutine interaction_text

   subroutine start_row(line, pos, first, second)
      !! Blank the line and lay down the two leading numbers
      ! TODO(mqc): neither `start_row` nor `append_orbital` is called from
      ! anywhere -- `print_quao_report` writes its rows with whole-line format
      ! statements instead -- so both, and the five `WIDTH_*` parameters that
      ! exist only for them, are dead.
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
      !! The cursor is carried rather than recovered with `len_trim`: a
      !! fixed-width field ending in blanks cannot be found again that way, and
      !! the columns after it would drift.
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
      ! orbital has three or more partners, rather than truncating the list.
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

end module mqc_czt_quao_report
