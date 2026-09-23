!! The quasi-atomic bonding tables as data rather than as printed lines
module mqc_quao_rows
   !! What a quasi-atomic bonding analysis found in one calculation, kept so it
   !! can travel with the result: over MPI to the coordinator, and into the
   !! JSON output. Backend-neutral, because `calculation_result_t` carries it.
   !!
   !! Two levels, and they are not equally robust. The **orbital pairs** are
   !! the rows of the report's bonds and delocalization tables, and depend on
   !! the orientation of the orbitals within each atom. The **atom pairs** are
   !! sums over every orbital on one atom against every orbital on another,
   !! which no rotation within an atom can change (Paper I eq (5.4)), so they
   !! hold even where the orientation is slow to settle.
   use pic_types, only: dp
   implicit none
   private

   public :: quao_rows_t
   public :: quao_type_name
   public :: QUAO_ROW_BOND, QUAO_ROW_DELOCALIZATION
   public :: QUAO_ROW_INTEGERS, QUAO_ROW_REALS
   public :: QUAO_TYPE_LPMOD, QUAO_TYPE_NONE, QUAO_TYPE_LONE_PAIR, QUAO_TYPE_RADICAL
   public :: QUAO_TYPE_SIGMA, QUAO_TYPE_PI, QUAO_TYPE_RDLP, QUAO_TYPE_RDNV
   public :: QUAO_TYPE_NVMOD, QUAO_TYPE_NV

   ! Orbital types, numbered as GAMESS numbers them so the two can be compared
   ! without a translation table.
   integer, parameter :: QUAO_TYPE_LPMOD = -1
      !! Lone-pair occupancy, but engaged in a strong bond. A lone pair that is
      !! not quite one.
   integer, parameter :: QUAO_TYPE_NONE = 0
      !! Unclassified. Reached only on a transition metal, where the sigma/pi
      !! test is skipped.
   integer, parameter :: QUAO_TYPE_LONE_PAIR = 1
   integer, parameter :: QUAO_TYPE_RADICAL = 2
   integer, parameter :: QUAO_TYPE_SIGMA = 3
   integer, parameter :: QUAO_TYPE_PI = 4
   integer, parameter :: QUAO_TYPE_RDLP = 6      !! Reduced lone pair
   integer, parameter :: QUAO_TYPE_RDNV = 7      !! Reduced, and nearly empty
   integer, parameter :: QUAO_TYPE_NVMOD = 8     !! Empty, but in a strong bond
   integer, parameter :: QUAO_TYPE_NV = 9        !! Empty

   integer, parameter :: QUAO_ROW_BOND = 1
      !! Both ends bonding orbitals naming each other's atom as partner
   integer, parameter :: QUAO_ROW_DELOCALIZATION = 2
      !! Anything else above threshold; end 1 is the donor

   integer, parameter :: QUAO_ROW_INTEGERS = 11
      !! Integers per row in `pack_integers`
   integer, parameter :: QUAO_ROW_REALS = 4
      !! Reals per row in `pack_reals`

   type :: quao_rows_t
      !! Orbital-pair rows and atom-pair sums from one bonding analysis
      integer :: n = 0
         !! Orbital-pair rows
      integer, allocatable :: kind(:)
         !! (n), `QUAO_ROW_BOND` or `QUAO_ROW_DELOCALIZATION`
      integer, allocatable :: orbital(:, :)
         !! (2, n), numbered as the molecular orbitals are, core included
      integer, allocatable :: atom(:, :)
         !! (2, n), 1-based atom of the calculation the analysis ran on
      integer, allocatable :: orbital_type(:, :)
         !! (2, n), one of the `QUAO_TYPE_*` codes
      integer, allocatable :: dominant_l(:, :)
         !! (2, n), 0 to 3
      integer, allocatable :: partner_atom(:, :)
         !! (2, n), 1-based atom of the orbital's first strong-bond partner,
         !! 0 when it has none
      real(dp), allocatable :: occupation(:, :)
         !! (2, n), electrons
      real(dp), allocatable :: bond_order(:)
         !! (n), the population bond order between the two orbitals
      real(dp), allocatable :: kinetic_bond_order(:)
         !! (n), kcal/mol, negative for a bonding interaction
      real(dp), allocatable :: atom_bond_index(:, :)
         !! (n_atoms, n_atoms), sum of squared bond orders over every orbital
         !! pair between two atoms; zero on the diagonal
      real(dp), allocatable :: atom_kinetic_bond_order(:, :)
         !! (n_atoms, n_atoms), kcal/mol, sum of kinetic bond orders over every
         !! orbital pair between two atoms; zero on the diagonal
      real(dp) :: threshold = 0.0_dp
         !! kcal/mol; delocalization rows weaker than this were not kept
      logical :: orientation_stalled = .false.
         !! The orientation stopped at its sweep limit; the orbital pairs are
         !! resolved only to the angle it was still turning by
   contains
      procedure :: destroy => rows_destroy
      procedure :: pack_integers => rows_pack_integers
      procedure :: pack_reals => rows_pack_reals
      procedure :: unpack => rows_unpack
   end type quao_rows_t

contains

   pure function quao_type_name(orbital_type, dominant_l) result(name)
      !! The label a report prints for an orbital type
      !!
      !! `dominant_l` only matters for the lone-pair types, where it is the
      !! difference between an s lone pair and a p one. Every type has a name,
      !! including the ones GAMESS leaves blank.
      integer, intent(in) :: orbital_type
      integer, intent(in) :: dominant_l
      character(len=10) :: name

      ! Lower case and spelled out rather than GAMESS's `SIGMA`, `PLP`,
      ! `NVMOD`, which are its internal codes: a blank column cannot be told
      ! from an orbital that genuinely has no type.
      select case (orbital_type)
      case (QUAO_TYPE_LPMOD)
         name = shell_letter(dominant_l)//"-lone/bnd"
      case (QUAO_TYPE_LONE_PAIR)
         name = shell_letter(dominant_l)//"-lone"
      case (QUAO_TYPE_RADICAL)
         name = "radical"
      case (QUAO_TYPE_SIGMA)
         name = "sigma"
      case (QUAO_TYPE_PI)
         name = "pi"
      case (QUAO_TYPE_RDLP)
         name = "part-lone"
      case (QUAO_TYPE_RDNV)
         name = "part-empty"
      case (QUAO_TYPE_NVMOD)
         name = "empty/bnd"
      case (QUAO_TYPE_NV)
         name = "empty"
      case default
         name = "unclassed"
      end select
   end function quao_type_name

   pure function shell_letter(l) result(letter)
      integer, intent(in) :: l
      character(len=1) :: letter

      select case (l)
      case (0)
         letter = "s"
      case (1)
         letter = "p"
      case (2)
         letter = "d"
      case (3)
         letter = "f"
      case default
         letter = "?"
      end select
   end function shell_letter

   subroutine rows_destroy(this)
      !! Release every array and return to the empty state
      class(quao_rows_t), intent(inout) :: this

      if (allocated(this%kind)) deallocate (this%kind)
      if (allocated(this%orbital)) deallocate (this%orbital)
      if (allocated(this%atom)) deallocate (this%atom)
      if (allocated(this%orbital_type)) deallocate (this%orbital_type)
      if (allocated(this%dominant_l)) deallocate (this%dominant_l)
      if (allocated(this%partner_atom)) deallocate (this%partner_atom)
      if (allocated(this%occupation)) deallocate (this%occupation)
      if (allocated(this%bond_order)) deallocate (this%bond_order)
      if (allocated(this%kinetic_bond_order)) deallocate (this%kinetic_bond_order)
      if (allocated(this%atom_bond_index)) deallocate (this%atom_bond_index)
      if (allocated(this%atom_kinetic_bond_order)) deallocate (this%atom_kinetic_bond_order)
      this%n = 0
      this%threshold = 0.0_dp
      this%orientation_stalled = .false.
   end subroutine rows_destroy

   function rows_pack_integers(this) result(flat)
      !! Every integer column, `QUAO_ROW_INTEGERS` per row, for one MPI message
      class(quao_rows_t), intent(in) :: this
      integer, allocatable :: flat(:)

      integer :: k, base

      allocate (flat(QUAO_ROW_INTEGERS*this%n))
      do k = 1, this%n
         base = QUAO_ROW_INTEGERS*(k - 1)
         flat(base + 1) = this%kind(k)
         flat(base + 2:base + 3) = this%orbital(:, k)
         flat(base + 4:base + 5) = this%atom(:, k)
         flat(base + 6:base + 7) = this%orbital_type(:, k)
         flat(base + 8:base + 9) = this%dominant_l(:, k)
         flat(base + 10:base + 11) = this%partner_atom(:, k)
      end do
   end function rows_pack_integers

   function rows_pack_reals(this) result(flat)
      !! Every per-row real column, `QUAO_ROW_REALS` per row
      !!
      !! The atom-pair matrices are not included; they go as 2-D arrays.
      class(quao_rows_t), intent(in) :: this
      real(dp), allocatable :: flat(:)

      integer :: k, base

      allocate (flat(QUAO_ROW_REALS*this%n))
      do k = 1, this%n
         base = QUAO_ROW_REALS*(k - 1)
         flat(base + 1:base + 2) = this%occupation(:, k)
         flat(base + 3) = this%bond_order(k)
         flat(base + 4) = this%kinetic_bond_order(k)
      end do
   end function rows_pack_reals

   subroutine rows_unpack(this, integers, reals)
      !! The inverse of `pack_integers` and `pack_reals`
      !!
      !! Sizes the row arrays from `integers`; leaves the atom-pair matrices,
      !! the threshold and the stall flag alone.
      class(quao_rows_t), intent(inout) :: this
      integer, intent(in) :: integers(:)
      real(dp), intent(in) :: reals(:)

      integer :: k, base, n

      n = size(integers)/QUAO_ROW_INTEGERS
      this%n = n
      allocate (this%kind(n), this%orbital(2, n), this%atom(2, n), &
                this%orbital_type(2, n), this%dominant_l(2, n), this%partner_atom(2, n), &
                this%occupation(2, n), this%bond_order(n), this%kinetic_bond_order(n))
      do k = 1, n
         base = QUAO_ROW_INTEGERS*(k - 1)
         this%kind(k) = integers(base + 1)
         this%orbital(:, k) = integers(base + 2:base + 3)
         this%atom(:, k) = integers(base + 4:base + 5)
         this%orbital_type(:, k) = integers(base + 6:base + 7)
         this%dominant_l(:, k) = integers(base + 8:base + 9)
         this%partner_atom(:, k) = integers(base + 10:base + 11)
         base = QUAO_ROW_REALS*(k - 1)
         this%occupation(:, k) = reals(base + 1:base + 2)
         this%bond_order(k) = reals(base + 3)
         this%kinetic_bond_order(k) = reals(base + 4)
      end do
   end subroutine rows_unpack

end module mqc_quao_rows
