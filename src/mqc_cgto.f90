!! Data structures for cartesian contracted Gaussian type orbitals (CGTOs)
module mqc_cgto
   !! Defines data structures for cartesian contracted Gaussian type orbitals (CGTOs)
   use pic_types, only: dp
   implicit none
   private

   public :: cgto_type, atomic_basis_type, molecular_basis_type
   public :: ANGULAR_FORM_UNSET, ANGULAR_FORM_SPHERICAL, ANGULAR_FORM_CARTESIAN

   ! Which angular form the basis is written in. Cartesian and spherical are
   ! the same thing for s and p -- one function and three, either way -- so a
   ! set that never goes above p says nothing about the convention, and that is
   ! a third state rather than a default. Collapsing it onto "spherical" too
   ! early is what would make hydrogen look like it disagreed with the oxygen
   ! it is bonded to in 6-31G*, where only the d shell is marked Cartesian.
   integer, parameter :: ANGULAR_FORM_UNSET = 0
      !! No shell above p, so the file says nothing either way
   integer, parameter :: ANGULAR_FORM_SPHERICAL = 1
      !! 5d, 7f: `gto_spherical` in Basis Set Exchange JSON
   integer, parameter :: ANGULAR_FORM_CARTESIAN = 2
      !! 6d, 10f: `gto_cartesian` in Basis Set Exchange JSON

   type :: cgto_type
      !! Contracted Gaussian type orbital (CGTO) data structure
      integer :: ang_mom
        !! Angular momentum quantum number (0=s, 1=p, 2=d, etc.)
      integer :: nfunc
        !! Number of primitive Gaussians in the contraction
      real(dp), allocatable :: exponents(:)
        !! Exponents (alpha values)
      real(dp), allocatable :: coefficients(:)
        !! Contraction coefficients
   contains
      procedure :: allocate_arrays => cgto_allocate_arrays
      procedure :: destroy => cgto_destroy
      procedure :: num_basis_functions => cgto_num_basis_functions
   end type cgto_type

   type :: atomic_basis_type
      !! Atomic basis set data structure
      character(len=:), allocatable :: element
      !! element symbol
      type(cgto_type), allocatable :: shells(:)
      !! array of contracted shells
      integer :: nshells = 0
      !! number of shells in type. Zero until `allocate_shells` runs, so an
      !! atom the basis file never covered reads as empty rather than as
      !! whatever was on the stack.
      integer :: angular_form = ANGULAR_FORM_UNSET
      !! Cartesian or spherical, from this element's `function_type` entries
   contains
      procedure :: allocate_shells => allocate_basis_shells
      procedure :: destroy => atomic_basis_destroy
      procedure :: num_basis_functions => atomic_basis_num_basis_functions
   end type atomic_basis_type

   type :: molecular_basis_type
      !! Molecular basis set data structure (assembled basis)
      type(atomic_basis_type), allocatable :: elements(:)
      !! array of atomic basis types
      integer :: nelements = 0
      !! total number of atoms/elements in a molecule
      integer :: angular_form = ANGULAR_FORM_UNSET
      !! Cartesian or spherical, agreed across every atom in the molecule.
      !! The integral layer picks one set of entry points from this, so it has
      !! to be one answer for the whole basis; the reader refuses a basis whose
      !! atoms disagree rather than choosing for them.
   contains
      procedure :: allocate_elements => basis_set_allocate_elements
      procedure :: destroy => basis_set_destroy
      procedure :: num_basis_functions => molecular_basis_num_basis_functions
      procedure :: is_cartesian => molecular_basis_is_cartesian
   end type molecular_basis_type

contains

   pure subroutine cgto_allocate_arrays(self, nfunc)
      !! Allocate arrays for exponents and coefficients in a CGTO
      class(cgto_type), intent(inout) :: self
      integer, intent(in) :: nfunc

      self%nfunc = nfunc
      allocate (self%exponents(nfunc))
      allocate (self%coefficients(nfunc))
   end subroutine cgto_allocate_arrays

   pure subroutine cgto_destroy(self)
      !! Clean up allocated memory in a CGTO
      class(cgto_type), intent(inout) :: self

      if (allocated(self%exponents)) deallocate (self%exponents)
      if (allocated(self%coefficients)) deallocate (self%coefficients)
      self%nfunc = 0
      self%ang_mom = 0
   end subroutine cgto_destroy

   pure subroutine allocate_basis_shells(self, nshells)
      !! Allocate array of shells in an atomic basis
      class(atomic_basis_type), intent(inout) :: self
      integer, intent(in) :: nshells

      self%nshells = nshells
      allocate (self%shells(nshells))
   end subroutine allocate_basis_shells

   pure subroutine atomic_basis_destroy(self)
      !! Clean up allocated memory in an atomic basis
      class(atomic_basis_type), intent(inout) :: self
      integer :: i

      if (allocated(self%shells)) then
         do i = 1, self%nshells
            call self%shells(i)%destroy()
         end do
         deallocate (self%shells)
      end if
      if (allocated(self%element)) deallocate (self%element)
      self%nshells = 0
      self%angular_form = ANGULAR_FORM_UNSET
   end subroutine atomic_basis_destroy

   pure subroutine basis_set_allocate_elements(self, nelements)
      !! Allocate array of atomic basis elements in a molecular basis set
      class(molecular_basis_type), intent(inout) :: self
      integer, intent(in) :: nelements

      self%nelements = nelements
      allocate (self%elements(nelements))

   end subroutine basis_set_allocate_elements

   pure subroutine basis_set_destroy(self)
      !! Clean up allocated memory in a molecular basis set
      class(molecular_basis_type), intent(inout) :: self
      integer :: i

      if (allocated(self%elements)) then
         do i = 1, self%nelements
            call self%elements(i)%destroy()
         end do
         deallocate (self%elements)
      end if

      self%nelements = 0
      self%angular_form = ANGULAR_FORM_UNSET
   end subroutine basis_set_destroy

   pure function molecular_basis_is_cartesian(self) result(cartesian)
      !! Whether the integrals over this basis must be built in the Cartesian form
      !!
      !! `ANGULAR_FORM_UNSET` answers false, and correctly: a basis with no
      !! shell above p has the same integrals either way, so the spherical
      !! entry points are as right as the Cartesian ones and are what everything
      !! else already uses.
      class(molecular_basis_type), intent(in) :: self
      logical :: cartesian

      cartesian = self%angular_form == ANGULAR_FORM_CARTESIAN
   end function molecular_basis_is_cartesian

   pure function cgto_num_basis_functions(self) result(nbf)
      !! Get number of basis functions in a shell (Cartesian)
      class(cgto_type), intent(in) :: self
      integer :: nbf

      ! Cartesian: (ang_mom+1)*(ang_mom+2)/2
      nbf = (self%ang_mom + 1)*(self%ang_mom + 2)/2
   end function cgto_num_basis_functions

   pure function atomic_basis_num_basis_functions(self) result(nbf)
      !! Get total number of basis functions for an atom
      class(atomic_basis_type), intent(in) :: self
      integer :: nbf
      integer :: ishell

      nbf = 0
      do ishell = 1, self%nshells
         nbf = nbf + self%shells(ishell)%num_basis_functions()
      end do
   end function atomic_basis_num_basis_functions

   pure function molecular_basis_num_basis_functions(self) result(nbf)
      !! Get total number of basis functions for the molecule
      class(molecular_basis_type), intent(in) :: self
      integer :: nbf
      integer :: iatom

      nbf = 0
      do iatom = 1, self%nelements
         nbf = nbf + self%elements(iatom)%num_basis_functions()
      end do
   end function molecular_basis_num_basis_functions

end module mqc_cgto
