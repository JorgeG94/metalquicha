module mqc_cgto
   use pic_types, only: dp
   implicit none
   private

   public :: cgto_type, atomic_basis_type, molecular_basis_type

   type :: cgto_type
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
      character(len=:), allocatable :: element
      !! element symbol
      type(cgto_type), allocatable :: shells(:)
      !! array of contracted shells
      integer :: nshells
      !! number of shells in type
   contains
      procedure :: allocate_shells => allocate_basis_shells
      procedure :: destroy => atomic_basis_destroy
      procedure :: num_basis_functions => atomic_basis_num_basis_functions
   end type atomic_basis_type

   type :: molecular_basis_type
      type(atomic_basis_type), allocatable :: elements(:)
      !! array of atomic basis types
      integer :: nelements
      !! total number of atoms/elements in a molecule
   contains
      procedure :: allocate_elements => basis_set_allocate_elements
      procedure :: destroy => basis_set_destroy
      procedure :: num_basis_functions => molecular_basis_num_basis_functions
   end type molecular_basis_type

contains

   pure subroutine cgto_allocate_arrays(self, nfunc)
      class(cgto_type), intent(inout) :: self
      integer, intent(in) :: nfunc

      self%nfunc = nfunc
      allocate (self%exponents(nfunc))
      allocate (self%coefficients(nfunc))
   end subroutine cgto_allocate_arrays

   pure subroutine cgto_destroy(self)
      class(cgto_type), intent(inout) :: self

      if (allocated(self%exponents)) deallocate (self%exponents)
      if (allocated(self%coefficients)) deallocate (self%coefficients)
      self%nfunc = 0
      self%ang_mom = 0
   end subroutine cgto_destroy

   pure subroutine allocate_basis_shells(self, nshells)
      class(atomic_basis_type), intent(inout) :: self
      integer, intent(in) :: nshells

      self%nshells = nshells
      allocate (self%shells(nshells))
   end subroutine allocate_basis_shells

   pure subroutine atomic_basis_destroy(self)
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
   end subroutine atomic_basis_destroy

   pure subroutine basis_set_allocate_elements(self, nelements)
      class(molecular_basis_type), intent(inout) :: self
      integer, intent(in) :: nelements

      self%nelements = nelements
      allocate (self%elements(nelements))

   end subroutine basis_set_allocate_elements

   pure subroutine basis_set_destroy(self)
      class(molecular_basis_type), intent(inout) :: self
      integer :: i

      if (allocated(self%elements)) then
         do i = 1, self%nelements
            call self%elements(i)%destroy()
         end do
         deallocate (self%elements)
      end if

      self%nelements = 0
   end subroutine basis_set_destroy

!> Get number of basis functions in a shell (Cartesian)
   pure function cgto_num_basis_functions(self) result(nbf)
      class(cgto_type), intent(in) :: self
      integer :: nbf

      ! Cartesian: (ang_mom+1)*(ang_mom+2)/2
      nbf = (self%ang_mom + 1)*(self%ang_mom + 2)/2
   end function cgto_num_basis_functions

!> Get total number of basis functions for an atom
   pure function atomic_basis_num_basis_functions(self) result(nbf)
      class(atomic_basis_type), intent(in) :: self
      integer :: nbf
      integer :: ishell

      nbf = 0
      do ishell = 1, self%nshells
         nbf = nbf + self%shells(ishell)%num_basis_functions()
      end do
   end function atomic_basis_num_basis_functions

!> Get total number of basis functions for the molecule
   pure function molecular_basis_num_basis_functions(self) result(nbf)
      class(molecular_basis_type), intent(in) :: self
      integer :: nbf
      integer :: iatom

      nbf = 0
      do iatom = 1, self%nelements
         nbf = nbf + self%elements(iatom)%num_basis_functions()
      end do
   end function molecular_basis_num_basis_functions

end module mqc_cgto
