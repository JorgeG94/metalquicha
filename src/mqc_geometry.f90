module mqc_geometry
   use pic_types, only: dp
   implicit none
   private

   public :: geometry_type

   ! Parameters
   integer, parameter :: MAX_ELEMENT_SYMBOL_LEN = 4

   type :: geometry_type
      integer :: natoms
      character(len=:), allocatable :: elements(:)
      real(dp), allocatable :: coords(:, :)  ! coords(3, natoms)
      character(len=:), allocatable :: comment
   contains
      procedure :: destroy => geometry_destroy
   end type geometry_type

contains

   subroutine geometry_destroy(this)
      class(geometry_type), intent(inout) :: this
      if (allocated(this%elements)) deallocate (this%elements)
      if (allocated(this%coords)) deallocate (this%coords)
      if (allocated(this%comment)) deallocate (this%comment)
      this%natoms = 0
   end subroutine geometry_destroy

end module mqc_geometry
