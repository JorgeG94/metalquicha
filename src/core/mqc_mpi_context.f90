!! MPI context type for error handling across the codebase
module mqc_mpi_context
   use pic_mpi_lib, only: comm_t, abort_comm
   use pic_logger, only: logger => global_logger

   implicit none
   private

   public :: mpi_context_t

   type :: mpi_context_t
      !! Bundles MPI communicators for passing through call chains
      !! Provides convenient error abort functionality
      type(comm_t) :: world_comm  !! Global MPI communicator
      type(comm_t) :: node_comm   !! Node-local MPI communicator
   contains
      procedure :: abort => mpi_context_abort
   end type mpi_context_t

contains

   subroutine mpi_context_abort(this, message, error_code)
      !! Abort all MPI processes with an error message
      class(mpi_context_t), intent(in) :: this
      character(len=*), intent(in) :: message
      integer, intent(in), optional :: error_code

      integer :: code

      code = 1
      if (present(error_code)) code = error_code

      call logger%error(message)
      call abort_comm(this%world_comm, code)
   end subroutine mpi_context_abort

end module mqc_mpi_context
