!! Resources container for calculation methods
module mqc_resources
   !! Container type for calculation resources - extensible for future needs
   use mqc_mpi_comms, only: mpi_comms_t
   implicit none
   private

   public :: resources_t

   type :: resources_t
      !! Container for calculation resources
      !!
      !! This type bundles all resources needed by calculation methods.
      !! Currently holds MPI communicators, but can be extended to include
      !! BLAS providers, basis readers, and other shared resources.
      type(mpi_comms_t) :: mpi_comms   !! MPI communicators
      ! Future extensions:
      ! type(blas_provider_t), pointer :: blas => null()
      ! type(basis_reader_t), pointer :: basis_reader => null()
      ! type(memory_pool_t), pointer :: memory => null()
   end type resources_t

end module mqc_resources
