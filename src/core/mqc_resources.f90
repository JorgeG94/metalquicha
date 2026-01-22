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

      ! Hardware info
      integer :: num_threads = 1       !! OpenMP threads
      integer :: num_gpus = 0          !! Available GPUs
      logical :: use_gpu = .false.     !! GPU acceleration flag

      ! Future extensions:
      ! type(blas_provider_t), pointer :: blas => null()
      ! type(basis_reader_t), pointer :: basis_reader => null()
      ! type(memory_pool_t), pointer :: memory => null()
   contains
      procedure :: init => resources_init
      procedure :: finalize => resources_finalize
   end type resources_t

contains

   subroutine resources_init(this)
      !! Initialize resources with default values
      class(resources_t), intent(inout) :: this

      this%num_threads = 1
      this%num_gpus = 0
      this%use_gpu = .false.
   end subroutine resources_init

   subroutine resources_finalize(this)
      !! Clean up resources
      class(resources_t), intent(inout) :: this

      ! Currently nothing to clean up, but this provides
      ! a hook for future resource cleanup (GPU handles, memory pools, etc.)
      this%num_threads = 1
      this%num_gpus = 0
      this%use_gpu = .false.
   end subroutine resources_finalize

end module mqc_resources
