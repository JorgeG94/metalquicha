!! MPI communicators container for the application
module mqc_mpi_comms
   !! Container type for MPI communicators - extensible for future parallelism patterns
   use pic_mpi_lib, only: comm_t
   implicit none
   private

   public :: mpi_comms_t

   type :: mpi_comms_t
      !! Container for MPI communicators
      !!
      !! This type bundles all MPI communicators needed by the application.
      !! Currently supports world and node communicators, but can be extended
      !! to support team-based parallelism patterns.
      type(comm_t) :: world_comm       !! Global MPI communicator
      type(comm_t) :: node_comm        !! Node-local communicator
      ! Future extensions:
      ! type(comm_t) :: team_comm           !! Intra-team communicator
      ! type(comm_t) :: team_worker_comm    !! Team workers (excluding leader)
      ! type(comm_t) :: inter_team_comm     !! Cross-team communicator
      ! type(comm_t) :: team_leaders_comm   !! Team leaders only
   end type mpi_comms_t

end module mqc_mpi_comms
