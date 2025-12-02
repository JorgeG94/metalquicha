!! contains MPI tags used in the MQC parallel implementation
module mqc_mpi_tags
   !! Module defining MPI communication tags for clarity and maintainability
   use pic_types, only: default_int
   implicit none
   private

   ! Local worker communication tags (shared memory within a node)
   integer(default_int), parameter, public :: TAG_WORKER_REQUEST = 200
      !! Worker requests work from node coordinator
   integer(default_int), parameter, public :: TAG_WORKER_FRAGMENT = 201
      !! Coordinator sends fragment data to worker
   integer(default_int), parameter, public :: TAG_WORKER_FINISH = 202
      !! Coordinator signals worker to finish
   integer(default_int), parameter, public :: TAG_WORKER_SCALAR_RESULT = 203
      !! Worker sends scalar results back to coordinator
   integer(default_int), parameter, public :: TAG_WORKER_MATRIX_RESULT = 204
      !! Worker sends matrix results back to coordinator

   ! Remote node communication tags (between nodes via world communicator)
   integer(default_int), parameter, public :: TAG_NODE_REQUEST = 300
      !! Node coordinator requests work from global coordinator
   integer(default_int), parameter, public :: TAG_NODE_FRAGMENT = 301
      !! Global coordinator sends fragment data to node coordinator
   integer(default_int), parameter, public :: TAG_NODE_FINISH = 302
      !! Global coordinator signals node coordinator to finish
   integer(default_int), parameter, public :: TAG_NODE_SCALAR_RESULT = 303
      !! Node coordinator sends results (fragment_idx + scalar) to global coordinator
   integer(default_int), parameter, public :: TAG_NODE_MATRIX_RESULT = 304
      !! Node coordinator sends matrix results to global coordinator

end module mqc_mpi_tags
