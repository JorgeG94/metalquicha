# Met'al q'uicha (metalquicha)

Met'al q'uicha, which I'll just write as metalquicha is a sample quantum chemistry backend 
with focus on using the [pic](https://github.com/JorgeG94/pic) library and its derivatives:
[pic-mpi](https://github.com/JorgeG94/pic-mpi) and [pic-blas](https://github.com/JorgeG94/pic-blas)
which are Fortran based implementations of commonly used routines such as sorting algorithms, 
array handling, strings, loggers, timers, etc. While pic-mpi and pic-blas provide modern 
Fortran interfaces to MPI and BLAS implementations in a portable way. Specifically, the MPI
library lets the user switch betwen the `mpi` and `mpi_f08` modules with ease. 

Metalquicha implements a naive backend for unfragmented and fragmented quantum chemistry
calculations. Currently, metalquicha uses [tblite](https://github.com/tblite/tblite) as 
its chemistry engine which performs energy calculations. 

The main purpose of this package is to showcase the ability of being able to write 
simple, powerful Fortran based programs that are able to access massively parallel 
ecosystems with ease. 

For example, the entire coordination scheme is just:

```fortran
   ! Run the test based on role
   if (world_comm%leader() .and. node_comm%leader()) then
      ! Global coordinator (rank 0, node leader on node 0)
      print *, "Rank 0: Acting as global coordinator"
      call global_coordinator(world_comm, node_comm, total_fragments, polymers, max_level, &
                                    node_leader_ranks, num_nodes, matrix_size)
   else if (node_comm%leader()) then
      ! Node coordinator (node leader on other nodes)
      print '(a,i0,a)', "Rank ", world_rank, ": Acting as node coordinator"
      call node_coordinator(world_comm, node_comm, max_level, matrix_size)
   else
      ! Worker
      print '(a,i0,a)', "Rank ", world_rank, ": Acting as worker"
      call node_worker(world_comm, node_comm, max_level, sys_geom)
   end if
```

