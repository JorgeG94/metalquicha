! LFortran 0.64.0: a derived-type PARAMETER whose value is a structure
! constructor crashes the compiler when a variable of a type that defaults to
! it is declared.
!
!   $ lfortran -c 07_struct_constructor_ice.f90
!   Internal Compiler Error: Unhandled exception
!   LCompilersException: visit_StructConstructor() not implemented
!
! The module itself compiles. The crash happens at the *use* site, on the
! `type(outer_t) :: v` declaration -- there is no structure constructor in the
! program unit that dies. That is exactly how it presents in this project:
! `libmetalquicha.a` archives (all 314 objects, cenzontle included) and then
! app/main.f90 dies on `type(resources_t) :: resources`, which reaches
! `MPI_COMM_NULL = MPI_Comm(0)` in pic-mpi's serial backend three types down.
!
! Initializing the component with an inline `inner_t(0)` instead of the named
! constant compiles. gfortran accepts both.
module ice
   implicit none

   type :: inner_t
      integer :: handle = 0
   end type inner_t

   type(inner_t), parameter :: INNER_NULL = inner_t(0)

   type :: outer_t
      type(inner_t) :: part = INNER_NULL
   end type outer_t
end module ice

program use_it
   use ice, only: outer_t
   implicit none
   type(outer_t) :: v
   print *, v%part%handle
end program use_it
