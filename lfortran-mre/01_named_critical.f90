! LFortran 0.64.0: a named OpenMP critical region is rejected.
!
!   $ lfortran --openmp -c 01_named_critical.f90
!   semantic error: The clause  is not supported for parallel sections
!    --> 01_named_critical.f90:14:4
!      |
!   14 |    !$omp critical (accumulate)
!
! gfortran -fopenmp compiles and runs this. Dropping the name `(accumulate)`
! removes THIS error but then hits MRE 02 instead, so the two are independent.
subroutine named_critical(total)
   implicit none
   integer, intent(inout) :: total
   !$omp critical (accumulate)
   total = total + 1
   !$omp end critical (accumulate)
end subroutine named_critical
