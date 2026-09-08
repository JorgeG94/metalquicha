! LFortran 0.64.0: referencing a variable from the enclosing scope inside an
! OpenMP region segfaults the compiler.
!
!   $ lfortran --openmp -c 02_outer_variable.f90
!   Segmentation fault (core dumped)     [exit 139]
!
! A *read* is enough -- no assignment, no critical, no reduction. Naming the
! variable explicitly (`shared(total)`) does not change it, and neither does
! `parallel do` in place of `parallel`. Remove the `print` and the same file
! compiles, which is the whole of the difference.
program outer_variable
   implicit none
   integer :: total
   total = 7
   !$omp parallel
   print *, total
   !$omp end parallel
end program outer_variable
