! Control for MRE 02: identical region, no enclosing-scope variable. Compiles.
program control
   implicit none
   !$omp parallel
   print *, "hi"
   !$omp end parallel
end program control
