! LFortran 0.64.0: the `default` clause is rejected on a worksharing directive.
!
!   $ lfortran --openmp -c 03_default_clause.f90
!   semantic error: The clause default is not supported for parallel sections
!
! This is the first error metalquicha's own sources hit (src/methods/ci/mqc_ci.f90).
subroutine default_clause(a, n)
   implicit none
   integer, intent(in) :: n
   real, intent(inout) :: a(n)
   integer :: i
   !$omp parallel do default(shared) private(i) schedule(static)
   do i = 1, n
      a(i) = a(i)*2.0
   end do
   !$omp end parallel do
end subroutine default_clause
