! LFortran 0.64.0: an OpenMP directive continued onto a second line is rejected.
!
!   $ lfortran --openmp -c 04_continued_directive.f90
!   semantic error: The clause & is not supported for parallel sections
!
! The `&` is being read as a clause name rather than as the line continuation
! it is, so the second `!$omp` line is never joined to the first. This is the
! blocker with the widest reach in metalquicha: 359 continued directives across
! 23 files, which is how a directive with more than two or three clauses is
! written everywhere in this codebase.
subroutine continued_directive(a, n)
   implicit none
   integer, intent(in) :: n
   real, intent(inout) :: a(n)
   integer :: i, j
   !$omp parallel do schedule(static) &
   !$omp   private(i, j)
   do i = 1, n
      j = i
      a(j) = a(j)*2.0
   end do
   !$omp end parallel do
end subroutine continued_directive
