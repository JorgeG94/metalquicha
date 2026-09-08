! LFortran 0.64.0: `default(none)` is rejected, same as `default(shared)`.
!
!   $ lfortran --openmp -c 05_default_none.f90
!   semantic error: The clause default is not supported for parallel sections
!
! Kept separate from MRE 03 because `default(shared)` is removable -- shared is
! already the default for `parallel`, so dropping it changes nothing. Dropping
! `default(none)` is NOT neutral: it turns off exactly the compile-time check
! that every variable was classified on purpose.
subroutine default_none(a, n)
   implicit none
   integer, intent(in) :: n
   real, intent(inout) :: a(n)
   integer :: i
   !$omp parallel do default(none) shared(a, n) private(i)
   do i = 1, n
      a(i) = a(i)*2.0
   end do
   !$omp end parallel do
end subroutine default_none
