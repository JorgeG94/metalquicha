! LFortran 0.64.0: a derived-type component whose name matches a module
! parameter cannot be default-initialized from that parameter.
!
!   $ lfortran -c 06_component_name_shadows_parameter.f90
!   semantic error: Initialization of `grid_level` must reduce to a compile
!   time constant.
!
! Fortran is case-insensitive, so `GRID_LEVEL` and `grid_level` are one name.
! LFortran appears to resolve the initializer to the component it is declaring
! rather than to the module parameter, so the value stops being constant.
! Rename the component (see the `ok_` type below) and it compiles, which is the
! whole of the difference. gfortran accepts both.
module shadow
   implicit none
   integer, parameter :: GRID_LEVEL = 1

   type :: bad_t
      integer :: grid_level = GRID_LEVEL
   end type bad_t

   type :: ok_t
      integer :: level = GRID_LEVEL
   end type ok_t
end module shadow
