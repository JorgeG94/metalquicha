module mqc_result_types
   use pic_types, only: dp
   implicit none
   private

   public :: calculation_result_t

   !> Container for calculation results
   type :: calculation_result_t
      real(dp) :: energy = 0.0_dp
         !! Total energy
      real(dp), allocatable :: gradient(:, :)
         !! Energy gradient (3, natoms)
      real(dp), allocatable :: hessian(:, :)
         !! Energy hessian (future)
      real(dp), allocatable :: dipole(:)
         !! Dipole moment (3)

      ! Flags for what's been computed
      logical :: has_energy = .false.
      logical :: has_gradient = .false.
      logical :: has_hessian = .false.
      logical :: has_dipole = .false.
   contains
      procedure :: destroy => result_destroy
      procedure :: reset => result_reset
   end type calculation_result_t

contains

   subroutine result_destroy(this)
      class(calculation_result_t), intent(inout) :: this
      if (allocated(this%gradient)) deallocate (this%gradient)
      if (allocated(this%hessian)) deallocate (this%hessian)
      if (allocated(this%dipole)) deallocate (this%dipole)
      call this%reset()
   end subroutine result_destroy

   subroutine result_reset(this)
      class(calculation_result_t), intent(inout) :: this
      this%energy = 0.0_dp
      this%has_energy = .false.
      this%has_gradient = .false.
      this%has_hessian = .false.
      this%has_dipole = .false.
   end subroutine result_reset

end module mqc_result_types
