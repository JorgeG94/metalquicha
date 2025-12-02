module mqc_result_types
   !! Quantum chemistry calculation result containers
   !!
   !! Defines data structures for storing and managing results from
   !! quantum chemistry calculations including energies, gradients, and properties.
   use pic_types, only: dp
   implicit none
   private

   public :: calculation_result_t  !! Main result container type

   type :: calculation_result_t
      !! Container for quantum chemistry calculation results
      !!
      !! Stores computed quantities from QC calculations with flags
      !! indicating which properties have been computed.
      real(dp) :: energy = 0.0_dp               !! Total electronic energy (Hartree)
      real(dp), allocatable :: gradient(:, :)   !! Energy gradient (3, natoms) (Hartree/Bohr)
      real(dp), allocatable :: hessian(:, :)    !! Energy hessian (future implementation)
      real(dp), allocatable :: dipole(:)        !! Dipole moment vector (3) (Debye)

      ! Computation status flags
      logical :: has_energy = .false.    !! Energy has been computed
      logical :: has_gradient = .false.  !! Gradient has been computed
      logical :: has_hessian = .false.   !! Hessian has been computed
      logical :: has_dipole = .false.    !! Dipole moment has been computed
   contains
      procedure :: destroy => result_destroy  !! Clean up allocated memory
      procedure :: reset => result_reset      !! Reset all values and flags
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
