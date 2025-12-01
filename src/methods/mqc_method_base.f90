module mqc_method_base
   use pic_types, only: dp
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   implicit none
   private

   public :: qc_method_t

   !> Abstract base type for all quantum chemistry methods
   type, abstract :: qc_method_t
   contains
      procedure(calc_energy_interface), deferred :: calc_energy
      procedure(calc_gradient_interface), deferred :: calc_gradient
   end type qc_method_t

   abstract interface
      !> Calculate energy for a fragment
      subroutine calc_energy_interface(this, fragment, result)
         import :: qc_method_t, calculation_result_t, physical_fragment_t
         class(qc_method_t), intent(in) :: this
         type(physical_fragment_t), intent(in) :: fragment
         type(calculation_result_t), intent(out) :: result
      end subroutine calc_energy_interface

      !> Calculate energy and gradient for a fragment
      subroutine calc_gradient_interface(this, fragment, result)
         import :: qc_method_t, calculation_result_t, physical_fragment_t
         class(qc_method_t), intent(in) :: this
         type(physical_fragment_t), intent(in) :: fragment
         type(calculation_result_t), intent(out) :: result
      end subroutine calc_gradient_interface
   end interface

end module mqc_method_base
