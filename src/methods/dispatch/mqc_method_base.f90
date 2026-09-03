!! Abstract base module for quantum chemistry method implementations
module mqc_method_base
   !! Common interface every quantum chemistry method implements: energy,
   !! gradient and Hessian for a single fragment.
   use pic_types, only: dp
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   implicit none
   private

   public :: qc_method_t  !! Abstract base type for all QC methods

   type, abstract :: qc_method_t
      !! Abstract base type for all quantum chemistry methods
   contains
      procedure(calc_energy_interface), deferred :: calc_energy    !! Energy calculation interface
      procedure(calc_gradient_interface), deferred :: calc_gradient  !! Gradient calculation interface
      procedure(calc_hessian_interface), deferred :: calc_hessian  !! Hessian calculation interface
   end type qc_method_t

   abstract interface
      subroutine calc_energy_interface(this, fragment, result)
         !! Electronic energy of a molecular fragment
         import :: qc_method_t, calculation_result_t, physical_fragment_t
         implicit none
         class(qc_method_t), intent(in) :: this      !! Method instance
         type(physical_fragment_t), intent(in) :: fragment  !! Molecular fragment
         type(calculation_result_t), intent(out) :: result  !! Calculation results
      end subroutine calc_energy_interface

      subroutine calc_gradient_interface(this, fragment, result)
         !! Electronic energy and nuclear gradient of a molecular fragment
         import :: qc_method_t, calculation_result_t, physical_fragment_t
         implicit none
         class(qc_method_t), intent(in) :: this      !! Method instance
         type(physical_fragment_t), intent(in) :: fragment  !! Molecular fragment
         type(calculation_result_t), intent(out) :: result
      end subroutine calc_gradient_interface

      subroutine calc_hessian_interface(this, fragment, result)
         !! Electronic energy, nuclear gradient and Hessian of a molecular fragment
         import :: qc_method_t, calculation_result_t, physical_fragment_t
         implicit none
         class(qc_method_t), intent(in) :: this      !! Method instance
         type(physical_fragment_t), intent(in) :: fragment  !! Molecular fragment
         type(calculation_result_t), intent(out) :: result
      end subroutine calc_hessian_interface
   end interface

end module mqc_method_base
