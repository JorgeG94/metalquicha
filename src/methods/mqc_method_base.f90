module mqc_method_base
   !! Abstract base module for quantum chemistry method implementations
   !!
   !! Defines the common interface that all quantum chemistry methods must implement,
   !! providing a unified API for energy and gradient calculations.
   use pic_types, only: dp
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   implicit none
   private

   public :: qc_method_t  !! Abstract base type for all QC methods

   type, abstract :: qc_method_t
      !! Abstract base type for all quantum chemistry methods
      !!
      !! Defines the required interface for energy and gradient calculations
      !! that must be implemented by all concrete method types (XTB, HF, etc.).
   contains
      procedure(calc_energy_interface), deferred :: calc_energy    !! Energy calculation interface
      procedure(calc_gradient_interface), deferred :: calc_gradient !! Gradient calculation interface
   end type qc_method_t

   abstract interface
      subroutine calc_energy_interface(this, fragment, result)
         !! Interface for energy-only calculations
         !!
         !! Computes the electronic energy for a molecular fragment
         !! using the specified quantum chemistry method.
         import :: qc_method_t, calculation_result_t, physical_fragment_t
         class(qc_method_t), intent(in) :: this      !! Method instance
         type(physical_fragment_t), intent(in) :: fragment  !! Molecular fragment
         type(calculation_result_t), intent(out) :: result  !! Calculation results
      end subroutine calc_energy_interface

      subroutine calc_gradient_interface(this, fragment, result)
         !! Interface for energy and gradient calculations
         !!
         !! Computes both electronic energy and nuclear gradients for a
         !! molecular fragment using the specified quantum chemistry method.
         import :: qc_method_t, calculation_result_t, physical_fragment_t
         class(qc_method_t), intent(in) :: this      !! Method instance
         type(physical_fragment_t), intent(in) :: fragment  !! Molecular fragment
         type(calculation_result_t), intent(out) :: result
      end subroutine calc_gradient_interface
   end interface

end module mqc_method_base
