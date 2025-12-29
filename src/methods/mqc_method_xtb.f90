!! Extended Tight-Binding (xTB) quantum chemistry method implementation
module mqc_method_xtb
   !! Provides GFN1-xTB and GFN2-xTB methods via the tblite library,
   !! implementing the abstract method interface for energy and gradient calculations.
   use pic_types, only: dp
   use mqc_method_base, only: qc_method_t
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t

   ! tblite imports (reuse from mqc_mbe)
   use mctc_env, only: wp, error_type
   use mctc_io, only: structure_type, new
   use tblite_context_type, only: context_type
   use tblite_wavefunction, only: wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only: xtb_calculator
   use tblite_xtb_gfn1, only: new_gfn1_calculator
   use tblite_xtb_gfn2, only: new_gfn2_calculator
   use tblite_xtb_singlepoint, only: xtb_singlepoint

   implicit none
   private

   public :: xtb_method_t  !! XTB method implementation type

   type, extends(qc_method_t) :: xtb_method_t
      !! Extended Tight-Binding (xTB) method implementation
      !!
      !! Concrete implementation of the abstract quantum chemistry method
      !! interface for GFN1-xTB and GFN2-xTB calculations via tblite.
      character(len=:), allocatable :: variant  !! XTB variant: "gfn1" or "gfn2"
      logical :: verbose = .false.              !! Print calculation details
      real(wp) :: accuracy = 0.01_wp            !! Numerical accuracy parameter
      real(wp) :: kt = 300.0_wp*3.166808578545117e-06_wp  !! Electronic temperature (300 K)
   contains
      procedure :: calc_energy => xtb_calc_energy      !! Energy-only calculation
      procedure :: calc_gradient => xtb_calc_gradient  !! Energy + gradient calculation
   end type xtb_method_t

contains

   subroutine xtb_calc_energy(this, fragment, result)
      !! Calculate electronic energy using Extended Tight-Binding (xTB) method
      class(xtb_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      ! tblite calculation variables
      type(error_type), allocatable :: error
      type(structure_type) :: mol
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: num(:)
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp) :: energy
      type(context_type) :: ctx
      integer :: verbosity

      if (this%verbose) then
         print *, "XTB: Calculating energy using ", this%variant
         print *, "XTB: Fragment has", fragment%n_atoms, "atoms"
         print *, "XTB: nelec =", fragment%nelec
         print *, "XTB: charge =", fragment%charge
      end if

      ! Convert fragment to tblite format
      allocate (num(fragment%n_atoms))
      allocate (xyz(3, fragment%n_atoms))

      num = fragment%element_numbers
      xyz = fragment%coordinates  ! Already in Bohr

      ! Create molecular structure
      ! charge is real(wp), multiplicity converted to uhf (unpaired electrons)
      call new(mol, num, xyz, charge=real(fragment%charge, wp), &
               uhf=fragment%multiplicity - 1)

      ! Select and create appropriate GFN calculator
      select case (this%variant)
      case ("gfn1")
         call new_gfn1_calculator(calc, mol, error)
      case ("gfn2")
         call new_gfn2_calculator(calc, mol, error)
      case default
         error stop "Unknown XTB variant: "//this%variant
      end select

      if (allocated(error)) then
         error stop "Failed to create XTB calculator"
      end if

      ! Create wavefunction and run single point calculation
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, this%kt)
      energy = 0.0_wp

      verbosity = merge(1, 0, this%verbose)
      call xtb_singlepoint(ctx, mol, calc, wfn, this%accuracy, energy, verbosity=verbosity)

      ! Store result (XTB is a semi-empirical method, store as SCF energy)
      result%energy%scf = real(energy, dp)
      result%has_energy = .true.

      if (this%verbose) then
         print *, "XTB: Energy =", result%energy%total()
      end if

      deallocate (num, xyz)

   end subroutine xtb_calc_energy

   subroutine xtb_calc_gradient(this, fragment, result)
      !! Calculate energy gradient using Extended Tight-Binding (xTB) method
      class(xtb_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      ! tblite calculation variables
      type(error_type), allocatable :: error
      type(structure_type) :: mol
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: num(:)
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp) :: energy
      type(context_type) :: ctx
      integer :: verbosity
      real(wp), allocatable :: gradient(:, :)
      real(wp), allocatable :: sigma(:, :)

      if (this%verbose) then
         print *, "XTB: Calculating gradient using ", this%variant
         print *, "XTB: Fragment has", fragment%n_atoms, "atoms"
         print *, "XTB: nelec =", fragment%nelec
         print *, "XTB: charge =", fragment%charge
      end if

      ! Convert fragment to tblite format
      allocate (num(fragment%n_atoms))
      allocate (xyz(3, fragment%n_atoms))

      num = fragment%element_numbers
      xyz = fragment%coordinates  ! Already in Bohr

      ! Create molecular structure
      call new(mol, num, xyz, charge=real(fragment%charge, wp), &
               uhf=fragment%multiplicity - 1)

      ! Select and create appropriate GFN calculator
      select case (this%variant)
      case ("gfn1")
         call new_gfn1_calculator(calc, mol, error)
      case ("gfn2")
         call new_gfn2_calculator(calc, mol, error)
      case default
         error stop "Unknown XTB variant: "//this%variant
      end select

      if (allocated(error)) then
         error stop "Failed to create XTB calculator"
      end if

      ! Allocate gradient and sigma arrays
      allocate (gradient(3, fragment%n_atoms))
      allocate (sigma(3, 3))

      ! Create wavefunction and run single point calculation with gradient
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, this%kt, grad=.true.)
      energy = 0.0_wp

      verbosity = merge(1, 0, this%verbose)
      call xtb_singlepoint(ctx, mol, calc, wfn, this%accuracy, energy, &
                           gradient=gradient, sigma=sigma, verbosity=verbosity)

      ! Store results (XTB is a semi-empirical method, store as SCF energy)
      result%energy%scf = real(energy, dp)
      result%has_energy = .true.

      ! Store gradient
      allocate (result%gradient(3, fragment%n_atoms))
      result%gradient = real(gradient, dp)
      result%has_gradient = .true.

      ! Store sigma (stress tensor)
      allocate (result%sigma(3, 3))
      result%sigma = real(sigma, dp)
      result%has_sigma = .true.

      if (this%verbose) then
         print *, "XTB: Energy =", result%energy%total()
         print *, "XTB: Gradient norm =", sqrt(sum(result%gradient**2))
         print *, "XTB: Gradient calculation complete"
      end if

      deallocate (num, xyz, gradient, sigma)

   end subroutine xtb_calc_gradient

end module mqc_method_xtb
