!! Extended Tight-Binding (xTB) quantum chemistry method implementation
module mqc_method_xtb
   !! Provides GFN1-xTB and GFN2-xTB methods via the tblite library,
   !! implementing the abstract method interface for energy and gradient calculations.
   use pic_types, only: dp
   use mqc_method_base, only: qc_method_t
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_error, only: ERROR_GENERIC, ERROR_VALIDATION
   use pic_logger, only: logger => global_logger

   ! tblite imports (reuse from mqc_mbe)
   use mctc_env, only: wp, error_type
   use mctc_io, only: structure_type, new
   use pic_timer, only: timer_type
   use tblite_context_type, only: context_type
   use tblite_wavefunction, only: wavefunction_type, new_wavefunction, get_molecular_dipole_moment
   use tblite_container, only: container_type
   use tblite_solvation, only: solvation_type, solvation_input, alpb_input, &
                               new_solvation, new_solvation_cds, new_solvation_shift, &
                               cds_input, shift_input, cpcm_input, cpcm_solvation
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
      ! Solvation settings (leave solvent unallocated for gas phase)
      character(len=:), allocatable :: solvent  !! Solvent name: "water", "ethanol", etc.
      character(len=:), allocatable :: solvation_model  !! "alpb" (default), "gbsa", or "cpcm"
      logical :: use_cds = .true.               !! Include non-polar CDS terms (not for CPCM)
      logical :: use_shift = .true.             !! Include solution state shift (not for CPCM)
      ! CPCM-specific settings
      real(wp) :: dielectric = -1.0_wp          !! Direct dielectric constant (-1 = use solvent lookup)
      integer :: cpcm_nang = 110                !! Number of angular points for CPCM cavity
      real(wp) :: cpcm_rscale = 1.0_wp          !! Radii scaling for CPCM cavity
   contains
      procedure :: calc_energy => xtb_calc_energy      !! Energy-only calculation
      procedure :: calc_gradient => xtb_calc_gradient  !! Energy + gradient calculation
      procedure :: calc_hessian => xtb_calc_hessian  !! Placeholder for Hessian calculation
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
      real(wp) :: dipole_wp(3)

      if (this%verbose) then
         print *, "XTB: Calculating energy using ", this%variant
         print *, "XTB: Fragment has", fragment%n_atoms, "atoms"
         print *, "XTB: nelec =", fragment%nelec
         print *, "XTB: charge =", fragment%charge
         if (allocated(this%solvent)) then
            if (allocated(this%solvation_model)) then
               print *, "XTB: Solvation: ", trim(this%solvation_model), " with solvent = ", this%solvent
            else
               print *, "XTB: Solvation: alpb with solvent = ", this%solvent
            end if
         else
            print *, "XTB: Solvation: none (gas phase)"
         end if
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
         call result%error%set(ERROR_VALIDATION, "Unknown XTB variant: "//this%variant)
         result%has_error = .true.
         return
      end select

      if (allocated(error)) then
         call result%error%set(ERROR_GENERIC, "Failed to create XTB calculator")
         result%has_error = .true.
         return
      end if

      ! Add solvation if configured (either solvent name or direct dielectric)
      if (allocated(this%solvent) .or. this%dielectric > 0.0_wp) then
         if (allocated(this%solvation_model)) then
            call add_solvation_to_calc(calc, mol, this%solvent, this%solvation_model, this%variant, &
                                       this%use_cds, this%use_shift, this%dielectric, &
                                       this%cpcm_nang, this%cpcm_rscale, error)
         else
            call add_solvation_to_calc(calc, mol, this%solvent, "alpb", this%variant, &
                                       this%use_cds, this%use_shift, this%dielectric, &
                                       this%cpcm_nang, this%cpcm_rscale, error)
         end if
         if (allocated(error)) then
            call result%error%set(ERROR_GENERIC, "Failed to add solvation: "//error%message)
            result%has_error = .true.
            return
         end if
      end if

      ! Create wavefunction and run single point calculation
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, this%kt)
      energy = 0.0_wp

      verbosity = merge(1, 0, this%verbose)
      call xtb_singlepoint(ctx, mol, calc, wfn, this%accuracy, energy, verbosity=verbosity)

      ! Compute molecular dipole moment from wavefunction
      dipole_wp(:) = matmul(mol%xyz, wfn%qat(:, 1)) + sum(wfn%dpat(:, :, 1), 2)

      ! Store result (XTB is a semi-empirical method, store as SCF energy)
      result%energy%scf = real(energy, dp)
      result%has_energy = .true.

      ! Store dipole moment
      allocate (result%dipole(3))
      result%dipole = real(dipole_wp, dp)
      result%has_dipole = .true.

      if (this%verbose) then
         print *, "XTB: Energy =", result%energy%total()
         print *, "XTB: Dipole (e*Bohr) =", result%dipole
         print *, "XTB: Dipole magnitude (Debye) =", norm2(result%dipole)*2.541746_dp
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
      real(wp) :: dipole_wp(3)

      if (this%verbose) then
         print *, "XTB: Calculating gradient using ", this%variant
         print *, "XTB: Fragment has", fragment%n_atoms, "atoms"
         print *, "XTB: nelec =", fragment%nelec
         print *, "XTB: charge =", fragment%charge
         if (allocated(this%solvent)) then
            if (allocated(this%solvation_model)) then
               print *, "XTB: Solvation: ", trim(this%solvation_model), " with solvent = ", this%solvent
            else
               print *, "XTB: Solvation: alpb with solvent = ", this%solvent
            end if
         else
            print *, "XTB: Solvation: none (gas phase)"
         end if
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
         call result%error%set(ERROR_VALIDATION, "Unknown XTB variant: "//this%variant)
         result%has_error = .true.
         return
      end select

      if (allocated(error)) then
         call result%error%set(ERROR_GENERIC, "Failed to create XTB calculator")
         result%has_error = .true.
         return
      end if

      ! Add solvation if configured (either solvent name or direct dielectric)
      if (allocated(this%solvent) .or. this%dielectric > 0.0_wp) then
         if (allocated(this%solvation_model)) then
            call add_solvation_to_calc(calc, mol, this%solvent, this%solvation_model, this%variant, &
                                       this%use_cds, this%use_shift, this%dielectric, &
                                       this%cpcm_nang, this%cpcm_rscale, error)
         else
            call add_solvation_to_calc(calc, mol, this%solvent, "alpb", this%variant, &
                                       this%use_cds, this%use_shift, this%dielectric, &
                                       this%cpcm_nang, this%cpcm_rscale, error)
         end if
         if (allocated(error)) then
            call result%error%set(ERROR_GENERIC, "Failed to add solvation: "//error%message)
            result%has_error = .true.
            return
         end if
      end if

      ! Allocate gradient and sigma arrays (initialize to zero)
      allocate (gradient(3, fragment%n_atoms))
      allocate (sigma(3, 3))
      gradient = 0.0_wp
      sigma = 0.0_wp

      ! Create wavefunction and run single point calculation with gradient
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, this%kt, grad=.true.)
      energy = 0.0_wp

      verbosity = merge(1, 0, this%verbose)
      call xtb_singlepoint(ctx, mol, calc, wfn, this%accuracy, energy, &
                           gradient=gradient, sigma=sigma, verbosity=verbosity)

      ! Compute molecular dipole moment from wavefunction
      dipole_wp(:) = matmul(mol%xyz, wfn%qat(:, 1)) + sum(wfn%dpat(:, :, 1), 2)

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

      ! Store dipole moment
      allocate (result%dipole(3))
      result%dipole = real(dipole_wp, dp)
      result%has_dipole = .true.

      if (this%verbose) then
         print *, "XTB: Energy =", result%energy%total()
         print *, "XTB: Gradient norm =", sqrt(sum(result%gradient**2))
         print *, "XTB: Dipole (e*Bohr) =", result%dipole
         print *, "XTB: Dipole magnitude (Debye) =", norm2(result%dipole)*2.541746_dp
         print *, "XTB: Gradient calculation complete"
      end if

      deallocate (num, xyz, gradient, sigma)

   end subroutine xtb_calc_gradient

   subroutine xtb_calc_hessian(this, fragment, result)
      !! Calculate Hessian using finite differences of gradients
      !!
      !! Since tblite does not natively support analytic Hessians, this routine
      !! computes the Hessian numerically via central finite differences:
      !!   H[i,j] = (grad_j(x_i + h) - grad_j(x_i - h)) / (2h)
      !!
      !! This requires 6N gradient calculations (forward and backward for each coordinate)
      use mqc_finite_differences, only: generate_perturbed_geometries, displaced_geometry_t, &
                                        finite_diff_hessian_from_gradients, finite_diff_dipole_derivatives, &
                                        DEFAULT_DISPLACEMENT
      use pic_logger, only: logger => global_logger
      use pic_io, only: to_char

      class(xtb_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      type(displaced_geometry_t), allocatable :: forward_geoms(:), backward_geoms(:)
      real(dp), allocatable :: forward_gradients(:, :, :)  ! (n_displacements, 3, n_atoms)
      real(dp), allocatable :: backward_gradients(:, :, :)  ! (n_displacements, 3, n_atoms)
      real(dp), allocatable :: forward_dipoles(:, :)   ! (n_displacements, 3) for IR intensities
      real(dp), allocatable :: backward_dipoles(:, :)  ! (n_displacements, 3) for IR intensities
      type(calculation_result_t) :: disp_result
      real(dp) :: displacement
      integer :: n_atoms, n_displacements, i
      logical :: compute_dipole_derivs

      n_atoms = fragment%n_atoms
      n_displacements = 3*n_atoms
      displacement = DEFAULT_DISPLACEMENT

      if (this%verbose) then
         call logger%info("XTB: Computing Hessian via finite differences")
         call logger%info("  Method: Central differences of gradients")
         call logger%info("  Atoms: "//to_char(n_atoms))
         call logger%info("  Gradient calculations needed: "//to_char(2*n_displacements))
         call logger%info("  Finite difference step size: "//to_char(displacement)//" Bohr")
      end if

      ! Generate all perturbed geometries
      call generate_perturbed_geometries(fragment, displacement, forward_geoms, backward_geoms)

      ! Allocate storage for gradients at displaced geometries
      allocate (forward_gradients(n_displacements, 3, n_atoms))
      allocate (backward_gradients(n_displacements, 3, n_atoms))

      ! Allocate storage for dipoles at displaced geometries (for IR intensities)
      allocate (forward_dipoles(n_displacements, 3))
      allocate (backward_dipoles(n_displacements, 3))
      forward_dipoles = 0.0_dp
      backward_dipoles = 0.0_dp
      compute_dipole_derivs = .true.  ! Will be set to false if any dipole is missing

      ! Compute gradients at all forward-displaced geometries
      if (this%verbose) then
         call logger%info("  Computing forward-displaced gradients...")
      end if
      do i = 1, n_displacements

         ! Forward
         call this%calc_gradient(forward_geoms(i)%geometry, disp_result)
         if (disp_result%has_error .or. .not. disp_result%has_gradient) then
            call result%error%set(ERROR_GENERIC, "Failed to compute gradient for forward displacement "//to_char(i))
            result%has_error = .true.
            call disp_result%destroy()
            return
         end if
         forward_gradients(i, :, :) = disp_result%gradient
         ! Capture dipole for IR intensity calculation
         if (disp_result%has_dipole) then
            forward_dipoles(i, :) = disp_result%dipole
         else
            compute_dipole_derivs = .false.
         end if
         call disp_result%destroy()

         ! Backward
         call this%calc_gradient(backward_geoms(i)%geometry, disp_result)
         if (disp_result%has_error .or. .not. disp_result%has_gradient) then
            call result%error%set(ERROR_GENERIC, "Failed to compute gradient for backward displacement "//to_char(i))
            result%has_error = .true.
            call disp_result%destroy()
            return
         end if
         backward_gradients(i, :, :) = disp_result%gradient
         ! Capture dipole for IR intensity calculation
         if (disp_result%has_dipole) then
            backward_dipoles(i, :) = disp_result%dipole
         else
            compute_dipole_derivs = .false.
         end if
         call disp_result%destroy()

      end do
      if (this%verbose) then
         call logger%info("  Forward and backward gradient calculations complete ")
      end if
      ! Compute Hessian from finite differences
      if (this%verbose) then
         call logger%info("  Assembling Hessian matrix...")
      end if
      call finite_diff_hessian_from_gradients(fragment, forward_gradients, backward_gradients, &
                                              displacement, result%hessian)

      ! Compute dipole derivatives for IR intensity calculation
      if (compute_dipole_derivs) then
         call finite_diff_dipole_derivatives(n_atoms, forward_dipoles, backward_dipoles, &
                                             displacement, result%dipole_derivatives)
         result%has_dipole_derivatives = .true.
         if (this%verbose) then
            call logger%info("  Dipole derivatives computed for IR intensities")
         end if
      end if

      ! Also compute energy and gradient at reference geometry for completeness
      call this%calc_gradient(fragment, disp_result)
      if (disp_result%has_error) then
         result%error = disp_result%error
         result%has_error = .true.
         call disp_result%destroy()
         return
      end if
      result%energy = disp_result%energy
      result%has_energy = disp_result%has_energy
      if (disp_result%has_gradient) then
         allocate (result%gradient(3, n_atoms))
         result%gradient = disp_result%gradient
         result%has_gradient = .true.
      end if
      call disp_result%destroy()

      result%has_hessian = .true.

      if (this%verbose) then
         call logger%info("  Hessian calculation complete")
      end if

      ! Cleanup
      deallocate (forward_gradients, backward_gradients)
      deallocate (forward_dipoles, backward_dipoles)
      do i = 1, n_displacements
         call forward_geoms(i)%destroy()
         call backward_geoms(i)%destroy()
      end do
      deallocate (forward_geoms, backward_geoms)

   end subroutine xtb_calc_hessian

   subroutine add_solvation_to_calc(calc, mol, solvent, solvation_model, method, use_cds, use_shift, &
                                    dielectric, cpcm_nang, cpcm_rscale, error)
      !! Add implicit solvation model to XTB calculator
      !!
      !! Adds ALPB, GBSA, or CPCM solvation. For ALPB/GBSA, optionally adds CDS and shift corrections.
      !! CPCM does not support CDS or shift corrections.
      type(xtb_calculator), intent(inout) :: calc
      type(structure_type), intent(in) :: mol
      character(len=*), intent(in) :: solvent           !! Solvent name (can be empty if dielectric > 0)
      character(len=*), intent(in) :: solvation_model   !! "alpb", "gbsa", or "cpcm"
      character(len=*), intent(in) :: method            !! "gfn1" or "gfn2"
      logical, intent(in) :: use_cds, use_shift
      real(wp), intent(in) :: dielectric                !! Direct dielectric constant (-1 = use solvent lookup)
      integer, intent(in) :: cpcm_nang                  !! Angular grid points for CPCM
      real(wp), intent(in) :: cpcm_rscale               !! Radii scaling for CPCM
      type(error_type), allocatable, intent(out) :: error

      type(solvation_input) :: solv_input
      class(container_type), allocatable :: cont
      class(solvation_type), allocatable :: solv
      logical :: use_alpb
      real(wp) :: eps

      ! Handle CPCM model separately
      if (trim(solvation_model) == 'cpcm') then
         ! CPCM does not support CDS or shift - silently skip them
         ! (use_cds and use_shift are ignored for CPCM)

         ! Get dielectric constant (direct value or lookup from solvent name)
         if (dielectric > 0.0_wp) then
            eps = dielectric
         else if (len_trim(solvent) > 0) then
            eps = get_solvent_dielectric(solvent)
            if (eps < 0.0_wp) then
               allocate (error)
               error%message = "Unknown solvent for CPCM dielectric lookup: "//trim(solvent)
               return
            end if
         else
            allocate (error)
            error%message = "CPCM requires either solvent name or dielectric constant"
            return
         end if

         ! Create CPCM solvation
         allocate (solv_input%cpcm)
         solv_input%cpcm%dielectric_const = eps
         solv_input%cpcm%nang = cpcm_nang
         solv_input%cpcm%rscale = cpcm_rscale

         call new_solvation(solv, mol, solv_input, error)
         if (allocated(error)) return
         call move_alloc(solv, cont)
         call calc%push_back(cont)

         return
      end if

      ! For ALPB/GBSA, we need a solvent name
      if (len_trim(solvent) == 0) then
         allocate (error)
         error%message = "ALPB/GBSA solvation requires a solvent name"
         return
      end if

      ! Determine if using ALPB or GBSA (GBSA = ALPB with alpb flag false)
      use_alpb = .true.
      if (trim(solvation_model) == 'gbsa') then
         use_alpb = .false.
      end if

      ! 1. Add ALPB/GBSA (polar electrostatic solvation)
      allocate (solv_input%alpb)
      solv_input%alpb%solvent = solvent
      solv_input%alpb%alpb = use_alpb

      call new_solvation(solv, mol, solv_input, error, method)
      if (allocated(error)) return
      call move_alloc(solv, cont)
      call calc%push_back(cont)

      deallocate (solv_input%alpb)

      ! 2. Add CDS (non-polar: cavity, dispersion, surface) if requested
      if (use_cds) then
         allocate (solv_input%cds)
         solv_input%cds%solvent = solvent

         call new_solvation_cds(solv, mol, solv_input, error, method)
         if (allocated(error)) return
         call move_alloc(solv, cont)
         call calc%push_back(cont)

         deallocate (solv_input%cds)
      end if

      ! 3. Add shift (solution state correction) if requested
      if (use_shift) then
         allocate (solv_input%shift)
         solv_input%shift%solvent = solvent

         call new_solvation_shift(solv, solv_input, error, method)
         if (allocated(error)) return
         call move_alloc(solv, cont)
         call calc%push_back(cont)
      end if
   end subroutine add_solvation_to_calc

   pure function get_solvent_dielectric(solvent_name) result(eps)
      !! Get dielectric constant for a named solvent
      !!
      !! Returns the static dielectric constant (relative permittivity) for common solvents.
      !! Returns -1.0 if the solvent is not found.
      character(len=*), intent(in) :: solvent_name
      real(wp) :: eps

      character(len=32) :: name_lower
      integer :: i

      ! Convert to lowercase for case-insensitive matching
      name_lower = solvent_name
      do i = 1, len_trim(name_lower)
         if (name_lower(i:i) >= 'A' .and. name_lower(i:i) <= 'Z') then
            name_lower(i:i) = char(ichar(name_lower(i:i)) + 32)
         end if
      end do

      select case (trim(name_lower))
         ! Water
      case ('water', 'h2o')
         eps = 78.4_wp
         ! Alcohols
      case ('methanol', 'ch3oh')
         eps = 32.7_wp
      case ('ethanol', 'c2h5oh')
         eps = 24.6_wp
      case ('1-propanol', 'propanol')
         eps = 20.1_wp
      case ('2-propanol', 'isopropanol')
         eps = 19.9_wp
      case ('1-butanol', 'butanol')
         eps = 17.5_wp
      case ('2-butanol')
         eps = 15.8_wp
      case ('1-octanol', 'octanol')
         eps = 9.9_wp
         ! Polar aprotic
      case ('acetone')
         eps = 20.7_wp
      case ('acetonitrile', 'ch3cn')
         eps = 37.5_wp
      case ('dmso', 'dimethylsulfoxide')
         eps = 46.7_wp
      case ('dmf', 'dimethylformamide')
         eps = 36.7_wp
      case ('thf', 'tetrahydrofuran')
         eps = 7.6_wp
      case ('formamide')
         eps = 109.5_wp
         ! Aromatics
      case ('benzene')
         eps = 2.3_wp
      case ('toluene')
         eps = 2.4_wp
      case ('pyridine')
         eps = 12.4_wp
      case ('aniline')
         eps = 6.9_wp
      case ('nitrobenzene')
         eps = 34.8_wp
      case ('chlorobenzene')
         eps = 5.6_wp
         ! Halogenated
      case ('chloroform', 'chcl3')
         eps = 4.8_wp
      case ('dichloromethane', 'ch2cl2', 'dcm')
         eps = 8.9_wp
      case ('carbon tetrachloride', 'ccl4')
         eps = 2.2_wp
         ! Ethers
      case ('diethylether', 'ether')
         eps = 4.3_wp
      case ('dioxane')
         eps = 2.2_wp
      case ('furan')
         eps = 2.9_wp
         ! Alkanes
      case ('pentane')
         eps = 1.8_wp
      case ('hexane', 'n-hexane')
         eps = 1.9_wp
      case ('cyclohexane')
         eps = 2.0_wp
      case ('heptane', 'n-heptane')
         eps = 1.9_wp
      case ('octane', 'n-octane')
         eps = 1.9_wp
      case ('decane')
         eps = 2.0_wp
      case ('hexadecane')
         eps = 2.0_wp
         ! Other
      case ('nitromethane')
         eps = 35.9_wp
      case ('cs2', 'carbondisulfide')
         eps = 2.6_wp
      case ('ethyl acetate', 'ethylacetate')
         eps = 6.0_wp
      case ('acetic acid', 'aceticacid')
         eps = 6.2_wp
      case ('formic acid', 'formicacid')
         eps = 51.1_wp
      case ('phenol')
         eps = 9.8_wp
      case ('woctanol')
         eps = 8.1_wp
         ! Infinite dielectric (conductor)
      case ('inf')
         eps = 1.0e10_wp
      case default
         eps = -1.0_wp  ! Unknown solvent
      end select
   end function get_solvent_dielectric

end module mqc_method_xtb
