program test_method_api
   use pic_types, only: dp
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_method_base, only: qc_method_t
   use mqc_method_hf, only: hf_method_t
   use mqc_method_xtb, only: xtb_method_t
   use mqc_result_types, only: calculation_result_t
   implicit none

   type(physical_fragment_t) :: h2_fragment
   type(hf_method_t) :: hf_method
   type(xtb_method_t) :: xtb_method
   type(calculation_result_t) :: result

   ! Create a simple H2 molecule
   call setup_h2_fragment(h2_fragment)

   print *, "============================================"
   print *, "Testing Method API with H2 molecule"
   print *, "============================================"
   print *, ""

   ! Test HF method
   print *, "--------------------------------------------"
   print *, "Testing HF method (dummy implementation)"
   print *, "--------------------------------------------"
   hf_method%options%verbose = .true.

   call hf_method%calc_energy(h2_fragment, result)
   print *, ""
   print *, "Result from HF:"
   print *, "  Energy:", result%energy%total()
   print *, "  Has energy:", result%has_energy
   print *, ""

   call result%destroy()

   ! Test XTB method
   print *, "--------------------------------------------"
   print *, "Testing XTB/GFN2 method (real implementation)"
   print *, "--------------------------------------------"
   xtb_method%variant = "gfn2"
   xtb_method%verbose = .true.

   call xtb_method%calc_energy(h2_fragment, result)
   print *, ""
   print *, "Result from XTB/GFN2:"
   print *, "  Energy:", result%energy%total()
   print *, "  Has energy:", result%has_energy
   print *, ""

   call result%destroy()

   ! Test gradient calculation
   print *, "--------------------------------------------"
   print *, "Testing gradient calculation with XTB/GFN2"
   print *, "--------------------------------------------"
   call xtb_method%calc_gradient(h2_fragment, result)
   print *, ""
   print *, "Result from XTB/GFN2 gradient:"
   print *, "  Energy:", result%energy%total()
   print *, "  Has gradient:", result%has_gradient
   if (result%has_gradient) then
      print *, "  Gradient shape:", shape(result%gradient)
   end if
   print *, ""

   call result%destroy()

   ! Cleanup
   call h2_fragment%destroy()

   print *, "============================================"
   print *, "Method API test completed successfully!"
   print *, "============================================"

contains

   subroutine setup_h2_fragment(frag)
      type(physical_fragment_t), intent(out) :: frag

      ! H2 molecule at ~0.74 Angstrom bond length
      ! Coordinates in Bohr
      frag%n_atoms = 2

      allocate (frag%element_numbers(2))
      allocate (frag%coordinates(3, 2))

      ! Both atoms are hydrogen (Z=1)
      frag%element_numbers = [1, 1]

      ! H-H along z-axis, centered at origin
      ! 0.74 Å = 1.4 Bohr approximately
      frag%coordinates(:, 1) = [0.0_dp, 0.0_dp, -0.7_dp]
      frag%coordinates(:, 2) = [0.0_dp, 0.0_dp, 0.7_dp]

      ! Electronic structure info
      frag%charge = 0
      frag%multiplicity = 1  ! Singlet

      ! Compute number of electrons
      call frag%compute_nelec()

      print *, "H2 fragment setup:"
      print *, "  n_atoms:", frag%n_atoms
      print *, "  charge:", frag%charge
      print *, "  multiplicity:", frag%multiplicity
      print *, "  nelec:", frag%nelec
      print *, ""

   end subroutine setup_h2_fragment

end program test_method_api
