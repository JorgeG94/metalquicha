!! Finite difference utilities for numerical derivatives
module mqc_finite_differences
   !! Provides utilities for generating perturbed geometries and computing
   !! numerical derivatives via finite differences (gradients, Hessians, etc.)
   use pic_types, only: dp
   use mqc_physical_fragment, only: physical_fragment_t
   implicit none
   private

   public :: generate_perturbed_geometries  !! Generate forward/backward displacements
   public :: displaced_geometry_t           !! Container for displaced geometry
   public :: finite_diff_hessian_from_gradients  !! Compute Hessian from gradient differences
   public :: finite_diff_dipole_derivatives  !! Compute dipole derivatives from dipole differences
   public :: copy_and_displace_geometry   !! Copy and displace geometry

   ! Default displacement step size (Bohr)
   real(dp), parameter, public :: DEFAULT_DISPLACEMENT = 0.005_dp  !! ~0.05 Angstrom

   type :: displaced_geometry_t
      !! Container for a single displaced geometry
      integer :: atom_index      !! Which atom was displaced (1-based)
      integer :: coordinate      !! Which coordinate was displaced (1=x, 2=y, 3=z)
      integer :: direction       !! +1 for forward, -1 for backward
      real(dp) :: displacement   !! Displacement magnitude in Bohr
      type(physical_fragment_t) :: geometry  !! The displaced geometry
   contains
      procedure :: destroy => displaced_geometry_destroy
   end type displaced_geometry_t

contains

   subroutine generate_perturbed_geometries(reference_geom, displacement, forward_geoms, backward_geoms)
      !! Generate all forward and backward displaced geometries for finite difference calculations
      !!
      !! For a system with N atoms, this generates:
      !!   - 3N forward-displaced geometries (+x, +y, +z for each atom)
      !!   - 3N backward-displaced geometries (-x, -y, -z for each atom)
      !!
      !! These can be used to compute:
      !!   - Gradient: from energies at ±displacement
      !!   - Hessian: from gradients at ±displacement
      !!
      !! Args:
      !!   reference_geom: The reference geometry to perturb
      !!   displacement: Step size in Bohr (typical: 0.001 Bohr)
      !!   forward_geoms: Output array of forward-displaced geometries (size: 3*n_atoms)
      !!   backward_geoms: Output array of backward-displaced geometries (size: 3*n_atoms)
      type(physical_fragment_t), intent(in) :: reference_geom
      real(dp), intent(in) :: displacement
      type(displaced_geometry_t), intent(out), allocatable :: forward_geoms(:)
      type(displaced_geometry_t), intent(out), allocatable :: backward_geoms(:)

      integer :: n_atoms, n_displacements
      integer :: iatom, icoord, idx
      integer :: i

      n_atoms = reference_geom%n_atoms
      n_displacements = 3*n_atoms  ! x, y, z for each atom

      allocate (forward_geoms(n_displacements))
      allocate (backward_geoms(n_displacements))

      ! Generate all displaced geometries
      idx = 0
      do iatom = 1, n_atoms
         do icoord = 1, 3  ! x, y, z
            idx = idx + 1

            ! Forward displacement (+h)
            forward_geoms(idx)%atom_index = iatom
            forward_geoms(idx)%coordinate = icoord
            forward_geoms(idx)%direction = +1
            forward_geoms(idx)%displacement = displacement
            call copy_and_displace_geometry(reference_geom, iatom, icoord, +displacement, &
                                            forward_geoms(idx)%geometry)

            ! Backward displacement (-h)
            backward_geoms(idx)%atom_index = iatom
            backward_geoms(idx)%coordinate = icoord
            backward_geoms(idx)%direction = -1
            backward_geoms(idx)%displacement = displacement
            call copy_and_displace_geometry(reference_geom, iatom, icoord, -displacement, &
                                            backward_geoms(idx)%geometry)
         end do
      end do

   end subroutine generate_perturbed_geometries

   subroutine copy_and_displace_geometry(reference_geom, atom_idx, coord_idx, displacement, displaced_geom)
      !! Create a copy of reference geometry with one coordinate displaced
      !!
      !! Args:
      !!   reference_geom: Original geometry to copy
      !!   atom_idx: Atom to displace (1-based)
      !!   coord_idx: Coordinate to displace (1=x, 2=y, 3=z)
      !!   displacement: Amount to displace in Bohr (positive or negative)
      !!   displaced_geom: Output displaced geometry
      type(physical_fragment_t), intent(in) :: reference_geom
      integer, intent(in) :: atom_idx, coord_idx
      real(dp), intent(in) :: displacement
      type(physical_fragment_t), intent(out) :: displaced_geom

      ! Copy basic properties
      displaced_geom%n_atoms = reference_geom%n_atoms
      displaced_geom%charge = reference_geom%charge
      displaced_geom%multiplicity = reference_geom%multiplicity
      displaced_geom%nelec = reference_geom%nelec
      displaced_geom%n_caps = reference_geom%n_caps

      ! Allocate and copy arrays
      allocate (displaced_geom%element_numbers(displaced_geom%n_atoms))
      allocate (displaced_geom%coordinates(3, displaced_geom%n_atoms))

      displaced_geom%element_numbers = reference_geom%element_numbers
      displaced_geom%coordinates = reference_geom%coordinates

      ! Copy hydrogen cap information if present
      if (reference_geom%n_caps > 0) then
         allocate (displaced_geom%cap_replaces_atom(displaced_geom%n_caps))
         displaced_geom%cap_replaces_atom = reference_geom%cap_replaces_atom
      end if

      ! Copy gradient redistribution mapping if present
      if (allocated(reference_geom%local_to_global)) then
         allocate (displaced_geom%local_to_global(size(reference_geom%local_to_global)))
         displaced_geom%local_to_global = reference_geom%local_to_global
      end if

      ! Apply displacement to specified coordinate
      displaced_geom%coordinates(coord_idx, atom_idx) = &
         displaced_geom%coordinates(coord_idx, atom_idx) + displacement

      ! Copy basis set if present (same basis, just different geometry)
      if (allocated(reference_geom%basis)) then
         ! Note: Basis set will need to be rebuilt with new coordinates
         ! For now, we don't copy it - it should be set up during calculation
      end if

   end subroutine copy_and_displace_geometry

   subroutine finite_diff_hessian_from_gradients(reference_geom, forward_gradients, backward_gradients, &
                                                 displacement, hessian)
      !! Compute Hessian matrix from finite differences of gradients
      !!
      !! Uses central finite differences: H_ij = (grad_i(+h) - grad_i(-h)) / (2h)
      !!
      !! Args:
      !!   reference_geom: Reference geometry (for dimensioning)
      !!   forward_gradients: Gradients at forward-displaced geometries (3*n_atoms, 3, n_atoms)
      !!   backward_gradients: Gradients at backward-displaced geometries (3*n_atoms, 3, n_atoms)
      !!   displacement: Step size used in Bohr
      !!   hessian: Output Hessian matrix (3*n_atoms, 3*n_atoms)
      type(physical_fragment_t), intent(in) :: reference_geom
      real(dp), intent(in) :: forward_gradients(:, :, :)   !! (n_displacements, 3, n_atoms)
      real(dp), intent(in) :: backward_gradients(:, :, :)  !! (n_displacements, 3, n_atoms)
      real(dp), intent(in) :: displacement
      real(dp), intent(out), allocatable :: hessian(:, :)  !! (3*n_atoms, 3*n_atoms)

      integer :: n_atoms, n_coords
      integer :: iatom, jatom, icoord, jcoord
      integer :: i_global, j_global
      integer :: disp_idx

      n_atoms = reference_geom%n_atoms
      n_coords = 3*n_atoms

      allocate (hessian(n_coords, n_coords))
      hessian = 0.0_dp

      ! Build Hessian using central differences
      ! H[i,j] = d²E/(dx_i dx_j) = (dE/dx_j at x_i+h - dE/dx_j at x_i-h) / (2h)

      disp_idx = 0
      do iatom = 1, n_atoms
         do icoord = 1, 3
            disp_idx = disp_idx + 1
            i_global = 3*(iatom - 1) + icoord

            ! For each displacement, compute derivatives of all gradient components
            do jatom = 1, n_atoms
               do jcoord = 1, 3
                  j_global = 3*(jatom - 1) + jcoord

                  ! Central difference: (grad_j(+h) - grad_j(-h)) / (2h)
                  hessian(i_global, j_global) = &
                     (forward_gradients(disp_idx, jcoord, jatom) - &
                      backward_gradients(disp_idx, jcoord, jatom))/(2.0_dp*displacement)
               end do
            end do
         end do
      end do

      ! Symmetrize the Hessian: H = (H + H^T) / 2
      ! This reduces numerical noise from finite differences
      do i_global = 1, n_coords
         do j_global = i_global + 1, n_coords
            hessian(i_global, j_global) = 0.5_dp*(hessian(i_global, j_global) + hessian(j_global, i_global))
            hessian(j_global, i_global) = hessian(i_global, j_global)
         end do
      end do

   end subroutine finite_diff_hessian_from_gradients

   subroutine displaced_geometry_destroy(this)
      !! Clean up memory for displaced geometry
      class(displaced_geometry_t), intent(inout) :: this
      call this%geometry%destroy()
   end subroutine displaced_geometry_destroy

   subroutine finite_diff_dipole_derivatives(n_atoms, forward_dipoles, backward_dipoles, &
                                             displacement, dipole_derivatives)
      !! Compute dipole moment derivatives via central finite differences
      !!
      !! Computes: d_mu_k/d_x_j = (mu_k(+h) - mu_k(-h)) / (2h)
      !!
      !! This is used for IR intensity calculations, where we need the dipole
      !! moment derivatives with respect to Cartesian nuclear coordinates.
      !!
      !! Args:
      !!   n_atoms: Number of atoms
      !!   forward_dipoles: Dipoles at forward-displaced geometries (3*n_atoms, 3)
      !!   backward_dipoles: Dipoles at backward-displaced geometries (3*n_atoms, 3)
      !!   displacement: Step size in Bohr
      !!   dipole_derivatives: Output derivatives (3, 3*n_atoms) in a.u.
      integer, intent(in) :: n_atoms
      real(dp), intent(in) :: forward_dipoles(:, :)   !! (3*n_atoms, 3)
      real(dp), intent(in) :: backward_dipoles(:, :)  !! (3*n_atoms, 3)
      real(dp), intent(in) :: displacement
      real(dp), intent(out), allocatable :: dipole_derivatives(:, :)  !! (3, 3*n_atoms)

      integer :: n_coords, j, k

      n_coords = 3*n_atoms
      allocate (dipole_derivatives(3, n_coords))

      ! Central difference: d_mu_k/d_x_j = (mu_k(x_j+h) - mu_k(x_j-h)) / (2h)
      do j = 1, n_coords       ! Cartesian coordinate index (which coordinate was displaced)
         do k = 1, 3           ! Dipole component (x, y, z)
            dipole_derivatives(k, j) = (forward_dipoles(j, k) - backward_dipoles(j, k))/ &
                                       (2.0_dp*displacement)
         end do
      end do

   end subroutine finite_diff_dipole_derivatives

end module mqc_finite_differences
