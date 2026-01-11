!! Vibrational frequency analysis from Hessian matrix
module mqc_vibrational_analysis
   !! Computes vibrational frequencies from the mass-weighted Hessian matrix.
   !! Uses LAPACK eigenvalue decomposition via pic-blas interfaces.
   use pic_types, only: dp
   use pic_lapack_interfaces, only: pic_syev, pic_gesvd
   use pic_logger, only:
   use mqc_elements, only: element_mass, element_number_to_symbol
   use pic_logger, only: logger => global_logger
   use mqc_physical_constants, only: AU_TO_CM1, AU_TO_MDYNE_ANG, AU_TO_KMMOL, AMU_TO_AU
   use mqc_thermochemistry, only: thermochemistry_result_t, compute_thermochemistry, print_thermochemistry
   implicit none
   private

   public :: compute_vibrational_frequencies
   public :: compute_vibrational_analysis
   public :: mass_weight_hessian
   public :: project_translation_rotation
   public :: compute_reduced_masses
   public :: compute_force_constants
   public :: compute_cartesian_displacements
   public :: compute_ir_intensities
   public :: print_vibrational_analysis

contains

   subroutine compute_vibrational_frequencies(hessian, element_numbers, frequencies, &
                                              eigenvalues_out, eigenvectors, &
                                              coordinates, project_trans_rot, &
                                              projected_hessian_out)
      !! Compute vibrational frequencies from the Hessian matrix.
      !!
      !! Algorithm:
      !! 1. Mass-weight the Hessian: H_mw = M^{-1/2} * H * M^{-1/2}
      !! 2. Optionally project out translation/rotation modes
      !! 3. Diagonalize H_mw to get eigenvalues
      !! 4. Convert eigenvalues to frequencies in cm⁻¹
      !!
      !! Negative eigenvalues produce negative frequencies (imaginary modes,
      !! indicating transition states or saddle points).
      real(dp), intent(in) :: hessian(:, :)
         !! Hessian matrix in Hartree/Bohr² (3*N x 3*N)
      integer, intent(in) :: element_numbers(:)
         !! Atomic numbers for each atom (N atoms)
      real(dp), allocatable, intent(out) :: frequencies(:)
         !! Vibrational frequencies in cm⁻¹ (3*N modes, or 3*N-6 if projected)
      real(dp), allocatable, intent(out), optional :: eigenvalues_out(:)
         !! Raw eigenvalues from diagonalization (Hartree/Bohr²/amu)
      real(dp), allocatable, intent(out), optional :: eigenvectors(:, :)
         !! Normal mode eigenvectors (3*N x 3*N), columns are modes
      real(dp), intent(in), optional :: coordinates(:, :)
         !! Atomic coordinates in Bohr (3, N) - required for projection
      logical, intent(in), optional :: project_trans_rot
         !! If true, project out translation/rotation modes (requires coordinates)
      real(dp), allocatable, intent(out), optional :: projected_hessian_out(:, :)
         !! Mass-weighted Hessian after trans/rot projection (before diagonalization)

      real(dp), allocatable :: mw_hessian(:, :)
      real(dp), allocatable :: eigenvalues(:)
      integer :: n_coords, info, i
      logical :: do_projection

      n_coords = size(hessian, 1)

      ! Check if projection is requested
      do_projection = .false.
      if (present(project_trans_rot)) then
         if (project_trans_rot) then
            if (.not. present(coordinates)) then
               ! Cannot project without coordinates - fall back to no projection
               call logger%warning("Missing coordinates, not projecting out tran/rot motions")
               do_projection = .false.
            else
               do_projection = .true.
            end if
         end if
      end if

      ! Mass-weight the Hessian
      call mass_weight_hessian(hessian, element_numbers, mw_hessian)

      ! Optionally project out translation/rotation modes
      if (do_projection) then
         call project_translation_rotation(mw_hessian, coordinates, element_numbers)
      end if

      ! Return projected Hessian if requested (before diagonalization destroys it)
      if (present(projected_hessian_out)) then
         allocate (projected_hessian_out(n_coords, n_coords))
         projected_hessian_out = mw_hessian
      end if

      ! Allocate eigenvalue storage
      allocate (eigenvalues(n_coords))

      ! Diagonalize the mass-weighted Hessian
      ! pic_syev overwrites mw_hessian with eigenvectors (if jobz='V', default)
      call pic_syev(mw_hessian, eigenvalues, info=info)

      if (info /= 0) then
         ! Eigenvalue decomposition failed
         call logger%error("Eigenvalue decomposition in vibrational frequencies failed")
         allocate (frequencies(n_coords))
         frequencies = 0.0_dp
         return
      end if

      ! Convert eigenvalues to frequencies in cm⁻¹
      allocate (frequencies(n_coords))
      do i = 1, n_coords
         if (eigenvalues(i) >= 0.0_dp) then
            ! Real frequency
            frequencies(i) = sqrt(eigenvalues(i)*AU_TO_CM1)
         else
            ! Imaginary frequency (negative eigenvalue) - report as negative
            frequencies(i) = -sqrt(abs(eigenvalues(i))*AU_TO_CM1)
         end if
      end do

      ! Return eigenvalues if requested
      if (present(eigenvalues_out)) then
         allocate (eigenvalues_out(n_coords))
         eigenvalues_out = eigenvalues
      end if

      ! Return eigenvectors if requested
      if (present(eigenvectors)) then
         allocate (eigenvectors(n_coords, n_coords))
         eigenvectors = mw_hessian
      end if

      deallocate (eigenvalues, mw_hessian)

   end subroutine compute_vibrational_frequencies

   subroutine compute_vibrational_analysis(hessian, element_numbers, frequencies, &
                                           reduced_masses, force_constants, &
                                           cartesian_displacements, &
                                           eigenvalues_out, eigenvectors_out, &
                                           coordinates, project_trans_rot, &
                                           force_constants_mdyne, &
                                           dipole_derivatives, ir_intensities)
      !! Perform complete vibrational analysis from Hessian matrix.
      !!
      !! This is a convenience wrapper that computes:
      !! - Vibrational frequencies in cm⁻¹
      !! - Reduced masses in amu
      !! - Force constants in Hartree/Bohr² (and optionally mdyne/Å)
      !! - Cartesian displacement vectors (normalized)
      !! - IR intensities in km/mol (if dipole_derivatives provided)
      !!
      !! Optionally projects out translation/rotation modes.
      real(dp), intent(in) :: hessian(:, :)
         !! Hessian matrix in Hartree/Bohr² (3*N x 3*N)
      integer, intent(in) :: element_numbers(:)
         !! Atomic numbers for each atom (N atoms)
      real(dp), allocatable, intent(out) :: frequencies(:)
         !! Vibrational frequencies in cm⁻¹
      real(dp), allocatable, intent(out) :: reduced_masses(:)
         !! Reduced masses in amu
      real(dp), allocatable, intent(out) :: force_constants(:)
         !! Force constants in Hartree/Bohr²
      real(dp), allocatable, intent(out) :: cartesian_displacements(:, :)
         !! Cartesian displacement vectors (3*N x 3*N)
      real(dp), allocatable, intent(out), optional :: eigenvalues_out(:)
         !! Raw eigenvalues from diagonalization
      real(dp), allocatable, intent(out), optional :: eigenvectors_out(:, :)
         !! Mass-weighted eigenvectors
      real(dp), intent(in), optional :: coordinates(:, :)
         !! Atomic coordinates in Bohr (3, N) - required for projection
      logical, intent(in), optional :: project_trans_rot
         !! If true, project out translation/rotation modes
      real(dp), allocatable, intent(out), optional :: force_constants_mdyne(:)
         !! Force constants in mdyne/Å
      real(dp), intent(in), optional :: dipole_derivatives(:, :)
         !! Cartesian dipole derivatives (3, 3*N) in a.u. for IR intensities
      real(dp), allocatable, intent(out), optional :: ir_intensities(:)
         !! IR intensities in km/mol

      real(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: eigenvectors(:, :)

      ! First compute frequencies and eigenvectors
      call compute_vibrational_frequencies(hessian, element_numbers, frequencies, &
                                           eigenvalues_out=eigenvalues, &
                                           eigenvectors=eigenvectors, &
                                           coordinates=coordinates, &
                                           project_trans_rot=project_trans_rot)

      ! Compute reduced masses from eigenvectors
      call compute_reduced_masses(eigenvectors, element_numbers, reduced_masses)

      ! Compute force constants from eigenvalues and reduced masses
      call compute_force_constants(eigenvalues, reduced_masses, force_constants, &
                                   force_constants_mdyne)

      ! Compute Cartesian displacements from eigenvectors
      call compute_cartesian_displacements(eigenvectors, element_numbers, &
                                           cartesian_displacements)

      ! Compute IR intensities if dipole derivatives are provided
      if (present(dipole_derivatives) .and. present(ir_intensities)) then
         call compute_ir_intensities(dipole_derivatives, eigenvectors, element_numbers, &
                                     ir_intensities)
      end if

      ! Optionally return eigenvalues and eigenvectors
      if (present(eigenvalues_out)) then
         allocate (eigenvalues_out(size(eigenvalues)))
         eigenvalues_out = eigenvalues
      end if
      if (present(eigenvectors_out)) then
         allocate (eigenvectors_out(size(eigenvectors, 1), size(eigenvectors, 2)))
         eigenvectors_out = eigenvectors
      end if

      deallocate (eigenvalues, eigenvectors)

   end subroutine compute_vibrational_analysis

   subroutine mass_weight_hessian(hessian, element_numbers, mw_hessian)
      !! Apply mass weighting to Hessian matrix.
      !!
      !! H_mw(i,j) = H(i,j) / sqrt(m_i * m_j)
      !!
      !! where m_i is the mass of the atom corresponding to coordinate i.
      !! Each atom contributes 3 coordinates (x, y, z).
      real(dp), intent(in) :: hessian(:, :)
         !! Input Hessian in Hartree/Bohr² (3*N x 3*N)
      integer, intent(in) :: element_numbers(:)
         !! Atomic numbers for each atom (N atoms)
      real(dp), allocatable, intent(out) :: mw_hessian(:, :)
         !! Mass-weighted Hessian (3*N x 3*N)

      real(dp), allocatable :: inv_sqrt_mass(:)
      integer :: n_atoms, n_coords, iatom, icoord, i, j
      real(dp) :: mass

      n_atoms = size(element_numbers)
      n_coords = 3*n_atoms

      ! Build inverse square root mass vector (each mass repeated 3x for x,y,z)
      allocate (inv_sqrt_mass(n_coords))
      do iatom = 1, n_atoms
         mass = element_mass(element_numbers(iatom))
         do icoord = 1, 3
            inv_sqrt_mass(3*(iatom - 1) + icoord) = 1.0_dp/sqrt(mass)
         end do
      end do

      ! Apply mass weighting: H_mw(i,j) = H(i,j) * inv_sqrt_mass(i) * inv_sqrt_mass(j)
      allocate (mw_hessian(n_coords, n_coords))
      do j = 1, n_coords
         do i = 1, n_coords
            mw_hessian(i, j) = hessian(i, j)*inv_sqrt_mass(i)*inv_sqrt_mass(j)
         end do
      end do

      deallocate (inv_sqrt_mass)

   end subroutine mass_weight_hessian

   subroutine project_translation_rotation(mw_hessian, coordinates, element_numbers)
      !! Project out translation and rotation modes from mass-weighted Hessian.
      !!
      !! Builds 6 vectors (3 translation + 3 rotation) in mass-weighted coordinates,
      !! orthonormalizes them using SVD, then projects them out:
      !!   H_proj = (I - D @ D^T) @ H @ (I - D @ D^T)
      !!
      !! This sets the 6 translation/rotation eigenvalues to exactly zero.
      real(dp), intent(inout) :: mw_hessian(:, :)
         !! Mass-weighted Hessian (modified in place)
      real(dp), intent(in) :: coordinates(:, :)
         !! Atomic coordinates in Bohr (3, N)
      integer, intent(in) :: element_numbers(:)
         !! Atomic numbers for each atom (N atoms)

      real(dp), allocatable :: D(:, :)        ! Translation/rotation vectors (3N, 6)
      real(dp), allocatable :: com(:)         ! Center of mass
      real(dp), allocatable :: r(:, :)        ! Coordinates relative to COM
      real(dp), allocatable :: sqrt_mass(:)   ! sqrt(mass) for each atom
      real(dp), allocatable :: S(:)           ! Singular values
      real(dp), allocatable :: U(:, :)        ! Left singular vectors
      real(dp), allocatable :: VT(:, :)       ! Right singular vectors (transposed)
      real(dp), allocatable :: D_orth(:, :)   ! Orthonormalized D vectors
      real(dp), allocatable :: proj(:, :)     ! Projector matrix
      real(dp), allocatable :: temp(:, :)     ! Temporary matrix
      real(dp) :: total_mass, mass, norm
      integer :: n_atoms, n_coords, iatom, i, j, k, n_modes, info
      integer :: idx

      n_atoms = size(element_numbers)
      n_coords = 3*n_atoms

      ! Allocate arrays
      allocate (D(n_coords, 6))
      allocate (com(3))
      allocate (r(3, n_atoms))
      allocate (sqrt_mass(n_atoms))

      ! Compute sqrt(mass) for each atom and total mass
      total_mass = 0.0_dp
      do iatom = 1, n_atoms
         mass = element_mass(element_numbers(iatom))
         sqrt_mass(iatom) = sqrt(mass)
         total_mass = total_mass + mass
      end do

      ! Compute center of mass
      com = 0.0_dp
      do iatom = 1, n_atoms
         mass = element_mass(element_numbers(iatom))
         com(:) = com(:) + mass*coordinates(:, iatom)
      end do
      com = com/total_mass

      ! Compute coordinates relative to center of mass
      do iatom = 1, n_atoms
         r(:, iatom) = coordinates(:, iatom) - com(:)
      end do

      ! Initialize D to zero
      D = 0.0_dp

      ! Build translation vectors (mass-weighted)
      ! D_trans_k: displacement along axis k, weighted by sqrt(mass)
      do iatom = 1, n_atoms
         idx = 3*(iatom - 1)
         ! Translation along x
         D(idx + 1, 1) = sqrt_mass(iatom)
         ! Translation along y
         D(idx + 2, 2) = sqrt_mass(iatom)
         ! Translation along z
         D(idx + 3, 3) = sqrt_mass(iatom)
      end do

      ! Build rotation vectors (mass-weighted)
      ! D_rot_k: rotation around axis k, proportional to r × e_k, weighted by sqrt(mass)
      do iatom = 1, n_atoms
         idx = 3*(iatom - 1)
         ! Rotation around x-axis: r × e_x = (0, r_z, -r_y)
         D(idx + 2, 4) = sqrt_mass(iatom)*r(3, iatom)
         D(idx + 3, 4) = -sqrt_mass(iatom)*r(2, iatom)
         ! Rotation around y-axis: r × e_y = (-r_z, 0, r_x)
         D(idx + 1, 5) = -sqrt_mass(iatom)*r(3, iatom)
         D(idx + 3, 5) = sqrt_mass(iatom)*r(1, iatom)
         ! Rotation around z-axis: r × e_z = (r_y, -r_x, 0)
         D(idx + 1, 6) = sqrt_mass(iatom)*r(2, iatom)
         D(idx + 2, 6) = -sqrt_mass(iatom)*r(1, iatom)
      end do

      ! Normalize each column of D
      do k = 1, 6
         norm = sqrt(sum(D(:, k)**2))
         if (norm > 1.0e-10_dp) then
            D(:, k) = D(:, k)/norm
         end if
      end do

      ! Orthonormalize D using SVD: D = U @ S @ VT
      ! The orthonormal basis is given by the columns of U corresponding to non-zero singular values
      allocate (S(6))
      allocate (U(n_coords, 6))
      allocate (VT(6, 6))

      ! pic_gesvd(A, S, U, VT, info) - A is input, U and VT are separate outputs
      call pic_gesvd(D, S, U, VT, info=info)

      ! Count non-zero singular values (determines number of modes to project)
      n_modes = 0
      do k = 1, 6
         if (S(k) > 1.0e-10_dp) n_modes = n_modes + 1
      end do

      ! Build orthonormalized D matrix from U (columns with non-zero singular values)
      allocate (D_orth(n_coords, n_modes))
      j = 0
      do k = 1, 6
         if (S(k) > 1.0e-10_dp) then
            j = j + 1
            D_orth(:, j) = U(:, k)
         end if
      end do

      ! Build projector: P = I - D_orth @ D_orth^T
      allocate (proj(n_coords, n_coords))
      proj = 0.0_dp
      do i = 1, n_coords
         proj(i, i) = 1.0_dp
      end do

      ! Subtract D_orth @ D_orth^T
      do i = 1, n_coords
         do j = 1, n_coords
            do k = 1, n_modes
               proj(i, j) = proj(i, j) - D_orth(i, k)*D_orth(j, k)
            end do
         end do
      end do

      ! Apply projection: H_proj = P @ H @ P
      allocate (temp(n_coords, n_coords))

      ! temp = H @ P
      do i = 1, n_coords
         do j = 1, n_coords
            temp(i, j) = 0.0_dp
            do k = 1, n_coords
               temp(i, j) = temp(i, j) + mw_hessian(i, k)*proj(k, j)
            end do
         end do
      end do

      ! H_proj = P @ temp
      do i = 1, n_coords
         do j = 1, n_coords
            mw_hessian(i, j) = 0.0_dp
            do k = 1, n_coords
               mw_hessian(i, j) = mw_hessian(i, j) + proj(i, k)*temp(k, j)
            end do
         end do
      end do

      ! Cleanup
      deallocate (D, com, r, sqrt_mass, S, U, VT, D_orth, proj, temp)

   end subroutine project_translation_rotation

   subroutine compute_reduced_masses(eigenvectors, element_numbers, reduced_masses)
      !! Compute reduced masses for each normal mode.
      !!
      !! The reduced mass μ_k for mode k is defined as:
      !!   μ_k = 1 / Σ_i (L_mw_{i,k}² / m_i)
      !!
      !! where L_mw is the mass-weighted eigenvector (normalized to 1).
      !! This formula arises from the relationship Q_k = Σ_i √m_i * x_i * L_mw_{i,k}
      !! and ensures that the harmonic oscillator relation ω² = k/μ holds.
      real(dp), intent(in) :: eigenvectors(:, :)
         !! Mass-weighted eigenvectors from diagonalization (3*N x 3*N)
         !! Columns are normal modes, assumed normalized (Σ_i L²_{i,k} = 1)
      integer, intent(in) :: element_numbers(:)
         !! Atomic numbers for each atom (N atoms)
      real(dp), allocatable, intent(out) :: reduced_masses(:)
         !! Reduced masses in amu (one per mode)

      integer :: n_atoms, n_coords, iatom, icoord, k, idx
      real(dp) :: mass, sum_over_mass

      n_atoms = size(element_numbers)
      n_coords = 3*n_atoms

      allocate (reduced_masses(n_coords))

      ! For each normal mode k
      do k = 1, n_coords
         sum_over_mass = 0.0_dp

         ! Sum over all 3N coordinates: Σ_i (L²_{i,k} / m_i)
         do iatom = 1, n_atoms
            mass = element_mass(element_numbers(iatom))
            do icoord = 1, 3
               idx = 3*(iatom - 1) + icoord
               sum_over_mass = sum_over_mass + eigenvectors(idx, k)**2/mass
            end do
         end do

         ! μ_k = 1 / Σ_i (L²_{i,k} / m_i)
         if (sum_over_mass > 1.0e-14_dp) then
            reduced_masses(k) = 1.0_dp/sum_over_mass
         else
            ! Near-zero contribution (e.g., trans/rot mode) - assign a large mass
            reduced_masses(k) = 1.0e10_dp
         end if
      end do

   end subroutine compute_reduced_masses

   subroutine compute_force_constants(eigenvalues, reduced_masses, force_constants, &
                                      force_constants_mdyne)
      !! Compute force constants for each normal mode.
      !!
      !! From the harmonic oscillator relation:
      !!   ω² = k/μ  →  k = ω² × μ = eigenvalue × μ
      !!
      !! Returns force constants in both atomic units (Hartree/Bohr²) and mdyne/Å.
      real(dp), intent(in) :: eigenvalues(:)
         !! Eigenvalues from mass-weighted Hessian diagonalization (1/amu)
      real(dp), intent(in) :: reduced_masses(:)
         !! Reduced masses in amu
      real(dp), allocatable, intent(out) :: force_constants(:)
         !! Force constants in atomic units (Hartree/Bohr²)
      real(dp), allocatable, intent(out), optional :: force_constants_mdyne(:)
         !! Force constants in mdyne/Å (common experimental unit)

      integer :: n_modes, k

      n_modes = size(eigenvalues)
      allocate (force_constants(n_modes))

      ! k = eigenvalue × μ (eigenvalue has units Hartree/(Bohr²·amu), μ in amu)
      ! So force_constant has units Hartree/Bohr²
      do k = 1, n_modes
         if (eigenvalues(k) >= 0.0_dp) then
            force_constants(k) = eigenvalues(k)*reduced_masses(k)
         else
            ! Imaginary frequency mode - report absolute value
            force_constants(k) = -abs(eigenvalues(k))*reduced_masses(k)
         end if
      end do

      ! Optionally convert to mdyne/Å
      if (present(force_constants_mdyne)) then
         allocate (force_constants_mdyne(n_modes))
         force_constants_mdyne = force_constants*AU_TO_MDYNE_ANG
      end if

   end subroutine compute_force_constants

   subroutine compute_cartesian_displacements(eigenvectors, element_numbers, &
                                              cartesian_displacements, normalize_max)
      !! Convert mass-weighted eigenvectors to Cartesian displacements.
      !!
      !! The Cartesian displacement for coordinate i in mode k is:
      !!   x_{i,k} = L_mw_{i,k} / √(m_i)
      !!
      !! The output can be normalized in two ways:
      !! - normalize_max=.true. (default): normalize so max|x| = 1 for each mode (Gaussian convention)
      !! - normalize_max=.false.: normalize so Σ_i x²_{i,k} = 1
      real(dp), intent(in) :: eigenvectors(:, :)
         !! Mass-weighted eigenvectors (3*N x 3*N)
      integer, intent(in) :: element_numbers(:)
         !! Atomic numbers for each atom (N atoms)
      real(dp), allocatable, intent(out) :: cartesian_displacements(:, :)
         !! Cartesian displacement vectors (3*N x 3*N), columns are modes
      logical, intent(in), optional :: normalize_max
         !! If true, normalize so max displacement = 1 (default: true)

      integer :: n_atoms, n_coords, iatom, icoord, k, idx
      real(dp) :: mass, inv_sqrt_mass, norm, max_disp
      logical :: use_max_norm

      n_atoms = size(element_numbers)
      n_coords = 3*n_atoms

      use_max_norm = .true.
      if (present(normalize_max)) use_max_norm = normalize_max

      allocate (cartesian_displacements(n_coords, n_coords))

      ! Convert from mass-weighted to Cartesian: x = L_mw / √m
      do k = 1, n_coords
         do iatom = 1, n_atoms
            mass = element_mass(element_numbers(iatom))
            inv_sqrt_mass = 1.0_dp/sqrt(mass)
            do icoord = 1, 3
               idx = 3*(iatom - 1) + icoord
               cartesian_displacements(idx, k) = eigenvectors(idx, k)*inv_sqrt_mass
            end do
         end do
      end do

      ! Normalize each mode
      do k = 1, n_coords
         if (use_max_norm) then
            ! Gaussian convention: normalize so max |displacement| = 1
            max_disp = maxval(abs(cartesian_displacements(:, k)))
            if (max_disp > 1.0e-14_dp) then
               cartesian_displacements(:, k) = cartesian_displacements(:, k)/max_disp
            end if
         else
            ! Standard normalization: Σ_i x²_{i,k} = 1
            norm = sqrt(sum(cartesian_displacements(:, k)**2))
            if (norm > 1.0e-14_dp) then
               cartesian_displacements(:, k) = cartesian_displacements(:, k)/norm
            end if
         end if
      end do

   end subroutine compute_cartesian_displacements

   subroutine compute_ir_intensities(dipole_derivatives, eigenvectors, element_numbers, ir_intensities)
      !! Compute IR intensities from dipole derivatives and normal modes.
      !!
      !! IR intensities are computed by transforming Cartesian dipole derivatives
      !! to normal mode coordinates and computing the squared magnitude.
      !!
      !! For each normal mode i:
      !!   trdip(k) = Σ_j dipd(k,j) * L(j,i) * 1/√m_j
      !!   IR(i) = AU_TO_KMMOL * (trdip(1)² + trdip(2)² + trdip(3)²)
      !!
      !! where:
      !!   dipd(k,j) = ∂μ_k/∂x_j (Cartesian dipole derivative)
      !!   L(j,i) = mass-weighted eigenvector component
      !!   m_j = atomic mass for coordinate j
      !!
      real(dp), intent(in) :: dipole_derivatives(:, :)
         !! Cartesian dipole derivatives (3, 3*N) in atomic units
      real(dp), intent(in) :: eigenvectors(:, :)
         !! Mass-weighted eigenvectors from Hessian diagonalization (3*N x 3*N)
      integer, intent(in) :: element_numbers(:)
         !! Atomic numbers for each atom (N atoms)
      real(dp), allocatable, intent(out) :: ir_intensities(:)
         !! IR intensities in km/mol (one per mode)

      integer :: n_atoms, n_coords, iatom, i, j, k
      real(dp) :: mass, inv_sqrt_mass, trdip(3)

      n_atoms = size(element_numbers)
      n_coords = 3*n_atoms

      allocate (ir_intensities(n_coords))

      ! For each normal mode i
      do i = 1, n_coords
         trdip = 0.0_dp

         ! Transform dipole derivative from Cartesian to normal mode coordinates
         ! trdip(k) = Σ_j dipd(k,j) * L(j,i) * amass_au(j)
         ! where amass_au(j) = 1/√(m_j in atomic units) = 1/√(m_amu * AMU_TO_AU)
         ! This matches xtb's formula in hessian.F90 lines 526-535
         do j = 1, n_coords
            iatom = (j - 1)/3 + 1
            mass = element_mass(element_numbers(iatom))
            ! Convert mass to atomic units (electron masses) before taking sqrt
            inv_sqrt_mass = 1.0_dp/sqrt(mass*AMU_TO_AU)

            do k = 1, 3  ! x, y, z components of dipole
               trdip(k) = trdip(k) + dipole_derivatives(k, j)*eigenvectors(j, i)*inv_sqrt_mass
            end do
         end do

         ! IR intensity = |dμ/dQ|² * conversion factor
         ir_intensities(i) = AU_TO_KMMOL*(trdip(1)**2 + trdip(2)**2 + trdip(3)**2)
      end do

   end subroutine compute_ir_intensities

   subroutine print_vibrational_analysis(frequencies, reduced_masses, force_constants, &
                                         cartesian_displacements, element_numbers, &
                                         force_constants_mdyne, print_displacements, &
                                         n_atoms, ir_intensities, &
                                         coordinates, electronic_energy)
      !! Print vibrational analysis results in a formatted table.
      !!
      !! Output format is similar to Gaussian, with frequencies grouped in columns.
      !! Optionally prints Cartesian displacement vectors for each mode.
      !! If coordinates and electronic_energy are provided, also computes and prints
      !! thermochemistry using the RRHO approximation.
      real(dp), intent(in) :: frequencies(:)
         !! Vibrational frequencies in cm⁻¹
      real(dp), intent(in) :: reduced_masses(:)
         !! Reduced masses in amu
      real(dp), intent(in) :: force_constants(:)
         !! Force constants in Hartree/Bohr² (or mdyne/Å if force_constants_mdyne provided)
      real(dp), intent(in) :: cartesian_displacements(:, :)
         !! Cartesian displacement vectors (3*N x 3*N)
      integer, intent(in) :: element_numbers(:)
         !! Atomic numbers for each atom
      real(dp), intent(in), optional :: force_constants_mdyne(:)
         !! Force constants in mdyne/Å (if provided, these are printed instead)
      logical, intent(in), optional :: print_displacements
         !! If true, print Cartesian displacement vectors (default: true)
      integer, intent(in), optional :: n_atoms
         !! Number of atoms (if not provided, derived from size of element_numbers)
      real(dp), intent(in), optional :: ir_intensities(:)
         !! IR intensities in km/mol
      real(dp), intent(in), optional :: coordinates(:, :)
         !! Atomic coordinates (3, n_atoms) in Bohr - needed for thermochemistry
      real(dp), intent(in), optional :: electronic_energy
         !! Electronic energy in Hartree - needed for thermochemistry

      integer :: n_modes, n_at, n_groups, igroup, imode, iatom, icoord, k
      integer :: mode_start, mode_end, modes_in_group
      logical :: do_print_disp
      character(len=512) :: line
      character(len=16) :: freq_str, mass_str, fc_str, ir_str
      character(len=2) :: elem_sym
      character(len=3) :: coord_label
      real(dp) :: fc_value

      n_modes = size(frequencies)
      if (present(n_atoms)) then
         n_at = n_atoms
      else
         n_at = size(element_numbers)
      end if

      do_print_disp = .true.
      if (present(print_displacements)) do_print_disp = print_displacements

      call logger%info(" ")
      call logger%info("============================================================")
      call logger%info("                  VIBRATIONAL ANALYSIS")
      call logger%info("============================================================")
      call logger%info(" ")

      ! Print in groups of 3 modes (like Gaussian)
      n_groups = (n_modes + 2)/3

      do igroup = 1, n_groups
         mode_start = (igroup - 1)*3 + 1
         mode_end = min(igroup*3, n_modes)
         modes_in_group = mode_end - mode_start + 1

         ! Mode numbers header
         line = "                    "
         do k = mode_start, mode_end
            write (freq_str, '(i12)') k
            line = trim(line)//freq_str
         end do
         call logger%info(trim(line))

         ! Frequencies (show "i" only for significant imaginary frequencies)
         line = " Frequencies --  "
         do k = mode_start, mode_end
            if (frequencies(k) < 0.0_dp .and. abs(frequencies(k)) > 10.0_dp) then
               ! Significant imaginary frequency - show with "i"
               write (freq_str, '(f12.4,a)') abs(frequencies(k)), "i"
            else
               ! Real or near-zero frequency
               write (freq_str, '(f12.4)') abs(frequencies(k))
            end if
            line = trim(line)//freq_str
         end do
         call logger%info(trim(line))

         ! Reduced masses
         line = " Red. masses --  "
         do k = mode_start, mode_end
            write (mass_str, '(f12.4)') reduced_masses(k)
            line = trim(line)//mass_str
         end do
         call logger%info(trim(line))

         ! Force constants
         if (present(force_constants_mdyne)) then
            line = " Frc consts  --  "
            do k = mode_start, mode_end
               write (fc_str, '(f12.4)') force_constants_mdyne(k)
               line = trim(line)//fc_str
            end do
         else
            line = " Frc consts  --  "
            do k = mode_start, mode_end
               write (fc_str, '(f12.6)') force_constants(k)
               line = trim(line)//fc_str
            end do
         end if
         call logger%info(trim(line))

         ! IR intensities (if provided)
         if (present(ir_intensities)) then
            line = " IR Intens  --  "
            do k = mode_start, mode_end
               write (ir_str, '(f12.4)') ir_intensities(k)
               line = trim(line)//ir_str
            end do
            call logger%info(trim(line))
         end if

         ! Cartesian displacements
         if (do_print_disp) then
          call logger%info(" Atom          X         Y         Z       X         Y         Z       X         Y         Z")

            do iatom = 1, n_at
               elem_sym = element_number_to_symbol(element_numbers(iatom))

               ! Build line with atom info and displacements for each mode
               write (line, '(i4,1x,a2)') iatom, elem_sym

               do k = mode_start, mode_end
                  do icoord = 1, 3
                     write (freq_str, '(f10.5)') cartesian_displacements(3*(iatom - 1) + icoord, k)
                     line = trim(line)//freq_str
                  end do
               end do
               call logger%info(trim(line))
            end do
         end if

         call logger%info(" ")
      end do

      ! Summary statistics
      call logger%info("------------------------------------------------------------")
      call logger%info(" Summary:")

      ! Count real vs imaginary frequencies
      block
         integer :: n_real, n_imag, n_zero
         real(dp) :: zero_thresh
         zero_thresh = 10.0_dp  ! frequencies below 10 cm⁻¹ considered "zero"

         n_real = 0
         n_imag = 0
         n_zero = 0
         do k = 1, n_modes
            if (abs(frequencies(k)) < zero_thresh) then
               n_zero = n_zero + 1
            else if (frequencies(k) < 0.0_dp) then
               n_imag = n_imag + 1
            else
               n_real = n_real + 1
            end if
         end do

         write (line, '(a,i5)') "   Total modes:              ", n_modes
         call logger%info(trim(line))
         write (line, '(a,i5)') "   Real frequencies:         ", n_real
         call logger%info(trim(line))
         write (line, '(a,i5)') "   Imaginary frequencies:    ", n_imag
         call logger%info(trim(line))
         write (line, '(a,i5)') "   Near-zero (trans/rot):    ", n_zero
         call logger%info(trim(line))

         if (n_imag > 0) then
            call logger%warning("System has imaginary frequencies - may be a transition state")
         end if
      end block

      call logger%info("============================================================")
      call logger%info(" ")

      ! Compute and print thermochemistry if coordinates and energy are provided
      if (present(coordinates) .and. present(electronic_energy)) then
         block
            type(thermochemistry_result_t) :: thermo_result
            call compute_thermochemistry(coordinates, element_numbers, frequencies, &
                                         n_at, n_modes, thermo_result)
            call print_thermochemistry(thermo_result, electronic_energy)
         end block
      end if

   end subroutine print_vibrational_analysis

end module mqc_vibrational_analysis
