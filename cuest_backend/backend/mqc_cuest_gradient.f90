!! Analytic nuclear gradients for the cuEST-backed SCF
module mqc_cuest_gradient
   !! Assembles the closed-shell SCF nuclear gradient from cuEST's derivative
   !! entry points.
   !!
   !! For a converged SCF the gradient is
   !!
   !!   dE/dR = dE_nuc/dR
   !!         + sum_uv D_uv dH_uv/dR        (kinetic + nuclear attraction)
   !!         - sum_uv W_uv dS_uv/dR        (Pulay term)
   !!         + dE_JK/dR                    (density-fitted Coulomb + exchange)
   !!         + dE_xc/dR                    (Kohn-Sham only)
   !!
   !! with D the total density and W = 2 sum_i eps_i C_ui C_vi the
   !! energy-weighted density. There is no response (CPHF) term: at
   !! convergence the orbitals are stationary, which is exactly what the Pulay
   !! term accounts for. That also means the gradient is only as good as the
   !! SCF convergence -- a loosely converged density gives a gradient wrong at
   !! roughly the square root of the energy error.
   !!
   !! One approximation is inherited from the engine: J and K are density
   !! fitted, so the gradient carries the derivative of the fitting error too.
   !!
   !! The XC term is complete. cuestXCDerivativeRKSCompute returns the whole
   !! derivative for a built-in functional, grid-weight (Becke partition)
   !! terms included -- the separate cuestXCGridDerivativeCompute exists for
   !! the user-supplied-functional path, where the caller evaluates the XC
   !! energy density itself and has to hand the per-grid-point contributions
   !! back. Verified numerically: validation/check_gradient agrees with
   !! central finite differences to 3e-8 Ha/Bohr for PBE0, the same noise
   !! floor as Hartree-Fock, which a missing weight term would not do.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t
   use mqc_cuest_integrals, only: cuest_system_t
   use mqc_cuest_scf, only: scf_result_t
   implicit none
   private

   public :: compute_scf_gradient  !! Total nuclear gradient, Hartree/Bohr

contains

   subroutine nuclear_repulsion_gradient(atomic_numbers, coordinates, gradient)
      !! Gradient of the point-charge nuclear repulsion, Hartree/Bohr
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)  !! (3, n_atoms), Bohr
      real(dp), intent(out) :: gradient(:, :)    !! (3, n_atoms)

      integer :: iatom, jatom
      real(dp) ::  distance, factor
      real(dp) :: separation(3)

      gradient = 0.0_dp
      do iatom = 1, size(atomic_numbers)
         do jatom = 1, size(atomic_numbers)
            if (iatom == jatom) cycle
            separation = coordinates(:, iatom) - coordinates(:, jatom)
            distance = norm2(separation)
            factor = -real(atomic_numbers(iatom), dp)*real(atomic_numbers(jatom), dp) &
                     /(distance**3)
            gradient(:, iatom) = gradient(:, iatom) + factor*separation
         end do
      end do
   end subroutine nuclear_repulsion_gradient

   subroutine energy_weighted_density(occupied, orbital_energies, n_occ, weighted, occupancy)
      !! W = f sum_i eps_i C_ui C_vi, the Pulay-term density
      !!
      !! `occupancy` is the electrons per orbital: 2 for a restricted
      !! calculation (the default), 1 for one spin channel of an unrestricted
      !! one, where the caller sums the two channels.
      real(dp), intent(in) :: occupied(:, :)         !! C_occ, (n_ao, n_occ)
      real(dp), intent(in) :: orbital_energies(:)    !! eps, at least n_occ long
      integer, intent(in) :: n_occ
      real(dp), intent(inout) :: weighted(:, :)      !! (n_ao, n_ao)
      real(dp), intent(in), optional :: occupancy

      real(dp), allocatable :: scaled(:, :)
      real(dp) :: factor
      integer :: i

      factor = 2.0_dp
      if (present(occupancy)) factor = occupancy

      weighted = 0.0_dp
      if (n_occ <= 0) return

      allocate (scaled(size(occupied, 1), n_occ))
      do i = 1, n_occ
         scaled(:, i) = occupied(:, i)*orbital_energies(i)
      end do

      call pic_gemm(scaled, occupied(:, 1:n_occ), weighted, transb="T", &
                    alpha=factor, beta=0.0_dp)
      deallocate (scaled)
   end subroutine energy_weighted_density

   subroutine compute_scf_gradient(system, atomic_numbers, coordinates, scf, gradient, error)
      !! Total analytic gradient for a converged closed-shell SCF
      type(cuest_system_t), intent(inout) :: system   !! Live cuEST objects, same geometry
      integer, intent(in) :: atomic_numbers(:)        !! Z per atom
      real(dp), intent(in) :: coordinates(:, :)       !! (3, n_atoms), Bohr
      type(scf_result_t), intent(in) :: scf           !! Converged SCF
      real(dp), allocatable, intent(out) :: gradient(:, :)  !! (3, n_atoms), Hartree/Bohr
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: term(:, :), weighted(:, :), half_density(:, :)
      real(dp), allocatable :: weighted_beta(:, :)
      integer :: n_atoms, n_ao

      n_atoms = size(atomic_numbers)
      n_ao = size(scf%density, 1)

      allocate (gradient(3, n_atoms), term(3, n_atoms))
      allocate (weighted(n_ao, n_ao), half_density(n_ao, n_ao))

      ! ---- nuclear repulsion (host arithmetic) ------------------------------
      call nuclear_repulsion_gradient(atomic_numbers, coordinates, gradient)

      ! ---- one-electron: kinetic and nuclear attraction ---------------------
      call system%gradient_kinetic(scf%density, term, error)
      if (.not. error%has_error()) gradient = gradient + term

      call system%gradient_potential(scf%density, term, error)
      if (.not. error%has_error()) gradient = gradient + term

      ! ---- Pulay term, which enters with a minus sign -----------------------
      !
      ! Unrestricted needs both spin channels, each with one electron per
      ! orbital rather than two.
      if (scf%unrestricted) then
         allocate (weighted_beta(n_ao, n_ao))
         call energy_weighted_density(scf%occupied, scf%orbital_energies, &
                                      scf%n_occupied, weighted, occupancy=1.0_dp)
         call energy_weighted_density(scf%occupied_beta, scf%orbital_energies_beta, &
                                      scf%n_occupied_beta, weighted_beta, occupancy=1.0_dp)
         weighted = weighted + weighted_beta
         deallocate (weighted_beta)
      else
         call energy_weighted_density(scf%occupied, scf%orbital_energies, scf%n_occupied, weighted)
      end if
      call system%gradient_overlap(weighted, term, error)
      if (.not. error%has_error()) gradient = gradient - term

      ! ---- two-electron ------------------------------------------------------
      if (scf%unrestricted) then
         ! The unrestricted call takes the TOTAL density with densityScale 0.5,
         ! where the restricted one takes the half density with scale 2.
         call system%gradient_two_electron_uks(scf%density, scf%occupied, scf%occupied_beta, &
                                               term, error)
      else
         ! cuEST's DF derivative wants the HALF density with densityScale = 2,
         ! not the total density the SCF carries around.
         half_density = 0.5_dp*scf%density
         call system%gradient_two_electron(half_density, scf%occupied, term, error)
      end if
      if (.not. error%has_error()) gradient = gradient + term

      ! ---- exchange-correlation (no-op for Hartree-Fock) --------------------
      if (scf%unrestricted) then
         call system%gradient_xc_uks(scf%occupied, scf%occupied_beta, term, error)
      else
         call system%gradient_xc(scf%occupied, term, error)
      end if
      if (.not. error%has_error()) gradient = gradient + term

      deallocate (term, weighted, half_density)

      if (error%has_error()) then
         call error%add_context("mqc_cuest_gradient:compute_scf_gradient")
         gradient = 0.0_dp
      end if
   end subroutine compute_scf_gradient

end module mqc_cuest_gradient
