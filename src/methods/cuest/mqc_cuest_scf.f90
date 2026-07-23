!! Restricted Hartree-Fock SCF driven by cuEST integrals
module mqc_cuest_scf
   !! Closed-shell SCF on the host, with every integral supplied by cuEST.
   !!
   !! cuEST provides S, T, V, J and K but no SCF machinery, so the
   !! orthogonalization, diagonalization, DIIS and energy evaluation all live
   !! here and run on the CPU through pic-blas/LAPACK. For fragment-sized
   !! molecules the O(N^3) host work is small next to the integrals.
   !!
   !! Exchange is density-fitted, which is not a choice: cuEST exposes no
   !! conventional four-index ERI path, so an auxiliary (JKFIT) basis is
   !! mandatory and the energies carry the usual RI fitting error.
   !!
   !! Fock matrix convention. cuEST returns
   !!   J[D]_uv = sum_ls D_ls (uv|ls)
   !!   K[C]_uv = sum_ls sum_i C_il C_is (ul|vs)
   !! where K already carries the functional's exact-exchange fraction, since
   !! that is baked into the DF plan rather than applied afterwards. With the
   !! *total* density D = 2 sum_i C_ui C_vi, the closed-shell Fock matrix is
   !!   F = H + J[D] - K[C_occ] + Vxc
   !! and
   !!   E_elec = sum_uv D_uv H_uv + 1/2 sum_uv D_uv J_uv
   !!            - 1/2 sum_uv D_uv K_uv + Exc
   !!
   !! This covers Hartree-Fock and Kohn-Sham with one expression: HF is the
   !! case Vxc = 0, Exc = 0 and exchange fraction 1. Note that Exc enters
   !! directly and is NOT 1/2 sum D Vxc -- writing the energy as
   !! 1/2 sum D (H + F), which is correct for pure HF, would double count the
   !! XC contribution.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cuest_integrals, only: cuest_system_t
   implicit none
   private

   public :: scf_result_t   !! Converged SCF quantities
   public :: run_rhf_scf    !! Drive a closed-shell SCF to convergence

   real(dp), parameter :: LINEAR_DEPENDENCE_TOL = 1.0e-7_dp
      !! Overlap eigenvalues below this are dropped from the orbital space

   type :: scf_result_t
      !! What an SCF run produces
      real(dp) :: electronic_energy = 0.0_dp  !! Electronic energy, Hartree
      real(dp) :: xc_energy = 0.0_dp          !! Exchange-correlation energy, Hartree
      real(dp) :: nuclear_repulsion = 0.0_dp  !! Nuclear repulsion, Hartree
      real(dp) :: total_energy = 0.0_dp       !! Sum of the two
      integer :: iterations = 0               !! Iterations actually taken
      logical :: converged = .false.          !! Whether both criteria were met
      real(dp), allocatable :: orbital_energies(:)  !! Eigenvalues, Hartree
      real(dp), allocatable :: density(:, :)        !! Total density (n_ao, n_ao)
      real(dp), allocatable :: occupied(:, :)       !! Occupied MOs (n_ao, n_occ)
      integer :: n_occupied = 0                     !! Doubly occupied orbital count
   end type scf_result_t

contains

   function nuclear_repulsion_energy(atomic_numbers, coordinates) result(energy)
      !! Point-charge nuclear repulsion, coordinates in Bohr
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)  !! (3, n_atoms)
      real(dp) :: energy

      integer :: iatom, jatom
      real(dp) :: distance

      energy = 0.0_dp
      do iatom = 1, size(atomic_numbers)
         do jatom = iatom + 1, size(atomic_numbers)
            distance = norm2(coordinates(:, iatom) - coordinates(:, jatom))
            energy = energy + real(atomic_numbers(iatom), dp) &
                     *real(atomic_numbers(jatom), dp)/distance
         end do
      end do
   end function nuclear_repulsion_energy

   subroutine build_orthogonalizer(overlap, transform, n_mo, error)
      !! Canonical orthogonalizer X = U s^(-1/2), with near-null modes dropped
      !!
      !! Canonical rather than symmetric orthogonalization so that a basis with
      !! near-linear dependence -- diffuse functions, or two fragments close
      !! together -- loses the offending modes instead of producing garbage.
      real(dp), intent(in) :: overlap(:, :)                  !! S, (n_ao, n_ao)
      real(dp), allocatable, intent(out) :: transform(:, :)  !! X, (n_ao, n_mo)
      integer, intent(out) :: n_mo                           !! Surviving orbitals
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: eigenvectors(:, :), eigenvalues(:)
      integer :: n_ao, i, kept, info

      n_ao = size(overlap, 1)
      allocate (eigenvectors(n_ao, n_ao), eigenvalues(n_ao))
      eigenvectors = overlap

      call pic_syev(eigenvectors, eigenvalues, jobz='V', uplo='U', info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "SCF: overlap matrix diagonalization failed")
         return
      end if

      ! pic_syev returns eigenvalues ascending, so the discarded ones lead.
      n_mo = count(eigenvalues > LINEAR_DEPENDENCE_TOL)
      if (n_mo == 0) then
         call error%set(ERROR_VALIDATION, "SCF: overlap matrix is singular")
         return
      end if

      allocate (transform(n_ao, n_mo))
      kept = 0
      do i = 1, n_ao
         if (eigenvalues(i) <= LINEAR_DEPENDENCE_TOL) cycle
         kept = kept + 1
         transform(:, kept) = eigenvectors(:, i)/sqrt(eigenvalues(i))
      end do
   end subroutine build_orthogonalizer

   subroutine diagonalize_fock(fock, transform, orbitals, energies, error)
      !! Solve FC = SCe in the orthogonal basis and back-transform
      real(dp), intent(in) :: fock(:, :)                    !! F, (n_ao, n_ao)
      real(dp), intent(in) :: transform(:, :)               !! X, (n_ao, n_mo)
      real(dp), allocatable, intent(out) :: orbitals(:, :)  !! C, (n_ao, n_mo)
      real(dp), allocatable, intent(out) :: energies(:)     !! (n_mo)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: fock_ortho(:, :), scratch(:, :)
      integer :: n_ao, n_mo, info

      n_ao = size(transform, 1)
      n_mo = size(transform, 2)
      allocate (fock_ortho(n_mo, n_mo), scratch(n_ao, n_mo), energies(n_mo))

      ! F' = X^T F X
      call pic_gemm(fock, transform, scratch)
      call pic_gemm(transform, scratch, fock_ortho, transa='T')

      call pic_syev(fock_ortho, energies, jobz='V', uplo='U', info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "SCF: Fock matrix diagonalization failed")
         return
      end if

      ! Back-transform: C = X C'
      allocate (orbitals(n_ao, n_mo))
      call pic_gemm(transform, fock_ortho, orbitals)
   end subroutine diagonalize_fock

   subroutine build_density(occupied, density)
      !! Total closed-shell density D = 2 C_occ C_occ^T
      real(dp), intent(in) :: occupied(:, :)     !! C_occ, (n_ao, n_occ)
      real(dp), intent(inout) :: density(:, :)   !! D, (n_ao, n_ao)

      call pic_gemm(occupied, occupied, density, transb='T', alpha=2.0_dp, beta=0.0_dp)
   end subroutine build_density

   subroutine solve_diis(b_matrix, coefficients, ok)
      !! Solve the small DIIS linear system by Gaussian elimination
      !!
      !! The system is at most (diis_size+1) square, so a dedicated LAPACK
      !! solve would cost more in dependencies than it saves in time.
      real(dp), intent(in) :: b_matrix(:, :)
      real(dp), allocatable, intent(out) :: coefficients(:)
      logical, intent(out) :: ok

      real(dp), allocatable :: augmented(:, :)
      integer :: n, i, j, pivot_row
      real(dp) :: pivot, factor

      n = size(b_matrix, 1)
      allocate (augmented(n, n + 1), coefficients(n))
      augmented(:, 1:n) = b_matrix
      augmented(:, n + 1) = 0.0_dp
      augmented(n, n + 1) = -1.0_dp
      ok = .true.

      do i = 1, n
         pivot_row = i
         do j = i + 1, n
            if (abs(augmented(j, i)) > abs(augmented(pivot_row, i))) pivot_row = j
         end do
         if (pivot_row /= i) then
            augmented([i, pivot_row], :) = augmented([pivot_row, i], :)
         end if

         pivot = augmented(i, i)
         if (abs(pivot) < 1.0e-14_dp) then
            ok = .false.
            return
         end if

         do j = i + 1, n
            factor = augmented(j, i)/pivot
            augmented(j, i:n + 1) = augmented(j, i:n + 1) - factor*augmented(i, i:n + 1)
         end do
      end do

      do i = n, 1, -1
         coefficients(i) = (augmented(i, n + 1) - &
                            sum(augmented(i, i + 1:n)*coefficients(i + 1:n)))/augmented(i, i)
      end do
   end subroutine solve_diis

   subroutine run_rhf_scf(system, atomic_numbers, coordinates, n_electrons, &
                          max_iterations, energy_tolerance, density_tolerance, &
                          use_diis, diis_size, verbose, result, error)
      !! Drive a closed-shell RHF calculation to convergence
      type(cuest_system_t), intent(inout) :: system   !! Live cuEST objects for this molecule
      integer, intent(in) :: atomic_numbers(:)        !! Z per atom
      real(dp), intent(in) :: coordinates(:, :)       !! (3, n_atoms), Bohr
      integer, intent(in) :: n_electrons              !! Total electrons, must be even
      integer, intent(in) :: max_iterations
      real(dp), intent(in) :: energy_tolerance        !! Convergence on |dE|
      real(dp), intent(in) :: density_tolerance       !! Convergence on the DIIS error norm
      logical, intent(in) :: use_diis
      integer, intent(in) :: diis_size                !! DIIS subspace size
      logical, intent(in) :: verbose                  !! Print the iteration table
      type(scf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: overlap(:, :), kinetic(:, :), potential(:, :)
      real(dp), allocatable :: core_hamiltonian(:, :), fock(:, :), density(:, :)
      real(dp), allocatable :: coulomb(:, :), exchange(:, :), xc_potential(:, :)
      real(dp), allocatable :: transform(:, :), orbitals(:, :), orbital_energies(:)
      real(dp), allocatable :: occupied(:, :), diis_error(:, :), error_ortho(:, :)
      real(dp), allocatable :: scratch(:, :), ortho_scratch(:, :)
      real(dp), allocatable :: fock_history(:, :, :), error_history(:, :, :)
      real(dp), allocatable :: b_matrix(:, :), diis_coefficients(:)
      integer :: n_ao, n_mo, n_occ, iteration, n_stored, i, j
      real(dp) :: electronic_energy, previous_energy, energy_change, error_norm
      real(dp) :: xc_energy
      logical :: diis_ok

      n_ao = int(system%n_ao)
      n_occ = n_electrons/2

      if (mod(n_electrons, 2) /= 0) then
         call error%set(ERROR_VALIDATION, &
                        "cuEST SCF: open-shell systems are not supported yet (odd electron count)")
         return
      end if
      if (n_occ > n_ao) then
         call error%set(ERROR_VALIDATION, &
                        "cuEST SCF: more occupied orbitals than basis functions")
         return
      end if

      ! ---- one-electron integrals and the core Hamiltonian -----------------
      allocate (overlap(n_ao, n_ao), kinetic(n_ao, n_ao), potential(n_ao, n_ao))
      call system%compute_overlap(overlap, error)
      call system%compute_kinetic(kinetic, error)
      call system%compute_potential(potential, error)
      if (error%has_error()) then
         call error%add_context("mqc_cuest_scf:one-electron integrals")
         return
      end if

      allocate (core_hamiltonian(n_ao, n_ao))
      core_hamiltonian = kinetic + potential

      call build_orthogonalizer(overlap, transform, n_mo, error)
      if (error%has_error()) return
      if (n_occ > n_mo) then
         call error%set(ERROR_VALIDATION, &
                        "cuEST SCF: linear dependence left too few orbitals to hold the electrons")
         return
      end if

      result%nuclear_repulsion = nuclear_repulsion_energy(atomic_numbers, coordinates)

      ! ---- core guess ------------------------------------------------------
      allocate (fock(n_ao, n_ao))
      fock = core_hamiltonian
      call diagonalize_fock(fock, transform, orbitals, orbital_energies, error)
      if (error%has_error()) return

      allocate (occupied(n_ao, max(n_occ, 1)), density(n_ao, n_ao))
      occupied(:, 1:n_occ) = orbitals(:, 1:n_occ)
      call build_density(occupied(:, 1:n_occ), density)

      allocate (coulomb(n_ao, n_ao), exchange(n_ao, n_ao), diis_error(n_ao, n_ao))
      allocate (xc_potential(n_ao, n_ao))
      xc_energy = 0.0_dp
      allocate (error_ortho(n_mo, n_mo), scratch(n_ao, n_ao), ortho_scratch(n_ao, n_mo))
      if (use_diis) then
         ! The error vectors live in the orthogonal basis, the Fock matrices in the AO basis.
         allocate (fock_history(n_ao, n_ao, diis_size), error_history(n_mo, n_mo, diis_size))
      end if

      previous_energy = 0.0_dp
      n_stored = 0
      result%converged = .false.

      if (verbose) then
         if (system%has_xc) then
            write (*, '(A)') "  cuEST RKS (density-fitted J/K, grid XC)"
         else
            write (*, '(A)') "  cuEST RHF (density-fitted J/K)"
         end if
         write (*, '(A,I0,A,I0,A,I0)') "    n_ao = ", n_ao, "   n_mo = ", n_mo, &
            "   n_occ = ", n_occ
         write (*, '(A)') "    iter            energy (Ha)          dE        DIIS error"
      end if

      ! ---- SCF iterations --------------------------------------------------
      do iteration = 1, max_iterations
         call system%compute_coulomb(density, coulomb, error)

         ! A pure functional has no exact exchange; cuEST would compute K and
         ! then scale it by zero, so its own header advises skipping the call.
         if (system%needs_exchange) then
            call system%compute_exchange(occupied(:, 1:n_occ), exchange, error)
         else
            exchange = 0.0_dp
         end if

         call system%compute_xc(occupied(:, 1:n_occ), xc_energy, xc_potential, error)
         if (error%has_error()) then
            call error%add_context("mqc_cuest_scf:Fock build")
            return
         end if

         fock = core_hamiltonian + coulomb - exchange + xc_potential
         electronic_energy = sum(density*core_hamiltonian) &
                             + 0.5_dp*sum(density*coulomb) &
                             - 0.5_dp*sum(density*exchange) &
                             + xc_energy
         energy_change = electronic_energy - previous_energy
         previous_energy = electronic_energy

         ! Commutator error FDS - SDF, zero at convergence, then projected into
         ! the orthogonal basis so its norm is a basis-independent measure.
         call pic_gemm(density, overlap, scratch)
         call pic_gemm(fock, scratch, diis_error)
         call pic_gemm(density, fock, scratch)
         call pic_gemm(overlap, scratch, diis_error, alpha=-1.0_dp, beta=1.0_dp)

         call pic_gemm(diis_error, transform, ortho_scratch)
         call pic_gemm(transform, ortho_scratch, error_ortho, transa='T')
         error_norm = sqrt(sum(error_ortho**2))

         if (verbose) then
            write (*, '(A,I5,F24.12,2ES14.4)') "    ", iteration, &
               electronic_energy + result%nuclear_repulsion, energy_change, error_norm
         end if

         if (iteration > 1 .and. abs(energy_change) < energy_tolerance &
             .and. error_norm < density_tolerance) then
            result%converged = .true.
            result%iterations = iteration
            exit
         end if

         ! ---- DIIS extrapolation -------------------------------------------
         if (use_diis) then
            if (n_stored < diis_size) then
               n_stored = n_stored + 1
            else
               fock_history(:, :, 1:diis_size - 1) = fock_history(:, :, 2:diis_size)
               error_history(:, :, 1:diis_size - 1) = error_history(:, :, 2:diis_size)
            end if
            fock_history(:, :, n_stored) = fock
            error_history(:, :, n_stored) = error_ortho

            if (n_stored > 1) then
               allocate (b_matrix(n_stored + 1, n_stored + 1))
               b_matrix = -1.0_dp
               b_matrix(n_stored + 1, n_stored + 1) = 0.0_dp
               do i = 1, n_stored
                  do j = 1, n_stored
                     b_matrix(i, j) = sum(error_history(:, :, i)*error_history(:, :, j))
                  end do
               end do

               call solve_diis(b_matrix, diis_coefficients, diis_ok)
               if (diis_ok) then
                  fock = 0.0_dp
                  do i = 1, n_stored
                     fock = fock + diis_coefficients(i)*fock_history(:, :, i)
                  end do
               end if
               deallocate (b_matrix, diis_coefficients)
            end if
         end if

         ! ---- new orbitals and density --------------------------------------
         deallocate (orbitals, orbital_energies)
         call diagonalize_fock(fock, transform, orbitals, orbital_energies, error)
         if (error%has_error()) return

         occupied(:, 1:n_occ) = orbitals(:, 1:n_occ)
         call build_density(occupied(:, 1:n_occ), density)
      end do

      if (.not. result%converged) result%iterations = max_iterations

      result%electronic_energy = electronic_energy
      result%xc_energy = xc_energy
      result%total_energy = electronic_energy + result%nuclear_repulsion
      result%orbital_energies = orbital_energies
      result%density = density
      ! The gradient needs the occupied orbitals and their energies to form
      ! the energy-weighted density, so hand them back rather than recomputing.
      result%occupied = occupied(:, 1:n_occ)
      result%n_occupied = n_occ

      if (.not. result%converged) then
         call error%set(ERROR_VALIDATION, "cuEST SCF did not converge in the iteration limit")
      end if
   end subroutine run_rhf_scf

end module mqc_cuest_scf
