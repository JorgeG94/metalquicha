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
   public :: run_uks_scf    !! Drive an unrestricted (open-shell) SCF
   public :: spin_occupations  !! Electron count + multiplicity -> alpha/beta counts

   integer, parameter, public :: SCF_GUESS_CORE = 0  !! F = H
   integer, parameter, public :: SCF_GUESS_GWH = 1   !! Generalized Wolfsberg-Helmholz
   integer, parameter, public :: SCF_GUESS_SAC = 2   !! Superposition of atomic coefficients

   real(dp), parameter :: GWH_K = 1.75_dp
      !! The Wolfsberg-Helmholz constant. 1.75 is the value in universal use.

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

      ! Unrestricted results. `occupied`/`density` above hold the alpha channel
      ! and the TOTAL density respectively when `unrestricted` is set, so the
      ! restricted consumers keep working unchanged.
      logical :: unrestricted = .false.
      real(dp), allocatable :: occupied_beta(:, :)
      real(dp), allocatable :: orbital_energies_beta(:)
      real(dp), allocatable :: density_beta(:, :)   !! Beta density alone
      integer :: n_occupied_beta = 0
      real(dp) :: spin_squared = 0.0_dp             !! <S^2>, for spin contamination
      real(dp) :: dipole(3) = 0.0_dp                !! Electric dipole, a.u.
      logical :: has_dipole = .false.
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

   subroutine build_guess_fock(guess, core_hamiltonian, overlap, guess_fock)
      !! Initial Fock matrix for the first diagonalization
      !!
      !! The core guess (F = H) ignores all electron repulsion, which leaves
      !! the orbital ordering badly wrong -- on a radical that can put the hole
      !! in the wrong orbital and converge tidily onto an excited state.
      !!
      !! GWH scales the off-diagonal elements by the overlap:
      !!
      !!   F_ij = 1/2 K S_ij (H_ii + H_jj),   F_ii = H_ii
      !!
      !! which is a crude stand-in for the screening the core guess omits, and
      !! costs nothing beyond matrices already in hand.
      integer, intent(in) :: guess
      real(dp), intent(in) :: core_hamiltonian(:, :), overlap(:, :)
      real(dp), intent(out) :: guess_fock(:, :)

      integer :: i, j, n

      n = size(core_hamiltonian, 1)
      select case (guess)
      case (SCF_GUESS_GWH)
         do j = 1, n
            do i = 1, n
               if (i == j) then
                  guess_fock(i, j) = core_hamiltonian(i, j)
               else
                  guess_fock(i, j) = 0.5_dp*GWH_K*overlap(i, j) &
                                     *(core_hamiltonian(i, i) + core_hamiltonian(j, j))
               end if
            end do
         end do
      case default
         guess_fock = core_hamiltonian
      end select
   end subroutine build_guess_fock

   subroutine sac_guess_fock(system, core_hamiltonian, guess_alpha, guess_beta, &
                             unrestricted, fock_a, fock_b, error)
      !! First Fock matrix from superposed atomic coefficients
      !!
      !! The guess coefficients carry the summed atomic occupations, which is
      !! generally not the molecule's own electron count -- they are only ever
      !! used to build this matrix, never occupied.
      !!
      !! Restricted case: with C_all the alpha and beta blocks side by side,
      !! D_total = C_all C_all^T, and since K is linear in the density,
      !! K[C_all] = K[D_total]. The restricted Fock wants -1/2 K[D_total],
      !! hence the explicit half. The XC side wants a coefficient matrix whose
      !! implied density is D_total under its own D = 2 C C^T convention, so
      !! C_all is scaled by 1/sqrt(2).
      type(cuest_system_t), intent(inout) :: system
      real(dp), intent(in) :: core_hamiltonian(:, :)
      real(dp), intent(in) :: guess_alpha(:, :), guess_beta(:, :)
      logical, intent(in) :: unrestricted
      real(dp), intent(inout) :: fock_a(:, :), fock_b(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: density_a(:, :), density_b(:, :), density_total(:, :)
      real(dp), allocatable :: coulomb(:, :), exch_a(:, :), exch_b(:, :)
      real(dp), allocatable :: vxc_a(:, :), vxc_b(:, :), combined(:, :)
      real(dp) :: xc_energy
      integer :: n_ao, n_a, n_b

      n_ao = size(core_hamiltonian, 1)
      n_a = size(guess_alpha, 2)
      n_b = size(guess_beta, 2)

      allocate (density_a(n_ao, n_ao), density_b(n_ao, n_ao), density_total(n_ao, n_ao))
      allocate (coulomb(n_ao, n_ao), exch_a(n_ao, n_ao), exch_b(n_ao, n_ao))
      allocate (vxc_a(n_ao, n_ao), vxc_b(n_ao, n_ao))

      call pic_gemm(guess_alpha, guess_alpha, density_a, transb='T', alpha=1.0_dp, beta=0.0_dp)
      call pic_gemm(guess_beta, guess_beta, density_b, transb='T', alpha=1.0_dp, beta=0.0_dp)
      density_total = density_a + density_b

      call system%compute_coulomb(density_total, coulomb, error)
      vxc_a = 0.0_dp
      vxc_b = 0.0_dp
      exch_a = 0.0_dp
      exch_b = 0.0_dp

      if (unrestricted) then
         if (system%needs_exchange) then
            call system%compute_exchange(guess_alpha, exch_a, error, n_occupied=n_a)
            call system%compute_exchange(guess_beta, exch_b, error, n_occupied=n_b)
         end if
         if (system%has_xc) then
            call system%compute_xc_uks(guess_alpha, guess_beta, xc_energy, vxc_a, vxc_b, error)
         end if
         if (error%has_error()) return
         fock_a = core_hamiltonian + coulomb - exch_a + vxc_a
         fock_b = core_hamiltonian + coulomb - exch_b + vxc_b
      else
         allocate (combined(n_ao, n_a + n_b))
         combined(:, 1:n_a) = guess_alpha
         combined(:, n_a + 1:n_a + n_b) = guess_beta
         if (system%needs_exchange) then
            call system%compute_exchange(combined, exch_a, error, n_occupied=n_a + n_b)
         end if
         if (system%has_xc) then
            combined = combined/sqrt(2.0_dp)
            call system%compute_xc(combined, xc_energy, vxc_a, error)
         end if
         deallocate (combined)
         if (error%has_error()) return
         fock_a = core_hamiltonian + coulomb - 0.5_dp*exch_a + vxc_a
         fock_b = fock_a
      end if

      deallocate (density_a, density_b, density_total, coulomb, exch_a, exch_b, vxc_a, vxc_b)
   end subroutine sac_guess_fock

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
                          use_diis, diis_size, verbose, result, error, guess, &
                          guess_alpha, guess_beta)
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
      integer, intent(in), optional :: guess   !! SCF_GUESS_*, default GWH
      real(dp), intent(in), optional :: guess_alpha(:, :), guess_beta(:, :)
         !! Superposed atomic coefficients, when guess is SCF_GUESS_SAC

      integer :: guess_type
      real(dp), allocatable :: fock_scratch(:, :)
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

      guess_type = SCF_GUESS_GWH
      if (present(guess)) guess_type = guess

      n_ao = int(system%n_ao)
      n_occ = n_electrons/2

      if (mod(n_electrons, 2) /= 0) then
         call error%set(ERROR_VALIDATION, &
                        "cuEST RHF: odd electron count -- this is the restricted path; "// &
                        "an open-shell system should have reached run_uks_scf")
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
      if (guess_type == SCF_GUESS_SAC .and. present(guess_alpha) .and. present(guess_beta)) then
         allocate (fock_scratch(n_ao, n_ao))
         call sac_guess_fock(system, core_hamiltonian, guess_alpha, guess_beta, &
                             .false., fock, fock_scratch, error)
         deallocate (fock_scratch)
         if (error%has_error()) return
      else
         call build_guess_fock(guess_type, core_hamiltonian, overlap, fock)
      end if
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
      call system%compute_dipole(density, result%dipole, error)
      result%has_dipole = .not. error%has_error()
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

   pure subroutine spin_occupations(n_electrons, multiplicity, n_alpha, n_beta, ok)
      !! Alpha and beta occupations from the electron count and multiplicity
      !!
      !! n_alpha - n_beta = multiplicity - 1 and n_alpha + n_beta = n_electrons,
      !! so both are fixed. `ok` is false when the two are inconsistent -- an
      !! even electron count with an even multiplicity, for instance.
      integer, intent(in) :: n_electrons, multiplicity
      integer, intent(out) :: n_alpha, n_beta
      logical, intent(out) :: ok

      integer :: unpaired

      unpaired = multiplicity - 1
      ok = (multiplicity >= 1) .and. (unpaired <= n_electrons) .and. &
           (mod(n_electrons - unpaired, 2) == 0)
      if (.not. ok) then
         n_alpha = 0
         n_beta = 0
         return
      end if
      n_beta = (n_electrons - unpaired)/2
      n_alpha = n_beta + unpaired
   end subroutine spin_occupations

   function spin_contamination(occ_alpha, occ_beta, overlap, n_alpha, n_beta) result(s_squared)
      !! <S^2> for a UHF/UKS determinant
      !!
      !!   <S^2> = S_z(S_z+1) + n_beta - sum_ij |<phi_i^a|phi_j^b>|^2
      !!
      !! The exact value for a pure spin state is S_z(S_z+1); the excess is
      !! spin contamination, and reporting it is the cheapest way to notice
      !! that an open-shell answer is not describing the state intended.
      real(dp), intent(in) :: occ_alpha(:, :), occ_beta(:, :), overlap(:, :)
      integer, intent(in) :: n_alpha, n_beta
      real(dp) :: s_squared

      real(dp), allocatable :: scratch(:, :), mo_overlap(:, :)
      real(dp) :: sz

      sz = 0.5_dp*real(n_alpha - n_beta, dp)
      s_squared = sz*(sz + 1.0_dp) + real(n_beta, dp)
      if (n_alpha == 0 .or. n_beta == 0) return

      allocate (scratch(size(overlap, 1), n_beta), mo_overlap(n_alpha, n_beta))
      call pic_gemm(overlap, occ_beta(:, 1:n_beta), scratch)
      call pic_gemm(occ_alpha(:, 1:n_alpha), scratch, mo_overlap, transa='T')
      s_squared = s_squared - sum(mo_overlap**2)
      deallocate (scratch, mo_overlap)
   end function spin_contamination

   subroutine run_uks_scf(system, atomic_numbers, coordinates, n_electrons, multiplicity, &
                          max_iterations, energy_tolerance, density_tolerance, &
                          use_diis, diis_size, verbose, result, error, guess, &
                          guess_alpha, guess_beta)
      !! Drive an unrestricted (open-shell) SCF to convergence
      !!
      !! Spin conventions, matching what cuEST's density-fitted routines expect:
      !!
      !!   D^a = C_a C_a^T,  D^b = C_b C_b^T,  D^t = D^a + D^b   (no factor 2)
      !!   F^a = H + J[D^t] - K[C_a] + Vxc_a
      !!   E   = tr(D^t H) + 1/2 tr(D^t J) - 1/2 (tr(D^a K_a) + tr(D^b K_b)) + Exc
      !!
      !! These reduce to the restricted expressions when D^a = D^b = D/2, which
      !! is the check that the factors are right.
      !!
      !! DIIS uses the two spin commutators stacked into one error vector, so a
      !! single extrapolation drives both channels.
      type(cuest_system_t), intent(inout) :: system
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)
      integer, intent(in) :: n_electrons
      integer, intent(in) :: multiplicity
      integer, intent(in) :: max_iterations
      real(dp), intent(in) :: energy_tolerance, density_tolerance
      logical, intent(in) :: use_diis
      integer, intent(in) :: diis_size
      logical, intent(in) :: verbose
      type(scf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: guess   !! SCF_GUESS_*, default GWH
      real(dp), intent(in), optional :: guess_alpha(:, :), guess_beta(:, :)
         !! Superposed atomic coefficients, when guess is SCF_GUESS_SAC

      integer :: guess_type
      real(dp), allocatable :: overlap(:, :), kinetic(:, :), potential(:, :)
      real(dp), allocatable :: core_hamiltonian(:, :), transform(:, :)
      real(dp), allocatable :: fock_a(:, :), fock_b(:, :)
      real(dp), allocatable :: density_a(:, :), density_b(:, :), density_total(:, :)
      real(dp), allocatable :: coulomb(:, :), exchange_a(:, :), exchange_b(:, :)
      real(dp), allocatable :: vxc_a(:, :), vxc_b(:, :)
      real(dp), allocatable :: orbitals_a(:, :), orbitals_b(:, :)
      real(dp), allocatable :: energies_a(:), energies_b(:)
      real(dp), allocatable :: occ_a(:, :), occ_b(:, :)
      real(dp), allocatable :: err_a(:, :), err_b(:, :), scratch(:, :), ortho_scratch(:, :)
      real(dp), allocatable :: fock_history(:, :, :, :), error_history(:, :, :, :)
      real(dp), allocatable :: b_matrix(:, :), diis_coefficients(:)
      integer :: n_ao, n_mo, n_alpha, n_beta, iteration, n_stored, i, j
      real(dp) :: electronic_energy, previous_energy, energy_change, error_norm, xc_energy
      logical :: diis_ok, occupations_ok

      guess_type = SCF_GUESS_GWH
      if (present(guess)) guess_type = guess

      n_ao = int(system%n_ao)
      call spin_occupations(n_electrons, multiplicity, n_alpha, n_beta, occupations_ok)
      if (.not. occupations_ok) then
         call error%set(ERROR_VALIDATION, "cuEST UKS: electron count and multiplicity are inconsistent")
         return
      end if

      allocate (overlap(n_ao, n_ao), kinetic(n_ao, n_ao), potential(n_ao, n_ao))
      call system%compute_overlap(overlap, error)
      call system%compute_kinetic(kinetic, error)
      call system%compute_potential(potential, error)
      if (error%has_error()) then
         call error%add_context("mqc_cuest_scf:one-electron integrals (UKS)")
         return
      end if

      allocate (core_hamiltonian(n_ao, n_ao))
      core_hamiltonian = kinetic + potential

      call build_orthogonalizer(overlap, transform, n_mo, error)
      if (error%has_error()) return
      if (n_alpha > n_mo) then
         call error%set(ERROR_VALIDATION, &
                        "cuEST UKS: linear dependence left too few orbitals for the alpha electrons")
         return
      end if

      result%nuclear_repulsion = nuclear_repulsion_energy(atomic_numbers, coordinates)

      allocate (fock_a(n_ao, n_ao), fock_b(n_ao, n_ao))
      if (guess_type == SCF_GUESS_SAC .and. present(guess_alpha) .and. present(guess_beta)) then
         call sac_guess_fock(system, core_hamiltonian, guess_alpha, guess_beta, &
                             .true., fock_a, fock_b, error)
         if (error%has_error()) return
      else
         call build_guess_fock(guess_type, core_hamiltonian, overlap, fock_a)
         fock_b = fock_a
      end if
      call diagonalize_fock(fock_a, transform, orbitals_a, energies_a, error)
      call diagonalize_fock(fock_b, transform, orbitals_b, energies_b, error)
      if (error%has_error()) return

      allocate (occ_a(n_ao, max(n_alpha, 1)), occ_b(n_ao, max(n_beta, 1)))
      ! An unoccupied channel is passed on as a zero column, so it must be
      ! zero rather than whatever was on the stack.
      occ_a = 0.0_dp
      occ_b = 0.0_dp
      allocate (density_a(n_ao, n_ao), density_b(n_ao, n_ao), density_total(n_ao, n_ao))
      call set_spin_density(orbitals_a, n_alpha, occ_a, density_a)
      call set_spin_density(orbitals_b, n_beta, occ_b, density_b)
      density_total = density_a + density_b

      allocate (coulomb(n_ao, n_ao), exchange_a(n_ao, n_ao), exchange_b(n_ao, n_ao))
      allocate (vxc_a(n_ao, n_ao), vxc_b(n_ao, n_ao))
      allocate (err_a(n_ao, n_ao), err_b(n_ao, n_ao))
      allocate (scratch(n_ao, n_ao), ortho_scratch(n_ao, n_mo))
      if (use_diis) then
         ! Index 1 of the leading dimension is alpha, index 2 beta: one DIIS
         ! problem spanning both channels.
         allocate (fock_history(n_ao, n_ao, 2, diis_size))
         allocate (error_history(n_mo, n_mo, 2, diis_size))
      end if
      xc_energy = 0.0_dp
      previous_energy = 0.0_dp
      n_stored = 0
      result%converged = .false.

      if (verbose) then
         if (system%has_xc) then
            write (*, '(A)') "  cuEST UKS (density-fitted J/K, grid XC)"
         else
            write (*, '(A)') "  cuEST UHF (density-fitted J/K)"
         end if
         write (*, '(A,I0,A,I0,A,I0,A,I0)') "    n_ao = ", n_ao, "   n_mo = ", n_mo, &
            "   n_alpha = ", n_alpha, "   n_beta = ", n_beta
         write (*, '(A)') "    iter            energy (Ha)          dE        DIIS error"
      end if

      do iteration = 1, max_iterations
         call system%compute_coulomb(density_total, coulomb, error)

         if (system%needs_exchange) then
            call system%compute_exchange(occ_a(:, 1:n_alpha), exchange_a, error, n_occupied=n_alpha)
            if (n_beta > 0) then
               call system%compute_exchange(occ_b(:, 1:n_beta), exchange_b, error, n_occupied=n_beta)
            else
               exchange_b = 0.0_dp
            end if
         else
            exchange_a = 0.0_dp
            exchange_b = 0.0_dp
         end if

         if (system%has_xc) then
            call system%compute_xc_uks(occ_a(:, 1:n_alpha), occ_b(:, 1:max(n_beta, 1)), &
                                       xc_energy, vxc_a, vxc_b, error)
         else
            vxc_a = 0.0_dp
            vxc_b = 0.0_dp
            xc_energy = 0.0_dp
         end if
         if (error%has_error()) then
            call error%add_context("mqc_cuest_scf:Fock build (UKS)")
            return
         end if

         fock_a = core_hamiltonian + coulomb - exchange_a + vxc_a
         fock_b = core_hamiltonian + coulomb - exchange_b + vxc_b

         electronic_energy = sum(density_total*core_hamiltonian) &
                             + 0.5_dp*sum(density_total*coulomb) &
                             - 0.5_dp*(sum(density_a*exchange_a) + sum(density_b*exchange_b)) &
                             + xc_energy
         energy_change = electronic_energy - previous_energy
         previous_energy = electronic_energy

         call commutator_error(fock_a, density_a, overlap, transform, scratch, ortho_scratch, err_a)
         call commutator_error(fock_b, density_b, overlap, transform, scratch, ortho_scratch, err_b)
         error_norm = sqrt(sum(err_a(1:n_mo, 1:n_mo)**2) + sum(err_b(1:n_mo, 1:n_mo)**2))

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

         if (use_diis) then
            if (n_stored < diis_size) then
               n_stored = n_stored + 1
            else
               fock_history(:, :, :, 1:diis_size - 1) = fock_history(:, :, :, 2:diis_size)
               error_history(:, :, :, 1:diis_size - 1) = error_history(:, :, :, 2:diis_size)
            end if
            fock_history(:, :, 1, n_stored) = fock_a
            fock_history(:, :, 2, n_stored) = fock_b
            error_history(:, :, 1, n_stored) = err_a(1:n_mo, 1:n_mo)
            error_history(:, :, 2, n_stored) = err_b(1:n_mo, 1:n_mo)

            if (n_stored > 1) then
               allocate (b_matrix(n_stored + 1, n_stored + 1))
               b_matrix = -1.0_dp
               b_matrix(n_stored + 1, n_stored + 1) = 0.0_dp
               do i = 1, n_stored
                  do j = 1, n_stored
                     b_matrix(i, j) = sum(error_history(:, :, :, i)*error_history(:, :, :, j))
                  end do
               end do
               call solve_diis(b_matrix, diis_coefficients, diis_ok)
               if (diis_ok) then
                  fock_a = 0.0_dp
                  fock_b = 0.0_dp
                  do i = 1, n_stored
                     fock_a = fock_a + diis_coefficients(i)*fock_history(:, :, 1, i)
                     fock_b = fock_b + diis_coefficients(i)*fock_history(:, :, 2, i)
                  end do
               end if
               deallocate (b_matrix, diis_coefficients)
            end if
         end if

         deallocate (orbitals_a, energies_a, orbitals_b, energies_b)
         call diagonalize_fock(fock_a, transform, orbitals_a, energies_a, error)
         call diagonalize_fock(fock_b, transform, orbitals_b, energies_b, error)
         if (error%has_error()) return

         call set_spin_density(orbitals_a, n_alpha, occ_a, density_a)
         call set_spin_density(orbitals_b, n_beta, occ_b, density_b)
         density_total = density_a + density_b
      end do

      if (.not. result%converged) result%iterations = max_iterations

      result%unrestricted = .true.
      call system%compute_dipole(density_total, result%dipole, error)
      result%has_dipole = .not. error%has_error()
      result%electronic_energy = electronic_energy
      result%xc_energy = xc_energy
      result%total_energy = electronic_energy + result%nuclear_repulsion
      result%orbital_energies = energies_a
      result%orbital_energies_beta = energies_b
      result%density = density_total
      result%density_beta = density_b
      result%occupied = occ_a(:, 1:max(n_alpha, 1))
      result%occupied_beta = occ_b(:, 1:max(n_beta, 1))
      result%n_occupied = n_alpha
      result%n_occupied_beta = n_beta
      result%spin_squared = spin_contamination(occ_a, occ_b, overlap, n_alpha, n_beta)

      if (verbose) then
         write (*, '(A,F12.6,A,F12.6,A)') "    <S^2> = ", result%spin_squared, &
            "   (exact ", 0.25_dp*real(n_alpha - n_beta, dp)*(real(n_alpha - n_beta, dp) + 2.0_dp), ")"
      end if

      if (.not. result%converged) then
         call error%set(ERROR_VALIDATION, "cuEST UKS did not converge in the iteration limit")
      end if
   end subroutine run_uks_scf

   subroutine set_spin_density(orbitals, n_occ, occupied, density)
      !! Take the lowest n_occ orbitals of one spin and form D = C C^T
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_occ
      real(dp), intent(inout) :: occupied(:, :)
      real(dp), intent(inout) :: density(:, :)

      density = 0.0_dp
      if (n_occ <= 0) return
      occupied(:, 1:n_occ) = orbitals(:, 1:n_occ)
      ! One electron per spin orbital, so no factor of two here.
      call pic_gemm(occupied(:, 1:n_occ), occupied(:, 1:n_occ), density, &
                    transb='T', alpha=1.0_dp, beta=0.0_dp)
   end subroutine set_spin_density

   subroutine commutator_error(fock, density, overlap, transform, scratch, ortho_scratch, err)
      !! X^T (FDS - SDF) X, the DIIS error for one spin channel
      real(dp), intent(in) :: fock(:, :), density(:, :), overlap(:, :), transform(:, :)
      real(dp), intent(inout) :: scratch(:, :), ortho_scratch(:, :)
      real(dp), intent(inout) :: err(:, :)

      integer :: n_mo

      n_mo = size(transform, 2)
      call pic_gemm(density, overlap, scratch)
      call pic_gemm(fock, scratch, err)
      call pic_gemm(density, fock, scratch)
      call pic_gemm(overlap, scratch, err, alpha=-1.0_dp, beta=1.0_dp)

      call pic_gemm(err, transform, ortho_scratch)
      call pic_gemm(transform, ortho_scratch, err(1:n_mo, 1:n_mo), transa='T')
   end subroutine commutator_error

end module mqc_cuest_scf
