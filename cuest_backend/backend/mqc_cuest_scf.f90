!! Restricted Hartree-Fock SCF driven by cuEST integrals
module mqc_cuest_scf
   !! Closed-shell SCF, with every integral supplied by cuEST.
   !!
   !! cuEST provides S, T, V, J and K but no SCF machinery, so the
   !! orthogonalization, diagonalization, DIIS and energy evaluation all live
   !! here.
   !!
   !! **The restricted SCF loop performs no n_ao^2 transfers.** J, K and Vxc are
   !! built into buffers of their own and assembled into a Fock matrix on the
   !! card; the commutator FDS - SDF and its projection into the orthogonal
   !! basis are cuBLAS; the DIIS history is device-side, so the extrapolation
   !! is a device-to-device combination; the Fock is diagonalized by cuSOLVER
   !! and the density rebuilt from the resulting coefficients without either
   !! ever crossing the bus. The loop closes on itself: what cuEST reads at the
   !! top of an iteration is what cuBLAS and cuSOLVER wrote at the bottom of
   !! the previous one.
   !!
   !! What does cross, per iteration, is scalars: three energy contractions,
   !! the DIIS error norm, one row of the DIIS overlap matrix, and cuSOLVER's
   !! convergence flag. Every one of them is a value the host needs in order to
   !! decide something -- to print, to converge, or to solve the small DIIS
   !! system -- rather than data in transit.
   !!
   !! The guess is uploaded once before the loop; the density, occupied
   !! orbitals and orbital energies are fetched once after it, for the dipole,
   !! the gradient and the printed frontier orbitals. The overlap
   !! diagonalization that builds the orthogonalizer stays host LAPACK: it runs
   !! once, outside the loop, and moving it would buy nothing per iteration.
   !!
   !! The unrestricted path below still round-trips every term. It is the same
   !! transformation twice over with a second spin channel, and has not been
   !! done yet; `commutator_error` and `mqc_diis` are its host reference, and
   !! remain the reference the device versions here were diffed against.
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
   use, intrinsic :: iso_c_binding, only: c_int, c_int64_t
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_diis, only: diis_state_t
   use mqc_diis_device, only: diis_device_t
   use mqc_cuest_integrals, only: cuest_system_t
   use mqc_cuest_context, only: cuest_context_t
   use mqc_cuest_runtime, only: device_sync, cublas_status_check, &
                                cusolver_status_check, copy_int_to_host, copy_to_host
   use cublas, only: cublasDgemm, cublasDnrm2, cublasDcopy, CUBLAS_OP_N, CUBLAS_OP_T, &
                     CUBLAS_FILL_MODE_UPPER
   use cusolver, only: cusolverDnDsyevd, cusolverDnDsyevd_bufferSize, &
                       CUSOLVER_EIG_MODE_VECTOR
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

      call pic_gemm(guess_alpha, guess_alpha, density_a, transb="T", alpha=1.0_dp, beta=0.0_dp)
      call pic_gemm(guess_beta, guess_beta, density_b, transb="T", alpha=1.0_dp, beta=0.0_dp)
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

      call pic_syev(eigenvectors, eigenvalues, jobz="V", uplo="U", info=info)
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
      call pic_gemm(transform, scratch, fock_ortho, transa="T")

      call pic_syev(fock_ortho, energies, jobz="V", uplo="U", info=info)
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

      call pic_gemm(occupied, occupied, density, transb="T", alpha=2.0_dp, beta=0.0_dp)
   end subroutine build_density

   subroutine run_rhf_scf(system, context, atomic_numbers, coordinates, n_electrons, &
                          max_iterations, energy_tolerance, density_tolerance, &
                          use_diis, diis_size, verbose, result, error, guess, &
                          guess_alpha, guess_beta)
      !! Drive a closed-shell RHF calculation to convergence
      type(cuest_system_t), intent(inout) :: system   !! Live cuEST objects for this molecule
      type(cuest_context_t), intent(inout) :: context
         !! Per-rank handle and pools. Needed because the DIIS history is
         !! device-resident and is drawn from the same pools as everything
         !! else, rather than being allocated afresh for each fragment.
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
      real(dp), allocatable :: transform(:, :), orbitals(:, :), orbital_energies(:)
      real(dp), allocatable :: occupied(:, :)
      type(diis_device_t) :: diis
      integer :: n_ao, n_mo, n_occ, iteration
      real(dp) :: electronic_energy, previous_energy, energy_change, error_norm
      real(dp) :: xc_energy, trace_h, trace_j, trace_k
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

      xc_energy = 0.0_dp
      if (use_diis) then
         ! The error vectors live in the orthogonal basis, the Fock matrices in the AO basis.
         call diis%init(context, diis_size, n_ao*n_ao, n_mo*n_mo, error)
         if (error%has_error()) return
      end if

      previous_energy = 0.0_dp
      result%converged = .false.

      if (verbose) then
         if (system%has_xc) then
            write (*, "(A)") "  cuEST RKS (density-fitted J/K, grid XC)"
         else
            write (*, "(A)") "  cuEST RHF (density-fitted J/K)"
         end if
         write (*, "(A,I0,A,I0,A,I0)") "    n_ao = ", n_ao, "   n_mo = ", n_mo, &
            "   n_occ = ", n_occ
         write (*, "(A)") "    iter            energy (Ha)          dE        DIIS error"
      end if

      ! Everything the loop reads goes up here, and nothing of size n_ao^2 goes
      ! either way again until it is over.
      !
      ! H, S and X are iteration-invariant. The density and the occupied
      ! coefficients are not, but they are the *guess*: from the first
      ! diagonalization onwards the device produces its own, so this is the
      ! only time they are uploaded rather than computed in place.
      call system%stage_to(system%d_core, core_hamiltonian, "core Hamiltonian", error)
      call system%stage_to(system%d_ovlp, overlap, "overlap matrix", error)
      call system%stage_to(system%d_transform, transform, "orthogonalizer", error)
      call system%stage_density(density, error)
      call system%stage_occupied(occupied(:, 1:n_occ), error)
      if (error%has_error()) then
         call error%add_context("mqc_cuest_scf:initial upload")
         return
      end if

      call solver_workspace(system, context, n_mo, error)
      if (error%has_error()) return

      ! The guess has served its purpose. Neither the host Fock nor the host
      ! orbitals are read again -- the loop below produces and consumes both on
      ! the device -- so drop them rather than carrying two dead n_ao^2 arrays
      ! for the rest of the run.
      deallocate (fock, orbitals)

      ! ---- SCF iterations --------------------------------------------------
      do iteration = 1, max_iterations
         call system%coulomb_device(system%d_matrix, system%d_j, error)

         ! A pure functional has no exact exchange; cuEST would compute K and
         ! then scale it by zero, so its own header advises skipping the call.
         ! Skipped here means d_k is never written, which is why the assembly
         ! below is told whether to add it rather than adding a zeroed buffer.
         if (system%needs_exchange) then
            call system%exchange_device(system%d_c_occ, system%d_k, error)
         end if

         call system%xc_device(system%d_c_occ, system%d_xc, xc_energy, error)

         ! The one synchronise in the iteration, standing between the cuEST
         ! integrals and the cuBLAS that consumes them: cuEST is free to queue
         ! its work on a stream of its own, and a cuBLAS call reading J before
         ! the Coulomb kernel has finished writing it would produce plausible
         ! numbers rather than an error.
         call device_sync("Fock terms", error)
         if (error%has_error()) then
            call error%add_context("mqc_cuest_scf:Fock build")
            return
         end if

         ! F = H + J - K + Vxc, assembled where the terms already are.
         call system%assemble_fock(system%needs_exchange, system%has_xc, error)

         ! The two energy traces come back as scalars instead of dragging J and
         ! K across with them.
         call system%matrix_dot(system%d_matrix, system%d_j, trace_j, "D.J", error)
         trace_k = 0.0_dp
         if (system%needs_exchange) then
            call system%matrix_dot(system%d_matrix, system%d_k, trace_k, "D.K", error)
         end if

         ! The energy is three device contractions and a scalar cuEST already
         ! handed back. No matrix is needed on the host to compute it.
         call system%matrix_dot(system%d_matrix, system%d_core, trace_h, "D.H", error)
         if (error%has_error()) then
            call error%add_context("mqc_cuest_scf:energy traces")
            return
         end if
         electronic_energy = trace_h &
                             + 0.5_dp*trace_j &
                             - 0.5_dp*trace_k &
                             + xc_energy
         energy_change = electronic_energy - previous_energy
         previous_energy = electronic_energy

         ! Commutator error FDS - SDF, zero at convergence, projected into the
         ! orthogonal basis so its norm is a basis-independent measure. All of
         ! it on the device, so the DIIS error vector is already where the
         ! history wants it.
         call commutator_error_device(system, n_ao, n_mo, error_norm, error)
         if (error%has_error()) then
            call error%add_context("mqc_cuest_scf:commutator")
            return
         end if

         if (verbose) then
            write (*, "(A,I5,F24.12,2ES14.4)") "    ", iteration, &
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
            ! Fock, error vector and history are all on the device: the push is
            ! device-to-device and the extrapolation writes back into d_fock in
            ! place. Nothing crosses but the overlap row, n_stored doubles.
            call diis%push(system%d_fock, system%d_error, error)
            call diis%extrapolate(system%d_fock, diis_ok, error)
            if (error%has_error()) then
               call error%add_context("mqc_cuest_scf:DIIS")
               return
            end if
         end if

         ! ---- new orbitals and density --------------------------------------
         ! Diagonalize, slice off the occupied block and rebuild the density,
         ! all on the card. This closes the loop: the density and coefficients
         ! the next iteration's J, K and Vxc read are already where cuEST wants
         ! them, so nothing of size n_ao^2 crosses the bus at all.
         call diagonalize_fock_device(system, n_ao, n_mo, n_occ, error)
         if (error%has_error()) then
            call error%add_context("mqc_cuest_scf:diagonalization")
            return
         end if
      end do

      if (.not. result%converged) result%iterations = max_iterations

      ! Bring the converged state back, once. Everything below this line is a
      ! consumer on the host -- the dipole, the gradient, the printed orbital
      ! energies -- so this is where the loop's residency ends.
      !
      ! On the converged iteration the loop exits before diagonalizing, so
      ! these hold the orbitals and density that produced the Fock just
      ! accepted, which is the same pairing the host path produced.
      call fetch_scf_state(system, n_ao, n_mo, n_occ, density, occupied, &
                           orbital_energies, error)
      if (error%has_error()) then
         call error%add_context("mqc_cuest_scf:converged state fetch")
         return
      end if

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

      if (verbose) then
         call logger%info("  "//frontier_orbital_text("HOMO", "LUMO", orbital_energies, n_occ))
      end if

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
      call pic_gemm(occ_alpha(:, 1:n_alpha), scratch, mo_overlap, transa="T")
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
      type(diis_state_t) :: diis
      real(dp), allocatable :: diis_fock(:), diis_err(:)
         !! Both spin channels packed end to end. fock_a and fock_b are separate
         !! allocations, so unlike the restricted case they cannot reach the flat
         !! history by sequence association.
      integer :: n_ao, n_mo, n_alpha, n_beta, iteration
      integer :: n_fock_spin, n_err_spin
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
         n_fock_spin = n_ao*n_ao
         n_err_spin = n_mo*n_mo
         call diis%init(diis_size, 2*n_fock_spin, 2*n_err_spin)
         allocate (diis_fock(2*n_fock_spin), diis_err(2*n_err_spin))
      end if
      xc_energy = 0.0_dp
      previous_energy = 0.0_dp
      result%converged = .false.

      if (verbose) then
         if (system%has_xc) then
            write (*, "(A)") "  cuEST UKS (density-fitted J/K, grid XC)"
         else
            write (*, "(A)") "  cuEST UHF (density-fitted J/K)"
         end if
         write (*, "(A,I0,A,I0,A,I0,A,I0)") "    n_ao = ", n_ao, "   n_mo = ", n_mo, &
            "   n_alpha = ", n_alpha, "   n_beta = ", n_beta
         write (*, "(A)") "    iter            energy (Ha)          dE        DIIS error"
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
            write (*, "(A,I5,F24.12,2ES14.4)") "    ", iteration, &
               electronic_energy + result%nuclear_repulsion, energy_change, error_norm
         end if

         if (iteration > 1 .and. abs(energy_change) < energy_tolerance &
             .and. error_norm < density_tolerance) then
            result%converged = .true.
            result%iterations = iteration
            exit
         end if

         if (use_diis) then
            ! Both spins stacked into one error vector, so a single
            ! extrapolation drives both channels.
            diis_fock(1:n_fock_spin) = reshape(fock_a, [n_fock_spin])
            diis_fock(n_fock_spin + 1:) = reshape(fock_b, [n_fock_spin])
            diis_err(1:n_err_spin) = reshape(err_a(1:n_mo, 1:n_mo), [n_err_spin])
            diis_err(n_err_spin + 1:) = reshape(err_b(1:n_mo, 1:n_mo), [n_err_spin])

            call diis%push(diis_fock, diis_err)
            call diis%extrapolate(diis_fock, diis_ok)
            if (diis_ok) then
               fock_a = reshape(diis_fock(1:n_fock_spin), [n_ao, n_ao])
               fock_b = reshape(diis_fock(n_fock_spin + 1:), [n_ao, n_ao])
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
         write (*, "(A,F12.6,A,F12.6,A)") "    <S^2> = ", result%spin_squared, &
            "   (exact ", 0.25_dp*real(n_alpha - n_beta, dp)*(real(n_alpha - n_beta, dp) + 2.0_dp), ")"
         call logger%info("  "//frontier_orbital_text("HOMO alpha", "LUMO alpha", &
                                                      result%orbital_energies, n_alpha))
         call logger%info("  "//frontier_orbital_text("HOMO beta", "LUMO beta", &
                                                      result%orbital_energies_beta, n_beta))
      end if

      if (.not. result%converged) then
         call error%set(ERROR_VALIDATION, "cuEST UKS did not converge in the iteration limit")
      end if
   end subroutine run_uks_scf

   function frontier_orbital_text(homo_label, lumo_label, orbital_energies, n_occ) result(text)
      !! "HOMO: ..., LUMO: ..." for one spin channel, omitting what does not exist
      !!
      !! Neither edge is guaranteed to be there. A spin channel with no
      !! electrons (the beta channel of a one-electron fragment) has no HOMO,
      !! and a saturated orbital space -- helium in a minimal basis, say -- has
      !! no virtual to call a LUMO. Indexing blindly reads past the array.
      character(len=*), intent(in) :: homo_label, lumo_label
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      character(len=:), allocatable :: text

      character(len=64) :: buffer

      text = ""
      if (n_occ >= 1 .and. n_occ <= size(orbital_energies)) then
         write (buffer, '(A,": ",F15.6)') homo_label, orbital_energies(n_occ)
         text = trim(buffer)
      end if

      if (n_occ + 1 >= 1 .and. n_occ + 1 <= size(orbital_energies)) then
         write (buffer, '(A,": ",F15.6)') lumo_label, orbital_energies(n_occ + 1)
         if (len(text) > 0) then
            text = text//", "//trim(buffer)
         else
            text = trim(buffer)
         end if
      end if

      if (len(text) == 0) text = "no frontier orbitals to report"
   end function frontier_orbital_text

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
                    transb="T", alpha=1.0_dp, beta=0.0_dp)
   end subroutine set_spin_density

   subroutine fetch_scf_state(system, n_ao, n_mo, n_occ, density, occupied, &
                              orbital_energies, error)
      !! Bring the density, occupied orbitals and orbital energies back once
      !!
      !! The SCF loop leaves all three on the device. Their host consumers --
      !! the dipole, the gradient's energy-weighted density, the printed
      !! frontier orbitals -- run after convergence, so one transfer each at
      !! the end serves them all.
      type(cuest_system_t), intent(inout) :: system
      integer, intent(in) :: n_ao, n_mo, n_occ
      real(dp), intent(inout) :: density(:, :)   !! (n_ao, n_ao)
      real(dp), intent(inout) :: occupied(:, :)  !! (n_ao, >= n_occ)
      real(dp), allocatable, intent(inout) :: orbital_energies(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: flat(:)

      if (error%has_error()) return

      call system%fetch(system%d_matrix, density, "converged density", error)
      if (error%has_error()) return

      if (allocated(orbital_energies)) deallocate (orbital_energies)
      allocate (orbital_energies(n_mo))
      call copy_to_host(orbital_energies, system%d_eigenvalues, "orbital energies", error)
      if (error%has_error()) return

      if (n_occ <= 0) return

      ! Column-major on both sides, so the occupied block comes back as one
      ! contiguous n_ao*n_occ run with no repacking.
      allocate (flat(n_ao*n_occ))
      call copy_to_host(flat, system%d_c_occ, "occupied orbitals", error)
      if (.not. error%has_error()) occupied(:, 1:n_occ) = reshape(flat, [n_ao, n_occ])
      deallocate (flat)
   end subroutine fetch_scf_state

   subroutine solver_workspace(system, context, n_mo, error)
      !! Size cuSOLVER's workspace for an n_mo eigenproblem and pool it
      !!
      !! The requirement depends on n_mo alone, not on the matrix, so this runs
      !! once before the SCF loop rather than once per iteration.
      type(cuest_system_t), intent(inout) :: system
      type(cuest_context_t), intent(inout) :: context
      integer, intent(in) :: n_mo
      type(error_t), intent(inout) :: error

      integer(c_int) :: lwork

      if (error%has_error()) return

      lwork = 0
      call cusolver_status_check(cusolverDnDsyevd_bufferSize(system%cusolver, &
                                                             CUSOLVER_EIG_MODE_VECTOR, &
                                                             CUBLAS_FILL_MODE_UPPER, &
                                                             int(n_mo, c_int), &
                                                             system%d_fock_ortho, &
                                                             int(n_mo, c_int), &
                                                             system%d_eigenvalues, lwork), &
                                 "cusolverDnDsyevd_bufferSize", error)
      if (error%has_error()) return

      call context%scratch_solver%ensure(int(max(lwork, 1), c_int64_t), &
                                         "eigensolver workspace", error)
      if (error%has_error()) return
      system%d_solver = context%scratch_solver%ptr
      system%solver_lwork = lwork
   end subroutine solver_workspace

   subroutine diagonalize_fock_device(system, n_ao, n_mo, n_occ, error)
      !! Solve FC = SCe on the device and rebuild the density from the result
      !!
      !! The device twin of `diagonalize_fock` plus `build_density`, which the
      !! unrestricted path still uses and which are the reference this was
      !! written against. Doing both here is what keeps the density on the
      !! card: C is produced on the device, the occupied block is sliced off it
      !! on the device, and D = 2 C_occ C_occ^T never touches the host.
      !!
      !! Slicing is free because both are column-major: the first n_occ columns
      !! of C are its first n_ao*n_occ elements, contiguous.
      type(cuest_system_t), intent(inout) :: system
      integer, intent(in) :: n_ao, n_mo, n_occ
      type(error_t), intent(inout) :: error

      integer(c_int) :: m, k, devinfo

      if (error%has_error()) return

      m = int(n_ao, c_int)
      k = int(n_mo, c_int)

      ! F' = X^T F X, through the same scratch the commutator used -- which is
      ! dead by now, the DIIS having already consumed its output.
      call cublas_status_check(cublasDgemm(system%cublas, CUBLAS_OP_N, CUBLAS_OP_N, &
                                           m, k, m, 1.0_dp, system%d_fock, m, &
                                           system%d_transform, m, 0.0_dp, system%d_work, m), &
                               "cublasDgemm(F X)", error)
      call cublas_status_check(cublasDgemm(system%cublas, CUBLAS_OP_T, CUBLAS_OP_N, &
                                           k, k, m, 1.0_dp, system%d_transform, m, &
                                           system%d_work, m, 0.0_dp, system%d_fock_ortho, k), &
                               "cublasDgemm(X^T F X)", error)
      if (error%has_error()) return

      ! Overwrites d_fock_ortho with the eigenvectors C', which is why F' gets
      ! a buffer of its own rather than sharing one with something still live.
      call cusolver_status_check(cusolverDnDsyevd(system%cusolver, &
                                                  CUSOLVER_EIG_MODE_VECTOR, &
                                                  CUBLAS_FILL_MODE_UPPER, k, &
                                                  system%d_fock_ortho, k, &
                                                  system%d_eigenvalues, system%d_solver, &
                                                  system%solver_lwork, system%d_devinfo), &
                                 "cusolverDnDsyevd", error)
      if (error%has_error()) return

      ! cuSOLVER reports non-convergence through this device int, not through
      ! the status above. Skipping the check would let a failed eigensolve pass
      ! for a successful one with meaningless orbitals.
      call copy_int_to_host(devinfo, system%d_devinfo, "eigensolver info", error)
      if (error%has_error()) return
      if (devinfo /= 0) then
         call error%set(ERROR_VALIDATION, &
                        "cuEST SCF: device Fock diagonalization failed (cusolverDnDsyevd info /= 0)")
         return
      end if

      ! Back-transform C = X C'
      call cublas_status_check(cublasDgemm(system%cublas, CUBLAS_OP_N, CUBLAS_OP_N, &
                                           m, k, k, 1.0_dp, system%d_transform, m, &
                                           system%d_fock_ortho, k, 0.0_dp, &
                                           system%d_orbitals, m), &
                               "cublasDgemm(X C')", error)

      if (n_occ <= 0) return

      ! The occupied block, and the density it implies. Two electrons per
      ! spatial orbital, hence the 2.
      call cublas_status_check(cublasDcopy(system%cublas, int(n_ao*n_occ, c_int), &
                                           system%d_orbitals, 1, system%d_c_occ, 1), &
                               "cublasDcopy(occupied orbitals)", error)
      call cublas_status_check(cublasDgemm(system%cublas, CUBLAS_OP_N, CUBLAS_OP_T, &
                                           m, m, int(n_occ, c_int), 2.0_dp, &
                                           system%d_c_occ, m, system%d_c_occ, m, &
                                           0.0_dp, system%d_matrix, m), &
                               "cublasDgemm(2 C_occ C_occ^T)", error)
   end subroutine diagonalize_fock_device

   subroutine commutator_error_device(system, n_ao, n_mo, error_norm, error)
      !! X^T (FDS - SDF) X and its norm, entirely on the device
      !!
      !! The device twin of `commutator_error` below, which the unrestricted
      !! path still uses and which is the reference this was written against.
      !!
      !! Layout. cuBLAS is column-major, which is Fortran's own convention, so
      !! the leading dimensions are the Fortran ones and nothing is transposed
      !! to compensate -- unlike cuEST's row-major matrices (see the header of
      !! mqc_cuest_integrals.f90). That difference is invisible for the
      !! symmetric matrices here (F, D, S) and does not arise for X, which this
      !! code uploaded itself from a Fortran (n_ao, n_mo) array.
      !!
      !! FDS - SDF is antisymmetric rather than symmetric, but it is formed and
      !! consumed on the device without ever crossing to the host, so the
      !! layout question never comes up for it.
      type(cuest_system_t), intent(inout) :: system
      integer, intent(in) :: n_ao, n_mo
      real(dp), intent(out) :: error_norm  !! ||X^T (FDS - SDF) X||_F
      type(error_t), intent(inout) :: error

      integer(c_int) :: m, k

      error_norm = 0.0_dp
      if (error%has_error()) return

      m = int(n_ao, c_int)
      k = int(n_mo, c_int)

      ! work = D S
      call cublas_status_check(cublasDgemm(system%cublas, CUBLAS_OP_N, CUBLAS_OP_N, &
                                           m, m, m, 1.0_dp, system%d_matrix, m, &
                                           system%d_ovlp, m, 0.0_dp, system%d_work, m), &
                               "cublasDgemm(D S)", error)
      ! commutator = F (D S)
      call cublas_status_check(cublasDgemm(system%cublas, CUBLAS_OP_N, CUBLAS_OP_N, &
                                           m, m, m, 1.0_dp, system%d_fock, m, &
                                           system%d_work, m, 0.0_dp, system%d_commutator, m), &
                               "cublasDgemm(F D S)", error)
      ! work = D F
      call cublas_status_check(cublasDgemm(system%cublas, CUBLAS_OP_N, CUBLAS_OP_N, &
                                           m, m, m, 1.0_dp, system%d_matrix, m, &
                                           system%d_fock, m, 0.0_dp, system%d_work, m), &
                               "cublasDgemm(D F)", error)
      ! commutator := -S (D F) + commutator
      call cublas_status_check(cublasDgemm(system%cublas, CUBLAS_OP_N, CUBLAS_OP_N, &
                                           m, m, m, -1.0_dp, system%d_ovlp, m, &
                                           system%d_work, m, 1.0_dp, system%d_commutator, m), &
                               "cublasDgemm(-S D F)", error)

      ! Project into the orthogonal basis, where the norm is a basis-independent
      ! measure. `work` is dead by now -- it last held D F, which the gemm above
      ! has already consumed -- so it carries the (n_ao, n_mo) intermediate.
      call cublas_status_check(cublasDgemm(system%cublas, CUBLAS_OP_N, CUBLAS_OP_N, &
                                           m, k, m, 1.0_dp, system%d_commutator, m, &
                                           system%d_transform, m, 0.0_dp, system%d_work, m), &
                               "cublasDgemm(E X)", error)
      call cublas_status_check(cublasDgemm(system%cublas, CUBLAS_OP_T, CUBLAS_OP_N, &
                                           k, k, m, 1.0_dp, system%d_transform, m, &
                                           system%d_work, m, 0.0_dp, system%d_error, k), &
                               "cublasDgemm(X^T E X)", error)
      if (error%has_error()) return

      ! Blocks, like every cuBLAS scalar result -- but the convergence test
      ! needs this value on the host anyway, so nothing is lost.
      call cublas_status_check(cublasDnrm2(system%cublas, int(n_mo*n_mo, c_int), &
                                           system%d_error, 1, error_norm), &
                               "cublasDnrm2(DIIS error)", error)
      if (error%has_error()) error_norm = 0.0_dp
   end subroutine commutator_error_device

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
      call pic_gemm(transform, ortho_scratch, err(1:n_mo, 1:n_mo), transa="T")
   end subroutine commutator_error

end module mqc_cuest_scf
