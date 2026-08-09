!! Closed-shell SCF on the CPU
module mqc_libcint_rhf
   !! Restricted Hartree-Fock over libcint integrals, in core.
   !!
   !! The same algorithm the cuEST SCF runs -- canonical orthogonalisation,
   !! DIIS on the FDS - SDF commutator, energy from the density -- written
   !! against host arrays instead of device pointers. It is not a shared
   !! implementation, and that is deliberate: the two agreeing is worth more
   !! than the duplication costs, because agreement between two independent
   !! codes is the only check the GPU path has ever had.
   !!
   !! `mqc_diis` is shared, though. That one is already backend-neutral and
   !! there is nothing to learn from writing it twice.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_diis, only: diis_state_t, diis_slot_of_age, diis_coefficients
   use mqc_libcint_integrals, only: libcint_molecule_t
   implicit none
   private

   public :: rhf_result_t
   public :: run_libcint_rhf

   type :: rhf_result_t
      !! What a converged closed-shell SCF leaves behind
      real(dp) :: energy = 0.0_dp              !! Total, including nuclear repulsion
      real(dp) :: electronic = 0.0_dp          !! Without it
      real(dp) :: nuclear_repulsion = 0.0_dp
      integer :: iterations = 0
      logical :: converged = .false.
      real(dp), allocatable :: orbital_energies(:)
      real(dp), allocatable :: orbitals(:, :)
      real(dp), allocatable :: density(:, :)
      integer :: n_occupied = 0
   end type rhf_result_t

contains

   subroutine run_libcint_rhf(mol, nelec, max_iter, energy_tol, density_tol, &
                              verbose, result, error)
      !! Drive a closed-shell SCF to convergence
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: nelec
      integer, intent(in) :: max_iter
      real(dp), intent(in) :: energy_tol, density_tol
      logical, intent(in) :: verbose
      type(rhf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: s(:, :), h(:, :), eri(:, :, :, :)
      real(dp), allocatable :: x(:, :), fock(:, :), density(:, :), density_old(:, :)
      real(dp), allocatable :: coeff(:, :), eigenvalues(:)
      real(dp) :: e_elec, e_old, de, drms
      integer :: n_ao, n_mo, n_occ, iter

      if (mod(nelec, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "RHF needs an even electron count; this "// &
                        "system has an odd one and wants an unrestricted method")
         return
      end if

      n_ao = mol%nao
      n_occ = nelec/2
      if (n_occ < 1) then
         call error%set(ERROR_VALIDATION, "RHF: no electrons to place")
         return
      end if

      call mol%overlap(s)
      call mol%core_hamiltonian(h)
      call mol%eris(eri)

      call build_orthogonalizer(s, x, n_mo, error)
      if (error%has_error()) return
      if (n_occ > n_mo) then
         call error%set(ERROR_VALIDATION, "RHF: more occupied orbitals than the basis "// &
                        "supports after near-null modes were dropped")
         return
      end if

      allocate (fock(n_ao, n_ao), density(n_ao, n_ao), density_old(n_ao, n_ao))

      ! Core guess: F = H. Crude, and enough for the small systems this backend
      ! is for; a better guess is a convergence question, not a correctness one.
      fock = h
      call diagonalize(fock, x, n_ao, n_mo, coeff, eigenvalues, error)
      if (error%has_error()) return
      call build_density(coeff, n_occ, density)

      e_old = 0.0_dp
      result%converged = .false.

      do iter = 1, max_iter
         density_old = density
         call build_fock(h, eri, density, fock)
         e_elec = electronic_energy(h, fock, density)

         call diagonalize(fock, x, n_ao, n_mo, coeff, eigenvalues, error)
         if (error%has_error()) return
         call build_density(coeff, n_occ, density)

         de = abs(e_elec - e_old)
         drms = sqrt(sum((density - density_old)**2)/real(n_ao*n_ao, dp))
         if (verbose) then
            write (*, "(a,i4,a,f20.12,a,es10.3,a,es10.3)") &
               "  iter ", iter, "  E = ", e_elec + mol%nuclear_repulsion(), &
               "  dE = ", de, "  dD = ", drms
         end if

         e_old = e_elec
         result%iterations = iter
         if (iter > 1 .and. de < energy_tol .and. drms < density_tol) then
            result%converged = .true.
            exit
         end if
      end do

      ! The energy that goes out is the one belonging to the density that
      ! satisfied the test, so it is recomputed from the final Fock rather
      ! than carried over from the loop.
      call build_fock(h, eri, density, fock)
      result%electronic = electronic_energy(h, fock, density)
      result%nuclear_repulsion = mol%nuclear_repulsion()
      result%energy = result%electronic + result%nuclear_repulsion
      result%n_occupied = n_occ
      call move_alloc(eigenvalues, result%orbital_energies)
      call move_alloc(coeff, result%orbitals)
      call move_alloc(density, result%density)
   end subroutine run_libcint_rhf

   subroutine build_orthogonalizer(overlap, transform, n_mo, error)
      !! Canonical orthogonaliser X = U s^(-1/2), near-null modes dropped
      !!
      !! Canonical rather than symmetric so a basis with near-linear
      !! dependence loses the offending combinations instead of amplifying
      !! them, which is the same choice the cuEST path makes.
      real(dp), intent(in) :: overlap(:, :)
      real(dp), allocatable, intent(out) :: transform(:, :)
      integer, intent(out) :: n_mo
      type(error_t), intent(inout) :: error

      real(dp), parameter :: NULL_THRESHOLD = 1.0e-7_dp
      real(dp), allocatable :: vectors(:, :), values(:)
      integer :: n_ao, i, kept, info

      n_ao = size(overlap, 1)
      allocate (vectors(n_ao, n_ao), values(n_ao))
      vectors = overlap
      call pic_syev(vectors, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "RHF: overlap diagonalisation failed")
         return
      end if

      kept = count(values > NULL_THRESHOLD)
      if (kept == 0) then
         call error%set(ERROR_VALIDATION, "RHF: the overlap matrix is singular")
         return
      end if

      allocate (transform(n_ao, kept))
      n_mo = 0
      do i = 1, n_ao
         if (values(i) <= NULL_THRESHOLD) cycle
         n_mo = n_mo + 1
         transform(:, n_mo) = vectors(:, i)/sqrt(values(i))
      end do
   end subroutine build_orthogonalizer

   subroutine diagonalize(fock, x, n_ao, n_mo, coeff, eigenvalues, error)
      !! F' = X^T F X, diagonalise, and bring the orbitals back: C = X C'
      real(dp), intent(in) :: fock(:, :), x(:, :)
      integer, intent(in) :: n_ao, n_mo
      real(dp), allocatable, intent(inout) :: coeff(:, :), eigenvalues(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: work(:, :), f_ortho(:, :)
      integer :: info

      allocate (work(n_ao, n_mo), f_ortho(n_mo, n_mo))
      call pic_gemm(fock, x, work)                       ! F X
      call pic_gemm(x, work, f_ortho, transa="T")        ! X^T F X

      if (allocated(eigenvalues)) deallocate (eigenvalues)
      allocate (eigenvalues(n_mo))
      call pic_syev(f_ortho, eigenvalues, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "RHF: Fock diagonalisation failed")
         return
      end if

      if (allocated(coeff)) deallocate (coeff)
      allocate (coeff(n_ao, n_mo))
      call pic_gemm(x, f_ortho, coeff)                   ! C = X C'
   end subroutine diagonalize

   pure subroutine build_density(coeff, n_occ, density)
      !! D = 2 C_occ C_occ^T, the closed-shell density
      real(dp), intent(in) :: coeff(:, :)
      integer, intent(in) :: n_occ
      real(dp), intent(out) :: density(:, :)

      integer :: mu, nu, i

      density = 0.0_dp
      do nu = 1, size(density, 2)
         do mu = 1, size(density, 1)
            do i = 1, n_occ
               density(mu, nu) = density(mu, nu) + 2.0_dp*coeff(mu, i)*coeff(nu, i)
            end do
         end do
      end do
   end subroutine build_density

   pure subroutine build_fock(h, eri, density, fock)
      !! F = H + J - K/2, with J and K contracted straight from the ERIs
      real(dp), intent(in) :: h(:, :), eri(:, :, :, :), density(:, :)
      real(dp), intent(out) :: fock(:, :)

      integer :: mu, nu, la, si, n

      n = size(h, 1)
      fock = h
      do nu = 1, n
         do mu = 1, n
            do si = 1, n
               do la = 1, n
                  ! (mu nu | la si) is Coulomb, (mu la | nu si) is exchange,
                  ! and the half is the closed-shell factor -- D already
                  ! carries the two electrons per orbital.
                  fock(mu, nu) = fock(mu, nu) + density(la, si) &
                                 *(eri(mu, nu, la, si) - 0.5_dp*eri(mu, la, nu, si))
               end do
            end do
         end do
      end do
   end subroutine build_fock

   pure function electronic_energy(h, fock, density) result(energy)
      !! E = 1/2 sum_uv D_uv (H_uv + F_uv)
      real(dp), intent(in) :: h(:, :), fock(:, :), density(:, :)
      real(dp) :: energy

      energy = 0.5_dp*sum(density*(h + fock))
   end function electronic_energy

end module mqc_libcint_rhf
