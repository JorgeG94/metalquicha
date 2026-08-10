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
   use mqc_diis, only: diis_state_t
   use mqc_libcint_direct, only: schwarz_bounds, build_fock_direct, direct_stats_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_df_tensor
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
                              verbose, result, error, aux, diis_vectors, in_core)
      !! Drive a closed-shell SCF to convergence
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: nelec
      integer, intent(in) :: max_iter
      real(dp), intent(in) :: energy_tol, density_tol
      logical, intent(in) :: verbose
      type(rhf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      type(libcint_molecule_t), intent(in), optional :: aux
         !! Auxiliary basis. Present means density-fitted J and K, which is
         !! what cuEST always does -- so a comparison against the GPU path
         !! needs this, and a comparison against an exact-ERI code needs it
         !! absent. Both are worth being able to ask for.
      integer, intent(in), optional :: diis_vectors
         !! Subspace size. Zero turns DIIS off, which is worth being able to
         !! do: the two paths should agree with and without it, and if they
         !! do not, the extrapolation is where to look.
      logical, intent(in), optional :: in_core
         !! Store every integral and contract from the tensor, instead of
         !! rebuilding the Fock matrix directly each iteration. Default is
         !! direct. This exists because the in-core path is the one validated
         !! against PySCF, so it is what the direct build is checked against --
         !! not because anything should run with it. It is n^4 in memory.

      integer :: diis_size
      logical :: use_in_core
      real(dp), allocatable :: bounds(:, :)
      type(direct_stats_t) :: stats

      real(dp), allocatable :: s(:, :), h(:, :), eri(:, :, :, :), bmat(:, :)
      real(dp), allocatable :: x(:, :), fock(:, :), density(:, :), density_old(:, :)
      real(dp), allocatable :: coeff(:, :), eigenvalues(:)
      real(dp), allocatable :: err(:, :), fock_flat(:)
      type(diis_state_t) :: diis
      logical :: extrapolated
      real(dp) :: e_elec, e_old, de, drms
      integer :: n_ao, n_mo, n_occ, iter

      if (mod(nelec, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "RHF needs an even electron count; this "// &
                        "system has an odd one and wants an unrestricted method")
         return
      end if

      diis_size = 8
      if (present(diis_vectors)) diis_size = diis_vectors
      use_in_core = .false.
      if (present(in_core)) use_in_core = in_core

      n_ao = mol%nao
      n_occ = nelec/2
      if (n_occ < 1) then
         call error%set(ERROR_VALIDATION, "RHF: no electrons to place")
         return
      end if

      ! Reported because this is what went unnoticed for as long as it did:
      ! 6-31G* read as spherical converges just as prettily as 6-31G* read
      ! Cartesian -- one basis function short and 1.4 mHartree out -- and
      ! neither the iterations nor the final energy look wrong on their own.
      ! The count and the angular form together are what name the basis.
      if (verbose) then
         if (mol%cartesian) then
            write (*, "(a,i0,a)") "  basis functions: ", n_ao, "  (Cartesian, 6d/10f)"
         else
            write (*, "(a,i0,a)") "  basis functions: ", n_ao, "  (spherical, 5d/7f)"
         end if
      end if

      call mol%overlap(s)
      call mol%core_hamiltonian(h)
      if (present(aux)) then
         call build_df_tensor(mol, aux, bmat, error)
         if (error%has_error()) return
      else if (use_in_core) then
         call mol%eris(eri)
      else
         ! The bounds depend on the basis and the geometry, not the density, so
         ! one set serves every iteration of the SCF.
         call schwarz_bounds(mol, bounds, error)
         if (error%has_error()) return
      end if

      call build_orthogonalizer(s, x, n_mo, error)
      if (error%has_error()) return
      if (n_occ > n_mo) then
         call error%set(ERROR_VALIDATION, "RHF: more occupied orbitals than the basis "// &
                        "supports after near-null modes were dropped")
         return
      end if

      allocate (fock(n_ao, n_ao), density(n_ao, n_ao), density_old(n_ao, n_ao))
      allocate (err(n_mo, n_mo), fock_flat(n_ao*n_ao))

      ! The error vector lives in the orthogonal basis, where it is n_mo
      ! square rather than n_ao -- that is the same shape the cuEST path uses
      ! and the reason both converge alike.
      call diis%init(diis_size, n_ao*n_ao, n_mo*n_mo)

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
         if (allocated(bmat)) then
            call build_fock_df(h, bmat, density, coeff, n_occ, fock)
         else if (allocated(eri)) then
            call build_fock(h, eri, density, fock)
         else
            call build_fock_direct(mol, h, density, bounds, fock, stats, error)
            if (error%has_error()) return
         end if

         ! The energy belongs to the Fock built from this density, so it is
         ! taken before extrapolation. A DIIS-mixed Fock is a convergence
         ! device, not a state anything is the energy of.
         e_elec = electronic_energy(h, fock, density)

         ! FDS - SDF, projected into the orthogonal basis. It vanishes exactly
         ! when F and D commute, which is what convergence means, so it is the
         ! quantity worth extrapolating against.
         call commutator(fock, density, s, x, err)
         fock_flat = reshape(fock, [n_ao*n_ao])
         call diis%push(fock_flat, reshape(err, [n_mo*n_mo]))
         call diis%extrapolate(fock_flat, extrapolated)
         if (extrapolated) fock = reshape(fock_flat, [n_ao, n_ao])

         call diagonalize(fock, x, n_ao, n_mo, coeff, eigenvalues, error)
         if (error%has_error()) return
         call build_density(coeff, n_occ, density)

         de = abs(e_elec - e_old)
         drms = sqrt(sum((density - density_old)**2)/real(n_ao*n_ao, dp))
         if (verbose) then
            write (*, "(a,i4,a,f20.12,a,es10.3,a,es10.3,a,i0)") &
               "  iter ", iter, "  E = ", e_elec + mol%nuclear_repulsion(), &
               "  dE = ", de, "  dD = ", drms, "  diis ", diis%count()
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
      if (allocated(bmat)) then
         call build_fock_df(h, bmat, density, coeff, n_occ, fock)
      else if (allocated(eri)) then
         call build_fock(h, eri, density, fock)
      else
         call build_fock_direct(mol, h, density, bounds, fock, stats, error)
         if (error%has_error()) return
      end if
      result%electronic = electronic_energy(h, fock, density)
      result%nuclear_repulsion = mol%nuclear_repulsion()
      result%energy = result%electronic + result%nuclear_repulsion
      result%n_occupied = n_occ
      call move_alloc(eigenvalues, result%orbital_energies)
      call move_alloc(coeff, result%orbitals)
      call move_alloc(density, result%density)
      call diis%destroy()
   end subroutine run_libcint_rhf

   subroutine commutator(fock, density, overlap, x, err)
      !! e = X^T (F D S - S D F) X, the DIIS error vector
      !!
      !! Projected into the orthogonal basis rather than left in the AO one:
      !! in the AO basis the commutator's size depends on the overlap's
      !! conditioning, so two bases describing the same molecule would
      !! converge differently for no physical reason.
      real(dp), intent(in) :: fock(:, :), density(:, :), overlap(:, :), x(:, :)
      real(dp), intent(out) :: err(:, :)

      real(dp), allocatable :: fd(:, :), fds(:, :), sd(:, :), sdf(:, :), work(:, :)
      integer :: n_ao, n_mo

      n_ao = size(fock, 1)
      n_mo = size(x, 2)
      allocate (fd(n_ao, n_ao), fds(n_ao, n_ao), sd(n_ao, n_ao), sdf(n_ao, n_ao))
      allocate (work(n_ao, n_mo))

      call pic_gemm(fock, density, fd)
      call pic_gemm(fd, overlap, fds)
      call pic_gemm(overlap, density, sd)
      call pic_gemm(sd, fock, sdf)
      fds = fds - sdf

      call pic_gemm(fds, x, work)
      call pic_gemm(x, work, err, transa="T")
   end subroutine commutator

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

   subroutine build_fock_df(h, b, density, coeff, n_occ, fock)
      !! F = H + J - K/2 from the fitted tensor rather than the exact ERIs
      !!
      !! Neither term ever forms a four-index object. Density fitting is not the
      !! opposite of a direct build -- it removes the four-index integrals
      !! altogether, so the direct/stored distinction does not apply to it. What
      !! is stored is B, which is n^2 by n_aux: 40 MB where the exact tensor is
      !! 800 MB, and n^3 rather than n^4.
      !!
      !! J is two contractions with the density: c_P = sum_uv B(uv,P) D_uv, then
      !! J_uv = sum_P B(uv,P) c_P.
      !!
      !! K goes through the occupied orbitals rather than the density. Writing
      !! D = 2 sum_i C_ui C_vi and substituting,
      !!
      !!    K_uv = sum_P sum_ls B(ul,P) D_ls B(sv,P)
      !!         = 2 sum_P sum_i W(u,i,P) W(v,i,P),   W^P = B^P C_occ
      !!
      !! which costs 2 n^2 n_occ per auxiliary function instead of the 2 n^3 of
      !! contracting the full density. The saving is a factor of n/n_occ, and it
      !! grows with basis size at fixed electron count -- which is the regime
      !! density fitting is used in. On water/cc-pVDZ that is already 4.8x.
      real(dp), intent(in) :: h(:, :), b(:, :), density(:, :)
      real(dp), intent(in) :: coeff(:, :)   !! MO coefficients; only the occupied block is read
      integer, intent(in) :: n_occ
      real(dp), intent(out) :: fock(:, :)

      real(dp), allocatable :: c(:), j(:, :), k(:, :), w(:, :), c_occ(:, :)
      integer :: n, naux, p

      n = size(h, 1)
      naux = size(b, 2)
      allocate (c(naux), j(n, n), k(n, n), w(n, n_occ), c_occ(n, n_occ))

      c_occ = coeff(:, 1:n_occ)

      do p = 1, naux
         c(p) = sum(reshape(b(:, p), [n, n])*density)
      end do

      j = 0.0_dp
      do p = 1, naux
         j = j + c(p)*reshape(b(:, p), [n, n])
      end do

      k = 0.0_dp
      do p = 1, naux
         ! `associate` gives a view of the column rather than a copy: b(:,p) is
         ! contiguous, so reshaping it is free as long as nothing forces a
         ! temporary, and going through gemm keeps it out of matmul's hands.
         associate (b_p => reshape(b(:, p), [n, n]))
            w = 0.0_dp
            call pic_gemm(b_p, c_occ, w)
            call pic_gemm(w, w, k, transb="T", alpha=2.0_dp, beta=1.0_dp)
         end associate
      end do

      fock = h + j - 0.5_dp*k
   end subroutine build_fock_df

   pure function electronic_energy(h, fock, density) result(energy)
      !! E = 1/2 sum_uv D_uv (H_uv + F_uv)
      real(dp), intent(in) :: h(:, :), fock(:, :), density(:, :)
      real(dp) :: energy

      energy = 0.5_dp*sum(density*(h + fock))
   end function electronic_energy

end module mqc_libcint_rhf
