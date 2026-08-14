!! Linear algebra every SCF backend needs, written once
module mqc_scf_common
   !! The parts of an SCF that do not depend on how the integrals were made.
   !!
   !! Orthogonalising the overlap, forming a density from occupied orbitals and
   !! measuring spin contamination are all pure linear algebra: the libcint and
   !! cuEST paths had line-for-line copies of each, down to the same tolerance
   !! declared under two different names. They live here so the two backends
   !! cannot drift apart on the numerics while agreeing on the physics.
   !!
   !! What is *not* here is anything that touches integrals, DIIS or the
   !! iteration itself -- those differ between the backends for real reasons.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: build_orthogonalizer
   public :: build_density_closed_shell
   public :: build_density_spin
   public :: spin_contamination
   public :: LINEAR_DEPENDENCE_TOL

   !> Overlap eigenvalues at or below this are dropped as linearly dependent.
   !>
   !> One number for both backends. It used to be `NULL_THRESHOLD` on the
   !> libcint side and `LINEAR_DEPENDENCE_TOL` on the cuEST side, both 1.0e-7,
   !> which meant the two paths agreed only by coincidence.
   real(dp), parameter :: LINEAR_DEPENDENCE_TOL = 1.0e-7_dp

contains

   subroutine build_orthogonalizer(overlap, transform, n_mo, error)
      !! Canonical orthogonaliser X = U s^(-1/2), near-null modes dropped
      !!
      !! Canonical rather than symmetric so a basis with near-linear dependence
      !! -- diffuse functions, or two fragments close together -- loses the
      !! offending combinations instead of amplifying them.
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

   subroutine build_density_closed_shell(coeff, n_occ, density)
      !! D = 2 C_occ C_occ^T, the closed-shell density
      real(dp), intent(in) :: coeff(:, :)        !! C, (n_ao, n_mo)
      integer, intent(in) :: n_occ               !! Doubly occupied orbitals
      real(dp), intent(inout) :: density(:, :)   !! D, (n_ao, n_ao)

      density = 0.0_dp
      if (n_occ <= 0) return
      call pic_gemm(coeff(:, 1:n_occ), coeff(:, 1:n_occ), density, &
                    transb="T", alpha=2.0_dp, beta=0.0_dp)
   end subroutine build_density_closed_shell

   subroutine build_density_spin(coeff, n_occ, density)
      !! D_sigma = C_occ C_occ^T, one electron per occupied orbital
      !!
      !! Deliberately a separate routine rather than the closed-shell one with a
      !! factor: that one carries the two electrons per spatial orbital, and
      !! reusing it here would double every spin density. That error converges,
      !! which is what makes it worth keeping the two names apart.
      real(dp), intent(in) :: coeff(:, :)        !! C_sigma, (n_ao, n_mo)
      integer, intent(in) :: n_occ               !! Occupied orbitals of this spin
      real(dp), intent(inout) :: density(:, :)   !! D_sigma, (n_ao, n_ao)

      density = 0.0_dp
      if (n_occ <= 0) return
      call pic_gemm(coeff(:, 1:n_occ), coeff(:, 1:n_occ), density, &
                    transb="T", alpha=1.0_dp, beta=0.0_dp)
   end subroutine build_density_spin

   function spin_contamination(occ_alpha, occ_beta, overlap, n_alpha, n_beta) result(s_squared)
      !! <S^2> for a UHF/UKS determinant
      !!
      !!   <S^2> = S_z(S_z+1) + n_beta - sum_ij |<phi_i^a|phi_j^b>|^2
      !!
      !! The exact value for a pure spin state is S_z(S_z+1); the excess is
      !! spin contamination, and reporting it is the cheapest way to notice
      !! that an open-shell answer is not describing the state intended. A UHF
      !! solution that has collapsed onto the restricted one, or onto a badly
      !! broken-symmetry state, says so here and nowhere else.
      real(dp), intent(in) :: occ_alpha(:, :), occ_beta(:, :), overlap(:, :)
      integer, intent(in) :: n_alpha, n_beta
      real(dp) :: s_squared

      real(dp), allocatable :: scratch(:, :), mo_overlap(:, :)
      real(dp) :: sz

      sz = 0.5_dp*real(n_alpha - n_beta, dp)
      s_squared = sz*(sz + 1.0_dp) + real(n_beta, dp)
      if (n_alpha == 0 .or. n_beta == 0) return

      ! S C_beta, then C_alpha^T (S C_beta): the alpha-beta MO overlap block.
      allocate (scratch(size(overlap, 1), n_beta), mo_overlap(n_alpha, n_beta))
      call pic_gemm(overlap, occ_beta(:, 1:n_beta), scratch)
      call pic_gemm(occ_alpha(:, 1:n_alpha), scratch, mo_overlap, transa="T")
      s_squared = s_squared - sum(mo_overlap**2)
      deallocate (scratch, mo_overlap)
   end function spin_contamination

end module mqc_scf_common
