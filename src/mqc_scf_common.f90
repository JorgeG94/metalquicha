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
   use pic_logger, only: logger => global_logger
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: build_orthogonalizer
   public :: report_linear_dependence
   public :: build_density_closed_shell
   public :: build_density_spin
   public :: spin_contamination
   public :: LINEAR_DEPENDENCE_TOL
   public :: GWH_K

   !> Overlap eigenvalues at or below this are dropped as linearly dependent.
   !>
   !> One number for both backends. It used to be `NULL_THRESHOLD` on the
   !> libcint side and `LINEAR_DEPENDENCE_TOL` on the cuEST side, both 1.0e-7,
   !> which meant the two paths agreed only by coincidence.
   real(dp), parameter :: LINEAR_DEPENDENCE_TOL = 1.0e-7_dp

   !> Overlap eigenvalues above the drop threshold but below this one are kept
   !> and said out loud.
   !>
   !> The zone between the two is the one that hurts. A mode below
   !> `LINEAR_DEPENDENCE_TOL` is discarded and the basis is smaller, which is
   !> a decision with a number attached; a mode just above it is *retained*,
   !> and X carries a 1/sqrt(s) of ten thousand or more into every iteration.
   !> Nothing fails. DIIS extrapolates along a direction that is mostly noise,
   !> convergence slows or stalls, and two runs of the same molecule in
   !> slightly different orientations can settle on different solutions.
   !> Diffuse functions on a large or compact system are how a basis gets
   !> here, which is precisely when nobody is watching for it.
   real(dp), parameter :: LINEAR_DEPENDENCE_WARN_TOL = 1.0e-5_dp
   public :: LINEAR_DEPENDENCE_WARN_TOL

   !> The empirical scale on the GWH off-diagonal. 1.75 is the value in
   !> universal use, and both backends have to start from the same matrix or
   !> comparing their iteration counts means nothing -- which is exactly why it
   !> is one constant here rather than one in each of them.
   real(dp), parameter :: GWH_K = 1.75_dp

contains

   subroutine build_orthogonalizer(overlap, transform, n_mo, error, &
                                   smallest_overlap, smallest_kept, threshold)
      !! Canonical orthogonaliser X = U s^(-1/2), near-null modes dropped
      !!
      !! Canonical rather than symmetric so a basis with near-linear dependence
      !! -- diffuse functions, or two fragments close together -- loses the
      !! offending combinations instead of amplifying them.
      real(dp), intent(in) :: overlap(:, :)                  !! S, (n_ao, n_ao)
      real(dp), allocatable, intent(out) :: transform(:, :)  !! X, (n_ao, n_mo)
      integer, intent(out) :: n_mo                           !! Surviving orbitals
      type(error_t), intent(inout) :: error
      real(dp), intent(out), optional :: smallest_overlap
         !! The smallest eigenvalue of S, dropped or not -- the conditioning of
         !! the basis as given, before this routine does anything about it.
      real(dp), intent(out), optional :: smallest_kept
         !! The smallest eigenvalue that survived, which is what X actually
         !! divides by and therefore what the SCF has to live with.
      real(dp), intent(in), optional :: threshold
         !! Overlap eigenvalues at or below this are dropped. Defaults to
         !! `LINEAR_DEPENDENCE_TOL`, and `keywords.scf.linear_dependence_threshold`
         !! is how a deck overrides it. A non-positive value is ignored rather
         !! than obeyed: zero would keep every mode including the numerically
         !! null ones, which is not a basis anyone means to ask for.
         !!
         !! Optional, all of them, so the call sites that only want X are
         !! unchanged -- there are five across two backends and only the ones
         !! that report need the numbers.

      real(dp), allocatable :: eigenvectors(:, :), eigenvalues(:)
      integer :: n_ao, i, kept, info
      real(dp) :: tol

      tol = LINEAR_DEPENDENCE_TOL
      if (present(threshold)) then
         if (threshold > 0.0_dp) tol = threshold
      end if

      if (present(smallest_overlap)) smallest_overlap = 0.0_dp
      if (present(smallest_kept)) smallest_kept = 0.0_dp

      n_ao = size(overlap, 1)
      allocate (eigenvectors(n_ao, n_ao), eigenvalues(n_ao))
      eigenvectors = overlap

      call pic_syev(eigenvectors, eigenvalues, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "SCF: overlap matrix diagonalization failed")
         return
      end if

      ! pic_syev returns eigenvalues ascending, so the discarded ones lead.
      if (present(smallest_overlap)) smallest_overlap = eigenvalues(1)
      n_mo = count(eigenvalues > tol)
      if (present(smallest_kept) .and. n_mo > 0) then
         smallest_kept = eigenvalues(n_ao - n_mo + 1)
      end if
      if (n_mo == 0) then
         call error%set(ERROR_VALIDATION, "SCF: overlap matrix is singular")
         return
      end if

      allocate (transform(n_ao, n_mo))
      kept = 0
      do i = 1, n_ao
         if (eigenvalues(i) <= tol) cycle
         kept = kept + 1
         transform(:, kept) = eigenvectors(:, i)/sqrt(eigenvalues(i))
      end do
   end subroutine build_orthogonalizer

   subroutine report_linear_dependence(n_ao, n_mo, smallest_overlap, smallest_kept, verbose, &
                                       threshold)
      !! Say what the orthogonaliser did to the basis, and whether to worry
      !!
      !! **The warnings are not gated on `verbose`, deliberately.** Dropping a
      !! basis function changes the space every method downstream works in:
      !! the SCF, the virtuals an MP2 or CC run correlates, the count in the
      !! output file. A quiet run that silently solves a different problem
      !! than the one asked for is the failure this exists to prevent, and
      !! `logger%info` on its own would hide it at exactly the verbosity most
      !! production runs use. Only the healthy-basis line is verbose-only.
      integer, intent(in) :: n_ao      !! Basis functions the basis set defines
      integer, intent(in) :: n_mo      !! Orbitals that survived orthogonalisation
      real(dp), intent(in) :: smallest_overlap
      real(dp), intent(in) :: smallest_kept
      logical, intent(in) :: verbose
      real(dp), intent(in), optional :: threshold
         !! The cutoff actually applied, which the message quotes. Passed
         !! rather than read from the parameter because a deck can move it,
         !! and a warning naming a cutoff the run did not use is worse than
         !! one that names none.

      character(len=160) :: line
      integer :: n_dropped
      real(dp) :: tol

      tol = LINEAR_DEPENDENCE_TOL
      if (present(threshold)) then
         if (threshold > 0.0_dp) tol = threshold
      end if

      n_dropped = n_ao - n_mo

      if (n_dropped > 0) then
         call logger%warning("")
         write (line, "(a,i0,a,i0,a)") "  LINEAR DEPENDENCE: ", n_dropped, " of ", n_ao, &
            " basis functions were dropped."
         call logger%warning(trim(line))
         write (line, "(a,es9.2,a,es9.2,a)") "     the overlap's smallest eigenvalue is ", &
            smallest_overlap, ", below the ", tol, " cutoff."
         call logger%warning(trim(line))
         write (line, "(a,i0,a)") "     the SCF and everything after it run in ", n_mo, &
            " orbitals, not "//trim(adjustl(int_text(n_ao)))//"."
         call logger%warning(trim(line))
         call logger%warning("     This is usually diffuse functions overlapping too "// &
                             "well to tell apart, and")
         call logger%warning("     dropping them is the right repair -- but the basis "// &
                             "is no longer the one")
         call logger%warning("     named in the input, so energies are not comparable "// &
                             "with a run that kept")
         call logger%warning("     them all. A smaller or less diffuse basis is the fix "// &
                             "if that matters.")
         call logger%warning("")
      else if (smallest_kept > 0.0_dp .and. &
               smallest_kept < max(LINEAR_DEPENDENCE_WARN_TOL, 100.0_dp*tol)) then
         ! Kept, and worth saying so. See LINEAR_DEPENDENCE_WARN_TOL.
         call logger%warning("")
         write (line, "(a,es9.2,a)") "  NEARLY LINEARLY DEPENDENT: the overlap's "// &
            "smallest eigenvalue is ", smallest_kept, "."
         call logger%warning(trim(line))
         call logger%warning("     Nothing was dropped, so the basis is intact, but X "// &
                             "carries a large")
         call logger%warning("     1/sqrt(s) into every iteration. Expect slow or "// &
                             "stalled convergence, and")
         call logger%warning("     treat a converged energy with care: an SCF this "// &
                             "ill-conditioned can settle")
         call logger%warning("     on different solutions from different guesses.")
         call logger%warning("")
      else if (verbose) then
         write (line, "(a,es9.2)") "  overlap: smallest eigenvalue ", smallest_overlap
         call logger%info(trim(line))
      end if
   end subroutine report_linear_dependence

   pure function int_text(value) result(text)
      !! An integer as a trimmed string, for the messages above
      integer, intent(in) :: value
      character(len=:), allocatable :: text
      character(len=12) :: buffer

      write (buffer, "(i0)") value
      text = trim(buffer)
   end function int_text

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
