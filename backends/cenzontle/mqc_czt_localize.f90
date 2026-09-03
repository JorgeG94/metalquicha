!! Orbital localization
module mqc_czt_localize
   !! Foster-Boys localization of the occupied orbitals, by Jacobi sweeps.
   !!
   !! Boys maximizes the sum of squared orbital centroids,
   !!
   !!     L = sum_i sum_k <i| r_k |i>^2
   !!
   !! which is the same thing as minimizing the sum of orbital spreads, since
   !! `sum_i <i|r^2|i>` is invariant under a rotation among the occupied
   !! orbitals. Only the centroids move, so the dipole integrals from
   !! `mqc_czt_multipole` are the whole input.
   !!
   !! **Why localize at all here:** an effective fragment potential is built from
   !! localized quantities -- multipoles distributed over bonds and lone pairs,
   !! one polarizability tensor per localized orbital, exchange repulsion from
   !! localized-orbital overlaps. GAMESS's MAKEFP defaults to Edmiston-Ruedenberg
   !! rather than Boys and that choice changes every distributed parameter, so a
   !! comparison against it has to name the localization; `LOCAL=BOYS` is
   !! settable there.
   !!
   !! **The rotation.** For a pair (i,j) rotated by gamma, and writing
   !! `u_k = (d_ii - d_jj)/2`, `v_k = d_ij`, `w_k = (d_ii + d_jj)/2` for each
   !! Cartesian component of the dipole matrix in the current basis,
   !!
   !!     d_i'i' = w + u cos2g + v sin2g,    d_j'j' = w - (u cos2g + v sin2g)
   !!
   !! so the pair's contribution is `2 sum_k w_k^2 + 2 sum_k p_k^2` with
   !! `p_k = u_k cos2g + v_k sin2g`. The first term does not move, and expanding
   !! the second gives a pure sinusoid in 4g:
   !!
   !!     sum_k p_k^2 = const + P cos4g + Q sin4g,
   !!     P = (sum_k u_k^2 - sum_k v_k^2)/2,   Q = sum_k u_k v_k
   !!
   !! maximized at `gamma = atan2(Q, P)/4`, with the gain `sqrt(P^2+Q^2) - P`,
   !! which is non-negative for every pair and is what makes the sweep monotone.
   !! The test asserts that rather than trusting the derivation: the sign
   !! conventions in the literature differ, and a wrong one still converges, to a
   !! stationary point that is not the maximum.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_multipole, only: multipole_matrices
   implicit none
   private

   public :: boys_localize

   ! Sweeps before giving up. Boys on a fragment converges in a handful.
   integer, parameter :: DEFAULT_MAX_SWEEPS = 200

   ! Convergence on the largest rotation angle in a sweep, radians. On the angle
   ! rather than on the functional: near the maximum the functional is quadratic
   ! in the angle, so a threshold on it stops while the orbitals are still moving
   ! at the square root of that threshold.
   real(dp), parameter :: DEFAULT_ANGLE_TOL = 1.0e-10_dp

contains

   subroutine boys_localize(mol, coefficients, n_occ, localized, centroids, error, &
                            max_sweeps, angle_tol, sweeps_taken, functional)
      !! Localize the occupied orbitals and report where they sit
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coefficients(:, :)   !! MO coefficients, (n_ao, n_mo)
      integer, intent(in) :: n_occ                 !! Occupied orbitals to localize
      real(dp), allocatable, intent(out) :: localized(:, :)
         !! The localized occupied orbitals, (n_ao, n_occ).
      real(dp), allocatable, intent(out) :: centroids(:, :)
         !! `<i| r |i>` per localized orbital, (3, n_occ), Bohr. These are the
         !! points a distributed polarizability is placed on.
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_sweeps
      real(dp), intent(in), optional :: angle_tol
      integer, intent(out), optional :: sweeps_taken
      real(dp), intent(out), optional :: functional
         !! The converged value of L, for a caller that wants to check it rose.

      real(dp), allocatable :: dip(:, :, :), d(:, :, :), rot(:, :), work(:, :)
      real(dp) :: u(3), v(3)
      real(dp) :: pp, qq, gamma, biggest, cg, sg, tol
      real(dp) :: di, dj
      integer :: k, i, j, m, sweep, limit
      character(len=16) :: text
      ! TODO(mqc): `rot` is declared and never used; the rotation is applied in
      ! place below.

      if (n_occ < 1 .or. n_occ > size(coefficients, 2)) then
         write (text, "(i0)") n_occ
         call error%set(ERROR_VALIDATION, "cannot localize "//trim(text)// &
                        " orbitals: that is not a subset of the ones supplied")
         return
      end if

      limit = DEFAULT_MAX_SWEEPS
      if (present(max_sweeps)) limit = max_sweeps
      tol = DEFAULT_ANGLE_TOL
      if (present(angle_tol)) tol = angle_tol

      ! The *functional* is origin-dependent but the localization is not:
      ! shifting the origin adds a constant to every centroid and rotates
      ! nothing, since `sum_i d_ii` is invariant.
      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, error)
      if (error%has_error()) return

      allocate (localized(mol%nao, n_occ))
      localized = coefficients(:, 1:n_occ)

      ! The dipole matrices in the occupied MO basis. Rotations act on these
      ! directly from here on, which keeps a sweep O(n_occ^2) rather than
      ! re-transforming from the AO basis for every pair.
      allocate (d(n_occ, n_occ, 3), work(mol%nao, n_occ))
      do k = 1, 3
         call pic_gemm(dip(:, :, k), localized, work)
         call pic_gemm(localized, work, d(:, :, k), transa="T")
      end do

      do sweep = 1, limit
         biggest = 0.0_dp
         do i = 1, n_occ - 1
            do j = i + 1, n_occ
               do k = 1, 3
                  u(k) = 0.5_dp*(d(i, i, k) - d(j, j, k))
                  v(k) = d(i, j, k)
               end do
               pp = 0.5_dp*(sum(u*u) - sum(v*v))
               qq = sum(u*v)
               ! Nothing to gain from a pair whose contribution is already
               ! stationary; atan2(0,0) is also undefined.
               if (abs(pp) < tiny(1.0_dp) .and. abs(qq) < tiny(1.0_dp)) cycle
               gamma = 0.25_dp*atan2(qq, pp)
               if (abs(gamma) <= tol) cycle
               biggest = max(biggest, abs(gamma))

               cg = cos(gamma)
               sg = sin(gamma)
               ! Rotate the two orbitals, and the two rows and columns of every
               ! dipole matrix that they index.
               do m = 1, mol%nao
                  di = localized(m, i)
                  dj = localized(m, j)
                  localized(m, i) = cg*di + sg*dj
                  localized(m, j) = -sg*di + cg*dj
               end do
               do k = 1, 3
                  do m = 1, n_occ
                     di = d(m, i, k)
                     dj = d(m, j, k)
                     d(m, i, k) = cg*di + sg*dj
                     d(m, j, k) = -sg*di + cg*dj
                  end do
                  do m = 1, n_occ
                     di = d(i, m, k)
                     dj = d(j, m, k)
                     d(i, m, k) = cg*di + sg*dj
                     d(j, m, k) = -sg*di + cg*dj
                  end do
               end do
            end do
         end do
         if (biggest <= tol) exit
      end do

      if (present(sweeps_taken)) sweeps_taken = min(sweep, limit)

      allocate (centroids(3, n_occ))
      do i = 1, n_occ
         do k = 1, 3
            centroids(k, i) = d(i, i, k)
         end do
      end do

      if (present(functional)) then
         functional = 0.0_dp
         do i = 1, n_occ
            functional = functional + sum(centroids(:, i)**2)
         end do
      end if

      deallocate (dip, d, work)
   end subroutine boys_localize

end module mqc_czt_localize
