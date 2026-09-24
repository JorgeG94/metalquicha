!! Orbital localization
module mqc_czt_localize
   !! Foster-Boys and Edmiston-Ruedenberg localization of the occupied
   !! orbitals, both by Jacobi sweeps of two-by-two rotations.
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
   !! **The Boys rotation.** For a pair (i,j) rotated by gamma, and writing
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
   !!
   !! **Edmiston-Ruedenberg** maximizes the orbital self-repulsion
   !!
   !!     D = sum_i (ii|ii)
   !!
   !! equivalently minimizes the interorbital Coulomb and exchange energy, since
   !! the sum of all of it is invariant. For a pair rotated as
   !! `i' = c i + s j`, `j' = -s i + c j` the change is exactly
   !! `A (1 - cos4g) + B sin4g`, with
   !!
   !!     A = (ij|ij) - [(ii|ii) + (jj|jj) - 2 (ii|jj)]/4,   B = (ii|ij) - (jj|ij)
   !!
   !! so the pair's maximum is at `gamma = atan2(B, -A)/4` and gains
   !! `A + sqrt(A^2 + B^2)`, again never negative. That is the Edmiston and
   !! Ruedenberg angle [Rev. Mod. Phys. 35, 457 (1963)] and GAMESS's `LOCROT`,
   !! whose `ALM` is -A. ER needs the occupied two-electron integrals rather
   !! than dipoles; `occupied_eri` makes them without the AO tensor, and each
   !! rotation carries them along, so no integral is recomputed in the sweeps.
   !!
   !! ER starts from the orbitals it is given, as GAMESS's does from the
   !! canonical ones (`LOCSET` and `LOCLIZ` in `local.src`: canonical order,
   !! identity transformation). It has several maxima, so the start and the
   !! sweep schedule decide which is reached; both follow GAMESS here, and a
   !! Boys start is offered for a caller that wants one. A pair whose functional
   !! is flat in its angle, `A = B = 0`, is left alone.
   use pic_types, only: dp, int64
   use pic_blas_interfaces, only: pic_gemm
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, two_electron_block, two_electron_optimizer, &
                                eri_shell_table_t, eri_shell_table, eri_schwarz_collapse, &
                                pair_index
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_multipole, only: multipole_matrices
   use libcint_fortran, only: libcint_del_optimizer
   implicit none
   private

   public :: boys_localize
   public :: er_localize
   public :: occupied_eri
   public :: LOCALIZER_BOYS
   public :: LOCALIZER_ER

   character(len=*), parameter :: LOCALIZER_BOYS = "boys"
      !! Foster-Boys, `boys_localize`
   character(len=*), parameter :: LOCALIZER_ER = "er"
      !! Edmiston-Ruedenberg, `er_localize`

   real(dp), parameter :: ER_SCREEN_TOL = 1.0e-12_dp
      !! Shell quartets whose Schwarz bound falls below this are left out of
      !! the occupied transform.
   integer, parameter :: ER_MAX_SWEEPS = 2000
      !! ER's sweep limit, ten times Boys'. Jacobi converges linearly, and
      !! slowly along the flat directions a pi system has: adenine in 6-31G*
      !! takes 354 sweeps to 1e-10. A sweep over the stored integrals is
      !! O(n_occ^5) and cheap -- those 354 cost half a second against 18 for
      !! the integrals -- so the limit is there to stop a cycle, not to save
      !! time.
   integer(int64), parameter :: ER_TRANSFORM_WORDS = 2_int64**27
      !! Doubles the half-transformed `(mn|kl)` may hold at once, 1 GiB. Past
      !! it the `kl` range is batched and the integrals are recomputed per
      !! batch.

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

   subroutine er_localize(mol, coefficients, n_occ, localized, centroids, error, &
                          max_sweeps, angle_tol, sweeps_taken, functional, guess, &
                          converged)
      !! Edmiston-Ruedenberg localization of the occupied orbitals, by Jacobi sweeps
      !!
      !! Maximizes the orbital self-repulsion `D = sum_i (ii|ii)`. The calling
      !! shape and the meaning of every argument are `boys_localize`'s, except
      !! that `functional` returns D, in Hartree.
      !!
      !! Cost: the occupied integrals `(ij|kl)` come from `occupied_eri`, which
      !! is O(n_ao^4 n_occ) time and needs no AO tensor; it is most of the time. They are then held
      !! packed, `(n_occ(n_occ+1)/2)^2` doubles -- 51 MB at seventy orbitals --
      !! and every pair rotation updates them in O(n_occ^3), so a sweep is
      !! O(n_occ^5).
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coefficients(:, :)   !! MO coefficients, (n_ao, n_mo)
      integer, intent(in) :: n_occ                 !! Occupied orbitals to localize
      real(dp), allocatable, intent(out) :: localized(:, :)
         !! The localized occupied orbitals, (n_ao, n_occ).
      real(dp), allocatable, intent(out) :: centroids(:, :)
         !! `<i| r |i>` per localized orbital, (3, n_occ), Bohr.
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_sweeps
         !! Sweeps allowed over all threshold stages together, `ER_MAX_SWEEPS`
         !! by default. Zero returns the starting orbitals.
      real(dp), intent(in), optional :: angle_tol
      integer, intent(out), optional :: sweeps_taken
      real(dp), intent(out), optional :: functional
         !! `sum_i (ii|ii)` of the orbitals returned, Hartree.
      character(len=*), intent(in), optional :: guess
         !! Where the sweeps start: "canonical" (the default, and GAMESS's) or
         !! "boys", the Foster-Boys orbitals.
      logical, intent(out), optional :: converged
         !! Whether the last stage's sweep rotated nothing, rather than the
         !! sweep limit stopping it.

      real(dp), allocatable :: v(:, :), dip(:, :, :), work(:, :), d(:, :), start(:, :)
      real(dp), allocatable :: boys_centroids(:, :)
      real(dp) :: tol, thr, a, b, gamma, cg, sg, di, dj
      integer :: limit, sweep, i, j, k, m, n_pair
      logical :: rotated, done
      character(len=16) :: text

      if (n_occ < 1 .or. n_occ > size(coefficients, 2)) then
         write (text, "(i0)") n_occ
         call error%set(ERROR_VALIDATION, "cannot localize "//trim(text)// &
                        " orbitals: that is not a subset of the ones supplied")
         return
      end if

      limit = ER_MAX_SWEEPS
      if (present(max_sweeps)) limit = max_sweeps
      tol = DEFAULT_ANGLE_TOL
      if (present(angle_tol)) tol = angle_tol

      allocate (start(size(coefficients, 1), n_occ))
      start = coefficients(:, 1:n_occ)
      if (present(guess)) then
         select case (trim(guess))
         case ("canonical")
         case ("boys")
            call boys_localize(mol, coefficients, n_occ, localized, boys_centroids, error)
            if (error%has_error()) return
            start = localized
            deallocate (localized)
         case default
            call error%set(ERROR_VALIDATION, "er_localize: unknown guess '"// &
                           trim(guess)//"'; expected 'canonical' or 'boys'")
            return
         end select
      end if

      call occupied_eri(mol, start, v, error)
      if (error%has_error()) return

      allocate (localized(mol%nao, n_occ))
      localized = start
      n_pair = n_occ*(n_occ + 1)/2

      ! GAMESS's schedule (`LOCROT`, local.src): a ladder of thresholds from
      ! 0.1 down by tens, sweeping at each until no rotation is as large as the
      ! threshold, and skipping the smaller ones meanwhile. Pair order and
      ! rotation sense are its too -- (j,i) with j > i, j' = c j + s i -- so
      ! that from the same canonical orbitals both follow the same path to the
      ! same one of ER's several maxima.
      sweep = 0
      thr = 0.1_dp
      done = .false.
      stages: do
         do
            if (sweep >= limit) exit stages
            sweep = sweep + 1
            rotated = .false.
            do j = 2, n_occ
               do i = 1, j - 1
                  call er_pair_terms(v, j, i, a, b)
                  ! A pair that is flat in its angle has nothing to give, and
                  ! atan2(0,0) is not defined.
                  if (abs(a) < tiny(1.0_dp) .and. abs(b) < tiny(1.0_dp)) cycle
                  gamma = 0.25_dp*atan2(b, -a)
                  if (abs(gamma) < thr) cycle
                  rotated = .true.
                  cg = cos(gamma)
                  sg = sin(gamma)
                  do m = 1, mol%nao
                     dj = localized(m, j)
                     di = localized(m, i)
                     localized(m, j) = cg*dj + sg*di
                     localized(m, i) = -sg*dj + cg*di
                  end do
                  call rotate_pair_tensor(v, n_occ, j, i, cg, sg)
               end do
            end do
            if (.not. rotated) exit
         end do
         if (thr <= tol) then
            done = .true.
            exit stages
         end if
         thr = max(0.1_dp*thr, tol)
      end do stages

      if (present(sweeps_taken)) sweeps_taken = sweep
      if (present(converged)) converged = done
      if (present(functional)) then
         functional = 0.0_dp
         do i = 1, n_occ
            k = pair_index(i, i)
            functional = functional + v(k, k)
         end do
      end if

      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, error)
      if (error%has_error()) return
      allocate (centroids(3, n_occ), work(mol%nao, n_occ), d(n_occ, n_occ))
      do k = 1, 3
         call pic_gemm(dip(:, :, k), localized, work)
         call pic_gemm(localized, work, d, transa="T")
         do i = 1, n_occ
            centroids(k, i) = d(i, i)
         end do
      end do
   end subroutine er_localize

   pure subroutine er_pair_terms(v, p, q, a, b)
      !! The two coefficients of one pair's ER functional
      !!
      !! Rotating p' = c p + s q, q' = -s p + c q by gamma changes
      !! `(pp|pp) + (qq|qq)` by exactly
      !!
      !!     A (1 - cos 4g) + B sin 4g,
      !!     A = (pq|pq) - [(pp|pp) + (qq|qq) - 2 (pp|qq)]/4,
      !!     B = (pp|pq) - (qq|pq),
      !!
      !! whose maximum is at `g = atan2(B, -A)/4` and gains `A + sqrt(A^2+B^2)`,
      !! never negative. This is GAMESS's `LOCROT` with its `ALM` equal to -A.
      real(dp), intent(in) :: v(:, :)
      integer, intent(in) :: p, q
      real(dp), intent(out) :: a, b

      integer :: pp, qq, pq

      pp = pair_index(p, p)
      qq = pair_index(q, q)
      pq = pair_index(p, q)
      a = v(pq, pq) - 0.25_dp*(v(pp, pp) + v(qq, qq) - 2.0_dp*v(pp, qq))
      b = v(pp, pq) - v(qq, pq)
   end subroutine er_pair_terms

   subroutine rotate_pair_tensor(v, n, p, q, c, s)
      !! Carry `(ij|kl)`, packed over both pairs, through p' = c p + s q, q' = -s p + c q
      !!
      !! The rotation acts on the pair index as a linear map R, and the tensor
      !! goes to R V R^T: the rows first, then the columns. Only pairs that hold
      !! p or q move -- 2n of them -- so the cost is O(n^3).
      real(dp), intent(inout) :: v(:, :)
      integer, intent(in) :: n, p, q
      real(dp), intent(in) :: c, s

      integer :: col, row, k, kp, kq, pp, qq, pq
      real(dp) :: x, y, z

      pp = pair_index(p, p)
      qq = pair_index(q, q)
      pq = pair_index(p, q)

      ! Rows.
      do col = 1, size(v, 2)
         do k = 1, n
            if (k == p .or. k == q) cycle
            kp = pair_index(k, p)
            kq = pair_index(k, q)
            x = v(kp, col)
            y = v(kq, col)
            v(kp, col) = c*x + s*y
            v(kq, col) = -s*x + c*y
         end do
         x = v(pp, col)
         y = v(qq, col)
         z = v(pq, col)
         v(pp, col) = c*c*x + s*s*y + 2.0_dp*c*s*z
         v(qq, col) = s*s*x + c*c*y - 2.0_dp*c*s*z
         v(pq, col) = -c*s*x + c*s*y + (c*c - s*s)*z
      end do

      ! Columns.
      do k = 1, n
         if (k == p .or. k == q) cycle
         kp = pair_index(k, p)
         kq = pair_index(k, q)
         do row = 1, size(v, 1)
            x = v(row, kp)
            y = v(row, kq)
            v(row, kp) = c*x + s*y
            v(row, kq) = -s*x + c*y
         end do
      end do
      do row = 1, size(v, 1)
         x = v(row, pp)
         y = v(row, qq)
         z = v(row, pq)
         v(row, pp) = c*c*x + s*s*y + 2.0_dp*c*s*z
         v(row, qq) = s*s*x + c*c*y - 2.0_dp*c*s*z
         v(row, pq) = -c*s*x + c*s*y + (c*c - s*s)*z
      end do
   end subroutine rotate_pair_tensor

   subroutine occupied_eri(mol, orbitals, v, error)
      !! `(ij|kl)` over a set of orbitals, both pairs packed, without the AO tensor
      !!
      !! `v(pair_index(i,j), pair_index(k,l))`, square and symmetric. Two
      !! passes, both GEMM-shaped:
      !!
      !! 1. One task per bra shell pair (M >= N). Its integrals `(mn|pq)` over
      !!    every ket pair are gathered, a column of packed ket pairs per `mn`, and
      !!    each column is unpacked and taken to `(mn|kl)` by `C^T G C`. That
      !!    lands in `h(pair_index(m,n), kl)`; tasks own disjoint rows.
      !! 2. One task per `kl`: the column of `h` unpacked and taken to
      !!    `(ij|kl)` the same way.
      !!
      !! Time O(n_ao^4 n) in the GEMMs, the integrals under Schwarz screening.
      !! Memory `n_ao(n_ao+1)/2 * n(n+1)/2` doubles for `h`, capped at
      !! `ER_TRANSFORM_WORDS` by splitting the `kl` range into batches that
      !! each repeat pass 1, plus `d_M d_N n_ao(n_ao+1)/2` per thread.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)   !! (n_ao, n)
      real(dp), allocatable, intent(out) :: v(:, :)
      type(error_t), intent(inout) :: error

      type(eri_shell_table_t) :: tab
      type(c_ptr) :: opt
      real(dp), allocatable :: split_bounds(:, :), q(:, :), c(:, :), h(:, :)
      real(dp), allocatable :: buf(:), g(:, :), sq(:, :), t(:, :), w(:, :)
      integer, allocatable :: pair_m(:), pair_n(:)
      real(dp) :: bound
      integer :: nao, n, n_pair, ao_pair, n_task, task, msh, nsh, psh, qsh
      integer :: dm, dn, dp_, dq, mo, no, po, qo, ret
      integer :: shls(4)
      integer :: i, j, k, l, a1, a2, a3, a4, mn, row, kl, kl0, kl1, width, batch

      nao = mol%nao
      n = size(orbitals, 2)
      n_pair = n*(n + 1)/2
      ao_pair = nao*(nao + 1)/2
      allocate (c(nao, n))
      c = orbitals

      call eri_shell_table(mol, tab)
      call schwarz_bounds(mol, split_bounds, error)
      if (error%has_error()) return
      call eri_schwarz_collapse(mol, split_bounds, q)
      bound = maxval(q)

      opt = c_null_ptr
      call two_electron_optimizer(mol%cartesian, opt, mol%atm, mol%natm, tab%bas, &
                                  tab%nbas, tab%env)

      n_task = tab%nbas*(tab%nbas + 1)/2
      allocate (pair_m(n_task), pair_n(n_task))
      task = 0
      do msh = 1, tab%nbas
         do nsh = 1, msh
            task = task + 1
            pair_m(task) = msh
            pair_n(task) = nsh
         end do
      end do

      width = int(min(int(n_pair, int64), max(1_int64, ER_TRANSFORM_WORDS/int(ao_pair, int64))))
      allocate (v(n_pair, n_pair), h(ao_pair, width))
      v = 0.0_dp

      do batch = 1, (n_pair + width - 1)/width
         kl0 = (batch - 1)*width + 1
         kl1 = min(batch*width, n_pair)
         h = 0.0_dp

         !$omp parallel default(none) &
         !$omp    shared(mol, tab, opt, q, bound, c, h, pair_m, pair_n, n_task, nao, n, &
         !$omp           ao_pair, kl0, kl1) &
         !$omp    private(task, msh, nsh, psh, qsh, dm, dn, dp_, dq, mo, no, po, qo, ret, &
         !$omp            shls, i, j, k, l, a1, a2, a3, a4, mn, row, kl, buf, g, sq, t, w)
         allocate (buf(tab%block_max**4), g(ao_pair, tab%block_max**2))
         allocate (sq(nao, nao), t(nao, n), w(n, n))
         !$omp do schedule(dynamic)
         do task = 1, n_task
            msh = pair_m(task)
            nsh = pair_n(task)
            if (q(msh, nsh)*bound < ER_SCREEN_TOL) cycle
            dm = tab%dims(msh)
            dn = tab%dims(nsh)
            mo = tab%offs(msh)
            no = tab%offs(nsh)

            ! (mn|pq) for this bra block and every ket pair, a column of packed
            ! pq per mn.
            g(:, 1:dm*dn) = 0.0_dp
            do psh = 1, tab%nbas
               dp_ = tab%dims(psh)
               po = tab%offs(psh)
               do qsh = 1, psh
                  if (q(msh, nsh)*q(psh, qsh) < ER_SCREEN_TOL) cycle
                  dq = tab%dims(qsh)
                  qo = tab%offs(qsh)
                  shls = [msh - 1, nsh - 1, psh - 1, qsh - 1]
                  ret = two_electron_block(mol%cartesian, buf, shls, mol%atm, mol%natm, &
                                           tab%bas, tab%nbas, tab%env, opt)
                  if (ret == 0) cycle
                  do a4 = 1, dq
                     do a3 = 1, dp_
                        k = pair_index(po + a3, qo + a4)
                        do a2 = 1, dn
                           do a1 = 1, dm
                              g(k, a1 + (a2 - 1)*dm) = &
                                 buf(a1 + dm*((a2 - 1) + dn*((a3 - 1) + dp_*(a4 - 1))))
                           end do
                        end do
                     end do
                  end do
               end do
            end do

            ! Each (mn) column to (mn|kl), into h. A diagonal block holds both
            ! (mn) and (nm); the one with m >= n is enough.
            do a2 = 1, dn
               do a1 = 1, dm
                  if (msh == nsh .and. a1 < a2) cycle
                  mn = a1 + (a2 - 1)*dm
                  row = pair_index(mo + a1, no + a2)
                  call unpack_pairs(g(:, mn), nao, sq)
                  call pic_gemm(sq, c, t)
                  call pic_gemm(c, t, w, transa="T")
                  do l = 1, n
                     do k = l, n
                        kl = pair_index(k, l)
                        if (kl < kl0 .or. kl > kl1) cycle
                        h(row, kl - kl0 + 1) = w(k, l)
                     end do
                  end do
               end do
            end do
         end do
         !$omp end do
         deallocate (buf, g, sq, t, w)
         !$omp end parallel

         ! Each (kl) column of h to (ij|kl).
         !$omp parallel default(none) shared(h, c, v, nao, n, kl0, kl1) &
         !$omp    private(kl, i, j, sq, t, w)
         allocate (sq(nao, nao), t(nao, n), w(n, n))
         !$omp do schedule(dynamic)
         do kl = kl0, kl1
            call unpack_pairs(h(:, kl - kl0 + 1), nao, sq)
            call pic_gemm(sq, c, t)
            call pic_gemm(c, t, w, transa="T")
            do j = 1, n
               do i = j, n
                  v(pair_index(i, j), kl) = w(i, j)
               end do
            end do
         end do
         !$omp end do
         deallocate (sq, t, w)
         !$omp end parallel
      end do

      call libcint_del_optimizer(opt)

      ! The two triangles came by different routes; make them one.
      v = 0.5_dp*(v + transpose(v))
   end subroutine occupied_eri

   pure subroutine unpack_pairs(packed, nao, square)
      !! A symmetric matrix from its lower triangle, packed as `pair_index` orders it
      real(dp), intent(in) :: packed(:)
      integer, intent(in) :: nao
      real(dp), intent(out) :: square(:, :)   !! (nao, nao)

      integer :: p, r, pr

      pr = 0
      do p = 1, nao
         do r = 1, p - 1
            pr = pr + 1
            square(p, r) = packed(pr)
            square(r, p) = packed(pr)
         end do
         pr = pr + 1
         square(p, p) = packed(pr)
      end do
   end subroutine unpack_pairs

end module mqc_czt_localize
