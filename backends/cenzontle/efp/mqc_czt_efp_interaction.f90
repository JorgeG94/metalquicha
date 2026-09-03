!! Interaction energies between effective fragments
module mqc_czt_efp_interaction
   !! What a fragment potential is *for*: given several of them placed in space,
   !! the energy of their interaction.
   !!
   !! Of the five EFP2 terms this carries electrostatics, polarization and the
   !! undamped dispersion `E6`. Exchange repulsion and charge transfer need
   !! integrals between *two fragments' basis sets*, which means building a
   !! combined molecule, and are absent.
   !!
   !! Every fragment's expansion points are concatenated into one set of arrays
   !! with a fragment index per point, so an energy is one loop over point pairs
   !! that skips same-fragment pairs.
   !!
   !! The signs and component orders were pinned against GAMESS rank by rank
   !! with `tools/efp_validation/zero_sections.py`, which zeroes all the
   !! multipoles but one and compares a dimer energy.
   use pic_types, only: dp
   use mqc_physical_constants, only: PI
   use pic_logger, only: logger => global_logger, verbose_level
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_efp_read, only: efp_fragment_t
   implicit none
   private

   public :: efp_system_t
   public :: build_efp_system
   public :: electrostatic_energy
   public :: dispersion_energy_e6
   public :: polarization_energy
   public :: CP_WEIGHT

   integer, parameter :: N_DIPOLE = 3
   integer, parameter :: N_QUADRUPOLE = 6
   integer, parameter :: N_OCTUPOLE = 10
      !! Components of each stored multipole, in the file's own order.

   integer, parameter :: QUAD_I(N_QUADRUPOLE) = [1, 2, 3, 1, 1, 2]
   integer, parameter :: QUAD_J(N_QUADRUPOLE) = [1, 2, 3, 2, 3, 3]
      !! Row and column of each stored quadrupole component. The file carries
      !! six numbers for a symmetric 3x3, ordered xx, yy, zz, xy, xz, yz.

   integer, parameter :: OCT_I(N_OCTUPOLE) = [1, 2, 3, 1, 1, 1, 2, 1, 2, 1]
   integer, parameter :: OCT_J(N_OCTUPOLE) = [1, 2, 3, 1, 1, 2, 2, 3, 3, 2]
   integer, parameter :: OCT_K(N_OCTUPOLE) = [1, 2, 3, 2, 3, 2, 3, 3, 3, 3]
      !! Indices of each stored octupole component. Ten numbers for a fully
      !! symmetric 3x3x3, ordered xxx, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz,
      !! xyz, which is GAMESS's order in `efelec.src`.

   integer, parameter :: N_FREQUENCIES = 12
   real(dp), parameter :: CP_WEIGHT(N_FREQUENCIES) = [ &
                          0.72086099022968040154e-02_dp, 0.17697067815034886394e-01_dp, &
                          0.30660908596251749739e-01_dp, 0.48381293256249884995e-01_dp, &
                          0.74878830420650517080e-01_dp, 0.11806515901361630228e+00_dp, &
                          0.19535413832209084204e+00_dp, 0.35055692324483221824e+00_dp, &
                          0.71577113554429568336e+00_dp, 1.81409759976323969729e+00_dp, &
                          6.97923445114870823247e+00_dp, 83.2480938829658453917e+00_dp]
      !! The twelve Casimir-Polder quadrature weights, as GAMESS carries them in
      !! `efdrvr.src`: 12-point Gauss-Legendre weights times the Jacobian of
      !! `nu = w0 (1 + t)/(1 - t)` at `w0 = 0.3`.

   type :: efp_system_t
      !! Several placed fragments, flattened into one set of point arrays
      integer :: n_fragments = 0
      integer :: n_points = 0
      integer, allocatable :: fragment_of(:)     !! (n_points), which fragment
      real(dp), allocatable :: points(:, :)      !! (3, n_points), Bohr
      real(dp), allocatable :: charge(:)         !! Total monopole per point
      real(dp), allocatable :: q_elec(:)         !! Electronic part, for screening
      real(dp), allocatable :: q_nuc(:)          !! Nuclear part, for screening
      real(dp), allocatable :: alpha(:)          !! Exponential damping exponent
      logical :: has_screening = .false.
      real(dp), allocatable :: offset(:, :)      !! (3, n_fragments), how each was placed
      real(dp), allocatable :: dipole(:, :)      !! (3, n_points)
      real(dp), allocatable :: quad(:, :, :)     !! (3, 3, n_points), symmetric
      real(dp), allocatable :: oct(:, :, :, :)   !! (3, 3, 3, n_points), symmetric
   contains
      procedure :: destroy => system_destroy
   end type efp_system_t

contains

   subroutine system_destroy(self)
      class(efp_system_t), intent(inout) :: self

      if (allocated(self%fragment_of)) deallocate (self%fragment_of)
      if (allocated(self%points)) deallocate (self%points)
      if (allocated(self%charge)) deallocate (self%charge)
      if (allocated(self%q_elec)) deallocate (self%q_elec)
      if (allocated(self%q_nuc)) deallocate (self%q_nuc)
      if (allocated(self%alpha)) deallocate (self%alpha)
      if (allocated(self%offset)) deallocate (self%offset)
      self%has_screening = .false.
      if (allocated(self%dipole)) deallocate (self%dipole)
      if (allocated(self%quad)) deallocate (self%quad)
      if (allocated(self%oct)) deallocate (self%oct)
      self%n_fragments = 0
      self%n_points = 0
   end subroutine system_destroy

   subroutine build_efp_system(fragments, translations, system, error)
      !! Flatten placed fragments into one point set
      !!
      !! **Translation only.** A fragment placed by rotation needs its multipoles
      !! rotated with it -- the dipole as a vector, the quadrupole and octupole
      !! as tensors of their rank -- which is not done here, so no rotation is
      !! accepted rather than one silently ignored.
      type(efp_fragment_t), intent(in) :: fragments(:)
      real(dp), intent(in) :: translations(:, :)   !! (3, n_fragments), Bohr
      type(efp_system_t), intent(out) :: system
      type(error_t), intent(inout) :: error

      integer :: n_frag, total, f, i, at, k

      n_frag = size(fragments)
      if (size(translations, 1) /= 3 .or. size(translations, 2) /= n_frag) then
         call error%set(ERROR_VALIDATION, "efp: translations must be (3, n_fragments)")
         return
      end if

      total = 0
      do f = 1, n_frag
         total = total + fragments(f)%n_points
      end do
      if (total == 0) then
         call error%set(ERROR_VALIDATION, "efp: no expansion points in any fragment")
         return
      end if

      system%n_fragments = n_frag
      system%n_points = total
      allocate (system%offset(3, n_frag))
      system%offset = translations
      allocate (system%fragment_of(total), system%points(3, total), &
                system%charge(total), system%q_elec(total), &
                system%q_nuc(total), system%alpha(total), &
                system%dipole(N_DIPOLE, total), &
                system%quad(3, 3, total), &
                system%oct(3, 3, 3, total))
      system%dipole = 0.0_dp
      system%alpha = 0.0_dp
      system%has_screening = .true.
      system%quad = 0.0_dp
      system%oct = 0.0_dp

      at = 0
      do f = 1, n_frag
         do i = 1, fragments(f)%n_points
            at = at + 1
            system%fragment_of(at) = f
            system%points(:, at) = fragments(f)%points(:, i) + translations(:, f)
            ! The monopole a multipole expansion sees is the total at the point.
            system%charge(at) = fragments(f)%q_elec(i) + fragments(f)%q_nuc(i)
            ! Kept apart as well as summed: the penetration correction screens
            ! the electronic cloud against the other point's electrons and
            ! nucleus by different amounts, so their sum is not enough there.
            system%q_elec(at) = fragments(f)%q_elec(i)
            system%q_nuc(at) = fragments(f)%q_nuc(i)
            if (fragments(f)%has_screen2) then
               system%alpha(at) = fragments(f)%screen2(i)
            else
               system%has_screening = .false.
            end if
            if (allocated(fragments(f)%dipole)) then
               system%dipole(:, at) = fragments(f)%dipole(:, i)
            end if
            if (allocated(fragments(f)%quadrupole)) then
               ! The file packs a symmetric 3x3 as six numbers, xx yy zz xy xz
               ! yz. Unpacked here so the inner loop sees plain matrix-vector
               ! products.
               do k = 1, N_QUADRUPOLE
                  system%quad(QUAD_I(k), QUAD_J(k), at) = fragments(f)%quadrupole(k, i)
                  system%quad(QUAD_J(k), QUAD_I(k), at) = fragments(f)%quadrupole(k, i)
               end do
            end if
            if (allocated(fragments(f)%octopole)) then
               ! Every permutation of each index triple carries the same
               ! value: the moment is fully symmetric and the file stores one of
               ! each set.
               do k = 1, N_OCTUPOLE
                  call spread_octupole(system%oct(:, :, :, at), OCT_I(k), OCT_J(k), &
                                       OCT_K(k), fragments(f)%octopole(k, i))
               end do
            end if
         end do
      end do
   end subroutine build_efp_system

   function electrostatic_energy(system, max_rank, screen) result(energy)
      !! The multipole interaction energy between different fragments
      !!
      !! `max_rank` truncates the expansion: 0 for charges only, 1 through the
      !! dipole, 2 the quadrupole, 3 the octupole.
      type(efp_system_t), intent(in) :: system
      integer, intent(in) :: max_rank
      logical, intent(in), optional :: screen
         !! Add the charge-penetration correction. Off by default, so the bare
         !! multipole sum -- what each rank was validated against -- stays
         !! reachable.
      real(dp) :: energy

      real(dp) :: r(3)
      real(dp) :: pair
      integer :: a, b
      logical :: screening

      screening = .false.
      if (present(screen)) screening = screen

      energy = 0.0_dp
      ! One loop over ordered pairs of points on different fragments. Each
      ! unordered pair is visited once, so nothing is double counted and no
      ! factor of a half appears anywhere.
      !
      ! Threaded over the outer index with a reduction on the one scalar. The
      ! pair count is triangular, so a static schedule would leave the threads
      ! that drew the high indices idle -- hence `guided`.
      !$omp parallel do default(none) private(a, b, r, pair) &
      !$omp shared(system, max_rank, screening) schedule(guided) &
      !$omp reduction(+:energy) if(system%n_points >= 64)
      do a = 1, system%n_points - 1
         do b = a + 1, system%n_points
            if (system%fragment_of(a) == system%fragment_of(b)) cycle
            r = system%points(:, b) - system%points(:, a)
            pair = pair_energy(system, a, b, r, max_rank)
            if (screening .and. system%has_screening) then
               pair = pair + penetration(system, a, b, r)
            end if
            energy = energy + pair
         end do
      end do
      !$omp end parallel do
   end function electrostatic_energy

   pure function pair_energy(system, a, b, r, max_rank) result(e)
      !! One point pair, through the requested rank
      !!
      !! `r` points from `a` to `b`, and every term is written in terms of that
      !! direction, so the antisymmetry between the two points -- a dipole on
      !! `a` sees the opposite field to one on `b` -- appears as the sign of `r`
      !! and nowhere else.
      type(efp_system_t), intent(in) :: system
      integer, intent(in) :: a, b
      real(dp), intent(in) :: r(3)
      integer, intent(in) :: max_rank
      real(dp) :: e

      real(dp) :: qa_r(3), qb_r(3)
      real(dp) :: dist, inv, inv3, inv5, inv7, inv9, r2
      real(dp) :: r_qa_r, r_qb_r, tr_qa, tr_qb
      real(dp) :: oa_tr(3), ob_tr(3)
      real(dp) :: oa_rrr, ob_rrr
      real(dp) :: da(3), db(3)
      real(dp) :: qa, qb
      real(dp) :: ra_dot_da, rb_dot_db

      dist = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
      inv = 1.0_dp/dist
      e = 0.0_dp

      qa = system%charge(a)
      qb = system%charge(b)

      ! --- charge-charge, rank 0 ------------------------------------------------
      e = e + qa*qb*inv
      if (max_rank < 1) return

      inv3 = inv*inv*inv
      da = system%dipole(:, a)
      db = system%dipole(:, b)
      ra_dot_da = dot_product(r, da)
      rb_dot_db = dot_product(r, db)

      ! --- charge-dipole -------------------------------------------------------
      !
      ! A dipole in a potential has energy `mu . grad phi`, so for a charge at
      ! `a` and a dipole at `b` that is `-q_a (mu_b . r) / R^3`, and the other
      ! way round it is `+q_b (mu_a . r) / R^3` because the separation is `-r`
      ! seen from `b`. GAMESS's stored dipoles are in the ordinary convention,
      ! so there is nothing to undo on the way in.
      e = e - qa*rb_dot_db*inv3 + qb*ra_dot_da*inv3

      ! --- dipole-dipole -------------------------------------------------------
      inv5 = inv3*inv*inv
      e = e + dot_product(da, db)*inv3 - 3.0_dp*ra_dot_da*rb_dot_db*inv5
      if (max_rank < 2) return

      ! --- rank two ------------------------------------------------------------
      !
      ! Written as closed-form contractions of the derivative tensors rather
      ! than as loops over their components. GAMESS's own source groups the
      ! prefactors differently, carrying a ninth on a differently normalized
      ! quadrupole expression; the physics is the same.
      qa_r = matmul(system%quad(:, :, a), r)
      qb_r = matmul(system%quad(:, :, b), r)
      r_qa_r = dot_product(r, qa_r)
      r_qb_r = dot_product(r, qb_r)
      tr_qa = system%quad(1, 1, a) + system%quad(2, 2, a) + system%quad(3, 3, a)
      tr_qb = system%quad(1, 1, b) + system%quad(2, 2, b) + system%quad(3, 3, b)
      r2 = dist*dist
      inv7 = inv5*inv*inv
      inv9 = inv7*inv*inv

      ! charge-quadrupole, `Q : T2` against each charge
      e = e + 0.5_dp*(qa*(3.0_dp*r_qb_r - r2*tr_qb) &
                      + qb*(3.0_dp*r_qa_r - r2*tr_qa))*inv5

      ! dipole-quadrupole, `mu . (Q : T3)`, antisymmetric between the two points
      ! because the tensor is of odd rank
      e = e + 0.5_dp*((15.0_dp*ra_dot_da*r_qb_r &
                       - 3.0_dp*r2*(ra_dot_da*tr_qb + 2.0_dp*dot_product(da, qb_r))) &
                      - (15.0_dp*rb_dot_db*r_qa_r &
                         - 3.0_dp*r2*(rb_dot_db*tr_qa + 2.0_dp*dot_product(db, qa_r)))) &
          *inv7

      ! quadrupole-quadrupole, `Q : T4 : Q`
      e = e + 0.25_dp*(105.0_dp*r_qa_r*r_qb_r &
                       - 15.0_dp*r2*(r_qa_r*tr_qb + tr_qa*r_qb_r &
                                     + 4.0_dp*dot_product(qa_r, qb_r)) &
                       + 3.0_dp*r2*r2*(tr_qa*tr_qb &
                                       + 2.0_dp*sum(system%quad(:, :, a) &
                                                    *system%quad(:, :, b))))*inv9
      if (max_rank < 3) return

      ! --- rank three: charge-octupole, and nothing else ------------------------
      !
      ! GAMESS's electrostatic total is ECC + ECD + EDD + ECQ + ECO + EDQ + EQQ,
      ! with no dipole-octupole, quadrupole-octupole or octupole-octupole term.
      ! That truncation is `FFELEC` in `efelec.src`, and it is matched here.
      !
      ! The raw file values serve with a single factor and no traceless
      ! projection: `T3` is itself traceless in any index pair. GAMESS projects
      ! and carries a fifteenth, which arrives here as a sixth on the
      ! unprojected moment.
      oa_rrr = triple_contract(system%oct(:, :, :, a), r)
      ob_rrr = triple_contract(system%oct(:, :, :, b), r)
      oa_tr = octupole_trace(system%oct(:, :, :, a))
      ob_tr = octupole_trace(system%oct(:, :, :, b))
      e = e + (1.0_dp/6.0_dp) &
          *(-qa*(15.0_dp*ob_rrr - 9.0_dp*r2*dot_product(r, ob_tr)) &
            + qb*(15.0_dp*oa_rrr - 9.0_dp*r2*dot_product(r, oa_tr)))*inv7
   end function pair_energy

   function polarization_energy(system, fragments, error, max_iter, tol) result(energy)
      !! Polarization: induced dipoles at the orbital centroids, solved together
      !!
      !! Each polarizable point carries a static polarizability tensor and sits in
      !! the field of every *other* fragment's multipoles. Its induced dipole is
      !! `mu = alpha F`, and because each induced dipole adds to the field the others
      !! see, the set is solved self-consistently. The energy is
      !!
      !!     E = -(1/2) sum_k mu_k . F_k
      !!
      !! against the **static** field -- from the permanent multipoles only, not
      !! including the induced dipoles' own contribution, which is what makes the
      !! half correct rather than double counting.
      !!
      !! **The field runs to the quadrupole and stops**, as GAMESS's does:
      !! charges, dipoles and quadrupoles contribute, octupoles do not.
      !!
      !! The polarizability tensor is not symmetric -- a localized-orbital
      !! polarizability has an antisymmetric part -- so `alpha F` is a genuine
      !! matrix-vector product and the order matters.
      type(efp_system_t), intent(in) :: system
      type(efp_fragment_t), intent(in) :: fragments(:)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      real(dp) :: energy

      integer, parameter :: DEFAULT_ITER = 200
      real(dp), parameter :: DEFAULT_TOL = 1.0e-12_dp
      real(dp), allocatable :: centres(:, :), pol(:, :, :), fstat(:, :)
      real(dp), allocatable :: mu(:, :), mu_new(:, :)
      real(dp) :: field(3)
      integer, allocatable :: owner(:)
      real(dp) :: r(3)
      real(dp) :: dist, inv3, inv5, change
      integer :: n_pol, f, i, k, at, iter, limit, j
      real(dp) :: use_tol
      integer :: log_level
      logical :: show_iter

      energy = 0.0_dp
      limit = DEFAULT_ITER
      if (present(max_iter)) limit = max_iter
      use_tol = DEFAULT_TOL
      if (present(tol)) use_tol = tol

      ! The caller (`run_efp`) is already rank-0 only, so no leader guard is
      ! needed around the per-iteration logging.
      call logger%configuration(level=log_level)
      show_iter = log_level >= verbose_level

      n_pol = 0
      do f = 1, size(fragments)
         if (fragments(f)%has_static_pol) n_pol = n_pol + fragments(f)%n_pol
      end do
      if (n_pol == 0) return

      ! Flattened exactly as the multipole points are, so the induction loop is
      ! over pairs of polarizable points with a fragment index to skip on.
      allocate (centres(3, n_pol), pol(3, 3, n_pol), owner(n_pol))
      allocate (fstat(3, n_pol), mu(3, n_pol), mu_new(3, n_pol))
      at = 0
      do f = 1, size(fragments)
         if (.not. fragments(f)%has_static_pol) cycle
         do i = 1, fragments(f)%n_pol
            at = at + 1
            centres(:, at) = fragments(f)%pol_points(:, i) + system%offset(:, f)
            pol(:, :, at) = fragments(f)%static_pol(:, :, i)
            owner(at) = f
         end do
      end do

      ! The static field, from the permanent multipoles of the other fragments.
      fstat = 0.0_dp
      ! Each polarizable point writes only its own column, so no reduction is
      ! needed here.
      !$omp parallel do default(none) private(i, k, r, dist, inv3, inv5) &
      !$omp shared(n_pol, system, centres, owner, fstat) schedule(static) &
      !$omp if(n_pol >= 32)
      do i = 1, n_pol
         do k = 1, system%n_points
            if (system%fragment_of(k) == owner(i)) cycle
            r = centres(:, i) - system%points(:, k)
            dist = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
            inv3 = 1.0_dp/(dist*dist*dist)
            inv5 = inv3/(dist*dist)
            ! Charge, then dipole, then quadrupole -- the quadrupole carrying
            ! the same half the charge-quadrupole energy does.
            fstat(:, i) = fstat(:, i) + system%charge(k)*r*inv3 &
                          + 3.0_dp*r*dot_product(system%dipole(:, k), r)*inv5 &
                          - system%dipole(:, k)*inv3 &
                          - 0.5_dp*quadrupole_field(system%quad(:, :, k), r)
         end do
      end do
      !$omp end parallel do

      mu = 0.0_dp
      change = huge(1.0_dp)
      do iter = 1, limit
         ! A Jacobi sweep: every new dipole is built from the previous
         ! iteration's set, so the points are independent within a sweep.
         change = 0.0_dp
         !$omp parallel do default(none) private(i, j, r, dist, inv3, inv5, field) &
         !$omp shared(n_pol, owner, centres, mu, mu_new, fstat, pol) &
         !$omp schedule(static) reduction(max:change) if(n_pol >= 32)
         do i = 1, n_pol
            field = fstat(:, i)
            do j = 1, n_pol
               if (owner(j) == owner(i)) cycle
               r = centres(:, i) - centres(:, j)
               dist = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
               inv3 = 1.0_dp/(dist*dist*dist)
               inv5 = inv3/(dist*dist)
               field = field + 3.0_dp*r*dot_product(mu(:, j), r)*inv5 &
                       - mu(:, j)*inv3
            end do
            mu_new(:, i) = matmul(pol(:, :, i), field)
            ! Convergence measured inside the same loop, as a max reduction,
            ! rather than by a serial pass over the whole array afterwards.
            change = max(change, maxval(abs(mu_new(:, i) - mu(:, i))))
         end do
         !$omp end parallel do
         mu = mu_new
         if (show_iter) then
            call logger%verbose("  efp induction "//to_char(iter)//": max change "// &
                                to_char(change))
         end if
         if (change < use_tol) exit
      end do
      if (change >= use_tol) then
         call error%set(ERROR_VALIDATION, "efp: the induced dipoles did not converge")
         return
      end if

      do i = 1, n_pol
         energy = energy - 0.5_dp*dot_product(mu(:, i), fstat(:, i))
      end do

      deallocate (centres, pol, owner, fstat, mu, mu_new)
   end function polarization_energy

   pure function quadrupole_field(quad, r) result(f)
      !! `Q_jk T3_ijk`, the shape the quadrupole's field takes
      real(dp), intent(in) :: quad(3, 3), r(3)
      real(dp) :: f(3)

      real(dp) :: qr(3)
      real(dp) :: dist, r2, r_q_r, tr_q, inv7

      dist = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
      r2 = dist*dist
      inv7 = 1.0_dp/(r2*r2*r2*dist)
      qr = matmul(quad, r)
      r_q_r = dot_product(r, qr)
      tr_q = quad(1, 1) + quad(2, 2) + quad(3, 3)
      f = -(15.0_dp*r*r_q_r - 3.0_dp*r2*(r*tr_q + 2.0_dp*qr))*inv7
   end function quadrupole_field

   function dispersion_energy_e6(system, fragments) result(energy)
      !! The `E6` dispersion energy, undamped
      !!
      !!     E6 = - sum_ij C6_ij / R_ij^6,   C6_ij = (3/pi) sum_k w_k a_i(k) a_j(k)
      !!
      !! summed over localized orbitals on different fragments, with `a` the
      !! isotropic dynamic polarizability -- a third of the tensor trace -- at
      !! each of the twelve Casimir-Polder frequencies. The form, the factor of
      !! `3/pi` and the weights are GAMESS's, from `efdrvr.src`.
      !!
      !! **Undamped, where GAMESS reports a damped number.** Its damping is
      !! built from overlaps between the two fragments' localized orbitals,
      !! which needs integrals over both fragments' basis sets and is not
      !! present here, so this runs about two parts in a thousand larger in
      !! magnitude than GAMESS's E6.
      !!
      !! Takes the fragments rather than the flattened system: dispersion sits
      !! on the localized orbital centroids, a different set of points from the
      !! multipole expansion.
      ! TODO(mqc): the quadrature index runs to the smaller of the two
      ! fragments' `n_freq`, which `mqc_czt_efp_read` takes from the file without
      ! checking it against `N_FREQUENCIES`. A potential written with more than
      ! twelve frequency blocks reads past the end of `CP_WEIGHT`.
      type(efp_system_t), intent(in) :: system
      type(efp_fragment_t), intent(in) :: fragments(:)
      real(dp) :: energy

      real(dp) :: sep(3)
      real(dp) :: c6, r6, dist
      integer :: fa, fb, ia, ib, k

      energy = 0.0_dp
      ! Threaded over the first fragment with a reduction, `guided` because the
      ! pair set is triangular.
      !$omp parallel do default(none) private(fa, fb, ia, ib, k, c6, sep, dist, r6) &
      !$omp shared(fragments, system) schedule(guided) reduction(+:energy) &
      !$omp if(size(fragments) >= 8)
      do fa = 1, size(fragments) - 1
         if (.not. fragments(fa)%has_dynamic) cycle
         do fb = fa + 1, size(fragments)
            if (.not. fragments(fb)%has_dynamic) cycle
            do ia = 1, fragments(fa)%n_lmo
               do ib = 1, fragments(fb)%n_lmo
                  c6 = 0.0_dp
                  do k = 1, min(fragments(fa)%n_freq, fragments(fb)%n_freq)
                     c6 = c6 + CP_WEIGHT(k) &
                          *isotropic(fragments(fa)%dyn_pol(:, :, ia, k)) &
                          *isotropic(fragments(fb)%dyn_pol(:, :, ib, k))
                  end do
                  c6 = c6*3.0_dp/PI
                  sep = fragments(fb)%centroids(:, ib) + system%offset(:, fb) &
                        - fragments(fa)%centroids(:, ia) - system%offset(:, fa)
                  dist = sqrt(sep(1)*sep(1) + sep(2)*sep(2) + sep(3)*sep(3))
                  r6 = dist**6
                  energy = energy - c6/r6
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end function dispersion_energy_e6

   pure function isotropic(tensor) result(a)
      !! A third of the trace, which is what the isotropic `C6` is built from
      real(dp), intent(in) :: tensor(3, 3)
      real(dp) :: a

      a = (tensor(1, 1) + tensor(2, 2) + tensor(3, 3))/3.0_dp
   end function isotropic

   pure function penetration(system, a, b, r) result(e)
      !! Charge-charge penetration between two points
      !!
      !! A multipole expansion overstates the repulsion where the two charge
      !! densities actually overlap, since real electron clouds penetrate each
      !! other and screen the nuclei they surround. The correction is
      !! exponential in `alpha R`, so it matters at hydrogen-bond distances and
      !! vanishes by a few Angstrom.
      !!
      !! Only the charge-charge term, from `EPENCHCH` in GAMESS's `efelec.src`.
      !! GAMESS's further corrections for the dipole, quadrupole and octupole
      !! sit behind its `HOCHPEN` flag, which a default run does not set.
      !!
      !! The screening *coefficient* is not read: every potential this program
      !! writes carries a coefficient of one at every point, where GAMESS treats
      !! a zero coefficient as unscreened by setting the exponent to 1e20.
      type(efp_system_t), intent(in) :: system
      integer, intent(in) :: a, b
      real(dp), intent(in) :: r(3)
      real(dp) :: e

      real(dp), parameter :: SAME = 1.0e-5_dp
      real(dp) :: dist, aa, ab, ap, bp, aa2, ab2
      real(dp) :: p0_e, p0_n1, p0_n2

      dist = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
      aa = system%alpha(a)
      ab = system%alpha(b)
      ap = aa*dist
      bp = ab*dist
      aa2 = aa*aa
      ab2 = ab*ab

      if (abs(aa - ab) < SAME) then
         ! The two-exponent form is singular when the exponents coincide; this is
         ! its limit, not a separate approximation.
         p0_e = -exp(-ap)*(1.0_dp + 0.5_dp*ap)
         p0_n1 = -exp(-ap)
         p0_n2 = -exp(-ap)
      else
         p0_e = -exp(-ap)*(ab2/(ab2 - aa2)) - exp(-bp)*(aa2/(aa2 - ab2))
         p0_n1 = -exp(-bp)
         p0_n2 = -exp(-ap)
      end if

      ! `a`'s electrons screened against `b`'s nucleus carry `a`'s own exponent,
      ! and the other way round carries `b`'s.
      e = (system%q_elec(a)*system%q_elec(b)*p0_e &
           + system%q_elec(a)*system%q_nuc(b)*p0_n2 &
           + system%q_elec(b)*system%q_nuc(a)*p0_n1)/dist
   end function penetration

   pure subroutine spread_octupole(oct, i, j, k, value)
      !! One stored component written into every permutation of its indices
      real(dp), intent(inout) :: oct(3, 3, 3)
      integer, intent(in) :: i, j, k
      real(dp), intent(in) :: value

      oct(i, j, k) = value
      oct(i, k, j) = value
      oct(j, i, k) = value
      oct(j, k, i) = value
      oct(k, i, j) = value
      oct(k, j, i) = value
   end subroutine spread_octupole

   pure function triple_contract(oct, r) result(v)
      !! `Omega_ijk r_i r_j r_k`
      real(dp), intent(in) :: oct(3, 3, 3), r(3)
      real(dp) :: v

      integer :: i, j, k

      v = 0.0_dp
      do k = 1, 3
         do j = 1, 3
            do i = 1, 3
               v = v + oct(i, j, k)*r(i)*r(j)*r(k)
            end do
         end do
      end do
   end function triple_contract

   pure function octupole_trace(oct) result(t)
      !! `Omega_ijj`, the vector left by tracing the last two indices
      real(dp), intent(in) :: oct(3, 3, 3)
      real(dp) :: t(3)

      integer :: i, j

      t = 0.0_dp
      do i = 1, 3
         do j = 1, 3
            t(i) = t(i) + oct(i, j, j)
         end do
      end do
   end function octupole_trace

end module mqc_czt_efp_interaction
