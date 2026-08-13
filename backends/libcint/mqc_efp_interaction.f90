!! Interaction energies between effective fragments
module mqc_efp_interaction
   !! What a fragment potential is *for*: given several of them placed in space,
   !! the energy of their interaction.
   !!
   !! **Electrostatics first, and only electrostatics so far.** An EFP2 interaction
   !! is five terms -- electrostatics, polarization, exchange repulsion, dispersion
   !! and charge transfer -- and they differ enormously in what they need.
   !! Electrostatics and dispersion are geometry and stored tensors; polarization
   !! needs a self-consistent solve for the induced dipoles; exchange repulsion and
   !! charge transfer need integrals between *two fragments' basis sets*, which
   !! means building a combined molecule. This module carries the first.
   !!
   !! **The layout is flat on purpose.** Every fragment's expansion points are
   !! concatenated into one set of arrays with a fragment index per point, so the
   !! energy is one loop over point pairs that skips same-fragment pairs. That is
   !! the shape a parallel version wants: no nested dispatch over fragments to
   !! unpick, no per-fragment allocation inside the hot loop, and a reduction over
   !! a single scalar. Making it parallel is adding a directive, not a rewrite.
   !!
   !! **Every convention here was measured against GAMESS, not assumed.** The
   !! ladder is in `tools/efp_validation/zero_sections.py`: zero all the multipoles
   !! but one rank, ask GAMESS for the electrostatic energy of a dimer, and compare.
   !! With monopoles alone there is nothing to get wrong and it agreed to 4e-10,
   !! which fixed the geometry -- plain translation, bond midpoints included, the
   !! charge being electronic plus nuclear. Each higher rank was pinned the same
   !! way. Signs and component orders in a multipole expansion are exactly the
   !! things that look right and are wrong.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_efp_read, only: efp_fragment_t
   implicit none
   private

   public :: efp_system_t
   public :: build_efp_system
   public :: electrostatic_energy

   !> Components of each stored multipole, in the file's own order.
   integer, parameter :: N_DIPOLE = 3
   integer, parameter :: N_QUADRUPOLE = 6
   integer, parameter :: N_OCTUPOLE = 10

   !> Row and column of each stored quadrupole component. The file carries six
   !> numbers for a symmetric 3x3, ordered xx, yy, zz, xy, xz, yz.
   integer, parameter :: QUAD_I(N_QUADRUPOLE) = [1, 2, 3, 1, 1, 2]
   integer, parameter :: QUAD_J(N_QUADRUPOLE) = [1, 2, 3, 2, 3, 3]

   !> Indices of each stored octupole component. Ten numbers for a fully symmetric
   !> 3x3x3, ordered xxx, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz, xyz -- which is
   !> GAMESS's order, read off how `efelec.src` unpacks `EFOCT` rather than guessed.
   integer, parameter :: OCT_I(N_OCTUPOLE) = [1, 2, 3, 1, 1, 1, 2, 1, 2, 1]
   integer, parameter :: OCT_J(N_OCTUPOLE) = [1, 2, 3, 1, 1, 2, 2, 3, 3, 2]
   integer, parameter :: OCT_K(N_OCTUPOLE) = [1, 2, 3, 2, 3, 2, 3, 3, 3, 3]

   type :: efp_system_t
      !! Several placed fragments, flattened into one set of point arrays
      integer :: n_fragments = 0
      integer :: n_points = 0
      integer, allocatable :: fragment_of(:)     !! (n_points), which fragment
      real(dp), allocatable :: points(:, :)      !! (3, n_points), Bohr
      real(dp), allocatable :: charge(:)         !! Total monopole per point
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
      if (allocated(self%dipole)) deallocate (self%dipole)
      if (allocated(self%quad)) deallocate (self%quad)
      if (allocated(self%oct)) deallocate (self%oct)
      self%n_fragments = 0
      self%n_points = 0
   end subroutine system_destroy

   subroutine build_efp_system(fragments, translations, system, error)
      !! Flatten placed fragments into one point set
      !!
      !! **Translation only, for now.** A fragment placed by rotation as well as
      !! translation needs its multipoles rotated with it -- the dipole as a vector,
      !! the quadrupole and octupole as tensors of their rank -- and that is a
      !! separate piece of work with its own conventions to pin down. A translated
      !! copy needs none of it, and is enough to establish every term of the
      !! interaction energy, which is what comes first. Rotating a fragment whose
      !! multipoles were not rotated with it would be silently wrong, so it is not
      !! offered rather than offered and approximate.
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
      allocate (system%fragment_of(total), system%points(3, total), &
                system%charge(total), system%dipole(N_DIPOLE, total), &
                system%quad(3, 3, total), &
                system%oct(3, 3, 3, total))
      system%dipole = 0.0_dp
      system%quad = 0.0_dp
      system%oct = 0.0_dp

      at = 0
      do f = 1, n_frag
         do i = 1, fragments(f)%n_points
            at = at + 1
            system%fragment_of(at) = f
            system%points(:, at) = fragments(f)%points(:, i) + translations(:, f)
            ! The monopole a multipole expansion sees is the total at the point:
            ! GAMESS stores the electronic and nuclear parts separately because
            ! they come from different places, not because they act differently.
            system%charge(at) = fragments(f)%q_elec(i) + fragments(f)%q_nuc(i)
            if (allocated(fragments(f)%dipole)) then
               system%dipole(:, at) = fragments(f)%dipole(:, i)
            end if
            if (allocated(fragments(f)%quadrupole)) then
               ! The file packs a symmetric 3x3 as six numbers, xx yy zz xy xz yz.
               ! Unpacked here rather than in the energy so the inner loop sees
               ! plain matrix-vector products.
               do k = 1, N_QUADRUPOLE
                  system%quad(QUAD_I(k), QUAD_J(k), at) = fragments(f)%quadrupole(k, i)
                  system%quad(QUAD_J(k), QUAD_I(k), at) = fragments(f)%quadrupole(k, i)
               end do
            end if
            if (allocated(fragments(f)%octopole)) then
               ! Every permutation of each index triple carries the same value: the
               ! moment is fully symmetric and the file stores one of each set.
               do k = 1, N_OCTUPOLE
                  call spread_octupole(system%oct(:, :, :, at), OCT_I(k), OCT_J(k), &
                                       OCT_K(k), fragments(f)%octopole(k, i))
               end do
            end if
         end do
      end do
   end subroutine build_efp_system

   function electrostatic_energy(system, max_rank) result(energy)
      !! The multipole interaction energy between different fragments
      !!
      !! Undamped: this is the bare multipole sum, with no charge-penetration
      !! screening. The screening term is a separate correction with its own fit,
      !! and keeping it separate is what lets each rank be checked on its own.
      !!
      !! `max_rank` truncates the expansion -- 0 for charges only, 1 through the
      !! dipole, 2 the quadrupole, 3 the octupole. Present because it is what the
      !! validation ladder needs, and useful in its own right: the higher ranks fall
      !! off fast and a screening study wants to switch them off.
      type(efp_system_t), intent(in) :: system
      integer, intent(in) :: max_rank
      real(dp) :: energy

      real(dp) :: r(3)
      real(dp) :: pair
      integer :: a, b

      energy = 0.0_dp
      ! One loop over ordered pairs of points on different fragments. Each
      ! unordered pair is visited once, so nothing is double counted and no factor
      ! of a half appears anywhere.
      !
      ! Ready to thread: the iterations are independent and the only shared write
      ! is the reduction. Not threaded yet -- correctness first, and a directive
      ! here would be measuring nothing until the higher ranks are in.
      do a = 1, system%n_points - 1
         do b = a + 1, system%n_points
            if (system%fragment_of(a) == system%fragment_of(b)) cycle
            r = system%points(:, b) - system%points(:, a)
            pair = pair_energy(system, a, b, r, max_rank)
            energy = energy + pair
         end do
      end do
   end function electrostatic_energy

   pure function pair_energy(system, a, b, r, max_rank) result(e)
      !! One point pair, through the requested rank
      !!
      !! `r` points from `a` to `b`. Every term below is written in terms of that
      !! direction, so the antisymmetry between the two points -- a dipole on `a`
      !! sees the opposite field to one on `b` -- appears as the sign of `r` and
      !! nowhere else.
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
      ! A dipole in a potential has energy `mu . grad phi`, so for a charge at `a`
      ! and a dipole at `b` that is `-q_a (mu_b . r) / R^3`, and the other way round
      ! it is `+q_b (mu_a . r) / R^3` because the separation is `-r` seen from `b`.
      !
      ! Both signs were checked against GAMESS rather than trusted: with only
      ! monopoles and dipoles left in the potential the two possibilities differ by
      ! 4e-03 Hartree on a water dimer, which is not a subtle discrepancy, and the
      ! first version of this line had it backwards. GAMESS's stored dipoles are in
      ! the ordinary convention -- there is nothing to undo on the way in.
      e = e - qa*rb_dot_db*inv3 + qb*ra_dot_da*inv3

      ! --- dipole-dipole -------------------------------------------------------
      inv5 = inv3*inv*inv
      e = e + dot_product(da, db)*inv3 - 3.0_dp*ra_dot_da*rb_dot_db*inv5
      if (max_rank < 2) return

      ! --- rank two ------------------------------------------------------------
      !
      ! Written as closed-form contractions of the derivative tensors rather than
      ! as loops over their components: the quadrupole-quadrupole term is a
      ! rank-four tensor with 81 entries, and contracting it literally costs 81
      ! multiply-adds per pair where the closed form costs two matrix-vector
      ! products. Correct either way, and this one stays cheap when the pair count
      ! grows.
      !
      ! The three factors -- a half, minus a half, a quarter -- were each solved
      ! for against GAMESS on a potential with everything but the relevant ranks
      ! zeroed, and each came out an exact rational to seven figures. GAMESS's own
      ! source groups them differently, carrying a ninth on a differently
      ! normalized quadrupole expression; the physics is the same and this grouping
      ! is the one that follows from the derivative tensors.
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
      ! GAMESS's electrostatic total is ECC + ECD + EDD + ECQ + ECO + EDQ + EQQ.
      ! There is no dipole-octupole, quadrupole-octupole or octupole-octupole term
      ! in it -- that is the truncation, read off `FFELEC` in `efelec.src` rather
      ! than inferred from the energy, which could not have distinguished a missing
      ! term from a small one.
      !
      ! The trace of the stored octupole cannot contribute, because `T3` is itself
      ! traceless in any index pair, which is why the raw file values serve with a
      ! single factor and no traceless projection. GAMESS projects and carries a
      ! fifteenth; the same physics arrives here as a sixth on the unprojected
      ! moment.
      oa_rrr = triple_contract(system%oct(:, :, :, a), r)
      ob_rrr = triple_contract(system%oct(:, :, :, b), r)
      oa_tr = octupole_trace(system%oct(:, :, :, a))
      ob_tr = octupole_trace(system%oct(:, :, :, b))
      e = e + (1.0_dp/6.0_dp) &
          *(-qa*(15.0_dp*ob_rrr - 9.0_dp*r2*dot_product(r, ob_tr)) &
            + qb*(15.0_dp*oa_rrr - 9.0_dp*r2*dot_product(r, oa_tr)))*inv7
   end function pair_energy

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

end module mqc_efp_interaction
