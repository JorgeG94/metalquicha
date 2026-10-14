!! A molecule spanning two fragments, and the overlaps between them
module mqc_czt_efp_pair
   !! Exchange repulsion, charge transfer and the damping on dispersion all come from
   !! the same thing: overlaps between *one fragment's* localized orbitals and
   !! *another's*. Those are integrals over both fragments' basis functions at once,
   !! so they need a molecule that spans the pair -- which is what this builds.
   !!
   !! The basis comes from the potential's own `PROJECTION BASIS SET`, recovered by
   !! `mqc_czt_efp_read` with GAMESS's primitive normalization divided back out. The
   !! *file's* basis, not a lookup of the name it was computed with: a potential is
   !! a self-contained object, and the shipped GAMESS library potentials do not all
   !! name a basis this program has.
   use pic_types, only: dp
   use mqc_physical_constants, only: PI
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cgto, only: molecular_basis_type, ANGULAR_FORM_CARTESIAN
   use mqc_elements, only: element_number_to_symbol
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_efp_read, only: efp_fragment_t, N_QUADRUPOLE
   use mqc_czt_cphf, only: casimir_polder_weights, N_CASIMIR_POLDER
   use mqc_czt_esp, only: esp_matrices, drinv_matrices, ddrinv_matrices
   use mqc_czt_efp_potential, only: from_gamess_ao_order, frozen_core
   implicit none
   private

   public :: fragment_basis
   public :: two_fragment_molecule
   public :: fragment_molecule
   public :: fragment_lmo
   public :: exchange_repulsion
   public :: dispersion_e6_damped
   public :: dispersion_e7_damped
   public :: dispersion_e8_damped
   public :: lmo_overlap
   public :: charge_transfer
   ! Slot counts of the two higher dispersion records, exposed because anything
   ! unpacking one needs the same number.
   public :: N_DQ_SLOTS, N_QQ_SLOTS

   real(dp), parameter :: DAMPING_MIN_OVERLAP = 1.0e-5_dp
      !! Orbital overlap at or below which dispersion is left undamped; the
      !! bound that justifies it is in `overlap_damping`.

   integer, parameter :: N_QUAD_PAIR = 9
      !! Cartesian index pairs of a quadrupole, `(j, k)` as one index.

   integer, parameter :: N_DQ_SLOTS = 27
      !! Slots in a `DIPOLE-QUADRUPOLE` record: a 3x3x3 tensor, last index fastest.

   integer, parameter :: N_QQ_SLOTS = 81
      !! Slots in the `LMOQQPOL` record: a 3x3x3x3 tensor, last index fastest.

contains

   subroutine fragment_basis(frag, basis, error)
      !! A basis object for one fragment, from the projection basis it carries
      !!
      !! Cartesian: the ordering maps in `mqc_czt_efp_potential` cover Cartesian s, p,
      !! d and f, and a spherical potential is refused at write time. **A GAMESS
      !! library potential read here would need that assumption revisited.**
      type(efp_fragment_t), intent(in) :: frag
      type(molecular_basis_type), intent(out) :: basis
      type(error_t), intent(inout) :: error

      integer :: natm, at, sh, k, count, first

      if (.not. frag%has_basis) then
         call error%set(ERROR_VALIDATION, "efp: this fragment carries no projection "// &
                        "basis, so no molecule can be built from it")
         return
      end if

      natm = frag%n_atoms
      call basis%allocate_elements(natm)
      do at = 1, natm
         ! How many shells this atom owns. The reader tags each shell with the
         ! atom it was read under, in file order.
         count = 0
         do sh = 1, frag%n_shells
            if (frag%shell_atom(sh) == at) count = count + 1
         end do
         if (count == 0) then
            call error%set(ERROR_VALIDATION, "efp: an atom of the projection basis "// &
                           "carries no shells")
            return
         end if
         ! The element, where the fragment knows it. The nearest integer to the
         ! charge is the same thing only for an atom presenting its whole
         ! nucleus, which a fragment cut across a bond does not.
         if (allocated(frag%element)) then
            basis%elements(at)%element = trim(element_number_to_symbol(frag%element(at)))
         else
            basis%elements(at)%element = trim(element_number_to_symbol(nint(frag%charge(at))))
         end if
         basis%elements(at)%angular_form = ANGULAR_FORM_CARTESIAN
         call basis%elements(at)%allocate_shells(count)
         count = 0
         do sh = 1, frag%n_shells
            if (frag%shell_atom(sh) /= at) cycle
            count = count + 1
            first = frag%shell_first(sh)
            call basis%elements(at)%shells(count)%allocate_arrays(frag%shell_nprim(sh))
            basis%elements(at)%shells(count)%ang_mom = frag%shell_l(sh)
            basis%elements(at)%shells(count)%nfunc = frag%shell_nprim(sh)
            do k = 1, frag%shell_nprim(sh)
               basis%elements(at)%shells(count)%exponents(k) = frag%prim_expo(first + k - 1)
               basis%elements(at)%shells(count)%coefficients(k) = frag%prim_coef(first + k - 1)
            end do
         end do
      end do
   end subroutine fragment_basis

   function charge_transfer(frag_a, frag_b, offset_a, offset_b, error) result(energy)
      !! Charge transfer between two fragments
      !!
      !! Li, Gordon and Jensen, J. Chem. Phys. 124, 214108 (2006): second-order
      !! perturbation theory for an electron moving from an occupied orbital of
      !! one fragment into a virtual of the other, with the other fragment's
      !! multipoles as the perturbing potential `V`. For a donor `A`, an
      !! acceptor `B`, `i` occupied on `A`, `n` virtual on `B`, `j` occupied on
      !! `B` and `m` running over all of `A`'s orbitals, the virtual is first
      !! made orthogonal to `A`,
      !!
      !!     |n~> = |n> - sum_m S_mn |m>,      <n~|n~> = 1 - sum_m S_mn^2
      !!
      !! and then
      !!
      !!     E(A->B) = sum_in  V_in~ ( V_in~ + sum_j S_ij T_n~j )
      !!                       / ( <n~|n~> (e_i - T_nn) )
      !!
      !!     E_CT = 2 E(A->B) + 2 E(B->A)
      !!
      !! with `S` the overlap, `T` the kinetic energy, `e_i` the donor's orbital
      !! energy and the two closed shells supplying the factors of two. The
      !! acceptor's energy enters through its kinetic diagonal `T_nn`.
      !!
      !! The orbitals are the file's `CTVEC` -- occupied *and* virtual, since this
      !! moves density into the latter -- and `e_i` is its `CTFOK`.
      !!
      !! `V` is the other fragment's multipole potential over these orbitals:
      !! charges, dipoles and quadrupoles, and no octupole. Ending the rank at
      !! the quadrupole matches GAMESS, whose answer an octupole section does
      !! not change.
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      type(error_t), intent(inout) :: error
      real(dp) :: energy

      real(dp), allocatable :: s_ao(:, :), t_ao(:, :), ct_a(:, :), ct_b(:, :)
      real(dp), allocatable :: v_from_b(:, :), v_from_a(:, :)
      type(czt_molecule_t) :: pair
      integer :: n_ao_a

      energy = 0.0_dp
      if (.not. (frag_a%has_ctvec .and. frag_b%has_ctvec)) then
         call error%set(ERROR_VALIDATION, "efp: charge transfer needs CTVEC from both "// &
                        "fragments")
         return
      end if
      if (.not. (frag_a%has_ctfok .and. frag_b%has_ctfok)) then
         call error%set(ERROR_VALIDATION, "efp: charge transfer needs CTFOK from both "// &
                        "fragments")
         return
      end if

      call two_fragment_molecule(frag_a, frag_b, offset_a, offset_b, pair, n_ao_a, error)
      if (error%has_error()) return
      call pair%overlap(s_ao)
      call pair%kinetic(t_ao)
      call padded_ctvec(frag_a, pair, 0, n_ao_a, ct_a, error)
      if (error%has_error()) return
      call padded_ctvec(frag_b, pair, n_ao_a, pair%nao - n_ao_a, ct_b, error)
      if (error%has_error()) return
      call multipole_potential(pair, frag_b, offset_b, v_from_b, error)
      if (error%has_error()) return
      call multipole_potential(pair, frag_a, offset_a, v_from_a, error)
      if (error%has_error()) return

      energy = 2.0_dp*one_direction(frag_a, frag_b, ct_a, ct_b, s_ao, t_ao, v_from_b) &
               + 2.0_dp*one_direction(frag_b, frag_a, ct_b, ct_a, s_ao, t_ao, v_from_a)

      call pair%destroy()
      deallocate (s_ao, t_ao, ct_a, ct_b, v_from_b, v_from_a)
   end function charge_transfer

   function one_direction(frag_a, frag_b, ct_a, ct_b, s_ao, t_ao, v_b) result(ct)
      !! `E(A->B)` of `charge_transfer`: `a`'s occupied orbitals into `b`'s virtuals
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: ct_a(:, :), ct_b(:, :)
      real(dp), intent(in) :: s_ao(:, :), t_ao(:, :), v_b(:, :)
      real(dp) :: ct

      real(dp), allocatable :: s_ab(:, :), t_ab(:, :), t_bb(:, :), v_ab(:, :), v_aa(:, :)
      real(dp), allocatable :: s_to_vir(:, :), coupling(:, :), kinetic(:, :)
      real(dp), allocatable :: t_vir_occ(:, :), norm_vir(:), t_vir_diag(:)
      integer :: n_occ_a, n_occ_b, n_mo_b, i, n

      n_occ_a = frag_a%n_occ_ct
      n_occ_b = frag_b%n_occ_ct
      n_mo_b = frag_b%n_mo_ct

      s_ab = project(ct_a, s_ao, ct_b)
      t_ab = project(ct_a, t_ao, ct_b)
      t_bb = project(ct_b, t_ao, ct_b)
      v_ab = project(ct_a, v_b, ct_b)
      v_aa = project(ct_a, v_b, ct_a)

      ! Overlaps of every orbital of a with b's virtuals: (n_mo_a, n_vir_b).
      s_to_vir = s_ab(:, n_occ_b + 1:n_mo_b)

      ! <i|V|n~>, (n_occ_a, n_vir_b).
      coupling = v_ab(1:n_occ_a, n_occ_b + 1:n_mo_b) &
                 - matmul(v_aa(1:n_occ_a, :), s_to_vir)

      ! <n~|T|j>, (n_vir_b, n_occ_b), then its overlap-weighted sum over j.
      t_vir_occ = t_bb(n_occ_b + 1:n_mo_b, 1:n_occ_b) &
                  - matmul(transpose(s_to_vir), t_ab(:, 1:n_occ_b))
      kinetic = matmul(s_ab(1:n_occ_a, 1:n_occ_b), transpose(t_vir_occ))

      norm_vir = 1.0_dp - sum(s_to_vir**2, dim=1)
      t_vir_diag = [(t_bb(n, n), n=n_occ_b + 1, n_mo_b)]

      ct = 0.0_dp
      do n = 1, n_mo_b - n_occ_b
         do i = 1, n_occ_a
            ct = ct + coupling(i, n)*(coupling(i, n) + kinetic(i, n)) &
                 /(norm_vir(n)*(frag_a%eps_occ(i) - t_vir_diag(n)))
         end do
      end do
   end function one_direction

   function project(left, matrix, right) result(out)
      !! `left^T matrix right`, both sets already in the pair's AO space
      real(dp), intent(in) :: left(:, :), matrix(:, :), right(:, :)
      real(dp), allocatable :: out(:, :)

      real(dp), allocatable :: work(:, :)

      allocate (work(size(matrix, 1), size(right, 2)))
      allocate (out(size(left, 2), size(right, 2)))
      call pic_gemm(matrix, right, work)
      call pic_gemm(left, work, out, transa="T")
      deallocate (work)
   end function project

   subroutine multipole_potential(pair, frag, offset, v, error)
      !! The potential of one fragment's multipoles over the pair's basis
      !!
      !! Charges, dipoles and quadrupoles, and **no octupole** -- one rank fewer
      !! than the electrostatic energy uses.
      !!
      !! Every rank carries a minus, because this is the potential energy of an
      !! *electron* in the multipole's field. The ranks are successive
      !! centre-gradients of one integral:
      !!
      !!     charge      -        q     <mu| 1/r_C |nu>
      !!     dipole      -        d_a   grad_a <mu| 1/r_C |nu>
      !!     quadrupole  - (1/3)  Q_ab  grad_a grad_b <mu| 1/r_C |nu>
      !!
      !! all gradients with respect to `C`, and the quadrupole one summed over all
      !! nine `ab` (Li, Gordon and Jensen, J. Chem. Phys. 124, 214108 (2006)).
      !!
      !! The dipole and quadrupole moments here are **electronic only**; the
      !! nucleus sits entirely in the monopole, which is how the potential stores
      !! them.
      type(czt_molecule_t), intent(in) :: pair
      type(efp_fragment_t), intent(in) :: frag
      real(dp), intent(in) :: offset(3)
      real(dp), allocatable, intent(out) :: v(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: grids(:, :), matrices(:, :, :), grad(:, :, :, :)
      real(dp), allocatable :: second(:, :, :, :, :)
      real(dp) :: quad(3, 3)
      integer :: k, x, y

      allocate (grids(3, frag%n_points))
      do k = 1, frag%n_points
         grids(:, k) = frag%points(:, k) + offset
      end do
      call esp_matrices(pair, grids, matrices, error)
      if (error%has_error()) return
      allocate (v(pair%nao, pair%nao))
      v = 0.0_dp
      do k = 1, frag%n_points
         ! Minus: this is the potential energy of an *electron* in the field of that
         ! monopole, and the electron's charge is negative.
         v = v - (frag%q_elec(k) + frag%q_nuc(k))*matrices(:, :, k)
      end do
      deallocate (matrices)

      ! The dipole rank. `drinv_matrices` is `grad_C` of the integral above, its
      ! sign fixed in `test_mqc_czt_esp` against a difference of `esp_matrices`:
      ! libcint's `nabla-rinv` naming does not say whether the gradient is with
      ! respect to `r` or to `C`.
      call drinv_matrices(pair, grids, grad, error)
      if (error%has_error()) return
      do k = 1, frag%n_points
         do x = 1, 3
            v = v - frag%dipole(x, k)*grad(:, :, x, k)
         end do
      end do
      deallocate (grad)

      ! The quadrupole rank, summed over all nine `ab` of a symmetric `Q`.
      !
      ! The `1/3` is the coefficient a *traceless* quadrupole carries against
      ! `grad grad (1/R)`: the expansion is `(1/2) M_ab grad_a grad_b` over the
      ! second moments, and substituting `Q = (3M - tr M)/2` turns the half into
      ! a third, because `delta_ab grad_a grad_b (1/R)` vanishes.
      call ddrinv_matrices(pair, grids, second, error)
      if (error%has_error()) return
      do k = 1, frag%n_points
         quad = traceless_quadrupole(frag%quadrupole(:, k))
         do y = 1, 3
            do x = 1, 3
               v = v - quad(x, y)*second(:, :, x, y, k)/3.0_dp
            end do
         end do
      end do
      deallocate (grids, second)
   end subroutine multipole_potential

   pure function traceless_quadrupole(stored) result(quad)
      !! The stored second moments as a traceless electric quadrupole
      !!
      !! A potential's `QUADRUPOLES` section holds **second moments, not
      !! traceless quadrupoles**, as a nonzero stored trace shows.
      !!
      !!     Q_aa = (3 M_aa - tr M) / 2        Q_ab = (3/2) M_ab,  a /= b
      !!
      !! The tracelessness is not cosmetic here. `grad_a grad_b (1/r_C)` is traceless
      !! away from `C` but carries a delta function at it, so a trace left in `Q` would
      !! contribute through that contact term rather than cancelling.
      real(dp), intent(in) :: stored(N_QUADRUPOLE)   !! `XX YY ZZ XY XZ YZ`
      real(dp) :: quad(3, 3)

      real(dp) :: trace

      trace = stored(1) + stored(2) + stored(3)
      quad(1, 1) = (3.0_dp*stored(1) - trace)*0.5_dp
      quad(2, 2) = (3.0_dp*stored(2) - trace)*0.5_dp
      quad(3, 3) = (3.0_dp*stored(3) - trace)*0.5_dp
      quad(1, 2) = 1.5_dp*stored(4)
      quad(1, 3) = 1.5_dp*stored(5)
      quad(2, 3) = 1.5_dp*stored(6)
      quad(2, 1) = quad(1, 2)
      quad(3, 1) = quad(1, 3)
      quad(3, 2) = quad(2, 3)
   end function traceless_quadrupole

   subroutine padded_ctvec(frag, pair, offset_ao, n_own, padded, error)
      !! `CTVEC` in our AO order, padded into the pair's space
      type(efp_fragment_t), intent(in) :: frag
      type(czt_molecule_t), intent(in) :: pair
      integer, intent(in) :: offset_ao, n_own
      real(dp), allocatable, intent(out) :: padded(:, :)
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: own
      real(dp), allocatable :: ct(:, :)

      call fragment_molecule(frag, [0.0_dp, 0.0_dp, 0.0_dp], own, error)
      if (error%has_error()) return
      call from_gamess_ao_order(own, frag%ctvec_gamess, ct, error)
      call own%destroy()
      if (error%has_error()) return
      allocate (padded(pair%nao, frag%n_mo_ct))
      padded = 0.0_dp
      padded(offset_ao + 1:offset_ao + n_own, :) = ct
      deallocate (ct)
   end subroutine padded_ctvec

   subroutine lmo_overlap(frag_a, frag_b, offset_a, offset_b, s_lmo, t_lmo, error)
      !! Overlap and kinetic energy between two fragments' localized orbitals
      !!
      !! `(n_lmo_proj_a, n_lmo_proj_b)`. Both exchange repulsion and the dispersion
      !! damping are built from this.
      ! TODO(mqc): `t_lmo` is only used by the exchange repulsion; all three
      ! dispersion routines allocate it, pay for the kinetic integrals and two
      ! gemms, and discard it.
      ! TODO(mqc): the shapes are `n_lmo_proj`, from the projection wavefunction,
      ! while the dispersion callers index them by `n_lmo`, from the dynamic
      ! polarizability block. Nothing checks that the two sections of a potential
      ! agree, so a file where they differ reads out of bounds.
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      real(dp), allocatable, intent(out) :: s_lmo(:, :), t_lmo(:, :)
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: pair
      real(dp), allocatable :: s_ao(:, :), t_ao(:, :), lmo_a(:, :), lmo_b(:, :)
      real(dp), allocatable :: work(:, :)
      integer :: n_ao_a

      call two_fragment_molecule(frag_a, frag_b, offset_a, offset_b, pair, n_ao_a, error)
      if (error%has_error()) return
      call pair%overlap(s_ao)
      call pair%kinetic(t_ao)
      call padded_lmo(frag_a, pair, 0, n_ao_a, lmo_a, error)
      if (error%has_error()) return
      call padded_lmo(frag_b, pair, n_ao_a, pair%nao - n_ao_a, lmo_b, error)
      if (error%has_error()) return

      allocate (work(pair%nao, frag_b%n_lmo_proj))
      allocate (s_lmo(frag_a%n_lmo_proj, frag_b%n_lmo_proj))
      allocate (t_lmo(frag_a%n_lmo_proj, frag_b%n_lmo_proj))
      call pic_gemm(s_ao, lmo_b, work)
      call pic_gemm(lmo_a, work, s_lmo, transa="T")
      call pic_gemm(t_ao, lmo_b, work)
      call pic_gemm(lmo_a, work, t_lmo, transa="T")

      call pair%destroy()
      deallocate (s_ao, t_ao, lmo_a, lmo_b, work)
   end subroutine lmo_overlap

   function dispersion_e6_damped(frag_a, frag_b, offset_a, offset_b, error) result(energy)
      !! `E6` with overlap-based damping, which is what GAMESS reports
      !!
      !! Each orbital pair is damped by `F6` of `overlap_damping`, a function of
      !! that pair's own overlap (Slipchenko and Gordon, Mol. Phys. 107, 999
      !! (2009)).
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      type(error_t), intent(inout) :: error
      real(dp) :: energy

      real(dp), allocatable :: s_lmo(:, :), t_lmo(:, :)
      real(dp) :: sep(3), weight(N_CASIMIR_POLDER)
      real(dp) :: c6, dist, f6, f7, f8
      integer :: i, j, k

      energy = 0.0_dp
      weight = casimir_polder_weights()
      if (.not. (frag_a%has_dynamic .and. frag_b%has_dynamic)) then
         call error%set(ERROR_VALIDATION, "efp: damped E6 needs the dynamic "// &
                        "polarizabilities of both fragments")
         return
      end if
      call lmo_overlap(frag_a, frag_b, offset_a, offset_b, s_lmo, t_lmo, error)
      if (error%has_error()) return

      do i = 1, frag_a%n_lmo
         do j = 1, frag_b%n_lmo
            c6 = 0.0_dp
            do k = 1, min(frag_a%n_freq, frag_b%n_freq)
               c6 = c6 + weight(k)*isotropic_alpha(frag_a%dyn_pol(:, :, i, k)) &
                    *isotropic_alpha(frag_b%dyn_pol(:, :, j, k))
            end do
            c6 = c6*3.0_dp/PI
            sep = frag_b%centroids(:, j) + offset_b - frag_a%centroids(:, i) - offset_a
            dist = sqrt(sep(1)*sep(1) + sep(2)*sep(2) + sep(3)*sep(3))

            call overlap_damping(s_lmo(i, j), f6, f7, f8)
            energy = energy - c6*f6/dist**6
         end do
      end do
      deallocate (s_lmo, t_lmo)
   end function dispersion_e6_damped

   pure subroutine overlap_damping(s, f6, f7, f8)
      !! Dispersion damping from an orbital pair's overlap
      !!
      !! Slipchenko and Gordon, Mol. Phys. 107, 999 (2009). Treating the two
      !! orbitals' overlap as that of two spherical exponentials puts the
      !! Tang-Toennies damping variable at `x = sqrt(-2 ln|S|)`, so that
      !! `exp(-x^2) = S^2`, and
      !!
      !!     F_n = 1 - S^2 sum_{k=0..n} x^k / k!,       n = 6, 7, 8
      !!
      !! **The series is in `x`, the square root**, so most of its terms carry
      !! half-integer powers of `-2 ln|S|`; an ordinary exponential series in
      !! `-2 ln|S|` would still damp, plausibly and wrongly. Each `F_n` drops one
      !! more positive term, so `F8 < F7 < F6 < 1` wherever they are evaluated.
      !!
      !! At `|S| <= DAMPING_MIN_OVERLAP` all three are one. There `x >= 4.80`
      !! and `S^2 <= 1e-10`, and `S^2 sum x^k/k!` only falls as `|S|` does, so
      !! `1 - F8` is below `1.2e-8` -- smaller than the precision dispersion is
      !! known to -- and the cutoff also keeps `ln|S|` away from zero overlap.
      real(dp), intent(in) :: s
      real(dp), intent(out) :: f6, f7, f8

      integer, parameter :: HIGHEST_ORDER = 8
      real(dp) :: x, term
      real(dp) :: partial(0:HIGHEST_ORDER)
         !! `partial(n) = sum_{k=0..n} x^k/k!`
      integer :: k

      f6 = 1.0_dp
      f7 = 1.0_dp
      f8 = 1.0_dp
      if (abs(s) <= DAMPING_MIN_OVERLAP) return

      x = sqrt(-2.0_dp*log(abs(s)))
      term = 1.0_dp
      partial(0) = 1.0_dp
      do k = 1, HIGHEST_ORDER
         term = term*x/real(k, dp)
         partial(k) = partial(k - 1) + term
      end do
      f6 = 1.0_dp - s*s*partial(6)
      f7 = 1.0_dp - s*s*partial(7)
      f8 = 1.0_dp - s*s*partial(8)
   end subroutine overlap_damping

   function dispersion_e8_damped(frag_a, frag_b, offset_a, offset_b, error) result(energy)
      !! `E8`, the isotropic one -- which is the one GAMESS prints
      !!
      !! GAMESS also forms an anisotropic `E8` but never prints it. What reaches
      !! `E8 DISPERSION ENERGY` is the isotropic form, computed whenever the
      !! potential carries an `LMOQQPOL` section at all:
      !!
      !!     C8 = (15/pi) sum_k w_k ( a_iso^A(k) A_QQ^B(k) + a_iso^B(k) A_QQ^A(k) )
      !!     E8 = - F8 C8 / R^8
      !!
      !! with `w_k` the Casimir-Polder quadrature weights.
      !!
      !! It therefore needs the trace of the dipole-dipole polarizability and the
      !! spherical average of the quadrupole-quadrupole one, and nothing else. Both
      !! averages are isotropic contractions and survive the rotation into the
      !! current frame unchanged, which is why nothing is rotated here.
      !!
      !! **Not `E6/3`.** That approximation is a separate option in GAMESS and
      !! does not run for a pair of file-based fragments; reconciling against it
      !! would fit a factor to the wrong quantity.
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      type(error_t), intent(inout) :: error
      real(dp) :: energy

      real(dp), allocatable :: s_lmo(:, :), t_lmo(:, :)
      real(dp) :: sep(3), weight(N_CASIMIR_POLDER)
      real(dp) :: c8, dist, f6, f7, f8
      integer :: i, j, k

      energy = 0.0_dp
      weight = casimir_polder_weights()
      if (.not. (frag_a%has_dynamic .and. frag_b%has_dynamic)) then
         call error%set(ERROR_VALIDATION, "efp: damped E8 needs the dynamic "// &
                        "polarizabilities of both fragments")
         return
      end if
      if (.not. (frag_a%has_quadquad .and. frag_b%has_quadquad)) then
         call error%set(ERROR_VALIDATION, "efp: damped E8 needs the "// &
                        "quadrupole-quadrupole polarizabilities of both fragments")
         return
      end if
      if (frag_a%n_quadquad /= N_QQ_SLOTS .or. frag_b%n_quadquad /= N_QQ_SLOTS) then
         call error%set(ERROR_VALIDATION, "efp: a quadrupole-quadrupole record "// &
                        "does not carry 81 values")
         return
      end if
      call lmo_overlap(frag_a, frag_b, offset_a, offset_b, s_lmo, t_lmo, error)
      if (error%has_error()) return

      do i = 1, frag_a%n_lmo
         do j = 1, frag_b%n_lmo
            c8 = 0.0_dp
            do k = 1, min(frag_a%n_freq, frag_b%n_freq)
               c8 = c8 + weight(k) &
                    *(isotropic_alpha(frag_a%dyn_pol(:, :, i, k)) &
                      *isotropic_quadquad(frag_b%quadquad(:, j, k)) &
                      + isotropic_alpha(frag_b%dyn_pol(:, :, j, k)) &
                      *isotropic_quadquad(frag_a%quadquad(:, i, k)))
            end do
            c8 = c8*15.0_dp/PI
            sep = frag_b%centroids(:, j) + offset_b - frag_a%centroids(:, i) - offset_a
            dist = sqrt(sep(1)*sep(1) + sep(2)*sep(2) + sep(3)*sep(3))

            call overlap_damping(s_lmo(i, j), f6, f7, f8)
            energy = energy - c8*f8/dist**8
         end do
      end do
      deallocate (s_lmo, t_lmo)
   end function dispersion_e8_damped

   function dispersion_e7_damped(frag_a, frag_b, offset_a, offset_b, error) result(energy)
      !! `E7`, the dipole-dipole/dipole-quadrupole cross term of the dispersion
      !!
      !! Xu, Zahariev and Gordon, J. Chem. Theory Comput. 10, 1576 (2014). The
      !! seventh-order term pairs a dipole-dipole polarizability on one orbital
      !! with a dipole-quadrupole one on the other, through one dipole-dipole and
      !! one dipole-quadrupole interaction tensor. Per orbital pair, over the
      !! Casimir-Polder frequencies `nu_k` with weights `w_k`,
      !!
      !!     E7 = -F7/(3 pi) sum_k w_k sum_abcde T2_ab T3_cde
      !!            [ alpha^A_ca(k) A^B_b,de(k) - alpha^B_eb(k) A^A_a,cd(k) ]
      !!
      !! with `T2 = grad grad (1/R)` and `T3 = -grad grad grad (1/R)` at the
      !! separation `R`, `alpha` as stored (`(dipole, field)`), `A_x,yz` a
      !! dipole-quadrupole tensor with its dipole index first and `F7` from
      !! `overlap_damping`. There is only one E7 -- no isotropic variant.
      !!
      !! **Three conventions here are sign- or transpose-critical**, and E7 is the
      !! first term that can see any of them; E6 and E8 reach the polarizabilities
      !! only through isotropic averages and the separation only through `R`.
      !!
      !! * *The separation runs A minus B*, A being this routine's first
      !!   fragment. E7 is odd in it: `T2` is even and `T3` is odd.
      !! * *`T3` is the negative of the third derivative.* With the A-minus-B
      !!   separation that is the pairing that gives E7 GAMESS's sign; the
      !!   textbook `+grad grad grad` gives the opposite one.
      !! * *Each fragment's `T2` slot takes the second (field) index of its
      !!   `alpha`, and a `T3` slot the first.* The tensor is not symmetric -- its
      !!   antisymmetric part is 12% of it here -- so the swap is load-bearing.
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      type(error_t), intent(inout) :: error
      real(dp) :: energy

      real(dp), allocatable :: s_lmo(:, :), t_lmo(:, :)
      real(dp) :: sep(3), t2(3, 3), t3(3, 3, 3), weight(N_CASIMIR_POLDER)
      real(dp) :: t3_first(3, N_QUAD_PAIR), t3_last(N_QUAD_PAIR, 3)
      real(dp) :: via_b(3, 3), via_a(3, 3)
      real(dp) :: dq_a(3, N_QUAD_PAIR), dq_b(3, N_QUAD_PAIR)
      real(dp) :: summed, f6, f7, f8
      integer :: i, j, k

      energy = 0.0_dp
      if (.not. (frag_a%has_dynamic .and. frag_b%has_dynamic)) then
         call error%set(ERROR_VALIDATION, "efp: damped E7 needs the dynamic "// &
                        "polarizabilities of both fragments")
         return
      end if
      if (.not. (frag_a%has_dipquad .and. frag_b%has_dipquad)) then
         call error%set(ERROR_VALIDATION, "efp: damped E7 needs the "// &
                        "dipole-quadrupole polarizabilities of both fragments")
         return
      end if
      if (frag_a%n_dipquad /= N_DQ_SLOTS .or. frag_b%n_dipquad /= N_DQ_SLOTS) then
         call error%set(ERROR_VALIDATION, "efp: a dipole-quadrupole record does "// &
                        "not carry 27 values")
         return
      end if
      call lmo_overlap(frag_a, frag_b, offset_a, offset_b, s_lmo, t_lmo, error)
      if (error%has_error()) return
      weight = casimir_polder_weights()

      do i = 1, frag_a%n_lmo
         do j = 1, frag_b%n_lmo
            sep = frag_a%centroids(:, i) + offset_a - frag_b%centroids(:, j) - offset_b
            call t_tensors(sep, t2, t3)
            ! T3 with one index apart and the other two as a column-major pair,
            ! in each of the two places that index is needed.
            t3_first = reshape(t3, [3, N_QUAD_PAIR])
            t3_last = reshape(t3, [N_QUAD_PAIR, 3])

            summed = 0.0_dp
            do k = 1, min(frag_a%n_freq, frag_b%n_freq)
               ! B's tensor through both interaction tensors, left open on the
               ! two indices A's alpha takes: via_b(a, c) =
               ! sum_bde T2_ab T3_cde A^B_b,de.
               dq_a = dq_pairs(frag_a%dipquad(:, i, k))
               dq_b = dq_pairs(frag_b%dipquad(:, j, k))
               via_b = matmul(t2, transpose(matmul(t3_first, transpose(dq_b))))
               ! A's likewise, open on B's two: via_a(b, e) =
               ! sum_acd T2_ab T3_cde A^A_a,cd.
               via_a = matmul(transpose(t2), matmul(dq_a, t3_last))
               summed = summed + weight(k) &
                        *(sum(frag_a%dyn_pol(:, :, i, k)*transpose(via_b)) &
                          - sum(frag_b%dyn_pol(:, :, j, k)*transpose(via_a)))
            end do

            call overlap_damping(s_lmo(i, j), f6, f7, f8)
            energy = energy - f7*summed/(3.0_dp*PI)
         end do
      end do
      deallocate (s_lmo, t_lmo)
   end function dispersion_e7_damped

   pure function dq_pairs(record) result(tensor)
      !! A `DIPOLE-QUADRUPOLE` record as a `(dipole, quadrupole pair)` matrix
      !!
      !! Column `j + 3 (k - 1)` holds quadrupole indices `(j, k)`, the same
      !! column-major pairing `reshape` gives a rank-3 `T3`.
      real(dp), intent(in) :: record(N_DQ_SLOTS)
      real(dp) :: tensor(3, N_QUAD_PAIR)

      integer :: dip, j, k

      do k = 1, 3
         do j = 1, 3
            do dip = 1, 3
               tensor(dip, j + 3*(k - 1)) = record(dq_slot(dip, j, k))
            end do
         end do
      end do
   end function dq_pairs

   pure subroutine t_tensors(sep, t2, t3)
      !! The rank-2 and rank-3 interaction tensors at separation `sep`
      !!
      !!     T2_ij  =   d_i d_j (1/R)  = ( 3 R_i R_j - R^2 delta_ij ) / R^5
      !!     T3_ijk = - d_i d_j d_k (1/R)
      !!            = ( 15 R_i R_j R_k
      !!                - 3 R^2 ( R_i delta_jk + R_j delta_ik + R_k delta_ij ) ) / R^7
      !!
      !! `T2` is Stone's (*The Theory of Intermolecular Forces*) and `T3` the
      !! negative of his; `dispersion_e7_damped` says why that sign is wanted.
      real(dp), intent(in) :: sep(3)
      real(dp), intent(out) :: t2(3, 3), t3(3, 3, 3)

      real(dp) :: eye(3, 3), rr(3, 3)
      real(dp) :: r2, r5, r7
      integer :: k

      eye = 0.0_dp
      do k = 1, 3
         eye(k, k) = 1.0_dp
      end do
      r2 = dot_product(sep, sep)
      r5 = r2*r2*sqrt(r2)
      r7 = r5*r2
      rr = outer(sep, sep)

      t2 = (3.0_dp*rr - r2*eye)/r5
      do k = 1, 3
         ! R_i delta_jk and R_j delta_ik are `sep` against the k-th unit
         ! vector, both ways round; R_k delta_ij is the identity.
         t3(:, :, k) = (15.0_dp*sep(k)*rr &
                        - 3.0_dp*r2*(outer(sep, eye(:, k)) + outer(eye(:, k), sep) &
                                     + sep(k)*eye))/r7
      end do
   end subroutine t_tensors

   pure function outer(x, y) result(xy)
      !! `x_i y_j`
      real(dp), intent(in) :: x(3), y(3)
      real(dp) :: xy(3, 3)

      xy = spread(x, 2, 3)*spread(y, 1, 3)
   end function outer

   pure function dq_slot(i, j, k) result(slot)
      !! Where `DQ(i,j,k)` sits in a `DIPOLE-QUADRUPOLE` record
      !!
      !! `i` is the dipole index and `(j,k)` the quadrupole pair. Row-major with the
      !! last index fastest, as the file format lays it out.
      integer, intent(in) :: i, j, k
      integer :: slot

      slot = (i - 1)*9 + (j - 1)*3 + k
   end function dq_slot

   pure function isotropic_quadquad(values) result(a)
      !! The spherical average of a quadrupole-quadrupole tensor
      !!
      !! The full rank-four tensor contracted against the isotropic projector
      !!
      !!     A_QQ = (1/5) sum_ijkl QQ(i,j,k,l)
      !!            [ (d_ik d_jl + d_il d_jk)/2 - d_ij d_kl / 3 ]
      !!
      !! which with the deltas resolved is the two sums below. The projector
      !! symmetrises in `(i,j) <-> (k,l)`, so whether the stored tensor is exactly
      !! symmetric under that swap does not affect this number.
      !!
      !! The 81 file slots are plain row-major with the last index fastest, unlike
      !! the nine dipole-dipole slots, which are permuted.
      real(dp), intent(in) :: values(N_QQ_SLOTS)
      real(dp) :: a

      real(dp) :: paired, traced
      integer :: i, j

      paired = 0.0_dp
      traced = 0.0_dp
      do i = 1, 3
         do j = 1, 3
            paired = paired + 0.5_dp*(values(qq_slot(i, j, i, j)) &
                                      + values(qq_slot(i, j, j, i)))
            traced = traced + values(qq_slot(i, i, j, j))
         end do
      end do
      a = (paired - traced/3.0_dp)/5.0_dp
   end function isotropic_quadquad

   pure function qq_slot(i, j, k, l) result(slot)
      !! Where `QQ(i,j,k,l)` sits in an `LMOQQPOL` record
      integer, intent(in) :: i, j, k, l
      integer :: slot

      slot = (i - 1)*27 + (j - 1)*9 + (k - 1)*3 + l
   end function qq_slot

   pure function isotropic_alpha(tensor) result(a)
      !! A third of the trace
      real(dp), intent(in) :: tensor(3, 3)
      real(dp) :: a

      a = (tensor(1, 1) + tensor(2, 2) + tensor(3, 3))/3.0_dp
   end function isotropic_alpha

   function exchange_repulsion(frag_a, frag_b, offset_a, offset_b, error) result(energy)
      !! Pauli exchange repulsion between two fragments
      !!
      !! Jensen and Gordon, Mol. Phys. 89, 1313 (1996) and J. Chem. Phys. 108,
      !! 4772 (1998): the exchange energy of two closed-shell fragments expanded
      !! to second order in the intermolecular overlap of their localized
      !! orbitals, `i` on A and `j` on B, with the two-electron integrals
      !! replaced by the overlap and the orbitals' centroids:
      !!
      !!     E = -4 sqrt(2/pi) sum_ij sqrt(-ln|S_ij|) S_ij^2 / R_ij
      !!         -2 sum_ij S_ij ( sum_k F^A_ik S_kj + sum_l F^B_jl S_il - 2 T_ij )
      !!         +2 sum_ij S_ij^2 ( V^B_i + V^A_j - 1/R_ij )
      !!
      !! `S` and `T` are the overlap and kinetic energy between the two sets,
      !! `F` each fragment's Fock matrix over its own orbitals, `R_ij` the
      !! distance between centroids, and `V^B_i` the electrostatic potential
      !! of B at `i`'s centroid -- B's valence nuclear charges and its doubly
      !! occupied orbitals, each orbital taken as two electrons at its centroid
      !! (`fragment_potential`).
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      type(error_t), intent(inout) :: error
      real(dp) :: energy

      real(dp), parameter :: LOG_MIN_OVERLAP = 1.0e-10_dp
         !! Overlap at or below which the `sqrt(-ln|S|)` term is dropped. The
         !! term tends to zero with `S`, but `ln 0` does not exist and an exact
         !! zero overlap is common, from symmetry or from orbitals too far apart
         !! for their overlap to be representable. At `|S| = 1e-10` the term is
         !! `4 sqrt(2/pi) sqrt(23) S^2 / R < 2e-19/R` Hartree, so nothing below
         !! the cutoff can reach the sum.
      type(czt_molecule_t) :: pair
      real(dp), allocatable :: s_ao(:, :), t_ao(:, :), lmo_a(:, :), lmo_b(:, :)
      real(dp), allocatable :: s_lmo(:, :), t_lmo(:, :), work(:, :)
      real(dp), allocatable :: cen_a(:, :), cen_b(:, :), v_at_a(:), v_at_b(:)
      real(dp), allocatable :: dist(:, :), fock_term(:, :)
      real(dp) :: overlap_log, kinetic_fock, potential_term
      integer :: n_ao_a, i, j, n_lmo_a, n_lmo_b

      energy = 0.0_dp
      if (.not. (frag_a%has_fock .and. frag_b%has_fock)) then
         call error%set(ERROR_VALIDATION, "efp: exchange repulsion needs the LMO Fock "// &
                        "matrix of both fragments")
         return
      end if

      call two_fragment_molecule(frag_a, frag_b, offset_a, offset_b, pair, n_ao_a, error)
      if (error%has_error()) return
      call pair%overlap(s_ao)
      call pair%kinetic(t_ao)

      ! Each fragment's orbitals live on its own functions, so they are padded into
      ! the pair's space rather than transformed with a rectangular block: the
      ! transform is then one gemm and the block boundary appears once, here.
      n_lmo_a = frag_a%n_lmo_proj
      n_lmo_b = frag_b%n_lmo_proj
      call padded_lmo(frag_a, pair, 0, n_ao_a, lmo_a, error)
      if (error%has_error()) return
      call padded_lmo(frag_b, pair, n_ao_a, pair%nao - n_ao_a, lmo_b, error)
      if (error%has_error()) return

      allocate (work(pair%nao, n_lmo_b), s_lmo(n_lmo_a, n_lmo_b), t_lmo(n_lmo_a, n_lmo_b))
      call pic_gemm(s_ao, lmo_b, work)
      call pic_gemm(lmo_a, work, s_lmo, transa="T")
      call pic_gemm(t_ao, lmo_b, work)
      call pic_gemm(lmo_a, work, t_lmo, transa="T")

      cen_a = frag_a%pol_points(:, 1:n_lmo_a) + spread(offset_a, 2, n_lmo_a)
      cen_b = frag_b%pol_points(:, 1:n_lmo_b) + spread(offset_b, 2, n_lmo_b)

      ! Each centroid's potential from the other fragment, once.
      allocate (v_at_a(n_lmo_a), v_at_b(n_lmo_b), dist(n_lmo_a, n_lmo_b))
      do i = 1, n_lmo_a
         v_at_a(i) = fragment_potential(cen_a(:, i), frag_b, offset_b, cen_b)
      end do
      do j = 1, n_lmo_b
         v_at_b(j) = fragment_potential(cen_b(:, j), frag_a, offset_a, cen_a)
      end do
      do j = 1, n_lmo_b
         do i = 1, n_lmo_a
            dist(i, j) = norm2(cen_a(:, i) - cen_b(:, j))
         end do
      end do

      ! (F^A S + S F^B^T - 2T)_ij, the bracket of the second term.
      fock_term = matmul(frag_a%fock_lmo(1:n_lmo_a, 1:n_lmo_a), s_lmo) &
                  + matmul(s_lmo, transpose(frag_b%fock_lmo(1:n_lmo_b, 1:n_lmo_b))) &
                  - 2.0_dp*t_lmo
      kinetic_fock = -2.0_dp*sum(s_lmo*fock_term)

      potential_term = 0.0_dp
      overlap_log = 0.0_dp
      do j = 1, n_lmo_b
         do i = 1, n_lmo_a
            potential_term = potential_term &
                             + s_lmo(i, j)**2*(v_at_a(i) + v_at_b(j) - 1.0_dp/dist(i, j))
            if (abs(s_lmo(i, j)) > LOG_MIN_OVERLAP) then
               overlap_log = overlap_log &
                             + sqrt(-log(abs(s_lmo(i, j))))*s_lmo(i, j)**2/dist(i, j)
            end if
         end do
      end do

      energy = -4.0_dp*sqrt(2.0_dp/PI)*overlap_log + kinetic_fock + 2.0_dp*potential_term

      call pair%destroy()
      deallocate (s_ao, t_ao, lmo_a, lmo_b, s_lmo, t_lmo, work, cen_a, cen_b, &
                  v_at_a, v_at_b, dist, fock_term)
   end function exchange_repulsion

   function fragment_potential(point, frag, offset, centroids) result(v)
      !! The electrostatic potential at `point` of a fragment's valence nuclei
      !! and doubly occupied orbitals, each orbital as two electrons at its
      !! centroid
      real(dp), intent(in) :: point(3)
      type(efp_fragment_t), intent(in) :: frag
      real(dp), intent(in) :: offset(3)
      real(dp), intent(in) :: centroids(:, :)   !! (3, n), already placed
      real(dp) :: v

      real(dp) :: valence(frag%n_atoms)
      integer :: at, k

      ! Valence charge, not the full nuclear one: the localized orbitals a
      ! potential carries are valence only, so the core electrons count as
      ! screening their own nucleus. The same number the PROJECTION BASIS SET
      ! atom header carries. Taken from the fragment where it knows it,
      ! because across a cut bond the core does not follow from the charge:
      ! the centre standing in for a detached atom presents +1 and holds no
      ! core, whatever element it is.
      if (allocated(frag%valence)) then
         valence = frag%valence(1:frag%n_atoms)
      else
         do at = 1, frag%n_atoms
            valence(at) = frag%charge(at) &
                          - 2.0_dp*real(frozen_core([nint(frag%charge(at))]), dp)
         end do
      end if

      v = 0.0_dp
      do at = 1, frag%n_atoms
         v = v - valence(at)/norm2(point - frag%points(:, at) - offset)
      end do
      do k = 1, size(centroids, 2)
         v = v + 2.0_dp/norm2(point - centroids(:, k))
      end do
   end function fragment_potential

   subroutine padded_lmo(frag, pair, offset_ao, n_own, padded, error)
      !! One fragment's orbitals in the pair's full AO space, zero elsewhere
      type(efp_fragment_t), intent(in) :: frag
      type(czt_molecule_t), intent(in) :: pair
      integer, intent(in) :: offset_ao, n_own
      real(dp), allocatable, intent(out) :: padded(:, :)
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: own
      real(dp), allocatable :: lmo(:, :)

      call fragment_molecule(frag, [0.0_dp, 0.0_dp, 0.0_dp], own, error)
      if (error%has_error()) return
      if (own%nao /= n_own) then
         call error%set(ERROR_VALIDATION, "efp: a fragment's own basis does not match "// &
                        "its block of the pair")
         call own%destroy()
         return
      end if
      call fragment_lmo(frag, own, lmo, error)
      call own%destroy()
      if (error%has_error()) return
      allocate (padded(pair%nao, frag%n_lmo_proj))
      padded = 0.0_dp
      padded(offset_ao + 1:offset_ao + n_own, :) = lmo
      deallocate (lmo)
   end subroutine padded_lmo

   subroutine fragment_molecule(frag, offset, mol, error)
      !! One fragment as a molecule of its own, from the basis it carries
      !!
      !! Also the reference the pair is checked against: the pair's diagonal block
      !! has to be this molecule's overlap.
      type(efp_fragment_t), intent(in) :: frag
      real(dp), intent(in) :: offset(3)
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: error

      type(molecular_basis_type) :: basis
      real(dp), allocatable :: coords(:, :)
      integer, allocatable :: z(:)
      integer :: at

      call fragment_basis(frag, basis, error)
      if (error%has_error()) return
      allocate (coords(3, frag%n_atoms), z(frag%n_atoms))
      do at = 1, frag%n_atoms
         coords(:, at) = frag%points(:, at) + offset
         z(at) = nint(frag%charge(at))
      end do
      call mol%build(z, coords, basis, error, force_cartesian=.true.)
      deallocate (coords, z)
   end subroutine fragment_molecule

   subroutine fragment_lmo(frag, mol, lmo, error)
      !! The fragment's localized orbitals in *our* AO order
      !!
      !! The file stores them in GAMESS's, so the d and f permutation and its
      !! normalizations are undone by `from_gamess_ao_order`.
      !!
      !! The check that this is right is that the orbitals come back orthonormal
      !! against the fragment's own overlap, `C^T S C = I`. Nothing weaker would
      !! notice a permutation applied in the wrong direction.
      type(efp_fragment_t), intent(in) :: frag
      type(czt_molecule_t), intent(in) :: mol
      real(dp), allocatable, intent(out) :: lmo(:, :)
      type(error_t), intent(inout) :: error

      if (.not. frag%has_lmo) then
         call error%set(ERROR_VALIDATION, "efp: this fragment carries no projection "// &
                        "wavefunction")
         return
      end if
      if (frag%nao_proj /= mol%nao) then
         call error%set(ERROR_VALIDATION, "efp: the projection wavefunction and the "// &
                        "molecule disagree on the number of basis functions")
         return
      end if
      call from_gamess_ao_order(mol, frag%lmo_gamess, lmo, error)
   end subroutine fragment_lmo

   subroutine two_fragment_molecule(frag_a, frag_b, offset_a, offset_b, mol, &
                                    n_ao_a, error)
      !! One molecule covering both fragments, placed
      !!
      !! The atoms of `a` come first and then those of `b`, so the overlap matrix
      !! this yields is block structured: the leading `n_ao_a` rows and columns are
      !! `a`'s own overlap, the trailing block is `b`'s, and the off-diagonal block
      !! is what the inter-fragment terms want. `n_ao_a` is returned so the caller
      !! can find it.
      !!
      !! Only the atoms are included. A fragment's expansion points also sit at bond
      !! midpoints, and those carry multipoles but no basis functions, which is why
      !! the count here is `n_atoms` and not `n_points`.
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      type(czt_molecule_t), intent(out) :: mol
      integer, intent(out) :: n_ao_a
      type(error_t), intent(inout) :: error

      type(molecular_basis_type) :: basis_a, basis_b, both
      real(dp), allocatable :: coords(:, :)
      integer, allocatable :: z(:)
      integer :: na, nb, at

      n_ao_a = 0
      call fragment_basis(frag_a, basis_a, error)
      if (error%has_error()) return
      call fragment_basis(frag_b, basis_b, error)
      if (error%has_error()) return

      na = frag_a%n_atoms
      nb = frag_b%n_atoms
      call both%allocate_elements(na + nb)
      do at = 1, na
         both%elements(at) = basis_a%elements(at)
      end do
      do at = 1, nb
         both%elements(na + at) = basis_b%elements(at)
      end do

      allocate (coords(3, na + nb), z(na + nb))
      do at = 1, na
         coords(:, at) = frag_a%points(:, at) + offset_a
         z(at) = nint(frag_a%charge(at))
      end do
      do at = 1, nb
         coords(:, na + at) = frag_b%points(:, at) + offset_b
         z(na + at) = nint(frag_b%charge(at))
      end do

      ! Cartesian explicitly, not left to the basis object to declare. Setting the
      ! angular form per element is not enough: the molecule asks the basis as a
      ! whole, which comes back spherical and gives five d functions where the
      ! potential has six.
      call mol%build(z, coords, both, error, force_cartesian=.true.)
      if (error%has_error()) then
         deallocate (coords, z)
         return
      end if

      ! Where a's functions stop, counted from the basis: the two fragments need not
      ! be the same species.
      n_ao_a = 0
      do at = 1, na
         n_ao_a = n_ao_a + basis_a%elements(at)%num_basis_functions()
      end do

      deallocate (coords, z)
   end subroutine two_fragment_molecule

end module mqc_czt_efp_pair
