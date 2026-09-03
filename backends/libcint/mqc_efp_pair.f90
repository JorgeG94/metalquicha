!! A molecule spanning two fragments, and the overlaps between them
module mqc_efp_pair
   !! Exchange repulsion, charge transfer and the damping on dispersion all come from
   !! the same thing: overlaps between *one fragment's* localized orbitals and
   !! *another's*. Those are integrals over both fragments' basis functions at once,
   !! so they need a molecule that spans the pair -- which is what this builds.
   !!
   !! The basis comes from the potential's own `PROJECTION BASIS SET`, recovered by
   !! `mqc_efp_read` with GAMESS's primitive normalization divided back out. The
   !! *file's* basis, not a lookup of the name it was computed with: a potential is
   !! a self-contained object, and the shipped GAMESS library potentials do not all
   !! name a basis this program has.
   use pic_types, only: dp
   use mqc_physical_constants, only: PI
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cgto, only: molecular_basis_type, ANGULAR_FORM_CARTESIAN
   use mqc_elements, only: element_number_to_symbol
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_efp_read, only: efp_fragment_t, N_QUADRUPOLE
   use mqc_efp_interaction, only: CP_WEIGHT
   use mqc_libcint_esp, only: esp_matrices, drinv_matrices, ddrinv_matrices
   use mqc_efp_potential, only: from_gamess_ao_order, frozen_core
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

   real(dp), parameter :: S_FLOOR = 1.0e-5_dp
      !! Overlap below which a pair's damping series is not evaluated at all,
      !! `efdrvr.src:4464`.

   integer, parameter :: N_DQ_SLOTS = 27
      !! Slots in a `DIPOLE-QUADRUPOLE` record: a 3x3x3 tensor, last index fastest.

   integer, parameter :: N_QQ_SLOTS = 81
      !! Slots in the `LMOQQPOL` record: a 3x3x3x3 tensor, last index fastest.

contains

   subroutine fragment_basis(frag, basis, error)
      !! A basis object for one fragment, from the projection basis it carries
      !!
      !! Cartesian: the ordering maps in `mqc_efp_potential` cover Cartesian s, p,
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
         basis%elements(at)%element = trim(element_number_to_symbol(nint(frag%charge(at))))
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
      !! Charge transfer, from `ECHTR` in `efchtr.src`
      !!
      !!     for i occupied on A, n virtual on B:
      !!       SAB2 = sum_M S_AB(M,n)^2                  M over A's occupied + virtual
      !!       W    = V_AB(i,n) - sum_M S_AB(M,n) V_AA(i,M)
      !!       STIN = sum_J ( T_BB(n,J) - sum_M S_AB(M,n) T_AB(M,J) ) S_AB(i,J)
      !!       CT  += W (W + STIN) / ( (1 - SAB2) (F_A(i) - T_BB(n,n)) )
      !!     E = 2 CT(A->B) + 2 CT(B->A)
      !!
      !! The orbitals are `CTVEC`'s -- occupied *and* virtual, since this moves
      !! density into the latter -- and `F_A(i)` is `CTFOK`. The denominator pairs
      !! an orbital energy with a *kinetic-energy diagonal* standing in for the
      !! virtual's; that is GAMESS's expression, not an approximation introduced
      !! here.
      !!
      !! `V` is the other fragment's multipole potential over these orbitals:
      !! charges, dipoles and quadrupoles, and no octupole. The rank ending at the
      !! quadrupole is GAMESS's answer being unchanged by an octupole section,
      !! not an omission.
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      type(error_t), intent(inout) :: error
      real(dp) :: energy

      real(dp), allocatable :: s_ao(:, :), t_ao(:, :), ct_a(:, :), ct_b(:, :)
      real(dp), allocatable :: v_from_b(:, :), v_from_a(:, :)
      type(libcint_molecule_t) :: pair
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
      !! Density moving from `a`'s occupied orbitals into `b`'s virtuals
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: ct_a(:, :), ct_b(:, :)
      real(dp), intent(in) :: s_ao(:, :), t_ao(:, :), v_b(:, :)
      real(dp) :: ct

      real(dp), allocatable :: sab(:, :), tab(:, :), tbb(:, :), vab(:, :), vaa(:, :)
      real(dp) :: sab2, w, stin, denom, inner
      integer :: n_a, n_b, occ_a, occ_b, i, n, m, j

      n_a = frag_a%n_mo_ct
      n_b = frag_b%n_mo_ct
      occ_a = frag_a%n_occ_ct
      occ_b = frag_b%n_occ_ct
      ct = 0.0_dp

      sab = project(ct_a, s_ao, ct_b)
      tab = project(ct_a, t_ao, ct_b)
      tbb = project(ct_b, t_ao, ct_b)
      vab = project(ct_a, v_b, ct_b)
      vaa = project(ct_a, v_b, ct_a)

      do i = 1, occ_a
         do n = occ_b + 1, n_b
            sab2 = 0.0_dp
            w = vab(i, n)
            do m = 1, n_a
               sab2 = sab2 + sab(m, n)*sab(m, n)
               w = w - sab(m, n)*vaa(i, m)
            end do
            stin = 0.0_dp
            do j = 1, occ_b
               inner = tbb(n, j)
               do m = 1, n_a
                  inner = inner - sab(m, n)*tab(m, j)
               end do
               stin = stin + inner*sab(i, j)
            end do
            denom = (1.0_dp - sab2)*(frag_a%eps_occ(i) - tbb(n, n))
            ct = ct + w*(w + stin)/denom
         end do
      end do
      deallocate (sab, tab, tbb, vab, vaa)
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
      !! `VEFP`, built by `EFCEF`, `EFDEF` and `EFQEF` in `efchtr.src`: charges,
      !! dipoles and quadrupoles, and **no octupole** -- one rank fewer than the
      !! electrostatic energy uses.
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
      !! nine `ab` -- `efchtr.src:1801-1807`'s three diagonal terms plus twice its
      !! three off-diagonal ones, written out.
      !!
      !! The dipole and quadrupole moments here are **electronic only**; the
      !! nucleus sits entirely in the monopole, which is how the potential stores
      !! them.
      type(libcint_molecule_t), intent(in) :: pair
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
      ! sign fixed in `test_mqc_libcint_esp` against a difference of `esp_matrices`:
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

      ! The quadrupole rank, `efchtr.src:1801-1807`. Its three diagonal terms and twice
      ! its three off-diagonal ones are the sum over all nine `ab` of a symmetric `Q`.
      !
      ! **The `1/3` is not in `efchtr.src`, and belongs here anyway.** It is the
      ! coefficient a *traceless* quadrupole carries against `grad grad (1/R)`: the
      ! expansion is `(1/2) M_ab grad_a grad_b` over the second moments, and
      ! substituting `Q = (3M - tr M)/2` turns the half into a third, because
      ! `delta_ab grad_a grad_b (1/R)` vanishes. GAMESS folds it into its Rys
      ! assembly instead of writing it at the contraction.
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
      !! The stored second moments as the electric quadrupole `EFQEF` contracts
      !!
      !! `efchtr.src:1557-1571`, which does this conversion on its own copy before
      !! using it -- so a potential's `QUADRUPOLES` section holds **second moments,
      !! not traceless quadrupoles**, as a nonzero stored trace shows.
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
      type(libcint_molecule_t), intent(in) :: pair
      integer, intent(in) :: offset_ao, n_own
      real(dp), allocatable, intent(out) :: padded(:, :)
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: own
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

      type(libcint_molecule_t) :: pair
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
      !! `E6` with GAMESS's overlap-based damping, which is what it reports
      !!
      !! GAMESS's default `IDISDMP = 1` damps each orbital pair by a Tang-Toennies
      !! series in that pair's own overlap (`efdrvr.src` around line 2734):
      !!
      !!     RB = -2 ln|S_ij|
      !!     F6 = 1 - S_ij^2 ( 1 + sqrt(RB) + RB/2 + RB^(3/2)/6
      !!                         + RB^2/24 + RB^(5/2)/120 + RB^3/720 )
      !!
      !! and `F6 = 1` where `|S_ij| <= 1e-5`, the series being meaningless once the
      !! orbitals do not overlap. **The series is in `sqrt(RB)`, not in `RB`**, so
      !! six of its seven terms would be wrong if read as an ordinary exponential
      !! expansion.
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      type(error_t), intent(inout) :: error
      real(dp) :: energy

      real(dp), allocatable :: s_lmo(:, :), t_lmo(:, :)
      real(dp) :: sep(3)
      real(dp) :: c6, dist, f6, f7, f8
      integer :: i, j, k

      energy = 0.0_dp
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
               c6 = c6 + CP_WEIGHT(k)*isotropic_alpha(frag_a%dyn_pol(:, :, i, k)) &
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

   pure subroutine overlap_damping(sab, f6, f7, f8)
      !! `Damping_for_Dispersion`'s overlap branch, `efdrvr.src:4462-4517`
      !!
      !! One routine produces all three damping factors, and they are one series
      !! truncated at three different orders:
      !!
      !!     F_n = 1 - S^2 sum_{k=0..K} x^k / k!,   x = sqrt(-2 ln|S|)
      !!
      !! with `K = 6, 7, 8`. **The variable is `sqrt(RB)`, not `RB`**, so six of
      !! F6's seven terms carry half-integer powers -- read as an ordinary
      !! exponential expansion the series still damps, plausibly and wrongly.
      !!
      !! `F8` is assembled as the source assembles it, by subtracting the two extra
      !! terms from `F6` rather than by re-summing; algebraically it is the `K = 8`
      !! sum. Each factor subtracts one more positive term than the last, so
      !! `F8 < F7 < F6 < 1` wherever the series is used at all.
      !!
      !! Below `|S| = 1e-5` all three are 1.
      real(dp), intent(in) :: sab
      real(dp), intent(out) :: f6, f7, f8

      real(dp) :: rb, s2

      f6 = 1.0_dp
      f7 = 1.0_dp
      f8 = 1.0_dp
      if (abs(sab) <= S_FLOOR) return

      rb = -2.0_dp*log(abs(sab))
      s2 = sab*sab
      f6 = 1.0_dp - s2*(1.0_dp + sqrt(rb) + rb/2.0_dp &
                        + rb**1.5_dp/6.0_dp + rb*rb/24.0_dp &
                        + rb**2.5_dp/120.0_dp + rb**3/720.0_dp)
      f7 = f6 - s2*rb**3.5_dp/5040.0_dp
      f8 = f7 - s2*rb**4/40320.0_dp
   end subroutine overlap_damping

   function dispersion_e8_damped(frag_a, frag_b, offset_a, offset_b, error) result(energy)
      !! `E8`, the isotropic one -- which is the one GAMESS prints
      !!
      !! `Disp8_LMOpol` (`efdrvr.src:4311`) computes two unrelated things, and the
      !! anisotropic one is never printed. What reaches `E8 DISPERSION ENERGY` is
      !! the isotropic form at `efdrvr.src:4401-4418`, computed whenever the
      !! potential carries an `LMOQQPOL` section at all:
      !!
      !!     C8 = sum_f (15/pi) FACT(f) ( a_iso^A A_QQ^B + a_iso^B A_QQ^A )
      !!     E8 = - F8 C8 / R^8
      !!
      !! It therefore needs the trace of the dipole-dipole polarizability and the
      !! spherical average of the quadrupole-quadrupole one, and nothing else. Both
      !! averages are isotropic contractions and survive the rotation into the
      !! current frame unchanged, which is why nothing is rotated here.
      !!
      !! **Not `E6/3`.** That approximation is `efdrvr.src:1917` and does not run
      !! for a pair of file-based fragments; reconciling against it would fit a
      !! factor to the wrong quantity.
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      type(error_t), intent(inout) :: error
      real(dp) :: energy

      real(dp), allocatable :: s_lmo(:, :), t_lmo(:, :)
      real(dp) :: sep(3)
      real(dp) :: c8, dist, f6, f7, f8
      integer :: i, j, k

      energy = 0.0_dp
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
               c8 = c8 + CP_WEIGHT(k) &
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
      !! `E7`, the dipole-dipole/dipole-quadrupole cross term
      !!
      !! `Disp7_LMOpol`, `efdrvr.src:3979-4308`, accumulating `efdrvr.src:4042-4049`:
      !!
      !!     DUM1 = DD_A(a,c) DQ_B(b,d,e)
      !!     DUM2 = DD_B(b,e) DQ_A(a,c,d)
      !!     E7 = F7 sum_f sum_abcde (-1/3pi) T2(a,b) T3(c,d,e) FACT(f) (DUM1 - DUM2)
      !!
      !! `a` and `b` are the two `T2` slots and `c,d,e` the three `T3` slots; `a`
      !! belongs to A and `b` to B. Each dipole-quadrupole tensor puts its dipole
      !! index on a `T2` slot and its quadrupole pair on `T3` slots. There is only
      !! one E7 -- no isotropic variant exists.
      !!
      !! **Three conventions here are sign- or transpose-critical**, and E7 is the
      !! first term that can see any of them; E6 and E8 reach the polarizabilities
      !! only through isotropic averages and the separation only through `R`.
      !!
      !! *The displacement runs A minus B* (`efdrvr.src:1724-1726`), A being this
      !! routine's first fragment. E7 is odd in `C`: `T2` is even and `T3` is odd,
      !! so the whole term changes sign with it.
      !!
      !! *`T3` carries a deliberate extra negative* (`efdrvr.src:3507`). The
      !! textbook `-grad grad grad 1/R` gives E7 the wrong sign; `t_tensors` builds
      !! the form GAMESS's routine actually hands over.
      !!
      !! *The dipole-dipole tensor arrives transposed.* GAMESS indexes it
      !! `(field, dipole)` and `mqc_efp_read` indexes it `(dipole, field)`, so
      !! GAMESS's `DYNDD_LMO_ROT(a,c)` is our `dyn_pol(c,a)`. The swap below is
      !! load-bearing: the tensor's antisymmetric part is 12% of it here.
      !!
      !! The dipole-quadrupole tensor needs no such care -- `mqc_efp_read` keeps it
      !! flat in file order, so GAMESS's own slot formula recovers GAMESS's tensor.
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      type(error_t), intent(inout) :: error
      real(dp) :: energy

      real(dp), allocatable :: s_lmo(:, :), t_lmo(:, :)
      real(dp) :: sep(3), t2(3, 3), t3(3, 3, 3)
      real(dp) :: pair, weighted, f6, f7, f8
      integer :: i, j, k, ia, ib, ic, id, ie

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

      do i = 1, frag_a%n_lmo
         do j = 1, frag_b%n_lmo
            ! A minus B, A being this routine's first fragment -- the role GAMESS
            ! gives the lower fragment index.
            sep = frag_a%centroids(:, i) + offset_a - frag_b%centroids(:, j) - offset_b
            call t_tensors(sep, t2, t3)

            pair = 0.0_dp
            do k = 1, min(frag_a%n_freq, frag_b%n_freq)
               weighted = 0.0_dp
               do ia = 1, 3
                  do ib = 1, 3
                     do ic = 1, 3
                        do id = 1, 3
                           do ie = 1, 3
                              weighted = weighted + t2(ia, ib)*t3(ic, id, ie) &
                                         *(frag_a%dyn_pol(ic, ia, i, k) &
                                           *frag_b%dipquad(dq_slot(ib, id, ie), j, k) &
                                           - frag_b%dyn_pol(ie, ib, j, k) &
                                           *frag_a%dipquad(dq_slot(ia, ic, id), i, k))
                           end do
                        end do
                     end do
                  end do
               end do
               pair = pair + CP_WEIGHT(k)*weighted
            end do
            pair = -pair/(3.0_dp*PI)

            call overlap_damping(s_lmo(i, j), f6, f7, f8)
            energy = energy + f7*pair
         end do
      end do
      deallocate (s_lmo, t_lmo)
   end function dispersion_e7_damped

   pure subroutine t_tensors(c, t2, t3)
      !! The rank-2 and rank-3 interaction tensors, `efdrvr.src:3449-3562`
      !!
      !!     T2(i,j)   = ( 3 C_i C_j - R^2 d_ij ) / R^5
      !!     T3(i,j,k) = ( 15 C_i C_j C_k
      !!                   - 3 R^2 ( C_i d_jk + C_j d_ik + C_k d_ij ) ) / R^7
      !!
      !! `T3`'s sign is GAMESS's, not the textbook's: `T_tensor_3` builds the
      !! negative of the form above and flips it in place at `efdrvr.src:3507`.
      !! What is written here is the net tensor its caller receives, and E7 is
      !! linear in `T3`.
      real(dp), intent(in) :: c(3)
      real(dp), intent(out) :: t2(3, 3), t3(3, 3, 3)

      real(dp) :: r2, r5, r7
      integer :: i, j, k

      r2 = c(1)*c(1) + c(2)*c(2) + c(3)*c(3)
      r5 = r2*r2*sqrt(r2)
      r7 = r5*r2

      do j = 1, 3
         do i = 1, 3
            t2(i, j) = (3.0_dp*c(i)*c(j) - r2*delta(i, j))/r5
         end do
      end do

      do k = 1, 3
         do j = 1, 3
            do i = 1, 3
               t3(i, j, k) = (15.0_dp*c(i)*c(j)*c(k) &
                              - 3.0_dp*r2*(c(i)*delta(j, k) + c(j)*delta(i, k) &
                                           + c(k)*delta(i, j)))/r7
            end do
         end do
      end do
   end subroutine t_tensors

   pure function delta(i, j) result(d)
      !! Kronecker delta, so the tensors above read like the source they came from
      integer, intent(in) :: i, j
      real(dp) :: d

      d = 0.0_dp
      if (i == j) d = 1.0_dp
   end function delta

   pure function dq_slot(i, j, k) result(slot)
      !! Where `DQ(i,j,k)` sits in a `DIPOLE-QUADRUPOLE` record
      !!
      !! `i` is the dipole index and `(j,k)` the quadrupole pair. Row-major with the
      !! last index fastest, in GAMESS's writer (`efinp.src:7635`) and reader
      !! (`efinp.src:12943-12949`) alike.
      integer, intent(in) :: i, j, k
      integer :: slot

      slot = (i - 1)*9 + (j - 1)*3 + k
   end function dq_slot

   pure function isotropic_quadquad(values) result(a)
      !! `DYNQQ_LMO_AVE`, the spherical average of a quadrupole-quadrupole tensor
      !!
      !! `efdrvr.src:1567-1572` contracts the full rank-four tensor against the
      !! isotropic projector
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
      !! From `EXREP` in GAMESS's `efpaul.src` -- the energy routine, not `ODM` in
      !! the same file, which is its gradient and carries factors that only mean
      !! something in a derivative. Three terms over pairs of localized orbitals on
      !! different fragments, with `S` and `T` in the localized-orbital basis:
      !!
      !!     XR1 = -2 sqrt(2/pi) sum_ij sqrt(-ln|S_ij|) S_ij^2 / R_ij
      !!     XR2 = -sum_ij S_ij ( sum_k F^A_ik S_kj + sum_l F^B_jl S_il - 2 T_ij )
      !!     XR3 = +sum_ij S_ij^2 ( V_i + V_j - 1/R_ij )
      !!     E   = 2 (XR1 + XR2 + XR3)
      !!
      !! `V_i` is the electrostatic potential at orbital `i`'s centroid from the *other*
      !! fragment: minus its nuclear charges over their distances, plus two over the
      !! distance to each of its orbital centroids. The factor of two outside, and the
      !! two inside `V`, are `ICOEFF` and `JCOEFF` at `MLSWTCH = 1`, which is the
      !! closed-shell case.
      !!
      !! The first term is skipped where `|S_ij| <= 1e-7`: its logarithm diverges as
      !! the overlap vanishes, and the `S^2` in front does not save it in floating
      !! point. GAMESS's separate screening of a whole fragment pair below `1e-6` is
      !! not copied here.
      type(efp_fragment_t), intent(in) :: frag_a, frag_b
      real(dp), intent(in) :: offset_a(3), offset_b(3)
      type(error_t), intent(inout) :: error
      real(dp) :: energy

      real(dp), parameter :: RT2PI = 0.7978845608028654_dp   !! sqrt(2/pi)
      real(dp), parameter :: S_FLOOR_LOCAL = 1.0e-7_dp
         !! Named apart from the module's `S_FLOOR` deliberately: that one is the
         !! 1e-5 damping cutoff from `efdrvr.src:4464`, this is a different and
         !! tighter threshold, and one name for both would make the module constant
         !! silently mean something else inside this procedure.
      type(libcint_molecule_t) :: pair
      real(dp), allocatable :: s_ao(:, :), t_ao(:, :), lmo_a(:, :), lmo_b(:, :)
      real(dp), allocatable :: s_lmo(:, :), t_lmo(:, :), work(:, :)
      real(dp), allocatable :: cen_a(:, :), cen_b(:, :), v_a(:), v_b(:)
      integer :: n_ao_a, na, nb, i, j, n_lmo_a, n_lmo_b
      real(dp) :: xr1, xr2, xr3, sij, sij2, tij, rij, fij, fji
      integer :: k, l

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

      na = frag_a%n_atoms
      nb = frag_b%n_atoms
      allocate (cen_a(3, n_lmo_a), cen_b(3, n_lmo_b))
      do i = 1, n_lmo_a
         cen_a(:, i) = frag_a%pol_points(:, i) + offset_a
      end do
      do j = 1, n_lmo_b
         cen_b(:, j) = frag_b%pol_points(:, j) + offset_b
      end do

      ! The potential at each centroid from the other fragment, built once.
      allocate (v_a(n_lmo_a), v_b(n_lmo_b))
      do i = 1, n_lmo_a
         v_a(i) = other_potential(cen_a(:, i), frag_b, offset_b, cen_b)
      end do
      do j = 1, n_lmo_b
         v_b(j) = other_potential(cen_b(:, j), frag_a, offset_a, cen_a)
      end do

      xr1 = 0.0_dp
      xr2 = 0.0_dp
      xr3 = 0.0_dp
      do j = 1, n_lmo_b
         do i = 1, n_lmo_a
            sij = s_lmo(i, j)
            sij2 = sij*sij
            tij = t_lmo(i, j)
            rij = sqrt(sum((cen_a(:, i) - cen_b(:, j))**2))

            if (abs(sij) > S_FLOOR_LOCAL) then
               xr1 = xr1 - 2.0_dp*RT2PI*sqrt(-log(abs(sij)))*sij2/rij
            end if

            fij = 0.0_dp
            do k = 1, n_lmo_a
               fij = fij + frag_a%fock_lmo(i, k)*s_lmo(k, j)
            end do
            fji = 0.0_dp
            do l = 1, n_lmo_b
               fji = fji + frag_b%fock_lmo(j, l)*s_lmo(i, l)
            end do
            xr2 = xr2 - sij*(fij + fji - 2.0_dp*tij)

            xr3 = xr3 + sij2*(v_a(i) + v_b(j) - 1.0_dp/rij)
         end do
      end do
      energy = 2.0_dp*(xr1 + xr2 + xr3)

      call pair%destroy()
      deallocate (s_ao, t_ao, lmo_a, lmo_b, s_lmo, t_lmo, work, cen_a, cen_b, v_a, v_b)
   end function exchange_repulsion

   function other_potential(point, frag, offset, centroids) result(v)
      !! The potential at `point` from `frag`: its nuclei, then its orbital centroids
      real(dp), intent(in) :: point(3)
      type(efp_fragment_t), intent(in) :: frag
      real(dp), intent(in) :: offset(3)
      real(dp), intent(in) :: centroids(:, :)
      real(dp) :: v

      integer :: k

      v = 0.0_dp
      do k = 1, frag%n_atoms
         ! Valence charge, not the full nuclear one: the localized orbitals a
         ! potential carries are valence only, so the core electrons count as
         ! screening their own nucleus. The same number the PROJECTION BASIS SET
         ! atom header carries.
         v = v - (frag%charge(k) - 2.0_dp*real(frozen_core([nint(frag%charge(k))]), dp)) &
             /sqrt(sum((point - frag%points(:, k) - offset)**2))
      end do
      do k = 1, size(centroids, 2)
         v = v + 2.0_dp/sqrt(sum((point - centroids(:, k))**2))
      end do
   end function other_potential

   subroutine padded_lmo(frag, pair, offset_ao, n_own, padded, error)
      !! One fragment's orbitals in the pair's full AO space, zero elsewhere
      type(efp_fragment_t), intent(in) :: frag
      type(libcint_molecule_t), intent(in) :: pair
      integer, intent(in) :: offset_ao, n_own
      real(dp), allocatable, intent(out) :: padded(:, :)
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: own
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
      type(libcint_molecule_t), intent(out) :: mol
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
      type(libcint_molecule_t), intent(in) :: mol
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
      type(libcint_molecule_t), intent(out) :: mol
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

end module mqc_efp_pair
