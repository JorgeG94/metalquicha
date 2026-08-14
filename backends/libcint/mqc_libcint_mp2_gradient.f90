!! The MP2 nuclear gradient over a closed-shell reference
module mqc_libcint_mp2_gradient
   !! **Why this is not the shape of every gradient before it.** Hartree-Fock,
   !! density fitting and Kohn-Sham are all variational: the energy is
   !! stationary with respect to the orbitals, so moving a nucleus changes the
   !! energy only through the integrals, and the orbital response never appears.
   !! MP2 is not stationary in the reference orbitals. Its derivative therefore
   !! carries a term for how the orbitals themselves relax, and that term is the
   !! solution of a linear system the size of the occupied-virtual space -- the
   !! Z-vector equation, solved once rather than once per displacement, which is
   !! the whole reason an analytic MP2 gradient is affordable at all.
   !!
   !! The assembly follows Handy and Schaefer's interchange, in the arrangement
   !! `pyscf.grad.mp2` uses, because that is what these numbers are checked
   !! against and a different but equivalent grouping would make a disagreement
   !! impossible to localise.
   !!
   !! **What is built, in order.**
   !!
   !! 1. `t2(i,j,a,b) = (ia|jb) / (e_i + e_j - e_a - e_b)`, the amplitudes.
   !! 2. `doo` and `dvv`, the occupied-occupied and virtual-virtual blocks of
   !!    the unrelaxed one-particle correction.
   !! 3. `gamma`, the two-particle density in the AO basis. This is the
   !!    expensive object -- `n_ao^4` -- and the reason this routine is for
   !!    validation-sized systems until the blocked version replaces it.
   !! 4. `Imat`, that two-particle density contracted with the *undifferentiated*
   !!    integrals, which is what makes the Lagrangian.
   !! 5. The Z-vector solve, through the response operator the CPHF module
   !!    already applies -- it is the same operator, and writing a second one
   !!    would mean two things to keep in agreement.
   !! 6. The contraction with the differentiated integrals.
   !!
   !! **Restricted, exact ERIs, no frozen core.** Each of those is a refusal
   !! rather than an approximation: an unrestricted reference needs separate
   !! alpha and beta amplitudes, a fitted one differentiates a different energy,
   !! and a frozen core adds occupied-frozen and virtual-frozen blocks to the
   !! relaxed density which are not built here.
   use pic_types, only: dp, int64
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim, atom_ao_blocks, &
                                    two_electron_block
   use mqc_libcint_mp2, only: transform_ovov
   use mqc_libcint_cphf, only: cphf_solve
   use mqc_libcint_rhf, only: build_fock
   use mqc_libcint_direct, only: build_fock_direct, schwarz_bounds, direct_stats_t
   use pic_blas_interfaces, only: pic_gemm
   use mqc_libcint_gradient, only: nuclear_repulsion_gradient, one_electron_deriv, &
                                   iprinv_deriv_at, DERIV_OVLP, DERIV_KIN, DERIV_NUC
   use libcint_fortran, only: libcint_2e_ip1_sph, libcint_2e_ip1_cart, &
                              libcint_2e_ip1_sph_optimizer, libcint_2e_ip1_cart_optimizer, &
                              libcint_del_optimizer
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
   implicit none
   private

   public :: libcint_mp2_gradient

   real(dp), parameter :: BLOCK_TARGET = 5.0e8_dp
      !! Bytes the blocked path aims to hold in its two live block arrays. Not a
      !! hard limit on the routine -- the half-transformed amplitudes and the
      !! one-particle quantities sit outside it -- but it is what sets the block
      !! size, and it is the number to raise on a machine with room.

   real(dp), parameter :: IN_CORE_LIMIT = 4.0e9_dp
      !! Bytes. The same ceiling `mqc_libcint_cphf` applies to its own stored
      !! tensor, for the same reason and deliberately not a different number.

contains

   subroutine libcint_mp2_gradient(mol, coeff, orbital_energies, n_occ, gradient, &
                                   error, n_frozen, block_bytes, force_blocked)
      !! dE(MP2)/dR for a closed-shell reference, in Hartree/Bohr
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)            !! C, (n_ao, n_mo)
      real(dp), intent(in) :: orbital_energies(:)    !! (n_mo), Hartree
      integer, intent(in) :: n_occ                   !! Doubly occupied count
      real(dp), allocatable, intent(out) :: gradient(:, :)   !! (3, natm)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_frozen
      real(dp), intent(in), optional :: block_bytes
         !! Byte target for the block sizes, in place of `BLOCK_TARGET`. Small
         !! enough, it forces several blocks of both the first index and the
         !! occupied one, which is what the check harness passes: a single
         !! block exercises neither loop nor the offsets they carry.
      logical, intent(in), optional :: force_blocked
         !! Take the blocked path whatever the size. Separate from the target
         !! on purpose -- the dense path blocks over the occupied index too,
         !! and a knob that forced the other path would leave that untested.

      real(dp), allocatable :: eri_packed(:, :), eri(:, :, :, :)
      real(dp), allocatable :: ovov(:, :, :, :), t2(:, :, :, :)
      real(dp), allocatable :: doo(:, :), dvv(:, :), dm1mo(:, :), zeta(:, :)
      real(dp), allocatable, target :: gamma_ao(:, :, :, :)
      real(dp), allocatable :: imat_ao(:, :), imat(:, :), im1(:, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :)
      real(dp), allocatable :: hf_density(:, :), dm1(:, :), dm1_total(:, :), dm1p(:, :)
      real(dp), allocatable :: veff(:, :), xvo(:, :, :), zvec(:, :, :)
      real(dp), allocatable :: overlap(:, :), work(:, :), p_occ(:, :), vhf_s1occ(:, :)
      real(dp), allocatable :: zero_h(:, :)
      real(dp), allocatable :: s1(:, :, :), h1(:, :, :), kin(:, :, :)
      real(dp), allocatable :: vrinv(:, :, :), hcore_a(:, :, :)
      real(dp), allocatable :: de2(:, :), vhf1(:, :, :, :)
      real(dp), allocatable :: im1_t(:, :), work_t(:, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: n_ao, n_mo, n_o, n_v, frozen, iatom, comp, p0, p1, p, q, i, a
      logical :: dense
      real(dp) :: target_bytes
      real(dp), allocatable :: bounds(:, :)

      n_ao = mol%nao
      n_mo = size(coeff, 2)
      n_o = n_occ
      n_v = n_mo - n_occ

      frozen = 0
      if (present(n_frozen)) frozen = n_frozen
      if (frozen /= 0) then
         call error%set(ERROR_VALIDATION, "the MP2 gradient does not implement a "// &
                        "frozen core: the relaxed density gains occupied-frozen and "// &
                        "virtual-frozen blocks, which are not built here. Run the "// &
                        "gradient with freeze_core off.")
         return
      end if
      if (n_v < 1 .or. n_o < 1) then
         call error%set(ERROR_VALIDATION, "the MP2 gradient needs both occupied and "// &
                        "virtual orbitals")
         return
      end if
      if (size(orbital_energies) /= n_mo) then
         call error%set(ERROR_VALIDATION, "the MP2 gradient: one orbital energy per "// &
                        "orbital")
         return
      end if

      ! Dense or blocked. Dense keeps two `n_ao^4` arrays live -- the integrals
      ! and the two-particle density -- and is worth taking whenever it fits,
      ! because one gemm over the whole first index beats the same gemm run
      ! block by block, and because the integrals get built once with their full
      ! eightfold symmetry instead of once per block with only the ket pair's.
      ! Blocked pays roughly four times the integral work for a memory ceiling
      ! that is a choice rather than the basis size.
      dense = 2.0_dp*real(n_ao, dp)**4*8.0_dp <= IN_CORE_LIMIT
      target_bytes = BLOCK_TARGET
      if (present(block_bytes)) target_bytes = block_bytes
      if (present(force_blocked)) then
         if (force_blocked) dense = .false.
      end if

      allocate (c_occ(n_ao, n_o), c_vir(n_ao, n_v))
      c_occ = coeff(:, 1:n_o)
      c_vir = coeff(:, n_o + 1:n_mo)

      ! ---- amplitudes -------------------------------------------------------
      call mol%eris_packed(eri_packed)
      call transform_ovov(eri_packed, coeff, frozen, n_occ, n_ao, n_mo, ovov)
      deallocate (eri_packed)

      call build_amplitudes(ovov, orbital_energies, n_o, n_v, n_occ, t2)

      ! ---- the unrelaxed one-particle correction ----------------------------
      call gamma1_intermediates(t2, n_o, n_v, doo, dvv)

      allocate (dm1mo(n_mo, n_mo))
      dm1mo = 0.0_dp
      dm1mo(1:n_o, 1:n_o) = doo + transpose(doo)
      dm1mo(n_o + 1:n_mo, n_o + 1:n_mo) = dvv + transpose(dvv)

      ! ---- the two-particle density, and what it contracts against ----------
      if (.not. dense) then
         ! The direct Fock builds below need these, and so does the response
         ! solve, which cannot be handed a tensor that was never formed.
         call schwarz_bounds(mol, bounds, error)
         if (error%has_error()) return
      else
         allocate (bounds(0, 0))
      end if

      allocate (imat_ao(n_ao, n_ao), de2(3, mol%natm), vhf1(n_ao, n_ao, 3, mol%natm))
      imat_ao = 0.0_dp
      de2 = 0.0_dp
      vhf1 = 0.0_dp

      allocate (hf_density(n_ao, n_ao))
      hf_density = 2.0_dp*matmul(c_occ, transpose(c_occ))

      if (dense) then
         call build_two_particle_density(t2, c_occ, c_vir, n_ao, n_o, n_v, &
                                         target_bytes, gamma_ao)
         call mol%eris(eri)
         call contract_gamma_eri(eri, gamma_ao, n_ao, imat_ao)
         call two_electron_mp2_terms(mol, gamma_ao, hf_density, de2, vhf1)
         deallocate (gamma_ao)
      else
         call blocked_two_electron_terms(mol, t2, c_occ, c_vir, n_ao, n_o, n_v, &
                                         target_bytes, hf_density, imat_ao, de2, &
                                         vhf1, error)
         if (error%has_error()) return
      end if
      deallocate (t2)

      ! ---- the Lagrangian, and the orbitals relaxing under it ---------------
      !
      ! `imat` in the MO basis is where the two-particle density enters the
      ! right-hand side; the `veff` term is where the unrelaxed one-particle
      ! correction does. Neither alone is the Lagrangian.
      call mol%overlap(overlap)
      allocate (imat(n_mo, n_mo), work(n_ao, n_ao))
      work = matmul(imat_ao, matmul(overlap, coeff(:, 1:n_mo)))
      imat = -matmul(transpose(coeff), work)
      deallocate (work)

      allocate (dm1(n_ao, n_ao))
      dm1 = matmul(coeff, matmul(dm1mo, transpose(coeff)))
      ! `build_fock` returns `H + J - K/2`; a zero core Hamiltonian leaves the
      ! two-electron part alone, which is what a response wants. Threaded over
      ! the target element, and already the hot routine of the CPHF solver.
      allocate (zero_h(n_ao, n_ao), veff(n_ao, n_ao))
      zero_h = 0.0_dp
      call two_electron_potential(dense, mol, eri, bounds, zero_h, dm1, veff, error)
      if (error%has_error()) return
      veff = 2.0_dp*veff

      allocate (xvo(n_v, n_o, 1))
      xvo(:, :, 1) = matmul(transpose(c_vir), matmul(veff, c_occ))
      do i = 1, n_o
         do a = 1, n_v
            xvo(a, i, 1) = xvo(a, i, 1) + imat(i, n_o + a) - imat(n_o + a, i)
         end do
      end do

      ! Tighter than the solver's default. The Z-vector residual enters the
      ! gradient directly rather than quadratically -- there is no variational
      ! principle to make it second order here -- so the default 1e-9 leaves
      ! ~1e-7 in the answer, which is the size of the terms being checked.
      ! In core, not screened, and over the tensor already in hand. The
      ! solver's default is integral-direct with a Schwarz bound, which is right
      ! where the response is one quantity among many; here the screening error
      ! would land in the gradient undiluted. Handing over the tensor built for
      ! the Lagrangian rather than letting it build its own saves a second
      ! identical pass over every integral -- 10% of this routine at cc-pVTZ.
      if (dense) then
         call cphf_solve(mol, coeff, orbital_energies, n_occ, response=zvec, &
                         error=error, mo_rhs=xvo, tol=1.0e-12_dp, max_iter=200, &
                         eri_in=eri)
      else
         ! Integral-direct, with the same Schwarz bound the Fock builds use.
         ! Measured against the stored path on water/cc-pVDZ, the screening
         ! moves the gradient by less than 1e-12 -- it is the reference's own
         ! Krylov solver, not screening, that put 4e-8 between the two codes.
         call cphf_solve(mol, coeff, orbital_energies, n_occ, response=zvec, &
                         error=error, mo_rhs=xvo, tol=1.0e-12_dp, max_iter=200)
      end if
      if (error%has_error()) return

      dm1mo(n_o + 1:n_mo, 1:n_o) = dm1mo(n_o + 1:n_mo, 1:n_o) + zvec(:, :, 1)
      dm1mo(1:n_o, n_o + 1:n_mo) = dm1mo(1:n_o, n_o + 1:n_mo) + transpose(zvec(:, :, 1))

      ! `imat` is not symmetric, and the occupied-virtual block is the one the
      ! Lagrangian defined; the virtual-occupied block is its transpose by
      ! construction rather than by computation.
      imat(n_o + 1:n_mo, 1:n_o) = transpose(imat(1:n_o, n_o + 1:n_mo))
      allocate (im1(n_ao, n_ao))
      im1 = matmul(coeff, matmul(imat, transpose(coeff)))

      ! ---- the energy-weighted densities ------------------------------------
      allocate (zeta(n_mo, n_mo))
      do q = 1, n_mo
         do p = 1, n_mo
            zeta(p, q) = 0.5_dp*(orbital_energies(p) + orbital_energies(q))
         end do
      end do
      ! An occupied orbital energy rather than an average wherever one index is
      ! occupied and the other virtual: that block multiplies the relaxation,
      ! which is a rotation of occupied orbitals into virtual ones.
      do i = 1, n_o
         do a = 1, n_v
            zeta(n_o + a, i) = orbital_energies(i)
            zeta(i, n_o + a) = orbital_energies(i)
         end do
      end do
      zeta = zeta*dm1mo
      allocate (work(n_ao, n_ao))
      work = matmul(coeff, matmul(zeta, transpose(coeff)))
      deallocate (zeta)

      dm1 = matmul(coeff, matmul(dm1mo, transpose(coeff)))

      allocate (p_occ(n_ao, n_ao))
      p_occ = matmul(c_occ, transpose(c_occ))
      call two_electron_potential(dense, mol, eri, bounds, zero_h, &
                                  dm1 + transpose(dm1), veff, error)
      if (error%has_error()) return
      allocate (vhf_s1occ(n_ao, n_ao))
      vhf_s1occ = matmul(p_occ, matmul(veff, p_occ))

      allocate (dm1p(n_ao, n_ao), dm1_total(n_ao, n_ao))
      dm1p = hf_density + 2.0_dp*dm1
      dm1_total = hf_density + dm1

      ! The Hartree-Fock energy-weighted density, on top of the correlation one
      ! already in `work`.
      do i = 1, n_o
         do q = 1, n_ao
            do p = 1, n_ao
               work(p, q) = work(p, q) + 2.0_dp*orbital_energies(i)*c_occ(p, i)*c_occ(q, i)
            end do
         end do
      end do

      ! ---- the one-electron derivatives -------------------------------------
      call one_electron_deriv(mol, s1, DERIV_OVLP)
      s1 = -s1
      call one_electron_deriv(mol, kin, DERIV_KIN)
      call one_electron_deriv(mol, h1, DERIV_NUC)
      h1 = -(kin + h1)
      deallocate (kin)

      allocate (gradient(3, mol%natm))
      ! Zeroed first: `nuclear_repulsion_gradient` accumulates.
      gradient = 0.0_dp
      call nuclear_repulsion_gradient(mol, gradient)
      gradient = gradient + de2

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)
      allocate (vrinv(n_ao, n_ao, 3), hcore_a(n_ao, n_ao, 3))
      allocate (im1_t(n_ao, n_ao), work_t(n_ao, n_ao))
      im1_t = transpose(im1)
      work_t = transpose(work)

      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)

         call iprinv_deriv_at(mol, iatom, vrinv)
         vrinv = -mol%charges(iatom)*vrinv
         if (counts(iatom) > 0) then
            vrinv(p0:p1, :, :) = vrinv(p0:p1, :, :) + h1(p0:p1, :, :)
         end if
         do comp = 1, 3
            hcore_a(:, :, comp) = vrinv(:, :, comp) + transpose(vrinv(:, :, comp))
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + sum(hcore_a(:, :, comp)*dm1_total)
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    - sum(vhf1(:, :, comp, iatom)*dm1p)
         end do

         if (counts(iatom) == 0) cycle

         do comp = 1, 3
            ! Both halves differentiate a basis function on this atom -- the
            ! bra in one, the ket in the other -- so both slice `s1` the same
            ! way and it is the *matrix* that is transposed. Slicing `s1` on
            ! the ket instead makes the two cancel, because the derivative of
            ! an overlap changes sign with which centre moves.
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + sum(s1(p0:p1, :, comp)*im1(p0:p1, :)) &
                                    + sum(s1(p0:p1, :, comp)*im1_t(p0:p1, :))
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    - sum(s1(p0:p1, :, comp)*work(p0:p1, :)) &
                                    - sum(s1(p0:p1, :, comp)*work_t(p0:p1, :))
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    - 2.0_dp*sum(s1(p0:p1, :, comp)*vhf_s1occ(p0:p1, :))
         end do
      end do

   end subroutine libcint_mp2_gradient

   subroutine build_amplitudes(ovov, orbital_energies, n_o, n_v, n_occ, t2)
      !! `t2(i,j,a,b) = (ia|jb) / (e_i + e_j - e_a - e_b)`
      real(dp), intent(in) :: ovov(:, :, :, :)       !! (i, a, j, b)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_o, n_v, n_occ
      real(dp), allocatable, intent(out) :: t2(:, :, :, :)   !! (i, j, a, b)

      integer :: i, j, a, b
      real(dp) :: denom

      allocate (t2(n_o, n_o, n_v, n_v))
      do b = 1, n_v
         do a = 1, n_v
            do j = 1, n_o
               do i = 1, n_o
                  denom = orbital_energies(i) + orbital_energies(j) &
                          - orbital_energies(n_occ + a) - orbital_energies(n_occ + b)
                  t2(i, j, a, b) = ovov(i, a, j, b)/denom
               end do
            end do
         end do
      end do
   end subroutine build_amplitudes

   subroutine gamma1_intermediates(t2, n_o, n_v, doo, dvv)
      !! The occupied-occupied and virtual-virtual blocks of the correction
      !!
      !! Both are the amplitudes contracted with themselves, once directly and
      !! once with the two virtual indices exchanged -- the same two spin cases
      !! the energy expression separates.
      real(dp), intent(in) :: t2(:, :, :, :)
      integer, intent(in) :: n_o, n_v
      real(dp), allocatable, intent(out) :: doo(:, :), dvv(:, :)

      integer :: i, j, k, a, b, c
      real(dp) :: acc

      allocate (doo(n_o, n_o), dvv(n_v, n_v))
      doo = 0.0_dp
      dvv = 0.0_dp

      do j = 1, n_o
         do i = 1, n_o
            acc = 0.0_dp
            do k = 1, n_o
               do b = 1, n_v
                  do a = 1, n_v
                     acc = acc + 2.0_dp*t2(k, i, a, b)*t2(k, j, a, b) &
                           - t2(k, i, a, b)*t2(k, j, b, a)
                  end do
               end do
            end do
            doo(i, j) = -acc
         end do
      end do

      do b = 1, n_v
         do a = 1, n_v
            acc = 0.0_dp
            do i = 1, n_o
               do j = 1, n_o
                  do c = 1, n_v
                     acc = acc + 2.0_dp*t2(i, j, c, a)*t2(i, j, c, b) &
                           - t2(i, j, c, a)*t2(i, j, b, c)
                  end do
               end do
            end do
            dvv(b, a) = acc
         end do
      end do
   end subroutine gamma1_intermediates

   subroutine two_electron_potential(dense, mol, eri, bounds, zero_h, density, &
                                     veff, error)
      !! `J - K/2`, from the stored tensor or from the integrals directly
      !!
      !! Which one is not a choice about accuracy -- the two agree to better
      !! than 1e-12 here -- but about whether a tensor exists to read. The
      !! blocked path never forms one.
      logical, intent(in) :: dense
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: eri(:, :, :, :), bounds(:, :)
      real(dp), intent(in) :: zero_h(:, :), density(:, :)
      real(dp), intent(out) :: veff(:, :)
      type(error_t), intent(inout) :: error

      type(direct_stats_t) :: stats

      if (dense) then
         call build_fock(zero_h, eri, density, veff)
      else
         call build_fock_direct(mol, zero_h, density, bounds, veff, stats, error)
      end if
   end subroutine two_electron_potential

   subroutine build_half_block(t2, c_vir, n_ao, n_o, n_v, j_lo, j_hi, half)
      !! `half(i, (q, r, j))` for one range of `j`: the amplitudes with both
      !! virtual indices in the AO basis and the energy's weighting applied
      !!
      !! `4 t2(i,j,a,b) - 2 t2(i,j,b,a)` is the coefficient the energy gives
      !! each integral.
      !!
      !! Whole, this is `n_o^2 n_ao^2` -- smaller than `n_ao^4` by
      !! `(n_o/n_ao)^2`, and still unbounded, which past a few hundred
      !! functions made it the footprint rather than the block target. `j` is
      !! the index to cut on because it is contracted last, against `C_occ`, so
      !! a range of it contributes an addend to the two-particle density rather
      !! than a slice of one. When the whole range fits, the caller asks for it
      !! in one call and this is what it always was.
      real(dp), intent(in) :: t2(:, :, :, :)
      real(dp), intent(in) :: c_vir(:, :)
      integer, intent(in) :: n_ao, n_o, n_v, j_lo, j_hi
      real(dp), allocatable, intent(out) :: half(:, :)

      real(dp), allocatable :: part(:, :, :, :), tmp(:, :)
      integer :: i, j, q, r, col, nj

      nj = j_hi - j_lo + 1

      ! part(i,j,p,q) = sum_ab t2(i,j,a,b) C_vir(p,a) C_vir(q,b)
      allocate (part(n_o, nj, n_ao, n_ao), tmp(n_v, n_ao))
      do j = 1, nj
         do i = 1, n_o
            tmp = matmul(t2(i, j_lo + j - 1, :, :), transpose(c_vir))
            part(i, j, :, :) = matmul(c_vir, tmp)
         end do
      end do
      deallocate (tmp)

      allocate (half(n_o, n_ao*n_ao*nj))
      !$omp parallel do collapse(3) default(none) &
      !$omp    shared(half, part, n_o, n_ao, nj) private(i, j, q, r, col)
      do j = 1, nj
         do r = 1, n_ao
            do q = 1, n_ao
               col = q + n_ao*(r - 1) + n_ao*n_ao*(j - 1)
               do i = 1, n_o
                  half(i, col) = 4.0_dp*part(i, j, q, r) - 2.0_dp*part(i, j, r, q)
               end do
            end do
         end do
      end do
      !$omp end parallel do
      deallocate (part)
   end subroutine build_half_block

   pure function occupied_per_block(n_ao, n_o, target_bytes) result(nj)
      !! How many occupied orbitals one pass over `j` may carry
      !!
      !! Sized on the two arrays that scale with it: the half-transformed
      !! amplitudes at `n_o n_ao^2` per `j`, and the transformation
      !! intermediate at `n_ao^3` per `j`. At least one, always -- a single
      !! occupied orbital is the smallest pass there is, and if that does not
      !! fit, nothing else here will either.
      integer, intent(in) :: n_ao, n_o
      real(dp), intent(in) :: target_bytes
      integer :: nj

      nj = int(target_bytes/((real(n_o, dp)*real(n_ao, dp)**2 &
                              + real(n_ao, dp)**3)*8.0_dp))
      nj = max(1, min(n_o, nj))
   end function occupied_per_block

   subroutine build_gamma_block(t2, c_occ, c_vir, n_ao, n_o, n_v, p_lo, p_hi, &
                                target_bytes, gamma_blk)
      !! The two-particle density for a range of its first index
      !!
      !! The same two contractions the dense path makes, with the first index
      !! restricted -- and the bra-pair term written out rather than obtained by
      !! transposing, because the transpose of a block is not a block: the
      !! `sum_i C(q,i) half(i,p,r,j)` term needs `p` restricted and `q` free,
      !! which is a different contraction rather than a rearrangement of the
      !! first one.
      !!
      !! Blocked over `j` as well, on the same target, so nothing here is
      !! `n_o^2 n_ao^2`. The price is that the half transform is redone for
      !! each block of the first index; it is `n_o^2 n_v^2 n_ao` against the
      !! `n_ao^4` of integrals a block already pays for, and it is what buys a
      !! footprint that does not grow with the basis.
      real(dp), intent(in) :: t2(:, :, :, :)
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :)
      integer, intent(in) :: n_ao, n_o, n_v, p_lo, p_hi
      real(dp), intent(in) :: target_bytes
      real(dp), allocatable, target, intent(out) :: gamma_blk(:, :, :, :)

      real(dp), allocatable, target :: xbuf(:)
      real(dp), allocatable :: half(:, :), hblk(:, :), tblk(:, :)
      real(dp), pointer :: x_first(:, :), x_second(:, :), gamma_flat(:, :)
      integer :: np, p, q, r, j, col, colb, j_lo, j_hi, nj, per_block

      np = p_hi - p_lo + 1
      per_block = occupied_per_block(n_ao, n_o, target_bytes)

      allocate (gamma_blk(np, n_ao, n_ao, n_ao))
      gamma_flat(1:np*n_ao*n_ao, 1:n_ao) => gamma_blk
      gamma_blk = 0.0_dp

      j_lo = 1
      do while (j_lo <= n_o)
         j_hi = min(n_o, j_lo + per_block - 1)
         nj = j_hi - j_lo + 1

         call build_half_block(t2, c_vir, n_ao, n_o, n_v, j_lo, j_hi, half)

         allocate (xbuf(int(np, kind=int64)*n_ao*n_ao*nj))
         x_first(1:np, 1:n_ao*n_ao*nj) => xbuf
         x_second(1:np*n_ao*n_ao, 1:nj) => xbuf

         ! term one: sum_i C(p,i) half(i, q, r, j), for p in the block
         call pic_gemm(c_occ(p_lo:p_hi, :), half, x_first)

         ! term two: sum_i C(q,i) half(i, p, r, j), same block of p. Gathered
         ! into `(i, (p r j))` first, because the block is a stride in `half`
         ! and BLAS wants it contiguous.
         allocate (hblk(n_o, np*n_ao*nj), tblk(n_ao, np*n_ao*nj))
         !$omp parallel do collapse(3) default(none) &
         !$omp    shared(hblk, half, n_o, n_ao, np, nj, p_lo) &
         !$omp    private(j, r, p, col, colb)
         do j = 1, nj
            do r = 1, n_ao
               do p = 1, np
                  colb = p + np*(r - 1) + np*n_ao*(j - 1)
                  col = (p_lo + p - 1) + n_ao*(r - 1) + n_ao*n_ao*(j - 1)
                  hblk(:, colb) = half(:, col)
               end do
            end do
         end do
         !$omp end parallel do
         deallocate (half)
         call pic_gemm(c_occ, hblk, tblk)
         deallocate (hblk)

         !$omp parallel do collapse(3) default(none) &
         !$omp    shared(x_first, tblk, n_ao, np, nj) private(j, r, p, q, col, colb)
         do j = 1, nj
            do r = 1, n_ao
               do q = 1, n_ao
                  col = q + n_ao*(r - 1) + n_ao*n_ao*(j - 1)
                  colb = np*(r - 1) + np*n_ao*(j - 1)
                  do p = 1, np
                     x_first(p, col) = x_first(p, col) + tblk(q, p + colb)
                  end do
               end do
            end do
         end do
         !$omp end parallel do
         deallocate (tblk)

         call pic_gemm(x_second, c_occ(:, j_lo:j_hi), gamma_flat, transb="T", &
                       beta=1.0_dp)
         deallocate (xbuf)
         j_lo = j_hi + 1
      end do
   end subroutine build_gamma_block

   subroutine build_two_particle_density(t2, c_occ, c_vir, n_ao, n_o, n_v, &
                                         target_bytes, gamma_ao)
      !! The non-separable two-particle density, in the AO basis
      !!
      !! **Two occupied transformations, as two matrix multiplies.** Written
      !! elementwise -- one small `matmul` per `(p,q,r,s)` -- this was 65% of
      !! the whole gradient, because `n_ao^4` is eleven million calls at
      !! cc-pVTZ and each one costs more in call overhead than in arithmetic.
      !! Contracted instead as `(n_ao x n_o) (n_o x n_ao^2 n_j)` and then
      !! `(n_ao^3 x n_j) (n_j x n_ao)`, the same flops sit in two gemms.
      !! Pointer bounds remapping is what lets the intermediate be read as
      !! `(p, q r j)` by the first and `(p q r, j)` by the second without a
      !! copy: the memory order is the same either way.
      !!
      !! `j` is contracted last, so a range of it gives an addend rather than a
      !! slice: the second gemm accumulates. One pass when the whole occupied
      !! range fits, which is the ordinary case.
      !!
      !! `n_ao^4` and dense in its result, which is the point: this is the path
      !! taken when the whole thing fits, where one gemm over the full first
      !! index beats the same gemm run block by block.
      !! `build_gamma_block` is the other path, and the driver picks on size.
      real(dp), intent(in) :: t2(:, :, :, :)
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :)
      integer, intent(in) :: n_ao, n_o, n_v
      real(dp), intent(in) :: target_bytes
      real(dp), allocatable, target, intent(out) :: gamma_ao(:, :, :, :)

      real(dp), allocatable, target :: xbuf(:)
      real(dp), allocatable :: half(:, :)
      real(dp), pointer :: x_first(:, :), x_second(:, :), gamma_flat(:, :)
      real(dp) :: acc
      integer :: p, q, r, s, j_lo, j_hi, nj, per_block

      per_block = occupied_per_block(n_ao, n_o, target_bytes)

      allocate (gamma_ao(n_ao, n_ao, n_ao, n_ao))
      gamma_flat(1:n_ao*n_ao*n_ao, 1:n_ao) => gamma_ao
      gamma_ao = 0.0_dp

      j_lo = 1
      do while (j_lo <= n_o)
         j_hi = min(n_o, j_lo + per_block - 1)
         nj = j_hi - j_lo + 1

         call build_half_block(t2, c_vir, n_ao, n_o, n_v, j_lo, j_hi, half)

         allocate (xbuf(int(n_ao, kind=int64)*n_ao*n_ao*nj))
         x_first(1:n_ao, 1:n_ao*n_ao*nj) => xbuf
         x_second(1:n_ao*n_ao*n_ao, 1:nj) => xbuf

         call pic_gemm(c_occ, half, x_first)
         deallocate (half)

         call pic_gemm(x_second, c_occ(:, j_lo:j_hi), gamma_flat, transb="T", &
                       beta=1.0_dp)
         deallocate (xbuf)
         j_lo = j_hi + 1
      end do

      ! The bra pair, symmetrised in place. This is the second of the two
      ! transformation terms -- `sum_i C(q,i) half(i,p,r,j)` is the first one
      ! with `p` and `q` exchanged -- so it needs no second contraction, only a
      ! transpose-add on each `(r,s)` slice.
      !
      ! **The ket pair is deliberately not symmetrised.** Both consumers
      ! contract this against something already symmetric in `(r,s)`: the
      ! integrals for the Lagrangian, and the differentiated integrals for the
      ! gradient, where the nabla sits on the bra and leaves `(rs|` alone. For a
      ! symmetric partner `sum B(r,s) X(s,r)` equals `sum B(r,s) X(r,s)`, so
      ! symmetrising here would compute a quantity that contracts to the same
      ! number -- at the cost of a second `n_ao^4` array and a strided pass over
      ! it. That pass, and the packing half it carried, cancel exactly against
      ! each other.
      !$omp parallel do collapse(2) default(none) &
      !$omp    shared(gamma_ao, n_ao) private(p, q, r, s, acc)
      do s = 1, n_ao
         do r = 1, n_ao
            do q = 1, n_ao
               do p = 1, q - 1
                  acc = gamma_ao(p, q, r, s) + gamma_ao(q, p, r, s)
                  gamma_ao(p, q, r, s) = acc
                  gamma_ao(q, p, r, s) = acc
               end do
               gamma_ao(q, q, r, s) = 2.0_dp*gamma_ao(q, q, r, s)
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine build_two_particle_density

   subroutine blocked_two_electron_terms(mol, t2, c_occ, c_vir, n_ao, n_o, n_v, &
                                         target_bytes, hf_density, imat_ao, de2, &
                                         vhf1, error)
      !! The Lagrangian and the gradient's two-electron term, a block at a time
      !!
      !! Neither the integrals nor the two-particle density is ever whole here.
      !! Both are built for one range of the first index, contracted, and
      !! discarded, so the peak is `BLOCK_TARGET` rather than `2 n_ao^4`.
      !!
      !! **What blocking costs.** The undifferentiated integrals are rebuilt for
      !! each block instead of being stored, and within a block only the ket
      !! pair's symmetry is available -- the first index is pinned to the block,
      !! which is what rules the other permutations out. That is roughly four
      !! times the integral work of the dense path's one symmetric build. It
      !! buys a ceiling set by choice instead of by the basis.
      !!
      !! Blocks are cut on shell boundaries, since a shell's functions share a
      !! quartet and cannot be split across two passes.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: t2(:, :, :, :), c_occ(:, :), c_vir(:, :)
      real(dp), intent(in) :: hf_density(:, :)
      integer, intent(in) :: n_ao, n_o, n_v
      real(dp), intent(in) :: target_bytes
      real(dp), intent(inout) :: imat_ao(:, :)
      real(dp), intent(inout) :: de2(:, :), vhf1(:, :, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable, target :: gamma_blk(:, :, :, :)
      real(dp), allocatable :: eri_blk(:, :, :, :), local(:, :)
      integer :: ish, ish_lo, ish_hi, p_lo, p_hi, np, per_block, r, s

      ! How many functions of the first index one block may carry, from the
      ! target and the two arrays that scale with it.
      per_block = max(1, int(target_bytes/(2.0_dp*real(n_ao, dp)**3*8.0_dp)))

      ish_lo = 1
      do while (ish_lo <= mol%nbas)
         ! Grow the block by whole shells until it would exceed the target.
         p_lo = mol%shell_offset(ish_lo) + 1
         ish_hi = ish_lo
         do ish = ish_lo, mol%nbas
            p_hi = mol%shell_offset(ish) + shell_dim(mol%cartesian, ish - 1, mol%bas)
            if (ish > ish_lo .and. p_hi - p_lo + 1 > per_block) exit
            ish_hi = ish
         end do
         p_hi = mol%shell_offset(ish_hi) + shell_dim(mol%cartesian, ish_hi - 1, mol%bas)
         np = p_hi - p_lo + 1

         call build_gamma_block(t2, c_occ, c_vir, n_ao, n_o, n_v, p_lo, p_hi, &
                                target_bytes, gamma_blk)

         allocate (eri_blk(np, n_ao, n_ao, n_ao))
         eri_blk = 0.0_dp
         call two_electron_mp2_terms(mol, gamma_blk, hf_density, de2, vhf1, &
                                     ish_lo=ish_lo, ish_hi=ish_hi, &
                                     p_offset=p_lo - 1, eri_blk=eri_blk)

         ! The Lagrangian, over this block of the contracted first index. Same
         ! product as the dense path, with `k` the block rather than `n_ao`.
         !$omp parallel default(none) shared(eri_blk, gamma_blk, imat_ao, n_ao) &
         !$omp    private(r, s, local)
         allocate (local(n_ao, n_ao))
         local = 0.0_dp
         !$omp do collapse(2) schedule(static)
         do s = 1, n_ao
            do r = 1, n_ao
               call pic_gemm(eri_blk(:, :, r, s), gamma_blk(:, :, r, s), local, &
                             transa="T", beta=1.0_dp)
            end do
         end do
         !$omp end do
         !$omp critical
         imat_ao = imat_ao + local
         !$omp end critical
         deallocate (local)
         !$omp end parallel

         deallocate (eri_blk, gamma_blk)
         ish_lo = ish_hi + 1
      end do

      if (error%has_error()) return
   end subroutine blocked_two_electron_terms

   subroutine contract_gamma_eri(eri, gamma_ao, n_ao, imat_ao)
      !! `Imat(p,q) = sum_{u,r,s} (u p|r s) gamma(u,q,r,s)`
      !!
      !! `n_ao^5` either way; the question is only whether those flops sit in a
      !! loop nest or in BLAS. Fixing the ket pair leaves an ordinary
      !! `(n_ao x n_ao) (n_ao x n_ao)` product of two contiguous slices -- the
      !! integrals are symmetric in `(u,p)`, so the slice can be read with
      !! either index first -- and the outer pair parallelises with no
      !! reduction beyond one accumulator per thread.
      real(dp), intent(in) :: eri(:, :, :, :), gamma_ao(:, :, :, :)
      integer, intent(in) :: n_ao
      real(dp), allocatable, intent(out) :: imat_ao(:, :)

      real(dp), allocatable :: local(:, :)
      integer :: r, s

      allocate (imat_ao(n_ao, n_ao))
      imat_ao = 0.0_dp

      !$omp parallel default(none) shared(eri, gamma_ao, imat_ao, n_ao) &
      !$omp    private(r, s, local)
      allocate (local(n_ao, n_ao))
      local = 0.0_dp
      !$omp do collapse(2) schedule(static)
      do s = 1, n_ao
         do r = 1, n_ao
            call pic_gemm(eri(:, :, r, s), gamma_ao(:, :, r, s), local, &
                          beta=1.0_dp)
         end do
      end do
      !$omp end do
      !$omp critical
      imat_ao = imat_ao + local
      !$omp end critical
      deallocate (local)
      !$omp end parallel
   end subroutine contract_gamma_eri

   subroutine two_electron_mp2_terms(mol, gamma_ao, hf_density, de2, vhf1, &
                                     ish_lo, ish_hi, p_offset, eri_blk)
      !! The differentiated integrals, contracted twice
      !!
      !! Once against the two-particle density, which gives a gradient
      !! contribution directly, and once against the reference density, which
      !! gives the per-atom matrices the relaxed density is contracted with
      !! afterwards. One quartet loop rather than two because the integrals are
      !! the expensive part and both contractions want the same ones.
      !!
      !! **One permutation survives, and it is used.** A differentiated quartet
      !! has none of the eightfold symmetry of an ordinary one: the nabla
      !! distinguishes the first index from the second and the bra from the
      !! ket. What is untouched is the ket pair, `(∇i j|k l) = (∇i j|l k)`, so
      !! this walks `l <= k` and writes the transposed contribution out rather
      !! than doubling -- the two-particle density is deliberately not
      !! symmetric in that pair, so `k` and `l` exchanged is a different
      !! contraction rather than the same number twice.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: gamma_ao(:, :, :, :), hf_density(:, :)
      real(dp), intent(inout) :: de2(:, :)        !! (3, natm), accumulated into
      real(dp), intent(inout) :: vhf1(:, :, :, :)  !! (nao, nao, 3, natm), likewise
      integer, intent(in), optional :: ish_lo, ish_hi
         !! Shells to differentiate over. Absent means all of them, which is the
         !! dense path; present is one block of the first index.
      integer, intent(in), optional :: p_offset
         !! What `gamma_ao`'s first index is offset by. Absent means zero.
      real(dp), intent(inout), optional, target :: eri_blk(:, :, :, :)
         !! Present, the undifferentiated integrals of the same quartets are
         !! stored here as they are computed. The blocked path needs them for
         !! the Lagrangian and cannot afford to keep the whole tensor; the dense
         !! path already has the tensor and leaves this absent, which is why the
         !! two share this loop rather than each having one.

      real(dp), allocatable :: buf(:)
      real(dp), allocatable :: de_local(:, :), vhf_local(:, :, :, :)
      integer, allocatable :: offsets(:), counts(:), shell_atom(:)
      type(c_ptr) :: opt
      integer :: shls(4)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: io, jo, ko, lo, i, j, k, l, comp, ret, mx, idx
      integer :: nao, nbas, natm, ia, first, last, off, gi
      real(dp) :: g, g0
      real(dp), allocatable :: eri_local(:)
      real(dp), pointer :: eri_out(:, :, :, :)
      logical :: want_eri

      nao = mol%nao
      nbas = mol%nbas
      natm = mol%natm
      mx = largest_shell(mol)

      allocate (offsets(natm), counts(natm))
      call atom_ao_blocks(mol, offsets, counts)
      allocate (shell_atom(nbas))
      do ish = 1, nbas
         io = mol%shell_offset(ish)
         shell_atom(ish) = 1
         do ia = 1, natm
            if (io >= offsets(ia) .and. io < offsets(ia) + counts(ia)) shell_atom(ish) = ia
         end do
      end do

      first = 1
      last = nbas
      if (present(ish_lo)) first = ish_lo
      if (present(ish_hi)) last = ish_hi
      off = 0
      if (present(p_offset)) off = p_offset
      ! Through a pointer rather than the dummy itself: an absent optional
      ! cannot appear in a data-sharing clause, and naming it there costs more
      ! than the branch it saves.
      want_eri = present(eri_blk)
      eri_out => null()
      if (want_eri) eri_out => eri_blk

      opt = c_null_ptr
      if (mol%cartesian) then
         call libcint_2e_ip1_cart_optimizer(opt, mol%atm, mol%natm, mol%bas, mol%nbas, mol%env)
      else
         call libcint_2e_ip1_sph_optimizer(opt, mol%atm, mol%natm, mol%bas, mol%nbas, mol%env)
      end if

      !$omp parallel default(none) &
      !$omp    shared(mol, gamma_ao, hf_density, de2, vhf1, opt, mx, nao, nbas, natm, &
      !$omp           shell_atom, first, last, off, want_eri, eri_out) &
      !$omp    private(ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo, &
      !$omp            i, j, k, l, comp, ret, idx, g, g0, gi, shls, ia, buf, &
      !$omp            de_local, vhf_local, eri_local)
      allocate (buf(mx**4*3))
      allocate (de_local(3, natm), vhf_local(nao, nao, 3, natm))
      allocate (eri_local(mx**4))
      de_local = 0.0_dp
      vhf_local = 0.0_dp

      !$omp do schedule(dynamic)
      do ish = first, last
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         ia = shell_atom(ish)
         do jsh = 1, nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            do ksh = 1, nbas
               dk = shell_dim(mol%cartesian, ksh - 1, mol%bas)
               ko = mol%shell_offset(ksh)
               do lsh = 1, ksh
                  dl = shell_dim(mol%cartesian, lsh - 1, mol%bas)
                  lo = mol%shell_offset(lsh)
                  shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]

                  if (mol%cartesian) then
                     ret = libcint_2e_ip1_cart(buf, shls, mol%atm, mol%natm, &
                                               mol%bas, nbas, mol%env, opt)
                  else
                     ret = libcint_2e_ip1_sph(buf, shls, mol%atm, mol%natm, &
                                              mol%bas, nbas, mol%env, opt)
                  end if
                  if (ret == 0) cycle

                  ! The same quartet undifferentiated, where the caller wants
                  ! it. Written straight into the block rather than accumulated,
                  ! because each quartet appears once in this loop.
                  if (want_eri) then
                     ret = two_electron_block(mol%cartesian, eri_local, shls, mol%atm, &
                                              mol%natm, mol%bas, nbas, mol%env)
                     if (ret /= 0) then
                        do l = 1, dl
                           do k = 1, dk
                              do j = 1, dj
                                 do i = 1, di
                                    ! libcint packs with the shells' own
                                    ! dimensions, not the buffer's, so this is
                                    ! indexed by hand rather than shaped.
                                    gi = i + di*(j - 1 + dj*(k - 1 + dk*(l - 1)))
                                    g0 = eri_local(gi)
                                    eri_out(io + i - off, jo + j, ko + k, lo + l) = g0
                                    if (lsh /= ksh) then
                                       eri_out(io + i - off, jo + j, lo + l, ko + k) = g0
                                    end if
                                 end do
                              end do
                           end do
                        end do
                     end if
                  end if

                  do comp = 1, 3
                     do l = 1, dl
                        do k = 1, dk
                           do j = 1, dj
                              do i = 1, di
                                 idx = i + di*(j - 1 + dj*(k - 1 + dk*(l - 1 + dl*(comp - 1))))
                                 g = buf(idx)

                                 de_local(comp, ia) = de_local(comp, ia) &
                                                      - 2.0_dp*g*gamma_ao(io + i - off, jo + j, ko + k, lo + l)

                                 vhf_local(ko + k, lo + l, comp, ia) = &
                                    vhf_local(ko + k, lo + l, comp, ia) &
                                    + g*hf_density(io + i, jo + j)
                                 vhf_local(ko + k, jo + j, comp, ia) = &
                                    vhf_local(ko + k, jo + j, comp, ia) &
                                    - 0.5_dp*g*hf_density(io + i, lo + l)
                                 vhf_local(io + i, jo + j, comp, ia) = &
                                    vhf_local(io + i, jo + j, comp, ia) &
                                    + g*hf_density(ko + k, lo + l)
                                 vhf_local(io + i, lo + l, comp, ia) = &
                                    vhf_local(io + i, lo + l, comp, ia) &
                                    - 0.5_dp*g*hf_density(jo + j, ko + k)

                                 ! The transposed ket pair: the same integral,
                                 ! and not the same contraction. Skipped on the
                                 ! diagonal, where it is the quartet already
                                 ! counted.
                                 if (lsh /= ksh) then
                                    de_local(comp, ia) = de_local(comp, ia) &
                                                         - 2.0_dp*g*gamma_ao(io + i - off, jo + j, lo + l, ko + k)
                                    vhf_local(lo + l, ko + k, comp, ia) = &
                                       vhf_local(lo + l, ko + k, comp, ia) &
                                       + g*hf_density(io + i, jo + j)
                                    vhf_local(lo + l, jo + j, comp, ia) = &
                                       vhf_local(lo + l, jo + j, comp, ia) &
                                       - 0.5_dp*g*hf_density(io + i, ko + k)
                                    vhf_local(io + i, jo + j, comp, ia) = &
                                       vhf_local(io + i, jo + j, comp, ia) &
                                       + g*hf_density(lo + l, ko + k)
                                    vhf_local(io + i, ko + k, comp, ia) = &
                                       vhf_local(io + i, ko + k, comp, ia) &
                                       - 0.5_dp*g*hf_density(jo + j, lo + l)
                                 end if
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end do

      !$omp critical
      de2 = de2 + de_local
      vhf1 = vhf1 + vhf_local
      !$omp end critical
      deallocate (buf, de_local, vhf_local, eri_local)
      !$omp end parallel

      call libcint_del_optimizer(opt)
   end subroutine two_electron_mp2_terms

   function largest_shell(mol) result(mx)
      !! The largest number of functions any one shell contributes
      type(libcint_molecule_t), intent(in) :: mol
      integer :: mx, ish

      mx = 1
      do ish = 1, mol%nbas
         mx = max(mx, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do
   end function largest_shell

end module mqc_libcint_mp2_gradient
