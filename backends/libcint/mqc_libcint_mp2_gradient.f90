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
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim, atom_ao_blocks
   use mqc_libcint_mp2, only: transform_ovov
   use mqc_libcint_cphf, only: cphf_solve
   use mqc_libcint_rhf, only: build_fock
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

   real(dp), parameter :: IN_CORE_LIMIT = 4.0e9_dp
      !! Bytes. The same ceiling `mqc_libcint_cphf` applies to its own stored
      !! tensor, for the same reason and deliberately not a different number.

contains

   subroutine libcint_mp2_gradient(mol, coeff, orbital_energies, n_occ, gradient, &
                                   error, n_frozen)
      !! dE(MP2)/dR for a closed-shell reference, in Hartree/Bohr
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)            !! C, (n_ao, n_mo)
      real(dp), intent(in) :: orbital_energies(:)    !! (n_mo), Hartree
      integer, intent(in) :: n_occ                   !! Doubly occupied count
      real(dp), allocatable, intent(out) :: gradient(:, :)   !! (3, natm)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_frozen

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
      character(len=32) :: text, gigabytes

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

      ! Refused rather than attempted. Two dense `n_ao^4` arrays are live at
      ! once here -- the integrals and the two-particle density -- so the
      ! ceiling arrives as an allocation failure or an out-of-memory kill, with
      ! nothing to say which array or how far over it went. Blocking the
      ! two-particle density over its first index is what lifts this; until
      ! then the limit is worth stating in the one place that knows it.
      if (2.0_dp*real(n_ao, dp)**4*8.0_dp > IN_CORE_LIMIT) then
         write (text, "(i0)") n_ao
         write (gigabytes, "(f0.1)") 2.0_dp*real(n_ao, dp)**4*8.0_dp/1.0e9_dp
         call error%set(ERROR_VALIDATION, "the MP2 gradient stores the integrals and "// &
                        "the two-particle density as dense n^4 arrays, which for "// &
                        trim(text)//" basis functions is "//trim(gigabytes)//" GB. "// &
                        "Blocking over the auxiliary index is not implemented, so this "// &
                        "is refused rather than attempted. Use a smaller basis, or "// &
                        "fragment the system.")
         return
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
      call build_two_particle_density(t2, c_occ, c_vir, n_ao, n_o, n_v, gamma_ao)

      call mol%eris(eri)
      call contract_gamma_eri(eri, gamma_ao, n_ao, imat_ao)

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

      allocate (hf_density(n_ao, n_ao))
      hf_density = 2.0_dp*matmul(c_occ, transpose(c_occ))

      allocate (dm1(n_ao, n_ao))
      dm1 = matmul(coeff, matmul(dm1mo, transpose(coeff)))
      ! `build_fock` returns `H + J - K/2`; a zero core Hamiltonian leaves the
      ! two-electron part alone, which is what a response wants. Threaded over
      ! the target element, and already the hot routine of the CPHF solver.
      allocate (zero_h(n_ao, n_ao), veff(n_ao, n_ao))
      zero_h = 0.0_dp
      call build_fock(zero_h, eri, dm1, veff)
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
      call cphf_solve(mol, coeff, orbital_energies, n_occ, response=zvec, &
                      error=error, mo_rhs=xvo, tol=1.0e-12_dp, max_iter=200, &
                      eri_in=eri)
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
      call build_fock(zero_h, eri, dm1 + transpose(dm1), veff)
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

      ! ---- the differentiated integrals -------------------------------------
      call two_electron_mp2_terms(mol, gamma_ao, hf_density, de2, vhf1)

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

   subroutine build_two_particle_density(t2, c_occ, c_vir, n_ao, n_o, n_v, gamma_ao)
      !! The non-separable two-particle density, in the AO basis
      !!
      !! `4 t2(i,j,a,b) - 2 t2(i,j,b,a)` is the coefficient the energy gives
      !! each integral, and the rest is four index transformations and the two
      !! symmetrisations that make the object contract correctly against
      !! integrals carrying the permutational symmetry of the ERIs.
      !!
      !! **Two occupied transformations, as two matrix multiplies.** Written
      !! elementwise -- one small `matmul` per `(p,q,r,s)` -- this was 65% of
      !! the whole gradient, because `n_ao^4` is eleven million calls at
      !! cc-pVTZ and each one costs more in call overhead than in arithmetic.
      !! Contracted instead as `(n_ao x n_o) (n_o x n_ao^2 n_o)` and then
      !! `(n_ao^3 x n_o) (n_o x n_ao)`, the same flops sit in two gemms.
      !! Pointer bounds remapping is what lets the intermediate be read as
      !! `(p, q r j)` by the first and `(p q r, j)` by the second without a
      !! copy: the memory order is the same either way.
      !!
      !! `n_ao^4` and dense. That is what limits this to validation-sized
      !! systems; blocking over the first index is the way out, and is a
      !! separate change.
      real(dp), intent(in) :: t2(:, :, :, :)
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :)
      integer, intent(in) :: n_ao, n_o, n_v
      real(dp), allocatable, target, intent(out) :: gamma_ao(:, :, :, :)

      real(dp), allocatable :: part(:, :, :, :), tmp(:, :), half(:, :)
      real(dp), allocatable, target :: xbuf(:)
      real(dp), pointer :: x_first(:, :), x_second(:, :), gamma_flat(:, :)
      real(dp) :: acc
      integer :: i, j, p, q, r, s, col

      ! part(i,j,p,q) = sum_ab t2(i,j,a,b) C_vir(p,a) C_vir(q,b)
      allocate (part(n_o, n_o, n_ao, n_ao), tmp(n_v, n_ao))
      do j = 1, n_o
         do i = 1, n_o
            tmp = matmul(t2(i, j, :, :), transpose(c_vir))
            part(i, j, :, :) = matmul(c_vir, tmp)
         end do
      end do
      deallocate (tmp)

      ! half(i, (q, r, j)), the combination the energy weights each integral by
      allocate (half(n_o, n_ao*n_ao*n_o))
      !$omp parallel do collapse(3) default(none) &
      !$omp    shared(half, part, n_o, n_ao) private(i, j, q, r, col)
      do j = 1, n_o
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

      allocate (xbuf(int(n_ao, kind=8)*n_ao*n_ao*n_o))
      x_first(1:n_ao, 1:n_ao*n_ao*n_o) => xbuf
      x_second(1:n_ao*n_ao*n_ao, 1:n_o) => xbuf

      call pic_gemm(c_occ, half, x_first)
      deallocate (half)

      allocate (gamma_ao(n_ao, n_ao, n_ao, n_ao))
      gamma_flat(1:n_ao*n_ao*n_ao, 1:n_ao) => gamma_ao
      call pic_gemm(x_second, c_occ, gamma_flat, transb="T")
      deallocate (xbuf)

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

   subroutine two_electron_mp2_terms(mol, gamma_ao, hf_density, de2, vhf1)
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
      real(dp), allocatable, intent(out) :: de2(:, :)        !! (3, natm)
      real(dp), allocatable, intent(out) :: vhf1(:, :, :, :)  !! (nao, nao, 3, natm)

      real(dp), allocatable :: buf(:)
      real(dp), allocatable :: de_local(:, :), vhf_local(:, :, :, :)
      integer, allocatable :: offsets(:), counts(:), shell_atom(:)
      type(c_ptr) :: opt
      integer :: shls(4)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: io, jo, ko, lo, i, j, k, l, comp, ret, mx, idx
      integer :: nao, nbas, natm, ia
      real(dp) :: g

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

      allocate (de2(3, natm), vhf1(nao, nao, 3, natm))
      de2 = 0.0_dp
      vhf1 = 0.0_dp

      opt = c_null_ptr
      if (mol%cartesian) then
         call libcint_2e_ip1_cart_optimizer(opt, mol%atm, mol%natm, mol%bas, mol%nbas, mol%env)
      else
         call libcint_2e_ip1_sph_optimizer(opt, mol%atm, mol%natm, mol%bas, mol%nbas, mol%env)
      end if

      !$omp parallel default(none) &
      !$omp    shared(mol, gamma_ao, hf_density, de2, vhf1, opt, mx, nao, nbas, natm, &
      !$omp           shell_atom) &
      !$omp    private(ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo, &
      !$omp            i, j, k, l, comp, ret, idx, g, shls, ia, buf, de_local, vhf_local)
      allocate (buf(mx**4*3))
      allocate (de_local(3, natm), vhf_local(nao, nao, 3, natm))
      de_local = 0.0_dp
      vhf_local = 0.0_dp

      !$omp do schedule(dynamic)
      do ish = 1, nbas
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

                  do comp = 1, 3
                     do l = 1, dl
                        do k = 1, dk
                           do j = 1, dj
                              do i = 1, di
                                 idx = i + di*(j - 1 + dj*(k - 1 + dk*(l - 1 + dl*(comp - 1))))
                                 g = buf(idx)

                                 de_local(comp, ia) = de_local(comp, ia) &
                                                      - 2.0_dp*g*gamma_ao(io + i, jo + j, ko + k, lo + l)

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
                                                         - 2.0_dp*g*gamma_ao(io + i, jo + j, lo + l, ko + k)
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
      deallocate (buf, de_local, vhf_local)
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
