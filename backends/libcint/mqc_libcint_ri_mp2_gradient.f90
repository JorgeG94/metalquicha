!! The RI-MP2 nuclear gradient over a closed-shell reference
module mqc_libcint_ri_mp2_gradient
   !! **What is differentiated.** An exact-ERI restricted Hartree-Fock
   !! reference with a *fitted* correlation energy -- which is what an `ri-mp2`
   !! deck asks for, the auxiliary basis naming the correlation fit and not the
   !! SCF. That is why the four-centre coupling tensor and the Fock derivative
   !! below are never fitted: approximating them would differentiate an energy
   !! nobody computed.
   !!
   !! The formulation is Stocks, Palethorpe and Barca, J. Chem. Theory Comput.
   !! 2024, 20, 2505, section 3. Only two things differ from the conventional
   !! MP2 gradient next door:
   !!
   !! * the two-particle density stays three-index, `Gamma^P_ia`, and contracts
   !!   against three- and two-centre derivative integrals rather than against
   !!   four-centre ones; and
   !! * the Lagrangian is built from `(pq|P)` instead of from a four-index
   !!   density.
   !!
   !! Everything one-particle -- the relaxed density, the energy-weighted
   !! density, the Z-vector, and the core-Hamiltonian, overlap and
   !! Hartree-Fock two-electron derivative terms -- is
   !! `mqc_libcint_mp2_gradient`'s, reused rather than rewritten.
   !!
   !! **Conventions that cost a day to establish**, none of them in the paper,
   !! all of them settled against the numpy prototype in
   !! `tools/cpu_validation/ri_mp2_gradient.py`:
   !!
   !! 1. The paper's `L'` and `L''` index the Lagrangian's own index first; the
   !!    assembly this feeds carries it second, so `imat` is their transpose.
   !!    The diagonals and the norms agree either way -- only the full matrix
   !!    comparison shows it.
   !! 2. The paper's `A_pqrs` contracted with a density is *twice* the response
   !!    operator the conventional assembly uses.
   !! 3. Terms 3 and 4 stayed wrong through three revisions while translational
   !!    invariance held at 1e-15. It is blind to every one of these.
   !!
   !! **Reproducibility under OpenMP.** This pipeline carries ~5e-11 of
   !! run-to-run scatter at more than one thread -- two 8-thread runs of
   !! `validation/check_ri_mp2_gradient` differ by that much, two
   !! single-threaded runs are bit-identical -- and the mechanism is known:
   !! the unordered `!$omp critical(mqc_direct_fock_accumulate)` merge of
   !! thread-local accumulators in `mqc_libcint_direct.f90`. Correctly
   !! synchronised, but whichever thread reaches the merge first adds first,
   !! so the summation order varies between runs. Measured, not guessed:
   !! `MKL_NUM_THREADS=1 OMP_NUM_THREADS=8` still varies while
   !! `OMP_NUM_THREADS=1 MKL_NUM_THREADS=8` is byte-identical, so it is our
   !! OpenMP and not threaded BLAS; and eight runs land on five contiguous
   !! values in the 11th digit, which is a handful of merge orders, not a
   !! race. The practical consequence: bit-identity claims about this
   !! gradient -- a refactor "changing nothing", a comparison "to rounding"
   !! -- are testable only at one thread.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, atom_ao_blocks, &
                                    build_df_mo_tensor, three_centre, two_centre, &
                                    metric_inverse_sqrt
   use mqc_libcint_cphf, only: cphf_solve, fitted_potential_general, IN_CORE_LIMIT
   use mqc_libcint_xc, only: xc_context_t, xc_kernel_apply
   use mqc_libcint_gradient, only: nuclear_repulsion_gradient, one_electron_deriv, &
                                   iprinv_deriv_at, three_centre_deriv, &
                                   two_centre_deriv, fitted_reference_gradient, &
                                   xc_potential_gradient, &
                                   DERIV_OVLP, DERIV_KIN, DERIV_NUC
   use mqc_libcint_direct, only: schwarz_bounds
   use mqc_libcint_mp2_gradient, only: gamma1_intermediates, two_electron_potential, &
                                       two_electron_mp2_terms
   implicit none
   private

   public :: libcint_ri_mp2_gradient

contains

   subroutine libcint_ri_mp2_gradient(mol, aux, coeff, orbital_energies, n_occ, &
                                      gradient, error, n_frozen, force_direct, &
                                      fitted_reference, xc, scf_density, pt2_scale)
      !! dE(RI-MP2)/dR for a closed-shell reference, in Hartree/Bohr
      type(libcint_molecule_t), intent(in) :: mol
      type(libcint_molecule_t), intent(in) :: aux   !! Auxiliary basis, same atoms
      real(dp), intent(in) :: coeff(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_frozen
         !! Core orbitals excluded from the correlation, exactly the
         !! conventional gradient's split: the fitted tensor, the amplitudes and
         !! both fitted densities span the active occupied space only, while the
         !! reference density, the response and every one-particle quantity keep
         !! the full occupied space.
      logical, intent(in), optional :: force_direct
         !! Recompute the reference integrals rather than storing them, whatever
         !! the size. The check harness passes it so both paths are exercised on
         !! a case small enough to take either -- a threshold nothing ever
         !! crosses in a test is a threshold nothing ever tests.
      logical, intent(in), optional :: fitted_reference
         !! The SCF underneath fitted its own J and K, so differentiate that
         !! rather than the exact reference.
         !!
         !! This is a different energy, not a cheaper route to the same one, and
         !! everything the reference touches moves with it: the response
         !! operator the Z-vector is solved against, the two potentials built
         !! from the relaxed density, and the two-electron derivative term,
         !! which stops being a four-centre contraction entirely. Absent is the
         !! exact-reference case an `ri-mp2` deck asks for -- the auxiliary basis
         !! naming only the correlation fit.
         !!
         !! One auxiliary basis serves both, which is what a deck can express:
         !! `model.aux_basis` is the only place a fitting set is named.
      type(xc_context_t), intent(inout), optional :: xc
         !! Present, this returns a *double hybrid's* perturbative term rather
         !! than an MP2 gradient -- the correlation alone, over a Kohn-Sham
         !! reference, scaled by `pt2_scale`. The same four differences the
         !! conventional gradient next door documents apply here:
         !!
         !! 1. Every response is the Kohn-Sham one: exact exchange scaled by the
         !!    functional's fraction, plus the exchange-correlation kernel.
         !! 2. The reference's own gradient is not assembled --
         !!    `libcint_scf_gradient` has already returned it.
         !! 3. The derivative of the exchange-correlation potential appears.
         !! 4. Everything is scaled once, at the end.
         !!
         !! What differs from that routine is what this one was for: the
         !! correlation is *fitted*, so the perturbative term costs `n^2 n_aux`
         !! rather than `n^4`. That is the whole reason a double hybrid gradient
         !! would come through here.
      real(dp), intent(in), optional :: scf_density(:, :)
         !! The converged reference density. Required with `xc`.
      real(dp), intent(in), optional :: pt2_scale
         !! The functional's PT2 coefficient, applied once and after the
         !! Z-vector solve -- scaling a linear system's right-hand side and
         !! scaling its solution are the same thing, and doing both is the error
         !! that hides at a coefficient of one.

      real(dp), allocatable :: bmat(:, :), bia_raw(:, :, :), metric(:, :), jm12(:, :)
      real(dp), allocatable :: ovov(:, :), xm(:, :)
      real(dp), allocatable :: t2(:, :, :, :)
      real(dp), allocatable :: gamma(:, :, :), gamma_pq(:, :), gamma_ao(:, :, :)
      real(dp), allocatable :: doo(:, :), dvv(:, :), dm1mo(:, :), zeta(:, :)
      real(dp), allocatable :: imat(:, :), im1(:, :), im1_t(:, :)
      real(dp), allocatable :: lag_o(:, :), lag_v(:, :), pq_p(:, :, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :), c_act(:, :)
      real(dp), allocatable :: hf_density(:, :), dm1(:, :), dm1_total(:, :), dm1p(:, :)
      real(dp), allocatable :: veff(:, :), xvo(:, :, :), zvec(:, :, :)
      real(dp), allocatable :: work(:, :), work_t(:, :), p_occ(:, :), vhf_s1occ(:, :)
      real(dp), allocatable :: s1(:, :, :), h1(:, :, :), kin(:, :, :), zero_h(:, :)
      real(dp), allocatable :: vrinv(:, :, :), hcore_a(:, :, :)
      real(dp), allocatable :: de2(:, :), vhf1(:, :, :, :)
      real(dp), allocatable, target :: eri(:, :, :, :)
      real(dp), allocatable :: bounds(:, :), gamma_null(:, :, :, :)
      real(dp), allocatable :: three_ao(:, :), ref_grad(:, :)
      real(dp), allocatable, target :: bref(:, :)
      real(dp), pointer :: bref_arg(:, :)
      real(dp), pointer :: eri_arg(:, :, :, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: n_ao, n_mo, n_o, n_oa, n_v, n_aux, n_ov, frozen
      integer :: i, j, a, b, p, q, iatom, comp, p0, p1
      real(dp) :: denom, kf, cscale
      logical :: dense, fitted, dh

      n_ao = mol%nao
      n_mo = size(coeff, 2)
      n_o = n_occ
      n_v = n_mo - n_occ

      frozen = 0
      if (present(n_frozen)) frozen = n_frozen
      if (frozen < 0 .or. frozen >= n_occ) then
         call error%set(ERROR_VALIDATION, "the RI-MP2 gradient's frozen core must "// &
                        "leave at least one occupied orbital to correlate")
         return
      end if
      n_oa = n_o - frozen
      if (n_v < 1 .or. n_o < 1) then
         call error%set(ERROR_VALIDATION, "the RI-MP2 gradient needs both occupied "// &
                        "and virtual orbitals")
         return
      end if

      ! Stored or recomputed, for the *reference* integrals only -- the ones the
      ! Z-vector's operator and the reference potential are built from. Those
      ! stay exact whatever the correlation is fitted with, because an `ri-mp2`
      ! energy fits nothing in its SCF, so the response has to be the exact
      ! reference's or it answers a question nobody asked.
      !
      ! What is fitted is not what makes this decision. Storing `n_ao^4` doubles
      ! would put the ceiling of a method chosen to avoid `n_ao^4` back where RI
      ! removed it, so past the limit the same integrals are recomputed instead
      ! and nothing about the answer changes.
      fitted = .false.
      if (present(fitted_reference)) fitted = fitted_reference

      ! The double hybrid's three settings, and the ordinary RI-MP2 values that
      ! make every expression below reduce to what it was.
      dh = present(xc)
      kf = 1.0_dp
      cscale = 1.0_dp
      if (dh) then
         if (.not. present(scf_density)) then
            call error%set(ERROR_VALIDATION, "a double hybrid gradient needs the "// &
                           "converged reference density: the exchange-correlation "// &
                           "kernel is evaluated at it and the potential differentiated "// &
                           "there")
            return
         end if
         kf = xc%exx_fraction
         if (present(pt2_scale)) cscale = pt2_scale
      end if

      dense = real(n_ao, dp)**4*8.0_dp <= IN_CORE_LIMIT
      if (present(force_direct)) then
         if (force_direct) dense = .false.
      end if
      ! With the reference fitted there are no four-centre integrals in the
      ! whole routine, so neither branch of the storage decision applies and
      ! both are left empty rather than chosen between.
      if (fitted) dense = .false.

      ! Two occupied spaces from here on, the conventional gradient's split:
      ! the full one, `c_occ`, is what the reference density, the response and
      ! every one-particle quantity run over; the active one, `c_act`, is all
      ! the fitted tensor and the amplitudes ever see. With no frozen core the
      ! two are the same columns.
      allocate (c_occ(n_ao, n_o), c_vir(n_ao, n_v), c_act(n_ao, n_oa))
      c_occ = coeff(:, 1:n_o)
      c_vir = coeff(:, n_o + 1:n_mo)
      c_act = coeff(:, frozen + 1:n_o)

      ! ---- the fitted integrals, and the amplitudes over them --------------
      ! `build_df_mo_tensor` lays its result out `(a, P, i)`, the shape the
      ! energy's gemms want. Everything below indexes `(i, a, P)` instead,
      ! matching the equations and the numpy reference, so it is repacked once
      ! here rather than transposed at every use.
      call build_df_mo_tensor(mol, aux, c_act, c_vir, bia_raw, error)
      if (error%has_error()) return
      n_aux = size(bia_raw, 2)
      n_ov = n_oa*n_v
      allocate (bmat(n_ov, n_aux))
      do p = 1, n_aux
         do i = 1, n_oa
            bmat(i:n_ov:n_oa, p) = bia_raw(:, p, i)
         end do
      end do
      deallocate (bia_raw)

      call two_centre(aux, metric)
      call metric_inverse_sqrt(metric, jm12, error)
      if (error%has_error()) return

      ! The raw three-centre integrals, built once. The Lagrangian below
      ! transforms them to `(pq|P)`, and a fitted reference needs them again --
      ! undoubled and unfitted -- for its own derivative term.
      call three_centre(mol, aux, three_ao)
      if (fitted) then
         ! `B(mu nu, P)`, the same tensor the SCF built. Formed from the
         ! integrals and metric already in hand rather than through
         ! `build_df_tensor`, which would repeat both.
         allocate (bref(n_ao*n_ao, n_aux))
         call pic_gemm(three_ao, jm12, bref)
      else
         allocate (bref(0, 0))
      end if

      ! `(ia|jb)_RI = sum_P B^P_ia B^P_jb` is one gemm once `B` is held as a
      ! matrix in the compound `(i,a)` index, which is the shape the repack
      ! above produces. Written as `n_o^2 n_v^2` separate dot products it is the
      ! same arithmetic against `B`'s longest stride, and it was the second most
      ! expensive step in the routine.
      allocate (ovov(n_ov, n_ov))
      call pic_gemm(bmat, bmat, ovov, transb="T")

      ! `i` and `j` count *active* occupied orbitals, which is why their
      ! energies are read at `frozen + i`: the fitted tensor already dropped
      ! the core, and the denominators have to drop the same orbitals or every
      ! amplitude divides by the wrong gap.
      allocate (t2(n_oa, n_oa, n_v, n_v))
      do b = 1, n_v
         do a = 1, n_v
            do j = 1, n_oa
               do i = 1, n_oa
                  denom = orbital_energies(frozen + i) + orbital_energies(frozen + j) &
                          - orbital_energies(n_occ + a) - orbital_energies(n_occ + b)
                  t2(i, j, a, b) = ovov(i + (a - 1)*n_oa, j + (b - 1)*n_oa)/denom
               end do
            end do
         end do
      end do
      deallocate (ovov)

      ! ---- the three- and two-index densities (eqs 8, 10) ------------------
      ! `2 X^ab_ij - X^ba_ij`, in the two compound indices, so that eq 8 is a
      ! gemm rather than the five-deep loop it reads as.
      allocate (xm(n_ov, n_ov))
      do b = 1, n_v
         do a = 1, n_v
            do j = 1, n_oa
               do i = 1, n_oa
                  xm(i + (a - 1)*n_oa, j + (b - 1)*n_oa) = &
                     2.0_dp*t2(i, j, a, b) - t2(i, j, b, a)
               end do
            end do
         end do
      end do

      call build_gamma(xm, bmat, jm12, n_oa, n_v, n_aux, gamma)
      call build_gamma_metric(gamma, bmat, jm12, n_oa, n_v, n_aux, gamma_pq)
      deallocate (xm)

      ! ---- the unrelaxed one-particle blocks -------------------------------
      call gamma1_intermediates(t2, n_oa, n_v, doo, dvv)
      deallocate (t2)

      ! The amplitude blocks land where the amplitudes live: `doo` in the
      ! *active* occupied rows and columns, `dvv` in the virtuals. The frozen
      ! rows and columns stay zero here -- their occupied-frozen block is not an
      ! amplitude quantity and is filled from the Lagrangian below, once the
      ! Lagrangian exists.
      allocate (dm1mo(n_mo, n_mo))
      dm1mo = 0.0_dp
      dm1mo(frozen + 1:n_o, frozen + 1:n_o) = doo + transpose(doo)
      dm1mo(n_o + 1:n_mo, n_o + 1:n_mo) = dvv + transpose(dvv)

      ! ---- the Lagrangian (eqs 13, 14) -------------------------------------
      !
      ! `(pq|P)` over every molecular orbital, so the three-centre integrals are
      ! transformed in full rather than only to the occupied-virtual block.
      allocate (pq_p(n_mo, n_mo, n_aux))
      call transform_three_centre(three_ao, coeff, n_ao, n_mo, n_aux, pq_p)

      ! `gamma`'s occupied index counts active orbitals, so the occupied
      ! Lagrangian has active rows and its `(iq|P)` reads at `frozen + i` --
      ! while the free index `q` still runs over every molecular orbital,
      ! which is where the occupied-frozen rows of `imat` come from.
      allocate (lag_o(n_oa, n_mo), lag_v(n_v, n_mo))
      lag_o = 0.0_dp
      lag_v = 0.0_dp
      do q = 1, n_mo
         do a = 1, n_v
            do i = 1, n_oa
               lag_o(i, q) = lag_o(i, q) &
                             + 2.0_dp*sum(gamma(:, i, a)*pq_p(q, n_o + a, :))
            end do
         end do
         do a = 1, n_v
            do i = 1, n_oa
               lag_v(a, q) = lag_v(a, q) &
                             + 2.0_dp*sum(gamma(:, i, a)*pq_p(frozen + i, q, :))
            end do
         end do
      end do
      deallocate (pq_p)

      ! Transposed into the slot: the paper indexes the Lagrangian's own index
      ! first, the assembly below carries it second. `imat`'s column is the
      ! fitted density's own MO index, so its frozen columns are zero --
      ! exactly what the conventional gradient's AO back-transform produces,
      ! where the two-particle density has no frozen component to land there.
      allocate (imat(n_mo, n_mo))
      imat = 0.0_dp
      imat(:, frozen + 1:n_o) = -transpose(lag_o)
      imat(:, n_o + 1:n_mo) = -transpose(lag_v)

      ! The occupied-frozen block of the relaxed density, the conventional
      ! gradient's construction verbatim. The amplitudes never touch a frozen
      ! orbital, but the Lagrangian does, and dividing its occupied-frozen
      ! elements by the orbital-energy gap is that rotation's own first-order
      ! response -- resolved directly because both orbitals are occupied and
      ! the coupled equations below span only occupied-virtual. It has to land
      ! before the `veff` build: the Z-vector's right-hand side is the response
      ! of the reference to the whole unrelaxed density, this block included.
      do i = frozen + 1, n_o
         do p = 1, frozen
            dm1mo(p, i) = imat(p, i)/(orbital_energies(p) - orbital_energies(i))
            dm1mo(i, p) = dm1mo(p, i)
         end do
      end do

      ! ---- the Z-vector ----------------------------------------------------
      allocate (hf_density(n_ao, n_ao))
      hf_density = 2.0_dp*matmul(c_occ, transpose(c_occ))

      if (dense) then
         call mol%eris(eri)
         allocate (bounds(0, 0))
      else
         allocate (eri(0, 0, 0, 0))
         call schwarz_bounds(mol, bounds, error)
         if (error%has_error()) return
      end if
      allocate (zero_h(n_ao, n_ao), veff(n_ao, n_ao))
      zero_h = 0.0_dp

      allocate (dm1(n_ao, n_ao))
      dm1 = matmul(coeff, matmul(dm1mo, transpose(coeff)))
      call reference_response(fitted, bref, dense, mol, eri, bounds, zero_h, dm1, &
                              veff, error, kf, dh, xc, scf_density)
      if (error%has_error()) return
      veff = 2.0_dp*veff

      allocate (xvo(n_v, n_o, 1))
      xvo(:, :, 1) = matmul(transpose(c_vir), matmul(veff, c_occ))
      do i = 1, n_o
         do a = 1, n_v
            xvo(a, i, 1) = xvo(a, i, 1) + imat(i, n_o + a) - imat(n_o + a, i)
         end do
      end do

      ! Tighter than the solver's default, for the reason the conventional
      ! gradient gives: the Z-vector residual enters the answer directly rather
      ! than quadratically. Handing over the tensor when there is one saves a
      ! second identical pass; when there is not, the solver goes direct with
      ! the same Schwarz bound the potentials above used.
      ! One integral source, never both: `cphf_solve` refuses a fitted tensor
      ! and an exact one together and is right to, since they are different
      ! operators. A fitted reference hands over `bmat` -- the operator the SCF
      ! actually used, because the orbitals are stationary for the fitted energy
      ! and for no other -- and withholds `eri`.
      !
      ! `xc` makes it the Kohn-Sham operator rather than the Hartree-Fock one.
      ! Omitting it there solves a different linear system, cleanly, and returns
      ! the answer to a question nobody asked.
      bref_arg => null()
      if (fitted) bref_arg => bref
      eri_arg => null()
      if (dense .and. .not. fitted) eri_arg => eri
      if (dh) then
         call cphf_solve(mol, coeff, orbital_energies, n_occ, response=zvec, &
                         error=error, mo_rhs=xvo, tol=1.0e-12_dp, max_iter=200, &
                         eri_in=eri_arg, bmat=bref_arg, xc=xc, density=scf_density)
      else
         call cphf_solve(mol, coeff, orbital_energies, n_occ, response=zvec, &
                         error=error, mo_rhs=xvo, tol=1.0e-12_dp, max_iter=200, &
                         eri_in=eri_arg, bmat=bref_arg)
      end if
      if (error%has_error()) return

      dm1mo(n_o + 1:n_mo, 1:n_o) = dm1mo(n_o + 1:n_mo, 1:n_o) + zvec(:, :, 1)
      dm1mo(1:n_o, n_o + 1:n_mo) = dm1mo(1:n_o, n_o + 1:n_mo) + transpose(zvec(:, :, 1))

      imat(n_o + 1:n_mo, 1:n_o) = transpose(imat(1:n_o, n_o + 1:n_mo))
      allocate (im1(n_ao, n_ao))
      im1 = matmul(coeff, matmul(imat, transpose(coeff)))

      ! ---- the energy-weighted density -------------------------------------
      allocate (zeta(n_mo, n_mo))
      do q = 1, n_mo
         do p = 1, n_mo
            zeta(p, q) = 0.5_dp*(orbital_energies(p) + orbital_energies(q))
         end do
      end do
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
      call reference_response(fitted, bref, dense, mol, eri, bounds, zero_h, &
                              dm1 + transpose(dm1), veff, error, kf, dh, xc, &
                              scf_density)
      if (error%has_error()) return
      allocate (vhf_s1occ(n_ao, n_ao))
      vhf_s1occ = matmul(p_occ, matmul(veff, p_occ))

      allocate (dm1p(n_ao, n_ao), dm1_total(n_ao, n_ao))
      ! A double hybrid leaves the reference's own contribution out of both:
      ! `libcint_scf_gradient` already returned it. Adding it here is the single
      ! most likely error in this routine -- the answer comes out roughly the
      ! reference gradient plus a small correction, right in shape and wrong in
      ! value, and only a finite difference catches it.
      if (dh) then
         dm1p = 2.0_dp*dm1
         dm1_total = dm1
      else
         dm1p = hf_density + 2.0_dp*dm1
         dm1_total = hf_density + dm1
      end if

      ! The Hartree-Fock energy-weighted density, on top of the correlation one
      ! already in `work`. A double hybrid's reference carries its own.
      if (.not. dh) then
         do i = 1, n_o
            do q = 1, n_ao
               do p = 1, n_ao
                  work(p, q) = work(p, q) &
                               + 2.0_dp*orbital_energies(i)*c_occ(p, i)*c_occ(q, i)
               end do
            end do
         end do
      end if

      ! ---- assembly --------------------------------------------------------
      allocate (gradient(3, mol%natm))
      gradient = 0.0_dp
      ! The nuclei repel each other in the reference energy, not in the
      ! correlation one, so a double hybrid's perturbative term has no such
      ! contribution.
      if (.not. dh) call nuclear_repulsion_gradient(mol, gradient)

      ! The derivative of the exchange-correlation potential, contracted with
      ! the relaxed correlation density. No counterpart in the Hartree-Fock
      ! assembly: there the reference operator is entirely two-electron, while a
      ! Kohn-Sham one keeps part of itself on a quadrature grid.
      if (dh) then
         call xc_potential_gradient(xc, mol, scf_density, dm1, gradient, error)
         if (error%has_error()) return
      end if

      allocate (de2(3, mol%natm), vhf1(n_ao, n_ao, 3, mol%natm))
      de2 = 0.0_dp
      vhf1 = 0.0_dp
      allocate (gamma_null(0, 0, 0, 0))
      if (fitted) then
         ! The reference's own two-electron gradient and its cross term with the
         ! relaxed density, both from the fitted integrals. The half is the same
         ! one the exact branch below carries implicitly: `two_electron_mp2_terms`
         ! builds two of the four differentiated positions and lets the
         ! integral's symmetry supply the rest, so what it hands back is already
         ! halved. This one is not, and says so.
         allocate (ref_grad(3, mol%natm))
         ref_grad = 0.0_dp
         call fitted_reference_gradient(mol, aux, three_ao, jm12, hf_density, &
                                        dm1p, ref_grad, k_scale=kf)
         gradient = gradient + 0.5_dp*ref_grad
      else
         ! Only the reference-density matrices: the fitted two-particle density
         ! contracts against three-centre derivatives instead, below.
         call two_electron_mp2_terms(mol, gamma_null, hf_density, de2, vhf1, &
                                     with_gamma=.false., k_scale=kf)
      end if

      call build_gamma_ao(gamma, c_act, c_vir, n_ao, n_oa, n_v, n_aux, gamma_ao)
      call fitted_two_particle_gradient(mol, aux, gamma_ao, gamma_pq, gradient)

      call one_electron_deriv(mol, s1, DERIV_OVLP)
      s1 = -s1
      call one_electron_deriv(mol, kin, DERIV_KIN)
      call one_electron_deriv(mol, h1, DERIV_NUC)
      h1 = -(kin + h1)
      deallocate (kin)

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
                                    + sum(hcore_a(:, :, comp)*dm1_total) &
                                    - sum(vhf1(:, :, comp, iatom)*dm1p)
         end do

         if (counts(iatom) == 0) cycle
         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + sum(s1(p0:p1, :, comp)*im1(p0:p1, :)) &
                                    + sum(s1(p0:p1, :, comp)*im1_t(p0:p1, :)) &
                                    - sum(s1(p0:p1, :, comp)*work(p0:p1, :)) &
                                    - sum(s1(p0:p1, :, comp)*work_t(p0:p1, :)) &
                                    - 2.0_dp*sum(s1(p0:p1, :, comp)*vhf_s1occ(p0:p1, :))
         end do
      end do

      ! Once, here, and nowhere else. The coefficient scales the correlation
      ! energy, so it scales every term derived from it exactly once -- and
      ! since the Z-vector solves a linear system, scaling its right-hand side
      ! and scaling its solution are the same operation. Doing both is the
      ! error, and it is invisible at a coefficient of one.
      if (dh) gradient = cscale*gradient

   end subroutine libcint_ri_mp2_gradient

   subroutine reference_response(fitted, b, dense, mol, eri, bounds, zero_h, &
                                 density, veff, error, k_scale, dh, xc, ref_density)
      !! The reference operator's response to a density, whichever reference
      !!
      !! Four combinations behind one call: fitted or exact integrals, and
      !! Hartree-Fock or Kohn-Sham. The kernel is an addend over whichever
      !! source ran rather than a fifth path -- `f_xc` is a grid quantity and
      !! knows nothing about where J and K came from.
      !!
      !! Whatever built the Z-vector's operator has to build this too: the
      !! Lagrangian's right-hand side and the linear system it is solved against
      !! belong to one reference, and mixing them converges to a plausible wrong
      !! Z-vector.
      logical, intent(in) :: fitted
      real(dp), intent(in) :: b(:, :)
      logical, intent(in) :: dense
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: eri(:, :, :, :), bounds(:, :)
      real(dp), intent(in) :: zero_h(:, :), density(:, :)
      real(dp), intent(out) :: veff(:, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in) :: k_scale
      logical, intent(in) :: dh
      type(xc_context_t), intent(inout), optional :: xc
      real(dp), intent(in), optional :: ref_density(:, :)

      if (fitted) then
         ! Not `build_fock_df`: `density` here is a relaxed correlation density,
         ! symmetric and indefinite, with no occupied orbitals to route exchange
         ! through.
         call fitted_potential_general(b, density, veff, k_scale=k_scale)
      else
         call two_electron_potential(dense, mol, eri, bounds, zero_h, density, veff, &
                                     error, k_scale=k_scale)
         if (error%has_error()) return
      end if
      if (.not. dh) return
      call xc_kernel_apply(xc, mol, ref_density, density, veff, error)
   end subroutine reference_response

   subroutine build_gamma(xm, bmat, jm12, n_o, n_v, n_aux, gamma)
      !! `Gamma^P_ia = sum_{bjQ} (2 X^ab_ij - X^ba_ij) B^Q_jb J^{-1/2}_PQ` (eq 8)
      !!
      !! Two gemms over the compound `(i,a)` index. The `sum_{jb}` is the first
      !! and the metric contraction is the second, and neither needs the rank-4
      !! form the equation is written in -- which is the whole reason `xm` and
      !! `bmat` are built as matrices upstream rather than reshaped here.
      real(dp), intent(in) :: xm(:, :)      !! `2X - X^swap` as `((i,a), (j,b))`
      real(dp), intent(in) :: bmat(:, :)    !! `B` as `((j,b), P)`
      real(dp), intent(in) :: jm12(:, :)
      integer, intent(in) :: n_o, n_v, n_aux
      real(dp), allocatable, intent(out) :: gamma(:, :, :)

      real(dp), allocatable :: half(:, :), gm(:, :)
      integer :: i, a, qq, n_ov

      n_ov = n_o*n_v
      allocate (half(n_ov, n_aux), gm(n_ov, n_aux))
      call pic_gemm(xm, bmat, half)
      call pic_gemm(half, jm12, gm, transb="T")
      deallocate (half)

      ! Back to `(P, i, a)`, the layout the Lagrangian and the AO
      ! back-transform read. One `n_ov` by `n_aux` transpose against two gemms
      ! of `n_ov^2 n_aux` and `n_ov n_aux^2`.
      allocate (gamma(n_aux, n_o, n_v))
      do a = 1, n_v
         do i = 1, n_o
            do qq = 1, n_aux
               gamma(qq, i, a) = gm(i + (a - 1)*n_o, qq)
            end do
         end do
      end do
      deallocate (gm)
   end subroutine build_gamma

   subroutine build_gamma_metric(gamma, bmat, jm12, n_o, n_v, n_aux, gamma_pq)
      !! `gamma_PQ = sum_{iaR} Gamma^P_ia B^R_ia J^{-1/2}_QR` (eq 10)
      !!
      !! Two gemms, in place of the `n_o n_v` rank-one updates this was: the
      !! `sum_ia` is an inner product over the compound index and the metric
      !! contraction follows it.
      real(dp), intent(in) :: gamma(:, :, :), bmat(:, :), jm12(:, :)
      integer, intent(in) :: n_o, n_v, n_aux
      real(dp), allocatable, intent(out) :: gamma_pq(:, :)

      real(dp), allocatable :: half(:, :), gflat(:, :)
      integer :: i, a, qq, n_ov

      n_ov = n_o*n_v
      allocate (gflat(n_ov, n_aux), half(n_aux, n_aux), gamma_pq(n_aux, n_aux))
      do a = 1, n_v
         do i = 1, n_o
            do qq = 1, n_aux
               gflat(i + (a - 1)*n_o, qq) = gamma(qq, i, a)
            end do
         end do
      end do
      call pic_gemm(gflat, bmat, half, transa="T")
      call pic_gemm(half, jm12, gamma_pq, transb="T")
      deallocate (half, gflat)
   end subroutine build_gamma_metric

   subroutine build_gamma_ao(gamma, c_occ, c_vir, n_ao, n_o, n_v, n_aux, gamma_ao)
      !! `Gamma^P_munu`, the three-index density back-transformed
      real(dp), intent(in) :: gamma(:, :, :), c_occ(:, :), c_vir(:, :)
      integer, intent(in) :: n_ao, n_o, n_v, n_aux
      real(dp), allocatable, intent(out) :: gamma_ao(:, :, :)

      integer :: pp

      allocate (gamma_ao(n_ao, n_ao, n_aux))
      !$omp parallel do default(none) &
      !$omp    shared(gamma_ao, gamma, c_occ, c_vir, n_aux) private(pp)
      do pp = 1, n_aux
         gamma_ao(:, :, pp) = matmul(c_occ, matmul(gamma(pp, :, :), transpose(c_vir)))
      end do
      !$omp end parallel do
   end subroutine build_gamma_ao

   subroutine transform_three_centre(three, coeff, n_ao, n_mo, n_aux, pq_p)
      !! `(mu nu|P)` to `(pq|P)`, both molecular indices
      real(dp), intent(in) :: three(:, :), coeff(:, :)
      integer, intent(in) :: n_ao, n_mo, n_aux
      real(dp), intent(out) :: pq_p(:, :, :)

      real(dp), allocatable :: block(:, :)
      integer :: pp

      allocate (block(n_ao, n_ao))
      !$omp parallel do default(none) &
      !$omp    shared(pq_p, three, coeff, n_ao, n_aux) private(pp, block)
      do pp = 1, n_aux
         block = reshape(three(:, pp), [n_ao, n_ao])
         pq_p(:, :, pp) = matmul(transpose(coeff), matmul(block, coeff))
      end do
      !$omp end parallel do
      deallocate (block)
   end subroutine transform_three_centre

   subroutine fitted_two_particle_gradient(mol, aux, gamma_ao, gamma_pq, gradient)
      !! The first two terms of eq 6, the only ones the fitting changes
      !!
      !! `4 sum Gamma^P_munu (P|munu)^xi` over both the orbital centres and the
      !! auxiliary one, then `-2 sum gamma_PQ (P|Q)^xi`. libcint's `ip`
      !! integrals differentiate the electronic coordinate, so the derivative
      !! with respect to a nucleus carrying that function is minus them.
      type(libcint_molecule_t), intent(in) :: mol, aux
      real(dp), intent(in) :: gamma_ao(:, :, :), gamma_pq(:, :)
      real(dp), intent(inout) :: gradient(:, :)

      real(dp), allocatable :: ip1(:, :, :, :), ip2(:, :, :, :), j1(:, :, :)
      real(dp), allocatable :: gt(:, :, :)
      integer, allocatable :: offsets(:), counts(:), aux_offsets(:), aux_counts(:)
      integer :: iatom, comp, p0, p1, q0, q1

      allocate (offsets(mol%natm), counts(mol%natm))
      allocate (aux_offsets(aux%natm), aux_counts(aux%natm))
      call atom_ao_blocks(mol, offsets, counts)
      call atom_ao_blocks(aux, aux_offsets, aux_counts)

      call three_centre_deriv(mol, aux, 1, ip1)
      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)
         if (counts(iatom) == 0) cycle
         ! Hoisted out of the component loop: it depends on the atom block and
         ! not on which Cartesian direction is being differentiated, so inside
         ! it was the same transpose and the same allocation done three times.
         gt = transpose_first_two(gamma_ao, p0, p1)
         do comp = 1, 3
            ! The bra index on this atom, and then the ket: `(mu nu|P)` is
            ! symmetric in mu and nu, so the second is the first transposed.
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    - 4.0_dp*sum(ip1(p0:p1, :, :, comp) &
                                                 *gamma_ao(p0:p1, :, :)) &
                                    - 4.0_dp*sum(ip1(p0:p1, :, :, comp)*gt)
         end do
         deallocate (gt)
      end do
      deallocate (ip1)

      call three_centre_deriv(mol, aux, 2, ip2)
      do iatom = 1, aux%natm
         q0 = aux_offsets(iatom) + 1
         q1 = aux_offsets(iatom) + aux_counts(iatom)
         if (aux_counts(iatom) == 0) cycle
         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    - 4.0_dp*sum(ip2(:, :, q0:q1, comp) &
                                                 *gamma_ao(:, :, q0:q1))
         end do
      end do
      deallocate (ip2)

      call two_centre_deriv(aux, j1)
      do iatom = 1, aux%natm
         q0 = aux_offsets(iatom) + 1
         q1 = aux_offsets(iatom) + aux_counts(iatom)
         if (aux_counts(iatom) == 0) cycle
         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + 2.0_dp*sum(j1(q0:q1, :, comp)*gamma_pq(q0:q1, :)) &
                                    + 2.0_dp*sum(j1(q0:q1, :, comp) &
                                                 *transpose(gamma_pq(:, q0:q1)))
         end do
      end do
      deallocate (j1)
   end subroutine fitted_two_particle_gradient

   function transpose_first_two(g, p0, p1) result(t)
      !! `g(nu, mu, P)` for `mu` in a block, as the ket-side contraction needs
      real(dp), intent(in) :: g(:, :, :)
      integer, intent(in) :: p0, p1
      real(dp), allocatable :: t(:, :, :)
      integer :: pp

      allocate (t(p1 - p0 + 1, size(g, 1), size(g, 3)))
      do pp = 1, size(g, 3)
         t(:, :, pp) = transpose(g(:, p0:p1, pp))
      end do
   end function transpose_first_two

end module mqc_libcint_ri_mp2_gradient
