!! Contracts an effective density against the derivative integrals
module mqc_czt_density_contraction
   !! `contract_effective_density` turns any `effective_density_t` into one
   !! derivative per column, knowing nothing of the method behind it:
   !!
   !!     g(x) = tr(W S^x) + tr(P h^x) + sum_pairs D_l . (c J^x - x K^x)(D_r)
   !!          + tr(Gamma (mn|ls)^x) + tr(A S^x_bra)
   !!
   !! with `h^x` the full core-Hamiltonian derivative, Hellmann-Feynman term
   !! included. Every pool density, and every column's Gamma, goes through one
   !! `two_electron_deriv_many` sweep.
   !!
   !! Only these traces are formed. The nuclear repulsion, and anything else
   !! that is not a trace against an effective density, is the caller's.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, atom_ao_blocks
   use mqc_czt_gradient, only: one_electron_deriv, iprinv_deriv_at, two_electron_deriv_many, &
                               DERIV_OVLP, DERIV_KIN, DERIV_NUC
   use mqc_czt_effective_density, only: effective_density_t, separable_pair_t
   implicit none
   private

   public :: contract_effective_density

   real(dp), parameter :: GAMMA_BLOCK_TARGET = 2.0e8_dp
      !! Bytes of AO Gamma held per first-index block when the caller does
      !! not choose a width

contains

   subroutine contract_effective_density(mol, provider, gradients, error, gamma_width)
      !! Every column of `provider` contracted against the derivative
      !! integrals of `mol`, in Hartree/Bohr
      type(czt_molecule_t), intent(in) :: mol
      class(effective_density_t), intent(in) :: provider
      real(dp), allocatable, intent(out) :: gradients(:, :, :)
         !! (3, natm, n_columns); no nuclear repulsion
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: gamma_width
         !! Widest first-index block of Gamma to ask for, in compressed AOs.
         !! Absent: sized to hold `GAMMA_BLOCK_TARGET` bytes for every column.

      real(dp), allocatable :: s1(:, :, :), kin(:, :, :), h1(:, :, :), vrinv(:, :, :)
      real(dp), allocatable :: hcore_a(:, :), p_all(:, :, :), matrix(:, :)
      real(dp), allocatable :: pool(:, :, :), vhfs(:, :, :, :)
      real(dp), allocatable :: jd(:, :, :, :), kd(:, :, :, :), gamma_grads(:, :, :)
      integer, allocatable :: offsets(:), counts(:), ao_map(:)
      type(separable_pair_t), allocatable :: pairs(:)
      type(error_t) :: fill_error
      integer :: n_ao, natm, n_col, n_pool, n_sig, width, ic, ia, comp, ip, l, r, p0, p1
      real(dp) :: c, x

      if (error%has_error()) return
      n_ao = mol%nao
      natm = mol%natm
      n_col = provider%n_columns()
      if (n_col < 1) then
         call error%set(ERROR_VALIDATION, "contract_effective_density: the effective "// &
                        "density has no columns")
         return
      end if
      allocate (gradients(3, natm, n_col))
      gradients = 0.0_dp

      allocate (offsets(natm), counts(natm))
      call atom_ao_blocks(mol, offsets, counts)

      ! The sign convention is libcint's: its `ip` integrals carry a nabla on
      ! the bra, and the derivative with respect to the atom the bra sits on is
      ! minus that.
      call one_electron_deriv(mol, s1, DERIV_OVLP)
      s1 = -s1
      call one_electron_deriv(mol, kin, DERIV_KIN)
      call one_electron_deriv(mol, h1, DERIV_NUC)
      h1 = -(kin + h1)
      deallocate (kin)

      ! W, and the antisymmetric overlap term where a column has one, need
      ! only `s1`; P is kept for the per-atom core-Hamiltonian pass below.
      allocate (p_all(n_ao, n_ao, n_col))
      do ic = 1, n_col
         call provider%energy_weighted(ic, matrix, error)
         if (.not. matches(matrix, n_ao, "W", error)) return
         call add_bra_overlap(-2.0_dp, matrix, ic)

         call provider%overlap_antisymmetric(ic, matrix, error)
         if (error%has_error()) return
         if (allocated(matrix)) then
            if (.not. matches(matrix, n_ao, "the antisymmetric overlap term", error)) return
            call add_bra_overlap(1.0_dp, matrix, ic)
         end if

         call provider%one_particle(ic, matrix, error)
         if (.not. matches(matrix, n_ao, "P", error)) return
         p_all(:, :, ic) = matrix
      end do

      ! dh/dR_A: the nucleus of A moving (the Hellmann-Feynman term, every
      ! AO pair) and the basis functions on A moving (A's rows of `h1`).
      allocate (vrinv(n_ao, n_ao, 3), hcore_a(n_ao, n_ao))
      do ia = 1, natm
         call iprinv_deriv_at(mol, ia, vrinv)
         vrinv = -mol%charges(ia)*vrinv
         if (counts(ia) > 0) then
            p0 = offsets(ia) + 1
            p1 = offsets(ia) + counts(ia)
            vrinv(p0:p1, :, :) = vrinv(p0:p1, :, :) + h1(p0:p1, :, :)
         end if
         do comp = 1, 3
            hcore_a = vrinv(:, :, comp) + transpose(vrinv(:, :, comp))
            do ic = 1, n_col
               gradients(comp, ia, ic) = gradients(comp, ia, ic) + sum(hcore_a*p_all(:, :, ic))
            end do
         end do
      end do
      deallocate (vrinv, hcore_a, p_all, h1)

      ! One sweep for every pool density and every column's Gamma.
      call provider%density_pool(pool, error)
      if (error%has_error()) return
      if (.not. allocated(pool)) then
         call error%set(ERROR_VALIDATION, "contract_effective_density: no density pool")
         return
      end if
      if (size(pool, 1) /= n_ao .or. size(pool, 2) /= n_ao) then
         call error%set(ERROR_VALIDATION, "contract_effective_density: the density pool "// &
                        "does not match this basis")
         return
      end if
      n_pool = size(pool, 3)

      call provider%gamma_ao_map(ao_map)
      n_sig = 0
      if (allocated(ao_map)) then
         if (size(ao_map) > 0) n_sig = max(0, maxval(ao_map))
      end if
      allocate (gamma_grads(3, natm, n_col))
      gamma_grads = 0.0_dp
      if (n_sig > 0) then
         if (size(ao_map) /= n_ao) then
            call error%set(ERROR_VALIDATION, "contract_effective_density: the Gamma map "// &
                           "does not match this basis")
            return
         end if
         width = max(1, int(GAMMA_BLOCK_TARGET/(real(n_sig, dp)**3*8.0_dp*real(n_col, dp))))
         if (present(gamma_width)) width = gamma_width
         call two_electron_deriv_many(mol, pool, vhfs, error, gamma_ao_map=ao_map, &
                                      gamma_width=width, gamma_fill=fill_gamma, &
                                      gamma_grads=gamma_grads, coulomb_derivs=jd, &
                                      exchange_derivs=kd)
         if (fill_error%has_error()) then
            call error%set(fill_error%get_code(), fill_error%get_message())
            return
         end if
      else
         call two_electron_deriv_many(mol, pool, vhfs, error, coulomb_derivs=jd, &
                                      exchange_derivs=kd)
      end if
      if (error%has_error()) return
      deallocate (vhfs)

      ! A pair is symmetric in its two densities: each sees the other's field,
      ! and the factor two counts the ket half of the bra-only derivatives.
      do ic = 1, n_col
         call provider%separable_pairs(ic, pairs, error)
         if (error%has_error()) return
         if (.not. allocated(pairs)) cycle
         do ip = 1, size(pairs)
            l = pairs(ip)%left
            r = pairs(ip)%right
            c = pairs(ip)%coulomb
            x = pairs(ip)%exchange
            if (l < 1 .or. l > n_pool .or. r < 1 .or. r > n_pool) then
               call error%set(ERROR_VALIDATION, "contract_effective_density: a separable "// &
                              "pair indexes outside the density pool")
               return
            end if
            do ia = 1, natm
               if (counts(ia) == 0) cycle
               p0 = offsets(ia) + 1
               p1 = offsets(ia) + counts(ia)
               do comp = 1, 3
                  gradients(comp, ia, ic) = gradients(comp, ia, ic) + 2.0_dp*( &
                                            sum(pool(p0:p1, :, l)*(c*jd(p0:p1, :, comp, r) &
                                                                   - x*kd(p0:p1, :, comp, r))) &
                                            + sum(pool(p0:p1, :, r)*(c*jd(p0:p1, :, comp, l) &
                                                                     - x*kd(p0:p1, :, comp, l))))
               end do
            end do
         end do
      end do

      gradients = gradients + gamma_grads

   contains

      subroutine add_bra_overlap(scale, a, ic)
         !! `scale` times `a` against the bra derivative of the overlap
         real(dp), intent(in) :: scale
         real(dp), intent(in) :: a(:, :)
         integer, intent(in) :: ic
         integer :: ja, jc, q0, q1
         do ja = 1, natm
            if (counts(ja) == 0) cycle
            q0 = offsets(ja) + 1
            q1 = offsets(ja) + counts(ja)
            do jc = 1, 3
               gradients(jc, ja, ic) = gradients(jc, ja, ic) + scale*sum(s1(q0:q1, :, jc)*a(q0:q1, :))
            end do
         end do
      end subroutine add_bra_overlap

      subroutine fill_gamma(p_lo, p_hi, gamma)
         !! The engine's Gamma callback, forwarded to the provider
         integer, intent(in) :: p_lo, p_hi
         real(dp), allocatable, intent(inout) :: gamma(:, :, :, :, :)
         call provider%gamma_block(p_lo, p_hi, gamma, fill_error)
      end subroutine fill_gamma

   end subroutine contract_effective_density

   function matches(matrix, n_ao, what, error) result(ok)
      !! Whether a provider's matrix came back, at (n_ao, n_ao); sets `error`
      !! if not
      real(dp), allocatable, intent(in) :: matrix(:, :)
      integer, intent(in) :: n_ao
      character(len=*), intent(in) :: what
      type(error_t), intent(inout) :: error
      logical :: ok

      ok = .false.
      if (error%has_error()) return
      if (.not. allocated(matrix)) then
         call error%set(ERROR_VALIDATION, "contract_effective_density: no "//what// &
                        " for a column")
         return
      end if
      if (size(matrix, 1) /= n_ao .or. size(matrix, 2) /= n_ao) then
         call error%set(ERROR_VALIDATION, "contract_effective_density: "//what// &
                        " does not match this basis")
         return
      end if
      ok = .true.
   end function matches

end module mqc_czt_density_contraction
