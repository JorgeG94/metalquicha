!! Analytic nuclear gradients on the CPU backend
module mqc_czt_gradient
   !! The derivative of a converged SCF energy with respect to the nuclear
   !! coordinates, for a restricted or unrestricted reference.
   !!
   !! The energy is stationary with respect to the orbitals, so their
   !! derivatives do not appear and nothing here solves a response equation.
   !! What is left is the derivative of the *integrals*, contracted with
   !! densities the SCF already produced.
   !!
   !! **The four terms.**
   !!
   !!   * nuclear repulsion, which is analytic and needs no integrals
   !!   * the core Hamiltonian, `Tr(D dH/dR)`
   !!   * the two-electron term, `Tr(D dG/dR)`
   !!   * the Pulay term, `-Tr(W dS/dR)`, from the basis functions moving with
   !!     their atoms. `W` is the energy-weighted density, not the density
   !!
   !! **libcint differentiates the bra, and only the bra.** Every `ip` entry
   !! point returns the derivative with respect to the first shell's centre;
   !! the ket's derivative is recovered from translational invariance rather
   !! than integrated for. That is what the transposes and the factors of two
   !! below are doing, and it is why the returned gradient sums to zero over
   !! atoms only if they are right -- which the tests check.
   use pic_types, only: dp
   use mqc_nuclear_repulsion, only: add_nuclear_repulsion_gradient
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_ecp, only: ecp_refuses_derivatives
   use mqc_czt_integrals, only: czt_molecule_t, shell_dim, max_block, atom_ao_blocks, &
                                build_df_shell_table, three_centre, two_centre, &
                                metric_inverse_sqrt, eri_shell_table_t, &
                                eri_shell_table, eri_schwarz_collapse
   use mqc_czt_ao, only: eval_ao_block, AO_POINT_BLOCK, AO_HESS_COMP, &
                         shell_extents, block_significant_aos, eval_rho
   use mqc_czt_xc, only: xc_context_t, xc_grid_lda_quantities, &
                         xc_grid_gga_quantities, xc_grid_kernel_quantities, &
                         ensure_nlc_grid
   use mqc_czt_vv10, only: vv10_nlc
   use mqc_dft_partition, only: becke_partition_derivatives
   use pic_blas_interfaces, only: pic_gemm
   use libcint_fortran, only: libcint_1e_ipovlp_sph, libcint_1e_ipovlp_cart, &
                              libcint_1e_ipkin_sph, libcint_1e_ipkin_cart, &
                              libcint_1e_ipnuc_sph, libcint_1e_ipnuc_cart, &
                              libcint_1e_iprinv_sph, libcint_1e_iprinv_cart, &
                              libcint_2e_ip1_sph, libcint_2e_ip1_cart, &
                              libcint_3c2e_ip1_sph, libcint_3c2e_ip1_cart, &
                              libcint_3c2e_ip2_sph, libcint_3c2e_ip2_cart, &
                              libcint_2c2e_ip1_sph, libcint_2c2e_ip1_cart, &
                              libcint_2e_ip1_sph_optimizer, libcint_2e_ip1_cart_optimizer, &
                              libcint_del_optimizer, LIBCINT_PTR_RINV_ORIG, &
                              LIBCINT_PTR_RANGE_OMEGA, LIBCINT_PTR_EXP, &
                              LIBCINT_NPRIM_OF, LIBCINT_ATOM_OF
   use mqc_czt_direct, only: schwarz_bounds, block_density_max, pair_degeneracy, &
                             pair_work_order
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
   implicit none
   private

   public :: czt_scf_gradient
   public :: nuclear_repulsion_gradient   !! Exposed for the tests
   public :: three_centre_deriv           !! Exposed for the tests
   public :: two_centre_deriv             !! Exposed for the tests
   ! The pieces an MP2 gradient assembles for itself. It cannot call
   ! `czt_scf_gradient` and add to the answer -- the relaxed density enters
   ! the one-electron and overlap terms in place of the reference one, not
   ! alongside it.
   public :: one_electron_deriv
   public :: two_electron_deriv
      !! Exposed for the MCSCF gradient, whose inactive and active blocks are
      !! separable and contract exactly as a closed-shell reference does.
   public :: two_electron_gradient
   public :: hellmann_feynman_gradient
   public :: iprinv_deriv_at
   public :: xc_potential_gradient
   public :: vv10_gradient_fixed_grid     !! The VV10 Hessian differences this
   public :: fitted_reference_gradient
   public :: DERIV_OVLP, DERIV_KIN, DERIV_NUC

   ! Which one-electron derivative `one_electron_deriv` should build.
   integer, parameter :: DERIV_OVLP = 1
   integer, parameter :: DERIV_KIN = 2
   integer, parameter :: DERIV_NUC = 3

   real(dp), parameter :: DERIV_SCREEN_TOL = 1.0e-10_dp
      !! Contribution below which a differentiated shell quartet is dropped.
      !!
      !! Set two decades tighter than the accuracy actually wanted, because the
      !! bra bound it multiplies is an estimate rather than a rigorous
      !! Cauchy-Schwarz one -- see `two_electron_deriv`.

   real(dp), parameter :: GRAD_SCREEN_TOL = 1.0e-12_dp
      !! Contribution below which `two_electron_gradient` skips one derivative
      !! integral.
      !!
      !! Two decades inside `DERIV_SCREEN_TOL` because the quantity it bounds
      !! is two decades smaller: the contribution of an integral to the
      !! gradient itself, weighted by a product of two density elements, where
      !! `two_electron_deriv` bounds an element of `vhf` weighted by one.
      !! Measured on a 74-atom 6-31G* Hartree-Fock gradient, this value
      !! reproduces the unscreened result to 1e-9 Hartree/Bohr, which is what
      !! the older loop gave at its own default; 1e-10 here would be 1e-7.

contains

   subroutine czt_scf_gradient(mol, density, density_beta, orbitals, orbitals_beta, &
                               orbital_energies, orbital_energies_beta, &
                               n_occupied, n_occupied_beta, gradient, error, aux, xc)
      !! dE/dR for a converged SCF, in Hartree/Bohr
      !!
      !! The beta arguments present mean unrestricted: the Coulomb field is
      !! built from the total density and exchange from each spin separately,
      !! which is the only place the two paths differ. Absent means closed
      !! shell, where the density already carries its factor of two.
      !!
      !! `aux` present means the SCF being differentiated fitted its J and K,
      !! and only the two-electron term changes. **Passing an auxiliary basis
      !! the SCF did not use differentiates a different energy**, and nothing
      !! here can detect it.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)          !! Alpha, or the total for RHF
      real(dp), intent(in), optional :: density_beta(:, :)
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in), optional :: orbitals_beta(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      real(dp), intent(in), optional :: orbital_energies_beta(:)
      integer, intent(in) :: n_occupied
      integer, intent(in), optional :: n_occupied_beta
      real(dp), allocatable, intent(out) :: gradient(:, :)   !! (3, natm)
      type(error_t), intent(inout) :: error
      type(czt_molecule_t), intent(in), optional :: aux
         !! The auxiliary basis the SCF fitted with. Absent means exact ERIs.
      type(xc_context_t), intent(inout), optional :: xc
         !! The functional the SCF used. Absent means Hartree-Fock. Present adds
         !! the exchange-correlation term and scales exact exchange by the
         !! fraction the functional asks for -- zero for a pure functional, so
         !! the Fock exchange derivative is not merely small but absent.

      real(dp), allocatable :: total_density(:, :), weighted(:, :)
      real(dp), allocatable :: s1(:, :, :), h1(:, :, :), kin(:, :, :)
      real(dp) :: exx
      integer, allocatable :: offsets(:), counts(:)
      integer :: nao, iatom, comp, p0, p1
      logical :: unrestricted
      logical :: fitted

      if (ecp_refuses_derivatives(mol%core_electrons, "nuclear gradient", error)) return

      unrestricted = present(density_beta)
      nao = mol%nao

      ! Hartree-Fock keeps all of its exchange; a functional keeps the fraction
      ! it declares.
      exx = 1.0_dp
      if (present(xc)) exx = xc%exx_fraction

      if (size(density, 1) /= nao) then
         call error%set(ERROR_VALIDATION, &
                        "czt_scf_gradient: the density does not match this basis")
         return
      end if

      allocate (gradient(3, mol%natm))
      gradient = 0.0_dp

      call nuclear_repulsion_gradient(mol, gradient)

      ! The total density drives the Coulomb field and the core Hamiltonian
      ! trace; exchange is per spin and is handled inside the two-electron
      ! contraction.
      allocate (total_density(nao, nao), source=density)
      if (unrestricted) total_density = total_density + density_beta

      call energy_weighted_density(orbitals, orbital_energies, n_occupied, &
                                   .not. unrestricted, weighted)
      if (unrestricted) then
         block
            real(dp), allocatable :: weighted_beta(:, :)
            call energy_weighted_density(orbitals_beta, orbital_energies_beta, &
                                         n_occupied_beta, .false., weighted_beta)
            weighted = weighted + weighted_beta
         end block
      end if

      ! The sign convention is libcint's: its `ip` integrals carry a nabla on
      ! the bra, and the derivative of the integral with respect to the atom
      ! the bra sits on is minus that.
      call one_electron_deriv(mol, s1, DERIV_OVLP)
      s1 = -s1
      call one_electron_deriv(mol, kin, DERIV_KIN)
      call one_electron_deriv(mol, h1, DERIV_NUC)
      h1 = -(kin + h1)
      deallocate (kin)

      ! The fitted two-electron gradient does not factor into a `vhf` matrix
      ! the way the exact one does -- it contracts three- and two-centre
      ! derivatives against densities of its own -- so it accumulates straight
      ! into `gradient` here and the per-atom loop below skips its term.
      fitted = present(aux)
      if (fitted) then
         if (unrestricted) then
            call df_two_electron_gradient(mol, aux, total_density, orbitals, n_occupied, &
                                          orbitals_beta=orbitals_beta, &
                                          n_occupied_beta=n_occupied_beta, &
                                          gradient=gradient, error=error, &
                                          exx_fraction=exx)
         else
            call df_two_electron_gradient(mol, aux, total_density, orbitals, n_occupied, &
                                          gradient=gradient, error=error, &
                                          exx_fraction=exx)
            if (error%has_error()) return
            ! The long-range exchange, differentiated against the kernel it was
            ! fitted with. Exchange only: `J` is complete from the pass above,
            ! and the attenuated tensor's Coulomb term belongs to no energy.
            if (present(xc)) then
               if (xc%range_separated) then
                  call df_two_electron_gradient(mol, aux, total_density, orbitals, &
                                                n_occupied, gradient=gradient, &
                                                error=error, exx_fraction=xc%rs_k_lr, &
                                                rs_omega=xc%rs_omega, with_coulomb=.false.)
               end if
            end if
         end if
         if (error%has_error()) return
      else
         ! The exact two-electron term accumulates straight into `gradient` as
         ! well, one pass over the unique shell quartets whichever way the
         ! spins go: the Coulomb field comes from the total density and
         ! exchange from each spin, or from the total carrying its factor of
         ! two at a closed shell.
         if (unrestricted) then
            call two_electron_gradient(mol, total_density, gradient, error, &
                                       density_alpha=density, density_beta=density_beta, &
                                       exx_fraction=exx)
         else
            call two_electron_gradient(mol, total_density, gradient, error, &
                                       exx_fraction=exx)
         end if
         if (error%has_error()) return
         ! The long-range half of a range-separated functional's exchange, at
         ! the screened omega. Exchange only: `J` is complete from the pass
         ! above, and the attenuated kernel carries no Coulomb term of its own.
         if (present(xc) .and. .not. unrestricted) then
            if (xc%range_separated) then
               call two_electron_gradient(mol, total_density, gradient, error, &
                                          exx_fraction=xc%rs_k_lr, omega=xc%rs_omega, &
                                          with_coulomb=.false.)
               if (error%has_error()) return
            end if
         end if
      end if

      ! dH/dR_A has two parts that look alike and are not. Moving atom A moves
      ! the *nucleus*, which every electron feels wherever its orbital sits:
      ! the Hellmann-Feynman term, with no basis-set derivative in it, and one
      ! pass over the shell pairs for every atom at once. Moving atom A also
      ! moves the basis functions centred on it, which is the `h1` block in the
      ! per-atom loop below.
      call hellmann_feynman_gradient(mol, total_density, gradient)

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)

      do iatom = 1, mol%natm
         if (counts(iatom) == 0) cycle
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)

         ! Twice, because the nabla was on the bra: the ket's contribution is
         ! the same number by the symmetry of the integrals, so it is counted
         ! rather than computed.
         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + 2.0_dp*sum(h1(p0:p1, :, comp)*total_density(p0:p1, :)) &
                                    - 2.0_dp*sum(s1(p0:p1, :, comp)*weighted(p0:p1, :))
         end do
      end do

      ! Last, and on top of everything above: the exchange-correlation term is
      ! an addition to the Hartree-Fock gradient, not a replacement for any part
      ! of it.
      if (present(xc)) then
         ! Both closed-shell paths build the second exchange derivative at the
         ! screened omega. Neither long-range pass is written in an unrestricted
         ! branch, and without this refusal an open-shell range-separated
         ! gradient would come back missing its dominant exchange term,
         ! converged and unflagged.
         if (xc%range_separated .and. unrestricted) then
            call error%set(ERROR_VALIDATION, "the gradient of a range-separated "// &
                           "functional needs a second exchange derivative at the "// &
                           "screened omega, and the open-shell path does not build "// &
                           "one. The closed-shell path does, fitted or exact. Take "// &
                           "the gradient at a closed shell, or at a functional that "// &
                           "is not range separated.")
            return
         end if
         if (unrestricted) then
            call xc_gradient(xc, mol, density, gradient, error, density_beta=density_beta)
         else
            ! The closed-shell density carries its factor of two, and
            ! `xc_grid_lda_quantities` builds an unpolarised rho from it, which
            ! is what the unpolarised functional expects.
            call xc_gradient(xc, mol, total_density, gradient, error)
         end if
         if (error%has_error()) return
      end if

   end subroutine czt_scf_gradient

   subroutine xc_gradient(ctx, mol, density, gradient, error, density_beta)
      !! The exchange-correlation contribution to dE/dR, accumulated in place
      !!
      !! Three terms, and the last two are the ones a first attempt leaves out.
      !!
      !! **The basis functions move.** rho is built from functions centred on
      !! atoms, so moving atom A changes rho wherever A's functions reach. This
      !! is the analogue of the one-electron Pulay term.
      !!
      !! **The grid points move.** The quadrature is atom-centred, so the points
      !! belonging to A travel with it. This contributes the same quantity as
      !! the first term but summed over A's *points* rather than A's
      !! *functions* -- over every basis function, not only the ones on A -- and
      !! with the opposite sign.
      !!
      !! **The partition weights move.** The Becke weight of every point depends
      !! on every nuclear position, so moving A reweights the whole grid,
      !! including points it does not own. Its absence leaves a gradient that
      !! looks entirely plausible and misses finite differences at about 1e-4,
      !! and it is the term translational invariance is sensitive to.
      !!
      !! LDA, GGA and meta-GGA, restricted or unrestricted; the non-local (VV10)
      !! term goes through `vv10_gradient` and is refused for an open shell.
      type(xc_context_t), intent(inout) :: ctx
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(inout) :: gradient(:, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: density_beta(:, :)

      real(dp), allocatable :: rho(:), exc(:), vrho(:), vtau(:)
      real(dp), allocatable :: rho_beta(:), vrho_beta(:)
      real(dp), allocatable :: gcoef(:, :), gcoef_beta(:, :)
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), allocatable :: dchi(:, :), dchi_beta(:, :)
      real(dp), allocatable :: dgchi(:, :, :), dgchi_beta(:, :, :)
      logical :: gga, mgga
      real(dp), allocatable :: dpart(:, :, :)
      real(dp), allocatable :: extents(:), d_sig(:, :), db_sig(:, :)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_list(:), ao_offset(:)
      integer, allocatable :: c_offsets(:), c_counts(:)
      real(dp) :: contrib(3)
      real(dp) :: wv, scale
      integer :: npts, nao, natm, g0, g1, nb, ig, gg, comp, ia, own
      integer :: n_sig, isig, jsig
      logical :: unrestricted

      unrestricted = present(density_beta)

      ! The non-local term, where the functional carries one. Refused for an
      ! open shell rather than omitted: omitting it survives every check a user
      ! can make -- the SCF converges, the energy is right, and only the forces
      ! are wrong.
      if (ctx%nlc_b /= 0.0_dp .or. ctx%nlc_c /= 0.0_dp) then
         if (unrestricted) then
            call error%set(ERROR_VALIDATION, "this functional carries a non-local "// &
                           "correlation term (VV10), whose nuclear gradient is "// &
                           "implemented for a closed shell only. Refused rather than "// &
                           "returning a gradient missing it.")
            return
         end if
         call vv10_gradient(ctx, mol, density, gradient, error)
         if (error%has_error()) return
      end if

      if (.not. ctx%active) return
      nao = mol%nao
      natm = mol%natm
      npts = ctx%grid%n_points

      gga = ctx%any_gga .or. ctx%any_mgga
      mgga = ctx%any_mgga

      if (mgga) then
         ! `vtau` is what makes this legal: `xc_grid_gga_quantities` refuses a
         ! meta-GGA outright when the caller cannot take it, because its own
         ! per-functional dispatch would otherwise fall through to the LDA
         ! branch. `density_beta` is passed where there is one, so that the
         ! meta-GGA refusal inside fires for the reason it was written rather
         ! than being caught later by the spin-consistency check.
         if (unrestricted) then
            call xc_grid_gga_quantities(ctx, mol, density, rho, exc, vrho, gcoef, error, &
                                        density_beta=density_beta, rho_beta=rho_beta, &
                                        vrho_beta=vrho_beta, grad_coeff_beta=gcoef_beta, &
                                        vtau=vtau)
         else
            call xc_grid_gga_quantities(ctx, mol, density, rho, exc, vrho, gcoef, error, &
                                        vtau=vtau)
         end if
      else if (gga) then
         if (unrestricted) then
            call xc_grid_gga_quantities(ctx, mol, density, rho, exc, vrho, gcoef, error, &
                                        density_beta=density_beta, rho_beta=rho_beta, &
                                        vrho_beta=vrho_beta, grad_coeff_beta=gcoef_beta)
         else
            call xc_grid_gga_quantities(ctx, mol, density, rho, exc, vrho, gcoef, error)
         end if
      else
         if (unrestricted) then
            call xc_grid_lda_quantities(ctx, mol, density, rho, exc, vrho, error, &
                                        density_beta=density_beta, rho_beta=rho_beta, &
                                        vrho_beta=vrho_beta)
         else
            call xc_grid_lda_quantities(ctx, mol, density, rho, exc, vrho, error)
         end if
      end if
      if (error%has_error()) return

      allocate (c_offsets(natm), c_counts(natm))
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(nao))
      allocate (d_sig(nao, nao))
      if (unrestricted) allocate (db_sig(nao, nao))
      call shell_extents(mol, ctx%screen_tol, extents)

      ! A closed-shell density already carries its factor of two, so the
      ! derivative of rho with respect to a moving basis function is 2*D*chi
      ! either way -- the two is the bra/ket pair, not the occupation.
      scale = 2.0_dp

      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         ! Which functions reach this block, and where each atom's kept ones
         ! begin in the compressed numbering.
         !
         ! **Both halves of every accumulation below run over that same
         ! compressed set** -- the per-atom ranges and the sum over every
         ! function -- and those ranges tile 1..n_sig without a gap. That is
         ! what keeps the two opposite-signed terms cancelling and the gradient
         ! translationally invariant; screening them against different sets
         ! would leave a small spurious net force rather than crash.
         call block_significant_aos(mol, ctx%grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig, &
                                    atom_offsets=c_offsets, atom_counts=c_counts)

         ! `n_sig == 0` is empty space that no basis function reaches, so the
         ! basis-function and grid-motion terms are zero there. The
         ! partition-weight term below is deliberately not skipped with them:
         ! it depends on the grid and the nuclei, not on the basis.
         if (n_sig > 0) then
            do jsig = 1, n_sig
               do isig = 1, n_sig
                  d_sig(isig, jsig) = density(ao_list(isig), ao_list(jsig))
               end do
            end do
            if (unrestricted) then
               do jsig = 1, n_sig
                  do isig = 1, n_sig
                     db_sig(isig, jsig) = density_beta(ao_list(isig), ao_list(jsig))
                  end do
               end do
            end if

            if (gga) then
               call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, &
                                  grad=ao_grad, hess=ao_hess, shell_mask=shell_mask, &
                                  ao_offset=ao_offset, n_ao_out=n_sig)
            else
               call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, &
                                  grad=ao_grad, shell_mask=shell_mask, &
                                  ao_offset=ao_offset, n_ao_out=n_sig)
            end if
            if (error%has_error()) return

            ! (D chi)_mu(g) = sum_nu D_mu,nu chi_nu(g), the partner every term
            ! below contracts the basis-function gradient against.
            if (allocated(dchi)) deallocate (dchi)
            allocate (dchi(nb, n_sig))
            call density_times_ao(ao, d_sig(1:n_sig, 1:n_sig), nb, n_sig, dchi)
            if (unrestricted) then
               if (allocated(dchi_beta)) deallocate (dchi_beta)
               allocate (dchi_beta(nb, n_sig))
               call density_times_ao(ao, db_sig(1:n_sig, 1:n_sig), nb, n_sig, dchi_beta)
            end if

            ! (D grad chi)_mu(g), which only the GGA term needs: it is the
            ! partner for the piece where the two first derivatives pair with
            ! each other.
            if (gga) then
               if (allocated(dgchi)) deallocate (dgchi)
               allocate (dgchi(nb, n_sig, 3))
               call density_times_ao_grad(ao_grad, d_sig(1:n_sig, 1:n_sig), nb, n_sig, &
                                          dgchi)
               if (unrestricted) then
                  if (allocated(dgchi_beta)) deallocate (dgchi_beta)
                  allocate (dgchi_beta(nb, n_sig, 3))
                  call density_times_ao_grad(ao_grad, db_sig(1:n_sig, 1:n_sig), nb, &
                                             n_sig, dgchi_beta)
               end if
            end if

            do ig = 1, nb
               gg = g0 + ig - 1
               own = ctx%grid%atom(gg)

               wv = ctx%grid%weights(gg)*vrho(gg)
               call accumulate_channel(ao_grad, dchi, ig, n_sig, wv, scale, &
                                       c_offsets, c_counts, natm, own, gradient)
               if (gga) then
                  call accumulate_gga_channel(ao_grad, ao_hess, dchi, dgchi, ig, n_sig, &
                                              ctx%grid%weights(gg)*gcoef(gg, :), scale, &
                                              c_offsets, c_counts, natm, own, gradient)
               end if
               if (mgga) then
                  call accumulate_mgga_channel(ao_grad, ao_hess, dgchi, ig, n_sig, &
                                               ctx%grid%weights(gg)*vtau(gg), scale, &
                                               c_offsets, c_counts, natm, own, gradient)
               end if
               if (unrestricted) then
                  wv = ctx%grid%weights(gg)*vrho_beta(gg)
                  call accumulate_channel(ao_grad, dchi_beta, ig, n_sig, wv, scale, &
                                          c_offsets, c_counts, natm, own, gradient)
                  if (gga) then
                     call accumulate_gga_channel(ao_grad, ao_hess, dchi_beta, &
                                                 dgchi_beta, ig, n_sig, &
                                                 ctx%grid%weights(gg)*gcoef_beta(gg, :), &
                                                 scale, c_offsets, c_counts, natm, own, &
                                                 gradient)
                  end if
               end if
            end do
         end if

         ! The partition-weight term. Differentiated blockwise for the same
         ! reason the density is: the full (3, natm, npoints) array is large and
         ! nothing needs it all at once.
         if (allocated(dpart)) deallocate (dpart)
         allocate (dpart(3, natm, nb))
         call becke_partition_derivatives(ctx%grid%coords(:, g0:g1), mol%coords, &
                                          ctx%grid%numbers, ctx%grid%atom(g0:g1), &
                                          ctx%grid%scheme, ctx%grid%adjust, dpart, error)
         if (error%has_error()) return

         do ig = 1, nb
            gg = g0 + ig - 1
            contrib(1) = ctx%grid%quad_weights(gg)*exc(gg)
            if (unrestricted) then
               contrib(1) = contrib(1)*(rho(gg) + rho_beta(gg))
            else
               contrib(1) = contrib(1)*rho(gg)
            end if
            do ia = 1, natm
               do comp = 1, 3
                  gradient(comp, ia) = gradient(comp, ia) + contrib(1)*dpart(comp, ia, ig)
               end do
            end do

            ! The partition weight of a point owned by A also changes because
            ! the *point* moves with A, not only because A moves. That second
            ! piece is what cancels the integrand's own grid-motion
            ! contribution above; leaving it out does not make the answer
            ! slightly worse but leaves that term uncancelled, which on water is
            ! an error of 0.77 Hartree/Bohr against a gradient of 0.06.
            !
            ! It needs no new derivative: displacing the point is the same as
            ! displacing every nucleus the other way, so the partition's
            ! gradient with respect to the point is minus the sum over atoms of
            ! what `becke_partition_derivatives` already returned.
            own = ctx%grid%atom(gg)
            do ia = 1, natm
               do comp = 1, 3
                  gradient(comp, own) = gradient(comp, own) &
                                        - contrib(1)*dpart(comp, ia, ig)
               end do
            end do
         end do
      end do
   end subroutine xc_gradient

   subroutine vv10_gradient(ctx, mol, density, gradient, error)
      !! VV10's contribution to dE/dR, accumulated in place
      type(xc_context_t), intent(inout) :: ctx
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(inout) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      call vv10_gradient_core(ctx, mol, density, gradient, error, .false.)
   end subroutine vv10_gradient

   subroutine vv10_gradient_fixed_grid(ctx, mol, density, gradient, error)
      !! VV10's dE/dR with the grid held fixed
      !!
      !! The first derivative a VV10 *Hessian* is the second derivative of, in
      !! the same sense `xc_gradient_fixed_grid` is for the semilocal term: the
      !! points do not move, the partition does not respond, and the kernel's
      !! dependence on where the points are drops out with them.
      !!
      !! **Not the physical gradient** -- `vv10_gradient` is, and the two differ
      !! by 4.5e-4 on water. This exists so a Hessian can be differenced against
      !! its own first derivative with the same approximation on both sides.
      type(xc_context_t), intent(inout) :: ctx
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(inout) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      call vv10_gradient_core(ctx, mol, density, gradient, error, .true.)
   end subroutine vv10_gradient_fixed_grid

   subroutine vv10_gradient_core(ctx, mol, density, gradient, error, fixed_grid)
      !! The two above, which differ only in whether the grid is allowed to move
      !!
      !! **The non-locality is already inside `vrho` and `vsigma`.** `E_nl` is a
      !! functional of rho and sigma, so
      !!
      !!     dE/dR = int dr [ dE/drho(r) drho(r)/dR + dE/dsigma(r) dsigma(r)/dR ]
      !!
      !! and the pair sum over the second grid is spent entirely on producing
      !! those two arrays. What is left is the ordinary GGA contraction, term
      !! for term, which is why this reads like `xc_gradient`'s GGA path with
      !! the quantities swapped.
      !!
      !! **The grid response is kept, unlike PySCF's**, which stops at the
      !! basis-function term: the semilocal gradient next door carries the
      !! moving points and moving partition weights, and a total made of one of
      !! each would be neither.
      !!
      !! Restricted only. The caller refuses an open shell before this is
      !! reached.
      type(xc_context_t), intent(inout) :: ctx
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(inout) :: gradient(:, :)
      type(error_t), intent(inout) :: error
      logical, intent(in) :: fixed_grid
         !! Skip the moving points and the responding partition, leaving the
         !! basis-function term alone.

      real(dp), allocatable :: rho(:), sigma(:), rho_grad(:, :)
      real(dp), allocatable :: exc(:), vrho(:), vsigma(:), gcoef(:, :)
      real(dp), allocatable :: dedw(:), fexp(:, :)
      real(dp), allocatable :: rho_blk(:), rho_grad_blk(:, :)
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), allocatable :: dchi(:, :), dgchi(:, :, :)
      real(dp), allocatable :: dpart(:, :, :)
      real(dp), allocatable :: extents(:), d_sig(:, :)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_list(:), ao_offset(:)
      integer, allocatable :: c_offsets(:), c_counts(:)
      real(dp) :: contrib, wv, scale
      integer :: npts, nao, natm, g0, g1, nb, ig, id, gg, comp, ia, own
      integer :: n_sig, isig, jsig

      if (ctx%nlc_b == 0.0_dp .and. ctx%nlc_c == 0.0_dp) return

      call ensure_nlc_grid(ctx, mol, error)
      if (error%has_error()) return
      npts = ctx%nlc_grid%n_points
      if (npts == 0) return

      nao = mol%nao
      natm = mol%natm
      scale = 2.0_dp

      ! Sweep one: rho and its gradient over the whole grid at once. The pair
      ! sum needs every point before it can produce any potential, so unlike the
      ! semilocal path this cannot be done a block at a time.
      allocate (rho(npts), sigma(npts), rho_grad(npts, 3))
      rho = 0.0_dp
      sigma = 0.0_dp
      rho_grad = 0.0_dp
      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1
         if (allocated(rho_blk)) deallocate (rho_blk, rho_grad_blk)
         allocate (rho_blk(nb), rho_grad_blk(nb, 3))
         call eval_ao_block(mol, ctx%nlc_grid%coords(:, g0:g1), ao, error, grad=ao_grad)
         if (error%has_error()) return
         call eval_rho(ao, density, rho_blk, ao_grad=ao_grad, rho_grad=rho_grad_blk)
         do ig = 1, nb
            rho(g0 + ig - 1) = rho_blk(ig)
            do id = 1, 3
               rho_grad(g0 + ig - 1, id) = rho_grad_blk(ig, id)
            end do
            sigma(g0 + ig - 1) = rho_grad_blk(ig, 1)**2 + rho_grad_blk(ig, 2)**2 &
                                 + rho_grad_blk(ig, 3)**2
         end do
      end do

      allocate (exc(npts), vrho(npts), vsigma(npts))
      if (fixed_grid) then
         allocate (dedw(0), fexp(3, 0))
         call vv10_nlc(ctx%nlc_b, ctx%nlc_c, ctx%nlc_grid%coords, rho, sigma, &
                       ctx%nlc_grid%coords, rho, sigma, ctx%nlc_grid%weights, &
                       exc, vrho, vsigma)
      else
         allocate (dedw(npts), fexp(3, npts))
         call vv10_nlc(ctx%nlc_b, ctx%nlc_c, ctx%nlc_grid%coords, rho, sigma, &
                       ctx%nlc_grid%coords, rho, sigma, ctx%nlc_grid%weights, &
                       exc, vrho, vsigma, dedw=dedw, fexp=fexp)
      end if

      ! dE/d(grad rho) = 2 vsigma grad rho, the chain rule the semilocal path
      ! applies to its own sigma, written the same way so the two cannot
      ! disagree about it.
      allocate (gcoef(npts, 3))
      do id = 1, 3
         do ig = 1, npts
            gcoef(ig, id) = 2.0_dp*vsigma(ig)*rho_grad(ig, id)
         end do
      end do

      allocate (c_offsets(natm), c_counts(natm))
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(nao))
      allocate (d_sig(nao, nao))
      call shell_extents(mol, ctx%screen_tol, extents)

      ! Sweep two: the three terms, on the grid the potential was built on.
      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         call block_significant_aos(mol, ctx%nlc_grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig, &
                                    atom_offsets=c_offsets, atom_counts=c_counts)

         if (n_sig > 0) then
            do jsig = 1, n_sig
               do isig = 1, n_sig
                  d_sig(isig, jsig) = density(ao_list(isig), ao_list(jsig))
               end do
            end do

            ! Second derivatives, because the sigma channel pairs a basis
            ! function's gradient with another derivative of it.
            call eval_ao_block(mol, ctx%nlc_grid%coords(:, g0:g1), ao, error, &
                               grad=ao_grad, hess=ao_hess, shell_mask=shell_mask, &
                               ao_offset=ao_offset, n_ao_out=n_sig)
            if (error%has_error()) return

            if (allocated(dchi)) deallocate (dchi)
            allocate (dchi(nb, n_sig))
            call density_times_ao(ao, d_sig(1:n_sig, 1:n_sig), nb, n_sig, dchi)

            if (allocated(dgchi)) deallocate (dgchi)
            allocate (dgchi(nb, n_sig, 3))
            call density_times_ao_grad(ao_grad, d_sig(1:n_sig, 1:n_sig), nb, n_sig, dgchi)

            do ig = 1, nb
               gg = g0 + ig - 1
               own = ctx%nlc_grid%atom(gg)
               wv = ctx%nlc_grid%weights(gg)*vrho(gg)
               call accumulate_channel(ao_grad, dchi, ig, n_sig, wv, scale, &
                                       c_offsets, c_counts, natm, own, gradient, &
                                       moving=.not. fixed_grid)
               call accumulate_gga_channel(ao_grad, ao_hess, dchi, dgchi, ig, n_sig, &
                                           ctx%nlc_grid%weights(gg)*gcoef(gg, :), scale, &
                                           c_offsets, c_counts, natm, own, gradient, &
                                           moving=.not. fixed_grid)
            end do
         end if

         if (fixed_grid) cycle

         ! The partition-weight term, and the point-motion term that cancels
         ! most of it, as `xc_gradient` writes them.
         !
         ! **Two things differ from the semilocal version and neither is
         ! optional.** The coefficient is `dE/dw` and not `rho*exc`, because a
         ! point's weight multiplies it twice over; and the points carry an
         ! explicit force of their own, because the kernel is a function of
         ! where they are. Together they are worth 4.5e-4 on water, and that
         ! does not shrink when the grid is refined.
         if (allocated(dpart)) deallocate (dpart)
         allocate (dpart(3, natm, nb))
         call becke_partition_derivatives(ctx%nlc_grid%coords(:, g0:g1), mol%coords, &
                                          ctx%nlc_grid%numbers, ctx%nlc_grid%atom(g0:g1), &
                                          ctx%nlc_grid%scheme, ctx%nlc_grid%adjust, &
                                          dpart, error)
         if (error%has_error()) return

         do ig = 1, nb
            gg = g0 + ig - 1
            own = ctx%nlc_grid%atom(gg)

            ! The point's own coordinate derivative, onto the atom it rides.
            do comp = 1, 3
               gradient(comp, own) = gradient(comp, own) &
                                     + ctx%nlc_grid%weights(gg)*fexp(comp, gg)
            end do

            contrib = ctx%nlc_grid%quad_weights(gg)*dedw(gg)
            do ia = 1, natm
               do comp = 1, 3
                  gradient(comp, ia) = gradient(comp, ia) + contrib*dpart(comp, ia, ig)
               end do
            end do
            do ia = 1, natm
               do comp = 1, 3
                  gradient(comp, own) = gradient(comp, own) &
                                        - contrib*dpart(comp, ia, ig)
               end do
            end do
         end do
      end do
   end subroutine vv10_gradient_core

   subroutine xc_potential_gradient(ctx, mol, density, pmat, gradient, error, fixed_grid)
      !! `d/dR Tr(P V_xc[D])`, with both densities held fixed, accumulated in place
      !!
      !! What a double hybrid's gradient needs and no gradient before it did:
      !! its perturbative term is not stationary, so eliminating the orbital
      !! response leaves the Z-vector contracted with the derivative of the
      !! *reference* operator, which for a Kohn-Sham reference contains `V_xc`.
      !! Omitting it gives a gradient of entirely plausible magnitude, short by
      !! a piece the same order as the Coulomb one beside it.
      !!
      !! **`P` is not a density in the usual sense.** It is the relaxed
      !! correlation density plus the Z-vector: symmetric but indefinite,
      !! integrating to zero rather than to the electron count, and nowhere near
      !! idempotent. The functional is evaluated at `D` alone and `P` appears
      !! only linearly, which is what makes this the derivative of a *linear
      !! form* rather than of an energy.
      !!
      !!     G(R) = int w(R) [ v_rho rho_P + 2 v_sigma grad rho . grad rho_P ]
      !!
      !! `dG/dR` has the same three sources the exchange-correlation gradient
      !! has -- the basis functions move, the grid points move, the partition
      !! weights move -- plus one: `v_rho` and `v_sigma` are functions of the
      !! reference density, which moves too. That is where the second
      !! derivatives enter.
      !!
      !! **Four channels, two of them the reference's**, collected by which
      !! object is differentiated:
      !!
      !!     d rho      : f_rr rho_P + 2 f_rs (grad rho . grad rho_P)
      !!     d grad rho : 2 f_rs rho_P grad rho
      !!                  + 4 f_ss (grad rho . grad rho_P) grad rho
      !!                  + 2 v_sigma grad rho_P
      !!     d rho_P    : v_rho
      !!     d grad rho_P: 2 v_sigma grad rho
      !!
      !! Each is a contraction `accumulate_channel` and `accumulate_gga_channel`
      !! already perform for the energy gradient; only the density they carry
      !! and the coefficient weighting it differ.
      !!
      !! LDA and GGA; a meta-GGA is refused upstream by
      !! `xc_grid_kernel_quantities`.
      type(xc_context_t), intent(inout) :: ctx
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)   !! The converged reference density
      real(dp), intent(in) :: pmat(:, :)      !! The density `V_xc` is traced against
      real(dp), intent(inout) :: gradient(:, :)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: fixed_grid
         !! Hold the quadrature fixed: omit both the partition-weight response
         !! and the grid points travelling with their owning atom.
         !!
         !! Absent or false gives the physical derivative, which is what every
         !! production caller wants. True gives the quantity a fixed-grid Hessian
         !! is the second derivative *of*, so the two can be differenced against
         !! each other with the same approximation on both sides --
         !! `xc_gradient_fixed_grid` is the same idea for the energy term and
         !! `vv10_gradient_fixed_grid` for the non-local one.
         !!
         !! Differencing the physical form against a fixed-grid Hessian instead
         !! disagrees by exactly the omitted term: 5.5e-04 on b3lyp/cc-pVDZ water
         !! at grid level 3, falling only to 6.5e-05 by level 5. Large enough to
         !! swamp a real error and small enough to be mistaken for one.

      real(dp), allocatable :: rho(:), rho_grad(:, :), vrho(:), vsigma(:)
      real(dp), allocatable :: frr(:), frs(:), fss(:)
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), allocatable :: dchi(:, :), dgchi(:, :, :)
      real(dp), allocatable :: pchi(:, :), pgchi(:, :, :)
      real(dp), allocatable :: rho_p(:), p_grad(:, :)
      real(dp), allocatable :: dpart(:, :, :)
      real(dp) :: wg_ref(3), wg_p(3)
      real(dp) :: gdotp, coef_rho, integrand, w
      real(dp), allocatable :: extents(:), d_sig(:, :), p_sig(:, :)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_list(:), ao_offset(:)
      integer, allocatable :: c_offsets(:), c_counts(:)
      integer :: npts, nao, natm, g0, g1, nb, ig, gg, id, ia, comp, own
      integer :: n_sig, isig, jsig
      logical :: gga
      logical :: hold_grid

      if (.not. ctx%active) return

      nao = mol%nao
      natm = mol%natm
      npts = ctx%grid%n_points
      gga = ctx%any_gga

      call xc_grid_kernel_quantities(ctx, mol, density, rho, rho_grad, vrho, &
                                     vsigma, frr, frs, fss, error)
      if (error%has_error()) return

      hold_grid = .false.
      if (present(fixed_grid)) hold_grid = fixed_grid

      allocate (c_offsets(natm), c_counts(natm))
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(nao))
      allocate (d_sig(nao, nao), p_sig(nao, nao))
      call shell_extents(mol, ctx%screen_tol, extents)

      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         ! Which functions reach this block. The screen pays for more here than
         ! in the energy gradient, because `eval_ao_block` is asked for second
         ! derivatives as well: ten components per function per point.
         call block_significant_aos(mol, ctx%grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig, &
                                    atom_offsets=c_offsets, atom_counts=c_counts)

         ! Allocated and zeroed *outside* the screening guard, because the
         ! partition-weight term further down reads them at every point of the
         ! block. In a block no basis function reaches they are legitimately
         ! zero, and what they must not be is whatever the previous block left
         ! in them -- which would also be the wrong length.
         if (allocated(rho_p)) deallocate (rho_p, p_grad)
         allocate (rho_p(nb), p_grad(nb, 3))
         rho_p = 0.0_dp
         p_grad = 0.0_dp

         if (n_sig > 0) then
            do jsig = 1, n_sig
               do isig = 1, n_sig
                  d_sig(isig, jsig) = density(ao_list(isig), ao_list(jsig))
                  p_sig(isig, jsig) = pmat(ao_list(isig), ao_list(jsig))
               end do
            end do

            ! Second derivatives of the basis functions, for the same reason the
            ! GGA energy gradient needs them: differentiating `grad chi` with
            ! respect to a nuclear coordinate is a second derivative.
            if (gga) then
               call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, &
                                  grad=ao_grad, hess=ao_hess, shell_mask=shell_mask, &
                                  ao_offset=ao_offset, n_ao_out=n_sig)
            else
               call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, &
                                  grad=ao_grad, shell_mask=shell_mask, &
                                  ao_offset=ao_offset, n_ao_out=n_sig)
            end if
            if (error%has_error()) return

            if (allocated(dchi)) deallocate (dchi, pchi)
            allocate (dchi(nb, n_sig), pchi(nb, n_sig))
            call density_times_ao(ao, d_sig(1:n_sig, 1:n_sig), nb, n_sig, dchi)
            call density_times_ao(ao, p_sig(1:n_sig, 1:n_sig), nb, n_sig, pchi)

            ! rho_P and its gradient, on the same footing as the reference's.
            ! Built from `pchi` directly rather than through `eval_rho`, since
            ! the partner matrices are wanted anyway.
            do ig = 1, nb
               rho_p(ig) = sum(pchi(ig, :)*ao(ig, :))
            end do
            if (gga) then
               if (allocated(dgchi)) deallocate (dgchi, pgchi)
               allocate (dgchi(nb, n_sig, 3), pgchi(nb, n_sig, 3))
               call density_times_ao_grad(ao_grad, d_sig(1:n_sig, 1:n_sig), nb, &
                                          n_sig, dgchi)
               call density_times_ao_grad(ao_grad, p_sig(1:n_sig, 1:n_sig), nb, &
                                          n_sig, pgchi)
               ! grad rho_P = 2 sum_uv P_uv chi_v grad chi_u -- the two is the
               ! symmetry of P, exactly as in `eval_rho`.
               do id = 1, 3
                  do ig = 1, nb
                     p_grad(ig, id) = 2.0_dp*sum(pchi(ig, :)*ao_grad(ig, :, id))
                  end do
               end do
            end if

            do ig = 1, nb
               gg = g0 + ig - 1
               own = ctx%grid%atom(gg)
               w = ctx%grid%weights(gg)

               gdotp = 0.0_dp
               if (gga) gdotp = sum(rho_grad(gg, :)*p_grad(ig, :))

               ! The reference density's channel: what multiplies d rho / dR.
               coef_rho = frr(gg)*rho_p(ig)
               if (gga) coef_rho = coef_rho + 2.0_dp*frs(gg)*gdotp
               call accumulate_channel(ao_grad, dchi, ig, n_sig, w*coef_rho, 2.0_dp, &
                                       c_offsets, c_counts, natm, own, gradient, &
                                       moving=.not. hold_grid)

               ! `P`'s own channel, weighted by the ordinary potential: the
               ! term that survives when the functional is a constant.
               call accumulate_channel(ao_grad, pchi, ig, n_sig, w*vrho(gg), 2.0_dp, &
                                       c_offsets, c_counts, natm, own, gradient, &
                                       moving=.not. hold_grid)

               if (gga) then
                  do id = 1, 3
                     wg_ref(id) = w*(2.0_dp*frs(gg)*rho_p(ig)*rho_grad(gg, id) &
                                     + 4.0_dp*fss(gg)*gdotp*rho_grad(gg, id) &
                                     + 2.0_dp*vsigma(gg)*p_grad(ig, id))
                     wg_p(id) = w*2.0_dp*vsigma(gg)*rho_grad(gg, id)
                  end do
                  call accumulate_gga_channel(ao_grad, ao_hess, dchi, dgchi, ig, &
                                              n_sig, wg_ref, 2.0_dp, c_offsets, &
                                              c_counts, natm, own, gradient, &
                                              moving=.not. hold_grid)
                  call accumulate_gga_channel(ao_grad, ao_hess, pchi, pgchi, ig, &
                                              n_sig, wg_p, 2.0_dp, c_offsets, &
                                              c_counts, natm, own, gradient, &
                                              moving=.not. hold_grid)
               end if
            end do
         end if

         ! The partition weights, on the integrand of this linear form rather
         ! than on the energy density. Same two pieces as the energy gradient:
         ! every nucleus reweights the whole grid, and a point owned by A also
         ! moves with A.
         ! The partition weights are one of the two halves of the grid
         ! response. The other is the points moving with their owning atom,
         ! which lives in the `moving=` forwards above rather than here --
         ! believing this block was the whole of it is what shipped a flag that
         ! held half the grid.
         if (hold_grid) cycle
         if (allocated(dpart)) deallocate (dpart)
         allocate (dpart(3, natm, nb))
         call becke_partition_derivatives(ctx%grid%coords(:, g0:g1), mol%coords, &
                                          ctx%grid%numbers, ctx%grid%atom(g0:g1), &
                                          ctx%grid%scheme, ctx%grid%adjust, dpart, error)
         if (error%has_error()) return

         do ig = 1, nb
            gg = g0 + ig - 1
            integrand = vrho(gg)*rho_p(ig)
            if (gga) integrand = integrand &
                                 + 2.0_dp*vsigma(gg)*sum(rho_grad(gg, :)*p_grad(ig, :))
            integrand = ctx%grid%quad_weights(gg)*integrand

            own = ctx%grid%atom(gg)
            do ia = 1, natm
               do comp = 1, 3
                  gradient(comp, ia) = gradient(comp, ia) + integrand*dpart(comp, ia, ig)
                  gradient(comp, own) = gradient(comp, own) - integrand*dpart(comp, ia, ig)
               end do
            end do
         end do
      end do
   end subroutine xc_potential_gradient

   subroutine fitted_reference_gradient(mol, aux, three, jm12, dm_a, dm_b, &
                                        gradient, error, k_scale)
      !! `d/dR Tr(G[A] B)` when `G` is built from fitted integrals
      !!
      !! The fitted counterpart of `two_electron_mp2_terms`. Fitted,
      !! `(mu nu|lam sig)` is a product of three-centre integrals and a metric,
      !! so its derivative has no four-centre term anywhere in it and this
      !! returns the two intermediates that contract against the three- and
      !! two-centre ones instead.
      !!
      !! **Coulomb.** `F_J = sum_PQ g^A_P J^-1_PQ g^B_Q` with
      !! `g^A_P = sum_uv A_uv (uv|P)`. Each `g` differentiates, and
      !! `d(J^-1) = -J^-1 (dJ) J^-1` supplies the metric term:
      !!
      !!     Gamma^P_uv  +=  A_uv rho^B_P + B_uv rho^A_P,   rho = J^-1 g
      !!     Omega_PQ    += -rho^A_P rho^B_Q
      !!
      !! **Exchange.** `F_K = -(k/2) sum_PQ J^-1_PQ Tr(B M^P A M^Q)`, writing
      !! `M^P_uv = (uv|P)`. Both `M` factors differentiate, and the two results
      !! are *transposes* of one another rather than equal -- which is the trap,
      !! because either choice leaves `Gamma` symmetric and only the assembled
      !! gradient can tell them apart. With `Y^P = sum_Q J^-1_PQ M^Q` and
      !! `S^P = B Y^P A`,
      !!
      !!     Gamma^P     += -(k/2) (S^P + S^P transposed)
      !!     Omega_RS    += +(k/2) sum_uv S^R_uv Y^S_uv
      !!
      !! The metric half is written through `S` and `Y` rather than through the
      !! `n_aux^2` intermediate `W_PQ = Tr(B M^P A M^Q)` it comes from, since
      !! `S` already carries one of the two applications of `J^-1`.
      !!
      !! **The caller halves this.** What an assembly wants is
      !! `(1/2) Tr(G^x[D] D) + Tr(G^x[D] P)`, which is half of
      !! `Tr(G^x[D](D + 2P))`. The exact path reaches the same half by building
      !! only two of the four differentiated positions, invisibly from its call
      !! site; this returns the honest derivative, so the factor is the caller's
      !! and can be read there.
      type(czt_molecule_t), intent(in) :: mol, aux
      real(dp), intent(in) :: three(:, :)   !! Raw `(mu nu|P)`, (n_ao^2, naux)
      real(dp), intent(in) :: jm12(:, :)    !! `J^(-1/2)`, as the energy fitted with
      real(dp), intent(in) :: dm_a(:, :)    !! The density inside the operator
      real(dp), intent(in) :: dm_b(:, :)    !! The density it is contracted against
      real(dp), intent(inout) :: gradient(:, :)   !! (3, natm), accumulated into
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: k_scale

      real(dp), allocatable :: jinv(:, :), y(:, :), s(:, :), gamma(:, :, :)
      real(dp), allocatable :: omega(:, :), gt(:, :, :)
      real(dp), allocatable :: ip1(:, :, :, :), ip2(:, :, :, :), j1(:, :, :)
      real(dp), allocatable :: rho_a(:), rho_b(:), work(:, :)
      integer, allocatable :: offsets(:), counts(:)
      integer, allocatable :: aux_offsets(:), aux_counts(:)
      integer :: n, naux, p, iatom, comp, p0, p1, q0, q1
      real(dp) :: kf

      if (error%has_error()) return

      ! The auxiliary loops below deposit into `gradient(:, iatom)` using the
      ! auxiliary basis's own atom numbering, so the two bases must sit on the
      ! same atoms in the same order. `df_two_electron_gradient` refuses the
      ! same mismatch for the same reason: a derivative lands on the wrong
      ! nucleus, which translational invariance catches only after the fact --
      ! and if the auxiliary basis has the *more* atoms, the deposit runs off
      ! the end of `gradient` entirely.
      if (aux%natm /= mol%natm) then
         call error%set(ERROR_VALIDATION, &
                        "density-fitted reference gradient: the auxiliary basis is "// &
                        "on a different number of atoms than the orbital basis")
         return
      end if
      ! Equal counts are not enough: two same-sized molecules with their atoms in
      ! a different order pass that and still put every auxiliary derivative on
      ! the wrong nucleus. Same comparison `merge_basis_sets` makes.
      if (maxval(abs(aux%coords - mol%coords)) > 1.0e-10_dp) then
         call error%set(ERROR_VALIDATION, &
                        "density-fitted reference gradient: the auxiliary basis is "// &
                        "on a different geometry than the orbital basis")
         return
      end if

      kf = 1.0_dp
      if (present(k_scale)) kf = k_scale

      n = mol%nao
      naux = aux%nao

      allocate (jinv(naux, naux))
      call pic_gemm(jm12, jm12, jinv)

      ! `Y^P = sum_Q J^-1_PQ M^Q`, wanted by both halves of the exchange term
      ! and by the two `rho` vectors, which are `<Y^P, D>`.
      allocate (y(n*n, naux))
      call pic_gemm(three, jinv, y)

      allocate (rho_a(naux), rho_b(naux))
      rho_a = matmul(reshape(dm_a, [n*n]), y)
      rho_b = matmul(reshape(dm_b, [n*n]), y)

      ! `S^P = B Y^P A`, one pair of n^3 products per auxiliary function.
      allocate (s(n*n, naux), work(n, n))
      !$omp parallel do default(none) shared(s, y, dm_a, dm_b, n, naux) private(p, work)
      do p = 1, naux
         work = matmul(dm_b, reshape(y(:, p), [n, n]))
         s(:, p) = reshape(matmul(work, dm_a), [n*n])
      end do
      !$omp end parallel do
      deallocate (work)

      allocate (gamma(n, n, naux))
      do p = 1, naux
         associate (sp => reshape(s(:, p), [n, n]))
            gamma(:, :, p) = dm_a*rho_b(p) + dm_b*rho_a(p) &
                             - 0.5_dp*kf*(sp + transpose(sp))
         end associate
      end do

      allocate (omega(naux, naux))
      call pic_gemm(s, y, omega, transa="T")
      omega = 0.5_dp*kf*omega
      do p = 1, naux
         omega(:, p) = omega(:, p) - rho_a*rho_b(p)
      end do
      omega = 0.5_dp*(omega + transpose(omega))

      ! ---- contraction, in the same conventions the correlation's own fitted
      ! terms use: libcint's `ip` integrals differentiate the electronic
      ! coordinate, so the derivative with respect to a nucleus carrying the
      ! function is minus them.
      ! TODO(mqc): the auxiliary loops below index `gradient` by the auxiliary
      ! basis's own atom numbering, which only matches the orbital basis's when
      ! the two sit on the same atoms in the same order. `df_two_electron_gradient`
      ! refuses a mismatch explicitly; this does not, and would attribute a
      ! derivative to the wrong nucleus.
      allocate (offsets(mol%natm), counts(mol%natm))
      allocate (aux_offsets(aux%natm), aux_counts(aux%natm))
      call atom_ao_blocks(mol, offsets, counts)
      call atom_ao_blocks(aux, aux_offsets, aux_counts)

      call three_centre_deriv(mol, aux, 1, ip1)
      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)
         if (counts(iatom) == 0) cycle
         ! `(mu nu|P)` is symmetric in mu and nu, so the ket derivative is the
         ! bra one against a transposed density rather than a second integral.
         gt = transposed_block(gamma, p0, p1)
         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    - sum(ip1(p0:p1, :, :, comp)*gamma(p0:p1, :, :)) &
                                    - sum(ip1(p0:p1, :, :, comp)*gt)
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
                                    - sum(ip2(:, :, q0:q1, comp)*gamma(:, :, q0:q1))
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
                                    - sum(j1(q0:q1, :, comp)*omega(q0:q1, :)) &
                                    - sum(j1(q0:q1, :, comp)*transpose(omega(:, q0:q1)))
         end do
      end do
      deallocate (j1)
   end subroutine fitted_reference_gradient

   function transposed_block(g, p0, p1) result(t)
      !! `g(nu, mu, P)` for `mu` in a block, as the ket-side contraction needs
      real(dp), intent(in) :: g(:, :, :)
      integer, intent(in) :: p0, p1
      real(dp), allocatable :: t(:, :, :)
      integer :: p

      allocate (t(p1 - p0 + 1, size(g, 1), size(g, 3)))
      do p = 1, size(g, 3)
         t(:, :, p) = transpose(g(:, p0:p1, p))
      end do
   end function transposed_block

   subroutine density_times_ao(ao, density, nb, nao, dchi)
      !! (D chi)_mu(g) = sum_nu D_mu,nu chi_nu(g)
      real(dp), intent(in) :: ao(:, :)
      real(dp), intent(in) :: density(:, :)
      integer, intent(in) :: nb, nao
      real(dp), intent(out) :: dchi(:, :)

      ! `dchi = ao D^T`. `D` is symmetric, so the transpose is free -- written
      ! anyway, so that this does not quietly depend on symmetry the day
      ! somebody passes an energy-weighted matrix that has none.
      call pic_gemm(ao, density, dchi, transb="T")
   end subroutine density_times_ao

   subroutine accumulate_channel(ao_grad, dchi, ig, nao, wv, scale, &
                                 offsets, counts, natm, own, gradient, moving)
      !! One spin channel's basis-function and grid-point terms at one point
      !!
      !! Both are the same contraction of `grad chi` against `D chi`, differing
      !! in which basis functions are summed and where the result lands: the
      !! basis-function term takes only the functions on atom A and lands on A,
      !! the grid-point term takes every function and lands on the atom owning
      !! the point. Their signs are opposite, which is what makes the two cancel
      !! when every atom moves together.
      real(dp), intent(in) :: ao_grad(:, :, :)
      real(dp), intent(in) :: dchi(:, :)
      integer, intent(in) :: ig, nao, natm, own
      real(dp), intent(in) :: wv, scale
      integer, intent(in) :: offsets(:), counts(:)
      real(dp), intent(inout) :: gradient(:, :)
      logical, intent(in), optional :: moving
         !! Whether the grid points travel with their owning atom. Default true,
         !! which is the physical gradient. False leaves the basis-function term
         !! alone, which is the *fixed-grid* first derivative a Hessian omitting
         !! the grid response is the second derivative of.

      integer :: ia, mu, comp, p0, p1
      real(dp) :: block_sum(3), total(3)
      logical :: points_move

      points_move = .true.
      if (present(moving)) points_move = moving

      total = 0.0_dp
      do ia = 1, natm
         p0 = offsets(ia) + 1
         p1 = offsets(ia) + counts(ia)
         block_sum = 0.0_dp
         do mu = p0, p1
            do comp = 1, 3
               block_sum(comp) = block_sum(comp) + ao_grad(ig, mu, comp)*dchi(ig, mu)
            end do
         end do
         do comp = 1, 3
            ! Moving atom A moves its functions; d chi / dR_A is minus the
            ! gradient with respect to the electron coordinate.
            gradient(comp, ia) = gradient(comp, ia) - scale*wv*block_sum(comp)
            total(comp) = total(comp) + block_sum(comp)
         end do
      end do

      ! The point travels with the atom that produced it, which is the same
      ! contraction over every function and with the opposite sign.
      if (points_move) then
         do comp = 1, 3
            gradient(comp, own) = gradient(comp, own) + scale*wv*total(comp)
         end do
      end if
   end subroutine accumulate_channel

   subroutine density_times_ao_grad(ao_grad, density, nb, nao, dgchi)
      !! (D grad chi)_mu,k(g) = sum_nu D_mu,nu d chi_nu / dx_k (g)
      !!
      !! The gradient counterpart of `density_times_ao`. Only the GGA term
      !! needs it, and only for the piece in which both derivatives are first
      !! derivatives -- one on the bra and one on the ket.
      real(dp), intent(in) :: ao_grad(:, :, :)
      real(dp), intent(in) :: density(:, :)
      integer, intent(in) :: nb, nao
      real(dp), intent(out) :: dgchi(:, :, :)

      integer :: k

      ! One gemm per Cartesian direction. `ao_grad(:, :, k)` and its target are
      ! contiguous, the direction being the last dimension, so each is a plain
      ! `(nb, nao) x (nao, nao)` and no copy is made. Written out, the innermost
      ! index walks `ao_grad(ig, nu, k)` along `nu`, striding by `nb` -- one
      ! useful element per cache line, and the single most expensive thing in a
      ! density functional gradient run.
      do k = 1, 3
         call pic_gemm(ao_grad(:, :, k), density, dgchi(:, :, k), transb="T")
      end do
   end subroutine density_times_ao_grad

   subroutine accumulate_mgga_channel(ao_grad, ao_hess, dgchi, ig, nao, wt, scale, &
                                      offsets, counts, natm, own, gradient)
      !! The kinetic-energy-density term at one grid point
      !!
      !! A meta-GGA reads `tau` as well, and moving a nucleus moves the
      !! functions it is built from:
      !!
      !!     tau = 1/2 sum_j sum_uv D_uv d_j chi_u d_j chi_v
      !!     dtau/dR_{A,k} = -sum_{mu in A} sum_j d2chi_mu/dx_k dx_j (D grad chi)_mu,j
      !!
      !! **The half is already gone by the time it is written that way.** `tau`
      !! carries libxc's 1/2, and differentiating pairs the moving function with
      !! either the bra or the ket, which are equal by the symmetry of D -- so
      !! the two cancel and the net coefficient is one, not one half and not
      !! two. It is applied through `scale` and an explicit half below, where
      !! the matching factors in `accumulate_xc_matrix` appear.
      !!
      !! Only one contraction, unlike the GGA term next door: `tau` pairs a
      !! first derivative with a first derivative, so differentiating it gives a
      !! second derivative against `(D grad chi)` and nothing against `chi`.
      real(dp), intent(in) :: ao_grad(:, :, :)
      real(dp), intent(in) :: ao_hess(:, :, :)
      real(dp), intent(in) :: dgchi(:, :, :)
      integer, intent(in) :: ig, nao, natm, own
      real(dp), intent(in) :: wt          !! weight times dE/dtau
      real(dp), intent(in) :: scale
      integer, intent(in) :: offsets(:), counts(:)
      real(dp), intent(inout) :: gradient(:, :)

      ! Which packed component holds d2/dx_k dx_j, as `accumulate_gga_channel`.
      integer, parameter :: HESS_INDEX(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

      integer :: ia, mu, j, k, p0, p1
      real(dp) :: block_sum(3)
      real(dp) :: total(3)
      real(dp) :: acc, factor

      factor = 0.5_dp*scale*wt
      total = 0.0_dp
      do ia = 1, natm
         p0 = offsets(ia) + 1
         p1 = offsets(ia) + counts(ia)
         block_sum = 0.0_dp
         do mu = p0, p1
            do k = 1, 3
               acc = 0.0_dp
               do j = 1, 3
                  acc = acc + ao_hess(ig, mu, HESS_INDEX(k, j))*dgchi(ig, mu, j)
               end do
               block_sum(k) = block_sum(k) + acc
            end do
         end do
         do k = 1, 3
            gradient(k, ia) = gradient(k, ia) - factor*block_sum(k)
            total(k) = total(k) + block_sum(k)
         end do
      end do

      ! The point travels with the atom that owns it, exactly as the other two.
      do k = 1, 3
         gradient(k, own) = gradient(k, own) + factor*total(k)
      end do
      ! TODO(mqc): dead statement, and `nao` is an unused dummy here as it is in
      ! `accumulate_channel` and `accumulate_gga_channel`. Drop the argument
      ! from all three rather than referencing it to keep a warning quiet.
      if (nao < 0) return
   end subroutine accumulate_mgga_channel

   subroutine accumulate_gga_channel(ao_grad, ao_hess, dchi, dgchi, ig, nao, wg, scale, &
                                     offsets, counts, natm, own, gradient, moving)
      !! One spin channel's density-gradient term at one grid point
      !!
      !! A GGA reads `grad rho` as well as `rho`, so moving a nucleus moves that
      !! too, and the chain rule needs
      !!
      !!     d(grad rho)_j / dR_{A,k}
      !!         = -2 sum_{mu in A, nu} D_mu,nu [ d2chi_mu/dx_k dx_j * chi_nu
      !!                                          + dchi_mu/dx_j * dchi_nu/dx_k ]
      !!
      !! **The two pieces are not one term doubled, and that is the trap.** The
      !! second derivative pairs with `chi_nu`; the two first derivatives pair
      !! with *each other*, one on the bra and one on the ket. A single doubled
      !! term gives a gradient of entirely plausible magnitude that misses
      !! finite differences, so the check that decides this routine is the
      !! numerical one and not translational invariance.
      !!
      !! `wg` is the weight times dE/d(grad rho), which
      !! `xc_grid_gga_quantities` returns already resolved over libxc's sigma
      !! channels, so nothing here knows whether the calculation was polarised.
      !!
      !! The second-derivative index follows `AO_HESS_COMP`'s packing: xx, xy,
      !! xz, yy, yz, zz, so element (k, j) of the symmetric matrix is looked up
      !! rather than stored twice.
      real(dp), intent(in) :: ao_grad(:, :, :)
      real(dp), intent(in) :: ao_hess(:, :, :)
      real(dp), intent(in) :: dchi(:, :)
      real(dp), intent(in) :: dgchi(:, :, :)
      integer, intent(in) :: ig, nao, natm, own
      logical, intent(in), optional :: moving
         !! As `accumulate_channel`: false leaves the basis-function term alone.
      real(dp), intent(in) :: wg(3)
      real(dp), intent(in) :: scale
      integer, intent(in) :: offsets(:), counts(:)
      real(dp), intent(inout) :: gradient(:, :)

      ! Which packed component holds d2/dx_k dx_j.
      integer, parameter :: HESS_INDEX(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

      integer :: ia, mu, j, k, p0, p1
      real(dp) :: block_sum(3)
      real(dp) :: total(3)
      logical :: points_move
      real(dp) :: acc

      points_move = .true.
      if (present(moving)) points_move = moving

      total = 0.0_dp
      do ia = 1, natm
         p0 = offsets(ia) + 1
         p1 = offsets(ia) + counts(ia)
         block_sum = 0.0_dp
         do mu = p0, p1
            do k = 1, 3
               acc = 0.0_dp
               do j = 1, 3
                  ! Note which index sits where. The derivative direction k is
                  ! on the *bra* in both pieces; the summed direction j rides
                  ! the second derivative in the first piece and the *ket* in
                  ! the second. Pairing j with the bra instead is wrong by about
                  ! six percent on water/PBE, and translational invariance does
                  ! not notice.
                  acc = acc + wg(j)*(ao_hess(ig, mu, HESS_INDEX(k, j))*dchi(ig, mu) &
                                     + ao_grad(ig, mu, k)*dgchi(ig, mu, j))
               end do
               block_sum(k) = block_sum(k) + acc
            end do
         end do
         do k = 1, 3
            gradient(k, ia) = gradient(k, ia) - scale*block_sum(k)
            total(k) = total(k) + block_sum(k)
         end do
      end do

      ! The point travels with the atom that owns it: the same contraction over
      ! every basis function, with the opposite sign, exactly as the LDA term
      ! next door.
      do k = 1, 3
         if (points_move) gradient(k, own) = gradient(k, own) + scale*total(k)
      end do
   end subroutine accumulate_gga_channel

   subroutine nuclear_repulsion_gradient(mol, gradient)
      !! dE_NN/dR, accumulated into `gradient`
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(inout) :: gradient(:, :)

      call add_nuclear_repulsion_gradient(mol%charges, mol%coords, gradient)
   end subroutine nuclear_repulsion_gradient

   subroutine energy_weighted_density(orbitals, energies, n_occ, closed_shell, weighted)
      !! W = sum_i occ_i eps_i C_i C_i^T
      !!
      !! The Pulay term contracts against this rather than the density, because
      !! what moves with the basis functions is the orthonormality constraint,
      !! whose Lagrange multipliers are the orbital energies.
      ! TODO(mqc): written as an n_ao^2 * n_occ triple loop where every other
      ! contraction in this module goes through `pic_gemm`.
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      logical, intent(in) :: closed_shell   !! Occupation two rather than one
      real(dp), allocatable, intent(out) :: weighted(:, :)

      integer :: nao, mu, nu, i
      real(dp) :: occupation, acc

      nao = size(orbitals, 1)
      allocate (weighted(nao, nao))
      weighted = 0.0_dp

      occupation = 1.0_dp
      if (closed_shell) occupation = 2.0_dp

      do nu = 1, nao
         do mu = 1, nao
            acc = 0.0_dp
            do i = 1, n_occ
               acc = acc + occupation*energies(i)*orbitals(mu, i)*orbitals(nu, i)
            end do
            weighted(mu, nu) = acc
         end do
      end do
   end subroutine energy_weighted_density

   subroutine one_electron_deriv(mol, matrix, which)
      !! A one-electron derivative matrix, as (nao, nao, 3)
      !!
      !! The same shell-pair loop as `one_electron` next door, with three
      !! components per block instead of one. libcint returns them with the
      !! component slowest, so a block is `buf(di, dj, 1:3)`.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), allocatable, intent(out) :: matrix(:, :, :)
      integer, intent(in) :: which

      real(dp), allocatable :: buf(:)
      integer :: shls(2)
      integer :: ish, jsh, di, dj, i, j, comp, ret, io, jo, mx

      allocate (matrix(mol%nao, mol%nao, 3))
      matrix = 0.0_dp

      ! Flat, and indexed by hand below. libcint packs a block with the *actual*
      ! shell dimensions, so a rank-3 buffer declared at the largest shell has
      ! the wrong strides for every smaller one -- invisible in a basis of s
      ! functions and wrong everywhere else.
      mx = max_block(mol)

      ! Threaded over the bra shell: every block a thread writes sits in its
      ! own rows, so nothing is shared but the molecule.
      !$omp parallel default(none) shared(mol, matrix, which, mx) &
      !$omp    private(ish, jsh, di, dj, io, jo, i, j, comp, ret, shls, buf)
      allocate (buf(mx*mx*3))
      !$omp do schedule(dynamic)
      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, mol%nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]

            if (mol%cartesian) then
               select case (which)
               case (DERIV_OVLP)
                  ret = libcint_1e_ipovlp_cart(buf, shls, mol%atm, mol%natm, &
                                               mol%bas, mol%nbas, mol%env)
               case (DERIV_KIN)
                  ret = libcint_1e_ipkin_cart(buf, shls, mol%atm, mol%natm, &
                                              mol%bas, mol%nbas, mol%env)
               case default
                  ret = libcint_1e_ipnuc_cart(buf, shls, mol%atm, mol%natm, &
                                              mol%bas, mol%nbas, mol%env)
               end select
            else
               select case (which)
               case (DERIV_OVLP)
                  ret = libcint_1e_ipovlp_sph(buf, shls, mol%atm, mol%natm, &
                                              mol%bas, mol%nbas, mol%env)
               case (DERIV_KIN)
                  ret = libcint_1e_ipkin_sph(buf, shls, mol%atm, mol%natm, &
                                             mol%bas, mol%nbas, mol%env)
               case default
                  ret = libcint_1e_ipnuc_sph(buf, shls, mol%atm, mol%natm, &
                                             mol%bas, mol%nbas, mol%env)
               end select
            end if

            if (ret == 0) cycle
            do comp = 1, 3
               do j = 1, dj
                  do i = 1, di
                     matrix(io + i, jo + j, comp) = buf(i + di*(j - 1 + dj*(comp - 1)))
                  end do
               end do
            end do
         end do
      end do
      !$omp end do
      deallocate (buf)
      !$omp end parallel
   end subroutine one_electron_deriv

   subroutine iprinv_deriv_at(mol, iatom, matrix)
      !! The 1/|r-R| derivative with the operator origin on atom `iatom`
      !!
      !! `env` is copied rather than modified, because the origin lives in it
      !! and the molecule is shared.
      type(czt_molecule_t), intent(in) :: mol
      integer, intent(in) :: iatom
      real(dp), intent(out) :: matrix(:, :, :)

      real(dp), allocatable :: env(:), buf(:)
      integer :: shls(2)
      integer :: ish, jsh, di, dj, i, j, comp, ret, io, jo, mx

      allocate (env(size(mol%env)), source=mol%env)
      env(LIBCINT_PTR_RINV_ORIG + 1:LIBCINT_PTR_RINV_ORIG + 3) = mol%coords(:, iatom)

      matrix = 0.0_dp
      mx = max_block(mol)

      !$omp parallel default(none) shared(mol, matrix, env, mx) &
      !$omp    private(ish, jsh, di, dj, io, jo, i, j, comp, ret, shls, buf)
      allocate (buf(mx*mx*3))
      !$omp do schedule(dynamic)
      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, mol%nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]

            if (mol%cartesian) then
               ret = libcint_1e_iprinv_cart(buf, shls, mol%atm, mol%natm, &
                                            mol%bas, mol%nbas, env)
            else
               ret = libcint_1e_iprinv_sph(buf, shls, mol%atm, mol%natm, &
                                           mol%bas, mol%nbas, env)
            end if

            if (ret == 0) cycle
            do comp = 1, 3
               do j = 1, dj
                  do i = 1, di
                     matrix(io + i, jo + j, comp) = buf(i + di*(j - 1 + dj*(comp - 1)))
                  end do
               end do
            end do
         end do
      end do
      !$omp end do
      deallocate (buf)
      !$omp end parallel
   end subroutine iprinv_deriv_at

   subroutine hellmann_feynman_gradient(mol, density, gradient)
      !! The Hellmann-Feynman term of dE/dR, added into `gradient`
      !!
      !! `-2 Z_A sum_ij D_ij (nabla i|1/|r-R_A||j)` for every atom A: the force
      !! of the electrons on a nucleus that moves while the basis stays put.
      !! Integration by parts turns the derivative of the operator into the
      !! derivative of the bra plus that of the ket, and `D` is symmetric, so
      !! the ket's half is the bra's counted twice.
      !!
      !! One pass over the shell pairs with the atoms inside, rather than a
      !! full `iprinv_deriv_at` matrix per atom: the integral count is the
      !! same and nothing of size `nao^2` is stored.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)       !! Symmetric, the total density
      real(dp), intent(inout) :: gradient(:, :)   !! (3, natm), added to

      real(dp), allocatable :: env(:), buf(:), g_local(:, :)
      real(dp) :: acc
      integer :: shls(2)
      integer :: ish, jsh, di, dj, i, j, comp, ret, io, jo, mx, iatom

      mx = max_block(mol)

      ! `env` is private because the operator origin lives in it and each
      ! thread moves it from atom to atom.
      !$omp parallel default(none) shared(mol, density, gradient, mx) &
      !$omp    private(ish, jsh, di, dj, io, jo, i, j, comp, ret, shls, buf, env, &
      !$omp            g_local, iatom, acc)
      allocate (buf(mx*mx*3), env(size(mol%env)), g_local(3, mol%natm))
      env = mol%env
      g_local = 0.0_dp
      !$omp do schedule(dynamic)
      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, mol%nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]

            do iatom = 1, mol%natm
               env(LIBCINT_PTR_RINV_ORIG + 1:LIBCINT_PTR_RINV_ORIG + 3) = mol%coords(:, iatom)
               if (mol%cartesian) then
                  ret = libcint_1e_iprinv_cart(buf, shls, mol%atm, mol%natm, &
                                               mol%bas, mol%nbas, env)
               else
                  ret = libcint_1e_iprinv_sph(buf, shls, mol%atm, mol%natm, &
                                              mol%bas, mol%nbas, env)
               end if
               if (ret == 0) cycle

               do comp = 1, 3
                  acc = 0.0_dp
                  do j = 1, dj
                     do i = 1, di
                        acc = acc + buf(i + di*(j - 1 + dj*(comp - 1)))*density(io + i, jo + j)
                     end do
                  end do
                  g_local(comp, iatom) = g_local(comp, iatom) - 2.0_dp*mol%charges(iatom)*acc
               end do
            end do
         end do
      end do
      !$omp end do
      !$omp critical
      gradient = gradient + g_local
      !$omp end critical
      deallocate (buf, env, g_local)
      !$omp end parallel
   end subroutine hellmann_feynman_gradient

   subroutine two_electron_deriv(mol, density, vhf, error, exchange_density, exx_fraction, &
                                 screen_tol, omega, with_coulomb)
      !! `J - K/2` built from the differentiated ERIs, as (nao, nao, 3)
      !!
      !! `density` drives the Coulomb field and is the total density in both the
      !! restricted and unrestricted cases. `exchange_density` is the one
      !! exchange is built from: absent for a closed shell, where the total
      !! density carries its factor of two and `J - K/2` is the right
      !! combination; present for one spin of an unrestricted reference, where
      !! exchange takes that spin alone and the coefficient is one, not a half.
      !!
      !! **One permutation survives, and it is used.** A differentiated quartet
      !! has none of the eightfold symmetry of an ordinary one -- the nabla
      !! distinguishes the first index from the second and the bra from the ket.
      !! What is left is the ket pair, `(∇i j|k l) = (∇i j|l k)`, so this walks
      !! `l <= k` and counts the off-diagonal pair twice: once as itself and
      !! once transposed, which is a different contraction and is written out
      !! rather than doubled.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: vhf(:, :, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: exchange_density(:, :)
      real(dp), intent(in), optional :: exx_fraction   !! Default one, for Hartree-Fock
      real(dp), intent(in), optional :: screen_tol
         !! Bound below which a shell quartet is skipped, on the *contribution*
         !! rather than the integral. Default `DERIV_SCREEN_TOL`.
      real(dp), intent(in), optional :: omega
         !! Switch the kernel to `erf(omega r)/r`, for the long-range exchange
         !! of a range-separated functional. It travels through `env` rather
         !! than a separate entry point, so an attenuated pass is this same loop
         !! over a modified copy.
      logical, intent(in), optional :: with_coulomb
         !! Default true. False builds exchange alone, which is what the
         !! long-range pass wants: `J` is complete from the full-range pass,
         !! and the attenuated kernel's Coulomb term belongs to no energy.

      real(dp), allocatable :: buf(:), vj(:, :, :), vk(:, :, :), env_use(:)
      real(dp), allocatable :: vj_local(:, :, :), vk_local(:, :, :)
      real(dp), allocatable :: exchange_from(:, :)
      real(dp), allocatable :: bounds(:, :), bq(:, :), bra_bound(:, :), dsh(:, :), esh(:, :)
      integer, allocatable :: dims(:), offs(:)
      real(dp) :: g, k_scale, tol, bra, est, amax, jx
      logical :: do_j
      type(c_ptr) :: opt
      type(eri_shell_table_t) :: tab
      integer :: shls(4)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: io, jo, ko, lo, i, j, k, l, comp, ret, mx, idx, nao, nbas
      integer :: nprim, ptr

      ! The shells the quartet loop runs over: the fused-sp view when the
      ! molecule carries one, its split shells otherwise -- the same table the
      ! direct Fock build takes. The AO order of a fused shell is its split
      ! shells' order, so the scatter into vj/vk needs nothing beyond the view's
      ! own offsets. libfint's `l_shell_grad_check` guards the stride an L
      ! shell's s and p coefficients are applied on, which is what makes the
      ! view usable for a multi-component integral at all.
      call eri_shell_table(mol, tab)
      mx = tab%block_max
      nao = mol%nao
      nbas = tab%nbas
      allocate (vj(nao, nao, 3), vk(nao, nao, 3))
      vj = 0.0_dp
      vk = 0.0_dp

      tol = DERIV_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol

      do_j = .true.
      if (present(with_coulomb)) do_j = with_coulomb
      jx = 0.0_dp
      if (do_j) jx = 1.0_dp

      ! A copy, because `tab%env` is shared and an attenuated pass must not
      ! leave the omega set behind for the next caller. The slot is one-based
      ! here and zero-based in libfint's headers, hence the `+ 1`; getting it
      ! wrong is silent, and the "attenuated" integrals come back full-range.
      env_use = tab%env
      if (present(omega)) env_use(LIBCINT_PTR_RANGE_OMEGA + 1) = omega

      allocate (exchange_from(nao, nao))
      if (present(exchange_density)) then
         exchange_from = exchange_density
         k_scale = 1.0_dp
      else
         exchange_from = density
         k_scale = 0.5_dp
      end if

      if (present(exx_fraction)) k_scale = k_scale*exx_fraction

      opt = c_null_ptr
      if (mol%cartesian) then
         call libcint_2e_ip1_cart_optimizer(opt, mol%atm, mol%natm, tab%bas, tab%nbas, env_use)
      else
         call libcint_2e_ip1_sph_optimizer(opt, mol%atm, mol%natm, tab%bas, tab%nbas, env_use)
      end if

      ! Shell dimensions and offsets once, rather than a `shell_dim` call at
      ! every level of a four-deep loop that runs millions of times.
      dims = tab%dims
      offs = tab%offs(1:nbas)

      ! Screening. Two factors go into the estimate:
      !
      !   * **The bra carries a derivative, so its Schwarz bound is not the
      !     ordinary one.** Differentiating a primitive brings down `2*alpha*r`
      !     and the extent of that primitive is `r ~ 1/sqrt(alpha)`, so the
      !     distribution grows by roughly `2*sqrt(alpha)`. Taking the largest
      !     exponent in the shell makes that an over-estimate, which is the
      !     direction a screen has to err in. It is **not** the rigorous
      !     Cauchy-Schwarz bound, which would need `(∇i j|∇i j)`, so the
      !     tolerance is two decades tighter than the accuracy wanted and the
      !     result is checked against an unscreened gradient in the tests.
      !   * **The density the quartet multiplies**, blocked to shells as
      !     `shell_density_max` does it next door.
      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) return
      ! Per split shell from `schwarz_bounds`, re-blocked onto the view's
      ! shells -- an identity copy when there is no view.
      call eri_schwarz_collapse(mol, bounds, bq)

      allocate (bra_bound(nbas, nbas))
      do ish = 1, nbas
         nprim = tab%bas(LIBCINT_NPRIM_OF, ish)
         ptr = tab%bas(LIBCINT_PTR_EXP, ish)
         amax = maxval(tab%env(ptr + 1:ptr + nprim))
         bra_bound(ish, :) = 2.0_dp*sqrt(amax)*bq(ish, :)
      end do

      allocate (dsh(nbas, nbas), esh(nbas, nbas))
      do ish = 1, nbas
         do jsh = 1, nbas
            dsh(ish, jsh) = maxval(abs(density(offs(ish) + 1:offs(ish) + dims(ish), &
                                               offs(jsh) + 1:offs(jsh) + dims(jsh))))
            esh(ish, jsh) = maxval(abs(exchange_from(offs(ish) + 1:offs(ish) + dims(ish), &
                                                     offs(jsh) + 1:offs(jsh) + dims(jsh))))
         end do
      end do

      ! Thread-local accumulators merged once, rather than atomics on a shared
      ! matrix: the inner update is two scattered writes per integral. The price
      ! is `2 * nao^2 * 3` doubles per thread, worth watching on a few thousand
      ! functions.
      !
      ! schedule(dynamic) because the `l <= k` triangle makes the work per `ish`
      ! uneven.
      !$omp parallel default(none) &
      !$omp    shared(mol, tab, density, exchange_from, opt, vj, vk, mx, nao, nbas, &
      !$omp           dims, offs, bq, bra_bound, dsh, esh, tol, env_use, do_j) &
      !$omp    private(ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo, &
      !$omp            i, j, k, l, comp, ret, idx, g, shls, buf, vj_local, vk_local, &
      !$omp            bra, est)
      allocate (buf(mx**4*3))
      allocate (vj_local(nao, nao, 3), vk_local(nao, nao, 3))
      vj_local = 0.0_dp
      vk_local = 0.0_dp

      !$omp do schedule(dynamic)
      do ish = 1, nbas
         di = dims(ish)
         io = offs(ish)
         do jsh = 1, nbas
            dj = dims(jsh)
            jo = offs(jsh)
            bra = bra_bound(ish, jsh)
            do ksh = 1, nbas
               dk = dims(ksh)
               ko = offs(ksh)
               do lsh = 1, ksh
                  dl = dims(lsh)
                  lo = offs(lsh)

                  ! The largest density element this quartet can reach: two
                  ! Coulomb, from the ket pair both ways round, and two
                  ! exchange, pairing the bra's second index with each of the
                  ! ket's -- exactly the four the inner loop reads.
                  if (do_j) then
                     est = bra*bq(ksh, lsh) &
                           *max(dsh(lsh, ksh), dsh(ksh, lsh), esh(jsh, ksh), esh(jsh, lsh))
                  else
                     est = bra*bq(ksh, lsh)*max(esh(jsh, ksh), esh(jsh, lsh))
                  end if
                  if (est < tol) cycle

                  shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]

                  if (mol%cartesian) then
                     ret = libcint_2e_ip1_cart(buf, shls, mol%atm, mol%natm, &
                                               tab%bas, nbas, env_use, opt)
                  else
                     ret = libcint_2e_ip1_sph(buf, shls, mol%atm, mol%natm, &
                                              tab%bas, nbas, env_use, opt)
                  end if
                  if (ret == 0) cycle

                  do comp = 1, 3
                     do l = 1, dl
                        do k = 1, dk
                           do j = 1, dj
                              do i = 1, di
                                 idx = i + di*(j - 1 + dj*(k - 1 + dk*(l - 1 + dl*(comp - 1))))
                                 g = buf(idx)
                                 vj_local(io + i, jo + j, comp) = vj_local(io + i, jo + j, comp) &
                                                                  + g*density(lo + l, ko + k)
                                 vk_local(io + i, lo + l, comp) = vk_local(io + i, lo + l, comp) &
                                                                  + g*exchange_from(jo + j, ko + k)
                                 ! The transposed ket pair. Skipped on the
                                 ! diagonal, where it is the same quartet.
                                 if (lsh /= ksh) then
                                    vj_local(io + i, jo + j, comp) = vj_local(io + i, jo + j, comp) &
                                                                     + g*density(ko + k, lo + l)
                                    vk_local(io + i, ko + k, comp) = vk_local(io + i, ko + k, comp) &
                                                                     + g*exchange_from(jo + j, lo + l)
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
      vj = vj + vj_local
      vk = vk + vk_local
      !$omp end critical
      deallocate (buf, vj_local, vk_local)
      !$omp end parallel

      call libcint_del_optimizer(opt)

      ! Minus, for the same reason the one-electron ones carry it: libcint
      ! returns the nabla on the bra, and the derivative with respect to the
      ! atom is its negative.
      allocate (vhf(nao, nao, 3))
      vhf = -(jx*vj - k_scale*vk)

      deallocate (bounds, bq, bra_bound, dsh, esh, dims, offs)
   end subroutine two_electron_deriv

   subroutine two_electron_gradient(mol, density, gradient, error, density_alpha, density_beta, &
                                    exx_fraction, screen_tol, omega, with_coulomb)
      !! The two-electron term of dE/dR, added into `gradient`
      !!
      !! `density` is the total density and drives the Coulomb term. Exchange
      !! comes from the two spin densities when both are present, and from the
      !! total alone when neither is -- a closed shell, where the total carries
      !! its factor of two and each spin is half of it.
      !!
      !! **Unique quartets, three integrals each.** A differentiated quartet
      !! `(nabla i j|kl)` has none of the eightfold symmetry of `(ij|kl)`, but
      !! the density it is contracted against has all of it, so the loop is the
      !! direct Fock build's -- unique shell quartets, each weighted by its
      !! degeneracy -- and inside one the derivative on each of the first three
      !! centres is `int2e_ip1` on a permuted quartet: `(nabla i j|kl)`,
      !! `(nabla j i|kl)` and `(nabla k l|ij)`. The fourth is translational
      !! invariance -- the four derivatives of an integral sum to zero -- so it
      !! is counted rather than computed. That is three integrals in eight where
      !! `two_electron_deriv`, which walks every ordered bra pair, takes four,
      !! and the result is the gradient itself rather than an `nao^2` matrix
      !! per component and thread.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)          !! The total density
      real(dp), intent(inout) :: gradient(:, :)      !! (3, natm), added to
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: density_alpha(:, :)
      real(dp), intent(in), optional :: density_beta(:, :)
         !! Both or neither. One spin alone is refused.
      real(dp), intent(in), optional :: exx_fraction   !! Default one, for Hartree-Fock
      real(dp), intent(in), optional :: screen_tol
         !! Bound below which one derivative integral's contribution to the
         !! gradient is skipped. Default `DERIV_SCREEN_TOL`.
      real(dp), intent(in), optional :: omega
         !! Switch the kernel to `erf(omega r)/r`, for the long-range exchange
         !! of a range-separated functional.
      logical, intent(in), optional :: with_coulomb
         !! Default true. False builds exchange alone, for the long-range pass.

      real(dp), allocatable :: da(:, :), db(:, :), env_use(:)
      real(dp), allocatable :: bounds(:, :), bq(:, :), dfac(:), dsh(:, :), dsa(:, :), dsb(:, :)
      real(dp), allocatable :: buf1(:), buf2(:), buf3(:), gam(:), gam2(:), gam3(:)
      real(dp), allocatable :: g_local(:, :)
      integer, allocatable :: dims(:), offs(:), atom_of(:), pair_i(:), pair_j(:), order(:)
      real(dp) :: f(3, 3), kx, jx, tol, coef, deg, schwarz, gmax, weight, amax
      logical :: do_j, unrestricted, have1, have2, have3, same_bra, same_pair
      type(c_ptr) :: opt
      type(eri_shell_table_t) :: tab
      integer :: shls(4)
      integer :: itask, ij, kl, npair, ipair, nbas, mx, n, nprim, ptr
      integer :: s1, s2, s3, s4, d1, d2, d3, d4, o1, o2, o3, o4
      integer :: a, b, c, d, ia, ib, ic, id, idx, comp, ret

      unrestricted = present(density_alpha)
      if (unrestricted .neqv. present(density_beta)) then
         call error%set(ERROR_VALIDATION, "two_electron_gradient: the spin densities "// &
                        "come as a pair; one alone is not an unrestricted reference")
         return
      end if

      call eri_shell_table(mol, tab)
      mx = tab%block_max
      nbas = tab%nbas
      dims = tab%dims
      offs = tab%offs(1:nbas)

      tol = GRAD_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol

      do_j = .true.
      if (present(with_coulomb)) do_j = with_coulomb
      jx = 0.0_dp
      if (do_j) jx = 1.0_dp

      ! The exchange part of the two-particle density, symmetrised over the
      ! quartet: `-(1/2) sum_spin (Ds_ik Ds_jl + Ds_il Ds_jk)`, times the
      ! fraction of exact exchange. A closed shell is the same expression with
      ! each spin half the total, which is why the total is split rather than
      ! given a coefficient of its own.
      kx = 0.5_dp
      if (present(exx_fraction)) kx = 0.5_dp*exx_fraction
      if (unrestricted) then
         da = density_alpha
         db = density_beta
      else
         da = 0.5_dp*density
         db = da
      end if

      ! A copy, because `tab%env` is shared and an attenuated pass must not
      ! leave the omega set behind for the next caller. The slot is one-based
      ! here and zero-based in libfint's headers, hence the `+ 1`; getting it
      ! wrong is silent, and the "attenuated" integrals come back full-range.
      env_use = tab%env
      if (present(omega)) env_use(LIBCINT_PTR_RANGE_OMEGA + 1) = omega

      opt = c_null_ptr
      if (mol%cartesian) then
         call libcint_2e_ip1_cart_optimizer(opt, mol%atm, mol%natm, tab%bas, tab%nbas, env_use)
      else
         call libcint_2e_ip1_sph_optimizer(opt, mol%atm, mol%natm, tab%bas, tab%nbas, env_use)
      end if

      ! Screening, as `two_electron_deriv` does it: the Schwarz bound of the
      ! undifferentiated quartet, grown by `2*sqrt(alpha_max)` for the shell
      ! the nabla lands on, times the largest two-particle density the block
      ! can hold. Not the rigorous Cauchy-Schwarz bound of the derivative
      ! distribution -- that was measured on a 74-atom 6-31G* case and removed
      ! three per cent of the quartets this one keeps, not worth its pass -- so
      ! the tolerance stays two decades inside the accuracy wanted, and the
      ! tests check the result against an unscreened gradient.
      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) return
      call eri_schwarz_collapse(mol, bounds, bq)

      allocate (dfac(nbas), atom_of(nbas))
      do s1 = 1, nbas
         nprim = tab%bas(LIBCINT_NPRIM_OF, s1)
         ptr = tab%bas(LIBCINT_PTR_EXP, s1)
         amax = maxval(tab%env(ptr + 1:ptr + nprim))
         dfac(s1) = 2.0_dp*sqrt(amax)
         atom_of(s1) = tab%bas(LIBCINT_ATOM_OF, s1) + 1
      end do

      call block_density_max(density, nbas, offs, dims, dsh)
      call block_density_max(da, nbas, offs, dims, dsa)
      call block_density_max(db, nbas, offs, dims, dsb)

      ! The quartet loop flattened onto shell pairs, most expensive first, so
      ! `schedule(dynamic)` does not finish with every other thread idle behind
      ! the last d-shell pair.
      npair = nbas*(nbas + 1)/2
      allocate (pair_i(npair), pair_j(npair))
      ipair = 0
      do s1 = 1, nbas
         do s2 = 1, s1
            ipair = ipair + 1
            pair_i(ipair) = s1
            pair_j(ipair) = s2
         end do
      end do
      call pair_work_order(pair_i, pair_j, dims, order)

      ! Each thread accumulates a (3, natm) of its own and adds it in once:
      ! the whole of what a quartet touches is four atoms.
      !$omp parallel default(none) &
      !$omp    shared(mol, tab, density, da, db, opt, mx, nbas, &
      !$omp           dims, offs, atom_of, bq, dfac, dsh, dsa, dsb, tol, env_use, jx, kx, &
      !$omp           pair_i, pair_j, order, npair, gradient) &
      !$omp    private(itask, ij, kl, s1, s2, s3, s4, d1, d2, d3, d4, o1, o2, o3, o4, &
      !$omp            a, b, c, d, ia, ib, ic, id, idx, comp, ret, n, shls, deg, schwarz, &
      !$omp            gmax, weight, coef, f, have1, have2, have3, same_bra, same_pair, &
      !$omp            buf1, buf2, buf3, gam, gam2, gam3, g_local)
      allocate (buf1(mx**4*3), buf2(mx**4*3), buf3(mx**4*3))
      allocate (gam(mx**4), gam2(mx**4), gam3(mx**4))
      allocate (g_local(3, mol%natm))
      g_local = 0.0_dp

      !$omp do schedule(dynamic)
      do itask = 1, npair
         ij = order(itask)
         s1 = pair_i(ij)
         s2 = pair_j(ij)
         d1 = dims(s1)
         o1 = offs(s1)
         d2 = dims(s2)
         o2 = offs(s2)
         same_bra = s1 == s2

         do kl = 1, ij
            s3 = pair_i(kl)
            s4 = pair_j(kl)
            d3 = dims(s3)
            o3 = offs(s3)
            d4 = dims(s4)
            o4 = offs(s4)
            same_pair = kl == ij

            ! The largest |Gamma| in the block, from the same density maxima
            ! the elements below are built from. `deg/2` is the weight the
            ! contribution carries, so it belongs in the bound.
            deg = pair_degeneracy(s1, s2, s3, s4)
            schwarz = bq(s1, s2)*bq(s3, s4)
            gmax = jx*dsh(s1, s2)*dsh(s3, s4) &
                   + kx*(dsa(s1, s3)*dsa(s2, s4) + dsa(s1, s4)*dsa(s2, s3) &
                         + dsb(s1, s3)*dsb(s2, s4) + dsb(s1, s4)*dsb(s2, s3))
            weight = 0.5_dp*deg*schwarz*gmax

            ! One decision per integral. A skipped integral is zero on its own
            ! centre and in the fourth centre's sum alike, so the three are
            ! independent. When two shells coincide the second integral is the
            ! first one read with its indices swapped, and needs no call.
            have1 = weight*dfac(s1) >= tol
            have2 = weight*dfac(s2) >= tol
            have3 = weight*dfac(s3) >= tol
            if (same_bra) have2 = have1
            if (same_pair) have3 = have1
            if (.not. (have1 .or. have2 .or. have3)) cycle

            n = d1*d2*d3*d4

            if (have1) then
               shls = [s1 - 1, s2 - 1, s3 - 1, s4 - 1]
               have1 = ip1_block(mol%cartesian, buf1, shls, mol%atm, mol%natm, &
                                 tab%bas, nbas, env_use, opt)
               if (same_bra) have2 = have1
               if (same_pair) have3 = have1
            end if
            if (have2 .and. .not. same_bra) then
               shls = [s2 - 1, s1 - 1, s3 - 1, s4 - 1]
               have2 = ip1_block(mol%cartesian, buf2, shls, mol%atm, mol%natm, &
                                 tab%bas, nbas, env_use, opt)
            end if
            if (have3 .and. .not. same_pair) then
               shls = [s3 - 1, s4 - 1, s1 - 1, s2 - 1]
               have3 = ip1_block(mol%cartesian, buf3, shls, mol%atm, mol%natm, &
                                 tab%bas, nbas, env_use, opt)
            end if
            if (.not. (have1 .or. have2 .or. have3)) cycle

            ! Gamma once per element, laid out three ways to match the three
            ! buffers: (a,b,c,d) for the nabla on a, (b,a,c,d) on b, (c,d,a,b)
            ! on c. The contraction is then a dot product per component.
            do d = 1, d4
               id = o4 + d
               do c = 1, d3
                  ic = o3 + c
                  do b = 1, d2
                     ib = o2 + b
                     do a = 1, d1
                        ia = o1 + a
                        idx = a + d1*(b - 1 + d2*(c - 1 + d3*(d - 1)))
                        gam(idx) = jx*density(ia, ib)*density(ic, id) &
                                   - kx*(da(ia, ic)*da(ib, id) + da(ia, id)*da(ib, ic) &
                                         + db(ia, ic)*db(ib, id) + db(ia, id)*db(ib, ic))
                        gam2(b + d2*(a - 1 + d1*(c - 1 + d3*(d - 1)))) = gam(idx)
                        gam3(c + d3*(d - 1 + d4*(a - 1 + d1*(b - 1)))) = gam(idx)
                     end do
                  end do
               end do
            end do

            f = 0.0_dp
            do comp = 1, 3
               if (have1) f(comp, 1) = dot_product(gam(1:n), buf1((comp - 1)*n + 1:comp*n))
               if (have2) then
                  if (same_bra) then
                     f(comp, 2) = dot_product(gam2(1:n), buf1((comp - 1)*n + 1:comp*n))
                  else
                     f(comp, 2) = dot_product(gam2(1:n), buf2((comp - 1)*n + 1:comp*n))
                  end if
               end if
               if (have3) then
                  if (same_pair) then
                     f(comp, 3) = dot_product(gam3(1:n), buf1((comp - 1)*n + 1:comp*n))
                  else
                     f(comp, 3) = dot_product(gam3(1:n), buf3((comp - 1)*n + 1:comp*n))
                  end if
               end if
            end do

            ! Minus, because libcint puts the nabla on the function and the
            ! derivative with respect to the atom it sits on is the negative;
            ! the fourth centre is minus the other three again.
            coef = 0.5_dp*deg
            g_local(:, atom_of(s1)) = g_local(:, atom_of(s1)) - coef*f(:, 1)
            g_local(:, atom_of(s2)) = g_local(:, atom_of(s2)) - coef*f(:, 2)
            g_local(:, atom_of(s3)) = g_local(:, atom_of(s3)) - coef*f(:, 3)
            g_local(:, atom_of(s4)) = g_local(:, atom_of(s4)) + coef*(f(:, 1) + f(:, 2) + f(:, 3))
         end do
      end do
      !$omp end do

      !$omp critical
      gradient = gradient + g_local
      !$omp end critical
      deallocate (buf1, buf2, buf3, gam, gam2, gam3, g_local)
      !$omp end parallel

      call libcint_del_optimizer(opt)
   end subroutine two_electron_gradient

   function ip1_block(cartesian, buf, shls, atm, natm, bas, nbas, env, opt) result(have)
      !! One `int2e_ip1` shell quartet; false when libcint found it zero
      logical, intent(in) :: cartesian
      real(dp), intent(inout), contiguous :: buf(:)
      integer, intent(in) :: shls(4)
      integer, intent(in), contiguous :: atm(:, :)
      integer, intent(in) :: natm
      integer, intent(in), contiguous :: bas(:, :)
      integer, intent(in) :: nbas
      real(dp), intent(in), contiguous :: env(:)
      type(c_ptr), intent(in) :: opt
      logical :: have

      integer :: ret

      if (cartesian) then
         ret = libcint_2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt)
      else
         ret = libcint_2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt)
      end if
      have = ret /= 0
   end function ip1_block

   subroutine df_two_electron_gradient(orb, aux, total_density, orbitals, n_occupied, &
                                       orbitals_beta, n_occupied_beta, gradient, error, &
                                       exx_fraction, rs_omega, with_coulomb)
      !! The two-electron gradient when J and K are density-fitted
      !!
      !! **What is being differentiated.** The fitted energy is
      !!
      !!    E_2e = (1/2) g^T J^-1 g  -  (exchange, below)
      !!
      !! with `g_P = sum_uv (P|uv) D_uv` and `J_PQ = (P|Q)`. Only the integrals
      !! depend on the geometry -- the SCF is stationary in its orbitals, so
      !! their derivatives are already accounted for by the Pulay term -- and
      !! `d(J^-1)/dx = -J^-1 J^x J^-1` turns the Coulomb derivative into
      !!
      !!    rho^T g^x  -  (1/2) rho^T J^x rho,      rho = J^-1 g
      !!
      !! Exchange has the same two shapes over the fully transformed
      !! `e^P_ij = sum_ul C_ui C_lj (ul|P)` and `f = J^-1 e`.
      !!
      !! **Everything collapses into two densities**, so the expensive
      !! derivative integrals are each contracted once:
      !!
      !!    Gamma^P_uv = rho_P D_uv - sum_spin c_s Z^{s,P}_uv
      !!    Omega_PQ   = -(1/2) rho_P rho_Q + sum_spin (c_s/2) sum_ij f^s_ij f^s_ij
      !!
      !! with `Z^{s,P}_uv = sum_ij C_ui f^{s,P}_ij C_vj`, and `c_s` the channel
      !! weight: two for the single closed-shell channel, one for each of the
      !! two unrestricted ones. That the ratio between the two terms is the
      !! same either way is why one routine covers both.
      !!
      !! **The pseudo-inverse must be the SCF's.** `J^-1` is formed as
      !! `J^-1/2 J^-1/2` from `metric_inverse_sqrt`, the same routine and the
      !! same eigenvalue threshold `build_df_tensor` uses. A fitting set is
      !! near-linearly-dependent by construction, so a `J^-1` keeping different
      !! modes would differentiate a slightly different energy from the one that
      !! converged.
      type(czt_molecule_t), intent(in) :: orb, aux
      real(dp), intent(in) :: total_density(:, :)
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_occupied
      real(dp), intent(in), optional :: orbitals_beta(:, :)
      integer, intent(in), optional :: n_occupied_beta
      real(dp), intent(inout) :: gradient(:, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: exx_fraction
         !! Fraction of exact exchange the fitted energy kept. Absent is one,
         !! which is Hartree-Fock; a hybrid passes its own and a pure
         !! functional passes zero.
      real(dp), intent(in), optional :: rs_omega
         !! Fit and differentiate the attenuated kernel `erf(omega r)/r`. Named
         !! `rs_omega` and not `omega` because `omega` is already the two-index
         !! density `Omega_PQ` in this routine.
         !!
         !! A range-separated functional calls this routine twice: once
         !! full-range for `J` and its short-range exchange fraction, and once
         !! attenuated with `with_coulomb` false for the long-range exchange.
      logical, intent(in), optional :: with_coulomb
         !! False drops the Coulomb terms -- both the `rho^T g^x` shape and the
         !! `-(1/2) rho^T J^x rho` one -- leaving exchange alone. Default true.
         !!
         !! It has to drop *both*. Keeping the metric term while dropping the
         !! other leaves a piece of a Coulomb energy that the second pass never
         !! contained, and it is the same size as the answer.

      real(dp), allocatable :: three(:, :), metric(:, :), half(:, :), jinv(:, :)
      real(dp), allocatable :: g(:), rho(:)
      real(dp), allocatable :: gamma(:, :, :), omega(:, :)
      real(dp), allocatable :: ip1(:, :, :, :), ip2(:, :, :, :), d2c(:, :, :)
      integer, allocatable :: orb_off(:), orb_cnt(:), aux_off(:), aux_cnt(:)
      integer :: nao, naux, iatom, comp, p0, p1, q0, q1, ip
      logical :: unrestricted, coulomb
      real(dp) :: kf

      nao = orb%nao
      naux = aux%nao
      unrestricted = present(orbitals_beta)

      ! The per-atom accumulation below indexes both bases by the same atom
      ! number, so a mismatch would attribute an auxiliary function's derivative
      ! to the wrong nucleus.
      if (aux%natm /= orb%natm) then
         call error%set(ERROR_VALIDATION, &
                        "density-fitted gradient: the auxiliary basis is on a "// &
                        "different set of atoms than the orbital basis")
         return
      end if

      coulomb = .true.
      if (present(with_coulomb)) coulomb = with_coulomb

      if (present(rs_omega)) then
         call three_centre(orb, aux, three, omega=rs_omega)
         call two_centre(aux, metric, omega=rs_omega)
      else
         call three_centre(orb, aux, three)
         call two_centre(aux, metric)
      end if
      call metric_inverse_sqrt(metric, half, error)
      if (error%has_error()) return

      ! J^-1 from the square root, for the reason in the header.
      allocate (jinv(naux, naux))
      jinv = 0.0_dp
      call pic_gemm(half, half, jinv)

      allocate (g(naux), rho(naux))
      do ip = 1, naux
         g(ip) = sum(reshape(three(:, ip), [nao, nao])*total_density)
      end do
      rho = matmul(jinv, g)
      if (.not. coulomb) rho = 0.0_dp

      allocate (gamma(nao, nao, naux), omega(naux, naux))
      do ip = 1, naux
         gamma(:, :, ip) = rho(ip)*total_density
      end do
      ! The exchange fraction rides the channel weight rather than being applied
      ! afterwards: both of exchange's shapes -- the one in `gamma` and the one
      ! in `omega` -- are linear in it, and `add_exchange_channel` builds them
      ! together.
      kf = 1.0_dp
      if (present(exx_fraction)) kf = exx_fraction

      ! Both Coulomb shapes go together when this is an exchange-only pass:
      ! `rho` is left zero above, so the `rho^T g^x` term vanishes with it, and
      ! this metric term has to be dropped here to match.
      if (coulomb) then
         omega = -0.5_dp*outer(rho, rho)
      else
         omega = 0.0_dp
      end if

      ! Skipped outright at zero rather than multiplied by it: a pure functional
      ! would otherwise pay the whole `e^P_ij` transform, the `J^-1 e` solve and
      ! the `Z^P` back-transform to add nothing.
      if (kf /= 0.0_dp) then
         if (unrestricted) then
            call add_exchange_channel(three, jinv, orbitals, n_occupied, kf, gamma, omega)
            call add_exchange_channel(three, jinv, orbitals_beta, n_occupied_beta, kf, &
                                      gamma, omega)
         else
            call add_exchange_channel(three, jinv, orbitals, n_occupied, 2.0_dp*kf, &
                                      gamma, omega)
         end if
      end if

      deallocate (three, metric, half, jinv, g, rho)

      if (present(rs_omega)) then
         call three_centre_deriv(orb, aux, 1, ip1, omega=rs_omega)
         call three_centre_deriv(orb, aux, 2, ip2, omega=rs_omega)
         call two_centre_deriv(aux, d2c, omega=rs_omega)
      else
         call three_centre_deriv(orb, aux, 1, ip1)
         call three_centre_deriv(orb, aux, 2, ip2)
         call two_centre_deriv(aux, d2c)
      end if

      allocate (orb_off(orb%natm), orb_cnt(orb%natm))
      allocate (aux_off(aux%natm), aux_cnt(aux%natm))
      call atom_ao_blocks(orb, orb_off, orb_cnt)
      call atom_ao_blocks(aux, aux_off, aux_cnt)

      do iatom = 1, orb%natm
         p0 = orb_off(iatom) + 1
         p1 = orb_off(iatom) + orb_cnt(iatom)
         q0 = aux_off(iatom) + 1
         q1 = aux_off(iatom) + aux_cnt(iatom)

         do comp = 1, 3
            ! Orbital indices, twice: `ip1` differentiates the first only, and
            ! the second contributes the transpose -- the same number, since
            ! both (uv|P) and Gamma^P are symmetric in u and v. Minus, because
            ! libcint's nabla is on the bra.
            if (orb_cnt(iatom) > 0) then
               gradient(comp, iatom) = gradient(comp, iatom) &
                                       - 2.0_dp*sum(ip1(p0:p1, :, :, comp)*gamma(p0:p1, :, :))
            end if

            ! The auxiliary index. Once, because P appears once.
            if (aux_cnt(iatom) > 0) then
               gradient(comp, iatom) = gradient(comp, iatom) &
                                       - sum(ip2(:, :, q0:q1, comp)*gamma(:, :, q0:q1))

               ! The metric. Twice, for the same reason as the orbital pair:
               ! (P|Q) is symmetric and so is Omega.
               gradient(comp, iatom) = gradient(comp, iatom) &
                                       - 2.0_dp*sum(d2c(q0:q1, :, comp)*omega(q0:q1, :))
            end if
         end do
      end do

   end subroutine df_two_electron_gradient

   subroutine add_exchange_channel(three, jinv, orbitals, n_occ, weight, gamma, omega)
      !! One spin channel's exchange contribution to Gamma and Omega
      !!
      !! `weight` is two for a closed shell, where the single channel carries
      !! both spins, and one for each channel of an unrestricted reference. The
      !! metric term takes half of it, which is what makes the restricted and
      !! unrestricted expressions agree when the two spins are identical.
      ! TODO(mqc): the `f = J^-1 e` and metric loops below are naux^2 * n_occ^2
      ! written out by hand, where every other contraction in this module goes
      ! through `pic_gemm`. They dominate a fitted exchange gradient.
      real(dp), intent(in) :: three(:, :)     !! (nao*nao, naux)
      real(dp), intent(in) :: jinv(:, :)      !! (naux, naux)
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: weight
      real(dp), intent(inout) :: gamma(:, :, :)   !! (nao, nao, naux)
      real(dp), intent(inout) :: omega(:, :)      !! (naux, naux)

      real(dp), allocatable :: c_occ(:, :), e(:, :, :), f(:, :, :)
      real(dp), allocatable :: tmp(:, :), block(:, :), zp(:, :)
      integer :: nao, naux, ip, iq

      nao = size(orbitals, 1)
      naux = size(jinv, 1)

      if (n_occ <= 0) return

      allocate (c_occ(nao, n_occ), source=orbitals(:, 1:n_occ))
      allocate (e(n_occ, n_occ, naux), f(n_occ, n_occ, naux))
      allocate (tmp(nao, n_occ), block(nao, nao), zp(nao, nao))

      ! e^P_ij = C^T (uv|P) C, one auxiliary function at a time.
      do ip = 1, naux
         block = reshape(three(:, ip), [nao, nao])
         tmp = 0.0_dp
         call pic_gemm(block, c_occ, tmp)
         e(:, :, ip) = 0.0_dp
         call pic_gemm(c_occ, tmp, e(:, :, ip), transa="T")
      end do

      ! f = J^-1 e, over the auxiliary index for every occupied pair.
      f = 0.0_dp
      do ip = 1, naux
         do iq = 1, naux
            if (jinv(ip, iq) == 0.0_dp) cycle
            f(:, :, ip) = f(:, :, ip) + jinv(ip, iq)*e(:, :, iq)
         end do
      end do

      ! Z^P = C f^P C^T, back in the AO basis, subtracted from Gamma.
      do ip = 1, naux
         tmp = 0.0_dp
         call pic_gemm(c_occ, f(:, :, ip), tmp)
         zp = 0.0_dp
         call pic_gemm(tmp, c_occ, zp, transb="T")
         gamma(:, :, ip) = gamma(:, :, ip) - weight*zp
      end do

      ! The metric term, sum_ij f^P_ij f^Q_ij.
      do iq = 1, naux
         do ip = 1, naux
            omega(ip, iq) = omega(ip, iq) &
                            + 0.5_dp*weight*sum(f(:, :, ip)*f(:, :, iq))
         end do
      end do

   end subroutine add_exchange_channel

   pure function outer(a, b) result(m)
      !! The outer product a b^T, which Fortran has no intrinsic for
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: m(size(a), size(b))

      integer :: i

      do i = 1, size(b)
         m(:, i) = a*b(i)
      end do
   end function outer

   subroutine three_centre_deriv(orb, aux, which, deriv, omega)
      !! d(mu nu | P) / dR, as (nao, nao, naux, 3)
      !!
      !! `which` selects whose centre is differentiated: 1 for the first orbital
      !! index, through int3c2e_ip1, and 2 for the auxiliary one, through
      !! int3c2e_ip2. The second orbital index has no entry point of its own --
      !! (mu nu|P) is symmetric in mu and nu, so its derivative is the transpose
      !! of the first.
      !!
      !! **No symmetry in the shell-pair loop**, unlike the energy's
      !! `three_centre` next door: ip1 differentiates mu and not nu, so the block
      !! for (ish, jsh) is not the transpose of (jsh, ish).
      type(czt_molecule_t), intent(in) :: orb, aux
      integer, intent(in) :: which
      real(dp), allocatable, intent(out) :: deriv(:, :, :, :)
      real(dp), intent(in), optional :: omega
         !! Attenuate to `erf(omega r)/r`, matching `three_centre`. A fitted
         !! range-separated gradient differentiates two fits, and each has to be
         !! differentiated with the kernel it was made from.

      integer, allocatable :: bas(:, :)
      real(dp), allocatable :: env(:), buf(:)
      integer :: shls(4)
      integer :: dummy, nbas_orb, nbas_aux
      integer :: ish, jsh, ksh, di, dj, dk, io, jo, ko, i, j, k, comp, ret, idx, mx

      nbas_orb = orb%nbas
      nbas_aux = aux%nbas
      call build_df_shell_table(orb, aux, bas, env, dummy)
      if (present(omega)) env(LIBCINT_PTR_RANGE_OMEGA + 1) = omega

      allocate (deriv(orb%nao, orb%nao, aux%nao, 3))
      deriv = 0.0_dp

      mx = max(max_block(orb), max_block(aux))

      ! Threaded over the first orbital shell. Each (ish, jsh, ksh) writes its
      ! own block of `deriv` and reads nothing another writes, so no reduction
      ! is needed; `buf` is the one genuinely private object.
      !$omp parallel default(none) &
      !$omp    shared(deriv, orb, aux, bas, env, dummy, nbas_orb, nbas_aux, mx, which) &
      !$omp    private(ish, jsh, ksh, di, dj, dk, io, jo, ko, i, j, k, comp, ret, idx, &
      !$omp            shls, buf)
      allocate (buf(mx**3*3))
      !$omp do schedule(dynamic)
      do ish = 1, nbas_orb
         di = shell_dim(orb%cartesian, ish - 1, bas)
         io = orb%shell_offset(ish)
         do jsh = 1, nbas_orb
            dj = shell_dim(orb%cartesian, jsh - 1, bas)
            jo = orb%shell_offset(jsh)
            do ksh = 1, nbas_aux
               dk = shell_dim(orb%cartesian, nbas_orb + ksh - 1, bas)
               ko = aux%shell_offset(ksh)
               shls = [ish - 1, jsh - 1, nbas_orb + ksh - 1, dummy - 1]

               if (orb%cartesian) then
                  if (which == 1) then
                     ret = libcint_3c2e_ip1_cart(buf, shls, orb%atm, orb%natm, bas, &
                                                 dummy, env)
                  else
                     ret = libcint_3c2e_ip2_cart(buf, shls, orb%atm, orb%natm, bas, &
                                                 dummy, env)
                  end if
               else
                  if (which == 1) then
                     ret = libcint_3c2e_ip1_sph(buf, shls, orb%atm, orb%natm, bas, &
                                                dummy, env)
                  else
                     ret = libcint_3c2e_ip2_sph(buf, shls, orb%atm, orb%natm, bas, &
                                                dummy, env)
                  end if
               end if
               if (ret == 0) cycle

               do comp = 1, 3
                  do k = 1, dk
                     do j = 1, dj
                        do i = 1, di
                           idx = i + di*(j - 1 + dj*(k - 1 + dk*(comp - 1)))
                           deriv(io + i, jo + j, ko + k, comp) = buf(idx)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end do
      deallocate (buf)
      !$omp end parallel
   end subroutine three_centre_deriv

   subroutine two_centre_deriv(aux, deriv, omega)
      !! d(P|Q) / dR with respect to P's centre, as (naux, naux, 3)
      type(czt_molecule_t), intent(in) :: aux
      real(dp), allocatable, intent(out) :: deriv(:, :, :)
      real(dp), intent(in), optional :: omega   !! As `two_centre`

      real(dp), allocatable :: buf(:), env_local(:)
      integer :: shls(4)
      integer :: ish, jsh, di, dj, io, jo, i, j, comp, ret, idx, mx

      allocate (deriv(aux%nao, aux%nao, 3))
      deriv = 0.0_dp

      ! A copy: `aux%env` is shared and must not keep the omega afterwards.
      env_local = aux%env
      if (present(omega)) env_local(LIBCINT_PTR_RANGE_OMEGA + 1) = omega

      mx = max_block(aux)
      allocate (buf(mx*mx*3))

      do ish = 1, aux%nbas
         di = shell_dim(aux%cartesian, ish - 1, aux%bas)
         io = aux%shell_offset(ish)
         do jsh = 1, aux%nbas
            dj = shell_dim(aux%cartesian, jsh - 1, aux%bas)
            jo = aux%shell_offset(jsh)
            shls = [ish - 1, jsh - 1, 0, 0]

            if (aux%cartesian) then
               ret = libcint_2c2e_ip1_cart(buf, shls, aux%atm, aux%natm, aux%bas, &
                                           aux%nbas, env_local)
            else
               ret = libcint_2c2e_ip1_sph(buf, shls, aux%atm, aux%natm, aux%bas, &
                                          aux%nbas, env_local)
            end if
            if (ret == 0) cycle

            do comp = 1, 3
               do j = 1, dj
                  do i = 1, di
                     idx = i + di*(j - 1 + dj*(comp - 1))
                     deriv(io + i, jo + j, comp) = buf(idx)
                  end do
               end do
            end do
         end do
      end do
   end subroutine two_centre_deriv

end module mqc_czt_gradient
