!! Analytic nuclear gradients on the CPU backend
module mqc_libcint_gradient
   !! The derivative of a converged SCF energy with respect to the nuclear
   !! coordinates, for a restricted or unrestricted reference.
   !!
   !! **Why a variational wavefunction makes this cheap.** The energy is
   !! stationary with respect to the orbitals, so their derivatives do not
   !! appear: nothing here solves a response equation. What is left is the
   !! derivative of the *integrals*, contracted with densities the SCF already
   !! produced. That is why Hartree-Fock gradients come before MP2 or coupled
   !! cluster ones, which are not stationary and do need a response solve.
   !!
   !! **The four terms.**
   !!
   !!   * nuclear repulsion, which is analytic and needs no integrals
   !!   * the core Hamiltonian, `Tr(D dH/dR)`
   !!   * the two-electron term, `Tr(D dG/dR)`
   !!   * the Pulay term, `-Tr(W dS/dR)`, from the basis functions moving with
   !!     their atoms. This is the one that has no analogue in a plane-wave
   !!     code, and `W` is the energy-weighted density rather than the density
   !!
   !! **libcint differentiates the bra, and only the bra.** Every `ip` entry
   !! point returns the derivative with respect to the first shell's centre;
   !! the ket's derivative is recovered from translational invariance rather
   !! than integrated for. That is what the transposes and the factors of two
   !! below are doing, and it is why the sum of the returned gradient over
   !! atoms is zero by construction only if they are right -- which is worth
   !! checking, and is checked in the tests.
   use pic_types, only: dp
   use mqc_nuclear_repulsion, only: add_nuclear_repulsion_gradient
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim, max_block, atom_ao_blocks, &
                                    build_df_shell_table, three_centre, two_centre, &
                                    metric_inverse_sqrt, eri_shell_table_t, &
                                    eri_shell_table, eri_schwarz_collapse
   use mqc_libcint_ao, only: eval_ao_block, AO_POINT_BLOCK, AO_HESS_COMP
   use mqc_libcint_xc, only: xc_context_t, xc_grid_lda_quantities, &
                             xc_grid_gga_quantities, xc_grid_kernel_quantities
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
                              LIBCINT_NPRIM_OF
   use mqc_libcint_direct, only: schwarz_bounds
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
   implicit none
   private

   public :: libcint_scf_gradient
   public :: nuclear_repulsion_gradient   !! Exposed for the tests
   public :: three_centre_deriv           !! Exposed for the tests
   public :: two_centre_deriv             !! Exposed for the tests
   ! The pieces an MP2 gradient assembles for itself. It cannot call
   ! `libcint_scf_gradient` and add to the answer -- the relaxed density enters
   ! the one-electron and overlap terms in place of the reference one, not
   ! alongside it -- so it needs the same derivative integrals rather than the
   ! same assembly.
   public :: one_electron_deriv
   public :: two_electron_deriv
      !! Exposed for the MCSCF gradient, whose inactive and active blocks are
      !! separable and contract exactly as a closed-shell reference does.
   public :: iprinv_deriv_at
   public :: xc_potential_gradient
   public :: fitted_reference_gradient
   public :: DERIV_OVLP, DERIV_KIN, DERIV_NUC

   ! Which one-electron derivative `one_electron_deriv` should build.
   integer, parameter :: DERIV_OVLP = 1
   integer, parameter :: DERIV_KIN = 2
   integer, parameter :: DERIV_NUC = 3

   real(dp), parameter :: DERIV_SCREEN_TOL = 1.0e-10_dp
      !! Contribution below which a differentiated shell quartet is dropped.
      !!
      !! Measured on a seventeen-atom closed shell in 6-31G, against the same
      !! gradient computed with no screening at all:
      !!
      !!     tolerance     time      error in |g|
      !!     none         16.7 s      --
      !!     1e-12        14.9 s      1e-11
      !!     1e-10        13.7 s      7e-11
      !!     1e-9         12.9 s      5e-9
      !!
      !! 1e-10 is the knee. Below it the error starts moving a decade per decade
      !! while the time barely does, and an optimiser converging on a 1e-4
      !! gradient norm has no use for the accuracy being bought. The margin
      !! matters here because the bra bound this multiplies is an estimate
      !! rather than a rigorous Cauchy-Schwarz one -- see `two_electron_deriv`.
      !!
      !! Note how little any of this saves. On a molecule this compact almost
      !! every quartet genuinely contributes, and the same is true of the
      !! undifferentiated build next door -- `build_fock_direct` screens 15 % at
      !! 1e-11. Screening is what makes a *large* system tractable; it is not
      !! where a small one's time goes, and it should not be tuned as though it
      !! were.

contains

   subroutine libcint_scf_gradient(mol, density, density_beta, orbitals, orbitals_beta, &
                                   orbital_energies, orbital_energies_beta, &
                                   n_occupied, n_occupied_beta, gradient, error, aux, xc)
      !! dE/dR for a converged SCF, in Hartree/Bohr
      !!
      !! The beta arguments are what makes this one routine rather than two.
      !! Present means unrestricted: the Coulomb field is built from the total
      !! density and exchange from each spin separately, which is the only
      !! place the two paths differ. Absent means closed shell, where the
      !! density already carries its factor of two.
      !!
      !! `aux` present means the SCF being differentiated fitted its J and K,
      !! and only the two-electron term changes: a density-fitted SCF is still
      !! variational, so the one-electron, Pulay and nuclear-repulsion terms
      !! are the same integrals contracted with the same densities. Passing an
      !! auxiliary basis that the SCF did not use would differentiate a
      !! different energy, which nothing here can detect -- it is the caller's
      !! job to pass the one it converged with.
      type(libcint_molecule_t), intent(in) :: mol
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
      type(libcint_molecule_t), intent(in), optional :: aux
         !! The auxiliary basis the SCF fitted with. Absent means exact ERIs.
      type(xc_context_t), intent(inout), optional :: xc
         !! The functional the SCF used. Absent means Hartree-Fock. Present adds
         !! the exchange-correlation term and scales exact exchange by the
         !! fraction the functional asks for -- zero for a pure functional, so
         !! the Fock exchange derivative is not merely small but absent.

      real(dp), allocatable :: total_density(:, :), weighted(:, :)
      real(dp), allocatable :: s1(:, :, :), h1(:, :, :), kin(:, :, :)
      real(dp), allocatable :: vhf(:, :, :), vhf_beta(:, :, :)
      real(dp), allocatable :: vrinv(:, :, :), hcore_a(:, :, :)
      real(dp) :: exx
      integer, allocatable :: offsets(:), counts(:)
      integer :: nao, iatom, comp, p0, p1
      logical :: unrestricted
      logical :: fitted

      unrestricted = present(density_beta)
      nao = mol%nao

      ! Hartree-Fock keeps all of its exchange; a functional keeps the fraction
      ! it declares. A range-separated one needs a second exchange derivative at
      ! a screened omega, which this does not build, so it is refused below
      ! rather than silently given its full-range part.
      exx = 1.0_dp
      if (present(xc)) exx = xc%exx_fraction

      if (size(density, 1) /= nao) then
         call error%set(ERROR_VALIDATION, &
                        "libcint_scf_gradient: the density does not match this basis")
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
      else if (unrestricted) then
         call two_electron_deriv(mol, total_density, vhf, error, exx_fraction=exx, &
                                 exchange_density=density)
         if (error%has_error()) return
         call two_electron_deriv(mol, total_density, vhf_beta, error, exx_fraction=exx, &
                                 exchange_density=density_beta)
         if (error%has_error()) return
      else
         call two_electron_deriv(mol, total_density, vhf, error, exx_fraction=exx)
         if (error%has_error()) return
      end if

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)

      allocate (vrinv(nao, nao, 3), hcore_a(nao, nao, 3))

      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)

         ! dH/dR_A has two parts that look alike and are not. Moving atom A
         ! moves the basis functions centred on it, which is the h1 block
         ! below; it also moves the *nucleus*, and every electron feels that
         ! wherever its orbital sits. The second is the Hellmann-Feynman term,
         ! it involves no basis-set derivative at all, and it is what `iprinv`
         ! is for -- differentiating the origin of 1/|r-R| rather than a
         ! basis-function centre.
         call iprinv_deriv_at(mol, iatom, vrinv)
         vrinv = -mol%charges(iatom)*vrinv
         if (counts(iatom) > 0) then
            vrinv(p0:p1, :, :) = vrinv(p0:p1, :, :) + h1(p0:p1, :, :)
         end if
         do comp = 1, 3
            hcore_a(:, :, comp) = vrinv(:, :, comp) + transpose(vrinv(:, :, comp))
         end do

         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + sum(hcore_a(:, :, comp)*total_density)
         end do

         if (counts(iatom) == 0) cycle

         ! Twice, because the nabla was on the bra: the ket's contribution is
         ! the same number by the symmetry of the integrals, so it is counted
         ! rather than computed.
         do comp = 1, 3
            if (.not. fitted) then
               if (unrestricted) then
                  gradient(comp, iatom) = gradient(comp, iatom) &
                                          + 2.0_dp*sum(vhf(p0:p1, :, comp)*density(p0:p1, :)) &
                                          + 2.0_dp*sum(vhf_beta(p0:p1, :, comp)*density_beta(p0:p1, :))
               else
                  gradient(comp, iatom) = gradient(comp, iatom) &
                                          + 2.0_dp*sum(vhf(p0:p1, :, comp)*total_density(p0:p1, :))
               end if
            end if
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    - 2.0_dp*sum(s1(p0:p1, :, comp)*weighted(p0:p1, :))
         end do
      end do

      ! Last, and on top of everything above: the exchange-correlation term is
      ! an addition to the Hartree-Fock gradient, not a replacement for any part
      ! of it.
      if (present(xc)) then
         ! The second exchange derivative at the screened omega now exists,
         ! but only on the fitted path -- `df_two_electron_gradient` takes an
         ! `rs_omega`, while `two_electron_deriv` on the exact path still builds
         ! one kernel per call and has no second pass wired to it.
         ! `.or. unrestricted`, and that half is not hypothetical-only: the
         ! fitted long-range pass below is written in the *restricted* branch
         ! alone. Without this, unrestricted + fitted + range separated falls
         ! through both the pass and the refusal and returns a gradient missing
         ! wB97X's dominant exchange term. It is unreachable today because the
         ! SCF refuses unrestricted fitting up front -- which is exactly the
         ! shape the fitted-exchange-fraction bug had before Kohn-Sham fitting
         ! became reachable and turned it live.
         if (xc%range_separated .and. (.not. fitted .or. unrestricted)) then
            call error%set(ERROR_VALIDATION, "the gradient of a range-separated "// &
                           "functional needs a second exchange derivative at the "// &
                           "screened omega. The density-fitted path builds it; the "// &
                           "exact-ERI one does not yet. Give an auxiliary basis, or "// &
                           "take the gradient at a functional that is not range "// &
                           "separated.")
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

   end subroutine libcint_scf_gradient

   subroutine xc_gradient(ctx, mol, density, gradient, error, density_beta)
      !! The exchange-correlation contribution to dE/dR, accumulated in place
      !!
      !! Three terms, and the last two are the ones a first attempt leaves out.
      !!
      !! **The basis functions move.** rho is built from functions centred on
      !! atoms, so moving atom A changes rho wherever A's functions reach. This
      !! is the term that looks like the Pulay term of the one-electron part and
      !! is written the same way.
      !!
      !! **The grid points move.** The quadrature is atom-centred, so the points
      !! belonging to A travel with it and the integrand is sampled at moved
      !! positions. This contributes the same quantity as the first term but
      !! summed over A's *points* rather than A's *functions* -- over every
      !! basis function, not only the ones on A -- and with the opposite sign.
      !!
      !! **The partition weights move.** The Becke weight of every point depends
      !! on every nuclear position, so moving A reweights the whole grid,
      !! including points it does not own. This is the term whose absence leaves
      !! a gradient that looks entirely plausible and misses finite differences
      !! at about 1e-4, and it is also the term translational invariance is
      !! sensitive to -- which is why that check is worth running on a DFT
      !! gradient even though it is free.
      !!
      !! LDA only. `xc_grid_lda_quantities` refuses anything else, because a GGA
      !! needs second derivatives of the basis functions and `eval_ao_block`
      !! produces first.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
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
      integer, allocatable :: offsets(:), counts(:)
      real(dp) :: contrib(3)
      real(dp) :: wv, scale
      integer :: npts, nao, natm, g0, g1, nb, ig, gg, mu, comp, ia, own
      logical :: unrestricted

      if (.not. ctx%active) return

      unrestricted = present(density_beta)
      nao = mol%nao
      natm = mol%natm
      npts = ctx%grid%n_points

      gga = ctx%any_gga .or. ctx%any_mgga
      mgga = ctx%any_mgga

      if (mgga) then
         ! `vtau` is what makes this legal: `xc_grid_gga_quantities` refuses a
         ! meta-GGA outright when the caller cannot take it, because its own
         ! per-functional dispatch would otherwise fall through to the LDA
         ! branch. Unrestricted is refused in there for a separate reason.
         ! `density_beta` where there is one, so the meta-GGA refusal inside
         ! fires for the reason it was written. Without it that routine sees a
         ! restricted call, its own guard never runs, and what refuses an
         ! unrestricted meta-GGA gradient is the spin-consistency check further
         ! down -- correct by accident, with a message about internal state
         ! rather than about the missing feature.
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

      allocate (offsets(natm), counts(natm))
      call atom_ao_blocks(mol, offsets, counts)

      ! A closed-shell density already carries its factor of two, so the
      ! derivative of rho with respect to a moving basis function is 2*D*chi
      ! either way -- the two is the bra/ket pair, not the occupation.
      scale = 2.0_dp

      do g0 = 1, npts, AO_POINT_BLOCK
         g1 = min(g0 + AO_POINT_BLOCK - 1, npts)
         nb = g1 - g0 + 1

         if (gga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, &
                               grad=ao_grad, hess=ao_hess)
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, grad=ao_grad)
         end if
         if (error%has_error()) return

         ! (D chi)_mu(g) = sum_nu D_mu,nu chi_nu(g), the partner every term
         ! below contracts the basis-function gradient against.
         if (allocated(dchi)) deallocate (dchi)
         allocate (dchi(nb, nao))
         call density_times_ao(ao, density, nb, nao, dchi)
         if (unrestricted) then
            if (allocated(dchi_beta)) deallocate (dchi_beta)
            allocate (dchi_beta(nb, nao))
            call density_times_ao(ao, density_beta, nb, nao, dchi_beta)
         end if

         ! (D grad chi)_mu(g), which only the GGA term needs: it is the partner
         ! for the piece where the two first derivatives pair with each other.
         if (gga) then
            if (allocated(dgchi)) deallocate (dgchi)
            allocate (dgchi(nb, nao, 3))
            call density_times_ao_grad(ao_grad, density, nb, nao, dgchi)
            if (unrestricted) then
               if (allocated(dgchi_beta)) deallocate (dgchi_beta)
               allocate (dgchi_beta(nb, nao, 3))
               call density_times_ao_grad(ao_grad, density_beta, nb, nao, dgchi_beta)
            end if
         end if

         do ig = 1, nb
            gg = g0 + ig - 1
            own = ctx%grid%atom(gg)

            wv = ctx%grid%weights(gg)*vrho(gg)
            call accumulate_channel(ao_grad, dchi, ig, nao, wv, scale, &
                                    offsets, counts, natm, own, gradient)
            if (gga) then
               call accumulate_gga_channel(ao_grad, ao_hess, dchi, dgchi, ig, nao, &
                                           ctx%grid%weights(gg)*gcoef(gg, :), scale, &
                                           offsets, counts, natm, own, gradient)
            end if
            if (mgga) then
               call accumulate_mgga_channel(ao_grad, ao_hess, dgchi, ig, nao, &
                                            ctx%grid%weights(gg)*vtau(gg), scale, &
                                            offsets, counts, natm, own, gradient)
            end if
            if (unrestricted) then
               wv = ctx%grid%weights(gg)*vrho_beta(gg)
               call accumulate_channel(ao_grad, dchi_beta, ig, nao, wv, scale, &
                                       offsets, counts, natm, own, gradient)
               if (gga) then
                  call accumulate_gga_channel(ao_grad, ao_hess, dchi_beta, dgchi_beta, &
                                              ig, nao, &
                                              ctx%grid%weights(gg)*gcoef_beta(gg, :), &
                                              scale, offsets, counts, natm, own, gradient)
               end if
            end if
         end do

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
            ! piece is what makes this term cancel the integrand's own
            ! grid-motion contribution above -- together they are the integral
            ! of grad(P f) over one atomic subgrid, which the quadrature
            ! evaluates as very nearly zero because P f decays. Leaving it out
            ! does not make the answer slightly worse: it leaves the
            ! grid-motion term uncancelled, which on water is an error of 0.77
            ! Hartree/Bohr against a gradient of 0.06.
            !
            ! It needs no new derivative. Displacing the point is the same as
            ! displacing every nucleus the other way, so the gradient of the
            ! partition with respect to the point is minus the sum over atoms
            ! of what `becke_partition_derivatives` already returned.
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

   subroutine xc_potential_gradient(ctx, mol, density, pmat, gradient, error)
      !! `d/dR Tr(P V_xc[D])`, with both densities held fixed, accumulated in place
      !!
      !! **Why a double hybrid needs this and no gradient before it did.** Every
      !! gradient up to here differentiates an energy that is stationary in the
      !! orbitals, so the reference operator's derivative never appears on its
      !! own. A double hybrid's perturbative term is not stationary: eliminating
      !! the orbital response leaves the Z-vector contracted with the derivative
      !! of the *reference* operator, and for a Kohn-Sham reference that operator
      !! contains `V_xc`. Omitting this term gives a gradient of entirely
      !! plausible magnitude -- it is not a correction, it is a missing piece of
      !! the same order as the Coulomb one beside it.
      !!
      !! **`P` is not a density in the usual sense.** It is the relaxed
      !! correlation density plus the Z-vector, so it is symmetric but
      !! indefinite, integrates to zero rather than to the electron count, and is
      !! nowhere near idempotent. Nothing here assumes otherwise; in particular
      !! the functional is evaluated at `D` alone and `P` only ever appears
      !! linearly, which is what makes this the derivative of a *linear form*
      !! rather than of an energy.
      !!
      !!     G(R) = int w(R) [ v_rho rho_P + 2 v_sigma grad rho . grad rho_P ]
      !!
      !! and `dG/dR` has the same three sources the exchange-correlation gradient
      !! has -- the basis functions move, the grid points move, the partition
      !! weights move -- with one addition: `v_rho` and `v_sigma` are themselves
      !! functions of the reference density, which moves too. That is where the
      !! second derivatives enter, and it is the part with no analogue next door.
      !!
      !! **Four channels, two of them the reference's.** Collecting terms by
      !! which object is being differentiated:
      !!
      !!     d rho      : f_rr rho_P + 2 f_rs (grad rho . grad rho_P)
      !!     d grad rho : 2 f_rs rho_P grad rho
      !!                  + 4 f_ss (grad rho . grad rho_P) grad rho
      !!                  + 2 v_sigma grad rho_P
      !!     d rho_P    : v_rho
      !!     d grad rho_P: 2 v_sigma grad rho
      !!
      !! Each is the same contraction `accumulate_channel` and
      !! `accumulate_gga_channel` already perform for the energy gradient -- what
      !! differs is the density they carry and the coefficient they are weighted
      !! by, neither of which those routines know about. So this adds no new
      !! arithmetic of the kind that is hard to get right, only new coefficients.
      !!
      !! LDA and GGA; a meta-GGA is refused upstream by
      !! `xc_grid_kernel_quantities`.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)   !! The converged reference density
      real(dp), intent(in) :: pmat(:, :)      !! The density `V_xc` is traced against
      real(dp), intent(inout) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: rho(:), rho_grad(:, :), vrho(:), vsigma(:)
      real(dp), allocatable :: frr(:), frs(:), fss(:)
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), allocatable :: dchi(:, :), dgchi(:, :, :)
      real(dp), allocatable :: pchi(:, :), pgchi(:, :, :)
      real(dp), allocatable :: rho_p(:), p_grad(:, :)
      real(dp), allocatable :: dpart(:, :, :)
      real(dp) :: wg_ref(3), wg_p(3)
      real(dp) :: gdotp, coef_rho, integrand, w
      integer, allocatable :: offsets(:), counts(:)
      integer :: npts, nao, natm, g0, g1, nb, ig, gg, id, ia, comp, own
      logical :: gga

      if (.not. ctx%active) return

      nao = mol%nao
      natm = mol%natm
      npts = ctx%grid%n_points
      gga = ctx%any_gga

      call xc_grid_kernel_quantities(ctx, mol, density, rho, rho_grad, vrho, &
                                     vsigma, frr, frs, fss, error)
      if (error%has_error()) return

      allocate (offsets(natm), counts(natm))
      call atom_ao_blocks(mol, offsets, counts)

      do g0 = 1, npts, AO_POINT_BLOCK
         g1 = min(g0 + AO_POINT_BLOCK - 1, npts)
         nb = g1 - g0 + 1

         ! Second derivatives of the basis functions, for the same reason the
         ! GGA energy gradient needs them: differentiating `grad chi` with
         ! respect to a nuclear coordinate is a second derivative.
         if (gga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, &
                               grad=ao_grad, hess=ao_hess)
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, grad=ao_grad)
         end if
         if (error%has_error()) return

         if (allocated(dchi)) deallocate (dchi, pchi)
         allocate (dchi(nb, nao), pchi(nb, nao))
         call density_times_ao(ao, density, nb, nao, dchi)
         call density_times_ao(ao, pmat, nb, nao, pchi)

         ! rho_P and its gradient, on the same footing as the reference's. Built
         ! from `pchi` directly rather than through `eval_rho`, since the partner
         ! matrices are wanted anyway and this saves a second pass.
         if (allocated(rho_p)) deallocate (rho_p, p_grad)
         allocate (rho_p(nb), p_grad(nb, 3))
         rho_p = 0.0_dp
         p_grad = 0.0_dp
         do ig = 1, nb
            rho_p(ig) = sum(pchi(ig, :)*ao(ig, :))
         end do
         if (gga) then
            if (allocated(dgchi)) deallocate (dgchi, pgchi)
            allocate (dgchi(nb, nao, 3), pgchi(nb, nao, 3))
            call density_times_ao_grad(ao_grad, density, nb, nao, dgchi)
            call density_times_ao_grad(ao_grad, pmat, nb, nao, pgchi)
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
            call accumulate_channel(ao_grad, dchi, ig, nao, w*coef_rho, 2.0_dp, &
                                    offsets, counts, natm, own, gradient)

            ! `P`'s own channel, weighted by the ordinary potential. This is the
            ! term that survives when the functional is a constant -- and the
            ! only one an LDA-shaped first attempt would write.
            call accumulate_channel(ao_grad, pchi, ig, nao, w*vrho(gg), 2.0_dp, &
                                    offsets, counts, natm, own, gradient)

            if (gga) then
               do id = 1, 3
                  wg_ref(id) = w*(2.0_dp*frs(gg)*rho_p(ig)*rho_grad(gg, id) &
                                  + 4.0_dp*fss(gg)*gdotp*rho_grad(gg, id) &
                                  + 2.0_dp*vsigma(gg)*p_grad(ig, id))
                  wg_p(id) = w*2.0_dp*vsigma(gg)*rho_grad(gg, id)
               end do
               call accumulate_gga_channel(ao_grad, ao_hess, dchi, dgchi, ig, nao, &
                                           wg_ref, 2.0_dp, offsets, counts, natm, &
                                           own, gradient)
               call accumulate_gga_channel(ao_grad, ao_hess, pchi, pgchi, ig, nao, &
                                           wg_p, 2.0_dp, offsets, counts, natm, &
                                           own, gradient)
            end if
         end do

         ! The partition weights, on the integrand of this linear form rather
         ! than on the energy density. Same two pieces as the energy gradient:
         ! every nucleus reweights the whole grid, and a point owned by A also
         ! moves with A.
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
                                        gradient, k_scale)
      !! `d/dR Tr(G[A] B)` when `G` is built from fitted integrals
      !!
      !! **What this replaces.** With exact integrals the reference part of a
      !! correlated gradient is a four-centre `int2e_ip1` contraction --
      !! `two_electron_mp2_terms` builds it, and the assembly contracts the
      !! result against the relaxed density. Fitted, the same quantity has no
      !! four-centre derivative anywhere in it: `(mu nu|lam sig)` is a product
      !! of three-centre integrals and a metric, and differentiating a product
      !! gives terms of each kind. So this returns the two intermediates that
      !! contract against those, in the shapes the correlation's own fitted
      !! terms already use.
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
      !! `n_aux^2` intermediate `W_PQ = Tr(B M^P A M^Q)` it comes from: applying
      !! `J^-1` on both sides of `W` is what the derivative wants anyway, and
      !! `S` already carries one of the two applications.
      !!
      !! **The caller halves this.** What an assembly wants is
      !! `(1/2) Tr(G^x[D] D) + Tr(G^x[D] P)` -- the reference's own two-electron
      !! gradient plus the cross term -- which is half of `Tr(G^x[D](D + 2P))`.
      !! The exact path reaches the same half by building only two of the four
      !! differentiated positions and letting the integral's symmetry supply the
      !! rest, which is invisible from its call site. Left in, the water
      !! gradient comes out thirty times too large with translational invariance
      !! holding at 1e-14 throughout, so this returns the honest derivative and
      !! the factor is written where it can be read.
      type(libcint_molecule_t), intent(in) :: mol, aux
      real(dp), intent(in) :: three(:, :)   !! Raw `(mu nu|P)`, (n_ao^2, naux)
      real(dp), intent(in) :: jm12(:, :)    !! `J^(-1/2)`, as the energy fitted with
      real(dp), intent(in) :: dm_a(:, :)    !! The density inside the operator
      real(dp), intent(in) :: dm_b(:, :)    !! The density it is contracted against
      real(dp), intent(inout) :: gradient(:, :)   !! (3, natm), accumulated into
      real(dp), intent(in), optional :: k_scale

      real(dp), allocatable :: jinv(:, :), y(:, :), s(:, :), gamma(:, :, :)
      real(dp), allocatable :: omega(:, :), gt(:, :, :)
      real(dp), allocatable :: ip1(:, :, :, :), ip2(:, :, :, :), j1(:, :, :)
      real(dp), allocatable :: rho_a(:), rho_b(:), work(:, :)
      integer, allocatable :: offsets(:), counts(:)
      integer, allocatable :: aux_offsets(:), aux_counts(:)
      integer :: n, naux, p, iatom, comp, p0, p1, q0, q1
      real(dp) :: kf

      kf = 1.0_dp
      if (present(k_scale)) kf = k_scale

      n = mol%nao
      naux = aux%nao

      allocate (jinv(naux, naux))
      call pic_gemm(jm12, jm12, jinv)

      ! `Y^P = sum_Q J^-1_PQ M^Q`, wanted by both halves of the exchange term
      ! and by the two `rho` vectors, which are `<Y^P, D>` rather than a second
      ! metric solve.
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

      ! As a gemm. `dchi = ao D^T`, and `D` is symmetric so the transpose is
      ! free -- it is written anyway, because the expression above is what the
      ! loop said and a routine that quietly depends on symmetry is one that
      ! breaks when somebody passes an energy-weighted matrix that is not.
      call pic_gemm(ao, density, dchi, transb="T")
   end subroutine density_times_ao

   subroutine accumulate_channel(ao_grad, dchi, ig, nao, wv, scale, &
                                 offsets, counts, natm, own, gradient)
      !! One spin channel's basis-function and grid-point terms at one point
      !!
      !! Both are the same contraction of `grad chi` against `D chi`. They
      !! differ in which basis functions are summed and where the result lands:
      !! the basis-function term takes only the functions on atom A and lands on
      !! A, the grid-point term takes every function and lands on the atom that
      !! owns the point. Their signs are opposite, which is what makes the two
      !! cancel when every atom moves together.
      real(dp), intent(in) :: ao_grad(:, :, :)
      real(dp), intent(in) :: dchi(:, :)
      integer, intent(in) :: ig, nao, natm, own
      real(dp), intent(in) :: wv, scale
      integer, intent(in) :: offsets(:), counts(:)
      real(dp), intent(inout) :: gradient(:, :)

      integer :: ia, mu, comp, p0, p1
      real(dp) :: block_sum(3), total(3)

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
      do comp = 1, 3
         gradient(comp, own) = gradient(comp, own) + scale*wv*total(comp)
      end do
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
      ! `(nb, nao) x (nao, nao)` and no copy is made.
      !
      ! This was the single most expensive routine in a density functional
      ! gradient run -- 382 s of CPU on twenty waters, against 684 s for the
      ! whole gradient phase -- and the reason is the loop it replaces rather
      ! than the arithmetic, which is unchanged. Written out, the innermost
      ! index walked `ao_grad(ig, nu, k)` along `nu`, which strides by `nb`:
      ! one useful element per cache line, and nothing a vector unit can hold.
      ! A gemm blocks the same multiply to fit the cache instead.
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
      !! **The half is already gone by the time it is written that way**, and
      !! that is the factor to be careful about. `tau` carries libxc's 1/2, and
      !! differentiating pairs the moving function with either the bra or the
      !! ket, which are equal by the symmetry of D -- so the two cancel and the
      !! net coefficient is one, not one half and not two. It is applied through
      !! `scale` and an explicit half below so that both appear where the
      !! corresponding factors in `accumulate_xc_matrix` do, rather than being
      !! silently folded into a single constant nobody can check.
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
      if (nao < 0) return
   end subroutine accumulate_mgga_channel

   subroutine accumulate_gga_channel(ao_grad, ao_hess, dchi, dgchi, ig, nao, wg, scale, &
                                     offsets, counts, natm, own, gradient)
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
      !! with *each other*, one on the bra and one on the ket. Writing this as a
      !! single doubled term gives a gradient of entirely plausible magnitude
      !! that misses finite differences -- which is why the check that decides
      !! this routine is the numerical one and not translational invariance.
      !!
      !! `wg` is the weight times dE/d(grad rho), the conjugate quantity
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
      real(dp), intent(in) :: wg(3)
      real(dp), intent(in) :: scale
      integer, intent(in) :: offsets(:), counts(:)
      real(dp), intent(inout) :: gradient(:, :)

      ! Which packed component holds d2/dx_k dx_j.
      integer, parameter :: HESS_INDEX(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

      integer :: ia, mu, j, k, p0, p1
      real(dp) :: block_sum(3)
      real(dp) :: total(3)
      real(dp) :: acc

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
                  ! the second. Pairing j with the bra instead -- the natural
                  ! misreading -- is wrong by about six percent on water/PBE,
                  ! and translational invariance does not notice.
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
         gradient(k, own) = gradient(k, own) + scale*total(k)
      end do
   end subroutine accumulate_gga_channel

   subroutine nuclear_repulsion_gradient(mol, gradient)
      !! dE_NN/dR, accumulated into `gradient`
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(inout) :: gradient(:, :)

      call add_nuclear_repulsion_gradient(mol%charges, mol%coords, gradient)
   end subroutine nuclear_repulsion_gradient

   subroutine energy_weighted_density(orbitals, energies, n_occ, closed_shell, weighted)
      !! W = sum_i occ_i eps_i C_i C_i^T
      !!
      !! The Pulay term contracts against this rather than the density,
      !! because what moves with the basis functions is the orthonormality
      !! constraint, and its Lagrange multipliers are the orbital energies.
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
      !! component slowest, so a block is buf(di, dj, 1:3).
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), allocatable, intent(out) :: matrix(:, :, :)
      integer, intent(in) :: which

      real(dp), allocatable :: buf(:)
      integer :: shls(2)
      integer :: ish, jsh, di, dj, i, j, comp, ret, io, jo, mx

      allocate (matrix(mol%nao, mol%nao, 3))
      matrix = 0.0_dp

      ! Flat, and indexed by hand below. libcint packs a block with the
      ! *actual* shell dimensions, so a rank-3 buffer declared at the largest
      ! shell has the wrong strides for every smaller one -- which is invisible
      ! in a basis of s functions and wrong everywhere else.
      mx = max_block(mol)
      allocate (buf(mx*mx*3))

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
   end subroutine one_electron_deriv

   subroutine iprinv_deriv_at(mol, iatom, matrix)
      !! The 1/|r-R| derivative with the operator origin on atom `iatom`
      !!
      !! `env` is copied rather than modified, because the origin lives in it
      !! and the molecule is shared. One copy per atom of an array that is a
      !! few kilobytes is not worth avoiding; mutating the caller's molecule
      !! and putting it back would be.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: iatom
      real(dp), intent(out) :: matrix(:, :, :)

      real(dp), allocatable :: env(:), buf(:)
      integer :: shls(2)
      integer :: ish, jsh, di, dj, i, j, comp, ret, io, jo, mx

      allocate (env(size(mol%env)), source=mol%env)
      env(LIBCINT_PTR_RINV_ORIG + 1:LIBCINT_PTR_RINV_ORIG + 3) = mol%coords(:, iatom)

      matrix = 0.0_dp
      mx = max_block(mol)
      allocate (buf(mx*mx*3))

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
   end subroutine iprinv_deriv_at

   subroutine two_electron_deriv(mol, density, vhf, error, exchange_density, exx_fraction, &
                                 screen_tol)
      !! `J - K/2` built from the differentiated ERIs, as (nao, nao, 3)
      !!
      !! `density` drives the Coulomb field and is the total density in both
      !! the restricted and unrestricted cases. `exchange_density` is the one
      !! exchange is built from: absent for a closed shell, where the total
      !! density carries its factor of two and `J - K/2` is the right
      !! combination; present for one spin of an unrestricted reference, where
      !! exchange takes that spin alone and the coefficient is one rather than
      !! a half.
      !!
      !! **One permutation survives, and it is used.** A differentiated quartet
      !! has none of the eightfold symmetry of an ordinary one: the nabla
      !! distinguishes the first index from the second and the bra from the
      !! ket. What is left untouched is the ket pair, `(∇i j|k l) = (∇i j|l k)`,
      !! so this walks `l <= k` and counts the off-diagonal pair twice -- once
      !! as itself and once transposed, which is not the same contraction and
      !! has to be written out rather than doubled.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: vhf(:, :, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: exchange_density(:, :)
      real(dp), intent(in), optional :: exx_fraction   !! Default one, for Hartree-Fock
      real(dp), intent(in), optional :: screen_tol
         !! Bound below which a shell quartet is skipped, on the *contribution*
         !! rather than the integral. Default `DERIV_SCREEN_TOL`.

      real(dp), allocatable :: buf(:), vj(:, :, :), vk(:, :, :)
      real(dp), allocatable :: vj_local(:, :, :), vk_local(:, :, :)
      real(dp), allocatable :: exchange_from(:, :)
      real(dp), allocatable :: bounds(:, :), bq(:, :), bra_bound(:, :), dsh(:, :), esh(:, :)
      integer, allocatable :: dims(:), offs(:)
      real(dp) :: g, k_scale, tol, bra, est, amax
      type(c_ptr) :: opt
      type(eri_shell_table_t) :: tab
      integer :: shls(4)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: io, jo, ko, lo, i, j, k, l, comp, ret, mx, idx, nao, nbas
      integer :: nprim, ptr

      ! The shells the quartet loop runs over: the fused-sp view when the
      ! molecule carries one, its split shells otherwise -- the same table the
      ! direct Fock build takes, for the same 6-31G-shaped reason. This was
      ! split-only for a while: libfint's int2e_ip1 applied an L shell's s and
      ! p coefficients on the wrong stride whenever an integral carries more
      ! than one tensor component, so the view was correct for int2e and
      ! garbage here. Fixed in libfint (l_shell_grad_check guards it), and the
      ! AO order of a fused shell is its split shells' order, so the scatter
      ! into vj/vk needs nothing beyond the view's own offsets.
      call eri_shell_table(mol, tab)
      mx = tab%block_max
      nao = mol%nao
      nbas = tab%nbas
      allocate (vj(nao, nao, 3), vk(nao, nao, 3))
      vj = 0.0_dp
      vk = 0.0_dp

      tol = DERIV_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol

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
         call libcint_2e_ip1_cart_optimizer(opt, mol%atm, mol%natm, tab%bas, tab%nbas, tab%env)
      else
         call libcint_2e_ip1_sph_optimizer(opt, mol%atm, mol%natm, tab%bas, tab%nbas, tab%env)
      end if

      ! Shell dimensions and offsets once, rather than a `shell_dim` call --
      ! which reaches into `bas` and recomputes an angular-momentum formula --
      ! at every level of a four-deep loop that runs millions of times.
      dims = tab%dims
      offs = tab%offs(1:nbas)

      ! Screening. The undifferentiated loop in `build_fock_direct` has had this
      ! from the start; this one had none at all and evaluated every quartet the
      ! basis admits, which is why a gradient cost eight Fock builds.
      !
      ! Two factors go into the estimate:
      !
      !   * **The bra carries a derivative, so its Schwarz bound is not the
      !     ordinary one.** Differentiating a primitive brings down `2*alpha*r`,
      !     and the extent of that primitive is `r ~ 1/sqrt(alpha)`, so the
      !     distribution grows by roughly `2*sqrt(alpha)`. Taking the largest
      !     exponent in the shell makes that an over- rather than an
      !     under-estimate, which is the direction a screen has to err in. It is
      !     not the rigorous Cauchy-Schwarz bound -- that needs `(∇i j|∇i j)`,
      !     and libcint's `int2e_ip1ip2` is not among the entry points bound
      !     here -- so the tolerance is set two decades tighter than the
      !     accuracy actually wanted, and the result is checked against an
      !     unscreened gradient in the tests.
      !   * **The density the quartet multiplies.** Blocked to shells, the same
      !     way `shell_density_max` does it next door. A quartet whose four
      !     density blocks are negligible contributes nothing however large its
      !     integral, and near a converged density most of them are.
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
      ! matrix: the inner update is two scattered writes per integral, and
      ! making those atomic costs more than the integral. The price is
      ! 2 * nao^2 * 3 doubles per thread, which is worth watching on a few
      ! thousand functions -- the same caveat the direct Fock build carries.
      !
      ! schedule(dynamic) because the l <= k triangle makes the work per ish
      ! uneven, and a static split leaves threads idle on the tail.
      !$omp parallel default(none) &
      !$omp    shared(mol, tab, density, exchange_from, opt, vj, vk, mx, nao, nbas, &
      !$omp           dims, offs, bq, bra_bound, dsh, esh, tol) &
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

                  ! The largest density element this quartet can reach. The
                  ! four entries are exactly the four the inner loop reads:
                  ! two Coulomb, from the ket pair both ways round, and two
                  ! exchange, pairing the bra's second index with each of the
                  ! ket's.
                  est = bra*bq(ksh, lsh) &
                        *max(dsh(lsh, ksh), dsh(ksh, lsh), esh(jsh, ksh), esh(jsh, lsh))
                  if (est < tol) cycle

                  shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]

                  if (mol%cartesian) then
                     ret = libcint_2e_ip1_cart(buf, shls, mol%atm, mol%natm, &
                                               tab%bas, nbas, tab%env, opt)
                  else
                     ret = libcint_2e_ip1_sph(buf, shls, mol%atm, mol%natm, &
                                              tab%bas, nbas, tab%env, opt)
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
      vhf = -(vj - k_scale*vk)

      deallocate (bounds, bq, bra_bound, dsh, esh, dims, offs)
   end subroutine two_electron_deriv

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
      !! near-linearly-dependent by construction, so a `J^-1` that kept
      !! different modes would be differentiating a slightly different energy
      !! than the one that converged -- which shows up as a gradient that
      !! disagrees with finite differences by a little, in a way that looks
      !! like a factor being wrong somewhere.
      type(libcint_molecule_t), intent(in) :: orb, aux
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
         !! density `Omega_PQ` in this routine -- the collision compiles in some
         !! orderings and means the wrong thing in all of them.
         !!
         !! A
         !! range-separated functional calls this routine twice: once
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

      ! The auxiliary basis sits on the same nuclei as the orbital one, and the
      ! per-atom accumulation below indexes both by the same atom number. A
      ! mismatch would silently attribute an auxiliary function's derivative to
      ! the wrong nucleus, which translational invariance would catch but only
      ! after the fact.
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
      ! The exchange fraction rides the channel weight rather than being
      ! applied afterwards, because both of exchange's shapes -- the one in
      ! `gamma` and the one in `omega` -- are linear in it and
      ! `add_exchange_channel` builds them together. A pure functional passes
      ! zero and the whole exchange assembly is skipped, which is most of the
      ! work in this routine rather than a small saving.
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

      ! Skipped outright at zero rather than multiplied by it. The header a
      ! few lines up claimed this and the code did not: a pure functional was
      ! paying the whole `e^P_ij` transform, the `J^-1 e` solve and the `Z^P`
      ! back-transform to add exactly nothing.
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
            ! Orbital indices. Twice: ip1 differentiates the first index only,
            ! and the second index contributes the transpose -- which is the
            ! same number because both (uv|P) and Gamma^P are symmetric in u
            ! and v. Minus, because libcint's nabla is on the bra.
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
      !! metric term takes half of it, which is the factor that makes the
      !! restricted and unrestricted expressions agree when the two spins are
      !! identical -- and the one worth checking by feeding a closed-shell
      !! system through the unrestricted path.
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
      !! `which` selects whose centre is differentiated: 1 for the first
      !! orbital index, through int3c2e_ip1, and 2 for the auxiliary one,
      !! through int3c2e_ip2. The second orbital index has no entry point of
      !! its own -- (mu nu|P) is symmetric in mu and nu, so its derivative is
      !! the transpose of the first, which is why nothing here integrates for
      !! it.
      !!
      !! No symmetry in the shell-pair loop, unlike the energy's `three_centre`
      !! next door: ip1 differentiates mu and not nu, so the block for (ish,
      !! jsh) is not the transpose of (jsh, ish).
      type(libcint_molecule_t), intent(in) :: orb, aux
      integer, intent(in) :: which
      real(dp), allocatable, intent(out) :: deriv(:, :, :, :)
      real(dp), intent(in), optional :: omega
         !! Attenuate to `erf(omega r)/r`, matching `three_centre`. A fitted
         !! range-separated gradient differentiates two fits, and the long-range
         !! one has to be differentiated with the kernel it was made from.

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
      ! own block of `deriv` and reads nothing another writes, so the only
      ! shared state is the output array and no reduction is needed. `buf` is
      ! allocated inside the region because libcint writes into it per quartet,
      ! which makes it the one genuinely private object here.
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
      type(libcint_molecule_t), intent(in) :: aux
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

end module mqc_libcint_gradient
