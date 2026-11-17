!! The RI-CCSD nuclear gradient over a closed-shell reference
module mqc_czt_ri_ccsd_gradient
   !! An exact-ERI restricted Hartree-Fock reference with a *fitted* CCSD
   !! correlation energy, which is what an `ri-ccsd` deck computes: the
   !! auxiliary basis names the correlation fit and not the SCF.
   !!
   !! The correlation energy is
   !!
   !!     E = sum_pq gamma(p,q) f(p,q) + 1/2 sum_pqrs Gamma(p,q,r,s) (pq|rs)_fit
   !!
   !! with the densities `mqc_czt_rcc_density` builds from the amplitudes and
   !! the Lambda amplitudes. Its gradient needs, of the two-particle density,
   !! only the contraction with the fit, `g3`; this module turns that and
   !! `gamma` into the four intermediates `relaxed_fitted_gradient` takes, and
   !! the relaxation and assembly are the RI-MP2 gradient's own.
   !!
   !! `tools/cpu_validation/ri_ccsd_gradient.py` is the same gradient in numpy,
   !! checked against finite differences of the energy, and is what the test
   !! holds this one to.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, three_centre, two_centre, &
                                metric_inverse_sqrt
   use mqc_czt_cphf, only: IN_CORE_LIMIT
   use mqc_czt_gemm_threads, only: gemm_over_columns
   use mqc_czt_ecp, only: ecp_refuses_derivatives
   use mqc_czt_rcc, only: rcc_state_t
   use mqc_czt_rcc_lambda, only: rcc_lambda_t
   use mqc_czt_rcc_density, only: rccsd_one_particle_density, &
                                  rccsd_fitted_two_particle_density
   use mqc_czt_ri_mp2_gradient, only: relaxation_input_t, relaxed_fitted_gradient
   implicit none
   private

   public :: czt_ri_ccsd_gradient

contains

   subroutine czt_ri_ccsd_gradient(mol, aux, coeff, orbital_energies, n_occ, frozen, &
                                   state, lambda, gradient, error, fit_response)
      !! dE(RI-CCSD)/dR for a closed-shell reference, in Hartree/Bohr
      !!
      !! The total gradient: the reference's and the correlation's together.
      type(czt_molecule_t), intent(in) :: mol
      type(czt_molecule_t), intent(in) :: aux   !! Auxiliary basis, same atoms
      real(dp), intent(in) :: coeff(:, :)       !! Every molecular orbital, frozen ones included
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      integer, intent(in) :: frozen             !! Core orbitals left out of the correlation
      type(rcc_state_t), intent(in) :: state    !! The converged amplitudes
      type(rcc_lambda_t), intent(in) :: lambda
      real(dp), allocatable, intent(out) :: gradient(:, :)   !! (3, n_atoms)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: fit_response
         !! Build the Z-vector's operator and the reference potentials from the
         !! fit. Absent is `.true.`, as in the RI-MP2 gradient and for its
         !! reason; `.false.` is the exact operator, which is what this energy's
         !! derivative strictly is.

      type(relaxation_input_t), target :: inp
      real(dp), allocatable :: metric(:, :), b_all(:, :, :), g3(:, :, :), gamma(:, :)
      real(dp), allocatable :: block(:, :), c_act(:, :), g_ao(:, :), overlap_lm(:, :)
      real(dp), allocatable :: half(:, :), gflat(:, :), bflat(:, :)
      integer :: n_ao, n_mo, n_act, n_aux, no, nv, k

      if (ecp_refuses_derivatives(mol%core_electrons, "RI-CCSD gradient", error)) return

      n_ao = mol%nao
      n_mo = size(coeff, 2)
      no = state%no
      nv = state%nv
      n_act = no + nv
      if (frozen + no /= n_occ .or. frozen + n_act /= n_mo) then
         call error%set(ERROR_VALIDATION, "the RI-CCSD gradient was handed amplitudes "// &
                        "over a different orbital space from its orbitals")
         return
      end if
      if (aux%natm /= mol%natm) then
         call error%set(ERROR_VALIDATION, "the RI-CCSD gradient's auxiliary basis is "// &
                        "on a different set of atoms from the orbital basis")
         return
      end if
      if (any(abs(aux%coords - mol%coords) > 1.0e-10_dp)) then
         call error%set(ERROR_VALIDATION, "the RI-CCSD gradient's auxiliary basis is "// &
                        "not on the orbital basis's atoms, in their order")
         return
      end if

      inp%frozen = frozen
      inp%fitted_reference = .false.
      inp%fit_response = .true.
      if (present(fit_response)) inp%fit_response = fit_response
      inp%dense = real(n_ao, dp)**4*8.0_dp <= IN_CORE_LIMIT .and. .not. inp%fit_response

      ! ---- the fit, over every molecular orbital ---------------------------
      ! Built here, with the symmetric factor, rather than taken from the
      ! energy calculation: that one used whichever factor was cheapest, and
      ! the densities below have to be in the factorisation the derivative
      ! integrals are contracted in.
      call three_centre(mol, aux, inp%three_ao)
      call two_centre(aux, metric)
      call metric_inverse_sqrt(metric, inp%jm12, error)
      if (error%has_error()) return
      deallocate (metric)
      n_aux = size(inp%jm12, 2)
      allocate (inp%bref(n_ao*n_ao, n_aux))
      call gemm_over_columns(inp%three_ao, inp%jm12, inp%bref)

      allocate (b_all(n_mo, n_mo, n_aux), block(n_ao, n_ao))
      do k = 1, n_aux
         block = reshape(inp%bref(:, k), [n_ao, n_ao])
         b_all(:, :, k) = matmul(transpose(coeff), matmul(block, coeff))
      end do
      deallocate (block)

      ! ---- the densities ---------------------------------------------------
      call rccsd_one_particle_density(state, lambda, gamma)
      call rccsd_fitted_two_particle_density(state, lambda, &
                                             b_all(frozen + 1:n_mo, frozen + 1:n_mo, :), &
                                             g3, error)
      if (error%has_error()) return

      ! The unrelaxed density, without its occupied-virtual block: the
      ! relaxation solves for the whole of that block, so whatever the
      ! amplitudes put there would be counted twice.
      allocate (inp%dm1mo(n_mo, n_mo), source=0.0_dp)
      inp%dm1mo(frozen + 1:n_occ, frozen + 1:n_occ) = gamma(1:no, 1:no)
      inp%dm1mo(n_occ + 1:n_mo, n_occ + 1:n_mo) = gamma(no + 1:n_act, no + 1:n_act)

      ! imat(t, p) = -sum_q sum_P b(t, q, P) g3(p, q, P), `t` over every
      ! orbital and `p`, `q` over the active ones.
      allocate (inp%imat(n_mo, n_mo), source=0.0_dp)
      do k = 1, n_aux
         call pic_gemm(b_all(:, frozen + 1:n_mo, k), g3(:, :, k), &
                       inp%imat(:, frozen + 1:n_mo), transb="T", alpha=-1.0_dp, beta=1.0_dp)
      end do

      ! A quarter of the density contracted twice with the fit, and once and
      ! carried to the atomic orbitals. Both leave the factor's own index for
      ! the auxiliary function's.
      allocate (gflat(n_act*n_act, n_aux), bflat(n_act*n_act, n_aux))
      gflat = reshape(g3, [n_act*n_act, n_aux])
      bflat = reshape(b_all(frozen + 1:n_mo, frozen + 1:n_mo, :), [n_act*n_act, n_aux])
      deallocate (b_all)
      allocate (overlap_lm(n_aux, n_aux), half(n_aux, n_aux), inp%gamma_pq(n_aux, n_aux))
      call pic_gemm(gflat, bflat, overlap_lm, transa="T")
      call pic_gemm(inp%jm12, overlap_lm, half)
      call pic_gemm(half, inp%jm12, inp%gamma_pq, transb="T")
      inp%gamma_pq = 0.25_dp*inp%gamma_pq
      deallocate (gflat, bflat, overlap_lm, half)

      allocate (c_act(n_ao, n_act), g_ao(n_ao*n_ao, n_aux), block(n_ao, n_ao))
      c_act = coeff(:, frozen + 1:n_mo)
      do k = 1, n_aux
         block = matmul(c_act, matmul(g3(:, :, k), transpose(c_act)))
         g_ao(:, k) = reshape(block, [n_ao*n_ao])
      end do
      deallocate (g3, block, c_act)
      allocate (half(n_ao*n_ao, n_aux))
      call gemm_over_columns(g_ao, inp%jm12, half, transb="T")
      deallocate (g_ao)
      allocate (inp%gamma_ao(n_ao, n_ao, n_aux))
      inp%gamma_ao = 0.25_dp*reshape(half, [n_ao, n_ao, n_aux])
      deallocate (half)

      call relaxed_fitted_gradient(mol, aux, coeff, orbital_energies, n_occ, inp, &
                                   gradient, error)
   end subroutine czt_ri_ccsd_gradient

end module mqc_czt_ri_ccsd_gradient
