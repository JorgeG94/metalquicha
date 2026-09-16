!! Stand-in for the stability analysis when the build has no OpenTrustRegion
!!
!! Beside the module it replaces rather than in src/methods/stubs: its
!! signature names this backend's own types, so it can only be compiled where
!! they are, and fpm globs src/ without ever compiling backends/.
module mqc_czt_stability
   !! Same name and same entry points as the real bridge, declining.
   !!
   !! `otr_available` is how a caller asks in advance, so the refusal can name
   !! the build option -- and so that a test can skip rather than fail.
   !!
   !! **This is not the stability analysis being unavailable.** Neither the
   !! Hessian nor the eigensolver that diagonalises it is optional:
   !! `mqc_czt_ov_hessian` needs nothing from OpenTrustRegion and
   !! `mqc_czt_native_stability` diagonalises it with this program's own
   !! Davidson, in every build. What is missing here is the *second*
   !! implementation, the one `keywords.scf.stability_engine: otr` asks for and
   !! the one the native path is cross-checked against.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_xc, only: xc_context_t
   use mqc_czt_ov_hessian, only: ov_hessian_t, stability_result_t
   implicit none
   private

   public :: scf_stability
   public :: stability_of_hessian
   public :: otr_available

   character(len=*), parameter :: REFUSAL = &
                                  "keywords.scf.stability_engine 'otr' needs OpenTrustRegion; build "// &
                                  "with -DMQC_ENABLE_OTR=ON, or leave the engine at 'native', which "// &
                                  "needs nothing"

contains

   pure function otr_available() result(available)
      !! No OpenTrustRegion means no eigensolver for the electronic Hessian.
      logical :: available

      available = .false.
   end function otr_available

   subroutine scf_stability(mol, orbitals, energies, n_occ, result, error, xc, &
                            reference, k_scale, rs_k_lr, rs_omega, bounds, &
                            conv_tol, max_iter, n_trial_vectors, seed)
      !! No-op stand-in: report the missing backend, compute nothing
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      type(stability_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional, target :: xc
      real(dp), intent(in), optional :: reference(:, :)
      real(dp), intent(in), optional :: k_scale
      real(dp), intent(in), optional :: rs_k_lr, rs_omega
      real(dp), intent(in), optional :: bounds(:, :)
      real(dp), intent(in), optional :: conv_tol
      integer, intent(in), optional :: max_iter
      integer, intent(in), optional :: n_trial_vectors
      integer, intent(in), optional :: seed

      result%ran = .false.
      call error%set(ERROR_VALIDATION, REFUSAL)

      ! The arguments are the contract, so they are referenced rather than
      ! renamed away; this is the shape every stub in this directory uses.
      associate (unused_mol => mol, unused_occ => n_occ)
      end associate
      if (size(orbitals) < 0 .or. size(energies) < 0) return
      if (present(xc) .and. present(reference)) return
      if (present(k_scale) .or. present(rs_k_lr) .or. present(rs_omega)) return
      if (present(bounds) .or. present(conv_tol)) return
      if (present(max_iter) .or. present(n_trial_vectors) .or. present(seed)) return
   end subroutine scf_stability

   subroutine stability_of_hessian(hessian, result, error, conv_tol, max_iter, &
                                   n_trial_vectors, seed)
      !! No-op stand-in: the Hessian can be applied, but not diagonalised
      type(ov_hessian_t), intent(inout), target :: hessian
      type(stability_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: conv_tol
      integer, intent(in), optional :: max_iter
      integer, intent(in), optional :: n_trial_vectors
      integer, intent(in), optional :: seed

      result%ran = .false.
      result%n_parameters = hessian%length()
      call error%set(ERROR_VALIDATION, REFUSAL)

      if (present(conv_tol)) return
      if (present(max_iter) .or. present(n_trial_vectors) .or. present(seed)) return
   end subroutine stability_of_hessian

end module mqc_czt_stability
