!! The analytic Hessian, for a build that cannot reach the integrals it needs
module mqc_libcint_hessian
   !! Stand-in compiled instead of `mqc_libcint_hessian.f90` when the CPU
   !! backend is libcint rather than libfint.
   !!
   !! **Not a refusal, an absence.** The second derivatives of the two-electron
   !! integrals come from libfint's native modules -- `cint_gen_hess` and the
   !! workspace type beside it -- and libcint exposes no such entry points
   !! through the compatibility layer this project links against. Rather than
   !! declare thirty `bind(C)` interfaces for a path that build does not need,
   !! `analytic_hessian_available` reports false and the caller does what it
   !! already does for an unrestricted or density-fitted reference: computes the
   !! Hessian by finite differences of analytic gradients, which is correct and
   !! merely slower.
   !!
   !! The two procedures below exist so `mqc_libcint_bridge` compiles against
   !! one module name either way. They are unreachable: the bridge tests
   !! availability before calling them, and they say so loudly rather than
   !! returning something plausible if that test is ever dropped.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t
   implicit none
   private

   public :: rhf_hessian
   public :: hessian_to_matrix
   public :: analytic_hessian_available

contains

   pure function analytic_hessian_available() result(available)
      !! Whether this build has an analytic Hessian at all
      logical :: available

      available = .false.
   end function analytic_hessian_available

   subroutine rhf_hessian(mol, atomic_numbers, density, orbitals, energies, n_occ, &
                          hess, error, max_iter, tol)
      !! Unreachable in this build; see the module note
      type(libcint_molecule_t), intent(in), target :: mol
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: hess(:, :, :, :)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol

      call error%set(ERROR_VALIDATION, "this build has no analytic Hessian: the "// &
                     "second-derivative integrals are libfint's and this is a "// &
                     "libcint build. Callers are meant to ask "// &
                     "analytic_hessian_available first and fall back to finite "// &
                     "differences, so reaching this is a bug in the caller.")
      ! Referenced so the arguments are used in every build, and so this cannot
      ! be mistaken for a routine that computed something.
      if (mol%natm < 0 .or. size(atomic_numbers) < 0 .or. size(density) < 0 &
          .or. size(orbitals) < 0 .or. size(energies) < 0 .or. n_occ < 0) return
      if (present(max_iter) .or. present(tol)) return
   end subroutine rhf_hessian

   subroutine hessian_to_matrix(hess, matrix)
      !! Unreachable in this build; see the module note
      real(dp), intent(in) :: hess(:, :, :, :)
      real(dp), allocatable, intent(out) :: matrix(:, :)

      allocate (matrix(3*size(hess, 3), 3*size(hess, 4)))
      matrix = 0.0_dp
   end subroutine hessian_to_matrix

end module mqc_libcint_hessian
