!! Contraction coefficient normalization for Gaussian basis shells
module mqc_basis_normalization
   !! Normalizes contracted Gaussian shell coefficients.
   !!
   !! Basis set files store coefficients for *unnormalized* primitives, but
   !! integral engines expect the contracted function to be normalized. This
   !! applies the standard two-stage normalization: a per-primitive factor,
   !! followed by an overall factor Q that normalizes the contracted function.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cgto, only: cgto_type
   use mqc_physical_constants, only: PI
   implicit none
   private

   public :: normalize_shell_coefficients  !! Normalize one shell's coefficients
   public :: normalized_coefficients       !! Normalize raw exponent/coefficient arrays

   integer, parameter :: MAX_ANG_MOM = 8  !! Highest L the closed form is valid for

contains

   subroutine normalized_coefficients(ang_mom, exponents, coefficients, &
                                      normalization, normalized, error)
      !! Normalized contraction coefficients for a single shell
      !!
      !! Stage 1 scales each primitive by its own normalization constant.
      !! Stage 2 applies the overall factor Q so that the self-overlap of the
      !! contracted function equals `normalization`.
      integer, intent(in) :: ang_mom              !! Angular momentum (0=s, 1=p, ...)
      real(dp), intent(in) :: exponents(:)        !! Primitive exponents, all > 0
      real(dp), intent(in) :: coefficients(:)     !! Raw contraction coefficients
      real(dp), intent(in) :: normalization       !! Target self-overlap, usually 1
      real(dp), intent(out) :: normalized(:)      !! Normalized coefficients
      type(error_t), intent(out) :: error

      real(dp) :: pi_three_halves, two_to_l, double_factorial, overlap, scale
      integer :: n_prim, i, j, k

      n_prim = size(exponents)

      if (any(exponents(1:n_prim) <= 0.0_dp)) then
         call error%set(ERROR_VALIDATION, "Basis normalization: exponents must be positive")
         return
      end if
      if (ang_mom < 0 .or. ang_mom > MAX_ANG_MOM) then
         call error%set(ERROR_VALIDATION, "Basis normalization: unsupported angular momentum")
         return
      end if
      if (normalization <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "Basis normalization: normalization must be positive")
         return
      end if

      pi_three_halves = PI**1.5_dp
      two_to_l = 2.0_dp**real(ang_mom, dp)

      ! (2L-1)!! -- the double factorial appearing in the Gaussian self-overlap
      double_factorial = 1.0_dp
      do k = 1, ang_mom
         double_factorial = double_factorial*real(2*k - 1, dp)
      end do

      ! Stage 1: per-primitive normalization
      do i = 1, n_prim
         normalized(i) = sqrt(two_to_l/(pi_three_halves*double_factorial) &
                              *(2.0_dp*exponents(i))**(real(ang_mom, dp) + 1.5_dp)) &
                         *coefficients(i)
      end do

      ! Stage 2: overall factor from the contracted self-overlap
      overlap = 0.0_dp
      do i = 1, n_prim
         do j = 1, n_prim
            overlap = overlap + (sqrt(4.0_dp*exponents(i)*exponents(j)) &
                                 /(exponents(i) + exponents(j)))**(real(ang_mom, dp) + 1.5_dp) &
                      *coefficients(i)*coefficients(j)
         end do
      end do
      scale = overlap**(-0.5_dp)*sqrt(normalization)

      do i = 1, n_prim
         normalized(i) = normalized(i)*scale
      end do
   end subroutine normalized_coefficients

   subroutine normalize_shell_coefficients(shell, normalized, error)
      !! Normalized coefficients for a `cgto_type` shell
      type(cgto_type), intent(in) :: shell            !! Shell as read from a basis file
      real(dp), allocatable, intent(out) :: normalized(:)  !! Normalized coefficients
      type(error_t), intent(out) :: error

      allocate (normalized(shell%nfunc))
      call normalized_coefficients(shell%ang_mom, &
                                   shell%exponents(1:shell%nfunc), &
                                   shell%coefficients(1:shell%nfunc), &
                                   1.0_dp, normalized, error)
      if (error%has_error()) then
         call error%add_context("mqc_basis_normalization:normalize_shell_coefficients")
      end if
   end subroutine normalize_shell_coefficients

end module mqc_basis_normalization
