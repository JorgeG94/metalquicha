!! What is read off a nuclear magnetic shielding tensor
module mqc_nmr_summary
   !! Kept apart from the code that computes the tensor so the output writer
   !! can summarise one without reaching into a backend.
   use pic_types, only: dp
   use mqc_physical_constants, only: PI
   implicit none
   private

   public :: shielding_summary

contains

   pure subroutine shielding_summary(tensor, isotropic, anisotropy)
      !! The two numbers a spectrum is read from
      !!
      !! The isotropic shielding is a third of the trace. The anisotropy is
      !! taken from the principal values of the symmetrised tensor, ordered so
      !! that `s33` is the one furthest from the isotropic value:
      !! `s33 - (s11 + s22)/2`.
      real(dp), intent(in) :: tensor(3, 3)   !! ppm
      real(dp), intent(out) :: isotropic, anisotropy

      real(dp) :: s(3, 3), p(3), shifted(3, 3)
      real(dp) :: q, r, phi, b, tmp
      integer :: i, far

      isotropic = (tensor(1, 1) + tensor(2, 2) + tensor(3, 3))/3.0_dp
      s = 0.5_dp*(tensor + transpose(tensor))
      ! Eigenvalues of a symmetric 3x3 in closed form.
      q = s(1, 2)**2 + s(1, 3)**2 + s(2, 3)**2
      if (q < tiny(1.0_dp)) then
         p = [s(1, 1), s(2, 2), s(3, 3)]
      else
         shifted = s
         do i = 1, 3
            shifted(i, i) = s(i, i) - isotropic
         end do
         b = sqrt((shifted(1, 1)**2 + shifted(2, 2)**2 + shifted(3, 3)**2 + 2.0_dp*q)/6.0_dp)
         shifted = shifted/b
         r = 0.5_dp*(shifted(1, 1)*(shifted(2, 2)*shifted(3, 3) - shifted(2, 3)**2) &
                     - shifted(1, 2)*(shifted(1, 2)*shifted(3, 3) - shifted(2, 3)*shifted(1, 3)) &
                     + shifted(1, 3)*(shifted(1, 2)*shifted(2, 3) - shifted(2, 2)*shifted(1, 3)))
         r = max(-1.0_dp, min(1.0_dp, r))
         phi = acos(r)/3.0_dp
         p(1) = isotropic + 2.0_dp*b*cos(phi)
         p(3) = isotropic + 2.0_dp*b*cos(phi + 2.0_dp*PI/3.0_dp)
         p(2) = 3.0_dp*isotropic - p(1) - p(3)
      end if
      far = maxloc(abs(p - isotropic), dim=1)
      tmp = p(far)
      p(far) = p(3)
      p(3) = tmp
      anisotropy = p(3) - 0.5_dp*(p(1) + p(2))
   end subroutine shielding_summary

end module mqc_nmr_summary
