!! Cartesian-to-spherical transform coefficients, taken from libcint
module mqc_libcint_ao_data
   !! libcint's own `g_trans_cart2sph` table, for angular momenta up to l = 4.
   !!
   !! Transcribed rather than re-derived, because the transform has to be
   !! *libcint's* and not merely a correct one: the integrals are built with these
   !! coefficients, and a basis function evaluated with a different-but-valid
   !! convention is a different basis. A permutation or a sign here converges just
   !! as prettily and gives a wrong energy.
   !!
   !! Three things make the C table unsafe to read casually, all of which cost an
   !! attempt before this was right:
   !!
   !!   * it carries **two p orderings** behind `#ifdef PYPZPX`, so a text scrape
   !!     picks up both and every later block shifts;
   !!   * the literals are bare integers where they are whole (`1`, `0`), which a
   !!     float-shaped pattern skips;
   !!   * **s and p normalisation is factored out** of the table into
   !!     `CINTcommon_fac_sp`, so the 1s and 0s below are not the whole transform.
   !!
   !! Validated on extraction: 19176 coefficients total, matching the
   !! sum over l of (2l+1)(l+1)(l+2)/2; l = 0 is 1; l = 1 is the px/py/pz
   !! identity; l = 2 reproduces the d solid harmonics derived independently.
   use pic_types, only: dp
   implicit none
   private

   public :: C2S_LMAX, c2s_block, common_fac_sp

   !! Highest angular momentum this table covers. Beyond it the caller must be
   !! refused rather than silently given wrong functions.
   integer, parameter :: C2S_LMAX = 4

   !! Row offsets into `C2S`, one per l. Each block is (2l+1) by (l+1)(l+2)/2,
   !! spherical row-major over Cartesian columns.
   integer, parameter :: N_OFFSETS = C2S_LMAX + 2
   integer, parameter :: C2S_OFFSET(0:N_OFFSETS - 1) = [0, 1, 10, 40, 110, 245]

   !! Total coefficients over l = 0..C2S_LMAX, i.e. sum of (2l+1)(l+1)(l+2)/2.
   integer, parameter :: N_C2S = 245

   real(dp), parameter :: C2S(N_C2S) = [ &
                          1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 1.0_dp, 0.0_dp, 1.0925484305920792_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          1.0925484305920792_dp, 0.0_dp, -0.31539156525252_dp, 0.0_dp, &
                          0.0_dp, -0.31539156525252_dp, 0.0_dp, 0.63078313050504_dp, &
                          0.0_dp, 0.0_dp, 1.0925484305920792_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.5462742152960396_dp, 0.0_dp, &
                          0.0_dp, -0.5462742152960396_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 1.7701307697799304_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, -0.5900435899266435_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 2.8906114426405543_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, -0.4570457994644657_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, -0.4570457994644657_dp, 0.0_dp, &
                          1.8281831978578629_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          -1.1195289977703462_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, -1.1195289977703462_dp, 0.0_dp, 0.7463526651802308_dp, &
                          -0.4570457994644657_dp, 0.0_dp, 0.0_dp, -0.4570457994644657_dp, &
                          0.0_dp, 1.8281831978578629_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          1.4453057213202771_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, -1.4453057213202771_dp, 0.0_dp, 0.0_dp, &
                          0.5900435899266435_dp, 0.0_dp, 0.0_dp, -1.7701307697799304_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 2.5033429417967046_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          -2.5033429417967046_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 5.310392309339791_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          -1.7701307697799304_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, -0.94617469575756_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, -0.94617469575756_dp, 0.0_dp, &
                          5.6770481745453605_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, -2.0071396306718676_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, -2.0071396306718676_dp, 0.0_dp, &
                          2.676186174229157_dp, 0.0_dp, 0.31735664074561293_dp, 0.0_dp, &
                          0.0_dp, 0.6347132814912259_dp, 0.0_dp, -2.5388531259649034_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.31735664074561293_dp, 0.0_dp, -2.5388531259649034_dp, 0.0_dp, &
                          0.8462843753216345_dp, 0.0_dp, 0.0_dp, -2.0071396306718676_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          -2.0071396306718676_dp, 0.0_dp, 2.676186174229157_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          -0.47308734787878_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 2.8385240872726802_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.47308734787878_dp, 0.0_dp, &
                          -2.8385240872726802_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 1.7701307697799304_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, -5.310392309339791_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.6258357354491761_dp, 0.0_dp, &
                          0.0_dp, -3.755014412695057_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.6258357354491761_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp]

contains

   pure function common_fac_sp(l) result(fac)
      !! The normalisation libcint keeps outside the transform table
      !!
      !! `CINTcommon_fac_sp` in g1e.c. Y_00 and Y_1m carry their constant here
      !! rather than in `C2S`, which is why the s and p blocks above are bare
      !! ones and zeros.
      integer, intent(in) :: l
      real(dp) :: fac

      select case (l)
      case (0)
         fac = 0.282094791773878143_dp
      case (1)
         fac = 0.488602511902919921_dp
      case default
         fac = 1.0_dp
      end select
   end function common_fac_sp

   pure subroutine c2s_block(l, block)
      !! The (2l+1) by n_cart transform for one angular momentum
      integer, intent(in) :: l
      real(dp), intent(out) :: block(:, :)

      integer :: n_sph, n_cart, m, c

      n_sph = 2*l + 1
      n_cart = (l + 1)*(l + 2)/2
      do c = 1, n_cart
         do m = 1, n_sph
            ! The C table is row-major over Cartesian columns.
            block(m, c) = C2S(C2S_OFFSET(l) + (m - 1)*n_cart + c)
         end do
      end do
   end subroutine c2s_block

end module mqc_libcint_ao_data
