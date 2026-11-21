!! The six independent components of a symmetric 3x3 tensor, in one order
module mqc_symmetric_tensor
   !! A polarizability, its derivatives and a quadrupole are symmetric, and
   !! are stored and sent as six numbers. This is the order everything that
   !! does so agrees on: `xx, yy, zz, xy, xz, yz`.
   use pic_types, only: dp
   implicit none
   private

   public :: pack_symmetric
   public :: N_SYMMETRIC, SYMMETRIC_ROW, SYMMETRIC_COL

   integer, parameter :: N_SYMMETRIC = 6
      !! Independent components of a symmetric 3x3 tensor
   integer, parameter :: SYMMETRIC_ROW(N_SYMMETRIC) = [1, 2, 3, 1, 1, 2]
      !! Row of each packed component
   integer, parameter :: SYMMETRIC_COL(N_SYMMETRIC) = [1, 2, 3, 2, 3, 3]
      !! Column of each packed component

contains

   pure function pack_symmetric(tensor) result(packed)
      !! The components `xx, yy, zz, xy, xz, yz` of a 3x3 tensor
      !!
      !! An off-diagonal component is the mean of the two elements it stands
      !! for, so a tensor that is symmetric only to a solver's tolerance is
      !! symmetrised rather than read from one triangle.
      real(dp), intent(in) :: tensor(3, 3)
      real(dp) :: packed(N_SYMMETRIC)

      integer :: k

      do k = 1, N_SYMMETRIC
         packed(k) = 0.5_dp*(tensor(SYMMETRIC_ROW(k), SYMMETRIC_COL(k)) &
                             + tensor(SYMMETRIC_COL(k), SYMMETRIC_ROW(k)))
      end do
   end function pack_symmetric

end module mqc_symmetric_tensor
