!! Second derivatives of the atomic orbitals, before anything is built on them
program check_ao_deriv2
   !! Prints chi, its gradient and its six unique second derivatives at a few
   !! scattered points, for comparison against PySCF's
   !! `mol.eval_gto('GTOval_sph_deriv2', coords)`.
   !!
   !! Checked at this level on purpose, and the reason is the same one that
   !! made the three-centre derivative integrals worth checking on their own: a
   !! normalisation or component-ordering mistake here does not announce itself
   !! in an assembled GGA gradient. It produces a number of the right magnitude
   !! that disagrees with finite differences for reasons that could be anywhere
   !! in the exchange-correlation assembly.
   !!
   !! cc-pVDZ on water, so d functions and therefore the spherical transform are
   !! genuinely exercised -- an s/p-only basis leaves the transform an identity
   !! and tests nothing about it. The points are scattered and off-axis for the
   !! same reason: a point on a symmetry axis makes whole components vanish.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_ao, only: eval_ao_block, AO_HESS_COMP
   implicit none

   integer, parameter :: N_DIM = 3
   integer, parameter :: N_POINT = 3
   type(czt_molecule_t) :: mol
   type(error_t) :: error
   real(dp), allocatable :: ao(:, :), grad(:, :, :), hess(:, :, :)
   real(dp) :: points(N_DIM, N_POINT)
   integer :: ig, mu
   character(len=4), parameter :: LABEL(AO_HESS_COMP) = ["xx  ", "xy  ", "xz  ", "yy  ", "yz  ", "zz  "]
   integer :: ih

   !! Which (j, k) each packed second-derivative component stands for. Must
   !! match the packing `eval_ao_block` documents: xx, xy, xz, yy, yz, zz.
   integer, parameter :: HESS_J(AO_HESS_COMP) = [1, 1, 1, 2, 2, 3]
   integer, parameter :: HESS_K(AO_HESS_COMP) = [1, 2, 3, 2, 3, 3]

   call build_czt_molecule([8, 1, 1], ["O", "H", "H"], &
                           reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                    0.0_dp, -1.4308_dp, 1.1078_dp, &
                                    0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                           "cc-pvdz", mol, error)
   if (error%has_error()) then
      write (*, "(a,a)") "basis failed: ", error%get_message()
      error stop 1
   end if

   points(:, 1) = [0.31_dp, -0.47_dp, 0.63_dp]
   points(:, 2) = [-0.85_dp, 1.02_dp, -0.29_dp]
   points(:, 3) = [1.73_dp, 0.58_dp, 2.11_dp]

   call eval_ao_block(mol, points, ao, error, grad=grad, hess=hess)
   if (error%has_error()) then
      write (*, "(a,a)") "evaluation failed: ", error%get_message()
      error stop 1
   end if

   write (*, "(a,i0,a,i0)") "nao = ", mol%nao, "   points = ", N_POINT

   ! A handful of basis functions rather than all of them: enough to cover an
   ! s, a p and a d without printing a wall of numbers nobody reads.
   do ig = 1, N_POINT
      write (*, "(a,i0,a,3f10.5)") "point ", ig, "  at ", points(:, ig)
      do mu = 1, min(mol%nao, 15), 3
         write (*, "(a,i3,a,f18.12)") "   ao ", mu, "  value ", ao(ig, mu)
         write (*, "(a,3f18.12)") "      grad  ", grad(ig, mu, 1:3)
         do ih = 1, AO_HESS_COMP
            write (*, "(a,a,f18.12)") "      d2", trim(LABEL(ih))//"  ", hess(ig, mu, ih)
         end do
      end do
   end do

   ! Two step sizes, and the ratio is the diagnostic. A central difference
   ! truncates at order h^2, so halving the step should quarter the deviation.
   ! If instead it barely moves, the disagreement is a real error in the
   ! second derivatives rather than the difference formula, and no choice of
   ! threshold would tell the two apart.
   write (*, "(a)") ""
   write (*, "(a)") "== second derivatives against finite differences of the first"
   call self_consistency(mol, points, 2.0e-4_dp)
   call self_consistency(mol, points, 1.0e-4_dp)
   call self_consistency(mol, points, 0.5e-4_dp)

contains

   subroutine self_consistency(mol, points, step)
      !! The second derivatives against finite differences of the first
      !!
      !! This is the check that decides whether the second derivatives are
      !! right, and it is deliberately free of PySCF. Comparing against PySCF
      !! measures two things at once -- whether the second derivatives are
      !! correct, and whether the two codes hold bit-identical basis data --
      !! and the AO *values* already disagree with it by about 2e-7 relative on
      !! this basis, which is the normalisation artefact this repository
      !! documents elsewhere. That floor would mask a real error of the same
      !! size.
      !!
      !! Differencing this program's own gradient has no such floor: any error
      !! in the second derivatives shows up against the first derivatives they
      !! are supposed to be the derivative of, in the same basis, with the same
      !! normalisation, whatever that normalisation happens to be.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: points(:, :)
      real(dp), intent(in) :: step

      real(dp) :: STEP_
      real(dp), allocatable :: shifted(:, :)
      real(dp), allocatable :: ao_p(:, :), ao_m(:, :)
      real(dp), allocatable :: grad_p(:, :, :), grad_m(:, :, :)
      real(dp) :: numeric, worst
      type(error_t) :: err
      integer :: k, ip, imu, ihc

      STEP_ = step
      worst = 0.0_dp
      allocate (shifted(size(points, 1), size(points, 2)))

      do k = 1, 3
         shifted = points
         shifted(k, :) = points(k, :) + STEP_
         call eval_ao_block(mol, shifted, ao_p, err, grad=grad_p)
         if (err%has_error()) then
            write (*, "(a)") "FAIL: displaced evaluation failed"
            error stop 1
         end if

         shifted = points
         shifted(k, :) = points(k, :) - STEP_
         call eval_ao_block(mol, shifted, ao_m, err, grad=grad_m)
         if (err%has_error()) then
            write (*, "(a)") "FAIL: displaced evaluation failed"
            error stop 1
         end if

         ! d/dx_k of d chi/d x_j is the (j, k) second derivative. Every packed
         ! component that mentions k gets tested this way, and the symmetric
         ! ones get tested twice from the two directions -- which is a check in
         ! itself, since a transposed index would disagree between them.
         do ihc = 1, AO_HESS_COMP
            if (HESS_J(ihc) /= k .and. HESS_K(ihc) /= k) cycle
            do ip = 1, size(points, 2)
               do imu = 1, mol%nao
                  if (HESS_J(ihc) == k) then
                     numeric = (grad_p(ip, imu, HESS_K(ihc)) &
                                - grad_m(ip, imu, HESS_K(ihc)))/(2.0_dp*STEP_)
                  else
                     numeric = (grad_p(ip, imu, HESS_J(ihc)) &
                                - grad_m(ip, imu, HESS_J(ihc)))/(2.0_dp*STEP_)
                  end if
                  worst = max(worst, abs(numeric - hess(ip, imu, ihc)))
               end do
            end do
         end do
         deallocate (ao_p, ao_m, grad_p, grad_m)
      end do

      write (*, "(a,es9.2,a,es12.4)") "   step ", STEP_, "   worst absolute deviation ", worst
      if (worst > 1.0e-5_dp) then
         write (*, "(a)") "   FAIL: the second derivatives do not differentiate the first"
         error stop 1
      end if
   end subroutine self_consistency

end program check_ao_deriv2
