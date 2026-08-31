!! Energy-based DIIS coefficient solvers: EDIIS and ADIIS
module mqc_ediis
   !! What DIIS cannot do, and why these exist.
   !!
   !! DIIS minimises the norm of a combination of *error* vectors and knows
   !! nothing about the energy. Its coefficients are unconstrained apart from
   !! summing to one, so they go negative and large and the extrapolated Fock
   !! matrix can lie far outside the densities it was built from. Near the
   !! solution that freedom is exactly what makes it fast. Far from it, the same
   !! freedom is what lets an SCF oscillate, stall, or settle onto a saddle
   !! point with a perfectly small residual.
   !!
   !! EDIIS and ADIIS replace the residual with an energy model and constrain
   !! the coefficients to the simplex -- non-negative and summing to one -- so
   !! the result is an *interpolation* between densities already visited. That
   !! is what makes them descend from a bad starting point, and also what makes
   !! them slow at the end: the true solution is generally outside the hull of
   !! what has been seen. So neither replaces DIIS. The pair is the method:
   !! energy-based while the error is large, DIIS once it is small.
   !!
   !! **EDIIS** (Kudin, Scuseria, Cances, JCP 116, 8255 (2002)):
   !!
   !!     E(c) = sum_i c_i E_i - 1/2 sum_ij c_i c_j <D_i - D_j | F_i - F_j>
   !!
   !! **ADIIS** (Hu, Yang, JCP 132, 054109 (2010)), built on the augmented
   !! Roothaan-Hall energy function, measured against the newest entry `n`:
   !!
   !!     E(c) = 2 sum_i c_i <D_i - D_n | F_n>
   !!            + sum_ij c_i c_j <D_i - D_n | F_j - F_n>
   !!
   !! Both reduce to a quadratic in `c` once the traces are collected into
   !! `df(i,j) = <D_i | F_j>`, which is the only thing either needs from the
   !! history -- so the caller keeps the matrices and this module never sees
   !! one. `df` is *not* symmetric: `<D_i|F_j>` and `<D_j|F_i>` differ off the
   !! diagonal, and the expansions below rely on that being carried through.
   !!
   !! Everything here is in age order, oldest first, matching what
   !! `diis_coefficients` returns, so a caller walks one history in one order.
   use pic_types, only: dp
   implicit none
   private

   public :: ediis_coefficients
   public :: adiis_coefficients
   public :: simplex_quadratic_min

   integer, parameter :: MAX_STEPS = 200
      !! BFGS iterations. PySCF stops at `tol = 1e-9`; a handful of variables
      !! reaches that in well under this, and the bound is a runaway guard.
      !! Iterations of the simplex minimisation. The subspace is a handful of
      !! vectors, so this is small work beside one Fock build; the bound exists
      !! to stop a pathological case spinning, not to ration effort.
   real(dp), parameter :: GRAD_FLOOR = 1.0e-10_dp
      !! Projected-gradient norm below which the minimum is taken as found.
      !! PySCF asks scipy for `tol = 1e-9`; this is the same order.
   real(dp), parameter :: STEP_FLOOR = 1.0e-14_dp
      !! Backtracking gives up below this and keeps the best point so far

contains

   subroutine ediis_coefficients(energies, df, n, coefficients, ok)
      !! EDIIS weights, oldest first
      !!
      !!     E(c) = sum_i c_i E_i - sum_ij c_i c_j <D_i - D_j | F_i - F_j>
      !!
      !! with the bracket expanded from `df` as
      !! `df_ii + df_jj - df_ij - df_ji`, so this is `b . c + c^T A c` with
      !! `b = E` and `A = -(df_ii + df_jj - df_ij - df_ji)`.
      !!
      !! **The published functional carries a 1/2 on that sum and this does
      !! not**, because `pyscf/scf/diis.py:ediis_minimize` does not: its
      !! `costf` is `c.es - c.df.c` against a `df` already symmetrised into the
      !! bracket above. The factor is not cosmetic -- it reweights the energies
      !! against the correction and moves where the minimum sits -- so one of
      !! the two had to be chosen rather than averaged. PySCF is chosen for the
      !! reason the VV10 port gives for following it term for term: it makes a
      !! disagreement a porting bug rather than a difference of convention, and
      !! there is a reference implementation to diff against.
      !!
      !! The quadratic term is *not* positive definite, so the minimisation is
      !! a descent method on the simplex rather than a solve.
      real(dp), intent(in) :: energies(:)   !! (n) total energy of each entry
      real(dp), intent(in) :: df(:, :)      !! (n, n) <D_i | F_j>, age order
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: coefficients(:)
      logical, intent(out) :: ok            !! False when the subspace is too small

      real(dp), allocatable :: a(:, :), b(:)
      integer :: i, j

      ok = .false.
      if (n < 2) return

      allocate (a(n, n), b(n))
      do i = 1, n
         b(i) = energies(i)
         do j = 1, n
            a(i, j) = -(df(i, i) + df(j, j) - df(i, j) - df(j, i))
         end do
      end do

      call simplex_quadratic_min(a, b, n, coefficients)
      deallocate (a, b)
      ok = .true.
   end subroutine ediis_coefficients

   subroutine adiis_coefficients(df, n, coefficients, ok)
      !! ADIIS weights, oldest first; the newest entry is age `n`
      !!
      !! Expanding against the newest entry:
      !!
      !!     E(c) = 2 sum_i c_i (df_in - df_nn)
      !!            + sum_ij c_i c_j (df_ij - df_in - df_nj + df_nn)
      !!
      !! so `b_i = 2 (df_in - df_nn)` and `A_ij = df_ij - df_in - df_nj + df_nn`.
      !! At `c = (0, ..., 1)` -- all weight on the newest -- both vanish, which
      !! is the sanity check on the signs: ADIIS measures how much better the
      !! model says a mixture is than simply keeping what was just built.
      real(dp), intent(in) :: df(:, :)      !! (n, n) <D_i | F_j>, age order
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: coefficients(:)
      logical, intent(out) :: ok

      real(dp), allocatable :: a(:, :), b(:)
      integer :: i, j

      ok = .false.
      if (n < 2) return

      allocate (a(n, n), b(n))
      do i = 1, n
         b(i) = 2.0_dp*(df(i, n) - df(n, n))
         do j = 1, n
            a(i, j) = df(i, j) - df(i, n) - df(n, j) + df(n, n)
         end do
      end do

      call simplex_quadratic_min(a, b, n, coefficients)
      deallocate (a, b)
      ok = .true.
   end subroutine adiis_coefficients

   subroutine simplex_quadratic_min(a, b, n, c)
      !! Minimise `b . c + c^T A c` over the simplex: c_i >= 0, sum(c) = 1
      !!
      !! **The constraint is removed rather than enforced.** Writing
      !! `c_i = x_i^2 / sum_j x_j^2` satisfies both conditions for any real `x`,
      !! which turns a constrained problem into an unconstrained one. This is
      !! the parameterisation `pyscf/scf/diis.py` uses, and the chain rule
      !! through it is the same:
      !!
      !!     dE/dx_n = (2 x_n / S) * (g_n - sum_k c_k g_k),   g = dE/dc
      !!
      !! The bracket is the gradient projected onto the simplex, so no separate
      !! projection appears -- it is what the parameterisation produces.
      !!
      !! **BFGS, not steepest descent.** `A` is not positive definite, and on a
      !! badly conditioned quadratic steepest descent stalls short of the
      !! minimum and returns coefficients that are merely nearby. That is
      !! invisible: the SCF still converges, via DIIS-like behaviour, just
      !! slower -- which is exactly the symptom that would be read as "the
      !! energy schemes are slow" rather than "the inner minimisation is not
      !! converged". PySCF uses BFGS at `tol = 1e-9`; this matches it, with a
      !! backtracking Armijo line search so an indefinite region cannot take an
      !! unbounded step.
      !!
      !! Started from a uniform `x`, as PySCF does. The obvious alternative --
      !! all weight on the newest entry -- is a different basin on a non-convex
      !! model, so it is not a free choice.
      real(dp), intent(in) :: a(:, :)   !! (n, n)
      real(dp), intent(in) :: b(:)      !! (n)
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: c(:)

      real(dp), allocatable :: x(:), gx(:), gx_new(:), p(:), sv(:), yv(:), hy(:)
      real(dp), allocatable :: hmat(:, :), x_new(:), c_new(:)
      real(dp) :: e_now, e_trial, step, sy, rho, slope, gnorm
      integer :: it, i, j, ls

      allocate (x(n), gx(n), gx_new(n), p(n), sv(n), yv(n), hy(n))
      allocate (hmat(n, n), x_new(n), c_new(n), c(n))

      x = 1.0_dp
      call weights(x, n, c)
      e_now = model(a, b, c, n)
      call x_gradient(a, b, c, x, n, gx)

      hmat = 0.0_dp
      do i = 1, n
         hmat(i, i) = 1.0_dp
      end do

      do it = 1, MAX_STEPS
         gnorm = sqrt(sum(gx*gx))
         if (gnorm < GRAD_FLOOR) exit

         p = -matmul(hmat, gx)
         slope = sum(p*gx)
         ! A non-descent direction means the inverse-Hessian estimate has gone
         ! indefinite, which an indefinite `A` permits. Reset to steepest
         ! descent rather than step uphill.
         if (slope >= 0.0_dp) then
            hmat = 0.0_dp
            do i = 1, n
               hmat(i, i) = 1.0_dp
            end do
            p = -gx
            slope = sum(p*gx)
         end if

         step = 1.0_dp
         e_trial = e_now
         do ls = 1, 60
            x_new = x + step*p
            call weights(x_new, n, c_new)
            e_trial = model(a, b, c_new, n)
            if (e_trial <= e_now + 1.0e-4_dp*step*slope) exit
            step = 0.5_dp*step
         end do
         if (step < STEP_FLOOR) exit

         call x_gradient(a, b, c_new, x_new, n, gx_new)
         sv = x_new - x
         yv = gx_new - gx
         sy = sum(sv*yv)

         x = x_new
         c = c_new
         e_now = e_trial
         gx = gx_new

         ! Skip the update when the curvature condition fails; keeping the last
         ! good inverse Hessian is safer than one built from a bad pair.
         if (sy > 1.0e-14_dp) then
            rho = 1.0_dp/sy
            hy = matmul(hmat, yv)
            do j = 1, n
               do i = 1, n
                  hmat(i, j) = hmat(i, j) - rho*(sv(i)*hy(j) + hy(i)*sv(j)) &
                               + rho*rho*sum(yv*hy)*sv(i)*sv(j) + rho*sv(i)*sv(j)
               end do
            end do
         end if
      end do

      deallocate (x, gx, gx_new, p, sv, yv, hy, hmat, x_new, c_new)
   end subroutine simplex_quadratic_min

   pure subroutine x_gradient(a, b, c, x, n, gx)
      !! `dE/dx`, the model gradient chained through the parameterisation
      real(dp), intent(in) :: a(:, :), b(:), c(:), x(:)
      integer, intent(in) :: n
      real(dp), intent(out) :: gx(:)

      real(dp), allocatable :: g(:)
      real(dp) :: s, avg
      integer :: i

      allocate (g(n))
      call model_gradient(a, b, c, n, g)
      avg = sum(c(1:n)*g(1:n))
      s = sum(x(1:n)*x(1:n))
      if (s <= 0.0_dp) s = 1.0_dp
      do i = 1, n
         gx(i) = 2.0_dp*x(i)/s*(g(i) - avg)
      end do
      deallocate (g)
   end subroutine x_gradient

   pure subroutine weights(x, n, c)
      !! The simplex parameterisation, `c_i = x_i^2 / sum_j x_j^2`
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: n
      real(dp), intent(out) :: c(:)
      real(dp) :: s

      s = sum(x(1:n)*x(1:n))
      if (s <= 0.0_dp) then
         c(1:n) = 1.0_dp/real(n, dp)
      else
         c(1:n) = x(1:n)*x(1:n)/s
      end if
   end subroutine weights

   pure function model(a, b, c, n) result(e)
      !! `b . c + c^T A c`
      real(dp), intent(in) :: a(:, :), b(:), c(:)
      integer, intent(in) :: n
      real(dp) :: e
      integer :: i, j

      e = sum(b(1:n)*c(1:n))
      do j = 1, n
         do i = 1, n
            e = e + c(i)*c(j)*a(i, j)
         end do
      end do
   end function model

   pure subroutine model_gradient(a, b, c, n, g)
      !! `dE/dc_k = b_k + sum_j c_j (A_kj + A_jk)`
      real(dp), intent(in) :: a(:, :), b(:), c(:)
      integer, intent(in) :: n
      real(dp), intent(out) :: g(:)
      integer :: k, j

      do k = 1, n
         g(k) = b(k)
         do j = 1, n
            g(k) = g(k) + c(j)*(a(k, j) + a(j, k))
         end do
      end do
   end subroutine model_gradient

end module mqc_ediis
