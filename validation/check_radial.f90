program check_radial
   !! Check the Treutler-Ahlrichs radial mesh against PySCF and against exact integrals
   !!
   !! Two independent checks. The first dumps nodes and weights for several
   !! elements and sizes so compare_radial.py can match them to PySCF's
   !! radi.treutler_ahlrichs. The second integrates functions whose spherical
   !! integrals are known in closed form, which catches a mesh that agrees with
   !! PySCF because both are wrong in the same way -- they would not be, but the
   !! check costs nothing and does not depend on PySCF being installed.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_dft_radial, only: treutler_ahlrichs_radial, radial_volume_weights, &
                             treutler_xi, bragg_radius
   use mqc_physical_constants, only: PI
   implicit none

   integer, parameter :: N_ELEMENTS = 5
   integer, parameter :: N_SIZES = 2
   integer, parameter :: ELEMENTS(N_ELEMENTS) = [1, 6, 8, 26, 79]
   integer, parameter :: SIZES(N_SIZES) = [25, 75]

   real(dp), allocatable :: r(:), dr(:), w(:)
   type(error_t) :: error
   integer :: unit, e, s, n_bad
   real(dp) :: integral, exact, worst_rel

   n_bad = 0
   worst_rel = 0.0_dp

   open (newunit=unit, file="radial_mesh.txt", status="replace", action="write")
   do e = 1, size(ELEMENTS)
      do s = 1, size(SIZES)
         call treutler_ahlrichs_radial(SIZES(s), ELEMENTS(e), r, dr, error)
         if (error%has_error()) then
            write (*, "(a,a)") "FAIL: ", error%get_message()
            n_bad = n_bad + 1
            cycle
         end if
         call dump(unit, ELEMENTS(e), SIZES(s), r, dr)
         deallocate (r, dr)
      end do
   end do
   close (unit)

   ! Two integrals with known values, testing different things.
   !
   ! exp(-r^2) is contained well inside the mesh, so it probes the Jacobian
   ! alone and should be exact to rounding.
   !
   ! exp(-r) is not contained: the M4 mapping at 200 points reaches only
   ! r = 17 Bohr, where exp(-r) is still 4e-8, so the error is the truncated
   ! tail rather than the quadrature. That is a property of the grid, not a
   ! defect -- PySCF's mesh gives the same 9.567e-8 -- and it is worth a test
   ! because it is the reason a diffuse density needs more radial points than
   ! its shape alone would suggest. So the assertion is on convergence: the
   ! error must fall as the mesh reaches further out.
   call treutler_ahlrichs_radial(200, 1, r, dr, error)
   allocate (w(size(r)))
   call radial_volume_weights(r, dr, w)
   integral = sum(w*exp(-r*r))
   call report("exp(-r^2), n=200  ", integral, PI**1.5_dp, worst_rel, n_bad, 1.0e-13_dp)
   deallocate (r, dr, w)

   call check_tail_convergence(n_bad)

   write (*, "(a,f6.3,a,f6.3)") "xi(H) = ", treutler_xi(1), "   xi(Au) = ", treutler_xi(79)
   write (*, "(a,f6.3,a,f6.3)") "bragg(H) = ", bragg_radius(1), "   bragg(U) = ", bragg_radius(92)
   write (*, "(a,es10.3)") "worst relative error: ", worst_rel

   if (n_bad > 0) then
      write (*, "(a)") "FAILED"
      stop 1
   end if
   write (*, "(a)") "radial OK -- now run compare_radial.py"

contains

   subroutine dump(unit, z, n, r, dr)
      integer, intent(in) :: unit, z, n
      real(dp), intent(in) :: r(:), dr(:)
      integer :: k
      do k = 1, size(r)
         write (unit, "(i0,1x,i0,2(1x,es24.16e3))") z, n, r(k), dr(k)
      end do
   end subroutine dump

   subroutine report(label, got, want, worst, n_bad, tol)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: got, want
      real(dp), intent(inout) :: worst
      integer, intent(inout) :: n_bad
      real(dp), intent(in) :: tol
      real(dp) :: rel

      rel = abs(got - want)/abs(want)
      worst = max(worst, rel)
      write (*, "(a,a,f20.14,a,f20.14,a,es10.3)") label, " = ", got, "  exact ", want, "  rel ", rel
      if (rel > tol) then
         write (*, "(a,es10.3)") "   FAIL: over tolerance ", tol
         n_bad = n_bad + 1
      end if
   end subroutine report

   subroutine check_tail_convergence(n_bad)
      !! exp(-r) must converge as the mesh reaches further out
      integer, intent(inout) :: n_bad

      integer, parameter :: MESHES(4) = [100, 200, 400, 800]
      real(dp), allocatable :: rr(:), ddr(:), ww(:)
      type(error_t) :: err
      real(dp) :: rel(size(MESHES))
      real(dp) :: exact
      integer :: k

      exact = 8.0_dp*PI
      do k = 1, size(MESHES)
         call treutler_ahlrichs_radial(MESHES(k), 1, rr, ddr, err)
         allocate (ww(size(rr)))
         call radial_volume_weights(rr, ddr, ww)
         rel(k) = abs(sum(ww*exp(-rr)) - exact)/exact
         write (*, "(a,i4,a,es10.3,a,f6.2)") "exp(-r)   n=", MESHES(k), &
            "  rel ", rel(k), "   mesh reaches r = ", maxval(rr)
         deallocate (rr, ddr, ww)
      end do

      do k = 2, size(MESHES)
         if (rel(k) >= rel(k - 1)) then
            write (*, "(a,i0)") "   FAIL: error did not fall going to n=", MESHES(k)
            n_bad = n_bad + 1
         end if
      end do
      if (rel(size(MESHES)) > 1.0e-8_dp) then
         write (*, "(a)") "   FAIL: exp(-r) not converged by n=800"
         n_bad = n_bad + 1
      end if
   end subroutine check_tail_convergence

end program check_radial
