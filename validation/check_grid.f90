program check_grid
   !! Integrate functions with known values on the assembled molecular grid
   !!
   !! This is the check that matters: the pieces are individually correct
   !! already, and what remains to go wrong is how they are combined -- a
   !! missing r^2, a 4*pi applied twice or not at all, the partition multiplied
   !! into the wrong points. Every one of those returns a smooth wrong number
   !! rather than failing, and every one of them is caught by integrating
   !! something whose answer is known.
   !!
   !! Two integrands, chosen for different reasons:
   !!
   !!   sum of Gaussians   smooth, so any error is in the assembly, not the mesh
   !!   sum of Slaters     cusped at each nucleus, which is what the radial mesh
   !!                      is built to resolve and where a product grid would
   !!                      struggle without the partition
   !!
   !! Both are sums over atoms, so each atomic sub-grid has to integrate parts
   !! of every other atom's function -- which is exactly what the partition is
   !! for and what a per-atom check would miss.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_dft_grid, only: dft_grid_t, build_dft_grid, grid_level_radial, grid_level_angular
   use mqc_dft_prune, only: PRUNE_NONE, PRUNE_NWCHEM, prune_scheme_name
   implicit none

   integer, parameter :: N_DIM = 3
   integer, parameter :: MAX_LEVEL = 5
   integer, parameter :: N_INTEGRANDS = 2
   real(dp), parameter :: PI = 3.14159265358979323846_dp
   real(dp), parameter :: GAUSS_ALPHA = 1.3_dp
   real(dp), parameter :: SLATER_ZETA = 1.7_dp

   real(dp) :: water(N_DIM, 3), auh(N_DIM, 2)
   integer :: n_bad

   water = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                    0.0_dp, -1.4308_dp, 1.1078_dp, &
                    0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3])
   auh = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                  0.0_dp, 0.0_dp, 2.8_dp], [N_DIM, 2])

   n_bad = 0

   call run_case("water", water, [8, 1, 1], PRUNE_NONE, n_bad)
   call run_case("auh  ", auh, [79, 1], PRUNE_NONE, n_bad)
   call run_case("water", water, [8, 1, 1], PRUNE_NWCHEM, n_bad)
   call run_case("auh  ", auh, [79, 1], PRUNE_NWCHEM, n_bad)

   if (n_bad > 0) then
      write (*, "(a)") "FAILED"
      stop 1
   end if
   write (*, "(a)") "grid integration OK -- now run compare_grid.py"

contains

   subroutine run_case(label, atom_coords, atomic_numbers, prune, n_bad)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: atom_coords(:, :)
      integer, intent(in) :: atomic_numbers(:)
      integer, intent(in) :: prune
      integer, intent(inout) :: n_bad

      type(dft_grid_t) :: grid
      type(error_t) :: error
      integer :: level, n_atoms, k
      real(dp) :: got_gauss, got_slater, exact
      real(dp) :: rel(MAX_LEVEL, N_INTEGRANDS)

      n_atoms = size(atomic_numbers)
      exact = real(n_atoms, dp)

      write (*, "(a)") ""
      write (*, "(a,a,a,i0,a,a)") "== ", label, " (", n_atoms, " atoms), pruning: ", &
         prune_scheme_name(prune)
      write (*, "(a)") "level   points      gaussians          slaters"

      do level = 1, MAX_LEVEL
         call build_dft_grid(atom_coords, atomic_numbers, grid, error, level=level, &
                             prune=prune)
         if (error%has_error()) then
            write (*, "(a,a)") "FAIL: ", error%get_message()
            n_bad = n_bad + 1
            return
         end if

         if (level == 3 .and. prune == PRUNE_NONE) call dump_grid(label, grid)
         if (level == 3 .and. prune == PRUNE_NWCHEM) call dump_grid("pruned_"//label, grid)

         got_gauss = integrate_gaussians(grid, atom_coords)
         got_slater = integrate_slaters(grid, atom_coords)

         write (*, "(i3,i10,2(3x,es15.6))") level, grid%n_points, &
            abs(got_gauss - exact)/exact, abs(got_slater - exact)/exact

         ! Tolerances are calibrated against what PySCF's identical grid
         ! delivers, not guessed: at level 3 it reaches 4.3e-9 on the Gaussians
         ! and 5.3e-10 on the Slaters. Level 3 is the production default and is
         ! not a converged grid, so demanding much more of it would be testing
         ! the wrong thing.
         rel(level, 1) = abs(got_gauss - exact)/exact
         rel(level, 2) = abs(got_slater - exact)/exact

         ! Level 3 must reach the same accuracy pruned as unpruned -- that is
         ! the entire claim pruning makes, and it holds for both molecules.
         if (level == 3) then
            call assert_below("level 3 gaussians", rel(level, 1), 1.0e-8_dp, n_bad)
            call assert_below("level 3 slaters  ", rel(level, 2), 1.0e-8_dp, n_bad)
         end if

         ! Refining only reaches machine precision unpruned. NWChem's inner two
         ! zones are fixed at 50 and 86 points whatever the level, so a pruned
         ! grid saturates instead of converging -- AuH's Slater integral even
         ! worsens from level 4 to 5, and PySCF's pruned grid does the identical
         ! thing to the same digits. Asserting convergence here would be
         ! asserting something the scheme does not promise.
         if (level == 5 .and. prune == PRUNE_NONE) then
            call assert_below("level 5 gaussians", rel(level, 1), 1.0e-10_dp, n_bad)
            call assert_below("level 5 slaters  ", rel(level, 2), 1.0e-10_dp, n_bad)
         end if
         if (level == 5 .and. prune /= PRUNE_NONE) then
            call assert_below("level 5 gaussians", rel(level, 1), 1.0e-8_dp, n_bad)
            call assert_below("level 5 slaters  ", rel(level, 2), 1.0e-8_dp, n_bad)
         end if

         call grid%destroy()
      end do

      ! A grid that is merely close at one level could be wrong in a way that
      ! cancels. Refining must improve it, every step, for both integrands --
      ! but only unpruned, for the reason given above.
      if (prune /= PRUNE_NONE) return
      do k = 2, MAX_LEVEL
         if (rel(k, 1) >= rel(k - 1, 1)) then
            write (*, "(a,i0)") "   FAIL: gaussians did not improve at level ", k
            n_bad = n_bad + 1
         end if
         if (rel(k, 2) >= rel(k - 1, 2)) then
            write (*, "(a,i0)") "   FAIL: slaters did not improve at level ", k
            n_bad = n_bad + 1
         end if
      end do
   end subroutine run_case

   subroutine assert_below(label, value, tol, n_bad)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: value, tol
      integer, intent(inout) :: n_bad

      if (value > tol) then
         write (*, "(a,a,a,es10.3,a,es10.3)") "   FAIL: ", label, " = ", value, " over ", tol
         n_bad = n_bad + 1
      end if
   end subroutine assert_below

   subroutine dump_grid(label, grid)
      !! Write the grid for the set comparison against PySCF
      character(len=*), intent(in) :: label
      type(dft_grid_t), intent(in) :: grid
      integer :: unit, k

      open (newunit=unit, file="grid_"//trim(adjustl(label))//".txt", &
            status="replace", action="write")
      do k = 1, grid%n_points
         write (unit, "(3(es24.16e3,1x),es24.16e3)") grid%coords(:, k), grid%weights(k)
      end do
      close (unit)
   end subroutine dump_grid

   function integrate_gaussians(grid, atom_coords) result(total)
      !! sum over atoms of (alpha/pi)^(3/2) exp(-alpha |r - R_A|^2), which
      !! integrates to the number of atoms
      type(dft_grid_t), intent(in) :: grid
      real(dp), intent(in) :: atom_coords(:, :)
      real(dp) :: total
      real(dp) :: f, d2, norm
      integer :: k, ia

      norm = (GAUSS_ALPHA/PI)**1.5_dp
      total = 0.0_dp
      do k = 1, grid%n_points
         f = 0.0_dp
         do ia = 1, size(atom_coords, 2)
            d2 = sum((grid%coords(:, k) - atom_coords(:, ia))**2)
            f = f + norm*exp(-GAUSS_ALPHA*d2)
         end do
         total = total + grid%weights(k)*f
      end do
   end function integrate_gaussians

   function integrate_slaters(grid, atom_coords) result(total)
      !! sum over atoms of (zeta^3/pi) exp(-2 zeta |r - R_A|), integrating to
      !! the number of atoms; cusped at every nucleus
      type(dft_grid_t), intent(in) :: grid
      real(dp), intent(in) :: atom_coords(:, :)
      real(dp) :: total
      real(dp) :: f, d, norm
      integer :: k, ia

      norm = SLATER_ZETA**3/PI
      total = 0.0_dp
      do k = 1, grid%n_points
         f = 0.0_dp
         do ia = 1, size(atom_coords, 2)
            d = norm2(grid%coords(:, k) - atom_coords(:, ia))
            f = f + norm*exp(-2.0_dp*SLATER_ZETA*d)
         end do
         total = total + grid%weights(k)*f
      end do
   end function integrate_slaters

end program check_grid
