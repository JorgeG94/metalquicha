module test_mqc_libcint_ao
   !! Pins basis-function evaluation against the integrals it has to agree with.
   !!
   !! libcint does not evaluate basis functions, so these values come from our own
   !! code and there is nothing in the integral library to check them against
   !! directly. What there is instead is an identity: the overlap matrix is the
   !! integral of a product of two basis functions, so integrating that product
   !! over the DFT grid must reproduce the analytic overlap. That single check
   !! covers the AO values, the normalisation, the Cartesian ordering, the
   !! spherical transform, *and* the grid weights and Becke partition together --
   !! and it needs no external reference, which is what makes it a unit test rather
   !! than a validation script.
   !!
   !! It cannot be exact. A quadrature is not an integral, and the residual is the
   !! grid's own error: at level 3 on water/cc-pVDZ it is around 1e-7, and a much
   !! smaller number would mean the test had stopped measuring anything. So the
   !! assertion is two-sided in spirit: small enough to prove the conventions,
   !! large enough to still be a quadrature.
   !!
   !! The elementwise comparison against PySCF's `eval_ao` lives in
   !! `validation/check_ao.f90` and its Python companion, because it needs PySCF.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_ao, only: eval_ao_block, eval_rho, max_ao_l, &
                             shell_extents, block_significant_aos, &
                             AO_HESS_COMP, AO_DERIV3_COMP
   use mqc_dft_grid, only: dft_grid_t, build_dft_grid
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   implicit none
   private
   public :: collect_mqc_libcint_ao_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_libcint_ao_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("numerical_overlap_matches_analytic", test_overlap_spherical), &
                  new_unittest("numerical_overlap_matches_analytic_cartesian", test_overlap_cart), &
                  new_unittest("grid_refinement_reduces_the_error", test_convergence), &
                  new_unittest("value_on_a_nucleus_is_finite", test_on_nucleus), &
                  new_unittest("high_angular_momentum_is_refused", test_l_limit), &
                  new_unittest("integrated_density_is_the_electron_count", test_rho), &
                  new_unittest("ao_gradients_match_finite_differences", test_ao_grad), &
                  new_unittest("shell_extent_really_bounds_the_shell", test_extent_bound), &
                  new_unittest("screening_off_keeps_every_function", test_screen_off), &
                  new_unittest("screened_block_matches_the_full_one", test_screen_matches), &
                  new_unittest("compressed_atom_ranges_are_contiguous", test_atom_ranges), &
                  new_unittest("ao_third_derivatives_match_differences", test_ao_deriv3) &
                  ]
   end subroutine collect_mqc_libcint_ao_tests

   subroutine water(mol, err, basis)
      type(libcint_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      character(len=*), intent(in) :: basis
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.10077199490609_dp*ANG, &
                   0.0_dp, 0.77250895271063_dp*ANG, -0.46780199741728_dp*ANG, &
                   0.0_dp, -0.77250895280218_dp*ANG, -0.46780199748881_dp*ANG], [3, 3])
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
   end subroutine water

   subroutine overlap_error(basis, level, worst, err)
      !! max |S_numerical - S_analytic| over the whole matrix
      character(len=*), intent(in) :: basis
      integer, intent(in) :: level
      real(dp), intent(out) :: worst
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
      type(dft_grid_t) :: grid
      real(dp), allocatable :: ao(:, :), s_ana(:, :), s_num(:, :)
      integer :: mu, nu, ig

      worst = huge(1.0_dp)
      call water(mol, err, basis)
      if (err%has_error()) return

      call build_dft_grid(mol%coords, [8, 1, 1], grid, err, level=level)
      if (err%has_error()) return

      call eval_ao_block(mol, grid%coords, ao, err)
      if (err%has_error()) return

      call mol%overlap(s_ana)
      allocate (s_num(mol%nao, mol%nao))
      s_num = 0.0_dp
      do nu = 1, mol%nao
         do mu = 1, mol%nao
            do ig = 1, grid%n_points
               s_num(mu, nu) = s_num(mu, nu) + grid%weights(ig)*ao(ig, mu)*ao(ig, nu)
            end do
         end do
      end do

      worst = maxval(abs(s_num - s_ana))
      call grid%destroy()
      call mol%destroy()
   end subroutine overlap_error

   subroutine test_atom_ranges(error)
      !! The per-atom compressed ranges partition the kept set, in order
      !!
      !! The gradient sums over one atom's range and lands the result on that
      !! atom, so a range that is off by one, overlaps its neighbour or leaves a
      !! gap attributes force to the wrong nucleus. That is a wrong gradient
      !! rather than a slow one, and it does not announce itself -- the energy
      !! is untouched and the forces still look like forces.
      !!
      !! Checked by reconstruction rather than by eye: walking the atoms in
      !! order and concatenating their ranges has to reproduce `ao_list`
      !! exactly, which fails on an overlap, a gap, or a wrong start.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: radius(:), pts(:, :)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_list(:), ao_offset(:), a_off(:), a_cnt(:)
      integer :: n_sig, iatom, i, seen, k
      logical :: ok

      call water(mol, err, "6-31g")
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return

      call shell_extents(mol, 1.0e-10_dp, radius)
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(mol%nao))
      allocate (a_off(mol%natm), a_cnt(mol%natm))

      ! Offset from every nucleus so the screen actually drops something --
      ! a block that keeps everything would not exercise the compression.
      allocate (pts(3, 3))
      do k = 1, 3
         pts(:, k) = mol%coords(:, 1) + [2.0_dp + 0.1_dp*k, 0.4_dp, 0.3_dp]
      end do

      call block_significant_aos(mol, pts, radius, shell_mask, ao_list, &
                                 ao_offset, n_sig, atom_offsets=a_off, &
                                 atom_counts=a_cnt)
      call check(error, n_sig > 0, "the block should keep something")
      if (allocated(error)) return

      call check(error, sum(a_cnt) == n_sig, &
                 "the per-atom counts should add up to the kept count")
      if (allocated(error)) return

      ! Walk the atoms in order; their ranges must lay out ao_list end to end.
      seen = 0
      ok = .true.
      do iatom = 1, mol%natm
         if (a_cnt(iatom) == 0) cycle
         if (a_off(iatom) /= seen) ok = .false.
         do i = 1, a_cnt(iatom)
            seen = seen + 1
            if (seen > n_sig) then
               ok = .false.
               exit
            end if
         end do
      end do
      call check(error, ok .and. seen == n_sig, &
                 "the per-atom ranges do not tile the compressed set in atom "// &
                 "order, so an atom's functions are not contiguous")
      if (allocated(error)) return

      ! And every kept function must belong to exactly one atom's range.
      do iatom = 1, mol%natm
         do i = 1, a_cnt(iatom)
            call check(error, ao_list(a_off(iatom) + i) >= 1 .and. &
                       ao_list(a_off(iatom) + i) <= mol%nao, &
                       "an atom range points outside the basis")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_atom_ranges

   subroutine test_extent_bound(error)
      !! Past every shell's radius, every basis function is below threshold
      !!
      !! The point of the test rather than of the radius: `shell_extents`
      !! returns something a truncation argument rests on, so it has to be
      !! checked against the functions themselves and not against the formula
      !! that produced it. An envelope that bounds one primitive instead of the
      !! contracted sum passes every energy comparison and is still not a bound.
      !!
      !! Tested against the largest radius in the molecule rather than shell by
      !! shell, since the shell-to-atom map lives behind libcint and this suite
      !! does not link it. Weaker, and still fails on a broken envelope: the
      !! bound has to hold for whichever shell reaches furthest.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: radius(:), ao(:, :), pts(:, :)
      real(dp), parameter :: TOL = 1.0e-10_dp
      integer :: iatom, k
      real(dp) :: far, worst, dir(3)

      call water(mol, err, "6-31g")
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return

      call shell_extents(mol, TOL, radius)
      call check(error, size(radius) == mol%nbas, "one radius per shell")
      if (allocated(error)) return
      call check(error, all(radius > 0.0_dp), "every radius should be positive")
      if (allocated(error)) return

      ! Far enough that every atom is past its furthest-reaching shell.
      far = maxval(radius)*1.001_dp
      do iatom = 1, mol%natm
         far = max(far, maxval(radius) + norm2(mol%coords(:, iatom)))
      end do

      allocate (pts(3, 6))
      do k = 1, 6
         dir = 0.0_dp
         dir(1 + mod(k - 1, 3)) = merge(1.0_dp, -1.0_dp, k <= 3)
         pts(:, k) = dir*far
      end do

      call eval_ao_block(mol, pts, ao, err)
      call check(error,.not. err%has_error(), "evaluating past every radius")
      if (allocated(error)) return

      worst = maxval(abs(ao))
      call check(error, worst <= TOL, &
                 "a basis function is above the threshold beyond every shell "// &
                 "radius, so the extent is not a bound")
   end subroutine test_extent_bound

   subroutine test_screen_off(error)
      !! A non-positive threshold keeps the whole basis
      !!
      !! The escape hatch a deck uses to measure what the screening is worth,
      !! so it has to genuinely disable it rather than merely widen it.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: radius(:), pts(:, :)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_list(:), ao_offset(:)
      integer :: n_sig

      call water(mol, err, "6-31g")
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return

      call shell_extents(mol, -1.0_dp, radius)
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(mol%nao))
      allocate (pts(3, 1))
      ! Far enough away that any real screen would drop everything.
      pts(:, 1) = mol%coords(:, 1) + [50.0_dp, 0.0_dp, 0.0_dp]
      call block_significant_aos(mol, pts, radius, shell_mask, ao_list, &
                                 ao_offset, n_sig)
      call check(error, n_sig == mol%nao, &
                 "screening off should keep every basis function")
   end subroutine test_screen_off

   subroutine test_screen_matches(error)
      !! The compressed block holds the same numbers as the full one
      !!
      !! Checks the plumbing rather than the physics: `ao_list` has to map the
      !! compressed columns back to the ones they came from. A transposed or
      !! off-by-one map still produces a plausible density and a wrong energy.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: radius(:), pts(:, :), full(:, :), part(:, :)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_list(:), ao_offset(:)
      integer :: n_sig, i, k
      real(dp) :: worst

      call water(mol, err, "6-31g")
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return

      call shell_extents(mol, 1.0e-10_dp, radius)
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(mol%nao))
      allocate (pts(3, 4))
      do k = 1, 4
         pts(:, k) = mol%coords(:, 1) + [0.3_dp*k, 0.1_dp, -0.2_dp]
      end do

      call block_significant_aos(mol, pts, radius, shell_mask, ao_list, &
                                 ao_offset, n_sig)
      call check(error, n_sig > 0 .and. n_sig <= mol%nao, "a sane kept count")
      if (allocated(error)) return

      call eval_ao_block(mol, pts, full, err)
      call check(error,.not. err%has_error(), "the full evaluation")
      if (allocated(error)) return
      call eval_ao_block(mol, pts, part, err, shell_mask=shell_mask, &
                         ao_offset=ao_offset, n_ao_out=n_sig)
      call check(error,.not. err%has_error(), "the compressed evaluation")
      if (allocated(error)) return

      call check(error, size(part, 2) == n_sig, "the compressed width")
      if (allocated(error)) return
      worst = 0.0_dp
      do i = 1, n_sig
         worst = max(worst, maxval(abs(part(:, i) - full(:, ao_list(i)))))
      end do
      call check(error, worst < 1.0e-14_dp, &
                 "a compressed column does not match the function ao_list "// &
                 "says it came from")
   end subroutine test_screen_matches

   subroutine test_overlap_spherical(error)
      !! Spherical basis with d functions: cc-pVDZ
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: worst

      call overlap_error("cc-pvdz", 3, worst, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      ! Small enough that the conventions must be right: a transposed transform or
      ! a doubled normalisation is a relative error of order one, not 1e-6.
      call check(error, worst < 1.0e-5_dp, "numerical overlap disagrees with analytic")
      if (allocated(error)) return
      ! And large enough that this is still a quadrature. If it ever comes out at
      ! 1e-14 the grid has been replaced by something exact and the test has
      ! stopped covering the grid.
      call check(error, worst > 1.0e-12_dp, &
                 "numerical overlap is suspiciously exact -- is the grid still a grid?")
   end subroutine test_overlap_spherical

   subroutine test_overlap_cart(error)
      !! Cartesian basis: 6-31G**, which is Cartesian from lithium up
      !!
      !! Separate case because the Cartesian path skips the spherical transform
      !! entirely, so it is the one that would still pass if that transform were
      !! deleted.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: worst

      call overlap_error("6-31g_st__st_", 3, worst, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, worst < 1.0e-5_dp, &
                 "Cartesian numerical overlap disagrees with analytic")
   end subroutine test_overlap_cart

   subroutine test_convergence(error)
      !! A finer grid must integrate better
      !!
      !! The sharpest statement available without an external number: whatever the
      !! residual is, it has to be *quadrature* error, and quadrature error falls
      !! when the quadrature improves. A constant offset -- a wrong normalisation,
      !! say -- would not move.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: coarse, fine

      call overlap_error("cc-pvdz", 1, coarse, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call overlap_error("cc-pvdz", 5, fine, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      call check(error, fine < coarse, &
                 "refining the grid did not reduce the overlap error, so the residual "// &
                 "is not quadrature error")
   end subroutine test_convergence

   subroutine test_on_nucleus(error)
      !! A point exactly on a nucleus evaluates, and to something positive
      !!
      !! Every Cartesian component there is x**0, which must be one rather than
      !! whatever a general power routine does with 0**0. Grids really do put
      !! points on nuclei.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: ao(:, :)
      real(dp) :: pts(3, 1)

      call water(mol, err, "cc-pvdz")
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      pts(:, 1) = mol%coords(:, 1)
      call eval_ao_block(mol, pts, ao, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      ! The oxygen 1s is the largest function anywhere and it peaks here.
      call check(error, ao(1, 1) > 1.0_dp, "the 1s value on its own nucleus should be large")
      if (allocated(error)) return
      call check(error, all(ao(1, :) == ao(1, :)), "a NaN appeared on the nucleus")
      call mol%destroy()
   end subroutine test_on_nucleus

   subroutine test_l_limit(error)
      !! The transform table's range is enforced rather than exceeded
      !!
      !! Beyond the tabulated angular momenta the functions would be silently
      !! wrong, which is the one outcome worse than an error.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(error_t) :: err

      call water(mol, err, "cc-pvtz")
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      ! cc-pVTZ has f functions on oxygen, which is l = 3 and inside the table.
      call check(error, max_ao_l(mol) >= 3, "cc-pVTZ should reach l = 3")
      call mol%destroy()
   end subroutine test_l_limit

   subroutine test_rho(error)
      !! sum_g w_g rho(r_g) must be the number of electrons
      !!
      !! An exact identity, and the one milestone-1 check that needs neither an
      !! external reference nor a functional. It covers the density assembly and
      !! the grid weights at once, and it is the check that would catch a factor
      !! of two in the restricted density -- an error that otherwise survives all
      !! the way to an exchange-correlation energy that is merely wrong.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(dft_grid_t) :: grid
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: ao(:, :), rho(:)
      real(dp) :: n_elec

      call water(mol, err, "cc-pvdz")
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      call run_libcint_rhf(mol, 10, 200, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err, &
                           in_core=.true.)
      call check(error,.not. err%has_error() .and. scf%converged, "water SCF must converge")
      if (allocated(error)) return

      call build_dft_grid(mol%coords, [8, 1, 1], grid, err, level=4)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      call eval_ao_block(mol, grid%coords, ao, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      call eval_rho(ao, scf%density, rho)
      n_elec = sum(grid%weights*rho)

      ! Ten electrons in water. The tolerance is the grid's, not the identity's.
      call check(error, abs(n_elec - 10.0_dp) < 1.0e-5_dp, &
                 "the integrated density is not the electron count")
      if (allocated(error)) return
      ! Positive everywhere it matters: a sign or transpose error shows here.
      call check(error, minval(rho) > -1.0e-10_dp, "the density went negative")

      call grid%destroy()
      call mol%destroy()
   end subroutine test_rho

   subroutine test_ao_grad(error)
      !! d chi / d r against a central difference of chi
      !!
      !! Worth its own test rather than being folded into a GGA energy: a wrong
      !! gradient and a wrong potential expression both move a GGA total, and
      !! separating them afterwards means reasoning about which. Finite differences
      !! settle the gradient on its own, with no functional in sight.
      !!
      !! Central differences at h = 1e-4 are good to about h^2 relative, so the
      !! tolerance is 1e-6 and not tighter -- a smaller bound would be testing the
      !! difference formula rather than the derivative.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: ao(:, :), grad(:, :, :), plus(:, :), minus(:, :)
      real(dp) :: pts(3, 3), shifted(3, 3), fd, worst
      real(dp), parameter :: H = 1.0e-4_dp
      integer :: id, ig, mu

      call water(mol, err, "cc-pvdz")
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      ! Off-axis and off-nucleus, so no component is zero by symmetry and the
      ! difference is not taken across a cusp.
      pts = reshape([0.31_dp, 0.17_dp, -0.23_dp, &
                     -0.44_dp, 0.62_dp, 0.11_dp, &
                     1.30_dp, -0.85_dp, 0.47_dp], [3, 3])

      call eval_ao_block(mol, pts, ao, err, grad=grad)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      worst = 0.0_dp
      do id = 1, 3
         shifted = pts
         shifted(id, :) = pts(id, :) + H
         call eval_ao_block(mol, shifted, plus, err)
         shifted(id, :) = pts(id, :) - H
         call eval_ao_block(mol, shifted, minus, err)
         do mu = 1, mol%nao
            do ig = 1, 3
               fd = (plus(ig, mu) - minus(ig, mu))/(2.0_dp*H)
               worst = max(worst, abs(fd - grad(ig, mu, id)))
            end do
         end do
      end do

      call check(error, worst < 1.0e-6_dp, "AO gradients disagree with finite differences")
      call mol%destroy()
   end subroutine test_ao_grad

   subroutine test_ao_deriv3(error)
      !! Third derivatives against central differences of the analytic second
      !!
      !! Differenced against the *second* derivatives rather than the values,
      !! for the same reason the partition cutoffs are: `hess` is already
      !! checked against finite differences of `grad`, so a disagreement here
      !! localises to the third derivative instead of being shared between two
      !! candidates. Every one of the ten packed components is reached, and on
      !! cc-pVDZ that means d functions -- where a transposed exponent in the
      !! angular part shows up and an s or p basis would not notice.
      !!
      !! The packing is the other thing being tested. `deriv3` claims
      !! `GTOval_sph_deriv3` order, and the mapping from a component index to
      !! its three Cartesian axes has to agree with the mapping used to
      !! differentiate -- a permutation there is silent on the diagonal
      !! components and wrong on the mixed ones.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: ao(:, :), hess(:, :, :), d3(:, :, :)
      real(dp), allocatable :: hplus(:, :, :), hminus(:, :, :), tmp(:, :)
      real(dp) :: pts(3, 3), shifted(3, 3), fd, worst
      real(dp), parameter :: H = 1.0e-4_dp
      ! Which packed second-derivative component each third-derivative
      ! component reduces to once one axis is differenced numerically, and
      ! which axis that is.
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
      integer, parameter :: DI(AO_DERIV3_COMP) = [1, 1, 1, 1, 1, 1, 2, 2, 2, 3]
      integer, parameter :: DJ(AO_DERIV3_COMP) = [1, 1, 1, 2, 2, 3, 2, 2, 3, 3]
      integer, parameter :: DK(AO_DERIV3_COMP) = [1, 2, 3, 2, 3, 3, 2, 3, 3, 3]
      integer :: ih, ig, mu, ax

      call water(mol, err, "cc-pvdz")
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      pts = reshape([0.31_dp, 0.17_dp, -0.23_dp, &
                     -0.44_dp, 0.62_dp, 0.11_dp, &
                     1.30_dp, -0.85_dp, 0.47_dp], [3, 3])

      call eval_ao_block(mol, pts, ao, err, hess=hess, deriv3=d3)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, size(d3, 3) == AO_DERIV3_COMP, "deriv3 has the wrong component count")
      if (allocated(error)) return

      worst = 0.0_dp
      do ih = 1, AO_DERIV3_COMP
         ! Difference along the first of the three axes and compare against the
         ! second derivative in the remaining two.
         ax = DI(ih)
         shifted = pts
         shifted(ax, :) = pts(ax, :) + H
         call eval_ao_block(mol, shifted, tmp, err, hess=hplus)
         shifted(ax, :) = pts(ax, :) - H
         call eval_ao_block(mol, shifted, tmp, err, hess=hminus)
         do mu = 1, mol%nao
            do ig = 1, 3
               fd = (hplus(ig, mu, PAIR(DJ(ih), DK(ih))) &
                     - hminus(ig, mu, PAIR(DJ(ih), DK(ih))))/(2.0_dp*H)
               worst = max(worst, abs(fd - d3(ig, mu, ih)))
            end do
         end do
      end do

      call check(error, worst < 1.0e-4_dp, &
                 "AO third derivatives disagree with differences of the second")
      call mol%destroy()
   end subroutine test_ao_deriv3

end module test_mqc_libcint_ao

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_libcint_ao, only: collect_mqc_libcint_ao_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_libcint_ao", collect_mqc_libcint_ao_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
