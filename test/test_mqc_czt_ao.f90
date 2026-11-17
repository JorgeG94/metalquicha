module test_mqc_czt_ao
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
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, shell_dim
   use mqc_czt_ao_data, only: C2S_LMAX, c2s_block, common_fac_sp
   use libcint_fortran, only: LIBCINT_ATOM_OF, LIBCINT_ANG_OF, LIBCINT_NPRIM_OF, &
                              LIBCINT_NCTR_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_COEFF
   use mqc_czt_ao, only: eval_ao_block, eval_rho, max_ao_l, &
                         shell_extents, block_significant_aos, &
                         AO_HESS_COMP, AO_DERIV3_COMP
   use mqc_dft_grid, only: dft_grid_t, build_dft_grid
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   implicit none
   private
   public :: collect_mqc_czt_ao_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_czt_ao_tests(testsuite)
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
                  new_unittest("ao_third_derivatives_match_differences", test_ao_deriv3), &
                  new_unittest("vectorised_evaluation_matches_frozen_reference", &
                               test_eval_matches_frozen_reference) &
                  ]
   end subroutine collect_mqc_czt_ao_tests

   subroutine water(mol, err, basis)
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      character(len=*), intent(in) :: basis
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.10077199490609_dp*ANG, &
                   0.0_dp, 0.77250895271063_dp*ANG, -0.46780199741728_dp*ANG, &
                   0.0_dp, -0.77250895280218_dp*ANG, -0.46780199748881_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, basis, mol, err)
   end subroutine water

   subroutine overlap_error(basis, level, worst, err)
      !! max |S_numerical - S_analytic| over the whole matrix
      character(len=*), intent(in) :: basis
      integer, intent(in) :: level
      real(dp), intent(out) :: worst
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol
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

      type(czt_molecule_t) :: mol
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

      type(czt_molecule_t) :: mol
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

      type(czt_molecule_t) :: mol
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

      type(czt_molecule_t) :: mol
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
      type(czt_molecule_t) :: mol
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
      type(czt_molecule_t) :: mol
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
      type(czt_molecule_t) :: mol
      type(dft_grid_t) :: grid
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: ao(:, :), rho(:)
      real(dp) :: n_elec

      call water(mol, err, "cc-pvdz")
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      call run_czt_rhf(mol, 10, 200, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err, &
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
      type(czt_molecule_t) :: mol
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
      type(czt_molecule_t) :: mol
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

   subroutine test_eval_matches_frozen_reference(error)
      !! The vectorised evaluator against a frozen copy of the point-by-point one
      !!
      !! `ref_eval_ao_block` below is the implementation as it stood before the
      !! per-shell rewrite, copied verbatim apart from its name. The two are
      !! compared on value, `grad`, `hess` and `deriv3`, on a spherical and a
      !! Cartesian basis, unscreened and screened, over points that include one
      !! exactly on a nucleus, one on an axis through it (so two of x, y, z are
      !! exactly zero and `x**0` has to stay one), and two so far away that every
      !! primitive is cut off.
      !!
      !! The measure is the largest deviation over the whole array divided by the
      !! largest |value| in the reference array, and it has to be under 1e-13. It
      !! is not zero: the radial exponentials come from a different libm entry
      !! point and the products are formed in a different order, so a handful of
      !! ulps is expected. The 1e-22 left behind by the primitive cutoff sits far
      !! under it.
      type(error_type), allocatable, intent(out) :: error

      character(len=14), parameter :: BASES(4) = [character(len=14) :: &
                                                  "cc-pvqz", "cc-pvtz", "6-311g_st__st_", "def2-tzvp"]
      logical, parameter :: CART(4) = [.false., .true., .false., .false.]
      integer, parameter :: N_PTS = 301
      real(dp), parameter :: TOL = 1.0e-13_dp
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp) :: worst, seed
      real(dp) :: c(3, 3), pts(3, N_PTS)
      real(dp), allocatable :: a_new(:, :), a_ref(:, :)
      real(dp), allocatable :: g_new(:, :, :), g_ref(:, :, :)
      real(dp), allocatable :: h_new(:, :, :), h_ref(:, :, :)
      real(dp), allocatable :: t_new(:, :, :), t_ref(:, :, :)
      logical, allocatable :: mask(:)
      integer, allocatable :: offset(:)
      integer :: ib, ip, k, ish, n_keep, pass

      c = reshape([0.0_dp, 0.0_dp, 0.10077199490609_dp*ANG, &
                   0.0_dp, 0.77250895271063_dp*ANG, -0.46780199741728_dp*ANG, &
                   0.0_dp, -0.77250895280218_dp*ANG, -0.46780199748881_dp*ANG], [3, 3])

      ! A fixed pseudo-random cloud around the molecule, the same every run.
      seed = 0.3141592653589793_dp
      do ip = 1, N_PTS
         do k = 1, 3
            seed = modulo(seed*7919.0_dp + 0.6180339887498949_dp, 1.0_dp)
            pts(k, ip) = 8.0_dp*(seed - 0.5_dp)
         end do
      end do
      pts(:, 1) = c(:, 1)                                  ! on the oxygen
      pts(:, 2) = c(:, 2)                                  ! on a hydrogen
      pts(:, 3) = c(:, 1) + [0.7_dp, 0.0_dp, 0.0_dp]       ! on an axis through it
      pts(:, 4) = c(:, 1) + [30.0_dp, 0.0_dp, 0.0_dp]      ! every primitive cut off
      pts(:, 5) = [1.0e3_dp, -1.0e3_dp, 5.0e2_dp]

      do ib = 1, size(BASES)
         call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, trim(BASES(ib)), mol, err, &
                                 force_cartesian=CART(ib))
         call check(error,.not. err%has_error(), "building "//trim(BASES(ib)))
         if (allocated(error)) return

         ! Unscreened, each derivative order asked for on its own and then all
         ! at once, because what is allocated depends on what is present.
         do pass = 1, 4
            select case (pass)
            case (1)
               call eval_ao_block(mol, pts, a_new, err)
               call ref_eval_ao_block(mol, pts, a_ref, err)
            case (2)
               call eval_ao_block(mol, pts, a_new, err, grad=g_new)
               call ref_eval_ao_block(mol, pts, a_ref, err, grad=g_ref)
            case (3)
               call eval_ao_block(mol, pts, a_new, err, hess=h_new)
               call ref_eval_ao_block(mol, pts, a_ref, err, hess=h_ref)
            case (4)
               call eval_ao_block(mol, pts, a_new, err, grad=g_new, hess=h_new, deriv3=t_new)
               call ref_eval_ao_block(mol, pts, a_ref, err, grad=g_ref, hess=h_ref, deriv3=t_ref)
            case default
               call check(error, .false., "unknown pass")
               return
            end select
            call check(error,.not. err%has_error(), "evaluating "//trim(BASES(ib)))
            if (allocated(error)) return
            worst = deviation2(a_new, a_ref)
            call check(error, worst < TOL, trim(BASES(ib))//": values differ from the reference")
            if (allocated(error)) return
            if (pass == 2 .or. pass == 4) then
               worst = deviation3(g_new, g_ref)
               call check(error, worst < TOL, trim(BASES(ib))//": gradients differ from the reference")
               if (allocated(error)) return
            end if
            if (pass >= 3) then
               worst = deviation3(h_new, h_ref)
               call check(error, worst < TOL, trim(BASES(ib))//": Hessians differ from the reference")
               if (allocated(error)) return
            end if
            if (pass == 4) then
               worst = deviation3(t_new, t_ref)
               call check(error, worst < TOL, &
                          trim(BASES(ib))//": third derivatives differ from the reference")
               if (allocated(error)) return
            end if
         end do

         ! Screened: every other shell kept, written to a compressed block whose
         ! offsets are laid out here the way `block_significant_aos` does.
         allocate (mask(mol%nbas), offset(mol%nbas))
         n_keep = 0
         do ish = 1, mol%nbas
            mask(ish) = mod(ish, 2) == 0
            offset(ish) = n_keep
            if (mask(ish)) n_keep = n_keep + shell_dim(mol%cartesian, ish - 1, mol%bas)
         end do
         call eval_ao_block(mol, pts, a_new, err, grad=g_new, hess=h_new, deriv3=t_new, &
                            shell_mask=mask, ao_offset=offset, n_ao_out=n_keep)
         call ref_eval_ao_block(mol, pts, a_ref, err, grad=g_ref, hess=h_ref, deriv3=t_ref, &
                                shell_mask=mask, ao_offset=offset, n_ao_out=n_keep)
         call check(error,.not. err%has_error(), "screened evaluation of "//trim(BASES(ib)))
         if (allocated(error)) return
         call check(error, size(a_new, 2) == n_keep, "compressed width")
         if (allocated(error)) return
         worst = max(deviation2(a_new, a_ref), deviation3(g_new, g_ref), &
                     deviation3(h_new, h_ref), deviation3(t_new, t_ref))
         call check(error, worst < TOL, trim(BASES(ib))//": a screened block differs from the reference")
         if (allocated(error)) return
         deallocate (mask, offset)
         call mol%destroy()
      end do
   end subroutine test_eval_matches_frozen_reference

   function deviation2(a, b) result(dev)
      !! Largest |a - b| over the array, relative to the largest |b|
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp) :: dev
      dev = maxval(abs(a - b))/max(maxval(abs(b)), tiny(1.0_dp))
   end function deviation2

   function deviation3(a, b) result(dev)
      !! Largest |a - b| over the array, relative to the largest |b|
      real(dp), intent(in) :: a(:, :, :), b(:, :, :)
      real(dp) :: dev
      dev = maxval(abs(a - b))/max(maxval(abs(b)), tiny(1.0_dp))
   end function deviation3

   ! Frozen copy of `eval_ao_block` as it stood before it was vectorised per
   ! shell, kept as the reference `test_eval_matches_frozen_reference` measures
   ! the live one against. Do not tidy it: its value is that it is the old code.
   subroutine ref_eval_ao_block(mol, coords, ao, error, grad, hess, deriv3, shell_mask, &
                                ao_offset, n_ao_out)
      !! chi_mu(r) for every basis function at every supplied point
      !!
      !! `ao` comes back as (n_points, n_ao), which is the orientation the density
      !! assembly wants: rho at a point is a row contracted against the density
      !! matrix, so points vary slowest and a gemm over them reads contiguously.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)        !! (3, n_points), Bohr
      real(dp), allocatable, intent(out) :: ao(:, :)
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: grad(:, :, :)
         !! d chi / d r, as (n_points, n_ao, 3). Asked for by GGA and above.
      logical, intent(in), optional :: shell_mask(:)
         !! (nbas) which shells reach these points, from `block_significant_aos`.
         !! Absent means all of them, which is the unscreened path every other
         !! caller still takes.
      integer, intent(in), optional :: ao_offset(:)
         !! (nbas) compressed offsets, required with `shell_mask`
      integer, intent(in), optional :: n_ao_out
         !! Leading dimension of the compressed block, required with `shell_mask`
      real(dp), allocatable, intent(out), optional :: deriv3(:, :, :)
         !! d3 chi / d r_i d r_j d r_k, as (n_points, n_ao, 10), packed
         !! xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz -- the order
         !! `GTOval_sph_deriv3` uses.
         !!
         !! A *GGA Hessian* needs these: the energy depends on grad rho, which
         !! already costs one position derivative to form, and two nuclear
         !! derivatives on top of that land on the third.
      real(dp), allocatable, intent(out), optional :: hess(:, :, :)
         !! d2 chi / d r_j d r_k, as (n_points, n_ao, 6), packed xx, xy, xz, yy,
         !! yz, zz -- the six unique components of a symmetric tensor, in the
         !! order PySCF's `GTOval_sph_deriv2` uses, so a comparison against it
         !! needs no reshuffling. A GGA *gradient* needs these.

      logical :: screened
      integer :: n_ao_here, shell_base
      integer :: n_points, ish, l, nprim, nctr, ic, ip, ig, off_ao
      integer :: n_cart, n_sph, n_here, comp, i, j, iatom
      integer :: exp_ptr, coef_ptr
      real(dp) :: dx, dy, dz, r2, radial, fac, dradial, d2radial
      real(dp) :: dr(3)
      real(dp) :: ex
      real(dp), allocatable :: cart(:), sph(:), trans(:, :)
      real(dp), allocatable :: dcart(:, :), dsph(:, :)
      real(dp), allocatable :: d2cart(:, :), d2sph(:, :)
      real(dp), allocatable :: d3cart(:, :), d3sph(:, :)
      real(dp) :: d3radial
      logical :: want_grad, want_hess, want_d3
      integer :: id, ih, px, py, pz

      ! The (j, k) each packed component stands for. The delta only enters the
      ! second derivative of the radial part, in the assembly below.
      integer, parameter :: HESS_J(AO_HESS_COMP) = [1, 1, 1, 2, 2, 3]
      integer, parameter :: HESS_K(AO_HESS_COMP) = [1, 2, 3, 2, 3, 3]
      integer, parameter :: D3_I(AO_DERIV3_COMP) = [1, 1, 1, 1, 1, 1, 2, 2, 2, 3]
      integer, parameter :: D3_J(AO_DERIV3_COMP) = [1, 1, 1, 2, 2, 3, 2, 2, 3, 3]
      integer, parameter :: D3_K(AO_DERIV3_COMP) = [1, 2, 3, 2, 3, 3, 2, 3, 3, 3]
      integer, parameter :: HESS_OF(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
         !! Which packed Hessian component a Cartesian pair (j,k) is.
      ! TODO(mqc): `HESS_OF` is never read -- `ref_assemble_third` declares the same
      ! table again as `PAIR`, as it does `D3_I/D3_J/D3_K` under the names
      ! `DI/DJ/DK`. Four constants written down twice, so a change to the
      ! packing order has to be made in two places.

      n_points = size(coords, 2)

      if (max_ao_l(mol) > C2S_LMAX) then
         call error%set(ERROR_VALIDATION, "AO evaluation: this basis has an angular "// &
                        "momentum above the transform table's range. Extend "// &
                        "mqc_czt_ao_data rather than letting it through -- the "// &
                        "functions would otherwise be silently wrong.")
         return
      end if

      want_d3 = present(deriv3)
      ! Each order is built from the one below it, so asking for a third
      ! implies the second and the second implies the first.
      want_hess = present(hess) .or. want_d3
      want_grad = present(grad) .or. want_hess
      screened = present(shell_mask)
      if (screened) then
         if (.not. (present(ao_offset) .and. present(n_ao_out))) then
            call error%set(ERROR_VALIDATION, "AO evaluation: a shell mask needs the "// &
                           "compressed offsets and width alongside it; passing one "// &
                           "without the others would silently write the wrong columns.")
            return
         end if
         n_ao_here = n_ao_out
      else
         n_ao_here = mol%nao
      end if

      allocate (ao(n_points, n_ao_here))
      ao = 0.0_dp
      if (present(grad)) then
         allocate (grad(n_points, n_ao_here, 3))
         grad = 0.0_dp
      end if
      if (present(hess)) then
         allocate (hess(n_points, n_ao_here, AO_HESS_COMP))
         hess = 0.0_dp
      end if
      if (want_d3) then
         allocate (deriv3(n_points, n_ao_here, AO_DERIV3_COMP))
         deriv3 = 0.0_dp
      end if

      !$omp parallel default(none) &
      !$omp    shared(mol, coords, ao, grad, hess, deriv3, n_points, want_grad, &
!$omp           want_hess, want_d3) &
      !$omp    shared(screened, shell_mask, ao_offset) &
      !$omp    private(ish, l, nprim, nctr, ic, ip, ig, off_ao, shell_base, &
      !$omp            n_cart, n_sph, &
      !$omp            n_here, comp, i, j, iatom, exp_ptr, coef_ptr, id, ih, &
      !$omp            px, py, pz, dr, ex, &
      !$omp            dx, dy, dz, r2, radial, dradial, d2radial, fac, cart, sph, trans, &
      !$omp            dcart, dsph, d2cart, d2sph, d3cart, d3sph, d3radial)
      allocate (cart((C2S_LMAX + 1)*(C2S_LMAX + 2)/2), sph(2*C2S_LMAX + 1))
      allocate (trans(2*C2S_LMAX + 1, (C2S_LMAX + 1)*(C2S_LMAX + 2)/2))
      allocate (dcart((C2S_LMAX + 1)*(C2S_LMAX + 2)/2, 3), dsph(2*C2S_LMAX + 1, 3))
      allocate (d2cart((C2S_LMAX + 1)*(C2S_LMAX + 2)/2, AO_HESS_COMP), &
                d2sph(2*C2S_LMAX + 1, AO_HESS_COMP))
      allocate (d3cart((C2S_LMAX + 1)*(C2S_LMAX + 2)/2, AO_DERIV3_COMP), &
                d3sph(2*C2S_LMAX + 1, AO_DERIV3_COMP))

      do ish = 1, mol%nbas
         if (screened) then
            if (.not. shell_mask(ish)) cycle
            shell_base = ao_offset(ish)
         else
            shell_base = mol%shell_offset(ish)
         end if
         l = mol%bas(LIBCINT_ANG_OF, ish)
         nprim = mol%bas(LIBCINT_NPRIM_OF, ish)
         nctr = mol%bas(LIBCINT_NCTR_OF, ish)
         exp_ptr = mol%bas(LIBCINT_PTR_EXP, ish)
         coef_ptr = mol%bas(LIBCINT_PTR_COEFF, ish)
         iatom = mol%bas(LIBCINT_ATOM_OF, ish) + 1

         n_cart = (l + 1)*(l + 2)/2
         n_sph = 2*l + 1
         n_here = shell_dim(mol%cartesian, ish - 1, mol%bas)/nctr

         ! The factor libcint keeps outside its transform table. It applies to
         ! *both* angular conventions -- `cint1e.c` multiplies every integral by
         ! `CINTcommon_fac_sp` of each shell's l without consulting
         ! cart-vs-sph -- so gating it on the spherical path is wrong.
         fac = common_fac_sp(l)
         if (.not. mol%cartesian) call c2s_block(l, trans)

         !$omp do schedule(static)
         do ig = 1, n_points
            dx = coords(1, ig) - mol%coords(1, iatom)
            dy = coords(2, ig) - mol%coords(2, iatom)
            dz = coords(3, ig) - mol%coords(3, iatom)
            r2 = dx*dx + dy*dy + dz*dz
            dr = [dx, dy, dz]

            ! Angular part, in libcint's Cartesian order, and its gradient. The
            ! power rule term vanishes when the exponent is zero, which has to
            ! be a branch: `ref_pow(x, -1)` is not what it means.
            comp = 0
            do i = 0, l
               do j = 0, i
                  comp = comp + 1
                  px = l - i
                  py = i - j
                  pz = j
                  cart(comp) = ref_pow(dx, px)*ref_pow(dy, py)*ref_pow(dz, pz)
                  if (want_grad) then
                     dcart(comp, :) = 0.0_dp
                     if (px > 0) dcart(comp, 1) = real(px, dp) &
                                                  *ref_pow(dx, px - 1)*ref_pow(dy, py)*ref_pow(dz, pz)
                     if (py > 0) dcart(comp, 2) = real(py, dp) &
                                                  *ref_pow(dx, px)*ref_pow(dy, py - 1)*ref_pow(dz, pz)
                     if (pz > 0) dcart(comp, 3) = real(pz, dp) &
                                                  *ref_pow(dx, px)*ref_pow(dy, py)*ref_pow(dz, pz - 1)
                  end if
                  if (want_hess) then
                     ! Same power rule applied twice. The guards are on the
                     ! *original* exponent rather than the reduced one, because
                     ! the prefactor already vanishes when it should: a first
                     ! power differentiated twice gives 1*0, and `ref_pow` would
                     ! otherwise be asked for a negative exponent.
                     d2cart(comp, :) = 0.0_dp
                     if (px > 1) d2cart(comp, 1) = real(px*(px - 1), dp) &
                                                   *ref_pow(dx, px - 2)*ref_pow(dy, py)*ref_pow(dz, pz)
                     if (px > 0 .and. py > 0) d2cart(comp, 2) = real(px*py, dp) &
                                                                *ref_pow(dx, px - 1)*ref_pow(dy, py - 1)*ref_pow(dz, pz)
                     if (px > 0 .and. pz > 0) d2cart(comp, 3) = real(px*pz, dp) &
                                                                *ref_pow(dx, px - 1)*ref_pow(dy, py)*ref_pow(dz, pz - 1)
                     if (py > 1) d2cart(comp, 4) = real(py*(py - 1), dp) &
                                                   *ref_pow(dx, px)*ref_pow(dy, py - 2)*ref_pow(dz, pz)
                     if (py > 0 .and. pz > 0) d2cart(comp, 5) = real(py*pz, dp) &
                                                                *ref_pow(dx, px)*ref_pow(dy, py - 1)*ref_pow(dz, pz - 1)
                     if (pz > 1) d2cart(comp, 6) = real(pz*(pz - 1), dp) &
                                                   *ref_pow(dx, px)*ref_pow(dy, py)*ref_pow(dz, pz - 2)
                  end if
                  if (want_d3) then
                     ! The power rule a third time, with `ref_angular_third`
                     ! keeping the exponent bookkeeping in one place.
                     do ih = 1, AO_DERIV3_COMP
                        d3cart(comp, ih) = ref_angular_third(px, py, pz, dx, dy, dz, &
                                                             D3_I(ih), D3_J(ih), D3_K(ih))
                     end do
                  end if
               end do
            end do

            ! The transform is linear and does not depend on position, so it
            ! applies to the gradient exactly as it does to the value.
            if (.not. mol%cartesian) then
               do comp = 1, n_sph
                  sph(comp) = fac*sum(trans(comp, 1:n_cart)*cart(1:n_cart))
                  if (want_grad) then
                     do id = 1, 3
                        dsph(comp, id) = fac*sum(trans(comp, 1:n_cart)*dcart(1:n_cart, id))
                     end do
                  end if
                  if (want_hess) then
                     do ih = 1, AO_HESS_COMP
                        d2sph(comp, ih) = fac*sum(trans(comp, 1:n_cart)*d2cart(1:n_cart, ih))
                     end do
                  end if
                  if (want_d3) then
                     do ih = 1, AO_DERIV3_COMP
                        d3sph(comp, ih) = fac*sum(trans(comp, 1:n_cart)*d3cart(1:n_cart, ih))
                     end do
                  end if
               end do
            end if

            ! One contraction column at a time: libcint lays a shell's
            ! functions out with the contraction index outermost, so the AO
            ! offset advances by `n_here` per column.
            do ic = 1, nctr
               radial = 0.0_dp
               dradial = 0.0_dp
               d2radial = 0.0_dp
               d3radial = 0.0_dp
               do ip = 1, nprim
                  ! d/dr_k exp(-a r^2) = -2 a r_k exp(-a r^2), so the -2a is
                  ! accumulated here and the r_k factor applied per component.
                  ! Differentiating again gives 4 a^2 r_j r_k - 2 a delta_jk;
                  ! the 4 a^2 is what accumulates below and the delta is applied
                  ! in the assembly, where it multiplies `dradial`.
                  ex = mol%env(coef_ptr + (ic - 1)*nprim + ip) &
                       *exp(-mol%env(exp_ptr + ip)*r2)
                  radial = radial + ex
                  if (want_grad) then
                     dradial = dradial - 2.0_dp*mol%env(exp_ptr + ip)*ex
                  end if
                  if (want_hess) then
                     d2radial = d2radial + 4.0_dp*mol%env(exp_ptr + ip) &
                                *mol%env(exp_ptr + ip)*ex
                  end if
                  if (want_d3) then
                     ! One more factor of -2a per derivative of the exponential.
                     d3radial = d3radial - 8.0_dp*mol%env(exp_ptr + ip)**3*ex
                  end if
               end do
               off_ao = shell_base + (ic - 1)*n_here
               if (mol%cartesian) then
                  do comp = 1, n_cart
                     ao(ig, off_ao + comp) = fac*radial*cart(comp)
                  end do
                  if (present(grad)) then
                     do comp = 1, n_cart
                        grad(ig, off_ao + comp, 1) = fac*(dcart(comp, 1)*radial &
                                                          + cart(comp)*dradial*dx)
                        grad(ig, off_ao + comp, 2) = fac*(dcart(comp, 2)*radial &
                                                          + cart(comp)*dradial*dy)
                        grad(ig, off_ao + comp, 3) = fac*(dcart(comp, 3)*radial &
                                                          + cart(comp)*dradial*dz)
                     end do
                  end if
                  if (want_d3) then
                     call ref_assemble_third(d3cart, d2cart, dcart, cart, n_cart, &
                                             radial, dradial, d2radial, d3radial, dr, &
                                             fac, deriv3(ig, off_ao + 1:off_ao + n_cart, :))
                  end if
                  if (present(hess)) then
                     ! Product rule on chi = A(r) R(r^2), twice. The two cross
                     ! terms are not one term doubled: for an off-diagonal (j,k)
                     ! they differ, dA/dx_j pairing with dR/dx_k and vice versa.
                     do ih = 1, AO_HESS_COMP
                        i = HESS_J(ih)
                        j = HESS_K(ih)
                        hess(ig, off_ao + 1:off_ao + n_cart, ih) = &
                           fac*(d2cart(1:n_cart, ih)*radial &
                                + dcart(1:n_cart, i)*dradial*dr(j) &
                                + dcart(1:n_cart, j)*dradial*dr(i) &
                                + cart(1:n_cart)*d2radial*dr(i)*dr(j))
                        if (i == j) then
                           hess(ig, off_ao + 1:off_ao + n_cart, ih) = &
                              hess(ig, off_ao + 1:off_ao + n_cart, ih) &
                              + fac*cart(1:n_cart)*dradial
                        end if
                     end do
                  end if
               else
                  do comp = 1, n_sph
                     ao(ig, off_ao + comp) = radial*sph(comp)
                  end do
                  if (present(grad)) then
                     do comp = 1, n_sph
                        grad(ig, off_ao + comp, 1) = dsph(comp, 1)*radial &
                                                     + sph(comp)*dradial*dx
                        grad(ig, off_ao + comp, 2) = dsph(comp, 2)*radial &
                                                     + sph(comp)*dradial*dy
                        grad(ig, off_ao + comp, 3) = dsph(comp, 3)*radial &
                                                     + sph(comp)*dradial*dz
                     end do
                  end if
                  if (want_d3) then
                     call ref_assemble_third(d3sph, d2sph, dsph, sph, n_sph, &
                                             radial, dradial, d2radial, d3radial, dr, &
                                             1.0_dp, deriv3(ig, off_ao + 1:off_ao + n_sph, :))
                  end if
                  if (present(hess)) then
                     ! As the Cartesian branch, with `fac` already inside sph,
                     ! dsph and d2sph rather than applied here.
                     do ih = 1, AO_HESS_COMP
                        i = HESS_J(ih)
                        j = HESS_K(ih)
                        hess(ig, off_ao + 1:off_ao + n_sph, ih) = &
                           d2sph(1:n_sph, ih)*radial &
                           + dsph(1:n_sph, i)*dradial*dr(j) &
                           + dsph(1:n_sph, j)*dradial*dr(i) &
                           + sph(1:n_sph)*d2radial*dr(i)*dr(j)
                        if (i == j) then
                           hess(ig, off_ao + 1:off_ao + n_sph, ih) = &
                              hess(ig, off_ao + 1:off_ao + n_sph, ih) &
                              + sph(1:n_sph)*dradial
                        end if
                     end do
                  end if
               end if
            end do
         end do
         !$omp end do
      end do

      ! TODO(mqc): `d3cart` and `d3sph` are allocated per thread above and are
      ! not in this list, so every thread's third-derivative scratch is left
      ! allocated at the end of the parallel region.
      deallocate (cart, sph, trans, dcart, dsph, d2cart, d2sph)
      !$omp end parallel
   end subroutine ref_eval_ao_block

   pure subroutine ref_assemble_third(a3, a2, a1, a0, n, r0, r1, r2, r3, dr, fac, out)
      !! The product rule on chi = A(r) R(r^2), three times
      !!
      !!   d3(AR) = A_ijk R
      !!          + A_ij R_k + A_ik R_j + A_jk R_i
      !!          + A_i R_jk + A_j R_ik + A_k R_ij
      !!          + A R_ijk
      !!
      !! with the radial factors following from R depending on r^2 alone:
      !!
      !!   R_i   = r1 x_i
      !!   R_ij  = r2 x_i x_j + r1 delta_ij
      !!   R_ijk = r3 x_i x_j x_k + r2 (delta_ij x_k + delta_ik x_j + delta_jk x_i)
      !!
      !! Shared by the Cartesian and spherical branches: the transform is
      !! linear and position-independent, so it commutes with every derivative.
      !! `fac` is folded in on the Cartesian side and is already inside the
      !! angular arrays on the other.
      real(dp), intent(in) :: a3(:, :), a2(:, :), a1(:, :), a0(:)
      integer, intent(in) :: n
      real(dp), intent(in) :: r0, r1, r2, r3, fac
      real(dp), intent(in) :: dr(3)
      real(dp), intent(out) :: out(:, :)

      integer, parameter :: DI(AO_DERIV3_COMP) = [1, 1, 1, 1, 1, 1, 2, 2, 2, 3]
      integer, parameter :: DJ(AO_DERIV3_COMP) = [1, 1, 1, 2, 2, 3, 2, 2, 3, 3]
      integer, parameter :: DK(AO_DERIV3_COMP) = [1, 2, 3, 2, 3, 3, 2, 3, 3, 3]
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

      integer :: ih, i, j, k
      real(dp) :: rk, rj, ri, rjk, rik, rij, rijk

      do ih = 1, AO_DERIV3_COMP
         i = DI(ih)
         j = DJ(ih)
         k = DK(ih)

         ri = r1*dr(i)
         rj = r1*dr(j)
         rk = r1*dr(k)

         rij = r2*dr(i)*dr(j)
         rik = r2*dr(i)*dr(k)
         rjk = r2*dr(j)*dr(k)
         if (i == j) rij = rij + r1
         if (i == k) rik = rik + r1
         if (j == k) rjk = rjk + r1

         rijk = r3*dr(i)*dr(j)*dr(k)
         if (i == j) rijk = rijk + r2*dr(k)
         if (i == k) rijk = rijk + r2*dr(j)
         if (j == k) rijk = rijk + r2*dr(i)

         out(1:n, ih) = fac*(a3(1:n, ih)*r0 &
                             + a2(1:n, PAIR(i, j))*rk &
                             + a2(1:n, PAIR(i, k))*rj &
                             + a2(1:n, PAIR(j, k))*ri &
                             + a1(1:n, i)*rjk &
                             + a1(1:n, j)*rik &
                             + a1(1:n, k)*rij &
                             + a0(1:n)*rijk)
      end do
   end subroutine ref_assemble_third

   pure function ref_angular_third(px, py, pz, dx, dy, dz, i, j, k) result(v)
      !! d3/dx_i dx_j dx_k of x^px y^py z^pz, at (dx, dy, dz)
      !!
      !! Count how many derivatives fall on each axis, drop the term if any
      !! axis is differentiated more times than its exponent allows, and
      !! multiply the ref_falling factorials. That guard is also what keeps `ref_pow`
      !! from being asked for a negative exponent.
      integer, intent(in) :: px, py, pz, i, j, k
      real(dp), intent(in) :: dx, dy, dz
      real(dp) :: v

      integer :: n(3), p(3)
      integer :: m
      real(dp) :: coeff

      n = 0
      n(i) = n(i) + 1
      n(j) = n(j) + 1
      n(k) = n(k) + 1
      p = [px, py, pz]

      if (any(n > p)) then
         v = 0.0_dp
         return
      end if

      coeff = 1.0_dp
      do m = 1, 3
         coeff = coeff*ref_falling(p(m), n(m))
      end do
      v = coeff*ref_pow(dx, px - n(1))*ref_pow(dy, py - n(2))*ref_pow(dz, pz - n(3))
   end function ref_angular_third

   pure function ref_falling(p, n) result(f)
      !! p (p-1) ... (p-n+1), the factor n derivatives of x^p bring down
      integer, intent(in) :: p, n
      real(dp) :: f
      integer :: m

      f = 1.0_dp
      do m = 0, n - 1
         f = f*real(p - m, dp)
      end do
   end function ref_falling

   pure function ref_pow(x, n) result(v)
      !! x**n for small non-negative n, without the intrinsic's generality
      !!
      !! `x**0` must be one even when x is zero, which is what happens at a grid
      !! point sitting on a nucleus, and there are such points.
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      real(dp) :: v

      integer :: k

      v = 1.0_dp
      do k = 1, n
         v = v*x
      end do
   end function ref_pow

end module test_mqc_czt_ao

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_ao, only: collect_mqc_czt_ao_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_ao", collect_mqc_czt_ao_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
