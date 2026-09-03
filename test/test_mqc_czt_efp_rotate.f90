!! Placing a fragment in an arbitrary orientation
module test_mqc_czt_efp_rotate
   !! The rotation machinery, checked on its own rather than only through an energy.
   !!
   !! `test_mqc_czt_efp_energy` asserts that the whole interaction is invariant under a
   !! rigid rotation, which is the statement that matters -- but it is one number at
   !! the end of a long chain, so when it breaks it does not say where. These are the
   !! pieces, each against a property that pins it without a reference value:
   !!
   !!   * `superpose` must recover a transform that was applied deliberately.
   !!   * `cartesian_rotation` must be a representation: `D(R1) D(R2) = D(R1 R2)`,
   !!     `D(I) = I`, and for `l = 1` it must *be* the rotation, since `p` functions
   !!     transform as vectors.
   !!   * `rotate_fragment` must preserve every rotational invariant the fragment
   !!     carries, and must leave the orbitals orthonormal in the rotated basis.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_czt_efp_potential, only: efp_potential_t, make_efp_potential, &
                                    write_efp_potential, from_gamess_ao_order
   use mqc_czt_efp_read, only: efp_fragment_t, read_efp_potential
   use mqc_czt_efp_rotate, only: superpose, rotate_fragment, cartesian_rotation
   use mqc_czt_efp_pair, only: fragment_molecule
   use mqc_czt_integrals, only: czt_molecule_t
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_czt_efp_rotate_tests

   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp

   !! Built once by `water_fragment` and copied thereafter: it is an SCF plus
   !! the response solves behind the polarizabilities, identical every time,
   !! and every test here wants one.
   type(efp_fragment_t), save :: cached_water
   logical, save :: cached_ready = .false.

contains

   subroutine collect_mqc_czt_efp_rotate_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("rotate_superpose_recovers_transform", test_superpose), &
                  new_unittest("rotate_superpose_refuses_degenerate", test_degenerate), &
                  new_unittest("rotate_cartesian_is_a_representation", test_representation), &
                  new_unittest("rotate_cartesian_p_is_the_rotation", test_p_shell), &
                  new_unittest("rotate_fragment_keeps_invariants", test_invariants), &
                  new_unittest("rotate_fragment_keeps_orbitals_orthonormal", test_orbitals) &
                  ]
   end subroutine collect_mqc_czt_efp_rotate_tests

   pure function sample_rotation(a, b, c) result(rot)
      !! Three successive axis rotations, so nothing is aligned with a coordinate axis
      real(dp), intent(in) :: a, b, c
      real(dp) :: rot(3, 3)

      real(dp) :: rx(3, 3), ry(3, 3), rz(3, 3)

      rx = reshape([1.0_dp, 0.0_dp, 0.0_dp, &
                    0.0_dp, cos(a), sin(a), &
                    0.0_dp, -sin(a), cos(a)], [3, 3])
      ry = reshape([cos(b), 0.0_dp, -sin(b), &
                    0.0_dp, 1.0_dp, 0.0_dp, &
                    sin(b), 0.0_dp, cos(b)], [3, 3])
      rz = reshape([cos(c), sin(c), 0.0_dp, &
                    -sin(c), cos(c), 0.0_dp, &
                    0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
      rot = matmul(rz, matmul(ry, rx))
   end function sample_rotation

   subroutine water_fragment(frag, err)
      type(efp_fragment_t), intent(out) :: frag
      type(error_t), intent(inout) :: err

      type(efp_potential_t) :: pot
      real(dp) :: c(3, 3)
      integer :: z(3)
      character(len=2) :: symbols(3)
      character(len=*), parameter :: path = "test_efp_rotate.efp"
      integer :: unit, stat

      if (cached_ready) then
         frag = cached_water
         return
      end if
      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.00000000000000_dp, 0.00000000009155_dp, 0.10077199490609_dp, &
                   0.00000000000000_dp, 0.77250895271063_dp, -0.46780199741728_dp, &
                   0.00000000000000_dp, -0.77250895280218_dp, -0.46780199748881_dp], &
                  [3, 3])*ANG
      call make_efp_potential(z, symbols, c, "6-31g*", "WATER", pot, err)
      if (err%has_error()) return
      call write_efp_potential(pot, path, err)
      if (err%has_error()) return
      call read_efp_potential(path, frag, err)
      call pot%destroy()
      open (newunit=unit, file=path, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete")
      if (.not. err%has_error()) then
         cached_water = frag
         cached_ready = .true.
      end if
   end subroutine water_fragment

   subroutine test_superpose(error)
      !! A transform applied on purpose comes back out
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp) :: own(3, 4), deck(3, 4), rot(3, 3), got(3, 3)
      real(dp) :: shift(3), translation(3), rmsd
      integer :: k

      ! Four points in general position -- not planar, so nothing is accidentally
      ! invariant under a reflection the frame construction should have excluded.
      own = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                     1.3_dp, 0.0_dp, 0.0_dp, &
                     0.2_dp, 1.1_dp, 0.0_dp, &
                     0.4_dp, 0.3_dp, 0.9_dp], [3, 4])
      rot = sample_rotation(0.41_dp, -0.87_dp, 1.27_dp)
      shift = [2.5_dp, -1.25_dp, 0.75_dp]
      do k = 1, 4
         deck(:, k) = matmul(rot, own(:, k)) + shift
      end do

      call superpose(own, deck, got, translation, rmsd, err)
      call check(error,.not. err%has_error(), "superpose failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, maxval(abs(got - rot)) < 1.0e-13_dp, &
                 "superpose did not recover the rotation")
      if (allocated(error)) return
      call check(error, maxval(abs(translation - shift)) < 1.0e-13_dp, &
                 "superpose did not recover the translation")
      if (allocated(error)) return
      ! Exact, because the placement is rigid by construction rather than fitted.
      call check(error, rmsd < 1.0e-13_dp, "a rigid placement should leave no residual")
      if (allocated(error)) return

      ! And the fourth atom, which took no part in building the frame, must land where
      ! the transform says -- that is what makes the rmsd a real check on the rest.
      call check(error, maxval(abs(matmul(got, own(:, 4)) + translation - deck(:, 4))) &
                 < 1.0e-13_dp, "an atom outside the frame triple was not placed")
   end subroutine test_superpose

   subroutine test_degenerate(error)
      !! Placements that fix no orientation are refused rather than guessed
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp) :: own(3, 3), deck(3, 3), rot(3, 3), translation(3), rmsd
      real(dp) :: two_own(3, 2), two_deck(3, 2)

      ! Two atoms leave a free rotation about their axis, and a fragment's tensors are
      ! not invariant under it, so there is nothing to choose.
      two_own = reshape([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], [3, 2])
      two_deck = two_own
      call superpose(two_own, two_deck, rot, translation, rmsd, err)
      call check(error, err%has_error(), "a two-atom fragment should be refused")
      if (allocated(error)) return
      call err%clear()

      ! Three collinear atoms fix an axis but not a rotation about it.
      own = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                     1.0_dp, 0.0_dp, 0.0_dp, &
                     2.0_dp, 0.0_dp, 0.0_dp], [3, 3])
      deck = own
      call superpose(own, deck, rot, translation, rmsd, err)
      call check(error, err%has_error(), "collinear atoms should be refused")
   end subroutine test_degenerate

   subroutine test_representation(error)
      !! `D` is a representation of the rotation group, for every shell it handles
      !!
      !! `D(R1) D(R2) = D(R1 R2)` is the property that fails if the monomial expansion
      !! has a transposed index or a missing multinomial factor, and it holds for no
      !! reason other than the transformation being right. Checked for `d` and `f` as
      !! well as `p`, since those are where the combinatorics actually does work.
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: r1(3, 3), r2(3, 3), r12(3, 3), ident(3, 3)
      real(dp), allocatable :: d1(:, :), d2(:, :), d12(:, :), di(:, :)
      integer :: l, k

      r1 = sample_rotation(0.31_dp, 0.92_dp, -0.44_dp)
      r2 = sample_rotation(-1.05_dp, 0.23_dp, 0.68_dp)
      r12 = matmul(r1, r2)
      ident = 0.0_dp
      do k = 1, 3
         ident(k, k) = 1.0_dp
      end do

      do l = 1, 4
         call cartesian_rotation(l, r1, d1)
         call cartesian_rotation(l, r2, d2)
         call cartesian_rotation(l, r12, d12)
         call cartesian_rotation(l, ident, di)

         call check(error, maxval(abs(matmul(d1, d2) - d12)) < 1.0e-11_dp, &
                    "D(R1) D(R2) /= D(R1 R2)")
         if (allocated(error)) then
            deallocate (d1, d2, d12, di)
            return
         end if
         call check(error, maxval(abs(di - identity(size(di, 1)))) < 1.0e-13_dp, &
                    "D(identity) is not the identity")
         if (allocated(error)) then
            deallocate (d1, d2, d12, di)
            return
         end if
         deallocate (d1, d2, d12, di)
      end do
   end subroutine test_representation

   pure function identity(n) result(m)
      integer, intent(in) :: n
      real(dp) :: m(n, n)
      integer :: k

      m = 0.0_dp
      do k = 1, n
         m(k, k) = 1.0_dp
      end do
   end function identity

   subroutine test_p_shell(error)
      !! A `p` shell transforms as a vector, so its block must be the rotation itself
      !!
      !! This is the one case where the answer is known independently of the
      !! expansion, which makes it the anchor for the convention: it fixes whether the
      !! matrix acts as `R` or as its transpose, and everything above `l = 1` inherits
      !! that from the same code path.
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: rot(3, 3)
      real(dp), allocatable :: d(:, :)

      rot = sample_rotation(0.77_dp, -0.35_dp, 1.4_dp)
      call cartesian_rotation(1, rot, d)
      call check(error, size(d, 1) == 3, "a p shell should have three components")
      if (allocated(error)) return
      call check(error, maxval(abs(d - rot)) < 1.0e-13_dp, &
                 "the p block is not the rotation itself")
      deallocate (d)
   end subroutine test_p_shell

   subroutine test_invariants(error)
      !! Rotating a fragment preserves every scalar built out of what it carries
      !!
      !! Traces and isotropic averages cannot change under a rotation, so they are
      !! free checks on the tensor ranks that no reference value is needed for -- and
      !! they are sensitive to exactly the mistakes that are easy to make, since a
      !! transposed or mis-packed tensor almost never preserves its own trace.
      type(error_type), allocatable, intent(out) :: error

      type(efp_fragment_t) :: plain, spun
      type(error_t) :: err
      real(dp) :: rot(3, 3), before, after, span_before, span_after
      integer :: k

      call water_fragment(plain, err)
      call water_fragment(spun, err)
      call check(error,.not. err%has_error(), "building the fragments failed: "//err%get_full_trace())
      if (allocated(error)) return

      rot = sample_rotation(0.53_dp, -1.17_dp, 0.61_dp)
      call rotate_fragment(spun, rot, err)
      call check(error,.not. err%has_error(), "rotating failed: "//err%get_full_trace())
      if (allocated(error)) return

      ! The quadrupole trace, summed over points: a rotational invariant of a rank-2
      ! tensor, and the thing a wrong packing destroys first.
      !
      ! `plain` and `spun` are two independent SCFs, so the invariant holds only to
      ! their agreement, not to machine precision. That is usually ~2e-15, but the
      ! threaded Fock build is not bit-reproducible and about one run in thirty
      ! lands ~1e-10 apart -- so the tolerance is 1e-8, well above that and far
      ! below the O(1) shift a genuinely wrong packing would make.
      before = sum(plain%quadrupole(1, :) + plain%quadrupole(2, :) + plain%quadrupole(3, :))
      after = sum(spun%quadrupole(1, :) + spun%quadrupole(2, :) + spun%quadrupole(3, :))
      call check(error, abs(after - before) < 1.0e-8_dp, &
                 "the quadrupole trace changed under rotation")
      if (allocated(error)) return

      ! The isotropic dynamic polarizability of every orbital at every frequency,
      ! which is what E6 and E8 consume.
      before = 0.0_dp
      after = 0.0_dp
      do k = 1, size(plain%dyn_pol, 3)
         before = before + sum([(plain%dyn_pol(1, 1, k, 1) + plain%dyn_pol(2, 2, k, 1) &
                                 + plain%dyn_pol(3, 3, k, 1))])
         after = after + sum([(spun%dyn_pol(1, 1, k, 1) + spun%dyn_pol(2, 2, k, 1) &
                               + spun%dyn_pol(3, 3, k, 1))])
      end do
      call check(error, abs(after - before) < 1.0e-12_dp, &
                 "the dynamic polarizability trace changed under rotation")
      if (allocated(error)) return

      ! Distances between expansion points: a rigid rotation moves none of them.
      span_before = 0.0_dp
      span_after = 0.0_dp
      do k = 1, plain%n_points
         span_before = span_before + sqrt(sum(plain%points(:, k)**2))
         span_after = span_after + sqrt(sum(spun%points(:, k)**2))
      end do
      call check(error, abs(span_after - span_before) < 1.0e-12_dp, &
                 "the fragment's own geometry was distorted")
      if (allocated(error)) return

      ! The dipole's magnitude, which pins that it rotated as a vector rather than
      ! being left alone -- and it must actually have moved, or this proves nothing.
      before = sqrt(sum(plain%dipole(:, 1)**2))
      after = sqrt(sum(spun%dipole(:, 1)**2))
      call check(error, abs(after - before) < 1.0e-12_dp, &
                 "the dipole magnitude changed under rotation")
      if (allocated(error)) return
      call check(error, maxval(abs(spun%dipole - plain%dipole)) > 1.0e-6_dp, &
                 "the dipole did not move at all, so nothing was tested")

      call plain%destroy()
      call spun%destroy()
   end subroutine test_invariants

   subroutine test_orbitals(error)
      !! The rotated orbitals are still orthonormal, in the rotated basis
      !!
      !! `C^T S C = I` is the property that catches a wrong Cartesian `d` block, and it
      !! catches it *specifically*: the multipole ranks can all be right while the
      !! orbitals are wrong, and then exchange repulsion, charge transfer and every
      !! damped dispersion term are wrong with nothing to object. Here the overlap is
      !! built from the fragment's own rotated geometry, so this asks whether the
      !! coefficients were transformed consistently with the basis moving.
      type(error_type), allocatable, intent(out) :: error

      type(efp_fragment_t) :: spun
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), c(:, :), sc(:, :), gram(:, :)
      real(dp) :: rot(3, 3), worst
      integer :: i, j

      call water_fragment(spun, err)
      call check(error,.not. err%has_error(), "building the fragment failed: "//err%get_full_trace())
      if (allocated(error)) return

      rot = sample_rotation(-0.64_dp, 1.02_dp, 0.29_dp)
      call rotate_fragment(spun, rot, err)
      call check(error,.not. err%has_error(), "rotating failed: "//err%get_full_trace())
      if (allocated(error)) return

      call fragment_molecule(spun, [0.0_dp, 0.0_dp, 0.0_dp], mol, err)
      call check(error,.not. err%has_error(), "building the molecule failed: "//err%get_full_trace())
      if (allocated(error)) return
      call mol%overlap(s)
      call from_gamess_ao_order(mol, spun%lmo_gamess, c, err)
      call check(error,.not. err%has_error(), "converting the orbitals failed: "//err%get_full_trace())
      if (allocated(error)) return

      allocate (sc(size(s, 1), size(c, 2)), gram(size(c, 2), size(c, 2)))
      call pic_gemm(s, c, sc)
      call pic_gemm(c, sc, gram, transa="T")

      worst = 0.0_dp
      do j = 1, size(gram, 2)
         do i = 1, size(gram, 1)
            if (i == j) then
               worst = max(worst, abs(gram(i, j) - 1.0_dp))
            else
               worst = max(worst, abs(gram(i, j)))
            end if
         end do
      end do
      ! The same 1e-6 the unrotated orbitals are held to: the tolerance is the
      ! potential's ten printed decimals, not the rotation.
      call check(error, worst < 1.0e-6_dp, &
                 "the rotated orbitals are not orthonormal in the rotated basis")

      call mol%destroy()
      call spun%destroy()
      deallocate (s, c, sc, gram)
   end subroutine test_orbitals

end module test_mqc_czt_efp_rotate

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_efp_rotate, only: collect_mqc_czt_efp_rotate_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_efp_rotate", collect_mqc_czt_efp_rotate_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
