!! Carrying a CI vector between orbital bases
module test_mqc_ci_transform
   !! The transformation is a sum of minor determinants, and the two ways it can
   !! be wrong are the two things tested here.
   !!
   !! **The magnitudes** are checked against minors worked out by hand, on
   !! spaces small enough to write down: two orbitals holding one electron,
   !! where the string transformation *is* the orbital rotation, and three
   !! holding two, where each element is a two-by-two determinant that fits on a
   !! line.
   !!
   !! **The signs** are the part that survives every plausible check but one. A
   !! transformation with a sign error is still orthogonal, still preserves the
   !! norm, and still gives an energy-invariant vector for a *symmetric*
   !! rotation. So `permutation_carries_the_fermionic_sign` uses a rotation that
   !! merely exchanges two orbitals, where the answer is a signed permutation
   !! and every sign is decidable by counting anticommutations.
   !!
   !! What is not tested here is the physics, because it does not live here: the
   !! statement that a transformed wave function has the energy it was solved at
   !! needs a Hamiltonian, and `no_sharing_analysis` asserts it at runtime on
   !! every calculation that takes the transform route.
   use, intrinsic :: iso_fortran_env, only: int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_ci_transform, only: string_rotation_matrix, transform_ci_vector, &
                               rotation_is_orthogonal, orthogonality_defect, &
                               scatter_restricted_ci
   use mqc_ormas_space, only: ormas_space_t, build_ormas_space, determinant_address, &
                              ormas_strings, describes_a_cas
   implicit none
   private

   public :: collect_mqc_ci_transform_tests

   real(dp), parameter :: TOL = 1.0e-12_dp

   !> A rotation by 30 degrees in a two-orbital space. With one electron the
   !> strings are the orbitals themselves, so the string transformation has to
   !> come back as this matrix unchanged.
   real(dp), parameter :: C30 = 0.86602540378443864676_dp   ! cos(30 degrees)
   real(dp), parameter :: S30 = 0.5_dp                      ! sin(30 degrees)
   real(dp), parameter :: ROT2(2, 2) = reshape([C30, S30, -S30, C30], [2, 2])

contains

   subroutine collect_mqc_ci_transform_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("one_electron_is_the_orbital_rotation", &
                               one_electron_is_the_orbital_rotation), &
                  new_unittest("two_electrons_are_two_by_two_minors", &
                               two_electrons_are_two_by_two_minors), &
                  new_unittest("permutation_carries_the_fermionic_sign", &
                               permutation_carries_the_fermionic_sign), &
                  new_unittest("the_identity_transforms_nothing", &
                               the_identity_transforms_nothing), &
                  new_unittest("string_matrix_is_orthogonal", &
                               string_matrix_is_orthogonal), &
                  new_unittest("a_round_trip_returns_the_vector", &
                               a_round_trip_returns_the_vector), &
                  new_unittest("a_filled_space_is_a_determinant", &
                               a_filled_space_is_a_determinant), &
                  new_unittest("a_complete_space_scatters_to_itself", &
                               a_complete_space_scatters_to_itself), &
                  new_unittest("a_restricted_space_leaves_holes", &
                               a_restricted_space_leaves_holes), &
                  new_unittest("a_scattered_vector_can_be_rotated", &
                               a_scattered_vector_can_be_rotated), &
                  new_unittest("the_defect_measures_the_miss", &
                               the_defect_measures_the_miss), &
                  new_unittest("refusals", refusals) &
                  ]
   end subroutine collect_mqc_ci_transform_tests

   subroutine one_electron_is_the_orbital_rotation(error)
      !! With one electron per spin the strings are the orbitals
      !!
      !! Every submatrix is one by one, so its determinant is the element
      !! itself and the string transformation is the orbital rotation. It fixes
      !! the *direction* of the convention: `matrix(I,J)` is row source, column
      !! target, and getting that backwards transposes everything below without
      !! disturbing any other property tested here.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: matrix(:, :)

      call string_rotation_matrix(ROT2, 1, matrix, err)
      call check(error,.not. err%has_error(), "the rotation was refused")
      if (allocated(error)) return

      call check(error, size(matrix, 1), 2, "there are two one-electron strings")
      if (allocated(error)) return
      call check(error, maxval(abs(matrix - ROT2)), 0.0_dp, thr=TOL, &
                 more="one electron does not reproduce the orbital rotation")
   end subroutine one_electron_is_the_orbital_rotation

   subroutine two_electrons_are_two_by_two_minors(error)
      !! Three orbitals, two electrons, against minors computed by hand
      !!
      !! The strings in address order are `{1,2}`, `{1,3}`, `{2,3}`, so
      !! `matrix(I,J)` is the two-by-two determinant of the rows `I` occupies
      !! against the columns `J` occupies. With a rotation that acts only in the
      !! `{1,2}` block and leaves orbital 3 alone, every one of those nine
      !! determinants is a product or a single element:
      !!
      !!     <12|12> = c*c - (-s)*s = 1        (the block is orthogonal)
      !!     <12|13> = c*0 - (-s)*0 = 0        (orbital 3 is untouched)
      !!     <13|13> = c*1 = c,  <13|23> = -s
      !!     <23|13> = s,        <23|23> = c
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: matrix(:, :)
      real(dp) :: rotation(3, 3), expected(3, 3)

      rotation = 0.0_dp
      rotation(1:2, 1:2) = ROT2
      rotation(3, 3) = 1.0_dp

      call string_rotation_matrix(rotation, 2, matrix, err)
      call check(error,.not. err%has_error(), "the rotation was refused")
      if (allocated(error)) return
      call check(error, size(matrix, 1), 3, "there are three two-electron strings")
      if (allocated(error)) return

      expected = reshape([1.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, C30, S30, &
                          0.0_dp, -S30, C30], [3, 3])
      call check(error, maxval(abs(matrix - expected)), 0.0_dp, thr=TOL, &
                 more="the two-electron minors are not the hand-computed ones")
   end subroutine two_electrons_are_two_by_two_minors

   subroutine permutation_carries_the_fermionic_sign(error)
      !! Exchanging two orbitals costs a sign, and the minors have to know
      !!
      !! The rotation swaps orbitals 1 and 2 and leaves 3 alone. Three orbitals,
      !! two electrons, strings `{1,2}`, `{1,3}`, `{2,3}`:
      !!
      !!     {1,2} -> {2,1}, which is {1,2} with the two operators exchanged,
      !!              so the coefficient is -1
      !!     {1,3} -> {2,3}, already ascending, so +1
      !!     {2,3} -> {1,3}, already ascending, so +1
      !!
      !! A transformation that dropped the sign would be a permutation matrix,
      !! orthogonal, norm-preserving and wrong. It is the one error that the
      !! guard inside `transform_ci_vector` cannot see.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: matrix(:, :)
      real(dp) :: swap(3, 3), expected(3, 3)

      swap = 0.0_dp
      swap(2, 1) = 1.0_dp     ! target orbital 1 is source orbital 2
      swap(1, 2) = 1.0_dp
      swap(3, 3) = 1.0_dp

      call string_rotation_matrix(swap, 2, matrix, err)
      call check(error,.not. err%has_error(), "the permutation was refused")
      if (allocated(error)) return

      expected = 0.0_dp
      expected(1, 1) = -1.0_dp
      expected(2, 3) = 1.0_dp
      expected(3, 2) = 1.0_dp
      call check(error, maxval(abs(matrix - expected)), 0.0_dp, thr=TOL, &
                 more="an orbital exchange did not produce the fermionic sign")
   end subroutine permutation_carries_the_fermionic_sign

   subroutine the_identity_transforms_nothing(error)
      !! No rotation, no change, in either the matrix or a vector through it
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: matrix(:, :), transformed(:, :)
      real(dp) :: identity(4, 4), ci(6, 6)
      integer :: i, j

      identity = 0.0_dp
      do i = 1, 4
         identity(i, i) = 1.0_dp
      end do

      call string_rotation_matrix(identity, 2, matrix, err)
      call check(error,.not. err%has_error(), "the identity was refused")
      if (allocated(error)) return
      call check(error, size(matrix, 1), 6, "C(4,2) is six strings")
      if (allocated(error)) return

      do j = 1, 6
         do i = 1, 6
            ci(i, j) = sin(real(5*i + 3*j, dp))
         end do
      end do
      ci = ci/sqrt(sum(ci**2))

      call transform_ci_vector(matrix, matrix, ci, transformed, err)
      call check(error,.not. err%has_error(), "the transformation was refused")
      if (allocated(error)) return
      call check(error, maxval(abs(transformed - ci)), 0.0_dp, thr=TOL, &
                 more="the identity changed the coefficients")
   end subroutine the_identity_transforms_nothing

   subroutine string_matrix_is_orthogonal(error)
      !! `M` inherits `T`'s orthogonality, for every electron count at once
      !!
      !! This is what makes the transformation a change of basis rather than an
      !! arbitrary linear map, and it holds for every number of electrons in the
      !! same space, so all four are checked against the one rotation.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: matrix(:, :)
      real(dp) :: rotation(5, 5)
      integer :: n_electrons

      call jacobi_rotation(rotation)
      call check(error, rotation_is_orthogonal(rotation), &
                 "the test's own rotation is not orthogonal")
      if (allocated(error)) return

      do n_electrons = 1, 4
         call string_rotation_matrix(rotation, n_electrons, matrix, err)
         call check(error,.not. err%has_error(), "the rotation was refused")
         if (allocated(error)) return
         call check(error, rotation_is_orthogonal(matrix), &
                    "the string transformation is not orthogonal")
         if (allocated(error)) return
         deallocate (matrix)
      end do
   end subroutine string_matrix_is_orthogonal

   subroutine a_round_trip_returns_the_vector(error)
      !! Transform by `T`, then by `T^T`, and get back what went in
      !!
      !! An independent check on the signs: the round trip is the identity only
      !! if every minor's sign agrees with its partner's, and a systematic sign
      !! error in one direction would have to be matched exactly in the other to
      !! survive this.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: forward(:, :), backward(:, :)
      real(dp), allocatable :: once(:, :), twice(:, :)
      real(dp) :: rotation(5, 5), ci(10, 10)
      integer :: i, j

      call jacobi_rotation(rotation)
      call string_rotation_matrix(rotation, 2, forward, err)
      call string_rotation_matrix(transpose(rotation), 2, backward, err)
      call check(error,.not. err%has_error(), "a rotation was refused")
      if (allocated(error)) return

      do j = 1, 10
         do i = 1, 10
            ci(i, j) = cos(real(7*i + 2*j, dp))
         end do
      end do
      ci = ci/sqrt(sum(ci**2))

      call transform_ci_vector(forward, forward, ci, once, err)
      call transform_ci_vector(backward, backward, once, twice, err)
      call check(error,.not. err%has_error(), "a transformation was refused")
      if (allocated(error)) return
      call check(error, maxval(abs(twice - ci)), 0.0_dp, thr=1.0e-11_dp, &
                 more="the round trip did not return the vector it started from")
   end subroutine a_round_trip_returns_the_vector

   subroutine a_filled_space_is_a_determinant(error)
      !! Fill every orbital and there is one string, worth `det(T)`
      !!
      !! The whole transformation collapses to a single number, which for an
      !! orthogonal rotation is plus or minus one. Worth pinning because it is
      !! the edge the loops are most likely to get wrong, and because a closed
      !! shell in a minimal space reaches it.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: matrix(:, :)
      real(dp) :: rotation(5, 5)

      call jacobi_rotation(rotation)
      call string_rotation_matrix(rotation, 5, matrix, err)
      call check(error,.not. err%has_error(), "the rotation was refused")
      if (allocated(error)) return
      call check(error, size(matrix, 1), 1, "a filled space has one string")
      if (allocated(error)) return
      call check(error, abs(matrix(1, 1)), 1.0_dp, thr=TOL, &
                 more="the determinant of an orthogonal rotation is plus or minus one")

      ! And no electrons at all is the empty determinant, which is the same
      ! object in any basis.
      deallocate (matrix)
      call string_rotation_matrix(rotation, 0, matrix, err)
      call check(error,.not. err%has_error(), "an empty space was refused")
      if (allocated(error)) return
      call check(error, matrix(1, 1), 1.0_dp, thr=TOL, &
                 more="the empty determinant is not itself in the other basis")
   end subroutine a_filled_space_is_a_determinant

   subroutine a_complete_space_scatters_to_itself(error)
      !! One subspace is a CAS, and writing it out must move nothing
      !!
      !! The degenerate case, and the one that pins the address arithmetic: with
      !! a single subspace every determinant survives, so the rectangle is full
      !! and every coefficient must land where `determinant_address` says it
      !! sits. Nothing here depends on the restriction being interesting.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(ormas_space_t) :: space
      real(dp), allocatable :: flat(:), ci(:, :)
      integer(int64), allocatable :: alpha(:), beta(:)
      integer :: ia, ib, n
      integer(int64) :: address

      call build_ormas_space([1], 4, 2, 2, [0], [4], space, err)
      call check(error,.not. err%has_error(), "the space was refused")
      if (allocated(error)) return
      call check(error, describes_a_cas(space), "one subspace should be a CAS")
      if (allocated(error)) return

      n = int(space%n_determinants)
      call check(error, n, 36, "C(4,2)^2 is thirty-six determinants")
      if (allocated(error)) return

      allocate (flat(n))
      do ia = 1, n
         flat(ia) = sin(real(ia, dp))
      end do
      flat = flat/sqrt(sum(flat**2))

      call scatter_restricted_ci(space, flat, ci, err)
      call check(error,.not. err%has_error(), "the scatter was refused")
      if (allocated(error)) return
      call check(error, size(ci, 1), 6, "there are six alpha strings")
      if (allocated(error)) return
      call check(error, sum(ci**2), 1.0_dp, thr=TOL, &
                 more="the scatter did not preserve the norm")
      if (allocated(error)) return

      call ormas_strings(space, alpha, beta, err)
      do ib = 1, size(beta)
         do ia = 1, size(alpha)
            address = determinant_address(space, ia, ib)
            call check(error, ci(ia, ib), flat(address), thr=TOL, &
                       more="a coefficient did not land at its own address")
            if (allocated(error)) return
         end do
      end do
      call space%destroy()
   end subroutine a_complete_space_scatters_to_itself

   subroutine a_restricted_space_leaves_holes(error)
      !! A genuine restriction: the excluded determinants come back zero
      !!
      !! Four orbitals split two and two, with at least three of the four
      !! electrons in the first subspace. That excludes every determinant that
      !! moves two electrons out, so the rectangle must be strictly sparser than
      !! it is wide, and the surviving count must be the space's own.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(ormas_space_t) :: space
      real(dp), allocatable :: flat(:), ci(:, :)
      integer :: n, filled

      call build_ormas_space([1, 3], 4, 2, 2, [3, 0], [4, 1], space, err)
      call check(error,.not. err%has_error(), "the space was refused")
      if (allocated(error)) return
      call check(error,.not. describes_a_cas(space), &
                 "this partition should not be a complete space")
      if (allocated(error)) return

      n = int(space%n_determinants)
      call check(error, n < 36, "a restricted space holds fewer than all 36")
      if (allocated(error)) return

      allocate (flat(n))
      flat = 1.0_dp/sqrt(real(n, dp))
      call scatter_restricted_ci(space, flat, ci, err)
      call check(error,.not. err%has_error(), "the scatter was refused")
      if (allocated(error)) return

      ! Every coefficient is the same nonzero value, so counting the nonzeros
      ! counts the determinants that were kept -- which must be the space's own
      ! count, no more and no fewer.
      filled = count(abs(ci) > 0.0_dp)
      call check(error, filled, n, &
                 "the rectangle holds a different number of coefficients than "// &
                 "the restricted space has determinants")
      if (allocated(error)) return
      call check(error, sum(ci**2), 1.0_dp, thr=TOL, &
                 more="the scatter did not preserve the norm")
      call space%destroy()
   end subroutine a_restricted_space_leaves_holes

   subroutine a_scattered_vector_can_be_rotated(error)
      !! The point of the exercise: rotate a restricted wave function
      !!
      !! A restricted space is not closed under rotation of its active
      !! orbitals, so this is the only way the rotation can be done at all --
      !! and the result is expected to be *denser* than what went in, because
      !! the rotated wave function does not live in the restricted space any
      !! more. The norm is what has to survive, not the sparsity.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(ormas_space_t) :: space
      real(dp), allocatable :: flat(:), ci(:, :), matrix(:, :), turned(:, :)
      real(dp) :: rotation(4, 4)
      integer :: n, i, j

      call build_ormas_space([1, 3], 4, 2, 2, [3, 0], [4, 1], space, err)
      n = int(space%n_determinants)
      allocate (flat(n))
      do i = 1, n
         flat(i) = cos(real(2*i, dp))
      end do
      flat = flat/sqrt(sum(flat**2))
      call scatter_restricted_ci(space, flat, ci, err)
      call check(error,.not. err%has_error(), "the scatter was refused")
      if (allocated(error)) return

      ! A rotation that mixes the two subspaces, which is exactly the kind the
      ! restriction is not invariant under.
      rotation = 0.0_dp
      do i = 1, 4
         rotation(i, i) = 1.0_dp
      end do
      rotation(2, 2) = C30
      rotation(3, 3) = C30
      rotation(2, 3) = -S30
      rotation(3, 2) = S30

      call string_rotation_matrix(rotation, 2, matrix, err)
      call transform_ci_vector(matrix, matrix, ci, turned, err)
      call check(error,.not. err%has_error(), "the rotation was refused")
      if (allocated(error)) return
      call check(error, sum(turned**2), 1.0_dp, thr=1.0e-11_dp, &
                 more="rotating the scattered vector did not preserve the norm")
      if (allocated(error)) return

      ! It should have left the restricted space: at least one determinant the
      ! restriction excluded now carries amplitude. If none does, the rotation
      ! did not mix the subspaces and this test is proving nothing.
      j = 0
      do i = 1, size(ci, 2)
         if (any(abs(ci(:, i)) == 0.0_dp .and. abs(turned(:, i)) > 1.0e-10_dp)) j = 1
      end do
      call check(error, j, 1, "the rotation should have carried amplitude outside "// &
                 "the restricted space, and did not")
      call space%destroy()
   end subroutine a_scattered_vector_can_be_rotated

   subroutine the_defect_measures_the_miss(error)
      !! The defect is a number, not only a verdict
      !!
      !! A caller deciding whether an offered wave function lives in this space
      !! needs to tell "the same space, reached by a different route" from "a
      !! different space of the same dimension". Both fail the boolean; only the
      !! magnitude separates them, and the separation is several orders.
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: rotation(5, 5), nudged(5, 5), oblong(5, 3)
      real(dp), parameter :: NUDGE = 1.0e-6_dp
      integer :: i

      call jacobi_rotation(rotation)
      call check(error, orthogonality_defect(rotation) < 1.0e-14_dp, &
                 "an orthogonal rotation should have essentially no defect")
      if (allocated(error)) return

      ! Stretched by `NUDGE` along one axis. Built from the identity rather than
      ! from `rotation` because then the answer is exact and derivable in a
      ! line: the only element of `R^T R - I` that moves is
      !
      !     (1 + NUDGE)^2 - 1 = 2*NUDGE + NUDGE^2
      !
      ! whereas perturbing a general rotation scales the defect by the
      ! magnitude of the element perturbed, which makes the expected value a
      ! property of the test's own matrix rather than of the routine.
      nudged = 0.0_dp
      do i = 1, 5
         nudged(i, i) = 1.0_dp
      end do
      nudged(1, 1) = 1.0_dp + NUDGE
      call check(error, orthogonality_defect(nudged), 2.0_dp*NUDGE + NUDGE**2, &
                 thr=1.0e-15_dp, &
                 more="a small perturbation should give a defect of its own size")
      if (allocated(error)) return
      call check(error,.not. rotation_is_orthogonal(nudged), &
                 "a defect above the tolerance is still not orthogonal")
      if (allocated(error)) return

      ! Not square at all: no defect is defined, and returning zero would read
      ! as a pass.
      oblong = 0.0_dp
      call check(error, orthogonality_defect(oblong) > 1.0e30_dp, &
                 "a non-square matrix has no orthogonality to measure")
   end subroutine the_defect_measures_the_miss

   subroutine refusals(error)
      !! What the module declines, and why each would otherwise go unnoticed
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: matrix(:, :), other(:, :), transformed(:, :)
      real(dp) :: stretched(3, 3), rotation(5, 5), ci(10, 10)

      ! A rotation that is not orthogonal does not relate two orthonormal
      ! determinant bases, so its minors transform nothing. Caught here rather
      ! than left to show up as an energy that moved.
      stretched = 0.0_dp
      stretched(1, 1) = 2.0_dp
      stretched(2, 2) = 1.0_dp
      stretched(3, 3) = 1.0_dp
      call string_rotation_matrix(stretched, 1, matrix, err)
      call check(error, err%has_error(), &
                 "a non-orthogonal rotation should be refused")
      if (allocated(error)) return
      call err%clear()

      ! A transformation built for a different active space. The shapes are the
      ! only evidence, and without this check the matrix products would either
      ! crash or silently use the wrong elements.
      call jacobi_rotation(rotation)
      call string_rotation_matrix(rotation, 2, matrix, err)
      call string_rotation_matrix(rotation, 1, other, err)
      call check(error,.not. err%has_error(), "a rotation was refused")
      if (allocated(error)) return
      ci = 0.0_dp
      ci(1, 1) = 1.0_dp
      call transform_ci_vector(other, other, ci, transformed, err)
      call check(error, err%has_error(), &
                 "a transformation of the wrong size should be refused")
      if (allocated(error)) return
      call err%clear()
   end subroutine refusals

   subroutine jacobi_rotation(rotation)
      !! A five-orbital orthogonal matrix that mixes every orbital with every
      !! other, built as a product of plane rotations so that it is orthogonal
      !! by construction rather than by a normalization that could hide an error
      real(dp), intent(out) :: rotation(5, 5)

      real(dp) :: c, s, row_i(5), row_j(5), angle
      integer :: i, j

      rotation = 0.0_dp
      do i = 1, 5
         rotation(i, i) = 1.0_dp
      end do
      do i = 1, 5
         do j = i + 1, 5
            angle = 0.1_dp*real(3*i + j, dp)
            c = cos(angle)
            s = sin(angle)
            row_i = c*rotation(i, :) + s*rotation(j, :)
            row_j = -s*rotation(i, :) + c*rotation(j, :)
            rotation(i, :) = row_i
            rotation(j, :) = row_j
         end do
      end do
   end subroutine jacobi_rotation

end module test_mqc_ci_transform

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_ci_transform, only: collect_mqc_ci_transform_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_ci_transform", collect_mqc_ci_transform_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
