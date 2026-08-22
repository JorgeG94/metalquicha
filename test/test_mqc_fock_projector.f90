!! Unit tests for the frozen-orbital Fock constraint
module test_mqc_fock_projector
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t
   use mqc_scf_common, only: build_orthogonalizer
   use mqc_fock_projector, only: fock_projector_t, build_frozen_basis
   implicit none

   !> The constraint is four GEMMs and an exact inverse, so what is left is
   !> rounding. But it is rounding on a matrix whose largest element is the
   !> level shift, and the back transform spreads that over every block -- so a
   !> comparison against zero has to be scaled by the magnitude in play, not
   !> stated in absolute Hartree. `ATOL` is that scaled tolerance.
   real(dp), parameter :: TOL = 1.0e-12_dp

   integer, parameter :: N = 6

   !> A modest level shift, and deliberately not the 1e6 the reference codes
   !> use. GAFO decouples the blocks rather than penalising them, so the shift
   !> only has to lift the frozen virtuals clear of the occupied manifold --
   !> it does not have to overwhelm a coupling, because there is no coupling
   !> left. See `projector_shift_sets_the_precision_floor` for what asking for
   !> more costs.
   real(dp), parameter :: SHIFT = 1.0e3_dp
   real(dp), parameter :: ATOL = TOL*SHIFT

   private
   public :: collect_mqc_fock_projector

contains

   subroutine collect_mqc_fock_projector(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("projector_nothing_frozen_is_bit_identical", test_inactive), &
                  new_unittest("projector_round_trip_returns_the_fock", test_round_trip), &
                  new_unittest("projector_cuts_every_frozen_coupling", test_decoupled), &
                  new_unittest("projector_frozen_virtuals_are_degenerate", test_degenerate), &
                  new_unittest("projector_leaves_the_unfrozen_block_alone", test_unfrozen_intact), &
                  new_unittest("projector_shift_sets_the_precision_floor", test_shift_floor), &
                  new_unittest("projector_refuses_blocks_that_do_not_nest", test_bad_blocks), &
                  new_unittest("frozen_basis_is_orthonormal_in_the_metric", test_basis_orthonormal), &
                  new_unittest("frozen_basis_spans_what_was_frozen", test_basis_spans), &
                  new_unittest("frozen_basis_keeps_the_blocks_apart", test_basis_blocks), &
                  new_unittest("frozen_basis_refuses_a_repeated_orbital", test_basis_repeat), &
                  new_unittest("frozen_basis_feeds_the_projector", test_basis_end_to_end) &
                  ]
   end subroutine collect_mqc_fock_projector

   subroutine model_problem(s, f, c, n_mo, err)
      !! A conditioned overlap, a symmetric Fock, and an S-orthonormal basis
      !!
      !! `S(i,j) = 0.7^|i-j|` is positive definite and well enough conditioned
      !! that no combination is dropped, so `n_mo == n_ao` and the round trip
      !! below is exact rather than exact-on-the-retained-space.
      real(dp), allocatable, intent(out) :: s(:, :), f(:, :), c(:, :)
      integer, intent(out) :: n_mo
      type(error_t), intent(inout) :: err
      integer :: i, j

      allocate (s(N, N), f(N, N))
      do i = 1, N
         do j = 1, N
            s(i, j) = 0.7_dp**abs(i - j)
            f(i, j) = sin(real(i*j, dp)) + cos(real(i + j, dp))
         end do
      end do
      f = 0.5_dp*(f + transpose(f))
      do i = 1, N
         f(i, i) = f(i, i) - real(i, dp)
      end do

      call build_orthogonalizer(s, c, n_mo, err)
   end subroutine model_problem

   function to_mo(c, f) result(f_mo)
      !! `C^T F C`, the view the constraint is stated in
      real(dp), intent(in) :: c(:, :), f(:, :)
      real(dp), allocatable :: f_mo(:, :), work(:, :)
      integer :: n_ao, n_mo

      n_ao = size(c, 1)
      n_mo = size(c, 2)
      allocate (work(n_ao, n_mo), f_mo(n_mo, n_mo))
      call pic_gemm(f, c, work)
      call pic_gemm(c, work, f_mo, transa="T")
   end function to_mo

   subroutine test_inactive(error)
      !! No frozen orbitals must not mean "the identity to roundoff"
      type(error_type), allocatable, intent(out) :: error
      type(fock_projector_t) :: proj
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), f(:, :), c(:, :), before(:, :)
      integer :: n_mo

      call model_problem(s, f, c, n_mo, err)
      before = f
      call proj%init(c, s, 0, 0, SHIFT, err)
      call check(error,.not. err%has_error(), "init rejected an empty partition")
      if (allocated(error)) return

      call proj%apply(f, err)
      call check(error,.not. err%has_error(), "apply failed with nothing frozen")
      if (allocated(error)) return
      call check(error, all(f == before), "a projector with nothing frozen moved the Fock")
   end subroutine test_inactive

   subroutine test_round_trip(error)
      !! Freeze everything as occupied: no block is zeroed, so only the
      !! transform runs and it must give the Fock back. This is the test that
      !! catches a transpose convention before any physics sits on top of it.
      type(error_type), allocatable, intent(out) :: error
      type(fock_projector_t) :: proj
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), f(:, :), c(:, :), before(:, :)
      integer :: n_mo

      call model_problem(s, f, c, n_mo, err)
      before = f
      call proj%init(c, s, n_mo, n_mo, SHIFT, err)
      call proj%apply(f, err)
      call check(error,.not. err%has_error(), "apply failed")
      if (allocated(error)) return
      call check(error, maxval(abs(f - before)) < TOL, &
                 "the round trip through the frozen basis did not return the Fock")
   end subroutine test_round_trip

   subroutine test_decoupled(error)
      !! Every frozen-to-anything-else element is gone
      type(error_type), allocatable, intent(out) :: error
      type(fock_projector_t) :: proj
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), f(:, :), c(:, :), f_mo(:, :)
      integer :: n_mo, nfo, nfr

      call model_problem(s, f, c, n_mo, err)
      nfo = 2
      nfr = 4
      call proj%init(c, s, nfo, nfr, SHIFT, err)
      call proj%apply(f, err)
      call check(error,.not. err%has_error(), "apply failed")
      if (allocated(error)) return

      f_mo = to_mo(c, f)
      call check(error, maxval(abs(f_mo(nfo + 1:, :nfo))) < ATOL, &
                 "the frozen-occupied block is still coupled to something")
      if (allocated(error)) return
      call check(error, maxval(abs(f_mo(:nfo, nfo + 1:))) < ATOL, &
                 "the constraint is not symmetric")
      if (allocated(error)) return
      call check(error, maxval(abs(f_mo(nfr + 1:, nfo + 1:nfr))) < ATOL, &
                 "a frozen virtual is still coupled to the unfrozen space")
   end subroutine test_decoupled

   subroutine test_degenerate(error)
      !! The frozen-virtual block comes out as `shift` times the identity
      type(error_type), allocatable, intent(out) :: error
      type(fock_projector_t) :: proj
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), f(:, :), c(:, :), f_mo(:, :)
      integer :: n_mo, nfo, nfr, i, j

      call model_problem(s, f, c, n_mo, err)
      nfo = 2
      nfr = 4
      call proj%init(c, s, nfo, nfr, SHIFT, err)
      call proj%apply(f, err)
      f_mo = to_mo(c, f)

      do i = nfo + 1, nfr
         do j = nfo + 1, nfr
            if (i == j) then
               call check(error, abs(f_mo(i, i) - SHIFT) < ATOL, &
                          "a frozen virtual is not held at the shift")
            else
               call check(error, abs(f_mo(i, j)) < ATOL, &
                          "the frozen-virtual block still rotates internally")
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine test_degenerate

   subroutine test_unfrozen_intact(error)
      !! What is variational stays variational, untouched
      type(error_type), allocatable, intent(out) :: error
      type(fock_projector_t) :: proj
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), f(:, :), c(:, :), before(:, :), after(:, :)
      integer :: n_mo, nfo, nfr

      call model_problem(s, f, c, n_mo, err)
      nfo = 2
      nfr = 4
      before = to_mo(c, f)
      call proj%init(c, s, nfo, nfr, SHIFT, err)
      call proj%apply(f, err)
      after = to_mo(c, f)

      call check(error, maxval(abs(after(nfr + 1:, nfr + 1:) - before(nfr + 1:, nfr + 1:))) < ATOL, &
                 "the unfrozen block was modified")
      if (allocated(error)) return
      call check(error, maxval(abs(after(:nfo, :nfo) - before(:nfo, :nfo))) < ATOL, &
                 "the frozen-occupied block was modified rather than decoupled")
   end subroutine test_unfrozen_intact

   subroutine test_shift_floor(error)
      !! What the level shift costs the part of the matrix nobody froze
      !!
      !! The constraint is applied in the frozen-orbital basis and carried back
      !! with `(SC) F_mo (SC)^T`, and that back transform mixes the shifted
      !! block into every element. So the unfrozen block -- where the physics
      !! is -- is clean only to about `shift * epsilon`. At the 1e6 the
      !! reference implementations use that is about 1e-10 Hartree sitting
      !! under a suite that validates energies at 1e-9, which is far too close
      !! to be comfortable.
      !!
      !! This is an argument specific to GAFO being right. A penalty scheme
      !! needs a huge constant, because it is fighting a coupling it never
      !! removed; decoupling first means a modest shift does the whole job, and
      !! the precision comes back.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: leak_modest, leak_huge

      call unfrozen_leak(1.0e3_dp, leak_modest)
      call unfrozen_leak(1.0e6_dp, leak_huge)

      call check(error, leak_modest < 1.0e-11_dp, &
                 "a modest shift is already contaminating the unfrozen block")
      if (allocated(error)) return
      call check(error, leak_huge < 1.0e-8_dp, &
                 "the contamination is worse than shift*epsilon explains")
      if (allocated(error)) return
      call check(error, leak_huge > leak_modest, &
                 "the contamination does not track the shift, so this analysis is wrong")
   end subroutine test_shift_floor

   subroutine unfrozen_leak(shift, leak)
      !! How far the unfrozen block moves when only frozen blocks were touched
      real(dp), intent(in) :: shift
      real(dp), intent(out) :: leak
      type(fock_projector_t) :: proj
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), f(:, :), c(:, :), before(:, :), after(:, :)
      integer :: n_mo, nfo, nfr

      call model_problem(s, f, c, n_mo, err)
      nfo = 2
      nfr = 4
      before = to_mo(c, f)
      call proj%init(c, s, nfo, nfr, shift, err)
      call proj%apply(f, err)
      after = to_mo(c, f)
      leak = maxval(abs(after(nfr + 1:, nfr + 1:) - before(nfr + 1:, nfr + 1:)))
   end subroutine unfrozen_leak

   subroutine raw_frozen(c, raw)
      !! Three orbitals to freeze, deliberately neither orthogonal nor normalised
      !!
      !! Built from an orthonormal basis and then mixed, which is what arrives
      !! from a localization: recognisable orbitals with no orthogonality left
      !! between them. The third overlaps the first, so the block projection has
      !! something to do.
      real(dp), intent(in) :: c(:, :)
      real(dp), allocatable, intent(out) :: raw(:, :)

      allocate (raw(size(c, 1), 3))
      raw(:, 1) = c(:, 1) + 0.4_dp*c(:, 2)
      raw(:, 2) = c(:, 2) - 0.3_dp*c(:, 3)
      raw(:, 3) = c(:, 4) + 0.2_dp*c(:, 1)
      raw = 1.7_dp*raw
   end subroutine raw_frozen

   function metric_residual(s, span, v) result(left)
      !! How much of `v` lies outside the columns of `span`, measured in `S`
      real(dp), intent(in) :: s(:, :), span(:, :), v(:)
      real(dp) :: left
      real(dp), allocatable :: sv(:), coeff(:), r(:), sr(:)
      integer :: n_ao, k, i

      n_ao = size(v)
      k = size(span, 2)
      allocate (sv(n_ao), coeff(k), r(n_ao), sr(n_ao))
      sv = matmul(s, v)
      do i = 1, k
         coeff(i) = dot_product(span(:, i), sv)
      end do
      r = v - matmul(span, coeff)
      sr = matmul(s, r)
      left = sqrt(abs(dot_product(r, sr)))
   end function metric_residual

   subroutine test_basis_orthonormal(error)
      !! `C^T S C = I` -- the property `apply` is built on
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), f(:, :), c(:, :), raw(:, :), basis(:, :), gram(:, :)
      integer :: n_mo, n_basis, i

      call model_problem(s, f, c, n_mo, err)
      call raw_frozen(c, raw)
      call build_frozen_basis(raw, 2, s, basis, n_basis, err)
      call check(error,.not. err%has_error(), "building the frozen basis failed")
      if (allocated(error)) return
      call check(error, n_basis, n_mo, "the frozen basis lost part of the space")
      if (allocated(error)) return

      gram = matmul(transpose(basis), matmul(s, basis))
      do i = 1, n_basis
         gram(i, i) = gram(i, i) - 1.0_dp
      end do
      call check(error, maxval(abs(gram)) < TOL, &
                 "the frozen basis is not orthonormal in the overlap metric")
   end subroutine test_basis_orthonormal

   subroutine test_basis_spans(error)
      !! Every orbital handed in still lies in the leading columns
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), f(:, :), c(:, :), raw(:, :), basis(:, :)
      integer :: n_mo, n_basis, i

      call model_problem(s, f, c, n_mo, err)
      call raw_frozen(c, raw)
      call build_frozen_basis(raw, 2, s, basis, n_basis, err)

      do i = 1, 3
         call check(error, metric_residual(s, basis(:, :3), raw(:, i)) < TOL, &
                    "a frozen orbital is not spanned by the frozen block")
         if (allocated(error)) return
      end do
   end subroutine test_basis_spans

   subroutine test_basis_blocks(error)
      !! The occupied orbitals stay inside the occupied block
      !!
      !! This is the one that fails if the frozen set is orthonormalised in a
      !! single pass: the blocks come out mixed, and then `apply` holds part of
      !! an occupied orbital at the level shift.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), f(:, :), c(:, :), raw(:, :), basis(:, :)
      integer :: n_mo, n_basis, i

      call model_problem(s, f, c, n_mo, err)
      call raw_frozen(c, raw)
      call build_frozen_basis(raw, 2, s, basis, n_basis, err)

      do i = 1, 2
         call check(error, metric_residual(s, basis(:, :2), raw(:, i)) < TOL, &
                    "a frozen occupied orbital leaked into the virtual block")
         if (allocated(error)) return
      end do
   end subroutine test_basis_blocks

   subroutine test_basis_repeat(error)
      !! The same orbital frozen twice is a caller's bug, not something to absorb
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), f(:, :), c(:, :), raw(:, :), basis(:, :)
      integer :: n_mo, n_basis

      call model_problem(s, f, c, n_mo, err)
      allocate (raw(N, 2))
      raw(:, 1) = c(:, 1)
      raw(:, 2) = c(:, 1)
      call build_frozen_basis(raw, 2, s, basis, n_basis, err)
      call check(error, err%has_error(), "a repeated frozen orbital was accepted")
   end subroutine test_basis_repeat

   subroutine test_basis_end_to_end(error)
      !! Localization output in, constrained Fock out
      type(error_type), allocatable, intent(out) :: error
      type(fock_projector_t) :: proj
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), f(:, :), c(:, :), raw(:, :), basis(:, :), f_mo(:, :)
      integer :: n_mo, n_basis

      call model_problem(s, f, c, n_mo, err)
      call raw_frozen(c, raw)
      call build_frozen_basis(raw, 2, s, basis, n_basis, err)
      call proj%init(basis, s, 2, 3, SHIFT, err)
      call proj%apply(f, err)
      call check(error,.not. err%has_error(), "the projector rejected a basis built for it")
      if (allocated(error)) return

      f_mo = to_mo(basis, f)
      call check(error, maxval(abs(f_mo(3:, :2))) < ATOL, &
                 "the frozen occupied block is still coupled")
      if (allocated(error)) return
      call check(error, abs(f_mo(3, 3) - SHIFT) < ATOL, &
                 "the frozen virtual is not held at the shift")
      if (allocated(error)) return
      call check(error, maxval(abs(f_mo(4:, 3))) < ATOL, &
                 "the frozen virtual is still coupled to the unfrozen space")
   end subroutine test_basis_end_to_end

   subroutine test_bad_blocks(error)
      !! The blocks have to nest, and being told otherwise is a caller's bug
      type(error_type), allocatable, intent(out) :: error
      type(fock_projector_t) :: proj
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), f(:, :), c(:, :)
      integer :: n_mo

      call model_problem(s, f, c, n_mo, err)
      call proj%init(c, s, 4, 2, SHIFT, err)
      call check(error, err%has_error(), "accepted n_frozen_occ greater than n_frozen")
      if (allocated(error)) return

      err = error_t()
      call proj%init(c, s, 0, n_mo + 1, SHIFT, err)
      call check(error, err%has_error(), "accepted more frozen orbitals than there are")
   end subroutine test_bad_blocks

end module test_mqc_fock_projector

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fock_projector, only: collect_mqc_fock_projector
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_fock_projector", collect_mqc_fock_projector) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
