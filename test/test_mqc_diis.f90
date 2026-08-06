!! Unit tests for DIIS subspace acceleration
module test_mqc_diis
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp, int64
   use mqc_diis, only: diis_state_t
   implicit none
   private

   public :: collect_mqc_diis

contains

   subroutine collect_mqc_diis(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("diis_no_extrapolation_below_two", test_below_two), &
                  new_unittest("diis_zero_error_reproduces_fock", test_zero_error), &
                  new_unittest("diis_coefficients_sum_to_one", test_coeffs_sum), &
                  new_unittest("diis_ring_evicts_oldest", test_ring_evicts), &
                  new_unittest("diis_matches_reference_algorithm", test_vs_reference), &
                  new_unittest("diis_cached_overlap_matches_direct", test_overlap_cache) &
                  ]
   end subroutine collect_mqc_diis

   subroutine test_below_two(error)
      !! A single stored vector cannot be extrapolated; fock must be untouched
      type(error_type), allocatable, intent(out) :: error
      type(diis_state_t) :: diis
      real(dp) :: fock(4), before(4)
      logical :: ok

      call diis%init(max_vectors=4, n_fock=4, n_error=4)
      fock = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      before = fock

      call diis%extrapolate(fock, ok)
      call check(error, .not. ok, "extrapolated with an empty subspace")
      if (allocated(error)) return
      call check(error, all(fock == before), "fock modified with an empty subspace")
      if (allocated(error)) return

      call diis%push(fock, [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])
      call diis%extrapolate(fock, ok)
      call check(error, .not. ok, "extrapolated with one stored vector")
      if (allocated(error)) return
      call check(error, all(fock == before), "fock modified with one stored vector")

      call diis%destroy()
   end subroutine test_below_two

   subroutine test_zero_error(error)
      !! With one exactly-converged entry among others, DIIS should land on it
      type(error_type), allocatable, intent(out) :: error
      type(diis_state_t) :: diis
      real(dp) :: fock(3), target_fock(3)
      logical :: ok

      call diis%init(max_vectors=4, n_fock=3, n_error=3)

      call diis%push([1.0_dp, 1.0_dp, 1.0_dp], [1.0_dp, 0.0_dp, 0.0_dp])
      call diis%push([2.0_dp, 2.0_dp, 2.0_dp], [0.0_dp, 1.0_dp, 0.0_dp])
      target_fock = [7.0_dp, 8.0_dp, 9.0_dp]
      call diis%push(target_fock, [1.0e-12_dp, 0.0_dp, 0.0_dp])

      fock = 0.0_dp
      call diis%extrapolate(fock, ok)
      call check(error, ok, "extrapolation failed")
      if (allocated(error)) return
      call check(error, maxval(abs(fock - target_fock)) < 1.0e-6_dp, &
                 "did not converge onto the near-zero-error entry")

      call diis%destroy()
   end subroutine test_zero_error

   subroutine test_coeffs_sum(error)
      !! The DIIS constraint is sum(c) = 1, so extrapolating a constant history
      !! must return that constant
      type(error_type), allocatable, intent(out) :: error
      type(diis_state_t) :: diis
      integer, parameter :: NF = 5, NPUSH = 4
      real(dp) :: fock(NF), e_in(NF, NPUSH)
      real(dp), parameter :: CONSTANT = 3.25_dp
      logical :: ok
      integer :: i

      ! Errors must be linearly independent or B is singular by construction --
      ! a fixed pattern like [i, 1/i, c, -i, c] spans only three dimensions no
      ! matter how many vectors you push.
      call fill_pseudorandom(e_in, 7)

      call diis%init(max_vectors=6, n_fock=NF, n_error=NF)
      do i = 1, NPUSH
         call diis%push(spread(CONSTANT, 1, NF), e_in(:, i))
      end do

      fock = 0.0_dp
      call diis%extrapolate(fock, ok)
      call check(error, ok, "extrapolation failed")
      if (allocated(error)) return
      call check(error, maxval(abs(fock - CONSTANT)) < 1.0e-10_dp, &
                 "coefficients do not sum to one")

      call diis%destroy()
   end subroutine test_coeffs_sum

   subroutine test_ring_evicts(error)
      !! Pushing past capacity must drop the oldest entry, not corrupt the ring
      !!
      !! Checked by equivalence rather than by bounding the result: DIIS
      !! coefficients are not convex, so an extrapolation may legitimately fall
      !! outside the range of the surviving entries. A history wrapped many
      !! times must instead be indistinguishable from one holding only the last
      !! max_vectors entries.
      type(error_type), allocatable, intent(out) :: error
      type(diis_state_t) :: wrapped, fresh
      integer, parameter :: NF = 6, NMAX = 3, NPUSH = 10
      real(dp) :: f_in(NF, NPUSH), e_in(NF, NPUSH)
      real(dp) :: fock_wrapped(NF), fock_fresh(NF)
      logical :: ok_wrapped, ok_fresh
      integer :: i

      call fill_pseudorandom(f_in, 17)
      call fill_pseudorandom(e_in, 23)

      call wrapped%init(max_vectors=NMAX, n_fock=NF, n_error=NF)
      do i = 1, NPUSH
         call wrapped%push(f_in(:, i), e_in(:, i))
      end do

      call fresh%init(max_vectors=NMAX, n_fock=NF, n_error=NF)
      do i = NPUSH - NMAX + 1, NPUSH
         call fresh%push(f_in(:, i), e_in(:, i))
      end do

      call check(error, wrapped%count() == NMAX, "subspace grew past max_vectors")
      if (allocated(error)) return

      fock_wrapped = 0.0_dp
      fock_fresh = 0.0_dp
      call wrapped%extrapolate(fock_wrapped, ok_wrapped)
      call fresh%extrapolate(fock_fresh, ok_fresh)

      call check(error, ok_wrapped .and. ok_fresh, "extrapolation failed")
      if (allocated(error)) return
      call check(error, all(fock_wrapped == fock_fresh), &
                 "wrapped history differs from the last max_vectors entries")

      call wrapped%destroy()
      call fresh%destroy()
   end subroutine test_ring_evicts

   subroutine test_vs_reference(error)
      !! Bit-for-bit against the shift-the-history, rebuild-B implementation the
      !! SCF driver used before the ring buffer and cached overlaps.
      type(error_type), allocatable, intent(out) :: error
      type(diis_state_t) :: diis
      integer, parameter :: NF = 12, NE = 12, NMAX = 4, NPUSH = 9
      real(dp) :: fock(NF), fock_ref(NF)
      real(dp) :: f_in(NF, NPUSH), e_in(NE, NPUSH)
      real(dp) :: hist_f(NF, NMAX), hist_e(NE, NMAX)
      integer :: n_stored, step, i
      logical :: ok, ok_ref

      call fill_pseudorandom(f_in, 11)
      call fill_pseudorandom(e_in, 29)

      call diis%init(max_vectors=NMAX, n_fock=NF, n_error=NE)
      n_stored = 0

      do step = 1, NPUSH
         call diis%push(f_in(:, step), e_in(:, step))

         ! Reference: shift down when full, append at the end
         if (n_stored < NMAX) then
            n_stored = n_stored + 1
         else
            hist_f(:, 1:NMAX - 1) = hist_f(:, 2:NMAX)
            hist_e(:, 1:NMAX - 1) = hist_e(:, 2:NMAX)
         end if
         hist_f(:, n_stored) = f_in(:, step)
         hist_e(:, n_stored) = e_in(:, step)

         fock = 0.0_dp
         call diis%extrapolate(fock, ok)
         call reference_extrapolate(hist_f, hist_e, n_stored, fock_ref, ok_ref)

         call check(error, ok .eqv. ok_ref, "solvability disagrees at step "//itoa(step))
         if (allocated(error)) return
         if (ok) then
            call check(error, all(fock == fock_ref), &
                       "extrapolated Fock differs from reference at step "//itoa(step))
            if (allocated(error)) return
         end if
      end do

      call diis%destroy()
   end subroutine test_vs_reference

   subroutine test_overlap_cache(error)
      !! The incrementally maintained overlap must equal a direct recomputation
      type(error_type), allocatable, intent(out) :: error
      type(diis_state_t) :: diis
      integer, parameter :: NE = 7, NMAX = 3, NPUSH = 8
      real(dp) :: e_in(NE, NPUSH), f_in(NE, NPUSH)
      real(dp) :: direct
      integer :: step, i, j, si, sj

      call fill_pseudorandom(e_in, 5)
      call fill_pseudorandom(f_in, 13)

      call diis%init(max_vectors=NMAX, n_fock=NE, n_error=NE)
      do step = 1, NPUSH
         call diis%push(f_in(:, step), e_in(:, step))
      end do

      do i = 1, diis%count()
         do j = 1, diis%count()
            si = modulo(diis%newest - diis%count() + i - 1, diis%max_vectors) + 1
            sj = modulo(diis%newest - diis%count() + j - 1, diis%max_vectors) + 1
            direct = sum(diis%error_history(:, si)*diis%error_history(:, sj))
            call check(error, diis%overlap(si, sj) == direct, &
                       "cached overlap differs from direct recomputation")
            if (allocated(error)) return
         end do
      end do

      call diis%destroy()
   end subroutine test_overlap_cache

   subroutine reference_extrapolate(hist_f, hist_e, n_stored, fock, ok)
      !! The pre-refactor algorithm, kept verbatim as the comparison target
      real(dp), intent(in) :: hist_f(:, :), hist_e(:, :)
      integer, intent(in) :: n_stored
      real(dp), intent(out) :: fock(:)
      logical, intent(out) :: ok

      real(dp), allocatable :: b_matrix(:, :), augmented(:, :), coefficients(:)
      integer :: n, i, j, pivot_row
      real(dp) :: pivot, factor

      ok = .false.
      fock = 0.0_dp
      if (n_stored < 2) return

      allocate (b_matrix(n_stored + 1, n_stored + 1))
      b_matrix = -1.0_dp
      b_matrix(n_stored + 1, n_stored + 1) = 0.0_dp
      do i = 1, n_stored
         do j = 1, n_stored
            b_matrix(i, j) = sum(hist_e(:, i)*hist_e(:, j))
         end do
      end do

      n = size(b_matrix, 1)
      allocate (augmented(n, n + 1), coefficients(n))
      augmented(:, 1:n) = b_matrix
      augmented(:, n + 1) = 0.0_dp
      augmented(n, n + 1) = -1.0_dp
      ok = .true.

      do i = 1, n
         pivot_row = i
         do j = i + 1, n
            if (abs(augmented(j, i)) > abs(augmented(pivot_row, i))) pivot_row = j
         end do
         if (pivot_row /= i) augmented([i, pivot_row], :) = augmented([pivot_row, i], :)
         pivot = augmented(i, i)
         if (abs(pivot) < 1.0e-14_dp) then
            ok = .false.
            return
         end if
         do j = i + 1, n
            factor = augmented(j, i)/pivot
            augmented(j, i:n + 1) = augmented(j, i:n + 1) - factor*augmented(i, i:n + 1)
         end do
      end do

      do i = n, 1, -1
         coefficients(i) = (augmented(i, n + 1) - &
                            sum(augmented(i, i + 1:n)*coefficients(i + 1:n)))/augmented(i, i)
      end do

      do i = 1, n_stored
         fock = fock + coefficients(i)*hist_f(:, i)
      end do
   end subroutine reference_extrapolate

   subroutine fill_pseudorandom(a, seed)
      !! Deterministic filler; a fixed stream keeps the reference comparison
      !! reproducible without depending on random_number's implementation
      real(dp), intent(out) :: a(:, :)
      integer, intent(in) :: seed
      integer :: i, j
      integer(int64) :: state
      integer(int64), parameter :: MULTIPLIER = 1103515245_int64
      integer(int64), parameter :: INCREMENT = 12345_int64
      integer(int64), parameter :: MODULUS = 2147483648_int64

      state = int(seed, int64)
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            state = modulo(MULTIPLIER*state + INCREMENT, MODULUS)
            a(i, j) = real(state, dp)/real(MODULUS, dp) - 0.5_dp
         end do
      end do
   end subroutine fill_pseudorandom

   pure function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=12) :: s
      write (s, "(i0)") i
   end function itoa

end module test_mqc_diis

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_diis, only: collect_mqc_diis
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_diis", collect_mqc_diis) &
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
