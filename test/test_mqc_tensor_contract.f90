module test_mqc_tensor_contract
   !! `contract` against the loops its index strings stand for.
   !!
   !! Every case fills its operands with numbers that have no symmetry at all
   !! and uses a different extent for every label, so a transposed index or a
   !! permutation applied backwards changes the answer and cannot hide behind
   !! two indices that happen to run over the same range.
   !!
   !! The strings are taken from the Lambda equations that are the reason the
   !! routine exists: a product that keeps the first operand's order, one that
   !! reorders the result, one that sums two labels, one that sums none, and
   !! the mixed-rank forms.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_tensor_contract, only: contract
   implicit none
   private
   public :: collect_mqc_tensor_contract_tests

   integer, parameter :: NI = 2, NJ = 3, NK = 4, NL = 5
   integer, parameter :: NA = 3, NB = 4, NC = 2, ND = 5
   real(dp), parameter :: TOL = 1.0e-12_dp

contains

   subroutine collect_mqc_tensor_contract_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("one_label_summed_order_kept", test_plain), &
                  new_unittest("result_reordered", test_reordered), &
                  new_unittest("two_labels_summed", test_two_summed), &
                  new_unittest("outer_product", test_outer), &
                  new_unittest("rank_four_to_rank_two", test_to_matrix), &
                  new_unittest("matrix_with_rank_four", test_matrix_first), &
                  new_unittest("adds_and_scales", test_accumulates), &
                  new_unittest("large_enough_to_be_threaded", test_large) &
                  ]
   end subroutine collect_mqc_tensor_contract_tests

   subroutine fill4(x, seed)
      !! Numbers with no symmetry between any two indices
      real(dp), intent(out) :: x(:, :, :, :)
      real(dp), intent(in) :: seed

      integer :: p, q, r, s

      do s = 1, size(x, 4)
         do r = 1, size(x, 3)
            do q = 1, size(x, 2)
               do p = 1, size(x, 1)
                  x(p, q, r, s) = sin(seed + 1.3_dp*real(p, dp) + 0.37_dp*real(q, dp)**2 &
                                      + 2.1_dp*real(r, dp) - 0.9_dp*real(s*p, dp))
               end do
            end do
         end do
      end do
   end subroutine fill4

   subroutine fill2(x, seed)
      real(dp), intent(out) :: x(:, :)
      real(dp), intent(in) :: seed

      integer :: p, q

      do q = 1, size(x, 2)
         do p = 1, size(x, 1)
            x(p, q) = cos(seed + 0.7_dp*real(p, dp) - 1.9_dp*real(q, dp) + 0.3_dp*real(p*q, dp))
         end do
      end do
   end subroutine fill2

   subroutine test_plain(error)
      !! 'ijac,cb->ijab'
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: l2(NI, NJ, NA, NC), v(NC, NB), got(NI, NJ, NA, NB), want(NI, NJ, NA, NB)
      integer :: i, j, a, b, c

      call fill4(l2, 0.1_dp)
      call fill2(v, 0.2_dp)
      want = 0.0_dp
      do b = 1, NB
         do c = 1, NC
            do a = 1, NA
               do j = 1, NJ
                  do i = 1, NI
                     want(i, j, a, b) = want(i, j, a, b) + l2(i, j, a, c)*v(c, b)
                  end do
               end do
            end do
         end do
      end do
      got = 0.0_dp
      call contract("ijac,cb->ijab", l2, v, got, 1.0_dp)
      call check(error, maxval(abs(got - want)) < TOL, "'ijac,cb->ijab' is wrong")
   end subroutine test_plain

   subroutine test_reordered(error)
      !! 'kalb,ijkl->ijab': the kept labels of the first operand end up last
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: ovov(NK, NA, NL, NB), x(NI, NJ, NK, NL)
      real(dp) :: got(NI, NJ, NA, NB), want(NI, NJ, NA, NB)
      integer :: i, j, k, l, a, b

      call fill4(ovov, 0.3_dp)
      call fill4(x, 0.4_dp)
      want = 0.0_dp
      do b = 1, NB
         do l = 1, NL
            do a = 1, NA
               do k = 1, NK
                  do j = 1, NJ
                     do i = 1, NI
                        want(i, j, a, b) = want(i, j, a, b) + ovov(k, a, l, b)*x(i, j, k, l)
                     end do
                  end do
               end do
            end do
         end do
      end do
      got = 0.0_dp
      call contract("kalb,ijkl->ijab", ovov, x, got, 1.0_dp)
      call check(error, maxval(abs(got - want)) < TOL, "'kalb,ijkl->ijab' is wrong")
   end subroutine test_reordered

   subroutine test_two_summed(error)
      !! 'ijcd,klcd->ijkl'
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: l2(NI, NJ, NC, ND), tau(NK, NL, NC, ND)
      real(dp) :: got(NI, NJ, NK, NL), want(NI, NJ, NK, NL)
      integer :: i, j, k, l, c, d

      call fill4(l2, 0.5_dp)
      call fill4(tau, 0.6_dp)
      want = 0.0_dp
      do d = 1, ND
         do c = 1, NC
            do l = 1, NL
               do k = 1, NK
                  do j = 1, NJ
                     do i = 1, NI
                        want(i, j, k, l) = want(i, j, k, l) + l2(i, j, c, d)*tau(k, l, c, d)
                     end do
                  end do
               end do
            end do
         end do
      end do
      got = 0.0_dp
      call contract("ijcd,klcd->ijkl", l2, tau, got, 1.0_dp)
      call check(error, maxval(abs(got - want)) < TOL, "'ijcd,klcd->ijkl' is wrong")
   end subroutine test_two_summed

   subroutine test_outer(error)
      !! 'ia,jb->ijab', with nothing summed
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: l1(NI, NA), v(NJ, NB), got(NI, NJ, NA, NB), want(NI, NJ, NA, NB)
      integer :: i, j, a, b

      call fill2(l1, 0.7_dp)
      call fill2(v, 0.8_dp)
      do b = 1, NB
         do a = 1, NA
            do j = 1, NJ
               do i = 1, NI
                  want(i, j, a, b) = l1(i, a)*v(j, b)
               end do
            end do
         end do
      end do
      got = 0.0_dp
      call contract("ia,jb->ijab", l1, v, got, 1.0_dp)
      call check(error, maxval(abs(got - want)) < TOL, "'ia,jb->ijab' is wrong")
   end subroutine test_outer

   subroutine test_to_matrix(error)
      !! 'ijbc,bacj->ia' and 'jiba,jb->ia'
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: l2(NI, NJ, NB, NC), w(NB, NA, NC, NJ), got(NI, NA), want(NI, NA)
      real(dp) :: m3(NJ, NI, NB, NA), t1(NJ, NB)
      integer :: i, j, a, b, c

      call fill4(l2, 0.9_dp)
      call fill4(w, 1.0_dp)
      want = 0.0_dp
      do c = 1, NC
         do b = 1, NB
            do j = 1, NJ
               do a = 1, NA
                  do i = 1, NI
                     want(i, a) = want(i, a) + l2(i, j, b, c)*w(b, a, c, j)
                  end do
               end do
            end do
         end do
      end do
      got = 0.0_dp
      call contract("ijbc,bacj->ia", l2, w, got, 1.0_dp)
      call check(error, maxval(abs(got - want)) < TOL, "'ijbc,bacj->ia' is wrong")
      if (allocated(error)) return

      call fill4(m3, 1.1_dp)
      call fill2(t1, 1.2_dp)
      want = 0.0_dp
      do b = 1, NB
         do j = 1, NJ
            do a = 1, NA
               do i = 1, NI
                  want(i, a) = want(i, a) + m3(j, i, b, a)*t1(j, b)
               end do
            end do
         end do
      end do
      got = 0.0_dp
      call contract("jiba,jb->ia", m3, t1, got, 1.0_dp)
      call check(error, maxval(abs(got - want)) < TOL, "'jiba,jb->ia' is wrong")
   end subroutine test_to_matrix

   subroutine test_matrix_first(error)
      !! 'ca,icjb->ijab' and 'kc,kjcb->jb'
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: m(NC, NA), ovov(NI, NC, NJ, NB), got(NI, NJ, NA, NB), want(NI, NJ, NA, NB)
      real(dp) :: l1(NK, NC), theta(NK, NJ, NC, NB), got2(NJ, NB), want2(NJ, NB)
      integer :: i, j, k, a, b, c

      call fill2(m, 1.3_dp)
      call fill4(ovov, 1.4_dp)
      want = 0.0_dp
      do b = 1, NB
         do j = 1, NJ
            do c = 1, NC
               do a = 1, NA
                  do i = 1, NI
                     want(i, j, a, b) = want(i, j, a, b) + m(c, a)*ovov(i, c, j, b)
                  end do
               end do
            end do
         end do
      end do
      got = 0.0_dp
      call contract("ca,icjb->ijab", m, ovov, got, 1.0_dp)
      call check(error, maxval(abs(got - want)) < TOL, "'ca,icjb->ijab' is wrong")
      if (allocated(error)) return

      call fill2(l1, 1.5_dp)
      call fill4(theta, 1.6_dp)
      want2 = 0.0_dp
      do b = 1, NB
         do c = 1, NC
            do j = 1, NJ
               do k = 1, NK
                  want2(j, b) = want2(j, b) + l1(k, c)*theta(k, j, c, b)
               end do
            end do
         end do
      end do
      got2 = 0.0_dp
      call contract("kc,kjcb->jb", l1, theta, got2, 1.0_dp)
      call check(error, maxval(abs(got2 - want2)) < TOL, "'kc,kjcb->jb' is wrong")
   end subroutine test_matrix_first

   subroutine test_accumulates(error)
      !! The result is added to, and scaled, rather than overwritten
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: a(NI, NK), b(NK, NJ), got(NI, NJ), once(NI, NJ)

      call fill2(a, 1.7_dp)
      call fill2(b, 1.8_dp)
      once = matmul(a, b)
      got = 1.0_dp
      call contract("ik,kj->ij", a, b, got, -0.5_dp)
      call contract("ik,kj->ij", a, b, got, 2.0_dp)
      call check(error, maxval(abs(got - (1.0_dp + 1.5_dp*once))) < TOL, &
                 "two calls did not add up to one and a half products on top of what was there")
   end subroutine test_accumulates

   subroutine test_large(error)
      !! 'kalb,ijkl->ijab' again, on operands big enough to cross the sizes at
      !! which the permutations and the product are split over threads
      !!
      !! The small cases all run on one thread, whatever the machine, so they
      !! say nothing about the split: a slab offset wrong by one, or a block of
      !! columns counted twice, shows only here.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: MI = 12, MJ = 13, MK = 37, ML = 41, MA = 39, MB = 43
      real(dp), allocatable :: ovov(:, :, :, :), x(:, :, :, :)
      real(dp), allocatable :: got(:, :, :, :), want(:, :, :, :)
      integer :: i, j, k, l, a, b

      allocate (ovov(MK, MA, ML, MB), x(MI, MJ, MK, ML))
      allocate (got(MI, MJ, MA, MB), want(MI, MJ, MA, MB))
      call fill4(ovov, 2.1_dp)
      call fill4(x, 2.2_dp)
      want = 0.0_dp
      do b = 1, MB
         do l = 1, ML
            do a = 1, MA
               do k = 1, MK
                  do j = 1, MJ
                     do i = 1, MI
                        want(i, j, a, b) = want(i, j, a, b) + ovov(k, a, l, b)*x(i, j, k, l)
                     end do
                  end do
               end do
            end do
         end do
      end do
      got = 0.0_dp
      call contract("kalb,ijkl->ijab", ovov, x, got, 1.0_dp)
      call check(error, maxval(abs(got - want)) < 1.0e-10_dp, &
                 "a contraction large enough to thread is wrong")
   end subroutine test_large

end module test_mqc_tensor_contract

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_tensor_contract, only: collect_mqc_tensor_contract_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_tensor_contract", collect_mqc_tensor_contract_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
