!! One- and two-particle density matrices of a CI wave function
module test_mqc_rdm
   !! The density matrices are the only thing the rest of MCSCF takes from the
   !! CI -- the orbital gradient, the generalised Fock and the energy are all
   !! contractions of them -- so an error here propagates into everything
   !! downstream while leaving the CI energy itself untouched. That is the
   !! reason for the shape of these tests: the CI energy is *not* evidence about
   !! them and cannot be used as such.
   !!
   !! Four kinds of check, deliberately overlapping, because the usual ways a
   !! two-particle density matrix is wrong are a transposed index pair and a
   !! factor of two, and each of them survives some of the tests below:
   !!
   !!   - traces, which catch normalisation but not ordering
   !!   - the `pq <-> rs` symmetry, which catches one transposition but is blind
   !!     to a uniform factor
   !!   - explicit elements against PySCF, which catch both but only where
   !!     they are looked at
   !!   - the energy rebuilt from the matrices against the CI eigenvalue, which
   !!     is the strongest of the four because it touches every element and
   !!     reaches the same number by completely different arithmetic
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_determinants, only: link_table_t, build_link_table
   use pic_lapack_interfaces, only: pic_syev
   use mqc_ci, only: absorb_one_electron, ci_hamiltonian
   use mqc_rdm, only: active_space_rdms, rdm_energy
   implicit none
   private

   public :: collect_mqc_rdm_tests

   integer, parameter :: NORB = 4
   integer, parameter :: NCHOL = 3

contains

   subroutine collect_mqc_rdm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("traces_count_electrons", test_traces), &
                  new_unittest("dm1_against_pyscf", test_dm1), &
                  new_unittest("dm2_against_pyscf", test_dm2), &
                  new_unittest("energy_rebuilt_from_rdms", test_energy), &
                  new_unittest("open_shell", test_open_shell), &
                  new_unittest("refusals", test_refusals) &
                  ]
   end subroutine collect_mqc_rdm_tests

   subroutine model_integrals(h1e, eri)
      real(dp), intent(out) :: h1e(NORB, NORB)
      real(dp), intent(out) :: eri(NORB, NORB, NORB, NORB)

      real(dp) :: b(NORB, NORB, NCHOL)
      integer :: p, q, r, s, l

      do q = 1, NORB
         do p = 1, NORB
            h1e(p, q) = -1.0_dp/real((p - 1) + (q - 1) + 2, dp)
         end do
      end do
      do l = 1, NCHOL
         do q = 1, NORB
            do p = 1, NORB
               b(p, q, l) = 1.0_dp/real((p - 1) + (q - 1) + (l - 1) + 3, dp)
            end do
         end do
      end do
      eri = 0.0_dp
      do s = 1, NORB
         do r = 1, NORB
            do q = 1, NORB
               do p = 1, NORB
                  do l = 1, NCHOL
                     eri(p, q, r, s) = eri(p, q, r, s) + b(p, q, l)*b(r, s, l)
                  end do
               end do
            end do
         end do
      end do
   end subroutine model_integrals

   subroutine ground_state(n_alpha, n_beta, h1e, eri, dm1, dm2, energy, err, ok)
      !! The exact ground state of the model, and its density matrices
      !!
      !! Diagonalised densely rather than with Davidson, deliberately. A density
      !! matrix is a property of a vector and has nothing to do with how the
      !! vector was found, so bringing the iterative solver in would couple
      !! these tests to it for no gain.
      !!
      !! It would also not work. This is the model from `test_mqc_ci.f90`, whose
      !! eigenvalues lie below every diagonal element and 3e-6 apart -- Davidson
      !! stalls on it, as `test_mqc_davidson.f90` explains at length. The model
      !! is kept because it is the one PySCF's reference density matrices were
      !! generated for, and because a dense solve does not care.
      integer, intent(in) :: n_alpha, n_beta
      real(dp), intent(out) :: h1e(NORB, NORB), eri(NORB, NORB, NORB, NORB)
      real(dp), allocatable, intent(out) :: dm1(:, :), dm2(:, :, :, :)
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      real(dp), allocatable :: folded(:, :), dense(:, :), values(:), vector(:, :)
      type(link_table_t) :: alpha, beta
      integer :: na, nb, info

      ok = .false.
      energy = 0.0_dp
      call model_integrals(h1e, eri)
      call absorb_one_electron(h1e, eri, n_alpha + n_beta, folded, err)
      call build_link_table(NORB, n_alpha, alpha, err)
      call build_link_table(NORB, n_beta, beta, err)
      if (err%has_error()) return
      call ci_hamiltonian(folded, alpha, beta, dense, err)
      if (err%has_error()) return

      na = alpha%n_strings
      nb = beta%n_strings
      allocate (values(na*nb))
      call pic_syev(dense, values, jobz="V", uplo="U", info=info)
      if (info /= 0) return
      energy = values(1)

      ! The sign of an eigenvector is arbitrary; the density matrices are
      ! quadratic in it, so it does not matter.
      allocate (vector(na, nb))
      vector = reshape(dense(:, 1), [na, nb])
      call active_space_rdms(vector, alpha, beta, dm1, dm2, err)
      ok = .not. err%has_error()
      call alpha%destroy()
      call beta%destroy()
   end subroutine ground_state

   subroutine test_traces(error)
      !! Both matrices count the electrons they describe
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: h1e(NORB, NORB), eri(NORB, NORB, NORB, NORB), energy
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      real(dp) :: trace, contracted
      logical :: ok
      integer :: p, q, r

      call ground_state(2, 2, h1e, eri, dm1, dm2, energy, err, ok)
      call check(error, ok, "the ground state should be found")
      if (allocated(error)) return

      trace = 0.0_dp
      do p = 1, NORB
         trace = trace + dm1(p, p)
      end do
      call check(error, trace, 4.0_dp, "the one-particle trace is the electron count", &
                 thr=1.0e-11_dp)
      if (allocated(error)) return

      call check(error, maxval(abs(dm1 - transpose(dm1))) < 1.0e-11_dp, &
                 "the one-particle matrix should be symmetric")
      if (allocated(error)) return

      ! sum_r d_pqrr = (N-1) D_pq, because sum_r E_rr is the number operator.
      ! This is the identity that pins the delta_qr correction: without it the
      ! contraction comes out N D_pq instead, which is a clean, plausible and
      ! completely wrong two-particle matrix.
      do q = 1, NORB
         do p = 1, NORB
            contracted = 0.0_dp
            do r = 1, NORB
               contracted = contracted + dm2(p, q, r, r)
            end do
            call check(error, abs(contracted - 3.0_dp*dm1(p, q)) < 1.0e-10_dp, &
                       "the partial trace of the two-particle matrix should be "// &
                       "(N-1) times the one-particle matrix")
            if (allocated(error)) return
         end do
      end do

      ! Exchanging the two electrons is a relabelling, not a different state.
      call check(error, maxval(abs(dm2 - reshape(dm2, [NORB, NORB, NORB, NORB], &
                                                 order=[3, 4, 1, 2]))) < 1.0e-11_dp, &
                 "the two-particle matrix should be symmetric under exchanging "// &
                 "the pq and rs pairs")
   end subroutine test_traces

   subroutine test_dm1(error)
      !! Every element of the one-particle matrix, against PySCF
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: h1e(NORB, NORB), eri(NORB, NORB, NORB, NORB), energy
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      logical :: ok
      integer :: p, q

      !> `fci.direct_spin1.make_rdm12` on the same model, row major.
      real(dp), parameter :: REF(16) = [ &
                             1.741672114860_dp, -0.175272798435_dp, -0.121591739916_dp, &
                             0.122535536906_dp, -0.175272798435_dp, 1.469591539197_dp, &
                             0.404123284177_dp, -0.534121281477_dp, -0.121591739916_dp, &
                             0.404123284177_dp, 0.240977509351_dp, 0.061678568088_dp, &
                             0.122535536906_dp, -0.534121281477_dp, 0.061678568088_dp, &
                             0.547758836591_dp]

      call ground_state(2, 2, h1e, eri, dm1, dm2, energy, err, ok)
      call check(error, ok, "the ground state should be found")
      if (allocated(error)) return

      do p = 1, NORB
         do q = 1, NORB
            call check(error, abs(dm1(p, q) - REF((p - 1)*NORB + q)) < 1.0e-9_dp, &
                       "one-particle element against PySCF")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_dm1

   subroutine test_dm2(error)
      !! Individual two-particle elements, against PySCF
      !!
      !! Three elements rather than all 256, chosen for what they distinguish:
      !! the fully diagonal one, one where the pairs are exchanged relative to
      !! each other, and one with four different indices, which is the only kind
      !! that can detect a transposition.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: h1e(NORB, NORB), eri(NORB, NORB, NORB, NORB), energy
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      logical :: ok

      call ground_state(2, 2, h1e, eri, dm1, dm2, energy, err, ok)
      call check(error, ok, "the ground state should be found")
      if (allocated(error)) return

      call check(error, dm2(1, 1, 1, 1), 1.622631009739_dp, "d(1,1,1,1)", thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, dm2(1, 2, 2, 1), -1.167262007598_dp, "d(1,2,2,1)", thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, dm2(2, 3, 4, 1), 0.128251555594_dp, "d(2,3,4,1)", thr=1.0e-9_dp)
   end subroutine test_dm2

   subroutine test_energy(error)
      !! The energy rebuilt from the density matrices is the CI eigenvalue
      !!
      !! The strongest check available on these matrices, and the one worth
      !! keeping if only one could be kept. Every element contributes; the two
      !! numbers are produced by entirely different arithmetic -- an iterative
      !! eigenvalue on one side, a contraction against integrals on the other --
      !! and a factor of two or a transposed index pair that leaves every trace
      !! identity above intact does not survive it.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: h1e(NORB, NORB), eri(NORB, NORB, NORB, NORB), energy
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      logical :: ok

      call ground_state(2, 2, h1e, eri, dm1, dm2, energy, err, ok)
      call check(error, ok, "the ground state should be found")
      if (allocated(error)) return
      call check(error, rdm_energy(h1e, eri, dm1, dm2), energy, &
                 "the energy from the density matrices should be the eigenvalue", &
                 thr=1.0e-10_dp)
      if (allocated(error)) return
      call check(error, energy, -1.058973188934_dp, &
                 "and both should be PySCF's", thr=1.0e-10_dp)
   end subroutine test_energy

   subroutine test_open_shell(error)
      !! Two alpha and one beta, where the spins are not interchangeable
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: h1e(NORB, NORB), eri(NORB, NORB, NORB, NORB), energy
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      real(dp) :: trace
      logical :: ok
      integer :: p

      call ground_state(2, 1, h1e, eri, dm1, dm2, energy, err, ok)
      call check(error, ok, "the ground state should be found")
      if (allocated(error)) return

      trace = 0.0_dp
      do p = 1, NORB
         trace = trace + dm1(p, p)
      end do
      call check(error, trace, 3.0_dp, "three electrons", thr=1.0e-11_dp)
      if (allocated(error)) return

      call check(error, dm1(1, 1), 1.632566590329_dp, "D(1,1) against PySCF", &
                 thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, dm1(2, 3), 0.314809213024_dp, "D(2,3) against PySCF", &
                 thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, dm2(2, 3, 4, 1), 0.042574561147_dp, "d(2,3,4,1)", thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, rdm_energy(h1e, eri, dm1, dm2), energy, &
                 "the energy should rebuild for an open shell too", thr=1.0e-10_dp)
   end subroutine test_open_shell

   subroutine test_refusals(error)
      !! Vectors that do not belong to the tables
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :), ci(:, :)

      call build_link_table(NORB, 2, alpha, err)
      call build_link_table(NORB, 2, beta, err)
      call check(error,.not. err%has_error(), "the tables should build")
      if (allocated(error)) return

      allocate (ci(alpha%n_strings + 1, beta%n_strings))
      ci = 0.0_dp
      call active_space_rdms(ci, alpha, beta, dm1, dm2, err)
      call check(error, err%has_error(), &
                 "a vector of the wrong shape should be refused")
      call err%clear()
      call alpha%destroy()
      call beta%destroy()
   end subroutine test_refusals

end module test_mqc_rdm

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_rdm, only: collect_mqc_rdm_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_rdm", collect_mqc_rdm_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
