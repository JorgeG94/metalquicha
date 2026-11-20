!! Intrinsic atomic and bond orbitals, against PySCF and against themselves
module test_mqc_czt_iao
   !! The pinned numbers are `pyscf.lo.iao` and `pyscf.lo.ibo` on water in
   !! 6-31G at this geometry, fed this repository's basis and PySCF's own
   !! `minao` -- the data `basis_sets/minao/minao.json` was written from.
   !!
   !! Beside them, what needs no other program:
   !!
   !!   * the intrinsic orbitals are orthonormal under the overlap, there is
   !!     one per minimal-basis function, and the occupied orbitals lie inside
   !!     their span -- the property that makes the charges add up exactly;
   !!   * localising rotates the occupied orbitals among themselves, so the
   !!     density must come out unchanged and the functional must not fall;
   !!   * every bond orbital's atomic shares sum to one.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_iao, only: build_iao, iao_charges, ibo_localize, ibo_populations
   implicit none
   private

   public :: collect_mqc_czt_iao_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   ! PySCF, water/6-31G: IAO charges, and the localisation functional
   ! `sum_i sum_A (Q_A^i)^p` at the converged bond orbitals.
   real(dp), parameter :: Q_REF(3) = [-0.69940123_dp, 0.34953775_dp, 0.34986348_dp]
   real(dp), parameter :: L4_REF = 3.4371735678_dp
   real(dp), parameter :: L2_REF = 4.1222858200_dp
   ! The oxygen share of each O-H bond orbital, both exponents.
   real(dp), parameter :: BOND_ON_O(2) = [0.674932_dp, 0.674769_dp]

contains

   subroutine collect_mqc_czt_iao_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("iaos_are_orthonormal_and_span_the_occupied_space", &
                               test_span), &
                  new_unittest("iao_charges_match_pyscf", test_charges), &
                  new_unittest("ibos_match_pyscf_and_keep_the_density", test_ibo), &
                  new_unittest("exponent_two_is_pipek_mezey", test_exponent_two), &
                  new_unittest("an_unknown_exponent_is_refused", test_bad_exponent) &
                  ]
   end subroutine collect_mqc_czt_iao_tests

   subroutine water(mol, scf, s, err)
      !! Water in 6-31G, converged tightly, with its overlap
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      real(dp), allocatable, intent(out) :: s(:, :)
      type(error_t), intent(inout) :: err
      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
      if (err%has_error()) return
      call mol%overlap(s)
      call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err, &
                       in_core=.true., grad_tol=1.0e-9_dp)
   end subroutine water

   pure function functional(populations, p) result(l)
      !! `sum_i sum_A (Q_A^i)^p`
      real(dp), intent(in) :: populations(:, :)
      integer, intent(in) :: p
      real(dp) :: l

      l = sum(populations**p)
   end function functional

   subroutine test_span(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), iao(:, :), metric(:, :), projected(:, :)
      integer, allocatable :: owner(:)
      integer :: i

      call water(mol, scf, s, err)
      call check(error,.not. err%has_error() .and. scf%converged, "water must converge")
      if (allocated(error)) return
      call build_iao(mol, scf%orbitals(:, 1:5), iao, owner, err)
      call check(error,.not. err%has_error(), "IAOs: "//err%get_full_trace())
      if (allocated(error)) return

      ! O carries 1s 2s 2p, each H one 1s: seven.
      call check(error, size(iao, 2), 7)
      if (allocated(error)) return
      call check(error, count(owner == 1) == 5 .and. count(owner == 2) == 1 &
                 .and. count(owner == 3) == 1, "five IAOs on O and one on each H")
      if (allocated(error)) return

      metric = matmul(transpose(iao), matmul(s, iao))
      do i = 1, 7
         metric(i, i) = metric(i, i) - 1.0_dp
      end do
      call check(error, maxval(abs(metric)) < 1.0e-10_dp, &
                 "the IAOs are not orthonormal under the overlap")
      if (allocated(error)) return

      ! A A^T S is the projector onto their span; the occupied orbitals must
      ! come through it untouched.
      projected = matmul(iao, matmul(transpose(iao), matmul(s, scf%orbitals(:, 1:5))))
      call check(error, maxval(abs(projected - scf%orbitals(:, 1:5))) < 1.0e-9_dp, &
                 "the occupied orbitals do not lie inside the IAO span")
      call mol%destroy()
   end subroutine test_span

   subroutine test_charges(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), iao(:, :), q(:)
      integer, allocatable :: owner(:)

      call water(mol, scf, s, err)
      call check(error,.not. err%has_error() .and. scf%converged, "water must converge")
      if (allocated(error)) return
      call build_iao(mol, scf%orbitals(:, 1:5), iao, owner, err)
      call iao_charges(mol, scf%density, iao, owner, q, err)
      call check(error,.not. err%has_error(), "IAO charges: "//err%get_full_trace())
      if (allocated(error)) return
      write (*, "(a,3f14.8)") "    IAO charges ", q
      call check(error, abs(sum(q)) < 1.0e-9_dp, "the IAO charges must sum to zero")
      if (allocated(error)) return
      call check(error, maxval(abs(q - Q_REF)) < 1.0e-6_dp, &
                 "the IAO charges do not match PySCF's")
      call mol%destroy()
   end subroutine test_charges

   subroutine test_ibo(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), iao(:, :), loc(:, :), pops(:, :), before(:, :)
      real(dp), allocatable :: d_before(:, :), d_after(:, :)
      real(dp) :: on_o(2)
      integer, allocatable :: owner(:)
      integer :: i, sweeps, n_bond

      call water(mol, scf, s, err)
      call check(error,.not. err%has_error() .and. scf%converged, "water must converge")
      if (allocated(error)) return
      call build_iao(mol, scf%orbitals(:, 1:5), iao, owner, err)
      call ibo_populations(mol, scf%orbitals(:, 1:5), iao, owner, before, err)
      call ibo_localize(mol, scf%orbitals(:, 1:5), iao, owner, loc, err, sweeps=sweeps)
      call check(error,.not. err%has_error(), "IBOs: "//err%get_full_trace())
      if (allocated(error)) return
      call ibo_populations(mol, loc, iao, owner, pops, err)
      write (*, "(a,i0,a,f14.10,a,f14.10)") "    ", sweeps, " sweeps; functional ", &
         functional(before, 4), " -> ", functional(pops, 4)

      call check(error, functional(pops, 4) >= functional(before, 4) - 1.0e-12_dp, &
                 "localising must not lower the functional")
      if (allocated(error)) return
      call check(error, abs(functional(pops, 4) - L4_REF) < 1.0e-7_dp, &
                 "the localisation functional does not match PySCF's")
      if (allocated(error)) return
      call check(error, maxval(abs(sum(pops, dim=1) - 1.0_dp)) < 1.0e-9_dp, &
                 "each bond orbital's atomic shares must sum to one")
      if (allocated(error)) return

      ! Three orbitals on oxygen alone, two O-H bonds.
      call check(error, count(pops(1, :) > 0.999_dp), 3)
      if (allocated(error)) return
      n_bond = 0
      do i = 1, 5
         if (pops(1, i) > 0.999_dp) cycle
         n_bond = n_bond + 1
         if (n_bond <= 2) on_o(n_bond) = pops(1, i)
      end do
      call check(error, n_bond, 2)
      if (allocated(error)) return
      call check(error, abs(maxval(on_o) - BOND_ON_O(1)) < 1.0e-5_dp .and. &
                 abs(minval(on_o) - BOND_ON_O(2)) < 1.0e-5_dp, &
                 "the oxygen share of the O-H bond orbitals does not match PySCF's")
      if (allocated(error)) return

      ! A rotation within the occupied space: the density is the same matrix.
      d_before = matmul(scf%orbitals(:, 1:5), transpose(scf%orbitals(:, 1:5)))
      d_after = matmul(loc, transpose(loc))
      call check(error, maxval(abs(d_after - d_before)) < 1.0e-10_dp, &
                 "localising changed the density")
      call mol%destroy()
   end subroutine test_ibo

   subroutine test_exponent_two(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), iao(:, :), loc(:, :), pops(:, :)
      integer, allocatable :: owner(:)

      call water(mol, scf, s, err)
      call check(error,.not. err%has_error() .and. scf%converged, "water must converge")
      if (allocated(error)) return
      call build_iao(mol, scf%orbitals(:, 1:5), iao, owner, err)
      call ibo_localize(mol, scf%orbitals(:, 1:5), iao, owner, loc, err, exponent=2)
      call ibo_populations(mol, loc, iao, owner, pops, err)
      call check(error,.not. err%has_error(), "IBOs: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, abs(functional(pops, 2) - L2_REF) < 1.0e-7_dp, &
                 "the exponent-2 functional does not match PySCF's")
      call mol%destroy()
   end subroutine test_exponent_two

   subroutine test_bad_exponent(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), iao(:, :), loc(:, :)
      integer, allocatable :: owner(:)

      call water(mol, scf, s, err)
      call build_iao(mol, scf%orbitals(:, 1:5), iao, owner, err)
      call check(error,.not. err%has_error(), "setup: "//err%get_full_trace())
      if (allocated(error)) return
      call ibo_localize(mol, scf%orbitals(:, 1:5), iao, owner, loc, err, exponent=3)
      call check(error, err%has_error(), "an exponent of three must be refused")
      call mol%destroy()
   end subroutine test_bad_exponent

end module test_mqc_czt_iao

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_iao, only: collect_mqc_czt_iao_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_iao", collect_mqc_czt_iao_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
