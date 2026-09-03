module test_mqc_mp2_hessian_ks_operator
   !! The identity the double-hybrid Hessian's rewrite rests on
   !!
   !! Every kernel contraction in `mqc_czt_mp2_hessian` has one of two
   !! shapes -- a matrix contracted against slots 1 and 3 of `L`, or against
   !! slots 2 and 4 -- and both are the reference operator applied to a
   !! symmetric generalized density:
   !!
   !!     sum_rs M(r,s) (L(r,p,s,q) + L(r,q,s,p))  =  2 V[M_sym]  in MO
   !!
   !! Two rather than four: the eight-fold symmetry of the integrals makes the
   !! two terms equal, so the pair is twice one application and not four times
   !! it. Writing `V = 2J - kK`, the same identity reads `4(J - kK/2)`, which is
   !! the same number said differently and the easier of the two to double by
   !! accident.
   !!
   !! That is what lets a double hybrid reach those sites by *applying* the
   !! Kohn-Sham operator rather than by materialising a four-index `f_xc`
   !! tensor, which would cost an `n_mo^4` array per perturbation. The whole
   !! shape of the remaining work depends on it, and until now it was checked
   !! only by a Python script in a scratch directory, against random integrals
   !! rather than against this code.
   !!
   !! **Run at three exchange fractions**, because the substitution has to hold
   !! for the hybrid weights a double hybrid uses and not only for the full
   !! exchange Hartree-Fock has. A rewrite that scaled exchange in one place and
   !! not the other passes at `k = 1` and fails at 0.53.
   !!
   !! **What this does not cover.** Only the two-electron half. Over a Kohn-Sham
   !! reference `L` also carries `f_xc`, and there is nothing to compare that
   !! against here: the whole point of the substitution is that no four-index
   !! `f_xc` tensor is ever built, so the kernel half is an apply on both sides
   !! of the identity and has no separate contraction to be checked against. It
   !! is pinned instead where it can be -- `xc_kernel_apply` against a density
   !! difference, one rung down in `test_mqc_xc_hessian`.
   !!
   !! The trial matrix is deliberately **not** symmetric. The contraction
   !! symmetrises one index pair itself, so its antisymmetric part must fall
   !! out entirely -- that is what makes an apply a legal substitution for
   !! sites whose generalized densities (`U^X`, `S^(X)`) are not symmetric.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available, &
                         xc_kernel_apply
   use mqc_czt_mp2_hessian, only: mp2_mo_eri_physicist
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_mp2_hessian_ks_operator_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape([ &
                                                0.0000000000_dp, 0.0000000000_dp, 0.2217589718_dp, &
                                                0.0000000000_dp, 1.4304281515_dp, -0.8870358873_dp, &
                                                0.0000000000_dp, -1.4304281515_dp, -0.8870358873_dp], [3, 3])

contains

   subroutine collect_mqc_mp2_hessian_ks_operator_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("hartree_fock_contraction_is_an_apply", test_hf_shape), &
                  new_unittest("hybrid_exchange_fractions_survive_the_rewrite", test_ks_shape), &
                  new_unittest("both_slot_pairings_agree", test_shapes_agree) &
                  ]
   end subroutine collect_mqc_mp2_hessian_ks_operator_tests

   subroutine test_hf_shape(error)
      !! `2J - K`, the operator the existing MP2 Hessian contracts against
      type(error_type), allocatable, intent(out) :: error

      call shape_against_apply("", 1.0_dp, error)
   end subroutine test_hf_shape

   subroutine test_ks_shape(error)
      !! The hybrid weights, where a half-applied scaling shows up
      !!
      !! 0.20 is B3LYP's and 0.53 is B2PLYP's. Neither is zero or one, which are
      !! the two values a `k_scale` dropped on one side of the identity would
      !! still reproduce.
      type(error_type), allocatable, intent(out) :: error

      call shape_against_apply("", 0.20_dp, error)
      if (allocated(error)) return
      call shape_against_apply("", 0.53_dp, error)
   end subroutine test_ks_shape

   subroutine test_shapes_agree(error)
      !! The two slot pairings are the same operator
      !!
      !! Sites in the Hessian contract slots 1 and 3 (`mp2_pair_rotation_augment`)
      !! or slots 2 and 4 (`mp2_perturbed_fock`). If those disagreed, one of the
      !! two families of rewrite would be wrong while the other passed.
      type(error_type), allocatable, intent(out) :: error

      call shape_against_apply("", 1.0_dp, error, compare_pairings=.true.)
   end subroutine test_shapes_agree

   subroutine shape_against_apply(functional, k_scale, error, compare_pairings)
      !! Contract `L` the way the Hessian does; compare against applying `V`
      character(len=*), intent(in) :: functional
      real(dp), intent(in) :: k_scale
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in), optional :: compare_pairings

      real(dp), parameter :: TOL = 1.0e-11_dp
         !! Two routes to one number through different arithmetic, so this is
         !! accumulation order and nothing else. There is no step and no grid
         !! difference: the kernel is evaluated once and used by both sides.
      type(czt_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: eri_mo(:, :, :, :), c(:, :), m(:, :), msym(:, :)
      real(dp), allocatable :: shape_a(:, :), shape_b(:, :), applied(:, :)
      real(dp), allocatable :: m_ao(:, :), v_ao(:, :), scratch(:, :)
      real(dp) :: acc, worst
      integer :: n, p, q, r, s, u, v
      logical :: pairings

      pairings = .false.
      if (present(compare_pairings)) pairings = compare_pairings

      ! Every setup failure below is a *test* failure, not a reason to stop
      ! quietly. Returning with `error` unallocated is what testdrive reads as a
      ! pass, so a molecule that would not build or an SCF that would not
      ! converge used to leave this suite green with no comparison made.
      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "the molecule did not build: "// &
                 err%get_message())
      if (allocated(error)) return
      if (len_trim(functional) > 0) then
         call xc_context_create(mol, functional, ctx, err, level=3)
         if (err%has_error()) call mol%destroy()
         call check(error,.not. err%has_error(), "the functional did not "// &
                    "resolve: "//err%get_message())
         if (allocated(error)) return
         call run_czt_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, xc=ctx)
      else
         call run_czt_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      end if
      if (err%has_error()) call mol%destroy()
      call check(error,.not. err%has_error(), "the reference did not converge: "// &
                 err%get_message())
      if (allocated(error)) return

      c = scf%orbitals
      n = size(c, 2)
      call mp2_mo_eri_physicist(mol, c, eri_mo, err)
      if (err%has_error()) call mol%destroy()
      call check(error,.not. err%has_error(), "the MO transform failed: "// &
                 err%get_message())
      if (allocated(error)) return

      ! Deliberately asymmetric: the contraction has to annihilate this part.
      allocate (m(n, n), msym(n, n))
      do q = 1, n
         do p = 1, n
            m(p, q) = 0.1_dp/(1.0_dp + real(abs(p - q), dp)) &
                      + 0.03_dp*real(p - q, dp)/real(n, dp)
         end do
      end do
      msym = 0.5_dp*(m + transpose(m))

      ! Shape A: contract slots 1 and 3, symmetrise the free pair.
      allocate (shape_a(n, n), shape_b(n, n))
      shape_a = 0.0_dp
      do q = 1, n
         do p = 1, n
            acc = 0.0_dp
            do s = 1, n
               do r = 1, n
                  acc = acc + m(r, s)*(l_element(eri_mo, k_scale, r, p, s, q) &
                                       + l_element(eri_mo, k_scale, r, q, s, p))
               end do
            end do
            shape_a(p, q) = acc
         end do
      end do

      ! Shape B: contract slots 2 and 4.
      shape_b = 0.0_dp
      do q = 1, n
         do p = 1, n
            acc = 0.0_dp
            do s = 1, n
               do r = 1, n
                  acc = acc + m(r, s)*(l_element(eri_mo, k_scale, p, r, q, s) &
                                       + l_element(eri_mo, k_scale, p, s, q, r))
               end do
            end do
            shape_b(p, q) = acc
         end do
      end do

      if (pairings) then
         call mol%destroy()
         worst = maxval(abs(shape_a - shape_b))
         call check(error, worst < TOL, "the two slot pairings should be one operator")
         return
      end if

      ! The apply: 4 V[M_sym], with V built the way the reference operator is.
      allocate (applied(n, n), m_ao(n, n), v_ao(n, n), scratch(n, n))
      applied = 0.0_dp
      do q = 1, n
         do p = 1, n
            acc = 0.0_dp
            do s = 1, n
               do r = 1, n
                  acc = acc + msym(r, s)*(2.0_dp*eri_mo(r, p, s, q) &
                                          - k_scale*eri_mo(r, p, q, s))
               end do
            end do
            applied(p, q) = 2.0_dp*acc
         end do
      end do

      call mol%destroy()

      worst = maxval(abs(shape_a - applied))
      call check(error, worst < TOL, &
                 "the Hessian's kernel contraction should equal 2 V[M_sym]")
   end subroutine shape_against_apply

   pure function l_element(eri_mo, k_scale, p, q, r, s) result(val)
      !! `L_pqrs = 2 <pq|rs> - k <pq|sr>`, the reference operator's weight
      real(dp), intent(in) :: eri_mo(:, :, :, :)
      real(dp), intent(in) :: k_scale
      integer, intent(in) :: p, q, r, s
      real(dp) :: val

      val = 2.0_dp*eri_mo(p, q, r, s) - k_scale*eri_mo(p, q, s, r)
   end function l_element

end module test_mqc_mp2_hessian_ks_operator

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mp2_hessian_ks_operator, only: collect_mqc_mp2_hessian_ks_operator_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [ &
                new_testsuite("mqc_mp2_hessian_ks_operator", &
                              collect_mqc_mp2_hessian_ks_operator_tests) &
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
