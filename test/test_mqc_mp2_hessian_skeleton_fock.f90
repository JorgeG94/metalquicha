module test_mqc_mp2_hessian_skeleton_fock
   !! `f^(X)` against a difference of the Fock matrix it is the derivative of
   !!
   !! `mp2_first_order_skeletons` says what it computes:
   !!
   !!     f^(X)_pq = C^T dh/dX C + sum_m [2 <pm|qm>^(X) - <pm|mq>^(X)]
   !!
   !! which is the nuclear derivative of `C^T (h + G[D]) C` with the coefficients
   !! and the density matrix both held fixed -- only the integrals move. So it
   !! can be differenced against the thing it differentiates: build the Fock
   !! matrix at displaced geometries from the *same* `C` and `D`, transform with
   !! the *same* `C`, and difference.
   !!
   !! **Why this exists before the code it guards.** A double-hybrid Hessian has
   !! to change this site: over a Kohn-Sham reference the exchange carries the
   !! functional's fraction and the exchange-correlation potential's derivative
   !! joins the sum. Nothing currently pins the term against anything external,
   !! so that edit would be unverifiable at the point it is made and would only
   !! surface much later, as a fraction-sized error in a Hessian that no symmetry
   !! check can reject. Pinning the Hartree-Fock behaviour first turns that edit
   !! into a change with a test either side of it.
   !!
   !! It is deliberately a check on the *skeleton* derivative alone. No orbital
   !! response enters here -- `C` never moves -- which is what makes the
   !! comparison exact rather than approximate, and is also why a fixed geometry
   !! step suffices.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_direct, only: build_fock_direct, schwarz_bounds, direct_stats_t
   use mqc_libcint_mp2_hessian, only: mp2_first_order_skeletons
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_mp2_hessian_skeleton_fock_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape([ &
                                                0.0000000000_dp, 0.0000000000_dp, 0.2217589718_dp, &
                                                0.0000000000_dp, 1.4304281515_dp, -0.8870358873_dp, &
                                                0.0000000000_dp, -1.4304281515_dp, -0.8870358873_dp], [3, 3])

contains

   subroutine collect_mqc_mp2_hessian_skeleton_fock_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("skeleton_fock_differences_the_fock", test_fx), &
                  new_unittest("skeleton_overlap_differences_the_overlap", test_sx) &
                  ]
   end subroutine collect_mqc_mp2_hessian_skeleton_fock_tests

   subroutine test_fx(error)
      !! Every entry of `f^(X)`, for every perturbation
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp), parameter :: TOL = 2.0e-8_dp
         !! Six times the measured worst entry, 3.38e-09 at this step, which is
         !! step error and nothing else: 8.54e-10 at 5e-5 and 1.35e-08 at 2e-4,
         !! ratios 3.95 and 4.00.
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: fx(:, :, :), sx(:, :, :), erix(:, :, :, :, :)
      real(dp), allocatable :: c(:, :), dens(:, :), coords(:, :)
      real(dp), allocatable :: fplus(:, :), fminus(:, :)
      real(dp) :: fd, worst
      integer :: n, n_occ, ia, comp, ix, p, q

      call reference(mol, scf, err)
      if (err%has_error()) return
      c = scf%orbitals
      dens = scf%density
      n = size(c, 1)
      n_occ = 5

      call mp2_first_order_skeletons(mol, c, n_occ, fx, sx, erix, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      allocate (coords(3, 3), fplus(n, n), fminus(n, n))
      worst = 0.0_dp
      ix = 0
      do ia = 1, 3
         do comp = 1, 3
            ix = ix + 1
            coords = WATER
            coords(comp, ia) = WATER(comp, ia) + H
            call mo_fock_at(coords, c, dens, fplus, err)
            if (err%has_error()) exit
            coords(comp, ia) = WATER(comp, ia) - H
            call mo_fock_at(coords, c, dens, fminus, err)
            if (err%has_error()) exit
            do q = 1, n
               do p = 1, n
                  fd = (fplus(p, q) - fminus(p, q))/(2.0_dp*H)
                  worst = max(worst, abs(fx(p, q, ix) - fd))
               end do
            end do
         end do
         if (err%has_error()) exit
      end do
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      call check(error, worst < TOL, &
                 "f^(X) should difference the Fock matrix it derives")
   end subroutine test_fx

   subroutine test_sx(error)
      !! `S^(X)` the same way, which costs nothing and rules out a shared fault
      !!
      !! If both `f^(X)` and `S^(X)` were transformed with the wrong
      !! coefficients, or indexed with the perturbation in the wrong slot, the
      !! Fock check alone could not tell that apart from a two-electron error.
      !! The overlap has no two-electron part at all, so it isolates the
      !! machinery they share.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp), parameter :: TOL = 1.2e-8_dp
         !! Six times the measured worst, 1.85e-09 at this step: 4.60e-10 at
         !! 5e-5 and 7.41e-09 at 2e-4, ratios 4.03 and 4.00.
         !!
         !! Not tighter than the Fock check, which is what a first guess assumed
         !! -- the overlap carries no two-electron accumulation, but the entries
         !! it is compared against are the same order, so the step error is too.
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: fx(:, :, :), sx(:, :, :), erix(:, :, :, :, :)
      real(dp), allocatable :: c(:, :), coords(:, :), splus(:, :), sminus(:, :)
      real(dp) :: fd, worst
      integer :: n, ia, comp, ix, p, q

      call reference(mol, scf, err)
      if (err%has_error()) return
      c = scf%orbitals
      n = size(c, 1)

      call mp2_first_order_skeletons(mol, c, 5, fx, sx, erix, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      allocate (coords(3, 3), splus(n, n), sminus(n, n))
      worst = 0.0_dp
      ix = 0
      do ia = 1, 3
         do comp = 1, 3
            ix = ix + 1
            coords = WATER
            coords(comp, ia) = WATER(comp, ia) + H
            call mo_overlap_at(coords, c, splus, err)
            if (err%has_error()) exit
            coords(comp, ia) = WATER(comp, ia) - H
            call mo_overlap_at(coords, c, sminus, err)
            if (err%has_error()) exit
            do q = 1, n
               do p = 1, n
                  fd = (splus(p, q) - sminus(p, q))/(2.0_dp*H)
                  worst = max(worst, abs(sx(p, q, ix) - fd))
               end do
            end do
         end do
         if (err%has_error()) exit
      end do
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      call check(error, worst < TOL, "S^(X) should difference the overlap")
   end subroutine test_sx

   subroutine reference(mol, scf, err)
      !! The converged reference, whose `C` and `D` every displacement reuses
      type(libcint_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
   end subroutine reference

   subroutine mo_fock_at(coords, c, dens, f_mo, err)
      !! `C^T (h + G[D]) C` at a displaced geometry, with `C` and `D` the
      !! reference's -- the quantity `f^(X)` is the derivative of
      real(dp), intent(in) :: coords(:, :), c(:, :), dens(:, :)
      real(dp), intent(out) :: f_mo(:, :)
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
      type(direct_stats_t) :: stats
      real(dp), allocatable :: h(:, :), bounds(:, :), f_ao(:, :)

      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "sto-3g", mol, err)
      if (err%has_error()) return
      allocate (f_ao(size(c, 1), size(c, 1)))
      call mol%core_hamiltonian(h)
      call schwarz_bounds(mol, bounds, err)
      if (.not. err%has_error()) then
         call build_fock_direct(mol, h, dens, bounds, f_ao, stats, err)
      end if
      call mol%destroy()
      if (err%has_error()) return
      f_mo = matmul(transpose(c), matmul(f_ao, c))
   end subroutine mo_fock_at

   subroutine mo_overlap_at(coords, c, s_mo, err)
      !! `C^T S C` at a displaced geometry, with the reference's `C`
      real(dp), intent(in) :: coords(:, :), c(:, :)
      real(dp), intent(out) :: s_mo(:, :)
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
      real(dp), allocatable :: s_ao(:, :)

      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "sto-3g", mol, err)
      if (err%has_error()) return
      call mol%overlap(s_ao)
      call mol%destroy()
      s_mo = matmul(transpose(c), matmul(s_ao, c))
   end subroutine mo_overlap_at

end module test_mqc_mp2_hessian_skeleton_fock

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mp2_hessian_skeleton_fock, only: collect_mqc_mp2_hessian_skeleton_fock_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [ &
                new_testsuite("mqc_mp2_hessian_skeleton_fock", &
                              collect_mqc_mp2_hessian_skeleton_fock_tests) &
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
