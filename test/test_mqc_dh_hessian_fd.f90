!! The double hybrid's perturbative Hessian, against its own gradient
module test_mqc_dh_hessian_fd
   !! Phase 5's gate: the term the analytic double-hybrid Hessian adds
   !!
   !! `mp2_correlation_hessian` handed an `xc` context returns the second
   !! derivative of a double hybrid's perturbative term over a Kohn-Sham
   !! reference, and `libcint_mp2_gradient` handed the same context returns the
   !! first derivative of that same term. One differences the other. Nothing in
   !! the ladder below this point tests the two together: the kernel identities
   !! pin single contractions and the skeleton checks pin single matrices, while
   !! this is the whole assembly -- response, Z-vector, relaxed density and all
   !! -- against a derivative with an independent reference behind it.
   !!
   !! **The grid is held fixed on both sides, and that is the whole design of
   !! this file.** `xc_context_create` builds an atom-centred grid, so rebuilding
   !! a context at each displaced geometry would move the grid with the nuclei
   !! and the finite difference would carry a grid-response term the analytic
   !! Hessian deliberately omits -- the same term, and the same trap, that
   !! `test_mqc_xc_hessian`'s `xc_at` documents one rung down. One context is
   !! built at the reference geometry and every displaced evaluation is handed
   !! it; only `mol` moves. That is also why `libcint_mp2_gradient` gained
   !! `fixed_grid`: without it the comparison lands about 1e-5 apart, which is
   !! too small to look like a bug and too large to be step noise, and reads as
   !! a grid-convergence problem rather than as two different derivatives.
   !!
   !! **The density is not held fixed**, unlike `xc_hessian`'s test: this is the
   !! full second derivative, so each displaced point re-converges its SCF. The
   !! orbital response is most of what is being checked.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_libcint_mp2_gradient, only: libcint_mp2_gradient
   use mqc_libcint_mp2_hessian, only: mp2_correlation_hessian
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   implicit none
   private

   public :: collect_mqc_dh_hessian_fd_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.000000_dp, 0.000000_dp, &
                           0.0_dp, 0.000000_dp, 1.814137_dp, &
                           0.0_dp, 1.756000_dp, -0.454300_dp], [3, 3])
   integer, parameter :: WATER_NELEC = 10

   !> B2PLYP: a GGA-based double hybrid, so the kernel ladder this exercises is
   !> the one that exists. A meta-GGA one would be refused by
   !> `xc_kernel2_apply` rather than silently approximated, and the bridge keeps
   !> it off the analytic path for that reason.
   character(len=*), parameter :: FUNCTIONAL = "b2plyp"

   !> Atom 1 (O), z -- the column `test_mqc_mp2_hessian_fd` perturbs, so a
   !> reader comparing the two files is looking at the same entry of the same
   !> molecule with only the reference changed.
   integer, parameter :: PERT_ATOM = 1
   integer, parameter :: PERT_CART = 3

   !> Bohr; the 7-point O(h^6) stencil's step, as one rung down.
   real(dp), parameter :: STEP = 2.0e-3_dp

   !> Twenty times the measured worst entry, 4.891e-10 with this step on this
   !> grid. Loose against the stencil's own reach and not against the
   !> derivative's: what the number absorbs is the two SCF convergences and the
   !> stencil's amplification of gradient-level noise, and a real disagreement
   !> between these two derivatives would sit two orders above it.
   real(dp), parameter :: TOL = 1.0e-8_dp

contains

   subroutine collect_mqc_dh_hessian_fd_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("perturbative_hessian_differences_its_own_gradient", &
                               column_against_fd), &
                  new_unittest("perturbative_hessian_column_sums_to_nothing", &
                               column_translates) &
                  ]
   end subroutine collect_mqc_dh_hessian_fd_tests

   !> The perturbative term's gradient at `coords`, on the grid in `ctx`
   subroutine pt2_gradient_at(ctx, coords, gradient, err, ok)
      type(xc_context_t), intent(inout) :: ctx
      real(dp), intent(in) :: coords(3, 3)
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "6-31g", mol, err)
      if (err%has_error()) return

      ! Tighter than the default for the reason one rung down gives: a loose
      ! SCF is the usual cause of a noisy finite-difference column and the
      ! cheapest thing to rule out before doubting the derivative.
      call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, &
                           .false., scf, err, xc=ctx)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if

      call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, &
                                WATER_NELEC/2, gradient, err, &
                                xc=ctx, scf_density=scf%density, &
                                pt2_scale=ctx%pt2_fraction, fixed_grid=.true.)
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine pt2_gradient_at

   !> `H[:, X] = d(perturbative gradient)/dX` by the 7-point O(h^6) stencil
   subroutine fd_column(ctx, column, err, ok)
      type(xc_context_t), intent(inout) :: ctx
      real(dp), intent(out) :: column(9)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      real(dp), parameter :: WEIGHT(6) = [-1.0_dp, 9.0_dp, -45.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]
      integer, parameter :: OFFSET(6) = [-3, -2, -1, 1, 2, 3]
      real(dp), allocatable :: gradient(:, :)
      real(dp) :: coords(3, 3), accumulated(3, 3)
      integer :: point

      ok = .false.
      accumulated = 0.0_dp
      do point = 1, 6
         coords = WATER
         coords(PERT_CART, PERT_ATOM) = coords(PERT_CART, PERT_ATOM) &
                                        + real(OFFSET(point), dp)*STEP
         call pt2_gradient_at(ctx, coords, gradient, err, ok)
         if (.not. ok) return
         accumulated = accumulated + WEIGHT(point)*gradient
      end do
      column = reshape(accumulated/(60.0_dp*STEP), [9])
      ok = .true.
   end subroutine fd_column

   !> The analytic column, from the same context and reference geometry
   subroutine analytic_column(ctx, column, err, ok)
      type(xc_context_t), intent(inout) :: ctx
      real(dp), intent(out) :: column(9)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: hess_corr(:, :, :, :), hess_ref(:, :, :, :)
      integer :: ia, a

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "6-31g", mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, &
                           .false., scf, err, xc=ctx)
      if (.not. err%has_error()) then
         ! `hess_ref` is Hartree-Fock-shaped and discarded by contract -- the
         ! bridge takes its whole reference block from `ks_hessian` instead.
         call mp2_correlation_hessian(mol, scf%orbitals, scf%orbital_energies, &
                                      scf%density, WATER_NELEC/2, 0, &
                                      hess_corr, hess_ref, err, &
                                      xc=ctx, pt2_scale=ctx%pt2_fraction)
      end if
      call mol%destroy()
      if (err%has_error()) return

      do ia = 1, 3
         do a = 1, 3
            column(3*(ia - 1) + a) = hess_corr(a, PERT_CART, ia, PERT_ATOM)
         end do
      end do
      ok = .true.
   end subroutine analytic_column

   subroutine column_against_fd(error)
      !! The gate this phase turns on
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(xc_context_t) :: ctx
      real(dp) :: fd(9), analytic(9)
      logical :: ok
      integer :: threads

      if (.not. xc_available()) return

      ! One thread, restored afterwards: the stencil divides by 60h = 0.12, so
      ! it multiplies gradient-level noise by about eight, and the unordered
      ! OpenMP critical merges do not sum in a fixed order. A gate this tight
      ! has to be handed a deterministic number.
      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call reference_context(ctx, err, ok)
      if (ok) call fd_column(ctx, fd, err, ok)
      if (ok) call analytic_column(ctx, analytic, err, ok)
      call omp_set_num_threads(threads)

      call check(error, ok, "the double-hybrid column did not evaluate: "// &
                 err%get_message())
      if (allocated(error)) return

      write (*, "(a, es10.3)") "        max |analytic - fd| = ", &
         maxval(abs(analytic - fd))
      call check(error, maxval(abs(analytic - fd)) < TOL, &
                 "the perturbative term's analytic Hessian does not difference "// &
                 "its own gradient")
   end subroutine column_against_fd

   subroutine column_translates(error)
      !! Translating the molecule cannot change its energy -- nearly
      !!
      !! Weak on its own: this repository has a Hessian wrong by 14
      !! Hartree/bohr that passed exactly this check. It is here for the one
      !! thing it catches cheaply, a term deposited on the wrong atom, and the
      !! tolerance is what makes it worth keeping rather than what makes it
      !! pass.
      !!
      !! **This column does not sum to zero, and should not.** The grid is
      !! atom-centred and held fixed, so translating the nuclei through it
      !! genuinely changes the quadrature -- the omitted grid response is
      !! exactly the residual. Measured: 2.219e-07 on the level-3 grid this
      !! case uses and 1.282e-08 on a level-5 one, a factor of seventeen for
      !! one refinement. That scaling is the evidence for calling it
      !! quadrature: a term deposited on the wrong atom is a property of the
      !! assembly and would not move with the grid at all, and it would be
      !! four orders larger.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(xc_context_t) :: ctx
      real(dp) :: analytic(9), net(3)
      logical :: ok
      integer :: a

      if (.not. xc_available()) return

      call reference_context(ctx, err, ok)
      if (ok) call analytic_column(ctx, analytic, err, ok)
      call check(error, ok, "the double-hybrid column did not evaluate: "// &
                 err%get_message())
      if (allocated(error)) return

      do a = 1, 3
         net(a) = analytic(a) + analytic(3 + a) + analytic(6 + a)
      end do
      write (*, "(a, es10.3)") "        max |net| = ", maxval(abs(net))
      call check(error, maxval(abs(net)) < 2.0e-6_dp, &
                 "the perturbative Hessian column does not sum to the grid "// &
                 "response alone")
   end subroutine column_translates

   !> The grid every evaluation in this file shares
   subroutine reference_context(ctx, err, ok)
      type(xc_context_t), intent(out) :: ctx
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "6-31g", mol, err)
      if (err%has_error()) return
      call xc_context_create(mol, FUNCTIONAL, ctx, err, level=3)
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine reference_context

end module test_mqc_dh_hessian_fd

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_dh_hessian_fd, only: collect_mqc_dh_hessian_fd_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [ &
                new_testsuite("mqc_dh_hessian_fd", collect_mqc_dh_hessian_fd_tests) &
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
