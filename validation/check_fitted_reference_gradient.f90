!! The derivative of a density-fitted two-electron operator, against finite difference
program check_fitted_reference_gradient
   !! **What is being checked, and why on its own.** A correlated gradient over
   !! a density-fitted reference needs
   !!
   !!     d/dR Tr(G[A] B),     G = J - k K/2 built from fitted integrals
   !!
   !! with both densities held fixed as matrices. With exact integrals that is a
   !! four-centre `int2e_ip1` contraction; fitted, it is three- and two-centre
   !! terms with intermediates of their own, and the exchange half of those is
   !! the part where a factor convention can go wrong without anything looking
   !! wrong.
   !!
   !! It is checked here rather than through the gradient it serves because by
   !! then it would sit behind amplitudes, a Lagrangian, a Z-vector solve and a
   !! relaxed density. This needs two matrices and nothing else.
   !!
   !! **Both routines are exercised, and that is the point.** The finite
   !! difference rebuilds the molecule, the auxiliary basis and the fitted
   !! tensor at each displaced geometry and evaluates `Tr(G[A] B)` through
   !! `fitted_potential_general`; the analytic side is
   !! `fitted_reference_gradient`. One is the derivative of the other, so a
   !! disagreement is in exactly one of the two.
   !!
   !! **Neither `A` nor `B` is a density.** They stand in for the reference
   !! density and the relaxed correlation density plus a Z-vector, and the
   !! second of those is symmetric, indefinite, and integrates to zero. A test
   !! that passed an SCF density twice would be blind to every place the code
   !! assumed idempotency -- which is the whole reason
   !! `fitted_potential_general` exists beside two builds that do assume it. So
   !! what goes in is a pair of fixed pseudo-random symmetric matrices,
   !! reproducible and resembling nothing that could be special-cased.
   !!
   !! The exchange fraction is varied too. At `k = 1` a term misplaced between
   !! the Coulomb and exchange halves can still cancel; at 0.5 it cannot.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, &
                                three_centre, two_centre, metric_inverse_sqrt
   use mqc_czt_cphf, only: fitted_potential_general
   use mqc_czt_gradient, only: fitted_reference_gradient
   use pic_blas_interfaces, only: pic_gemm
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none

   integer, parameter :: N_DIM = 3
   real(dp), parameter :: STEP = 1.0e-3_dp
   real(dp), parameter :: TOL = 1.0e-8_dp
      !! The integrals are analytic and there is no quadrature anywhere in this,
      !! so a central difference at 1e-3 is limited by its own `h^2` truncation
      !! and not by noise. Measured deviations sit two orders below this.
   integer :: n_bad

   n_bad = 0
   write (*, "(a)") "== the derivative of a fitted two-electron operator"

   call check_case("H2O / sto-3g, cc-pvdz-rifit, k=1", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", "cc-pvdz-rifit", 1.0_dp, n_bad)

   ! A hybrid's fraction: the Coulomb and exchange halves scale differently, so
   ! a term attributed to the wrong one shows up here and not above.
   call check_case("H2O / sto-3g, cc-pvdz-rifit, k=0.5", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", "cc-pvdz-rifit", 0.5_dp, n_bad)

   ! Asymmetric, and with p functions on two centres: the exchange term pairs
   ! `A` on one side of `Y^P` and `B` on the other, and swapping them survives a
   ! symmetric molecule.
   call check_case("HCN / sto-3g, cc-pvdz-rifit", [1, 6, 7], ["H", "C", "N"], &
                   reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                   "sto-3g", "cc-pvdz-rifit", 1.0_dp, n_bad)

   write (*, "(a)") ""
   if (n_bad == 0) then
      write (*, "(a)") "all fitted reference gradient checks passed"
   else
      write (*, "(a,i0,a)") "FAILED: ", n_bad, " case(s)"
      stop 1
   end if

contains

   subroutine check_case(label, numbers, symbols, coords, basis, aux_basis, &
                         k_scale, n_bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      real(dp), intent(in) :: k_scale
      integer, intent(inout) :: n_bad

      type(czt_molecule_t) :: mol, aux
      type(error_t) :: error
      real(dp), allocatable :: three(:, :), metric(:, :), jm12(:, :)
      real(dp), allocatable :: dm_a(:, :), dm_b(:, :)
      real(dp), allocatable :: analytic(:, :), numeric(:, :), shifted(:, :)
      real(dp) :: plus, minus, worst
      integer :: natm, ia, comp, nao

      natm = size(numbers)
      write (*, "(a)") ""
      write (*, "(a,a)") "== ", label
      flush (output_unit)

      call build_czt_molecule(numbers, symbols, coords, basis, mol, error)
      if (fail(error, n_bad)) return
      call build_czt_molecule(numbers, symbols, coords, aux_basis, aux, error)
      if (fail(error, n_bad)) return

      nao = mol%nao
      call pseudo_random_symmetric(nao, 0.5_dp, dm_a)
      call pseudo_random_symmetric(nao, 0.25_dp, dm_b)

      call three_centre(mol, aux, three)
      call two_centre(aux, metric)
      call metric_inverse_sqrt(metric, jm12, error)
      if (fail(error, n_bad)) return

      allocate (analytic(3, natm))
      analytic = 0.0_dp
      call fitted_reference_gradient(mol, aux, three, jm12, dm_a, dm_b, analytic, &
                                     error, k_scale=k_scale)
      call aux%destroy()
      call mol%destroy()

      allocate (numeric(3, natm), shifted(N_DIM, natm))
      do ia = 1, natm
         do comp = 1, 3
            shifted = coords
            shifted(comp, ia) = coords(comp, ia) + STEP
            call trace_at(numbers, symbols, shifted, basis, aux_basis, dm_a, dm_b, &
                          k_scale, plus, error)
            if (fail(error, n_bad)) return
            shifted(comp, ia) = coords(comp, ia) - STEP
            call trace_at(numbers, symbols, shifted, basis, aux_basis, dm_a, dm_b, &
                          k_scale, minus, error)
            if (fail(error, n_bad)) return
            numeric(comp, ia) = (plus - minus)/(2.0_dp*STEP)
         end do
      end do

      write (*, "(a)") "   atom      analytic (x, y, z)                        "// &
         "finite difference"
      do ia = 1, natm
         write (*, "(i6,3f14.9,a,3f14.9)") ia, analytic(:, ia), "   ", numeric(:, ia)
      end do
      worst = maxval(abs(analytic - numeric))
      write (*, "(a,es14.4)") "   largest deviation: ", worst
      flush (output_unit)

      if (worst > TOL) then
         write (*, "(a)") "  FAIL: analytic and finite difference disagree"
         n_bad = n_bad + 1
      end if

      ! Translational invariance. Worth little on its own -- it held at 1e-14
      ! while this same quantity was out by a factor of thirty -- but a term
      ! landing on the wrong atom fails it while the finite difference above
      ! still passes, so it costs nothing and covers a different mistake.
      if (maxval(abs(sum(analytic, dim=2))) > 1.0e-9_dp) then
         write (*, "(a,es14.4)") "  FAIL: does not sum to zero over atoms: ", &
            maxval(abs(sum(analytic, dim=2)))
         n_bad = n_bad + 1
      end if
   end subroutine check_case

   subroutine trace_at(numbers, symbols, coords, basis, aux_basis, dm_a, dm_b, &
                       k_scale, value, error)
      !! `Tr(G[A] B)` at one geometry, with both matrices carried in
      !!
      !! The molecule, the auxiliary basis and the fitted tensor are rebuilt;
      !! the densities are not. That is what makes this the *skeleton*
      !! derivative -- the one the analytic routine computes -- rather than the
      !! derivative of a converged quantity.
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      real(dp), intent(in) :: dm_a(:, :), dm_b(:, :)
      real(dp), intent(in) :: k_scale
      real(dp), intent(out) :: value
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: mol, aux
      real(dp), allocatable :: three(:, :), metric(:, :), jm12(:, :), b(:, :)
      real(dp), allocatable :: g(:, :)

      value = 0.0_dp
      call build_czt_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) return
      call build_czt_molecule(numbers, symbols, coords, aux_basis, aux, error)
      if (error%has_error()) return

      call three_centre(mol, aux, three)
      call two_centre(aux, metric)
      call metric_inverse_sqrt(metric, jm12, error)
      if (error%has_error()) return

      allocate (b(mol%nao*mol%nao, aux%nao))
      call pic_gemm(three, jm12, b)

      allocate (g(mol%nao, mol%nao))
      call fitted_potential_general(b, dm_a, g, k_scale=k_scale)
      value = sum(g*dm_b)

      call aux%destroy()
      call mol%destroy()
   end subroutine trace_at

   subroutine pseudo_random_symmetric(n, seed, p)
      !! A fixed symmetric matrix that is not a density
      !!
      !! Indefinite by construction, so nothing downstream can quietly rely on
      !! either argument being one: the real second argument is a reference
      !! density plus twice a relaxed correlation density, and the real first is
      !! contracted through a build that must not assume idempotency. Generated
      !! from a recurrence rather than `random_number` so a failure reproduces.
      integer, intent(in) :: n
      real(dp), intent(in) :: seed
      real(dp), allocatable, intent(out) :: p(:, :)

      integer :: i, j
      real(dp) :: x

      allocate (p(n, n))
      x = seed
      do j = 1, n
         do i = 1, j
            ! A cheap irrational-rotation sequence: equidistributed, and with no
            ! period a matrix of this size could line up with.
            x = mod(x + 0.6180339887498949_dp, 1.0_dp)
            p(i, j) = 0.1_dp*(x - 0.5_dp)
            p(j, i) = p(i, j)
         end do
      end do
   end subroutine pseudo_random_symmetric

   function fail(error, n_bad) result(bad)
      type(error_t), intent(inout) :: error
      integer, intent(inout) :: n_bad
      logical :: bad

      bad = error%has_error()
      if (bad) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
      end if
   end function fail

end program check_fitted_reference_gradient
