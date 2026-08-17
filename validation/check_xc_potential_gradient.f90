!! The derivative of the exchange-correlation potential, against finite difference
program check_xc_potential_gradient
   !! **What is being checked, and why on its own.** A double hybrid gradient
   !! eliminates its orbital response through a Z-vector, and what the Z-vector
   !! contracts with is the derivative of the *reference* operator. For a
   !! Kohn-Sham reference that operator contains `V_xc`, so the assembly needs
   !!
   !!     d/dR Tr(P V_xc[D])
   !!
   !! with both densities held fixed as matrices. That is a term of the same
   !! order as the Coulomb one beside it, and leaving it out gives a gradient
   !! with nothing obviously wrong about it.
   !!
   !! It is checked here rather than through the double hybrid because by then
   !! it would sit behind amplitudes, a Lagrangian, a Z-vector solve and a
   !! relaxed density, and a disagreement could be any of them. This needs a
   !! converged density and nothing else.
   !!
   !! **The finite difference has to move the grid.** `Tr(P V_xc[D])` depends on
   !! the nuclear positions through the basis functions, through the atom-centred
   !! quadrature, and through the Becke weights -- so the displaced evaluation
   !! rebuilds the molecule and its grid, and reuses only the two density
   !! *matrices*. Holding the grid fixed would check a different quantity and
   !! would agree with a routine missing two of its three terms.
   !!
   !! **`P` is deliberately not a density.** It stands in for the relaxed
   !! correlation density plus the Z-vector, which is symmetric, indefinite and
   !! integrates to zero. A test that passed `D` twice would be blind to every
   !! place the code assumed a positive or idempotent second argument, so what
   !! goes in here is a fixed pseudo-random symmetric matrix -- reproducible, and
   !! resembling nothing the routine could special-case.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_gradient, only: xc_potential_gradient
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_available, &
                             xc_add_potential
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none

   integer, parameter :: N_DIM = 3
   real(dp), parameter :: STEP = 1.0e-3_dp
      !! Central difference. Smaller starts to show the grid's own noise -- the
      !! quadrature is not a smooth function of the geometry at the level of the
      !! twelfth decimal -- and larger shows the fourth-order term.
   real(dp), parameter :: TOL = 2.0e-7_dp
   integer :: n_bad

   n_bad = 0
   write (*, "(a)") "== the derivative of the exchange-correlation potential"

   if (.not. xc_available()) then
      write (*, "(a)") "   skipped: this build has no libxc"
      stop 0
   end if

   ! LDA first: the second-derivative channel and the weight term with no
   ! gradient machinery on top of them.
   call check_case("H2O / sto-3g, SVWN", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", 10, "svwn", n_bad)

   ! A pure GGA, where all four channels are live and the two that carry
   ! `v_sigma` cannot hide behind the ones that carry `v_rho`.
   call check_case("H2O / sto-3g, BLYP", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", 10, "blyp", n_bad)

   ! Asymmetric, and in a basis with p functions on two centres: a term that
   ! pairs the wrong index of the second derivative survives a symmetric
   ! molecule and does not survive this one.
   call check_case("HCN / sto-3g, BLYP", [1, 6, 7], ["H", "C", "N"], &
                   reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                   "sto-3g", 14, "blyp", n_bad)

   ! B2PLYP's own semilocal part, at the weight it carries there. The
   ! coefficient is what the double hybrid will actually differentiate, and a
   ! functional weight dropped somewhere in the per-component loop shows up
   ! here and nowhere else.
   call check_case("H2O / sto-3g, B2PLYP semilocal", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", 10, "b2plyp", n_bad)

   write (*, "(a)") ""
   if (n_bad == 0) then
      write (*, "(a)") "all exchange-correlation potential gradient checks passed"
   else
      write (*, "(a,i0,a)") "FAILED: ", n_bad, " case(s)"
      stop 1
   end if

contains

   subroutine check_case(label, numbers, symbols, coords, basis, nelec, functional, n_bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nelec
      character(len=*), intent(in) :: functional
      integer, intent(inout) :: n_bad

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: xc
      type(error_t) :: error
      real(dp), allocatable :: analytic(:, :), numeric(:, :)
      real(dp), allocatable :: pmat(:, :), dens(:, :)
      real(dp), allocatable :: shifted(:, :)
      real(dp) :: plus, minus, worst
      integer :: natm, ia, comp, nao

      natm = size(numbers)
      write (*, "(a)") ""
      write (*, "(a,a)") "== ", label
      flush (output_unit)

      call build_libcint_molecule(numbers, symbols, coords, basis, mol, error)
      if (fail(error, n_bad)) return
      call xc_context_create(mol, functional, xc, error, polarized=.false.)
      if (fail(error, n_bad)) return
      call run_libcint_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, error, xc=xc)
      if (fail(error, n_bad)) return

      nao = mol%nao
      allocate (dens(nao, nao))
      dens = scf%density
      call pseudo_random_symmetric(nao, pmat)

      allocate (analytic(3, natm))
      analytic = 0.0_dp
      call xc_potential_gradient(xc, mol, dens, pmat, analytic, error)
      if (fail(error, n_bad)) return
      call xc%destroy()
      call mol%destroy()

      allocate (numeric(3, natm), shifted(N_DIM, natm))
      do ia = 1, natm
         do comp = 1, 3
            shifted = coords
            shifted(comp, ia) = coords(comp, ia) + STEP
            call trace_at(numbers, symbols, shifted, basis, functional, dens, pmat, &
                          plus, error)
            if (fail(error, n_bad)) return
            shifted(comp, ia) = coords(comp, ia) - STEP
            call trace_at(numbers, symbols, shifted, basis, functional, dens, pmat, &
                          minus, error)
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

      ! Translational invariance. Worth little on its own -- it held at 1e-15
      ! through three wrong versions of the RI-MP2 gradient -- but a term landing
      ! on the wrong atom fails it while the finite difference above still
      ! passes, so it costs nothing and covers a different mistake.
      if (maxval(abs(sum(analytic, dim=2))) > 1.0e-8_dp) then
         write (*, "(a,es14.4)") "  FAIL: does not sum to zero over atoms: ", &
            maxval(abs(sum(analytic, dim=2)))
         n_bad = n_bad + 1
      end if
   end subroutine check_case

   subroutine trace_at(numbers, symbols, coords, basis, functional, dens, pmat, &
                       value, error)
      !! `Tr(P V_xc[D])` at one geometry, with both matrices carried in
      !!
      !! The molecule and the grid are rebuilt; the densities are not. That is
      !! what makes this the *skeleton* derivative -- the one the analytic
      !! routine computes -- rather than the derivative of a converged quantity.
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, functional
      real(dp), intent(in) :: dens(:, :), pmat(:, :)
      real(dp), intent(out) :: value
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      type(xc_context_t) :: xc
      real(dp), allocatable :: v_xc(:, :)
      real(dp) :: e_xc, n_elec

      value = 0.0_dp
      call build_libcint_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) return
      call xc_context_create(mol, functional, xc, error, polarized=.false.)
      if (error%has_error()) return

      allocate (v_xc(mol%nao, mol%nao))
      call xc_add_potential(xc, mol, dens, v_xc, e_xc, n_elec, error)
      if (error%has_error()) return
      value = sum(pmat*v_xc)

      call xc%destroy()
      call mol%destroy()
   end subroutine trace_at

   subroutine pseudo_random_symmetric(n, p)
      !! A fixed symmetric matrix that is not a density
      !!
      !! Indefinite and traceless-ish by construction, so nothing downstream can
      !! quietly rely on `P` being a density: the real second argument is the
      !! relaxed correlation density plus a Z-vector, which is none of the things
      !! an SCF density is. Generated from a recurrence rather than
      !! `random_number` so a failure reproduces.
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: p(:, :)

      integer :: i, j
      real(dp) :: x

      allocate (p(n, n))
      x = 0.5_dp
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

end program check_xc_potential_gradient
