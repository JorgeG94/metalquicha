!! Second derivatives of the three- and two-centre fitting integrals
module test_mqc_hess_ints_df
   !! `hess_3c2c_contract` against differences of the first-derivative routines
   !!
   !! The reference is built from nothing the routine uses: `three_centre_deriv`
   !! and `two_centre_deriv`, the first derivatives the RI-MP2 gradient already
   !! runs on, contracted with the same density at displaced geometries and
   !! differenced with a five-point stencil. The density is held fixed, so what
   !! is differenced is exactly the function whose second derivative the
   !! routine claims to be.
   !!
   !! Water in cc-pVDZ with cc-pVDZ-RIFIT, so d shells on the orbital side and
   !! up to g on the fitting side, in both angular forms. Three densities, each
   !! for a reason:
   !!
   !!   * random and **not** symmetric in `mu nu`, as a correlation density is not:
   !!     the routine has to fold `g(ij)` and `g(ji)` itself, and a symmetric
   !!     density would let it get that wrong unnoticed;
   !!   * random and symmetric, where every orbit term is equal and a wrong
   !!     weight shows as a factor of two rather than as noise;
   !!   * built from SCF orbitals, `C_occ (.) C_vir^T` with smooth auxiliary
   !!     weights, which is the shape an RI-MP2 density actually has;
   !!   * the two terms separately, because a sign or a swap that cancels
   !!     between the three-centre and two-centre sums survives their total.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, atom_ao_blocks
   use mqc_czt_gradient, only: three_centre_deriv, two_centre_deriv
   use mqc_czt_hess_ints, only: hess_3c2c_contract
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   implicit none
   private

   public :: collect_mqc_hess_ints_df_tests

   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.05_dp, 0.02_dp, -0.03_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.04_dp, 1.4308_dp, 1.1078_dp], [3, 3])
   integer, parameter :: Z(3) = [8, 1, 1]
   character(len=2), parameter :: SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: STEP = 1.0e-3_dp
   real(dp), parameter :: TOL = 1.0e-9_dp
      !! Relative to the largest entry. Thirty times the worst measured (3.2e-11,
      !! Cartesian); the stencil's truncation error at this step is far below it
      !! and what it absorbs is the roundoff of a differenced sum.

contains

   subroutine collect_mqc_hess_ints_df_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("random_density_spherical", random_spherical), &
                  new_unittest("random_density_cartesian", random_cartesian), &
                  new_unittest("orbital_density_spherical", orbital_density), &
                  new_unittest("each_term_alone", each_term_alone), &
                  new_unittest("translations_leave_it_alone", translations), &
                  new_unittest("a_bad_shape_is_refused", bad_shape) &
                  ]
   end subroutine collect_mqc_hess_ints_df_tests

   subroutine build_pair(geo, cartesian, orb, aux, err)
      real(dp), intent(in) :: geo(3, 3)
      logical, intent(in) :: cartesian
      type(czt_molecule_t), intent(out) :: orb, aux
      type(error_t), intent(inout) :: err

      call build_czt_molecule(Z, SYM, geo, "cc-pvdz", orb, err, force_cartesian=cartesian)
      if (err%has_error()) return
      call build_czt_molecule(Z, SYM, geo, "cc-pvdz-rifit", aux, err, force_cartesian=cartesian)
   end subroutine build_pair

   subroutine owner_of(mol, owner)
      !! Owning atom per AO, 1-based
      type(czt_molecule_t), intent(in) :: mol
      integer, allocatable, intent(out) :: owner(:)

      integer, allocatable :: offsets(:), counts(:)
      integer :: c

      allocate (offsets(mol%natm), counts(mol%natm), owner(mol%nao))
      call atom_ao_blocks(mol, offsets, counts)
      do c = 1, mol%natm
         owner(offsets(c) + 1:offsets(c) + counts(c)) = c
      end do
   end subroutine owner_of

   subroutine contraction(geo, cartesian, g3, g2, vec, err)
      !! `sum g3 . d(mu nu|P)/dX + sum g2 . d(P|Q)/dX` for every nuclear
      !! coordinate `X = 3*(A-1) + c`, from the first-derivative routines
      !!
      !! `three_centre_deriv` and `two_centre_deriv` differentiate with respect
      !! to an electron coordinate, which is minus the nuclear one, and give the
      !! first function's centre only: the second orbital index is the transpose
      !! and the auxiliary index has its own entry point.
      real(dp), intent(in) :: geo(3, 3)
      logical, intent(in) :: cartesian
      real(dp), intent(in) :: g3(:, :, :), g2(:, :)
      real(dp), intent(out) :: vec(9)
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: orb, aux
      real(dp), allocatable :: ip1(:, :, :, :), ip2(:, :, :, :), ip2c(:, :, :)
      integer, allocatable :: own_o(:), own_a(:)
      integer :: mu, nu, p, q, c, atom

      vec = 0.0_dp
      call build_pair(geo, cartesian, orb, aux, err)
      if (err%has_error()) return
      call three_centre_deriv(orb, aux, 1, ip1)
      call three_centre_deriv(orb, aux, 2, ip2)
      call two_centre_deriv(aux, ip2c)
      call owner_of(orb, own_o)
      call owner_of(aux, own_a)

      do c = 1, 3
         do p = 1, aux%nao
            do nu = 1, orb%nao
               do mu = 1, orb%nao
                  atom = own_o(mu)
                  vec(3*(atom - 1) + c) = vec(3*(atom - 1) + c) &
                                          - g3(mu, nu, p)*ip1(mu, nu, p, c)
                  atom = own_o(nu)
                  vec(3*(atom - 1) + c) = vec(3*(atom - 1) + c) &
                                          - g3(mu, nu, p)*ip1(nu, mu, p, c)
                  atom = own_a(p)
                  vec(3*(atom - 1) + c) = vec(3*(atom - 1) + c) &
                                          - g3(mu, nu, p)*ip2(mu, nu, p, c)
               end do
            end do
         end do
         do q = 1, aux%nao
            do p = 1, aux%nao
               atom = own_a(p)
               vec(3*(atom - 1) + c) = vec(3*(atom - 1) + c) - g2(p, q)*ip2c(p, q, c)
               atom = own_a(q)
               vec(3*(atom - 1) + c) = vec(3*(atom - 1) + c) - g2(p, q)*ip2c(q, p, c)
            end do
         end do
      end do
      call orb%destroy()
      call aux%destroy()
   end subroutine contraction

   subroutine differenced(cartesian, g3, g2, hfd, err)
      !! `hfd(X, Y)`: five-point derivative of `contraction` along `Y`
      logical, intent(in) :: cartesian
      real(dp), intent(in) :: g3(:, :, :), g2(:, :)
      real(dp), intent(out) :: hfd(9, 9)
      type(error_t), intent(inout) :: err

      real(dp) :: geo(3, 3), v(9, 4)
      real(dp), parameter :: SHIFT(4) = [2.0_dp, 1.0_dp, -1.0_dp, -2.0_dp]
      real(dp), parameter :: WEIGHT(4) = [-1.0_dp, 8.0_dp, -8.0_dp, 1.0_dp]
      integer :: y, k, atom, comp

      do y = 1, 9
         atom = (y - 1)/3 + 1
         comp = y - 3*(atom - 1)
         do k = 1, 4
            geo = WATER
            geo(comp, atom) = geo(comp, atom) + SHIFT(k)*STEP
            call contraction(geo, cartesian, g3, g2, v(:, k), err)
            if (err%has_error()) return
         end do
         hfd(:, y) = (WEIGHT(1)*v(:, 1) + WEIGHT(2)*v(:, 2) + WEIGHT(3)*v(:, 3) &
                      + WEIGHT(4)*v(:, 4))/(12.0_dp*STEP)
      end do
   end subroutine differenced

   subroutine analytic(cartesian, g3, g2, hess, err, with3, with2)
      logical, intent(in) :: cartesian
      real(dp), intent(in) :: g3(:, :, :), g2(:, :)
      real(dp), intent(out) :: hess(9, 9)
      type(error_t), intent(inout) :: err
      logical, intent(in), optional :: with3, with2

      type(czt_molecule_t) :: orb, aux
      real(dp), allocatable :: h4(:, :, :, :)
      logical :: do3, do2
      integer :: a, b, ia, ib

      do3 = .true.
      do2 = .true.
      if (present(with3)) do3 = with3
      if (present(with2)) do2 = with2
      call build_pair(WATER, cartesian, orb, aux, err)
      if (err%has_error()) return
      allocate (h4(3, 3, 3, 3))
      h4 = 0.0_dp
      if (do3 .and. do2) then
         call hess_3c2c_contract(orb, aux, h4, err, gamma3=g3, gamma2=g2)
      else if (do3) then
         call hess_3c2c_contract(orb, aux, h4, err, gamma3=g3)
      else
         call hess_3c2c_contract(orb, aux, h4, err, gamma2=g2)
      end if
      do ib = 1, 3
         do b = 1, 3
            do ia = 1, 3
               do a = 1, 3
                  hess(3*(ia - 1) + a, 3*(ib - 1) + b) = h4(a, b, ia, ib)
               end do
            end do
         end do
      end do
      call orb%destroy()
      call aux%destroy()
   end subroutine analytic

   subroutine random_densities(n_ao, n_aux, symmetric, g3, g2)
      integer, intent(in) :: n_ao, n_aux
      logical, intent(in) :: symmetric
      real(dp), allocatable, intent(out) :: g3(:, :, :), g2(:, :)

      real(dp), allocatable :: t3(:, :, :), t2(:, :)

      allocate (g3(n_ao, n_ao, n_aux), g2(n_aux, n_aux), t3(n_ao, n_ao, n_aux), t2(n_aux, n_aux))
      call random_seed()
      call random_number(g3)
      call random_number(g2)
      g3 = g3 - 0.5_dp
      g2 = g2 - 0.5_dp
      if (symmetric) then
         t3 = reshape(g3, shape(g3), order=[2, 1, 3])
         g3 = 0.5_dp*(g3 + t3)
      end if
      t2 = transpose(g2)
      g2 = 0.5_dp*(g2 + t2)
   end subroutine random_densities

   subroutine compare(error, label, h_an, h_fd)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: h_an(9, 9), h_fd(9, 9)

      real(dp) :: scale, dev

      scale = max(1.0_dp, maxval(abs(h_fd)))
      dev = maxval(abs(h_an - h_fd))
      call check(error, dev < TOL*scale, label//": analytic against differenced first "// &
                 "derivatives; largest entry deviates")
      if (allocated(error)) then
         write (*, "(a,a,es10.3,a,es10.3)") label, " dev ", dev, "  scale ", scale
      end if
   end subroutine compare

   subroutine one_density_case(error, cartesian, symmetric, label)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: cartesian, symmetric
      character(len=*), intent(in) :: label

      type(czt_molecule_t) :: orb, aux
      type(error_t) :: err
      real(dp), allocatable :: g3(:, :, :), g2(:, :)
      real(dp) :: h_an(9, 9), h_fd(9, 9)

      call build_pair(WATER, cartesian, orb, aux, err)
      call check(error,.not. err%has_error(), "water and its fitting set must build")
      if (allocated(error)) return
      call random_densities(orb%nao, aux%nao, symmetric, g3, g2)
      call orb%destroy()
      call aux%destroy()

      call analytic(cartesian, g3, g2, h_an, err)
      call differenced(cartesian, g3, g2, h_fd, err)
      call check(error,.not. err%has_error(), "the reference must run: "//err%get_full_trace())
      if (allocated(error)) return
      call compare(error, label, h_an, h_fd)
   end subroutine one_density_case

   subroutine random_spherical(error)
      type(error_type), allocatable, intent(out) :: error

      call one_density_case(error, .false., .false., "spherical, unsymmetric density")
      if (allocated(error)) return
      call one_density_case(error, .false., .true., "spherical, symmetric density")
   end subroutine random_spherical

   subroutine random_cartesian(error)
      !! The `_cart` kernels: a different set of entry points, not a flag
      type(error_type), allocatable, intent(out) :: error

      call one_density_case(error, .true., .false., "Cartesian, unsymmetric density")
   end subroutine random_cartesian

   subroutine orbital_density(error)
      !! `Gamma^P_{mu nu} = sum_ia C_{mu i} C_{nu a} w^P_{ia}`, `gamma_PQ` smooth
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: orb, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: g3(:, :, :), g2(:, :)
      real(dp) :: h_an(9, 9), h_fd(9, 9), w
      integer :: i, a, p, q, nocc

      nocc = 5
      call build_pair(WATER, .false., orb, aux, err)
      call check(error,.not. err%has_error(), "water and its fitting set must build")
      if (allocated(error)) return
      call run_czt_rhf(orb, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error(), "the SCF must converge: "//err%get_full_trace())
      if (allocated(error)) return

      allocate (g3(orb%nao, orb%nao, aux%nao), g2(aux%nao, aux%nao))
      g3 = 0.0_dp
      do p = 1, aux%nao
         do a = nocc + 1, orb%nao
            do i = 1, nocc
               w = cos(0.3_dp*p + 0.7_dp*i - 1.1_dp*a)/(1.0_dp + 0.1_dp*p)
               g3(:, :, p) = g3(:, :, p) + w*spread(scf%orbitals(:, i), 2, orb%nao) &
                             *spread(scf%orbitals(:, a), 1, orb%nao)
            end do
         end do
      end do
      do q = 1, aux%nao
         do p = 1, aux%nao
            g2(p, q) = exp(-abs(real(p - q, dp))/6.0_dp)
         end do
      end do
      call orb%destroy()
      call aux%destroy()

      call analytic(.false., g3, g2, h_an, err)
      call differenced(.false., g3, g2, h_fd, err)
      call check(error,.not. err%has_error(), "the reference must run: "//err%get_full_trace())
      if (allocated(error)) return
      call compare(error, "orbital-built density", h_an, h_fd)
   end subroutine orbital_density

   subroutine each_term_alone(error)
      !! The three- and two-centre sums separately, each against its own
      !! reference, so that neither can hide the other's error
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: orb, aux
      type(error_t) :: err
      real(dp), allocatable :: g3(:, :, :), g2(:, :), z3(:, :, :), z2(:, :)
      real(dp) :: h_an(9, 9), h_fd(9, 9)

      call build_pair(WATER, .false., orb, aux, err)
      call check(error,.not. err%has_error(), "water and its fitting set must build")
      if (allocated(error)) return
      call random_densities(orb%nao, aux%nao, .false., g3, g2)
      allocate (z3(orb%nao, orb%nao, aux%nao), z2(aux%nao, aux%nao))
      z3 = 0.0_dp
      z2 = 0.0_dp
      call orb%destroy()
      call aux%destroy()

      call analytic(.false., g3, z2, h_an, err, with3=.true., with2=.false.)
      call differenced(.false., g3, z2, h_fd, err)
      call compare(error, "three-centre term alone", h_an, h_fd)
      if (allocated(error)) return

      call analytic(.false., z3, g2, h_an, err, with3=.false., with2=.true.)
      call differenced(.false., z3, g2, h_fd, err)
      call check(error,.not. err%has_error(), "the reference must run: "//err%get_full_trace())
      if (allocated(error)) return
      call compare(error, "two-centre term alone", h_an, h_fd)
   end subroutine each_term_alone

   subroutine translations(error)
      !! Moving every atom together moves nothing, so each row of the Hessian
      !! sums to zero over the second atom. This needs no reference at all, and
      !! it is what a missing translational partner or a dropped transpose
      !! placement breaks first.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: orb, aux
      type(error_t) :: err
      real(dp), allocatable :: g3(:, :, :), g2(:, :), h4(:, :, :, :)
      real(dp) :: worst, row
      integer :: a, b, ia

      call build_pair(WATER, .false., orb, aux, err)
      call check(error,.not. err%has_error(), "water and its fitting set must build")
      if (allocated(error)) return
      call random_densities(orb%nao, aux%nao, .false., g3, g2)
      allocate (h4(3, 3, 3, 3))
      h4 = 0.0_dp
      call hess_3c2c_contract(orb, aux, h4, err, gamma3=g3, gamma2=g2)
      call check(error,.not. err%has_error(), "the contraction must run")
      if (allocated(error)) return

      worst = 0.0_dp
      do ia = 1, 3
         do a = 1, 3
            do b = 1, 3
               row = sum(h4(a, b, ia, :))
               worst = max(worst, abs(row))
            end do
         end do
      end do
      call check(error, worst < 1.0e-9_dp, "translating all atoms must change nothing")
      if (allocated(error)) return

      worst = 0.0_dp
      do ia = 1, 3
         do b = 1, 3
            do a = 1, 3
               worst = max(worst, maxval(abs(h4(a, b, ia, :) - h4(b, a, :, ia))))
            end do
         end do
      end do
      call check(error, worst < 1.0e-9_dp, "the Hessian must be symmetric")
      call orb%destroy()
      call aux%destroy()
   end subroutine translations

   subroutine bad_shape(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: orb, aux
      type(error_t) :: err
      real(dp), allocatable :: h4(:, :, :, :), g3(:, :, :)

      call build_pair(WATER, .false., orb, aux, err)
      call check(error,.not. err%has_error(), "water and its fitting set must build")
      if (allocated(error)) return
      allocate (h4(3, 3, 3, 3), g3(orb%nao, orb%nao, aux%nao - 1))
      h4 = 0.0_dp
      g3 = 0.0_dp
      call hess_3c2c_contract(orb, aux, h4, err, gamma3=g3)
      call check(error, err%has_error(), "a density of the wrong shape must be refused")
      call orb%destroy()
      call aux%destroy()
   end subroutine bad_shape

end module test_mqc_hess_ints_df

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_hess_ints_df, only: collect_mqc_hess_ints_df_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_hess_ints_df", collect_mqc_hess_ints_df_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
