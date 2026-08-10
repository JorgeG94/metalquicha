program check_direct
   !! Check the direct Fock build against the in-core one, elementwise
   !!
   !! The in-core path is slow but validated against PySCF, which makes it the
   !! right oracle for this: the direct build must reproduce it to machine
   !! precision, on the same density, for every element.
   !!
   !! The density used is deliberately not a converged one. A converged density
   !! is close to idempotent and nearly diagonal in the MO basis, which can hide
   !! errors in off-diagonal Fock blocks -- exactly where a degeneracy-factor
   !! mistake shows up. A dense pseudo-random symmetric density exercises every
   !! block equally.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_direct, only: schwarz_bounds, build_fock_direct, direct_stats_t, &
                                 DEFAULT_SCREEN_TOL
   implicit none

   integer, parameter :: N_DIM = 3
   integer, parameter :: N_TOL = 3
   real(dp), parameter :: SCREEN_TOLS(N_TOL) = [1.0e-14_dp, 1.0e-11_dp, 1.0e-6_dp]

   integer :: n_bad

   n_bad = 0

   ! A single water is compact: every shell pair overlaps every other, so
   ! nothing screens and the test is purely of the symmetry bookkeeping.
   call run_case("water", ["O", "H", "H"], &
                 reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, -1.4308_dp, 1.1078_dp, &
                          0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                 [8, 1, 1], "cc-pvdz", .false., n_bad)

   ! Four waters strung out along z, 12 Bohr apart. Screening only pays when
   ! the system is larger than the range of the basis functions, which is the
   ! whole reason direct SCF beats storing integrals -- so it has to be tested
   ! on something extended, not on a molecule that fits inside one orbital.
   call run_case("chain4", [character(len=2) :: "O", "H", "H", "O", "H", "H", &
                            "O", "H", "H", "O", "H", "H"], &
                 water_chain(4), &
                 [8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1], "cc-pvdz", .true., n_bad)

   if (n_bad > 0) then
      write (*, "(a)") "FAILED"
      stop 1
   end if
   write (*, "(a)") "direct Fock matches the in-core build"

contains

   function water_chain(n_units) result(coords)
      !! `n_units` waters along z, far enough apart that screening has work to do
      integer, intent(in) :: n_units
      real(dp), allocatable :: coords(:, :)
      real(dp), parameter :: SPACING = 12.0_dp
      integer :: u, base

      allocate (coords(N_DIM, 3*n_units))
      do u = 1, n_units
         base = 3*(u - 1)
         coords(:, base + 1) = [0.0_dp, 0.0_dp, SPACING*real(u - 1, dp)]
         coords(:, base + 2) = [0.0_dp, -1.4308_dp, SPACING*real(u - 1, dp) + 1.1078_dp]
         coords(:, base + 3) = [0.0_dp, 1.4308_dp, SPACING*real(u - 1, dp) + 1.1078_dp]
      end do
   end function water_chain

   subroutine run_case(label, symbols, coords, atomic_numbers, basis, expect_screening, n_bad)
      character(len=*), intent(in) :: label
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: basis
      logical, intent(in) :: expect_screening  !! Whether the system is extended enough
      integer, intent(inout) :: n_bad

      type(libcint_molecule_t) :: mol
      type(direct_stats_t) :: stats
      type(error_t) :: error
      real(dp), allocatable :: h(:, :), eri(:, :, :, :), density(:, :)
      real(dp), allocatable :: fock_ref(:, :), fock_dir(:, :), bounds(:, :)
      real(dp) :: t0, t1, t_incore, t_direct, worst
      integer :: n, k

      call build_libcint_molecule(atomic_numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL building molecule: ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      n = mol%nao
      write (*, "(a,a,a,a,a,i0,a,i0,a)") "== ", label, "/", basis, ": ", n, &
         " basis functions, ", mol%nbas, " shells"

      allocate (h(n, n), density(n, n), fock_ref(n, n), fock_dir(n, n))
      call mol%core_hamiltonian(h)
      call pseudo_density(density)

      ! Reference: the in-core build, n^4 tensor and quadruple loop.
      call cpu_time(t0)
      call mol%eris(eri)
      call build_fock_incore(h, eri, density, fock_ref)
      call cpu_time(t1)
      t_incore = t1 - t0

      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL Schwarz bounds: ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      write (*, "(a)") "  tol        computed    screened   %screened   time(s)   max |dF|"
      do k = 1, N_TOL
         call cpu_time(t0)
         call build_fock_direct(mol, h, density, bounds, fock_dir, stats, error, &
                                screen_tol=SCREEN_TOLS(k))
         call cpu_time(t1)
         t_direct = t1 - t0

         if (error%has_error()) then
            write (*, "(a,a)") "FAIL direct build: ", error%get_message()
            n_bad = n_bad + 1
            return
         end if

         worst = maxval(abs(fock_dir - fock_ref))
         write (*, "(a,es8.1,i12,i12,f11.1,f10.4,es12.3)") "  ", SCREEN_TOLS(k), &
            stats%quartets_computed, stats%quartets_screened, &
            100.0_dp*stats%screened_fraction(), t_direct, worst

         ! At a screening tolerance of 1e-14 the two builds must agree to
         ! rounding: nothing meaningful has been discarded. The loose tolerance
         ! is expected to differ, and is here to show screening is actually
         ! doing something rather than being silently inert.
         if (k == 1 .and. worst > 1.0e-12_dp) then
            write (*, "(a,es12.3)") "   FAIL: direct and in-core disagree by ", worst
            n_bad = n_bad + 1
         end if
         ! Only assert this where the geometry gives screening something to do.
         if (k == N_TOL .and. expect_screening .and. stats%screened_fraction() < 0.1_dp) then
            write (*, "(a,f6.1,a)") "   FAIL: loose screening removed only ", &
               100.0_dp*stats%screened_fraction(), "% on an extended system"
            n_bad = n_bad + 1
         end if
      end do

      write (*, "(a,f8.4,a,f8.4,a,f6.2,a)") "  in-core ", t_incore, " s vs direct ", &
         t_direct, " s  (", t_incore/max(t_direct, 1.0e-9_dp), "x)"
      write (*, "(a,i0,a,f9.3,a)") "  in-core tensor: ", n, "^4 = ", &
         real(n, dp)**4*8.0_dp/1.0e6_dp, " MB"

      call mol%destroy()
   end subroutine run_case

   subroutine pseudo_density(density)
      !! A dense symmetric matrix standing in for a density
      !!
      !! Deterministic, so a failure reproduces. Not idempotent and not close to
      !! a real density -- that is the point, see the header.
      real(dp), intent(out) :: density(:, :)
      integer :: i, j, n

      n = size(density, 1)
      do j = 1, n
         do i = 1, n
            density(i, j) = sin(real(i*7 + j*13, dp))*exp(-0.05_dp*abs(i - j))
         end do
      end do
      density = 0.5_dp*(density + transpose(density))
   end subroutine pseudo_density

   pure subroutine build_fock_incore(h, eri, density, fock)
      !! The reference build, copied from mqc_libcint_rhf
      real(dp), intent(in) :: h(:, :), eri(:, :, :, :), density(:, :)
      real(dp), intent(out) :: fock(:, :)

      integer :: mu, nu, la, si, n

      n = size(h, 1)
      fock = h
      do nu = 1, n
         do mu = 1, n
            do si = 1, n
               do la = 1, n
                  fock(mu, nu) = fock(mu, nu) + density(la, si) &
                                 *(eri(mu, nu, la, si) - 0.5_dp*eri(mu, la, nu, si))
               end do
            end do
         end do
      end do
   end subroutine build_fock_incore

end program check_direct
