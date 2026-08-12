!! Coupled-perturbed Hartree-Fock
module test_mqc_libcint_cphf
   !! The analytical polarizability against the field it is the derivative of.
   !!
   !! This is the rare response property that can be checked to convergence
   !! without any external reference, because the derivative it computes is one
   !! the SCF can be made to take numerically: apply a field, converge, read off
   !! the dipole, difference. So these tests are here rather than in
   !! `validation/`, which needs PySCF -- an agreement between an analytical
   !! derivative and a numerical one is worth as much as an agreement with
   !! another program, and it runs on every commit.
   !!
   !! **The dipole is differenced, not the energy.** Both would work, but
   !! `alpha = -d2E/dF2` needs a second difference, whose error grows as the SCF
   !! threshold over the step squared -- forcing a large step, whose truncation
   !! error is then the limit. Differencing the dipole is one derivative instead
   !! of two, so a 1e-3 field and an SCF converged to 1e-12 leave the check
   !! limited by the 1e-6 truncation of the central difference, which is six
   !! digits of agreement rather than three.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_multipole, only: multipole_matrices
   use mqc_libcint_cphf, only: cphf_solve, static_polarizability
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_libcint_cphf_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

   !> Field strength for the numerical derivative, atomic units.
   !>
   !> Small enough that the cubic term (the hyperpolarizability) contributes
   !> below the tolerance, large enough that the SCF threshold does not. Water
   !> has a first hyperpolarizability of order 10 a.u., so at 1e-3 the cubic
   !> contamination of a central difference is around 1e-5 relative -- and it
   !> cancels to leading order anyway, being an even-order effect in a central
   !> difference of the dipole.
   real(dp), parameter :: FIELD_STRENGTH = 1.0e-3_dp

contains

   subroutine collect_mqc_libcint_cphf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("polarizability_matches_a_finite_field", &
                               test_finite_field), &
                  new_unittest("polarizability_is_symmetric", test_symmetry), &
                  new_unittest("coupling_shifts_the_polarizability", &
                               test_coupling_matters), &
                  new_unittest("response_is_orthogonal_to_the_occupied_space", &
                               test_origin_independence), &
                  new_unittest("cphf_refuses_a_closed_shell_it_cannot_solve", &
                               test_refusals) &
                  ]
   end subroutine collect_mqc_libcint_cphf_tests

   subroutine water(mol, scf, err, field)
      !! Water in STO-3G, optionally in a uniform electric field
      type(libcint_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err
      real(dp), intent(in), optional :: field(3)

      real(dp) :: c(3, 3)
      real(dp), allocatable :: dip(:, :, :), h_field(:, :)
      integer :: k

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_libcint_molecule([8, 1, 1], ["O ", "H ", "H "], c, "sto-3g", mol, err)
      if (err%has_error()) return

      if (present(field)) then
         call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, err)
         if (err%has_error()) return
         ! The electronic dipole operator is -r, so a field F perturbs the
         ! Hamiltonian by +r.F. Getting this sign wrong would flip the numerical
         ! derivative and turn the test into a mirror of itself, which is why the
         ! sign of alpha is also asserted below: a positive definite result is
         ! not something a sign error can fake.
         allocate (h_field(mol%nao, mol%nao))
         h_field = 0.0_dp
         do k = 1, 3
            h_field = h_field + field(k)*dip(:, :, k)
         end do
         call run_libcint_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, &
                              in_core=.true., h_extra=h_field)
         deallocate (dip, h_field)
      else
         call run_libcint_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, &
                              in_core=.true.)
      end if
   end subroutine water

   subroutine dipole_in_field(field, mu, err)
      !! Converge water in a field and return its total dipole, in a.u.
      real(dp), intent(in) :: field(3)
      real(dp), intent(out) :: mu(3)
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: dip(:, :, :)
      integer :: k, iatom

      mu = 0.0_dp
      call water(mol, scf, err, field=field)
      if (err%has_error()) return
      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, err)
      if (err%has_error()) return

      do k = 1, 3
         mu(k) = -sum(scf%density*dip(:, :, k))
         do iatom = 1, mol%natm
            mu(k) = mu(k) + mol%charges(iatom)*mol%coords(k, iatom)
         end do
      end do
      call mol%destroy()
      deallocate (dip)
   end subroutine dipole_in_field

   subroutine test_finite_field(error)
      !! alpha_kl = d mu_k / d F_l, by central difference
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp) :: alpha(3, 3), numeric(3, 3), field(3), plus(3), minus(3)
      integer :: k, l

      call water(mol, scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         call check(error, .false., "the reference SCF failed")
         return
      end if
      call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                 scf%n_occupied, alpha, err)
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "CPHF failed: "//err%get_message())
         return
      end if

      do l = 1, 3
         field = 0.0_dp
         field(l) = FIELD_STRENGTH
         call dipole_in_field(field, plus, err)
         field(l) = -FIELD_STRENGTH
         call dipole_in_field(field, minus, err)
         if (err%has_error()) then
            call check(error, .false., "a field SCF failed: "//err%get_message())
            return
         end if
         do k = 1, 3
            numeric(k, l) = (plus(k) - minus(k))/(2.0_dp*FIELD_STRENGTH)
         end do
      end do

      ! 1e-6 absolute on numbers of order 1 to 10: the central difference is
      ! O(FIELD_STRENGTH^2) truncated, so this is the truncation error and not slack.
      do l = 1, 3
         do k = 1, 3
            call check(error, alpha(k, l), numeric(k, l), &
                       "analytical and finite-field polarizabilities disagree", &
                       thr=2.0e-6_dp)
            if (allocated(error)) return
         end do
      end do

      ! And a polarizability is positive definite, so no sign convention can
      ! have cancelled between the two calculations above.
      do k = 1, 3
         call check(error, alpha(k, k) > 0.0_dp, &
                    "a diagonal polarizability came out negative")
         if (allocated(error)) return
      end do
   end subroutine test_finite_field

   subroutine test_symmetry(error)
      !! alpha is a second derivative of one scalar, so it is symmetric
      !!
      !! Not automatic: it comes out of contracting three separate CG solutions
      !! against three separate perturbation matrices, and a transposed dipole
      !! block or a mismatched index in that contraction would break it while
      !! leaving the diagonal -- and therefore the isotropic average that a user
      !! looks at -- entirely plausible.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp) :: alpha(3, 3)

      call water(mol, scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         call check(error, .false., "the reference SCF failed")
         return
      end if
      call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                 scf%n_occupied, alpha, err)
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "CPHF failed: "//err%get_message())
         return
      end if

      call check(error, alpha(1, 2), alpha(2, 1), "alpha is not symmetric", thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, alpha(1, 3), alpha(3, 1), "alpha is not symmetric", thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, alpha(2, 3), alpha(3, 2), "alpha is not symmetric", thr=1.0e-9_dp)
      if (allocated(error)) return

      ! Water lies in the xz plane here, so y is the out-of-plane direction and
      ! must be the least polarizable. A tensor that came back in the wrong
      ! frame would pass every symmetry check above and fail this.
      call check(error, alpha(2, 2) < alpha(1, 1), &
                 "the out-of-plane direction is not the least polarizable")
      if (allocated(error)) return
      call check(error, alpha(2, 2) < alpha(3, 3), &
                 "the out-of-plane direction is not the least polarizable")
   end subroutine test_symmetry

   subroutine test_coupling_matters(error)
      !! The coupling moves the polarizability by tens of percent
      !!
      !! The uncoupled polarizability is the sum-over-states expression with bare
      !! orbital-energy denominators, `4 sum_ai h_ai^2 / (eps_a - eps_i)`, which
      !! is what solving nothing at all would give. Comparing against it is the
      !! only way to show the CG solve does work: a solver that returned its
      !! preconditioned right-hand side would pass symmetry, pass origin
      !! independence, and be wrong by tens of percent. So the assertion is on the
      !! *size* of the difference.
      !!
      !! **The direction is asserted too, but it is a property of this system and
      !! not a law.** Coupling raises the polarizability here -- by 20% in x, 29%
      !! in z and 90% in the nearly rigid out-of-plane direction. The reason is
      !! visible in the diagonal of the response operator,
      !!
      !!     A_ai,ai = (eps_a - eps_i) + 3(ai|ai) - (aa|ii)
      !!
      !! where `(aa|ii)`, the Coulomb repulsion between the occupied and virtual
      !! densities, is the largest of the three and enters negatively. The
      !! operator therefore sits *below* its own uncoupled diagonal, its inverse
      !! above, and the polarizability above -- which is the familiar statement
      !! that uncoupled Hartree-Fock underestimates polarizabilities. Whether the
      !! exchange or the Coulomb term wins is basis dependent, so this inequality
      !! is a regression guard on water in STO-3G rather than something to expect
      !! everywhere. What fixes the *sign* of the two-electron term for good is
      !! the finite-field test, not this one.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: dip(:, :, :), h(:, :, :), work(:, :)
      real(dp) :: coupled(3, 3), uncoupled, shift, gap
      integer :: k, a, i, n_occ, n_vir

      call water(mol, scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         call check(error, .false., "the reference SCF failed")
         return
      end if
      n_occ = scf%n_occupied
      n_vir = size(scf%orbitals, 2) - n_occ

      call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                 n_occ, coupled, err)
      if (.not. err%has_error()) then
         call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, err)
      end if
      if (err%has_error()) then
         call check(error, .false., "CPHF failed: "//err%get_message())
         call mol%destroy()
         return
      end if

      ! The occupied-virtual dipole blocks, for the uncoupled sum.
      allocate (h(n_vir, n_occ, 3), work(mol%nao, n_occ))
      do k = 1, 3
         work = matmul(dip(:, :, k), scf%orbitals(:, 1:n_occ))
         h(:, :, k) = matmul(transpose(scf%orbitals(:, n_occ + 1:)), work)
      end do
      call mol%destroy()

      shift = 0.0_dp
      do k = 1, 3
         uncoupled = 0.0_dp
         do i = 1, n_occ
            do a = 1, n_vir
               gap = scf%orbital_energies(n_occ + a) - scf%orbital_energies(i)
               uncoupled = uncoupled + 4.0_dp*h(a, i, k)**2/gap
            end do
         end do
         call check(error, coupled(k, k) > uncoupled, &
                    "coupling did not raise the polarizability on this system, so "// &
                    "the two-electron response has changed sign or scale")
         if (allocated(error)) return
         shift = max(shift, (coupled(k, k) - uncoupled)/uncoupled)
      end do

      call check(error, shift > 0.05_dp, &
                 "the coupled and uncoupled polarizabilities are suspiciously "// &
                 "close, so the CG solve may not be doing anything")

      deallocate (dip, h, work)
   end subroutine test_coupling_matters

   subroutine test_origin_independence(error)
      !! Shifting the dipole origin cannot move the polarizability
      !!
      !! Because the perturbation changes by a multiple of the overlap and the
      !! occupied-virtual block of the overlap vanishes. Testing it tests that
      !! claim, and it is the cheapest available check that the occupied-virtual
      !! transform is really occupied against virtual: if the block were built
      !! from overlapping index ranges, the overlap contribution would survive
      !! and a shifted origin would change the answer.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: dip(:, :, :), shifted(:, :, :), u(:, :, :), v(:, :, :)
      real(dp), allocatable :: ovl(:, :)
      real(dp) :: worst
      integer :: k

      call water(mol, scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         call check(error, .false., "the reference SCF failed")
         return
      end if

      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, err)
      call mol%overlap(ovl)
      allocate (shifted(mol%nao, mol%nao, 3))
      ! r - R for a displacement R, which is what a shifted origin produces.
      do k = 1, 3
         shifted(:, :, k) = dip(:, :, k) - real(k, dp)*ovl
      end do

      call cphf_solve(mol, scf%orbitals, scf%orbital_energies, scf%n_occupied, &
                      dip, u, err)
      if (.not. err%has_error()) then
         call cphf_solve(mol, scf%orbitals, scf%orbital_energies, scf%n_occupied, &
                         shifted, v, err)
      end if
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "CPHF failed: "//err%get_message())
         return
      end if

      worst = maxval(abs(u - v))
      call check(error, worst < 1.0e-9_dp, &
                 "shifting the dipole origin changed the orbital response, so the "// &
                 "occupied-virtual block is not what it claims to be")

      deallocate (dip, shifted, u, v, ovl)
   end subroutine test_origin_independence

   subroutine test_refusals(error)
      !! Bad input is refused rather than answered
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: dip(:, :, :), u(:, :, :)

      call water(mol, scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         call check(error, .false., "the reference SCF failed")
         return
      end if
      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, err)

      ! No virtual space to rotate into.
      call cphf_solve(mol, scf%orbitals, scf%orbital_energies, &
                      size(scf%orbitals, 2), dip, u, err)
      call check(error, err%has_error(), "CPHF accepted a fully occupied basis")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      err = error_t()

      ! Too few iterations to converge, with a threshold that means it.
      call cphf_solve(mol, scf%orbitals, scf%orbital_energies, scf%n_occupied, &
                      dip, u, err, max_iter=1, tol=1.0e-14_dp)
      call check(error, err%has_error(), &
                 "CPHF returned an unconverged response without saying so")

      call mol%destroy()
      deallocate (dip)
   end subroutine test_refusals

end module test_mqc_libcint_cphf

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_libcint_cphf, only: collect_mqc_libcint_cphf_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_libcint_cphf", collect_mqc_libcint_cphf_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
