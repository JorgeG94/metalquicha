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
   use mqc_libcint_direct, only: build_fock_direct, build_fock_direct_many, &
                                 schwarz_bounds, direct_stats_t
   use mqc_libcint_cphf, only: cphf_solve, static_polarizability, distributed_polarizability, &
                               build_hessian, build_hessian_mo, static_response_dense
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
                               test_refusals), &
                  new_unittest("distributed_tensors_sum_to_the_total", &
                               test_distributed_sum_rule), &
                  new_unittest("distributed_tensors_are_individually_asymmetric", &
                               test_distributed_asymmetry), &
                  new_unittest("direct_and_stored_integrals_agree", test_direct), &
                  new_unittest("a_batch_of_densities_equals_one_at_a_time", &
                               test_density_batch), &
                  new_unittest("the_transformed_hessian_equals_the_probed_one", &
                               test_hessian_transform), &
                  new_unittest("the_dense_static_solve_equals_the_iterative_one", &
                               test_static_dense) &
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

   subroutine test_distributed_sum_rule(error)
      !! The per-orbital tensors add up to the molecule's total polarizability
      !!
      !! Exact by construction rather than approximately true, because the
      !! decomposition only inserts `W W^T` into a sum that was already there. So
      !! this is asserted at solver precision, not at chemical precision -- and
      !! that sharpness is the point: it catches a transposed `W`, a mismatched
      !! occupied index, or a localization that quietly dropped an orbital, none
      !! of which would move the total enough to notice any other way.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: tensors(:, :, :), centroids(:, :)
      real(dp) :: total(3, 3), alpha(3, 3)
      integer :: i, k, l

      call water(mol, scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         call check(error, .false., "the reference SCF failed")
         return
      end if

      call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                 scf%n_occupied, alpha, err)
      if (.not. err%has_error()) then
         call distributed_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                         scf%n_occupied, tensors, centroids, err)
      end if
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "failed: "//err%get_message())
         return
      end if

      call check(error, size(tensors, 3), scf%n_occupied, &
                 "one tensor per occupied orbital was expected with no core excluded")
      if (allocated(error)) return

      total = 0.0_dp
      do i = 1, size(tensors, 3)
         total = total + tensors(:, :, i)
      end do
      do l = 1, 3
         do k = 1, 3
            call check(error, total(k, l), alpha(k, l), &
                       "the distributed tensors do not sum to the total "// &
                       "polarizability", thr=1.0e-10_dp)
            if (allocated(error)) return
         end do
      end do

      ! And with the core excluded the sum must fall *short* of the total, by the
      ! core's own contribution -- small, and diagonal, since a 1s is spherical.
      ! GAMESS's water potential excludes it, so this is the arrangement a
      ! reference comparison actually uses.
      deallocate (tensors, centroids)
      call water(mol, scf, err)
      call distributed_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                      scf%n_occupied, tensors, centroids, err, &
                                      n_core=1)
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "failed with a core: "//err%get_message())
         return
      end if
      call check(error, size(tensors, 3), scf%n_occupied - 1, &
                 "excluding one core orbital did not remove one tensor")
      if (allocated(error)) return

      total = 0.0_dp
      do i = 1, size(tensors, 3)
         total = total + tensors(:, :, i)
      end do

      ! Per component the core can only lower a diagonal or leave it be, so assert
      ! that none *grows*. Not a strict `total < alpha`: water's 1s is spherical
      ! and, for this orientation, decoupled from y, so it contributes nothing to
      ! alpha_yy and that diagonal returns equal to the total to the last bit. A
      ! strict comparison there is deciding a one-ULP difference on the sign of the
      ! summation error, which flips between toolchains -- the flake this once was.
      do k = 1, 3
         call check(error, total(k, k) <= alpha(k, k) + 1.0e-10_dp, &
                    "excluding the core increased a diagonal polarizability")
         if (allocated(error)) return
      end do

      ! The reduction still has to appear somewhere, or n_core did nothing. It
      ! lands in the trace, which falls by the core's whole contribution -- about
      ! 6e-4 here, far above the per-component rounding -- so assert that robustly.
      call check(error, total(1, 1) + total(2, 2) + total(3, 3) < &
                 alpha(1, 1) + alpha(2, 2) + alpha(3, 3), &
                 "excluding the core did not reduce the total polarizability")
      if (allocated(error)) return
      call check(error, abs(total(1, 2) - alpha(1, 2)) < 1.0e-4_dp, &
                 "excluding the core moved an off-diagonal component, but a "// &
                 "spherical core should contribute only to the diagonal")

      deallocate (tensors, centroids)
   end subroutine test_distributed_sum_rule

   subroutine test_distributed_asymmetry(error)
      !! Each tensor is asymmetric; the *complete* sum is symmetric
      !!
      !! The contraction pairs perturbation `k` with response `l`, which for one
      !! orbital are different objects. GAMESS's water reference shows the same
      !! thing -- 0.17 per O-H bond tensor -- so a well-meaning
      !! `0.5*(a + transpose(a))` here would agree with the total and disagree with
      !! every reference file. This test exists to make that temptation fail.
      !!
      !! **Which sum is symmetric matters, and only the complete one is.** Writing
      !! the partial sum over retained orbitals as `4 h^k P A^-1 h^l` with `P` the
      !! projector onto them, `P` and `A^-1` do not commute, so excluding the core
      !! leaves a small genuine asymmetry -- 3e-7 here. It is not an error and
      !! cannot be tightened away: GAMESS's own water potential carries 4.6e-6 in
      !! the same place. Only `P = 1` gives back the exactly symmetric total, so
      !! that is where the sharp assertion goes.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: tensors(:, :, :), centroids(:, :)
      real(dp) :: total(3, 3), worst_single, sum_asym, valence_asym
      integer :: i

      call water(mol, scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         call check(error, .false., "the reference SCF failed")
         return
      end if
      call distributed_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                      scf%n_occupied, tensors, centroids, err)
      if (err%has_error()) then
         call check(error, .false., "failed: "//err%get_message())
         call mol%destroy()
         return
      end if

      total = 0.0_dp
      worst_single = 0.0_dp
      do i = 1, size(tensors, 3)
         total = total + tensors(:, :, i)
         worst_single = max(worst_single, &
                            maxval(abs(tensors(:, :, i) - transpose(tensors(:, :, i)))))
      end do
      sum_asym = maxval(abs(total - transpose(total)))

      call check(error, worst_single > 1.0e-3_dp, &
                 "the per-orbital tensors came back symmetric, so either they "// &
                 "were averaged or the occupied index was contracted wrongly")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      call check(error, sum_asym < 1.0e-10_dp, &
                 "the sum over every occupied orbital is not symmetric, and that "// &
                 "one must be exactly so")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      ! Excluding the core leaves a small asymmetry, by the projector argument
      ! above. Bounded rather than required to vanish -- and required to be *there*,
      ! since a zero would mean the core had been frozen out of the response too,
      ! which is a different convention from GAMESS's and would shift every tensor.
      deallocate (tensors, centroids)
      call distributed_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                      scf%n_occupied, tensors, centroids, err, &
                                      n_core=1)
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "failed with a core: "//err%get_message())
         return
      end if
      total = 0.0_dp
      do i = 1, size(tensors, 3)
         total = total + tensors(:, :, i)
      end do
      valence_asym = maxval(abs(total - transpose(total)))
      call check(error, valence_asym < 1.0e-5_dp, &
                 "the valence-only sum is more asymmetric than a projector "// &
                 "artefact can explain")
      if (allocated(error)) return
      call check(error, valence_asym > 1.0e-12_dp, &
                 "the valence-only sum is exactly symmetric, which means the core "// &
                 "was frozen out of the response rather than only out of the "// &
                 "distribution -- not the convention GAMESS uses")
      if (allocated(error)) return

      ! Every centroid must sit inside the molecule -- a crude bound, but it
      ! catches a centroid array that got out of step with its tensors.
      do i = 1, size(centroids, 2)
         call check(error, norm2(centroids(:, i)) < 5.0_dp, &
                    "a localized orbital centroid is nowhere near the molecule")
         if (allocated(error)) return
      end do

      deallocate (tensors, centroids)
   end subroutine test_distributed_asymmetry

   subroutine test_direct(error)
      !! Recomputing the integrals gives the same response as storing them
      !!
      !! The direct build is what will actually run: the stored tensor is `n_ao^4`,
      !! 2 GB for uracil in 6-31G* where the response equations want tens of MB, so
      !! storage and not the solver is what decides how large a fragment can be
      !! done. The in-core path stays because it is the one checked against a dense
      !! exact solve, which makes it the reference for this comparison -- the same
      !! arrangement `run_libcint_rhf` uses.
      !!
      !! Screening is why this is not a formality. The direct build drops shell
      !! quartets by a Schwarz bound computed from the basis, and a response density
      !! is not a converged SCF density -- it has a different magnitude and a
      !! different sparsity -- so a screening threshold that is invisible in an SCF
      !! can be visible here.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp) :: stored(3, 3), recomputed(3, 3)
      integer :: k, l

      call water(mol, scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         call check(error, .false., "the reference SCF failed")
         return
      end if

      call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                 scf%n_occupied, stored, err, in_core=.true.)
      if (.not. err%has_error()) then
         call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                    scf%n_occupied, recomputed, err, in_core=.false.)
      end if
      call mol%destroy()
      if (err%has_error()) then
         call check(error, .false., "CPHF failed: "//err%get_message())
         return
      end if

      do l = 1, 3
         do k = 1, 3
            call check(error, recomputed(k, l), stored(k, l), &
                       "the direct and stored integral paths disagree", &
                       thr=1.0e-9_dp)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_direct

   subroutine test_density_batch(error)
      !! Contracting many densities in one integral pass matches doing them singly
      !!
      !! The amortization the frequency-dependent equations need: an integral is
      !! computed once and reused across every density in hand. The risk is that
      !! batching changes an answer -- through a screening decision that depended on
      !! the density, or an accumulator indexed with the set offset wrong -- and
      !! both would be invisible in an SCF, which only ever passes one.
      !!
      !! Screening here is a Schwarz bound over the basis and so is density
      !! independent, which is what makes the two routes comparable at all; this test
      !! is what would notice if that ever stopped being true.
      !!
      !! Agreement is to 1e-12 rather than exact because the threads reduce into the
      !! accumulator in whatever order they finish, so the sums are differently
      !! ordered in the two routes -- the same reason the SCF has an OpenMP-order
      !! sensitivity at the last digits.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(direct_stats_t) :: stats
      real(dp), allocatable :: bounds(:, :), zero_h(:, :), dens(:, :, :)
      real(dp), allocatable :: batched(:, :, :), single(:, :)
      real(dp) :: worst
      integer :: n, i, nset

      nset = 4
      call water(mol, scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         call check(error, .false., "the reference SCF failed")
         return
      end if
      n = mol%nao
      call schwarz_bounds(mol, bounds, err)
      if (err%has_error()) then
         call check(error, .false., "Schwarz bounds failed")
         call mol%destroy()
         return
      end if

      allocate (zero_h(n, n), dens(n, n, nset), single(n, n))
      zero_h = 0.0_dp
      ! A converged density, then progressively scaled and shifted copies, so the
      ! sets are genuinely different and none is a multiple of another.
      do i = 1, nset
         dens(:, :, i) = scf%density/real(i, dp)
         dens(1, 1, i) = dens(1, 1, i) + 0.1_dp*real(i, dp)
         dens(:, :, i) = 0.5_dp*(dens(:, :, i) + transpose(dens(:, :, i)))
      end do

      call build_fock_direct_many(mol, zero_h, dens, bounds, batched, stats, err)
      if (err%has_error()) then
         call check(error, .false., "the batched build failed: "//err%get_message())
         call mol%destroy()
         return
      end if

      worst = 0.0_dp
      do i = 1, nset
         call build_fock_direct(mol, zero_h, dens(:, :, i), bounds, single, stats, err)
         if (err%has_error()) then
            call check(error, .false., "a single build failed: "//err%get_message())
            call mol%destroy()
            return
         end if
         worst = max(worst, maxval(abs(single - batched(:, :, i))))
      end do
      call mol%destroy()

      call check(error, worst < 1.0e-12_dp, &
                 "a batch of densities does not reproduce the same densities "// &
                 "contracted one at a time")
      if (allocated(error)) return

      ! And the batch must not have quietly collapsed to one set.
      call check(error, size(batched, 3), nset, "the batch lost a density")

      deallocate (bounds, zero_h, dens, batched, single)
   end subroutine test_density_batch

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

   subroutine test_hessian_transform(error)
      !! The transformed Hessian is the one the column build probes for
      !!
      !! `build_hessian` gets column `j` by applying the operator to the unit
      !! vector `e_j` -- one Fock build per column, `n_ov` of them. `build_hessian_mo`
      !! never forms a column: it transforms the integrals into the two MO blocks
      !! the operators are made of and assembles them. Same operator, and on a
      !! ten-atom fragment the difference between thirty-four seconds and two.
      !!
      !! What makes that substitution safe is only that the two agree, so this is
      !! the test that has to exist for it. Nothing downstream would catch a
      !! disagreement cheaply: the Hessian is reached through a makefp run, where
      !! it reappears as a dynamic polarizability with no independent value to
      !! check against.
      !!
      !! Element by element rather than by any norm, because the failure this
      !! guards against is a permuted index -- `(aj|ib)` read where `(ab|ij)`
      !! belongs -- which moves a few entries a long way and leaves a norm over
      !! four thousand of them looking respectable.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: eri(:, :, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), allocatable :: plus_col(:, :), minus_col(:, :)
      real(dp), allocatable :: plus_mo(:, :), minus_mo(:, :)
      integer :: n_ao, n_mo, n_occ, n_vir, a, i

      call water(mol, scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         call check(error, .false., "the reference SCF failed")
         return
      end if

      n_ao = mol%nao
      n_mo = size(scf%orbitals, 2)
      n_occ = scf%n_occupied
      n_vir = n_mo - n_occ
      allocate (c_occ(n_ao, n_occ), c_vir(n_ao, n_vir), gaps(n_vir, n_occ))
      c_occ = scf%orbitals(:, 1:n_occ)
      c_vir = scf%orbitals(:, n_occ + 1:n_mo)
      do i = 1, n_occ
         do a = 1, n_vir
            gaps(a, i) = scf%orbital_energies(n_occ + a) - scf%orbital_energies(i)
         end do
      end do

      call mol%eris(eri)
      allocate (zero_h(n_ao, n_ao), bounds(0, 0))
      zero_h = 0.0_dp

      ! The column build, on the stored tensor so that neither route involves the
      ! direct build's thread-order reduction and any difference found is the
      ! transformation's own.
      call build_hessian(mol, .false., eri, bounds, zero_h, c_occ, c_vir, gaps, &
                         plus_col, minus_col, 8, err)
      call check(error,.not. err%has_error(), "the column build failed: "//err%get_message())
      if (allocated(error)) return

      call build_hessian_mo(mol, eri, c_occ, c_vir, gaps, plus_mo, minus_mo, err)
      call check(error,.not. err%has_error(), "the transform failed: "//err%get_message())
      if (allocated(error)) return

      call check(error, size(plus_mo, 1) == size(plus_col, 1) .and. &
                 size(plus_mo, 2) == size(plus_col, 2), "the operators are different shapes")
      if (allocated(error)) return

      call check(error, maxval(abs(plus_mo - plus_col)) < 1.0e-10_dp, &
                 "(A+B) disagrees between the transform and the column build")
      if (allocated(error)) return
      call check(error, maxval(abs(minus_mo - minus_col)) < 1.0e-10_dp, &
                 "(A-B) disagrees between the transform and the column build")
      if (allocated(error)) return

      ! `(A-B)` is diagonal-dominant and symmetric, and the transform reads it out
      ! of a permuted index; a transpose slipped in there would still match the
      ! column build only if the column build had the same slip.
      call check(error, maxval(abs(minus_mo - transpose(minus_mo))) < 1.0e-12_dp, &
                 "(A-B) should be symmetric")

      call mol%destroy()
   end subroutine test_hessian_transform

   subroutine test_static_dense(error)
      !! Factorizing (A+B) gives what iterating against it gives
      !!
      !! The static response is `(A+B) U = -h`, and `cphf_solve` reaches it by
      !! conjugate gradients -- a Fock build per iteration per perturbation. Once
      !! the dynamic blocks have built `(A+B)` for their own solves the same
      !! equation is one factorization away, which is forty seconds of a hundred
      !! on a tripeptide.
      !!
      !! Worth a test because the substitution is safe only if the two really are
      !! the same equation, and there is a nearby one that is not: the
      !! zero-frequency member of the dynamic family is reached by multiplying
      !! through by `(A-B)`, and agrees with this only to about 1e-4. Confusing
      !! the two would move every static polarizability in a potential without
      !! anything failing.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: dip(:, :, :), u_iter(:, :, :), u_dense(:, :, :)
      real(dp), allocatable :: eri(:, :, :, :), c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), allocatable :: aplus(:, :), aminus(:, :)
      integer :: n_ao, n_mo, n_occ, n_vir, a, i
      real(dp) :: worst, scale

      call water(mol, scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         call check(error, .false., "the reference SCF failed")
         return
      end if
      n_ao = mol%nao
      n_mo = size(scf%orbitals, 2)
      n_occ = scf%n_occupied
      n_vir = n_mo - n_occ

      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, err)
      call cphf_solve(mol, scf%orbitals, scf%orbital_energies, n_occ, dip, u_iter, &
                      err, tol=1.0e-13_dp, in_core=.true.)
      call check(error, .not. err%has_error(), "the iterative solve failed: "//err%get_message())
      if (allocated(error)) return

      allocate (c_occ(n_ao, n_occ), c_vir(n_ao, n_vir), gaps(n_vir, n_occ))
      c_occ = scf%orbitals(:, 1:n_occ)
      c_vir = scf%orbitals(:, n_occ + 1:n_mo)
      do i = 1, n_occ
         do a = 1, n_vir
            gaps(a, i) = scf%orbital_energies(n_occ + a) - scf%orbital_energies(i)
         end do
      end do
      call mol%eris(eri)
      call build_hessian_mo(mol, eri, c_occ, c_vir, gaps, aplus, aminus, err)
      call check(error, .not. err%has_error(), "the Hessian build failed: "//err%get_message())
      if (allocated(error)) return

      call static_response_dense(aplus, dip, c_occ, c_vir, u_dense, err)
      call check(error, .not. err%has_error(), "the dense solve failed: "//err%get_message())
      if (allocated(error)) return

      worst = maxval(abs(u_iter - u_dense))
      scale = maxval(abs(u_iter))
      call check(error, worst <= 1.0e-9_dp*max(scale, 1.0_dp), &
                 "the dense and iterative static responses disagree")

      call mol%destroy()
   end subroutine test_static_dense

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
