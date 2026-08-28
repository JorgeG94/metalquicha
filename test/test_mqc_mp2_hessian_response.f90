!! The first-order skeletons and orbital-response carriers, checked internally
module test_mqc_mp2_hessian_response
   !! Units 1.4-1.6 of the MP2 Hessian ladder (`mp2-hessian-phased-plan.md`):
   !! the per-perturbation skeletons `f^(X)`/`S^(X)`/`<pq|rs>^(X)`, the
   !! skeleton Lagrangian `I'^(X)` with its carriers, and the orbital-response
   !! term of the pass-2 assembly.
   !!
   !! **Why these checks and not a pycc comparison.** The cross-code gates
   !! (pycc's `fX`/`SX`/`erix`, `Ip`/`Xx`/`I2x`, and the isolated `resp_orb`,
   !! symmetric and asymmetric geometry) were run when these units landed and
   !! their residuals are recorded in the commits; they are deliberately not
   !! stored as tests, for the reason `test_mqc_hess_ints`' header gives -- a
   !! stored comparison against an external library pins our layout to that
   !! library's conventions rather than to anything true. What stays in-tree
   !! is what must not drift silently:
   !!
   !! * rigid translation moves no integral, so each skeleton summed over
   !!   atoms cancels -- earned across the four-position permutation assembly
   !!   and the operator-centre terms, not built in;
   !! * the physicist tensor's pair symmetries, which the assembly deposits
   !!   from different integral orderings;
   !! * the transcribed Lagrangian against the gradient's own
   !!   `energy_weighted_ao` -- two of our routines, two decompositions of the
   !!   same matrix;
   !! * the pair-rotation augmentation is the exact identity when nothing is
   !!   frozen and the gauge is non-canonical -- the guard structure pycc's
   !!   Phase 2 branch relies on;
   !! * the orbital-response term inherits the translational zero through the
   !!   CPHF solve.
   !!
   !! One thread throughout, as everywhere on this ladder.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_mp2_gradient, only: libcint_mp2_gradient, build_amplitudes
   use mqc_libcint_mp2, only: transform_ovov
   use mqc_libcint_mp2_hessian, only: mp2_first_order_skeletons, mp2_mo_eri_physicist, &
                                      mp2_cumulant_2pdm, mp2_mo_lagrangian, &
                                      mp2_skeleton_lagrangian, mp2_pair_rotation_augment, &
                                      mp2_orbital_response_term, mp2_full_u, &
                                      mp2_perturbed_fock, mp2_perturbed_eri, &
                                      mp2_perturbed_t2
   use mqc_libcint_hessian, only: make_h1_atom, overlap_deriv_atom, solve_mo1_batch
   use mqc_libcint_hess_ints, only: eri_ip1_block
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   implicit none
   private

   public :: collect_mqc_mp2_hessian_response_tests

   !> The ladder's pinned case, deliberately at the **asymmetric** geometry:
   !> pycc record a ket-swap bug that the symmetric one masked, and nothing
   !> here depends on the residual planarity.
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.000000_dp, 0.000000_dp, &
                           0.031_dp, -0.024_dp, 1.814137_dp, &
                           0.0_dp, 1.756000_dp, -0.454300_dp], [3, 3])
   integer, parameter :: WATER_NELEC = 10

contains

   subroutine collect_mqc_mp2_hessian_response_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("first_order_skeletons_translate_to_nothing", skeletons_translate), &
                  new_unittest("derivative_eri_keeps_its_pair_symmetries", eri_symmetries), &
                  new_unittest("lagrangian_agrees_with_the_gradient", lagrangian_agrees), &
                  new_unittest("pair_augmentation_is_the_allelectron_identity", augment_identity), &
                  new_unittest("orbital_response_term_translates_to_nothing", response_translates), &
                  new_unittest("perturbed_amplitudes_translate_to_nothing", amplitudes_translate) &
                  ]
   end subroutine collect_mqc_mp2_hessian_response_tests

   !> One SCF and the stacked first-order skeletons, shared by every test.
   subroutine stage_at(mol, scf, fx, sx, erix, err)
      type(libcint_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      real(dp), allocatable, intent(out) :: fx(:, :, :), sx(:, :, :)
      real(dp), allocatable, intent(out) :: erix(:, :, :, :, :)
      type(error_t), intent(inout) :: err

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "6-31g", mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, .false., &
                           scf, err)
      if (err%has_error()) return
      call mp2_first_order_skeletons(mol, scf%orbitals, WATER_NELEC/2, &
                                     fx, sx, erix, err)
   end subroutine stage_at

   !> Rigid translation moves no integral, so every skeleton summed over the
   !> atoms of one Cartesian component cancels. Nothing imposes this: the
   !> per-atom pieces mix basis-centre and operator-centre derivatives, and
   !> the two-electron assembly deposits each atom from different index
   !> permutations. Elements reach ~30 (the core Hamiltonian), so 1e-12 is
   !> thirteen digits of cancellation.
   subroutine skeletons_translate(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: fx(:, :, :), sx(:, :, :), erix(:, :, :, :, :)
      real(dp) :: worst_f, worst_s, worst_e
      integer :: threads, comp, a, x, natm

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call stage_at(mol, scf, fx, sx, erix, err)
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), &
                 "the skeletons did not evaluate: "//err%get_message())
      if (allocated(error)) return

      natm = mol%natm
      worst_f = 0.0_dp
      worst_s = 0.0_dp
      worst_e = 0.0_dp
      do comp = 1, 3
         block
            real(dp), allocatable :: sum_f(:, :), sum_s(:, :), sum_e(:, :, :, :)
            allocate (sum_f, sum_s, mold=fx(:, :, 1))
            allocate (sum_e, mold=erix(:, :, :, :, 1))
            sum_f = 0.0_dp
            sum_s = 0.0_dp
            sum_e = 0.0_dp
            do a = 1, natm
               x = 3*(a - 1) + comp
               sum_f = sum_f + fx(:, :, x)
               sum_s = sum_s + sx(:, :, x)
               sum_e = sum_e + erix(:, :, :, :, x)
            end do
            worst_f = max(worst_f, maxval(abs(sum_f)))
            worst_s = max(worst_s, maxval(abs(sum_s)))
            worst_e = max(worst_e, maxval(abs(sum_e)))
         end block
      end do
      write (*, "(a, 3es10.2)") "        translation residuals f/S/eri = ", &
         worst_f, worst_s, worst_e
      call check(error, max(worst_f, worst_s, worst_e) < 1.0e-12_dp, &
                 "a first-order skeleton does not cancel under rigid translation")
      call mol%destroy()
   end subroutine skeletons_translate

   !> `<pq|rs>^(X) = <qp|sr>^(X) = <rq|ps>^(X)`: the derivative inherits the
   !> integral's permutational symmetry, but the assembly does not -- each
   !> relation pairs elements deposited from different `ip1` orderings and
   !> different owner tests, so this checks the permutation map.
   subroutine eri_symmetries(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: fx(:, :, :), sx(:, :, :), erix(:, :, :, :, :)
      real(dp) :: worst_x, worst_b, worst_h
      integer :: threads, x, p, q, r, s, n_mo

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call stage_at(mol, scf, fx, sx, erix, err)
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), &
                 "the skeletons did not evaluate: "//err%get_message())
      if (allocated(error)) return

      n_mo = size(fx, 1)
      worst_x = 0.0_dp
      worst_b = 0.0_dp
      worst_h = 0.0_dp
      do x = 1, size(erix, 5)
         do s = 1, n_mo
            do r = 1, n_mo
               do q = 1, n_mo
                  do p = 1, n_mo
                     worst_x = max(worst_x, abs(erix(p, q, r, s, x) &
                                                - erix(q, p, s, r, x)))
                     worst_b = max(worst_b, abs(erix(p, q, r, s, x) &
                                                - erix(r, q, p, s, x)))
                  end do
               end do
            end do
         end do
         do q = 1, n_mo
            do p = 1, n_mo
               worst_h = max(worst_h, abs(fx(p, q, x) - fx(q, p, x)), &
                             abs(sx(p, q, x) - sx(q, p, x)))
            end do
         end do
      end do
      write (*, "(a, 3es10.2)") "        symmetry residuals swap/bra/1e = ", &
         worst_x, worst_b, worst_h
      call check(error, max(worst_x, worst_b, worst_h) < 1.0e-13_dp, &
                 "a first-derivative tensor lost a permutational symmetry")
      call mol%destroy()
   end subroutine eri_symmetries

   !> The transcribed generalized-Fock Lagrangian against the gradient's
   !> `energy_weighted_ao`: `C I' C^T` must equal the correlation share of the
   !> matrix multiplying `dS/dR`, which is `(W_total + 2 W_ref) / 2` with
   !> `W_ref = 2 sum_i eps_i C_i C_i^T` (conventions note, s.4a). Two routines
   !> that never share a line, one matrix.
   subroutine lagrangian_agrees(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: gradient(:, :), dm1mo(:, :), w_ao(:, :)
      real(dp), allocatable :: eri_packed(:, :), ovov(:, :, :, :), t2(:, :, :, :)
      real(dp), allocatable :: eri_mo(:, :, :, :), gam(:, :, :, :), imat(:, :)
      real(dp), allocatable :: wref(:, :), wcorr(:, :), from_imat(:, :)
      real(dp) :: worst
      integer :: threads, n_ao, n_mo, n_o, i

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "6-31g", mol, err)
      call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, .false., &
                           scf, err)
      n_o = WATER_NELEC/2
      call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, n_o, &
                                gradient, err, n_frozen=0, relaxed_density_mo=dm1mo, &
                                energy_weighted_ao=w_ao)
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), &
                 "the gradient did not evaluate: "//err%get_message())
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      n_ao = mol%nao
      n_mo = size(scf%orbitals, 2)
      call mol%eris_packed(eri_packed)
      call transform_ovov(eri_packed, scf%orbitals, 0, n_o, n_ao, n_mo, ovov)
      call build_amplitudes(ovov, scf%orbital_energies, 0, n_o, n_mo - n_o, n_o, t2)
      call mp2_mo_eri_physicist(mol, scf%orbitals, eri_mo, err)
      call mp2_cumulant_2pdm(t2, 0, n_o, n_mo, gam)
      call mp2_mo_lagrangian(eri_mo, scf%orbital_energies, dm1mo, gam, n_o, imat)

      allocate (wref(n_ao, n_ao))
      wref = 0.0_dp
      do i = 1, n_o
         wref = wref + 2.0_dp*scf%orbital_energies(i) &
                *matmul(reshape(scf%orbitals(:, i), [n_ao, 1]), &
                        reshape(scf%orbitals(:, i), [1, n_ao]))
      end do
      allocate (wcorr(n_ao, n_ao))
      wcorr = 0.5_dp*(w_ao + 2.0_dp*wref)
      allocate (from_imat(n_ao, n_ao))
      from_imat = matmul(scf%orbitals, matmul(imat, transpose(scf%orbitals)))

      worst = maxval(abs(from_imat - wcorr))
      write (*, "(a, es10.2)") "        max |C I' C^T - W_corr| = ", worst
      ! Tight, not solver-limited: both sides consume the same gradient
      ! call's Z-vector, so the identity is arithmetic. Measured 3.2e-14.
      call check(error, worst < 1.0e-12_dp, &
                 "the transcribed Lagrangian and the gradient's energy-weighted "// &
                 "density disagree")
      call mol%destroy()
   end subroutine lagrangian_agrees

   !> All-electron and non-canonical, the pair-rotation augmentation must be
   !> the exact identity -- both rotation classes are guarded off, so the
   !> carriers copy through untouched and `P^(X)` is zero. This pins the guard
   !> structure Phase 2's frozen core will open, bitwise.
   subroutine augment_identity(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: gradient(:, :), dm1mo(:, :), w_ao(:, :)
      real(dp), allocatable :: eri_packed(:, :), ovov(:, :, :, :), t2(:, :, :, :)
      real(dp), allocatable :: fx(:, :, :), sx(:, :, :), erix(:, :, :, :, :)
      real(dp), allocatable :: eri_mo(:, :, :, :), gam(:, :, :, :), imat(:, :)
      real(dp), allocatable :: l_mo(:, :, :, :), ip(:, :), xov(:, :), i2(:, :)
      real(dp), allocatable :: xt(:, :), it(:, :), pf(:, :)
      real(dp) :: worst
      integer :: threads, n_ao, n_mo, n_o, x

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call stage_at(mol, scf, fx, sx, erix, err)
      n_o = WATER_NELEC/2
      call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, n_o, &
                                gradient, err, n_frozen=0, relaxed_density_mo=dm1mo, &
                                energy_weighted_ao=w_ao)
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), &
                 "the setup did not evaluate: "//err%get_message())
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      n_ao = mol%nao
      n_mo = size(scf%orbitals, 2)
      call mol%eris_packed(eri_packed)
      call transform_ovov(eri_packed, scf%orbitals, 0, n_o, n_ao, n_mo, ovov)
      call build_amplitudes(ovov, scf%orbital_energies, 0, n_o, n_mo - n_o, n_o, t2)
      call mp2_mo_eri_physicist(mol, scf%orbitals, eri_mo, err)
      call mp2_cumulant_2pdm(t2, 0, n_o, n_mo, gam)
      call mp2_mo_lagrangian(eri_mo, scf%orbital_energies, dm1mo, gam, n_o, imat)

      allocate (l_mo(n_mo, n_mo, n_mo, n_mo))
      block
         integer :: p, q, r, s
         do s = 1, n_mo
            do r = 1, n_mo
               do q = 1, n_mo
                  do p = 1, n_mo
                     l_mo(p, q, r, s) = 2.0_dp*eri_mo(p, q, r, s) - eri_mo(p, q, s, r)
                  end do
               end do
            end do
         end do
      end block

      worst = 0.0_dp
      do x = 1, size(fx, 3)
         call mp2_skeleton_lagrangian(fx(:, :, x), sx(:, :, x), erix(:, :, :, :, x), &
                                      dm1mo, gam, imat, n_o, ip, xov, i2)
         call mp2_pair_rotation_augment(ip, xov, i2, l_mo, scf%orbital_energies, &
                                        n_o, 0, .false., xt, it, pf)
         worst = max(worst, maxval(abs(xt - xov)), maxval(abs(it - i2)), &
                     maxval(abs(pf)))
         ! The driver is antisymmetric by construction; its occupied-virtual
         ! block is the Z-vector RHS every response consumer reads.
         worst = max(worst, maxval(abs(xov + transpose(xov))))
         deallocate (ip, xov, i2, xt, it, pf)
      end do
      write (*, "(a, es10.2)") "        max augmentation residual = ", worst
      call check(error, worst == 0.0_dp, &
                 "the all-electron non-canonical augmentation is not the identity")
      call mol%destroy()
   end subroutine augment_identity

   !> `2 U^Y X^(X) + S^(Y) I''^(X)` summed over the atoms of either
   !> perturbation cancels: translating every atom leaves the orbitals and the
   !> overlap alone, so both `U^Y` and the skeleton carriers sum to zero, and
   !> the term inherits it -- through the CPHF solve on the `Y` side, which is
   !> why the bound is the solver's, not machine epsilon.
   subroutine response_translates(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: gradient(:, :), dm1mo(:, :), w_ao(:, :)
      real(dp), allocatable :: eri_packed(:, :), ovov(:, :, :, :), t2(:, :, :, :)
      real(dp), allocatable :: fx(:, :, :), sx(:, :, :), erix(:, :, :, :, :)
      real(dp), allocatable :: eri_mo(:, :, :, :), gam(:, :, :, :), imat(:, :)
      real(dp), allocatable :: ip(:, :), xov(:, :), i2(:, :)
      real(dp), allocatable :: xov_st(:, :, :), i2_st(:, :, :)
      real(dp), allocatable :: eip1(:, :, :, :, :), h1(:, :, :, :), s1(:, :, :, :)
      real(dp), allocatable :: h1a(:, :, :), s1a(:, :, :), mo1(:, :, :, :)
      real(dp), allocatable :: orb(:, :)
      real(dp) :: worst, acc
      integer :: threads, n_ao, n_mo, n_o, natm, x, a, comp, iy, ix

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call stage_at(mol, scf, fx, sx, erix, err)
      n_o = WATER_NELEC/2
      call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, n_o, &
                                gradient, err, n_frozen=0, relaxed_density_mo=dm1mo, &
                                energy_weighted_ao=w_ao)
      if (.not. err%has_error()) then
         n_ao = mol%nao
         n_mo = size(scf%orbitals, 2)
         natm = mol%natm
         call mol%eris_packed(eri_packed)
         call transform_ovov(eri_packed, scf%orbitals, 0, n_o, n_ao, n_mo, ovov)
         call build_amplitudes(ovov, scf%orbital_energies, 0, n_o, n_mo - n_o, &
                               n_o, t2)
         call mp2_mo_eri_physicist(mol, scf%orbitals, eri_mo, err)
         call mp2_cumulant_2pdm(t2, 0, n_o, n_mo, gam)
         call mp2_mo_lagrangian(eri_mo, scf%orbital_energies, dm1mo, gam, n_o, imat)

         allocate (xov_st(n_mo, n_mo, 3*natm), i2_st(n_mo, n_mo, 3*natm))
         do x = 1, 3*natm
            call mp2_skeleton_lagrangian(fx(:, :, x), sx(:, :, x), &
                                         erix(:, :, :, :, x), dm1mo, gam, imat, &
                                         n_o, ip, xov, i2)
            xov_st(:, :, x) = xov
            i2_st(:, :, x) = i2
            deallocate (ip, xov, i2)
         end do

         call eri_ip1_block(mol, eip1, err)
         allocate (h1(n_ao, n_ao, 3, natm), s1(n_ao, n_ao, 3, natm))
         do a = 1, natm
            call make_h1_atom(mol, scf%density, eip1, a, h1a, err)
            call overlap_deriv_atom(mol, a, s1a, err)
            h1(:, :, :, a) = h1a
            s1(:, :, :, a) = s1a
            deallocate (h1a, s1a)
         end do
         call solve_mo1_batch(mol, scf%orbitals, scf%orbital_energies, n_o, h1, &
                              s1, mo1, err, max_iter=200, tol=1.0e-13_dp)
      end if
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), &
                 "the response setup did not evaluate: "//err%get_message())
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call mp2_orbital_response_term(mo1, sx, fx, xov_st, i2_st, orb)

      worst = 0.0_dp
      do comp = 1, 3
         do ix = 1, 3*natm
            acc = 0.0_dp
            do a = 1, natm
               acc = acc + orb(ix, 3*(a - 1) + comp)
            end do
            worst = max(worst, abs(acc))
         end do
         do iy = 1, 3*natm
            acc = 0.0_dp
            do a = 1, natm
               acc = acc + orb(3*(a - 1) + comp, iy)
            end do
            worst = max(worst, abs(acc))
         end do
      end do
      write (*, "(a, es10.2)") "        max translation residual = ", worst
      ! Measured 2.2e-16: the translated right-hand side is already zero, so
      ! the solve never gets to spend its tolerance. Room left for neither.
      call check(error, worst < 1.0e-12_dp, &
                 "the orbital-response term does not cancel under rigid translation")
      call mol%destroy()
   end subroutine response_translates

   !> Rigid translation leaves the orbitals, and so the amplitudes, exactly
   !> where they were -- so the perturbed amplitudes summed over the atoms of
   !> one Cartesian component cancel. This is earned across the whole of Unit
   !> 1.7: the skeleton derivatives, the full `U^Y` rotation of the `nmo^4`
   !> integrals, the perturbed Fock with its response fold, and the closed-form
   !> divide all have to agree about their conventions for the sum to vanish,
   !> and it runs through the CPHF solve on every perturbation, which is why
   !> the bound is the solver's rather than machine epsilon.
   !>
   !> What this cannot see is a gauge error -- a perturbed quantity wrong by an
   !> orbital rotation translates to zero just as well (the plan's Unit 1.7
   !> warning). The cross-code gate against pycc's `dt2` dump was run when the
   !> unit landed (sym 6.6e-12, asym 3.4e-11, one thread) and lives in the
   !> commit message, per `test_mqc_hess_ints`' rule on external comparisons.
   subroutine amplitudes_translate(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: fx(:, :, :), sx(:, :, :), erix(:, :, :, :, :)
      real(dp), allocatable :: eri_packed(:, :), ovov(:, :, :, :), t2(:, :, :, :)
      real(dp), allocatable :: eri_mo(:, :, :, :), l_mo(:, :, :, :)
      real(dp), allocatable :: eip1(:, :, :, :, :), h1(:, :, :, :), s1(:, :, :, :)
      real(dp), allocatable :: h1a(:, :, :), s1a(:, :, :), mo1(:, :, :, :)
      real(dp), allocatable :: u(:, :), df(:, :), deri(:, :, :, :), ta(:, :, :, :)
      real(dp), allocatable :: sum_t(:, :, :, :, :)
      real(dp) :: worst
      integer :: threads, n_ao, n_mo, n_o, natm, a, x, comp, p, q, r, s

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call stage_at(mol, scf, fx, sx, erix, err)
      n_o = WATER_NELEC/2
      if (.not. err%has_error()) then
         n_ao = mol%nao
         n_mo = size(scf%orbitals, 2)
         natm = mol%natm
         call mol%eris_packed(eri_packed)
         call transform_ovov(eri_packed, scf%orbitals, 0, n_o, n_ao, n_mo, ovov)
         call build_amplitudes(ovov, scf%orbital_energies, 0, n_o, n_mo - n_o, &
                               n_o, t2)
         call mp2_mo_eri_physicist(mol, scf%orbitals, eri_mo, err)
         allocate (l_mo(n_mo, n_mo, n_mo, n_mo))
         do s = 1, n_mo
            do r = 1, n_mo
               do q = 1, n_mo
                  do p = 1, n_mo
                     l_mo(p, q, r, s) = 2.0_dp*eri_mo(p, q, r, s) - eri_mo(p, q, s, r)
                  end do
               end do
            end do
         end do
         call eri_ip1_block(mol, eip1, err)
         allocate (h1(n_ao, n_ao, 3, natm), s1(n_ao, n_ao, 3, natm))
         do a = 1, natm
            call make_h1_atom(mol, scf%density, eip1, a, h1a, err)
            call overlap_deriv_atom(mol, a, s1a, err)
            h1(:, :, :, a) = h1a
            s1(:, :, :, a) = s1a
            deallocate (h1a, s1a)
         end do
         call solve_mo1_batch(mol, scf%orbitals, scf%orbital_energies, n_o, h1, &
                              s1, mo1, err, max_iter=200, tol=1.0e-13_dp)
      end if
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), &
                 "the perturbed-amplitude setup did not evaluate: "//err%get_message())
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      allocate (sum_t(n_o, n_o, n_mo - n_o, n_mo - n_o, 3))
      sum_t = 0.0_dp
      do x = 1, 3*natm
         a = (x - 1)/3 + 1
         comp = x - 3*(a - 1)
         call mp2_full_u(mo1(:, :, comp, a), sx(:, :, x), n_o, u)
         call mp2_perturbed_fock(fx(:, :, x), sx(:, :, x), u, l_mo, &
                                 scf%orbital_energies, n_o, df)
         call mp2_perturbed_eri(erix(:, :, :, :, x), u, eri_mo, deri)
         call mp2_perturbed_t2(deri, df, t2, scf%orbital_energies, 0, n_o, ta)
         sum_t(:, :, :, :, comp) = sum_t(:, :, :, :, comp) + ta
         deallocate (u, df, deri, ta)
      end do
      worst = maxval(abs(sum_t))
      write (*, "(a, es10.2)") "        max translation residual = ", worst
      ! Measured 6.0e-16 at one thread: as with the orbital-response term,
      ! the translated right-hand side is already zero, so the solve never
      ! gets to spend its tolerance. Room left for neither.
      call check(error, worst < 1.0e-12_dp, &
                 "the perturbed amplitudes do not cancel under rigid translation")
      call mol%destroy()
   end subroutine amplitudes_translate

end module test_mqc_mp2_hessian_response

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mp2_hessian_response, only: collect_mqc_mp2_hessian_response_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_mp2_hessian_response", &
                               collect_mqc_mp2_hessian_response_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
