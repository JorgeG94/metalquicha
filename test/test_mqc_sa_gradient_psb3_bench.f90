!! Phase-5 benchmark: SA-CASSCF gradients of every root on PSB3
module test_mqc_sa_gradient_psb3_bench
   !! The timing/profiling harness `SA_CASSCF_GRADIENT_PLAN.md` phase 5 asks
   !! for, before any fusion is designed. PSB3 (penta-2,4-dieniminium cation,
   !! `tools/sa_casscf/psb3.xyz`), SA-N-CAS(6,6)/6-31G* (Cartesian), N = 2, 3, 4.
   !!
   !! **Active-space selection.** `run_czt_casscf` takes the active window as
   !! orbitals `n_inactive+1 .. n_inactive+n_active` of whatever is handed to
   !! it, so the six pi orbitals have to be moved there by hand: RHF's
   !! canonical ordering interleaves pi and sigma orbitals near the frontier,
   !! it does not hand them over contiguously. PSB3 is built planar with every
   !! atom at `z = 0` (`tools/sa_casscf/psb3.xyz`), so reflection through the
   !! molecular plane (`z -> -z`) is an exact symmetry of the molecule and the
   !! basis, and it partitions every Cartesian AO into an even class (`s`,
   !! `px`, `py`, `dxx`, `dxy`, `dyy`, `dzz`, ... -- an even power of `z` in
   !! `x^(l-i) y^(i-j) z^j`, libcint's own Cartesian component order,
   !! `mqc_czt_ao.f90`) and an odd one (`pz`, `dxz`, `dyz`, ... an odd power).
   !! `S` and `F` are then exactly block-diagonal between the two classes, so a
   !! converged closed-shell RHF's own orbitals are pure symmetry eigenvectors:
   !! an orbital's AO coefficients are either entirely on the odd class (a pi
   !! orbital) or entirely on the even one (sigma/lone-pair/whatever), with no
   !! metric needed to tell them apart -- `pi_weight`, the raw sum of squared
   !! coefficients on the odd AOs, comes back at machine-precision 0 or an
   !! O(1) number, never in between. `classify_pi_orbitals` builds that odd/even
   !! split straight from `mol%bas`/`shell_offset` (no PySCF cross-reference
   !! needed: the Cartesian component order is this code's own convention,
   !! stated in `mqc_czt_ao.f90`), then `select_active_space` picks the three
   !! occupied and three virtual pi orbitals closest to the Fermi level and
   !! reorders the orbital columns so they land contiguously at
   !! `n_inactive+1 .. n_inactive+6`.
   use, intrinsic :: iso_fortran_env, only: int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use pic_io, only: to_char
   use mqc_physical_constants, only: BOHR_TO_ANGSTROM
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, shell_dim
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mcscf, only: run_czt_casscf, casscf_result_t
   use mqc_czt_sa_gradient, only: czt_sa_casscf_gradient, czt_sa_casscf_gradients
   use mqc_czt_sa_hessian, only: sa_hessian_t, build_sa_hessian, destroy_sa_hessian, &
                                 sa_hessian_n_param
   use libcint_fortran, only: LIBCINT_ANG_OF
   implicit none
   private

   public :: collect_mqc_sa_gradient_psb3_bench_tests

   integer, parameter :: N_ATOMS = 14
   integer, parameter :: PSB3_Z(N_ATOMS) = [6, 6, 6, 6, 6, 7, 1, 1, 1, 1, 1, 1, 1, 1]
   character(len=2), parameter :: PSB3_SYM(N_ATOMS) = ["C ", "C ", "C ", "C ", "C ", "N ", &
                                                       "H ", "H ", "H ", "H ", "H ", "H ", &
                                                       "H ", "H "]
   ! tools/sa_casscf/psb3.xyz, Angstrom -- RHF/6-31G* optimum (PySCF geometric).
   real(dp), parameter :: PSB3_XYZ(3, N_ATOMS) = reshape([ &
                                                         0.00370489_dp, -0.05103774_dp, -0.00000000_dp, &
                                                         1.11900869_dp, 0.67603172_dp, -0.00000000_dp, &
                                                         2.41564022_dp, 0.04610394_dp, -0.00000000_dp, &
                                                         3.59754549_dp, 0.71053179_dp, -0.00000000_dp, &
                                                         4.80239601_dp, -0.02131695_dp, -0.00000000_dp, &
                                                         5.98559633_dp, 0.50262760_dp, 0.00000000_dp, &
                                                         -0.96513394_dp, 0.41185645_dp, -0.00000000_dp, &
                                                         0.02310333_dp, -1.12642786_dp, 0.00000000_dp, &
                                                         1.07780722_dp, 1.75013103_dp, -0.00000000_dp, &
                                                         2.42630449_dp, -1.03175504_dp, -0.00000000_dp, &
                                                         3.63381571_dp, 1.78500997_dp, -0.00000000_dp, &
                                                         4.76799191_dp, -1.09647567_dp, -0.00000000_dp, &
                                                         6.80855581_dp, -0.06379908_dp, -0.00000000_dp, &
                                                         6.12403059_dp, 1.49351982_dp, 0.00000000_dp], [3, N_ATOMS])

   integer, parameter :: N_ELECTRONS = 44   !! sum(Z) - charge(+1)
   integer, parameter :: N_OCC = 22
   integer, parameter :: N_ACTIVE = 6
   integer, parameter :: N_INACTIVE = N_OCC - 3   !! 3 of the 22 occupied MOs are pi
   integer, parameter :: N_ALPHA = 3, N_BETA = 3
   real(dp), parameter :: PI_THRESHOLD = 1.0e-6_dp
   integer, parameter :: N_REPEAT = 3   !! Timing repeats; report the min
   integer, parameter :: MACRO_ITERS = 2
      !! Macro-iterations per SA solve. `run_czt_casscf` builds the explicit
      !! orbital Hessian every macro-iteration, which for PSB3/6-31G* (about
      !! 2100 rotations) costs minutes; the gradient's cost does not depend on
      !! how converged the orbitals are, so the timing point is a capped run.

contains

   subroutine collect_mqc_sa_gradient_psb3_bench_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("psb3_active_space_is_the_pi_system", test_active_space), &
                  new_unittest("psb3_sa_gradient_timing", test_timing) &
                  ]
   end subroutine collect_mqc_sa_gradient_psb3_bench_tests

   subroutine build_psb3(mol, error)
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: error

      call build_czt_molecule(PSB3_Z, PSB3_SYM, PSB3_XYZ/BOHR_TO_ANGSTROM, "6-31g_st_", &
                              mol, error)
   end subroutine build_psb3

   subroutine classify_pi_orbitals(mol, orbitals, pi_weight)
      !! `pi_weight(k)`: the raw sum of squared AO coefficients of orbital `k`
      !! over every AO whose Cartesian exponent carries an odd power of `z` --
      !! see the module docstring. `PSB3_XYZ` has `z = 0` for every atom, which
      !! is what makes this an exact (not approximate) classification.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), allocatable, intent(out) :: pi_weight(:)

      logical, allocatable :: odd(:)
      integer :: ish, l, d, off, comp, i, j, p, n_mo, k

      allocate (odd(mol%nao))
      odd = .false.
      do ish = 1, mol%nbas
         l = mol%bas(LIBCINT_ANG_OF, ish)
         d = shell_dim(mol%cartesian, ish - 1, mol%bas)
         off = mol%shell_offset(ish)
         comp = 0
         do i = 0, l
            do j = 0, i
               comp = comp + 1
               p = off + comp
               ! x^(l-i) y^(i-j) z^j: the z power is j.
               odd(p) = (mod(j, 2) == 1)
            end do
         end do
      end do

      n_mo = size(orbitals, 2)
      allocate (pi_weight(n_mo))
      do k = 1, n_mo
         pi_weight(k) = sum(orbitals(:, k)**2, mask=odd)
      end do
      deallocate (odd)
   end subroutine classify_pi_orbitals

   subroutine select_active_space(mol, orbitals_in, orbitals_out, n_pi_occ, n_pi_virt, error)
      !! Reorders columns so the three occupied and three virtual pi orbitals
      !! (by `classify_pi_orbitals`) land at `N_INACTIVE+1 .. N_INACTIVE+6`,
      !! occupied ones first, everything else keeping its relative order.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals_in(:, :)
      real(dp), allocatable, intent(out) :: orbitals_out(:, :)
      integer, intent(out) :: n_pi_occ, n_pi_virt
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: pi_weight(:)
      logical, allocatable :: is_pi(:)
      integer, allocatable :: order(:)
      integer :: n_mo, k, n_inact_fill, n_act_fill, n_virt_fill

      call classify_pi_orbitals(mol, orbitals_in, pi_weight)
      n_mo = size(orbitals_in, 2)
      allocate (is_pi(n_mo))
      is_pi = pi_weight > PI_THRESHOLD

      ! Every d shell's `dxz`/`dyz` component is also odd under the mirror, so
      ! the odd-AO subspace is bigger than the six pi/pi* valence combinations
      ! -- it also spans a batch of much higher-lying virtual combinations
      ! built mostly from polarization functions. `n_pi_occ` still comes back
      ! at exactly 3 (nothing occupied is built from a d function alone), but
      ! counting every `is_pi` virtual overcounts; the chemical pi* triad is
      ! the three *lowest* of them; canonical RHF orbitals are energy-ordered,
      ! so "lowest" is "first in the virtual block".
      n_pi_occ = count(is_pi(1:N_OCC))
      n_pi_virt = min(3, count(is_pi(N_OCC + 1:n_mo)))
      if (n_pi_occ /= 3 .or. count(is_pi(N_OCC + 1:n_mo)) < 3) then
         call error%set(ERROR_VALIDATION, "psb3 bench: expected 3 occupied and at least 3 "// &
                        "virtual pi orbitals, found "//to_char(n_pi_occ)//" and "// &
                        to_char(count(is_pi(N_OCC + 1:n_mo))))
         return
      end if
      n_pi_virt = 3

      allocate (order(n_mo))
      n_inact_fill = 0
      n_act_fill = 0
      n_virt_fill = 0
      ! Occupied pi orbitals first (active window's first half), then the
      ! three lowest-energy virtual pi orbitals (second half), matching how
      ! `active_space_rdms`/`run_czt_casscf` read the active window -- order
      ! among the six does not change the physics (any orthonormal active
      ! basis spans the same space and the CASSCF CI solve re-diagonalises
      ! it), only readability.
      do k = 1, N_OCC
         if (is_pi(k)) then
            order(N_INACTIVE + n_act_fill + 1) = k
            n_act_fill = n_act_fill + 1
         else
            n_inact_fill = n_inact_fill + 1
            order(n_inact_fill) = k
         end if
      end do
      do k = N_OCC + 1, n_mo
         if (is_pi(k) .and. n_virt_fill < 3) then
            order(N_INACTIVE + n_act_fill + 1) = k
            n_act_fill = n_act_fill + 1
            n_virt_fill = n_virt_fill + 1
         else
            order(N_INACTIVE + N_ACTIVE + (k - N_OCC - n_virt_fill)) = k
         end if
      end do

      allocate (orbitals_out(size(orbitals_in, 1), n_mo))
      do k = 1, n_mo
         orbitals_out(:, k) = orbitals_in(:, order(k))
      end do
      deallocate (pi_weight, is_pi, order)
   end subroutine select_active_space

   subroutine test_active_space(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: active_orbitals(:, :)
      integer :: n_pi_occ, n_pi_virt

      call build_psb3(mol, err)
      call check(error,.not. err%has_error(), "psb3 molecule should build")
      if (allocated(error)) return

      call run_czt_rhf(mol, N_ELECTRONS, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "psb3 RHF should converge")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      print '(a,i0,a,f18.10)', "# psb3: nao = ", mol%nao, ", RHF energy = ", scf%energy

      block
         real(dp), allocatable :: pi_weight(:)
         integer :: k
         call classify_pi_orbitals(mol, scf%orbitals, pi_weight)
         do k = 15, 30
            print '(a,i0,a,es14.4)', "# orbital ", k, " pi_weight = ", pi_weight(k)
         end do
      end block

      call select_active_space(mol, scf%orbitals, active_orbitals, n_pi_occ, n_pi_virt, err)
      call check(error,.not. err%has_error(), "psb3: the active space should be the pi system")
      call check(error, n_pi_occ, 3)
      call check(error, n_pi_virt, 3)

      deallocate (active_orbitals)
      call mol%destroy()
   end subroutine test_active_space

   subroutine converge_sa(mol, guess, n_states, weights, result, error)
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: guess(:, :)
      integer, intent(in) :: n_states
      real(dp), intent(in) :: weights(:)
      type(casscf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error

      integer(int64) :: c0, c1, rate

      call system_clock(c0, rate)
      call run_czt_casscf(mol, guess, N_INACTIVE, N_ACTIVE, N_ALPHA, N_BETA, result, error, &
                          max_iterations=MACRO_ITERS, gradient_tol=1.0e-4_dp, &
                          n_states=n_states, weights=weights)
      call system_clock(c1)
      print '(a,i0,a,i0,a,f10.2,a,es10.2)', "# SA-", n_states, ": ", result%iterations, &
         " macro-iterations in ", real(c1 - c0, dp)/real(rate, dp), &
         " s, orbital gradient ", result%gradient_norm
   end subroutine converge_sa

   subroutine time_fused(mol, result, n_roots, best_seconds, error)
      !! One `czt_sa_casscf_gradients` call for roots `1..n_roots`, `N_REPEAT`
      !! times, reporting the fastest
      type(czt_molecule_t), intent(in) :: mol
      type(casscf_result_t), intent(in) :: result
      integer, intent(in) :: n_roots
      real(dp), intent(out) :: best_seconds
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: gradients(:, :, :), weights(:)
      integer(int64) :: c0, c1, rate
      integer :: rep, root, n_states

      n_states = size(result%energies)
      allocate (weights(n_states))
      weights = 1.0_dp/real(n_states, dp)

      best_seconds = huge(1.0_dp)
      do rep = 1, N_REPEAT
         call system_clock(c0, rate)
         call czt_sa_casscf_gradients(mol, result%orbitals, N_INACTIVE, N_ACTIVE, N_ALPHA, &
                                      N_BETA, result%ci_vectors, result%energies, weights, &
                                      [(root, root=1, n_roots)], gradients, error)
         if (error%has_error()) return
         call system_clock(c1)
         best_seconds = min(best_seconds, real(c1 - c0, dp)/real(rate, dp))
      end do
   end subroutine time_fused

   subroutine time_unfused(mol, result, n_roots, best_seconds, error)
      !! `n_roots` separate `czt_sa_casscf_gradient` calls -- today's phase-4
      !! path -- `N_REPEAT` times, reporting the fastest
      type(czt_molecule_t), intent(in) :: mol
      type(casscf_result_t), intent(in) :: result
      integer, intent(in) :: n_roots
      real(dp), intent(out) :: best_seconds
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: gradient(:, :), weights(:)
      integer(int64) :: c0, c1, rate
      real(dp) :: elapsed
      integer :: rep, root, n_states

      n_states = size(result%energies)
      allocate (weights(n_states))
      weights = 1.0_dp/real(n_states, dp)

      best_seconds = huge(1.0_dp)
      do rep = 1, N_REPEAT
         call system_clock(c0, rate)
         do root = 1, n_roots
            call czt_sa_casscf_gradient(mol, result%orbitals, N_INACTIVE, N_ACTIVE, N_ALPHA, &
                                        N_BETA, result%ci_vectors, result%energies, weights, &
                                        root, gradient, error)
            if (error%has_error()) return
         end do
         call system_clock(c1)
         elapsed = real(c1 - c0, dp)/real(rate, dp)
         best_seconds = min(best_seconds, elapsed)
      end do
   end subroutine time_unfused

   subroutine time_breakdown(mol, result, error)
      !! One root's pieces, timed individually with the same public routines
      !! `sa_casscf_gradient_general` calls -- the profile Step 1 asks for,
      !! before any fusion. `build_sa_hessian` alone (paid once by the fused
      !! path per N roots, N times by the unfused one) is the headline number.
      use pic_blas_interfaces, only: pic_gemm
      use mqc_czt_mcscf, only: mcscf_fock_t, generalized_fock, orbital_gradient
      use mqc_czt_mcscf_gradient, only: czt_mcscf_gradient
      use mqc_czt_sa_gradient, only: sa_zvector_solve, orbital_response_gradient, &
                                     ci_response_gradient
      use mqc_rdm, only: active_space_rdms
      type(czt_molecule_t), intent(in) :: mol
      type(casscf_result_t), intent(in) :: result
      type(error_t), intent(inout) :: error

      type(sa_hessian_t) :: state
      real(dp), allocatable :: weights(:), dm1_i(:, :), dm2_i(:, :, :, :)
      real(dp), allocatable :: grad_full(:, :), rhs(:, :), x(:, :)
      real(dp), allocatable :: kappa_bar(:, :), xbar(:, :, :)
      real(dp), allocatable :: base_g(:, :), orb_g(:, :), ci_g(:, :)
      type(mcscf_fock_t) :: fock_i
      integer(int64) :: c0, c1, rate
      real(dp) :: t_hessian, t_rhs, t_zvector, t_base, t_orb, t_ci
      integer :: n_states, n_mo, na, nb, n_param, l, j, seg0, seg1, iterations
      real(dp) :: residual

      n_states = size(result%energies)
      allocate (weights(n_states))
      weights = 1.0_dp/real(n_states, dp)
      n_mo = size(result%orbitals, 2)

      call system_clock(c0, rate)
      call build_sa_hessian(mol, result%orbitals, N_INACTIVE, N_ACTIVE, N_ALPHA, N_BETA, &
                            result%ci_vectors, result%energies, weights, state, error)
      call system_clock(c1)
      t_hessian = real(c1 - c0, dp)/real(rate, dp)
      if (error%has_error()) return

      call system_clock(c0, rate)
      call active_space_rdms(result%ci_vectors(:, :, 1), state%alpha, state%beta, dm1_i, &
                             dm2_i, error)
      if (error%has_error()) then
         call destroy_sa_hessian(state)
         return
      end if
      call generalized_fock(mol, result%orbitals, N_INACTIVE, N_ACTIVE, dm1_i, dm2_i, &
                            fock_i, error)
      call orbital_gradient(fock_i, N_INACTIVE, N_ACTIVE, grad_full)
      n_param = sa_hessian_n_param(state)
      allocate (rhs(n_param, 1))
      rhs = 0.0_dp
      do l = 1, state%n_rot
         rhs(l, 1) = -grad_full(state%rows(l), state%cols(l))
      end do
      call system_clock(c1)
      t_rhs = real(c1 - c0, dp)/real(rate, dp)

      call system_clock(c0, rate)
      call sa_zvector_solve(state, rhs, x, iterations, residual, 1.0e-10_dp, 200, error)
      call system_clock(c1)
      t_zvector = real(c1 - c0, dp)/real(rate, dp)
      if (error%has_error()) then
         call destroy_sa_hessian(state)
         return
      end if

      allocate (kappa_bar(n_mo, n_mo))
      kappa_bar = 0.0_dp
      do l = 1, state%n_rot
         kappa_bar(state%rows(l), state%cols(l)) = x(l, 1)
         kappa_bar(state%cols(l), state%rows(l)) = -x(l, 1)
      end do
      na = state%alpha%n_strings
      nb = state%beta%n_strings
      allocate (xbar(na, nb, n_states))
      do j = 1, n_states
         seg0 = state%n_rot + (j - 1)*state%n_det + 1
         seg1 = state%n_rot + j*state%n_det
         xbar(:, :, j) = reshape(x(seg0:seg1, 1), [na, nb])
      end do

      call system_clock(c0, rate)
      call czt_mcscf_gradient(mol, result%orbitals, N_INACTIVE, N_ACTIVE, dm1_i, dm2_i, &
                              base_g, error)
      call system_clock(c1)
      t_base = real(c1 - c0, dp)/real(rate, dp)
      if (error%has_error()) then
         call destroy_sa_hessian(state)
         return
      end if

      call system_clock(c0, rate)
      call orbital_response_gradient(mol, result%orbitals, N_INACTIVE, N_ACTIVE, state, &
                                     kappa_bar, orb_g, error)
      call system_clock(c1)
      t_orb = real(c1 - c0, dp)/real(rate, dp)
      if (error%has_error()) then
         call destroy_sa_hessian(state)
         return
      end if

      call system_clock(c0, rate)
      call ci_response_gradient(mol, result%orbitals, N_INACTIVE, N_ACTIVE, state, xbar, &
                                weights, ci_g, error)
      call system_clock(c1)
      t_ci = real(c1 - c0, dp)/real(rate, dp)

      print '(a,i0,a)', "# --- per-root breakdown (n_states=", n_states, ") ---"
      print '(a,f10.4,a)', "#   build_sa_hessian : ", t_hessian, " s"
      print '(a,f10.4,a)', "#   RHS build        : ", t_rhs, " s"
      print '(a,f10.4,a,i0,a,es10.2)', "#   Z-vector solve   : ", t_zvector, " s (", &
         iterations, " iters, residual ", residual
      print '(a,f10.4,a)', "#   base term        : ", t_base, " s"
      print '(a,f10.4,a)', "#   orbital response : ", t_orb, " s"
      print '(a,f10.4,a)', "#   CI response      : ", t_ci, " s"

      call destroy_sa_hessian(state)
   end subroutine time_breakdown

   subroutine test_timing(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: active_orbitals(:, :)
      type(casscf_result_t) :: sa2, sa3, sa4
      real(dp) :: w2(2), w3(3), w4(4)
      real(dp) :: t1, t2, t3, t4, f1, f2, f3, f4
      integer :: n_pi_occ, n_pi_virt

      call build_psb3(mol, err)
      call check(error,.not. err%has_error(), "psb3 molecule should build")
      if (allocated(error)) return

      call run_czt_rhf(mol, N_ELECTRONS, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "psb3 RHF should converge")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call select_active_space(mol, scf%orbitals, active_orbitals, n_pi_occ, n_pi_virt, err)
      call check(error,.not. err%has_error(), "psb3: active space selection should succeed")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      w2 = [0.5_dp, 0.5_dp]
      w3 = [1.0_dp, 1.0_dp, 1.0_dp]/3.0_dp
      w4 = [0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp]

      print '(a)', "# psb3 SA-CASSCF: converging SA-2, SA-3, SA-4 CAS(6,6)/6-31G*(cart)"
      print '(a,i0,a)', "# (", MACRO_ITERS, " macro-iterations per solve: a timing point, "// &
         "not a converged energy)"
      call converge_sa(mol, active_orbitals, 2, w2, sa2, err)
      call check(error,.not. err%has_error(), "SA-2 should not error")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      print '(a,l1,a,i0)', "# SA-2 converged = ", sa2%converged, ", iterations = ", &
         sa2%iterations
      print '(a,*(f18.10,1x))', "# SA-2 energies: ", sa2%energies

      call converge_sa(mol, sa2%orbitals, 3, w3, sa3, err)
      call check(error,.not. err%has_error(), "SA-3 should not error")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      print '(a,l1,a,i0)', "# SA-3 converged = ", sa3%converged, ", iterations = ", &
         sa3%iterations
      print '(a,*(f18.10,1x))', "# SA-3 energies: ", sa3%energies

      call converge_sa(mol, sa3%orbitals, 4, w4, sa4, err)
      call check(error,.not. err%has_error(), "SA-4 should not error")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      print '(a,l1,a,i0)', "# SA-4 converged = ", sa4%converged, ", iterations = ", &
         sa4%iterations
      print '(a,*(f18.10,1x))', "# SA-4 energies: ", sa4%energies

      print '(a)', "# --- profile: one root's pieces, SA-3 (representative) ---"
      call time_breakdown(mol, sa3, err)
      call check(error,.not. err%has_error(), "breakdown should not error")

      call time_unfused(mol, sa2, 1, t1, err)
      if (.not. err%has_error()) call time_unfused(mol, sa2, 2, t2, err)
      if (.not. err%has_error()) call time_unfused(mol, sa3, 3, t3, err)
      if (.not. err%has_error()) call time_unfused(mol, sa4, 4, t4, err)
      call check(error,.not. err%has_error(), "the unfused timings should not error")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      call time_fused(mol, sa2, 1, f1, err)
      if (.not. err%has_error()) call time_fused(mol, sa2, 2, f2, err)
      if (.not. err%has_error()) call time_fused(mol, sa3, 3, f3, err)
      if (.not. err%has_error()) call time_fused(mol, sa4, 4, f4, err)
      call check(error,.not. err%has_error(), "the fused timings should not error")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      print '(a)', "#  N  states   unfused (s)  cost(N)/cost(1)   fused (s)  cost(N)/cost(1)"
      print '(a,i2,a,i2,2(f14.3,f15.3))', "# ", 1, "  SA-", 2, t1, 1.0_dp, f1, 1.0_dp
      print '(a,i2,a,i2,2(f14.3,f15.3))', "# ", 2, "  SA-", 2, t2, t2/t1, f2, f2/f1
      print '(a,i2,a,i2,2(f14.3,f15.3))', "# ", 3, "  SA-", 3, t3, t3/t1, f3, f3/f1
      print '(a,i2,a,i2,2(f14.3,f15.3))', "# ", 4, "  SA-", 4, t4, t4/t1, f4, f4/f1

      deallocate (active_orbitals)
      call mol%destroy()
   end subroutine test_timing

end module test_mqc_sa_gradient_psb3_bench

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_sa_gradient_psb3_bench, only: collect_mqc_sa_gradient_psb3_bench_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_sa_gradient_psb3_bench", &
                               collect_mqc_sa_gradient_psb3_bench_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
