!! Transition density matrices between two CI vectors
module test_mqc_transition_rdm
   !! `transition_rdms` generalises `active_space_rdms` to a bra different from
   !! the ket. Five checks, each blind to a different way it could be wrong:
   !!
   !!   - `bra == ket` must reproduce `active_space_rdms` -- the routine this
   !!     one generalises, and the one every CASSCF/CASCI energy in this tree
   !!     already depends on
   !!   - the polarisation identity, exact operator algebra with no reference
   !!     code: it pins the overall normalisation and the cross term together
   !!   - the transpose/pair symmetries, again exact algebra, which catch an
   !!     index swapped between the bra and the ket
   !!   - traces, which count electrons and (for two orthogonal CASCI roots)
   !!     must vanish
   !!   - the energy contraction against two distinct CASCI roots, which must
   !!     reproduce the off-diagonal (zero) and diagonal Hamiltonian matrix
   !!     elements
   !!   - explicit elements against `pyscf.fci.direct_spin1.trans_rdm12`
   !!
   !! The model is `test_mqc_rdm.f90`'s: a Cholesky-style `(pq|rs)` so the
   !! two-electron tensor has the full eightfold symmetry `absorb_one_electron`
   !! relies on, and CAS(4,4) throughout so every one-particle index pattern
   !! (diagonal, off-diagonal, four distinct indices) appears in a 6x6
   !! determinant space.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use pic_lapack_interfaces, only: pic_syev
   use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   use mqc_error, only: error_t
   use mqc_determinants, only: link_table_t, build_link_table
   use mqc_ci, only: absorb_one_electron, ci_hamiltonian
   use mqc_rdm, only: active_space_rdms, transition_rdms, rdm_energy
   implicit none
   private

   public :: collect_mqc_transition_rdm_tests

   integer, parameter :: NORB = 4
   integer, parameter :: NCHOL = 3
   integer, parameter :: NALPHA = 2
   integer, parameter :: NBETA = 2

contains

   subroutine collect_mqc_transition_rdm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("bra_equals_ket_is_active_space_rdms", test_bra_equals_ket), &
                  new_unittest("polarisation_identity", test_polarisation), &
                  new_unittest("transpose_and_pair_symmetry", test_symmetries), &
                  new_unittest("traces", test_traces), &
                  new_unittest("energy_against_two_roots", test_energy), &
                  new_unittest("elements_against_pyscf", test_pyscf), &
                  new_unittest("refusals", test_refusals) &
                  ]
   end subroutine collect_mqc_transition_rdm_tests

   subroutine model_integrals(h1e, eri)
      !! The model from `test_mqc_rdm.f90`, repeated rather than shared so
      !! that a change made for one test file cannot silently move the other's
      !! references.
      real(dp), intent(out) :: h1e(NORB, NORB)
      real(dp), intent(out) :: eri(NORB, NORB, NORB, NORB)

      real(dp) :: b(NORB, NORB, NCHOL)
      integer :: p, q, r, s, l

      do q = 1, NORB
         do p = 1, NORB
            h1e(p, q) = -1.0_dp/real((p - 1) + (q - 1) + 2, dp)
         end do
      end do
      do l = 1, NCHOL
         do q = 1, NORB
            do p = 1, NORB
               b(p, q, l) = 1.0_dp/real((p - 1) + (q - 1) + (l - 1) + 3, dp)
            end do
         end do
      end do
      eri = 0.0_dp
      do s = 1, NORB
         do r = 1, NORB
            do q = 1, NORB
               do p = 1, NORB
                  do l = 1, NCHOL
                     eri(p, q, r, s) = eri(p, q, r, s) + b(p, q, l)*b(r, s, l)
                  end do
               end do
            end do
         end do
      end do
   end subroutine model_integrals

   subroutine formula_vector(na, nb, ca, cb, phase, vector)
      !! A deterministic, arbitrary CI vector: `sin(ca*ia + cb*ib + phase)` on
      !! 1-based string addresses. Not normalised or symmetric under anything
      !! -- these are meant to exercise every `(p,q,r,s)` pattern, not to be a
      !! physical state.
      integer, intent(in) :: na, nb
      real(dp), intent(in) :: ca, cb, phase
      real(dp), intent(out) :: vector(na, nb)

      integer :: ia, ib

      do ib = 1, nb
         do ia = 1, na
            vector(ia, ib) = sin(ca*real(ia, dp) + cb*real(ib, dp) + phase)
         end do
      end do
   end subroutine formula_vector

   subroutine build_tables(alpha, beta, err, ok)
      type(link_table_t), intent(out) :: alpha, beta
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      call build_link_table(NORB, NALPHA, alpha, err)
      call build_link_table(NORB, NBETA, beta, err)
      ok = .not. err%has_error()
   end subroutine build_tables

   subroutine two_casci_roots(h1e, eri, alpha, beta, root1, root2, e1, e2, err, ok)
      !! The two lowest eigenvectors of the model's dense CI Hamiltonian, in
      !! (n_alpha_strings, n_beta_strings) shape -- exact algebra, not an
      !! iterative solve, so they are orthogonal to machine precision and
      !! their energies are exact eigenvalues of `h1e`/`eri`.
      real(dp), intent(in) :: h1e(NORB, NORB), eri(NORB, NORB, NORB, NORB)
      type(link_table_t), intent(in) :: alpha, beta
      real(dp), allocatable, intent(out) :: root1(:, :), root2(:, :)
      real(dp), intent(out) :: e1, e2
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      real(dp), allocatable :: folded(:, :), dense(:, :), values(:)
      integer :: na, nb, info

      ok = .false.
      call absorb_one_electron(h1e, eri, NALPHA + NBETA, folded, err)
      if (err%has_error()) return
      call ci_hamiltonian(folded, alpha, beta, dense, err)
      if (err%has_error()) return

      na = alpha%n_strings
      nb = beta%n_strings
      allocate (values(na*nb))
      call pic_syev(dense, values, jobz="V", uplo="U", info=info)
      if (info /= 0) return
      e1 = values(1)
      e2 = values(2)

      allocate (root1(na, nb), root2(na, nb))
      root1 = reshape(dense(:, 1), [na, nb])
      root2 = reshape(dense(:, 2), [na, nb])
      ok = .true.
   end subroutine two_casci_roots

   subroutine test_bra_equals_ket(error)
      !! `transition_rdms(ci, ci, ...)` must equal `active_space_rdms(ci, ...)`
      !! exactly. The two routines schedule their OpenMP reduction
      !! independently, and an unordered `critical` merge is not guaranteed to
      !! land in the same order between two separate parallel regions even
      !! within one process (see `openmp-merge-order-is-not-a-race`), so this
      !! is pinned to one thread, the only regime bit-identity is testable in.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(link_table_t) :: alpha, beta
      real(dp) :: h1e(NORB, NORB), eri(NORB, NORB, NORB, NORB)
      real(dp), allocatable :: ci(:, :)
      real(dp), allocatable :: dm1a(:, :), dm2a(:, :, :, :)
      real(dp), allocatable :: dm1t(:, :), dm2t(:, :, :, :)
      logical :: ok
      integer :: threads

      call model_integrals(h1e, eri)
      call build_tables(alpha, beta, err, ok)
      call check(error, ok, "the excitation tables should build")
      if (allocated(error)) return

      allocate (ci(alpha%n_strings, beta%n_strings))
      call formula_vector(alpha%n_strings, beta%n_strings, 0.3_dp, 0.7_dp, 0.0_dp, ci)

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call active_space_rdms(ci, alpha, beta, dm1a, dm2a, err)
      if (.not. err%has_error()) call transition_rdms(ci, ci, alpha, beta, dm1t, dm2t, err)
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), "both routines should succeed")
      if (allocated(error)) then
         call alpha%destroy()
         call beta%destroy()
         return
      end if

      call check(error, all(dm1a == dm1t), &
                 "dm1 should be bit-identical to active_space_rdms")
      if (.not. allocated(error)) &
         call check(error, all(dm2a == dm2t), &
                    "dm2 should be bit-identical to active_space_rdms")
      call alpha%destroy()
      call beta%destroy()
   end subroutine test_bra_equals_ket

   subroutine test_polarisation(error)
      !! `T(a,b) + T(b,a) = [Gamma(a+b) - Gamma(a-b)] / 2`, exact operator
      !! algebra for any real `a`, `b` and any bilinear form: no reference
      !! numbers, no CASCI, just the definition.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(link_table_t) :: alpha, beta
      real(dp) :: h1e(NORB, NORB), eri(NORB, NORB, NORB, NORB)
      real(dp), allocatable :: a(:, :), b(:, :), apb(:, :), amb(:, :)
      real(dp), allocatable :: dm1ab(:, :), dm2ab(:, :, :, :)
      real(dp), allocatable :: dm1ba(:, :), dm2ba(:, :, :, :)
      real(dp), allocatable :: dm1p(:, :), dm2p(:, :, :, :)
      real(dp), allocatable :: dm1m(:, :), dm2m(:, :, :, :)
      real(dp), allocatable :: lhs1(:, :), rhs1(:, :)
      real(dp), allocatable :: lhs2(:, :, :, :), rhs2(:, :, :, :)
      logical :: ok

      call model_integrals(h1e, eri)
      call build_tables(alpha, beta, err, ok)
      call check(error, ok, "the excitation tables should build")
      if (allocated(error)) return

      allocate (a(alpha%n_strings, beta%n_strings), b(alpha%n_strings, beta%n_strings))
      allocate (apb(alpha%n_strings, beta%n_strings), amb(alpha%n_strings, beta%n_strings))
      call formula_vector(alpha%n_strings, beta%n_strings, 0.3_dp, 0.7_dp, 0.0_dp, a)
      call formula_vector(alpha%n_strings, beta%n_strings, 0.5_dp, -0.2_dp, 0.9_dp, b)
      apb = a + b
      amb = a - b

      call transition_rdms(a, b, alpha, beta, dm1ab, dm2ab, err)
      if (.not. err%has_error()) call transition_rdms(b, a, alpha, beta, dm1ba, dm2ba, err)
      if (.not. err%has_error()) call active_space_rdms(apb, alpha, beta, dm1p, dm2p, err)
      if (.not. err%has_error()) call active_space_rdms(amb, alpha, beta, dm1m, dm2m, err)
      call check(error,.not. err%has_error(), "every density build should succeed")
      if (allocated(error)) then
         call alpha%destroy()
         call beta%destroy()
         return
      end if

      lhs1 = dm1ab + dm1ba
      rhs1 = 0.5_dp*(dm1p - dm1m)
      call check(error, maxval(abs(lhs1 - rhs1)) < 1.0e-12_dp, &
                 "the polarisation identity should hold for dm1")
      if (.not. allocated(error)) then
         lhs2 = dm2ab + dm2ba
         rhs2 = 0.5_dp*(dm2p - dm2m)
         call check(error, maxval(abs(lhs2 - rhs2)) < 1.0e-12_dp, &
                    "the polarisation identity should hold for dm2")
      end if
      call alpha%destroy()
      call beta%destroy()
   end subroutine test_polarisation

   subroutine test_symmetries(error)
      !! Three exact identities, none of them needing bra == ket:
      !!
      !!   - dm1(bra,ket) = transpose(dm1(ket,bra))
      !!   - dm2(bra,ket)(p,q,r,s) = dm2(ket,bra)(q,p,s,r)
      !!   - dm2(bra,ket)(p,q,r,s) = dm2(bra,ket)(r,s,p,q), the ordinary
      !!     pq<->rs pair symmetry `active_space_rdms` also has
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: bra(:, :), ket(:, :)
      real(dp), allocatable :: dm1bk(:, :), dm2bk(:, :, :, :)
      real(dp), allocatable :: dm1kb(:, :), dm2kb(:, :, :, :)
      logical :: ok
      integer :: p, q, r, s

      call build_tables(alpha, beta, err, ok)
      call check(error, ok, "the excitation tables should build")
      if (allocated(error)) return

      allocate (bra(alpha%n_strings, beta%n_strings), ket(alpha%n_strings, beta%n_strings))
      call formula_vector(alpha%n_strings, beta%n_strings, 0.3_dp, 0.7_dp, 0.0_dp, bra)
      call formula_vector(alpha%n_strings, beta%n_strings, 0.5_dp, -0.2_dp, 0.9_dp, ket)

      call transition_rdms(bra, ket, alpha, beta, dm1bk, dm2bk, err)
      if (.not. err%has_error()) &
         call transition_rdms(ket, bra, alpha, beta, dm1kb, dm2kb, err)
      call check(error,.not. err%has_error(), "both orderings should succeed")
      if (allocated(error)) then
         call alpha%destroy()
         call beta%destroy()
         return
      end if

      call check(error, maxval(abs(dm1bk - transpose(dm1kb))) < 1.0e-12_dp, &
                 "dm1(bra,ket) should be the transpose of dm1(ket,bra)")
      if (allocated(error)) then
         call alpha%destroy()
         call beta%destroy()
         return
      end if

      do s = 1, NORB
         do r = 1, NORB
            do q = 1, NORB
               do p = 1, NORB
                  if (abs(dm2bk(p, q, r, s) - dm2kb(q, p, s, r)) > 1.0e-12_dp) then
                     call check(error, .false., &
                                "dm2(bra,ket)(p,q,r,s) should equal "// &
                                "dm2(ket,bra)(q,p,s,r)")
                     call alpha%destroy()
                     call beta%destroy()
                     return
                  end if
               end do
            end do
         end do
      end do

      call check(error, maxval(abs(dm2bk - reshape(dm2bk, [NORB, NORB, NORB, NORB], &
                                                   order=[3, 4, 1, 2]))) < 1.0e-12_dp, &
                 "dm2(bra,ket) should be symmetric under exchanging the pq and rs pairs")
      call alpha%destroy()
      call beta%destroy()
   end subroutine test_symmetries

   subroutine test_traces(error)
      !! `sum_p dm1(p,p) = N <bra|ket>` and `sum_pr dm2(p,p,r,r) = N(N-1)
      !! <bra|ket>`, both because `sum_p E_pp` is the (constant, N) number
      !! operator on this fixed-particle-count space. Checked twice: once with
      !! an arbitrary, non-orthogonal `(bra,ket)` pair, where the identity is
      !! genuinely exercised (a sign error in `N` would still read as zero
      !! against an orthogonal pair), and once with two orthogonal CASCI roots,
      !! where both sides collapse to zero.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(link_table_t) :: alpha, beta
      real(dp) :: h1e(NORB, NORB), eri(NORB, NORB, NORB, NORB)
      real(dp), allocatable :: bra(:, :), ket(:, :)
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      real(dp), allocatable :: root1(:, :), root2(:, :)
      real(dp) :: overlap, trace1, trace2, e1, e2
      logical :: ok
      integer :: p, r
      integer, parameter :: NTOT = NALPHA + NBETA

      call model_integrals(h1e, eri)
      call build_tables(alpha, beta, err, ok)
      call check(error, ok, "the excitation tables should build")
      if (allocated(error)) return

      allocate (bra(alpha%n_strings, beta%n_strings), ket(alpha%n_strings, beta%n_strings))
      call formula_vector(alpha%n_strings, beta%n_strings, 0.3_dp, 0.7_dp, 0.0_dp, bra)
      call formula_vector(alpha%n_strings, beta%n_strings, 0.5_dp, -0.2_dp, 0.9_dp, ket)
      overlap = sum(bra*ket)

      call transition_rdms(bra, ket, alpha, beta, dm1, dm2, err)
      call check(error,.not. err%has_error(), "the transition build should succeed")
      if (allocated(error)) then
         call alpha%destroy()
         call beta%destroy()
         return
      end if

      trace1 = 0.0_dp
      do p = 1, NORB
         trace1 = trace1 + dm1(p, p)
      end do
      call check(error, trace1, real(NTOT, dp)*overlap, &
                 "sum_p dm1(p,p) should be N times the overlap", thr=1.0e-11_dp)
      if (allocated(error)) then
         call alpha%destroy()
         call beta%destroy()
         return
      end if

      trace2 = 0.0_dp
      do r = 1, NORB
         do p = 1, NORB
            trace2 = trace2 + dm2(p, p, r, r)
         end do
      end do
      call check(error, trace2, real(NTOT*(NTOT - 1), dp)*overlap, &
                 "sum_pr dm2(p,p,r,r) should be N(N-1) times the overlap", thr=1.0e-10_dp)
      if (allocated(error)) then
         call alpha%destroy()
         call beta%destroy()
         return
      end if

      call two_casci_roots(h1e, eri, alpha, beta, root1, root2, e1, e2, err, ok)
      call check(error, ok, "the two CASCI roots should be found")
      if (.not. allocated(error)) then
         call transition_rdms(root1, root2, alpha, beta, dm1, dm2, err)
         call check(error,.not. err%has_error(), "the orthogonal-root build should succeed")
      end if
      if (allocated(error)) then
         call alpha%destroy()
         call beta%destroy()
         return
      end if

      trace1 = 0.0_dp
      do p = 1, NORB
         trace1 = trace1 + dm1(p, p)
      end do
      call check(error, abs(trace1) < 1.0e-10_dp, &
                 "orthogonal roots should give a zero one-particle trace")
      if (.not. allocated(error)) then
         trace2 = 0.0_dp
         do r = 1, NORB
            do p = 1, NORB
               trace2 = trace2 + dm2(p, p, r, r)
            end do
         end do
         call check(error, abs(trace2) < 1.0e-9_dp, &
                    "orthogonal roots should give a zero two-particle trace")
      end if
      call alpha%destroy()
      call beta%destroy()
   end subroutine test_traces

   subroutine test_energy(error)
      !! `sum_pq h_pq dm1(p,q) + (1/2) sum_pqrs (pq|rs) dm2(p,q,r,s)` is the
      !! matrix element `<bra|H_active|ket>` of the active-space Hamiltonian
      !! -- `rdm_energy` computes exactly this contraction, unmodified, fed a
      !! transition rather than a same-state pair. For two distinct
      !! eigenvectors of the model's dense CI Hamiltonian this is an
      !! off-diagonal element and must be zero; fed one root against itself it
      !! must reproduce that root's eigenvalue.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(link_table_t) :: alpha, beta
      real(dp) :: h1e(NORB, NORB), eri(NORB, NORB, NORB, NORB)
      real(dp), allocatable :: root1(:, :), root2(:, :)
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      real(dp) :: e1, e2
      logical :: ok

      call model_integrals(h1e, eri)
      call build_tables(alpha, beta, err, ok)
      call check(error, ok, "the excitation tables should build")
      if (allocated(error)) return

      call two_casci_roots(h1e, eri, alpha, beta, root1, root2, e1, e2, err, ok)
      call check(error, ok, "the two CASCI roots should be found")
      if (allocated(error)) then
         call alpha%destroy()
         call beta%destroy()
         return
      end if

      call transition_rdms(root1, root2, alpha, beta, dm1, dm2, err)
      call check(error,.not. err%has_error(), "the off-diagonal transition build should succeed")
      if (.not. allocated(error)) &
         call check(error, abs(rdm_energy(h1e, eri, dm1, dm2)) < 1.0e-9_dp, &
                    "the off-diagonal energy contraction should vanish")
      if (allocated(error)) then
         call alpha%destroy()
         call beta%destroy()
         return
      end if

      call transition_rdms(root1, root1, alpha, beta, dm1, dm2, err)
      call check(error,.not. err%has_error(), "the diagonal transition build should succeed")
      if (.not. allocated(error)) &
         call check(error, rdm_energy(h1e, eri, dm1, dm2), e1, &
                    "the diagonal energy contraction should reproduce the root's eigenvalue", &
                    thr=1.0e-10_dp)
      call alpha%destroy()
      call beta%destroy()
   end subroutine test_energy

   subroutine test_pyscf(error)
      !! Explicit elements against `pyscf.fci.direct_spin1.trans_rdm12` on the
      !! same model and the same bra/ket formulas, generated by
      !! `tools/sa_casscf/trans_rdm_ref.py`.
      !!
      !! PySCF's own convention, read from `trans_rdm12`'s docstring, is
      !! `1pdm[p,q] = <q^dagger p>` and `2pdm[p,q,r,s] = <p^dagger r^dagger s
      !! q>` (0-based). The two-particle formula is the standard
      !! `a_p^dagger a_r^dagger a_s a_q = E_pq E_rs - delta_qr E_ps` identity,
      !! so `2pdm` lines up with this module's `dm2` index for index, no
      !! permutation. The one-particle formula is `<q^dagger p> = <E_qp>`,
      !! so `dm1(p,q)` here is PySCF's `1pdm[q-1,p-1]` (0-based, transposed) --
      !! confirmed against PySCF numerically on a toy one-electron case before
      !! trusting it here, since a transition dm1 need not be symmetric and
      !! the transpose would otherwise pass unnoticed.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: bra(:, :), ket(:, :)
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      logical :: ok

      call build_tables(alpha, beta, err, ok)
      call check(error, ok, "the excitation tables should build")
      if (allocated(error)) return

      allocate (bra(alpha%n_strings, beta%n_strings), ket(alpha%n_strings, beta%n_strings))
      call formula_vector(alpha%n_strings, beta%n_strings, 0.3_dp, 0.7_dp, 0.0_dp, bra)
      call formula_vector(alpha%n_strings, beta%n_strings, 0.5_dp, -0.2_dp, 0.9_dp, ket)

      call transition_rdms(bra, ket, alpha, beta, dm1, dm2, err)
      call check(error,.not. err%has_error(), "the transition build should succeed")
      if (allocated(error)) then
         call alpha%destroy()
         call beta%destroy()
         return
      end if

      call check(error, dm1(1, 1), 4.157307965617898_dp, "dm1(1,1)", thr=1.0e-12_dp)
      if (.not. allocated(error)) &
         call check(error, dm1(2, 2), -0.884781587483814_dp, "dm1(2,2)", thr=1.0e-12_dp)
      if (.not. allocated(error)) &
         call check(error, dm1(1, 2), -1.1852039606306641_dp, "dm1(1,2)", thr=1.0e-12_dp)
      if (.not. allocated(error)) &
         call check(error, dm1(2, 1), -4.2392123459681885_dp, "dm1(2,1)", thr=1.0e-12_dp)
      if (.not. allocated(error)) &
         call check(error, dm1(3, 4), 0.52738992006519_dp, "dm1(3,4)", thr=1.0e-12_dp)
      if (.not. allocated(error)) &
         call check(error, dm1(4, 1), 7.3894339450089985_dp, "dm1(4,1)", thr=1.0e-12_dp)
      if (allocated(error)) then
         call alpha%destroy()
         call beta%destroy()
         return
      end if

      call check(error, dm2(1, 1, 1, 1), 6.673378612038242_dp, "dm2(1,1,1,1)", thr=1.0e-11_dp)
      if (.not. allocated(error)) &
         call check(error, dm2(2, 3, 4, 1), 4.2765409538507395_dp, "dm2(2,3,4,1)", &
                    thr=1.0e-11_dp)
      if (.not. allocated(error)) &
         call check(error, dm2(1, 2, 2, 1), -5.634495964874194_dp, "dm2(1,2,2,1)", &
                    thr=1.0e-11_dp)
      if (.not. allocated(error)) &
         call check(error, dm2(3, 1, 2, 4), 0.2802969603785896_dp, "dm2(3,1,2,4)", &
                    thr=1.0e-11_dp)
      if (.not. allocated(error)) &
         call check(error, dm2(4, 4, 1, 1), -8.720551249300046_dp, "dm2(4,4,1,1)", &
                    thr=1.0e-11_dp)
      if (.not. allocated(error)) &
         call check(error, dm2(2, 1, 3, 4), -1.1451098460568194_dp, "dm2(2,1,3,4)", &
                    thr=1.0e-11_dp)
      call alpha%destroy()
      call beta%destroy()
   end subroutine test_pyscf

   subroutine test_refusals(error)
      !! Vectors that do not belong to the tables, in either slot
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :), bra(:, :), ket(:, :)
      logical :: ok

      call build_tables(alpha, beta, err, ok)
      call check(error, ok, "the excitation tables should build")
      if (allocated(error)) return

      allocate (bra(alpha%n_strings + 1, beta%n_strings), ket(alpha%n_strings, beta%n_strings))
      bra = 0.0_dp
      ket = 0.0_dp
      call transition_rdms(bra, ket, alpha, beta, dm1, dm2, err)
      call check(error, err%has_error(), "a mis-shaped bra should be refused")
      call err%clear()
      deallocate (bra, ket)

      allocate (bra(alpha%n_strings, beta%n_strings), ket(alpha%n_strings, beta%n_strings + 1))
      bra = 0.0_dp
      ket = 0.0_dp
      call transition_rdms(bra, ket, alpha, beta, dm1, dm2, err)
      call check(error, err%has_error(), "a mis-shaped ket should be refused")
      call err%clear()

      call alpha%destroy()
      call beta%destroy()
   end subroutine test_refusals

end module test_mqc_transition_rdm

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_transition_rdm, only: collect_mqc_transition_rdm_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_transition_rdm", collect_mqc_transition_rdm_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
