!! SAPT0's molecules: the dimer basis, and ghosts in it
module test_mqc_sapt
   !! M1 of the SAPT ladder. What is asserted here is everything that can be
   !! checked before a single term is computed, and each assertion corresponds to
   !! a mistake that is otherwise silent:
   !!
   !!   * a ghost that kept its nuclear charge -- caught by the charge sum
   !!   * a ghost that lost its basis functions -- caught by the AO count
   !!   * a monomer built in its own atom order rather than the dimer's -- caught
   !!     by `Tr(D S)`, which is the occupied count only when the AO orderings
   !!     agree. This one is not hypothetical; the PySCF prototype was written
   !!     that way first and every AO quantity was quietly permuted.
   !!   * ghosting that changed the nuclear repulsion of the real atoms -- caught
   !!     by comparing against the monomer alone
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_sapt, only: sapt_molecules_t, build_sapt_molecules, &
                       sapt_cache_t, build_sapt_cache, sapt_elst10, sapt_exch10_s2, sapt_exch10, &
                       sapt_induction, sapt_terms_t, sapt_disp20, sapt_exch_disp20
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_sapt_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_sapt_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("sapt_dimer_basis_and_ghosts", test_ghosts), &
                  new_unittest("sapt_monomer_ao_ordering", test_ordering), &
                  new_unittest("sapt_bsse_is_negative", test_bsse), &
                  new_unittest("sapt_elst10", test_elst10), &
                  new_unittest("sapt_exch10_s2", test_exch10_s2), &
                  new_unittest("sapt_exch10_sinf", test_exch10), &
                  new_unittest("sapt_induction", test_induction), &
                  new_unittest("sapt_disp20", test_disp20), &
                  new_unittest("sapt_exch_disp20", test_exch_disp20) &
                  ]
   end subroutine collect_mqc_sapt_tests

   subroutine water_dimer(mols, err)
      !! Two waters 3.0 Angstrom apart along x -- the geometry every EFP number
      !! in this repository is pinned to, so the two decompositions can one day
      !! be read against each other on a case where the EFP side is exact.
      type(sapt_molecules_t), intent(out) :: mols
      type(error_t), intent(inout) :: err

      integer :: z(3)
      character(len=2) :: sym(3)
      real(dp) :: a(3, 3), b(3, 3)

      z = [8, 1, 1]
      sym = ["O ", "H ", "H "]
      a = reshape([0.0_dp, 0.0_dp, 0.10077199_dp, &
                   0.0_dp, 0.77250895_dp, -0.46780200_dp, &
                   0.0_dp, -0.77250895_dp, -0.46780200_dp], [3, 3])*ANG
      b = a
      b(1, :) = b(1, :) + 3.0_dp*ANG
      call build_sapt_molecules(z, sym, a, z, sym, b, "6-31g", mols, err)
   end subroutine water_dimer

   subroutine test_ghosts(error)
      !! Both monomers span the dimer basis, and only their own atoms are charged
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(error_t) :: err

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return

      call check(error, mols%mono_a%nao == mols%dimer%nao, &
                 "monomer A must span the whole dimer basis")
      if (allocated(error)) return
      call check(error, mols%mono_b%nao == mols%dimer%nao, &
                 "monomer B must span the whole dimer basis")
      if (allocated(error)) return

      ! Charges: A's three atoms carry 10, B's ghosts carry nothing.
      call check(error, abs(sum(mols%mono_a%charges(1:3)) - 10.0_dp) < 1.0e-12_dp, &
                 "monomer A's own atoms lost their nuclear charge")
      if (allocated(error)) return
      call check(error, abs(sum(mols%mono_a%charges(4:6))) < 1.0e-12_dp, &
                 "monomer A's ghosts kept a nuclear charge")
      if (allocated(error)) return
      call check(error, abs(sum(mols%mono_b%charges(1:3))) < 1.0e-12_dp, &
                 "monomer B's ghosts kept a nuclear charge")
      if (allocated(error)) return

      ! And ghosting must not disturb the real atoms' mutual repulsion.
      call check(error, mols%mono_a%nuclear_repulsion() > 0.0_dp, &
                 "monomer A should still have a nuclear repulsion")
      if (allocated(error)) return
      call check(error, mols%dimer%nuclear_repulsion() > &
                 mols%mono_a%nuclear_repulsion() + mols%mono_b%nuclear_repulsion(), &
                 "the dimer's nuclear repulsion must exceed the monomers' sum")

      call mols%destroy()
   end subroutine test_ghosts

   subroutine test_ordering(error)
      !! `Tr(D S) = n_occ`, which holds only if the AO orderings agree
      !!
      !! `D` comes from a monomer's SCF and `S` from the dimer. If the monomer
      !! were built in its own atom order the two would be permuted relative to
      !! each other, the trace would come out wrong, and every SAPT term after it
      !! would be quietly meaningless.
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(error_t) :: err
      type(rhf_result_t) :: scf
      real(dp), allocatable :: s(:, :), d(:, :), ds(:, :)
      real(dp) :: trace
      integer :: i, nocc

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return

      call run_libcint_rhf(mols%mono_a, mols%n_elec_a, 200, 1.0e-11_dp, 1.0e-9_dp, &
                           .false., scf, err, in_core=.true.)
      call check(error,.not. err%has_error() .and. scf%converged, &
                 "monomer A's SCF failed in the dimer basis")
      if (allocated(error)) return

      nocc = mols%n_elec_a/2
      call mols%dimer%overlap(s)
      allocate (d(mols%dimer%nao, mols%dimer%nao))
      allocate (ds(mols%dimer%nao, mols%dimer%nao))
      call pic_gemm(scf%orbitals(:, 1:nocc), scf%orbitals(:, 1:nocc), d, transb="T")
      call pic_gemm(d, s, ds)
      trace = 0.0_dp
      do i = 1, mols%dimer%nao
         trace = trace + ds(i, i)
      end do

      call check(error, abs(trace - real(nocc, dp)) < 1.0e-8_dp, &
                 "Tr(D S) is not the occupied count -- the monomer's AO ordering "// &
                 "does not match the dimer's")

      deallocate (s, d, ds)
      call mols%destroy()
   end subroutine test_ordering

   subroutine test_bsse(error)
      !! The dimer basis must lower a monomer's energy, and by a sane amount
      !!
      !! No wrong ghost has this property: one that kept its nuclear charge gives
      !! a wildly lower energy, and one that lost its basis functions gives
      !! exactly the monomer's own energy.
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(libcint_molecule_t) :: alone
      type(rhf_result_t) :: dcbs, own
      type(error_t) :: err
      integer :: z(3)
      character(len=2) :: sym(3)
      real(dp) :: a(3, 3), bsse

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return

      z = [8, 1, 1]
      sym = ["O ", "H ", "H "]
      a = reshape([0.0_dp, 0.0_dp, 0.10077199_dp, &
                   0.0_dp, 0.77250895_dp, -0.46780200_dp, &
                   0.0_dp, -0.77250895_dp, -0.46780200_dp], [3, 3])*ANG
      call build_libcint_molecule(z, sym, a, "6-31g", alone, err)
      call check(error,.not. err%has_error(), "building the lone monomer failed")
      if (allocated(error)) return

      call run_libcint_rhf(mols%mono_a, mols%n_elec_a, 200, 1.0e-11_dp, 1.0e-9_dp, &
                           .false., dcbs, err, in_core=.true.)
      if (.not. err%has_error()) then
         call run_libcint_rhf(alone, mols%n_elec_a, 200, 1.0e-11_dp, 1.0e-9_dp, &
                              .false., own, err, in_core=.true.)
      end if
      call check(error,.not. err%has_error(), "an SCF failed")
      if (allocated(error)) return

      bsse = dcbs%energy - own%energy
      call check(error, bsse < 0.0_dp, &
                 "the dimer basis must lower the monomer's energy")
      if (allocated(error)) return
      ! Loose upper bound: a ghost that kept its charge is worth whole Hartrees,
      ! not milli-Hartrees, so this separates the two cases by orders of magnitude.
      call check(error, bsse > -0.1_dp, &
                 "the basis-set superposition error is far too large, which is "// &
                 "what a ghost that kept its nuclear charge looks like")

      call alone%destroy()
      call mols%destroy()
   end subroutine test_bsse

   subroutine test_elst10(error)
      !! `Elst10,r` against the conventional four-index reference
      !!
      !! The reference is `validation/check_sapt0.py`, which computes the same
      !! quantity in PySCF with exact integrals. psi4 cannot serve here: its
      !! closed-shell SAPT is density-fitted by construction, so its number is a
      !! different quantity in the same basis and agrees only to ~1e-6.
      !!
      !! Regenerate with:  python validation/check_sapt0.py --json <path>
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(error_t) :: err
      real(dp) :: e

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return

      e = sapt_elst10(cache)
      ! Conventional four-index, water dimer 3.0 Angstrom along x, 6-31G.
      call check(error, e, 0.006384228637_dp, thr=1.0e-9_dp, &
                 message="Elst10,r disagrees with the PySCF reference")

      call cache%destroy()
      call mols%destroy()
   end subroutine test_elst10

   subroutine test_exch10_s2(error)
      !! `Exch10(S^2)` against the conventional four-index reference
      !!
      !! The term that exposed psi4's `vector_dot` being elementwise rather than
      !! a trace: two of its operands are non-symmetric, and taking `Tr(X Y)`
      !! there puts the answer 55% high while leaving every shape valid.
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(error_t) :: err
      real(dp) :: e

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return

      e = sapt_exch10_s2(cache)
      call check(error, e, 0.003235393978_dp, thr=1.0e-9_dp, &
                 message="Exch10(S^2) disagrees with the PySCF reference")
      if (allocated(error)) return
      ! First-order exchange is repulsive. A sign slip here still produces a
      ! plausible-looking interaction energy once it is summed with the rest.
      call check(error, e > 0.0_dp, "first-order exchange must be repulsive")

      call cache%destroy()
      call mols%destroy()
   end subroutine test_exch10_s2

   subroutine test_exch10(error)
      !! `Exch10` at S^inf -- the form the SAPT0 total actually uses
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(error_t) :: err
      real(dp) :: e, e_s2

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return

      e = sapt_exch10(cache, err)
      call check(error,.not. err%has_error(), "Exch10 failed")
      if (allocated(error)) return
      call check(error, e, 0.003241500408_dp, thr=1.0e-9_dp, &
                 message="Exch10 disagrees with the PySCF reference")
      if (allocated(error)) return

      ! The S^2 form is the leading term of this one, so at a separation where
      ! the overlap is small the two must be close and S^inf the larger. They
      ! part company as the monomers approach -- 19% at 2 Angstrom on LiH.
      e_s2 = sapt_exch10_s2(cache)
      call check(error, e > e_s2, "S^inf exchange should exceed the S^2 form here")
      if (allocated(error)) return
      call check(error, abs(e - e_s2) < 0.05_dp*abs(e), &
                 "S^2 and S^inf differ by far more than the overlap warrants")

      call cache%destroy()
      call mols%destroy()
   end subroutine test_exch10

   subroutine test_induction(error)
      !! `Ind20` and `Exch-Ind20`, uncoupled and coupled
      !!
      !! The uncoupled pair is asserted as well as the coupled one, because it
      !! needs no response solver: if it passes and the coupled one does not, the
      !! fault is in `cphf_solve`'s conventions rather than in the contraction,
      !! and the two are otherwise indistinguishable from one failing number.
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(sapt_terms_t) :: t
      type(error_t) :: err

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return

      call sapt_induction(mols, cache, t, err)
      call check(error,.not. err%has_error(), "induction failed")
      if (allocated(error)) return

      call check(error, t%ind20_u, -0.000954737579_dp, thr=1.0e-9_dp, &
                 message="Ind20,u disagrees with the PySCF reference")
      if (allocated(error)) return
      call check(error, t%exch_ind20_u, 0.000805018637_dp, thr=1.0e-9_dp, &
                 message="Exch-Ind20,u disagrees with the PySCF reference")
      if (allocated(error)) return
      call check(error, t%ind20_r, -0.001134897260_dp, thr=1.0e-8_dp, &
                 message="Ind20,r disagrees with the PySCF reference")
      if (allocated(error)) return
      call check(error, t%exch_ind20_r, 0.000949574150_dp, thr=1.0e-8_dp, &
                 message="Exch-Ind20,r disagrees with the PySCF reference")
      if (allocated(error)) return

      ! Induction is stabilising and its exchange counterpart repulsive, and
      ! letting the monomer relax can only lower the energy further.
      call check(error, t%ind20_r < t%ind20_u, &
                 "the response must be more stabilising than the uncoupled form")
      if (allocated(error)) return
      call check(error, t%exch_ind20_r > 0.0_dp, &
                 "exchange-induction must be repulsive")

      call cache%destroy()
      call mols%destroy()
   end subroutine test_induction

   subroutine test_disp20(error)
      !! `Disp20` against the conventional four-index reference
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(error_t) :: err
      real(dp) :: e

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return

      e = sapt_disp20(cache)
      call check(error, e, -0.000317465444_dp, thr=1.0e-9_dp, &
                 message="Disp20 disagrees with the PySCF reference")
      if (allocated(error)) return
      ! Dispersion is attractive. The denominator is negative throughout, so a
      ! positive answer means the occupied and virtual energies were swapped.
      call check(error, e < 0.0_dp, "dispersion must be attractive")

      call cache%destroy()
      call mols%destroy()
   end subroutine test_disp20

   subroutine test_exch_disp20(error)
      !! `Exch-Disp20` against the conventional four-index reference
      type(error_type), allocatable, intent(out) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_cache_t) :: cache
      type(error_t) :: err
      real(dp) :: e, d

      call water_dimer(mols, err)
      call check(error,.not. err%has_error(), "building the molecules failed")
      if (allocated(error)) return
      call build_sapt_cache(mols, cache, err)
      call check(error,.not. err%has_error(), "the SAPT cache failed")
      if (allocated(error)) return

      e = sapt_exch_disp20(cache)
      call check(error, e, 0.000117840542_dp, thr=1.0e-9_dp, &
                 message="Exch-Disp20 disagrees with the PySCF reference")
      if (allocated(error)) return

      ! Repulsive, and a sizeable fraction of dispersion -- 10-25% across psi4's
      ! own tests, always opposite in sign. A term that quietly evaluated to zero
      ! would leave a total that is wrong in a systematic, plausible direction.
      d = sapt_disp20(cache)
      call check(error, e > 0.0_dp, "exchange-dispersion must be repulsive")
      if (allocated(error)) return
      call check(error, abs(e) > 0.05_dp*abs(d) .and. abs(e) < 0.6_dp*abs(d), &
                 "exchange-dispersion is not a plausible fraction of dispersion")

      call cache%destroy()
      call mols%destroy()
   end subroutine test_exch_disp20

end module test_mqc_sapt

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_sapt, only: collect_mqc_sapt_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_sapt", collect_mqc_sapt_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
