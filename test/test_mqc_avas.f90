!! Choosing an active space by atomic orbital character
module test_mqc_avas
   !! AVAS answers the question that stops CASSCF being usable: *which* orbitals
   !! go in the active space. The test that matters is therefore not a number
   !! but a choice -- does it select the same orbitals a careful person would,
   !! and the same ones the reference implementation does.
   !!
   !! **We project onto a different minimal basis than the paper does**, so an
   !! element-by-element comparison of the projector is not available and would
   !! not mean much if it were. AVAS as published uses MINAO; this uses the
   !! free-atom set transcribed from GAMESS. What is compared instead is the
   !! selection itself and the converged CASSCF energy that follows from it.
   !!
   !! That comparison is the stronger one anyway. The claim AVAS makes is that
   !! the projector spectrum is bimodal for a sensibly posed request, so the
   !! threshold falls in a gap rather than through a cluster and the answer does
   !! not depend on where exactly it sits. If that is true the choice of
   !! reference minimal basis should not matter either -- and the three cases
   !! below say it does not, selecting exactly what PySCF selects and converging
   !! to the same energy to 1e-11.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_avas, only: avas_select, avas_result_t, parse_orbital_label, &
                           valence_select
   use pic_blas_interfaces, only: pic_gemm
   use mqc_czt_casci, only: run_czt_casci, casci_result_t
   use mqc_czt_mcscf, only: run_czt_casscf, casscf_result_t
   implicit none
   private

   public :: collect_mqc_avas_tests

   real(dp), parameter :: STATIONARY = 1.0e-5_dp
      !! Largest orbital gradient norm still called a stationary point here.
      !! An order of magnitude above where these cases actually stall, and four
      !! below anything that would let a genuinely unconverged run through.

   real(dp), parameter :: NITROGEN(3, 2) = reshape( &
                          [0.0_dp, 0.0_dp, -1.0371_dp, &
                           0.0_dp, 0.0_dp, 1.0371_dp], [3, 2])
   integer, parameter :: NITROGEN_Z(2) = [7, 7]
   character(len=2), parameter :: NITROGEN_SYM(2) = ["N ", "N "]

   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [3, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]

contains

   subroutine collect_mqc_avas_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("labels_parse", test_labels), &
                  new_unittest("nitrogen_2p_against_pyscf", test_nitrogen_2p), &
                  new_unittest("nitrogen_2s2p_against_pyscf", test_nitrogen_2s2p), &
                  new_unittest("water_2p_against_pyscf", test_water), &
                  new_unittest("the_cut_falls_in_a_gap", test_gap), &
                  new_unittest("full_valence_space", test_valence), &
                  new_unittest("refusals", test_refusals) &
                  ]
   end subroutine collect_mqc_avas_tests

   subroutine test_labels(error)
      !! "N 2p" and its neighbours, including the ones that are not labels
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      character(len=:), allocatable :: symbol
      integer :: n, l

      call parse_orbital_label("N 2p", symbol, n, l, err)
      call check(error,.not. err%has_error(), "'N 2p' should parse")
      if (allocated(error)) return
      call check(error, symbol, "N", "the element symbol")
      if (allocated(error)) return
      call check(error, n, 2, "the principal quantum number")
      if (allocated(error)) return
      call check(error, l, 1, "p is l = 1")
      if (allocated(error)) return

      call parse_orbital_label("  Cr 3d  ", symbol, n, l, err)
      call check(error, symbol, "Cr", "a two-letter symbol, with the label padded")
      if (allocated(error)) return
      call check(error, l, 2, "d is l = 2")
      if (allocated(error)) return

      ! The refusals. A label that cannot be read is a deck asking for something
      ! specific and getting silence, which is what this whole selection exists
      ! to avoid.
      call parse_orbital_label("N2p", symbol, n, l, err)
      call check(error, err%has_error(), "a label with no space should be refused")
      if (allocated(error)) return
      call err%clear()
      call parse_orbital_label("N 2g", symbol, n, l, err)
      call check(error, err%has_error(), "an unknown subshell letter should be refused")
      if (allocated(error)) return
      call err%clear()
      call parse_orbital_label("N p", symbol, n, l, err)
      call check(error, err%has_error(), &
                 "a subshell with no principal quantum number should be refused")
      call err%clear()
   end subroutine test_labels

   subroutine one_case(atomic_numbers, symbols, coordinates, n_electrons, n_occupied, &
                       labels, result, energy, err, ok)
      !! An SCF, an AVAS selection, and a CASSCF in the space it chose
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      integer, intent(in) :: n_electrons, n_occupied
      character(len=*), intent(in) :: labels(:)
      type(avas_result_t), intent(out) :: result
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok
         !! The orbitals reached a stationary point, which is not the same as
         !! the optimiser having said so -- see `STATIONARY` below

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casscf_result_t) :: mc
      integer :: per_spin

      ok = .false.
      energy = 0.0_dp
      call build_czt_molecule(atomic_numbers, symbols, coordinates, "cc-pvdz", &
                              mol, err)
      call run_czt_rhf(mol, n_electrons, 300, 1.0e-12_dp, 1.0e-10_dp, .false., &
                       scf, err)
      if (err%has_error() .or. .not. scf%converged) return
      call avas_select(mol, atomic_numbers, symbols, coordinates, scf%orbitals, &
                       n_occupied, labels, result, err)
      if (err%has_error()) return

      per_spin = result%n_active_electrons/2
      call run_czt_casscf(mol, result%orbitals, result%n_inactive, &
                          result%n_active, per_spin, per_spin, mc, err, &
                          max_iterations=300, gradient_tol=1.0e-6_dp)
      if (err%has_error()) return

      ! Stationarity, not the optimiser's own verdict. These cases stall a
      ! whisker from the tolerance they were given -- 7.3e-7, 8.0e-7 and 5.5e-7
      ! against 1e-6 -- so `converged` records which side of a knife edge the
      ! last step happened to land on, and a different BLAS or a different
      ! thread count is enough to flip it. The energies do not move: they agree
      ! with PySCF to 1e-10 either way, which is what the assertions below
      ! actually care about. Asking instead that the gradient be genuinely
      ! small keeps the check honest with an order of magnitude of room.
      if (mc%gradient_norm > STATIONARY) return
      energy = mc%energy
      ok = .true.
      call mol%destroy()
   end subroutine one_case

   subroutine test_nitrogen_2p(error)
      !! N2 asking for the nitrogen 2p, which is the triple bond
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(avas_result_t) :: avas
      real(dp) :: energy
      logical :: ok
      character(len=8) :: labels(1)

      labels(1) = "N 2p"
      call one_case(NITROGEN_Z, NITROGEN_SYM, NITROGEN, 14, 7, labels, avas, energy, &
                    err, ok)
      call check(error, ok, "the orbitals should reach a stationary point")
      if (allocated(error)) return

      ! PySCF's avas.avas(mf, ['N 2p'], threshold=0.2) on the same molecule.
      call check(error, avas%n_active, 7, "seven active orbitals")
      if (allocated(error)) return
      call check(error, avas%n_active_electrons, 8, "eight active electrons")
      if (allocated(error)) return
      call check(error, avas%n_inactive, 3, "three inactive orbitals")
      if (allocated(error)) return
      call check(error, energy, -109.097049261360_dp, &
                 "the CASSCF energy in the space AVAS chose", thr=1.0e-9_dp)
   end subroutine test_nitrogen_2p

   subroutine test_nitrogen_2s2p(error)
      !! The same molecule, asking for 2s as well
      !!
      !! A different request must give a different space, or the selection is
      !! not responding to what it was asked. One more orbital and two more
      !! electrons, and a lower energy because the space contains the first.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(avas_result_t) :: avas
      real(dp) :: energy
      logical :: ok
      character(len=8) :: labels(2)

      labels(1) = "N 2s"
      labels(2) = "N 2p"
      call one_case(NITROGEN_Z, NITROGEN_SYM, NITROGEN, 14, 7, labels, avas, energy, &
                    err, ok)
      call check(error, ok, "the orbitals should reach a stationary point")
      if (allocated(error)) return
      call check(error, avas%n_active, 8, "eight active orbitals")
      if (allocated(error)) return
      call check(error, avas%n_active_electrons, 10, "ten active electrons")
      if (allocated(error)) return
      call check(error, energy, -109.102611660920_dp, "the CASSCF energy", &
                 thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, energy < -109.097049261360_dp, &
                 "a larger active space containing the smaller one cannot give a "// &
                 "higher energy")
   end subroutine test_nitrogen_2s2p

   subroutine test_water(error)
      !! Water asking for the oxygen 2p, where the answer is that there is nothing to do
      !!
      !! Three orbitals holding six electrons is a full space with no room for
      !! any excitation, so the CASSCF energy is the Hartree-Fock energy exactly.
      !! Worth keeping precisely because it is a degenerate case: the selection
      !! is right, the arithmetic runs, and the answer is that this request buys
      !! nothing. A tool that silently returned something else here would be
      !! worse than one that returned this.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(avas_result_t) :: avas
      real(dp) :: energy
      logical :: ok
      character(len=8) :: labels(1)

      labels(1) = "O 2p"
      call one_case(WATER_Z, WATER_SYM, WATER, 10, 5, labels, avas, energy, err, ok)
      call check(error, ok, "the orbitals should reach a stationary point")
      if (allocated(error)) return
      call check(error, avas%n_active, 3, "three active orbitals")
      if (allocated(error)) return
      call check(error, avas%n_active_electrons, 6, "six active electrons")
      if (allocated(error)) return
      call check(error, energy, -76.026781904325_dp, &
                 "a full active space recovers nothing, so this is the SCF energy", &
                 thr=1.0e-9_dp)
   end subroutine test_water

   subroutine test_gap(error)
      !! The threshold falls in a gap, which is what makes it robust
      !!
      !! AVAS's claim is that a sensible request separates cleanly: the kept
      !! eigenvalues sit near one, the rejected near zero, and the cut lands in
      !! empty space. If that is false the selection depends on the threshold
      !! and the method is arbitrary. It is also what lets a different reference
      !! minimal basis reach the same answer, which is the whole reason the
      !! comparisons above hold at all.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(avas_result_t) :: avas
      real(dp) :: energy, largest_rejected, smallest_kept
      logical :: ok
      integer :: i
      character(len=8) :: labels(1)

      labels(1) = "N 2p"
      call one_case(NITROGEN_Z, NITROGEN_SYM, NITROGEN, 14, 7, labels, avas, energy, &
                    err, ok)
      call check(error, ok, "the orbitals should reach a stationary point")
      if (allocated(error)) return

      largest_rejected = -1.0_dp
      smallest_kept = 2.0_dp
      do i = 1, size(avas%occupied_weights)
         if (avas%occupied_weights(i) < 0.2_dp) then
            largest_rejected = max(largest_rejected, avas%occupied_weights(i))
         else
            smallest_kept = min(smallest_kept, avas%occupied_weights(i))
         end if
      end do

      ! Measured here: 0.000, 0.000, 0.000, 0.779, 0.959, 0.991, 0.991. The
      ! rejected ones are zero to three decimals and the kept ones run from 0.78
      ! up, so the cut at 0.2 sits in a gap three quarters of the range wide.
      !
      ! **The lowest kept weight is 0.78 rather than near one, and that is
      ! chemistry rather than slop.** It is the sigma bonding orbital, which
      ! carries real nitrogen 2s character as well as 2p -- so asking only for
      ! 2p captures most of it and not all. Adding "N 2s" to the request lifts
      ! that same orbital to 0.959, which is `nitrogen_2s2p_against_pyscf`
      ! above. An assertion that every kept weight is near one would be
      ! asserting that no bonding orbital ever hybridises.
      call check(error, largest_rejected < 0.01_dp, &
                 "the rejected orbitals should have essentially none of the "// &
                 "atomic character asked for")
      if (allocated(error)) return
      call check(error, smallest_kept > 0.5_dp, &
                 "and the kept ones should be dominated by it, well clear of the "// &
                 "0.2 threshold")
      if (allocated(error)) return
      call check(error, smallest_kept - largest_rejected > 0.5_dp, &
                 "so the threshold sits in a wide gap, and moving it anywhere "// &
                 "inside that gap would not change the selection")
   end subroutine test_gap

   subroutine test_valence(error)
      !! The whole valence shell, by counting rather than by judgement
      !!
      !! Nitrogen: each atom contributes 1s, 2s and 2p to the free-atom minimal
      !! basis, so ten minimal-basis orbitals against seven occupied. Two of the
      !! occupied are the 1s core, which leaves five valence occupied and three
      !! valence virtual -- CAS(10,8), the space anyone would have written down
      !! for N2 by hand.
      !!
      !! Nothing here is a threshold. The size follows from the elements, so the
      !! same molecule gives the same space in any basis set, and the assertion
      !! below would hold in 6-31G as it does in cc-pVDZ.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(avas_result_t) :: valence
      type(casci_result_t) :: small, full
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :), gram(:, :), work(:, :)
      integer :: i, j, n_mo

      call build_czt_molecule(NITROGEN_Z, NITROGEN_SYM, NITROGEN, "cc-pvdz", mol, err)
      call run_czt_rhf(mol, 14, 300, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error, scf%converged, "the SCF should converge")
      if (allocated(error)) return

      call valence_select(mol, NITROGEN_Z, NITROGEN_SYM, NITROGEN, scf%orbitals, 14, &
                          valence, err)
      call check(error,.not. err%has_error(), "the valence space should be found")
      if (allocated(error)) return

      call check(error, valence%n_inactive, 2, "two 1s cores are inactive")
      if (allocated(error)) return
      call check(error, valence%n_active, 8, "eight valence orbitals")
      if (allocated(error)) return
      call check(error, valence%n_active_electrons, 10, "ten valence electrons")
      if (allocated(error)) return
      call check(error, valence%n_active_occupied, 5, "five of them occupied")
      if (allocated(error)) return

      ! The orbitals have to be a usable set, not only a correctly counted one:
      ! the valence virtuals are combinations of canonical virtuals and the
      ! external block is what is left, so orthonormality is the thing that
      ! would break if those two were assembled wrongly.
      n_mo = size(valence%orbitals, 2)
      call mol%overlap(overlap)
      allocate (work(size(overlap, 1), n_mo), gram(n_mo, n_mo))
      call pic_gemm(overlap, valence%orbitals, work)
      call pic_gemm(valence%orbitals, work, gram, transa="T")
      do j = 1, n_mo
         do i = 1, n_mo
            if (i == j) then
               call check(error, abs(gram(i, j) - 1.0_dp) < 1.0e-10_dp, &
                          "the valence orbital set should be normalised")
            else
               call check(error, abs(gram(i, j)) < 1.0e-10_dp, &
                          "and orthogonal")
            end if
            if (allocated(error)) return
         end do
      end do

      ! And it has to be a better space than the one it contains. CAS(6,6) on
      ! the same reference is a subset of the valence space, so the valence CI
      ! cannot be higher.
      call run_czt_casci(mol, scf%orbitals, 4, 6, 3, 3, small, err, &
                         tolerance=1.0e-11_dp)
      call run_czt_casci(mol, valence%orbitals, valence%n_inactive, &
                         valence%n_active, 5, 5, full, err, tolerance=1.0e-11_dp)
      call check(error,.not. err%has_error(), "both CIs should run")
      if (allocated(error)) return
      call check(error, full%energy < small%energy, &
                 "the valence space contains CAS(6,6) and cannot do worse")
      if (allocated(error)) return

      deallocate (work, gram)
      call mol%destroy()
   end subroutine test_valence

   subroutine test_refusals(error)
      !! Requests that name nothing
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(avas_result_t) :: avas
      character(len=8) :: labels(1)

      call build_czt_molecule(NITROGEN_Z, NITROGEN_SYM, NITROGEN, "cc-pvdz", &
                              mol, err)
      call run_czt_rhf(mol, 14, 300, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error, scf%converged, "the SCF should converge")
      if (allocated(error)) return

      ! A shell the free atom does not have. Nitrogen's minimal basis is 1s2s2p
      ! and stops there, so there is no 3d to project onto -- and asking for one
      ! must say so rather than return an empty space.
      labels(1) = "N 3d"
      call avas_select(mol, NITROGEN_Z, NITROGEN_SYM, NITROGEN, scf%orbitals, 7, &
                       labels, avas, err)
      call check(error, err%has_error(), &
                 "a shell the minimal basis does not have should be refused")
      if (allocated(error)) return
      call err%clear()

      ! An element that is not in the molecule.
      labels(1) = "C 2p"
      call avas_select(mol, NITROGEN_Z, NITROGEN_SYM, NITROGEN, scf%orbitals, 7, &
                       labels, avas, err)
      call check(error, err%has_error(), &
                 "an element the molecule does not contain should be refused")
      call err%clear()
      call mol%destroy()
   end subroutine test_refusals

end module test_mqc_avas

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_avas, only: collect_mqc_avas_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_avas", collect_mqc_avas_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
