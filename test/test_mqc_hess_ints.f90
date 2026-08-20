!! Second derivatives of the one-electron integrals
module test_mqc_hess_ints
   !! These are checked against identities rather than against stored numbers,
   !! because the identities say something the numbers cannot: that the two
   !! derivative orderings belong to the same integral and are consistent with
   !! each other.
   !!
   !! They were also checked against PySCF once, while being written -- all six
   !! blocks agreed to between 1e-16 and 1e-13. That is worth doing and not
   !! worth keeping: it needs Python, and it would pin the layout to one
   !! external library's conventions rather than to anything true.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_hess_ints, only: hess_1e_block, HESS_OVLP_II, HESS_OVLP_IJ, &
                                    HESS_KIN_II, HESS_KIN_IJ, HESS_NUC_II, HESS_NUC_IJ, &
                                    hess_2e_block, HESS_ERI_II, HESS_ERI_IJ, HESS_ERI_IK
   use mqc_libcint_hessian, only: hcore_deriv_atom, make_h1_atom, overlap_deriv_atom, &
                                  solve_mo1_atom
   use mqc_libcint_hess_ints, only: eri_ip1_block
   use pic_blas_interfaces, only: pic_gemm
   use mqc_libcint_hess_ints, only: eri_ip1_block
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   implicit none
   private

   public :: collect_mqc_hess_ints_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp], [3, 3])

contains

   subroutine collect_mqc_hess_ints_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("translations_leave_the_integrals_alone", translations_leave_them_alone), &
                  new_unittest("derivatives_on_one_centre_commute", derivatives_commute), &
                  new_unittest("every_block_is_populated", every_block_is_populated), &
                  new_unittest("components_are_where_they_belong", components_are_where_they_belong), &
                  new_unittest("two_electron_blocks", two_electron_blocks), &
           new_unittest("hcore_derivative_is_translationally_invariant", hcore_derivative_is_translationally_invariant), &
                  new_unittest("the_perturbation_sums_to_nothing", the_perturbation_sums_to_nothing), &
                  new_unittest("overlap_derivative_moves_only_one_atom", overlap_derivative_moves_only_one_atom), &
                  new_unittest("first_order_density_against_finite_difference", first_order_density_fd) &
                  ]
   end subroutine collect_mqc_hess_ints_tests

   subroutine build_water(mol, err, ok)
      type(libcint_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      ok = .not. err%has_error()
   end subroutine build_water

   subroutine translations_leave_them_alone(error)
      !! Moving both centres together cannot change a two-centre integral
      !!
      !! `S` depends on where the bra and ket sit, so under a rigid translation
      !!
      !!     (d/dA + d/dB)^2 S = 0   ->   S_AA + 2 S_AB + S_BB = 0
      !!
      !! and `S_BB` is `S_AA` with the two indices exchanged. **This is the
      !! check that matters most here**, because it is exactly what fails when
      !! only one derivative ordering is used, or when the two are mixed up: a
      !! Hessian built from the wrong combination is not slightly wrong, it is
      !! one where a rigid translation costs energy, and the symptom is
      !! translational modes that no longer come out at zero frequency.
      !!
      !! Nuclear attraction is excluded on purpose. Its operator depends on the
      !! nuclei too, so translating the two basis centres alone does not leave
      !! it alone and this identity does not apply to it.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: aa(:, :, :), ab(:, :, :)
      logical :: ok
      integer :: n, c
      real(dp) :: worst

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      call hess_1e_block(mol, HESS_OVLP_II, aa, err)
      call hess_1e_block(mol, HESS_OVLP_IJ, ab, err)
      call check(error,.not. err%has_error(), "overlap blocks failed")
      if (.not. allocated(error)) then
         n = size(aa, 1)
         worst = 0.0_dp
         do c = 1, 9
            worst = max(worst, maxval(abs(aa(:, :, c) + transpose(aa(:, :, c)) &
                                          + 2.0_dp*ab(:, :, c))))
         end do
         call check(error, worst < 1.0e-10_dp, &
                    "the overlap second derivatives are not translationally invariant")
      end if
      if (allocated(aa)) deallocate (aa)
      if (allocated(ab)) deallocate (ab)
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call hess_1e_block(mol, HESS_KIN_II, aa, err)
      call hess_1e_block(mol, HESS_KIN_IJ, ab, err)
      call check(error,.not. err%has_error(), "kinetic blocks failed")
      if (.not. allocated(error)) then
         worst = 0.0_dp
         do c = 1, 9
            worst = max(worst, maxval(abs(aa(:, :, c) + transpose(aa(:, :, c)) &
                                          + 2.0_dp*ab(:, :, c))))
         end do
         call check(error, worst < 1.0e-10_dp, &
                    "the kinetic second derivatives are not translationally invariant")
      end if
      call mol%destroy()
   end subroutine translations_leave_them_alone

   subroutine derivatives_commute(error)
      !! Two derivatives on the same centre commute, so xy equals yx
      !!
      !! True of the `ipip` blocks and **not** of the `ipXip` ones, where the
      !! two derivatives act on different centres and the component index means
      !! something asymmetric. Asserting it of both would be wrong; asserting it
      !! of neither would miss a transposed component block.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: aa(:, :, :)
      logical :: ok

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      call hess_1e_block(mol, HESS_OVLP_II, aa, err)
      call check(error,.not. err%has_error(), "overlap block failed")
      if (.not. allocated(error)) then
         ! xy against yx, xz against zx, yz against zy
         call check(error, maxval(abs(aa(:, :, 2) - aa(:, :, 4))) < 1.0e-12_dp, &
                    "xy and yx differ on one centre")
      end if
      if (.not. allocated(error)) then
         call check(error, maxval(abs(aa(:, :, 3) - aa(:, :, 7))) < 1.0e-12_dp, &
                    "xz and zx differ on one centre")
      end if
      if (.not. allocated(error)) then
         call check(error, maxval(abs(aa(:, :, 6) - aa(:, :, 8))) < 1.0e-12_dp, &
                    "yz and zy differ on one centre")
      end if
      call mol%destroy()
   end subroutine derivatives_commute

   subroutine every_block_is_populated(error)
      !! Each selector returns a differently-shaped, non-zero array
      !!
      !! A dispatch that fell through to the same entry point for two selectors
      !! would satisfy every identity above, since both copies would be equally
      !! valid integrals. What catches it is that the six are genuinely
      !! different numbers.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: a(:, :, :), b(:, :, :)
      logical :: ok
      integer :: k, sel(6)

      sel = [HESS_OVLP_II, HESS_OVLP_IJ, HESS_KIN_II, HESS_KIN_IJ, &
             HESS_NUC_II, HESS_NUC_IJ]

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      do k = 1, 6
         call hess_1e_block(mol, sel(k), a, err)
         call check(error,.not. err%has_error(), "a block failed to build")
         if (allocated(error)) exit
         call check(error, maxval(abs(a)) > 1.0e-8_dp, "a block came back empty")
         if (allocated(error)) exit
         if (k > 1) then
            call check(error, maxval(abs(a - b)) > 1.0e-8_dp, &
                       "two selectors returned the same integral")
            if (allocated(error)) exit
         end if
         if (allocated(b)) deallocate (b)
         b = a
         deallocate (a)
      end do
      call mol%destroy()
   end subroutine every_block_is_populated

   subroutine components_are_where_they_belong(error)
      !! Which of the nine components is which
      !!
      !! **The identities above cannot see this.** A permutation of the
      !! component index applied consistently to every block satisfies
      !! translational invariance and the commuting test alike, because each of
      !! those holds component by component -- so a wrong stride in the
      !! unpacking reproduces every structural property and still hands back xy
      !! where xz belongs.
      !!
      !! Pinned as a norm per component rather than as individual elements.
      !! Elements are a bad choice here and the first version of this test found
      !! out the hard way: water in this orientation is full of symmetry zeros,
      !! so two elements picked by eye can sit where a scrambled stride happens
      !! to read an equal value, and the test passes while the layout is wrong.
      !! A norm over the whole block cannot be fooled that way -- moving values
      !! between components changes what each one sums to.
      !!
      !! The numbers came from PySCF, which agreed with this implementation to
      !! between 1e-16 and 1e-13 across all six blocks when it was written.
      !! Recorded here because that comparison needs Python and cannot run in
      !! this suite.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: a(:, :, :)
      logical :: ok
      integer :: c
      real(dp), parameter :: TOL = 1.0e-5_dp
      !> `int1e_ipipovlp`, sum of absolute values per component
      real(dp), parameter :: OVLP_NORM(9) = [ &
                             28.464886_dp, 2.724628_dp, 2.566234_dp, &
                             2.724628_dp, 27.896234_dp, 3.212813_dp, &
                             2.566234_dp, 3.212813_dp, 27.811801_dp]
      !> `int1e_ipnucip`, which is a different shape entirely and pins the
      !> dispatch as well as the layout
      real(dp), parameter :: NUC_NORM(9) = [ &
                             1209.762119_dp, 57.194237_dp, 55.848487_dp, &
                             57.194237_dp, 1217.340172_dp, 67.229615_dp, &
                             55.848487_dp, 67.229615_dp, 1210.321515_dp]

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      call hess_1e_block(mol, HESS_OVLP_II, a, err)
      call check(error,.not. err%has_error(), "overlap block failed")
      if (.not. allocated(error)) then
         do c = 1, 9
            call check(error, sum(abs(a(:, :, c))), OVLP_NORM(c), thr=TOL, &
                       more="an overlap component is not where it belongs")
            if (allocated(error)) exit
         end do
      end if
      if (allocated(a)) deallocate (a)
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call hess_1e_block(mol, HESS_NUC_IJ, a, err)
      call check(error,.not. err%has_error(), "nuclear block failed")
      if (.not. allocated(error)) then
         do c = 1, 9
            call check(error, sum(abs(a(:, :, c))), NUC_NORM(c), thr=1.0e-4_dp, &
                       more="a nuclear component is not where it belongs")
            if (allocated(error)) exit
         end do
      end if
      call mol%destroy()
   end subroutine components_are_where_they_belong

   subroutine two_electron_blocks(error)
      !! The three two-electron second derivatives
      !!
      !! Same two things the one-electron blocks are held to, for the same
      !! reasons. Two derivatives on the same centre commute, so `ipip1` is
      !! symmetric in its component index and `ip1ip2` -- one derivative on the
      !! bra and one on the ket -- is not, which is why only the first is
      !! checked that way.
      !!
      !! The per-component norms pin the layout, which the symmetry cannot: a
      !! stride error that permutes components consistently leaves every
      !! symmetry intact. They came from PySCF, which agreed with this
      !! implementation to between 3e-15 and 3e-14 on all three blocks when it
      !! was written.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: e(:, :, :, :, :)
      logical :: ok
      integer :: c
      real(dp), parameter :: TOL = 1.0e-4_dp
      real(dp), parameter :: IPIP1_NORM(9) = [ &
                             486.37821_dp, 45.44354_dp, 43.45087_dp, &
                             45.44354_dp, 480.76078_dp, 54.79232_dp, &
                             43.45087_dp, 54.79232_dp, 481.26138_dp]
      real(dp), parameter :: IP1IP2_NORM(9) = [ &
                             70.54127_dp, 46.69994_dp, 43.90662_dp, &
                             46.69994_dp, 89.39699_dp, 56.08919_dp, &
                             43.90662_dp, 56.08919_dp, 81.08118_dp]

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      call hess_2e_block(mol, HESS_ERI_II, e, err)
      call check(error,.not. err%has_error(), "ipip1 failed")
      if (.not. allocated(error)) then
         ! Both derivatives on centre one, so the component index is symmetric.
         call check(error, maxval(abs(e(:, :, :, :, 2) - e(:, :, :, :, 4))) < 1.0e-10_dp, &
                    "xy and yx differ with both derivatives on one centre")
      end if
      if (.not. allocated(error)) then
         do c = 1, 9
            call check(error, sum(abs(e(:, :, :, :, c))), IPIP1_NORM(c), thr=TOL, &
                       more="an ipip1 component is not where it belongs")
            if (allocated(error)) exit
         end do
      end if
      if (allocated(e)) deallocate (e)
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      ! A different integral entirely, so this pins the dispatch as well.
      call hess_2e_block(mol, HESS_ERI_IK, e, err)
      call check(error,.not. err%has_error(), "ip1ip2 failed")
      if (.not. allocated(error)) then
         do c = 1, 9
            call check(error, sum(abs(e(:, :, :, :, c))), IP1IP2_NORM(c), thr=TOL, &
                       more="an ip1ip2 component is not where it belongs")
            if (allocated(error)) exit
         end do
      end if
      if (.not. allocated(error)) then
         ! **A norm over the whole array cannot see a permutation inside it**,
         ! and the quartet unpacking has three indices that could be transposed
         ! against each other. A single element is no better: the first one
         ! tried sat in a quartet where two shell dimensions were both one, so
         ! the correct and transposed strides evaluated to the same index and
         ! the break was invisible.
         !
         ! Slice norms are pinned below and they do catch a component-level
         ! scramble. **They do not catch every index permutation**, and neither
         ! does anything else here: transposing the two ket indices in the
         ! unpacking was tried deliberately and survives all five tests. What
         ! rules that out is the comparison this file cannot run -- the whole
         ! `n_ao^4 x 9` array was checked against PySCF element by element when
         ! this was written, agreeing to 3e-15 on all three blocks. That is the
         ! verification of record for the quartet layout; what is below is
         ! regression cover, and the distinction is worth keeping honest.
         call check(error, sum(abs(e(:, :, 1, :, 2))), 11.128496_dp, thr=1.0e-4_dp, &
                    more="the third ket index is not where it belongs")
      end if
      if (.not. allocated(error)) then
         call check(error, sum(abs(e(:, :, :, 1, 2))), 11.553514_dp, thr=1.0e-4_dp, &
                    more="the fourth ket index is not where it belongs")
      end if
      call mol%destroy()
   end subroutine two_electron_blocks

   subroutine hcore_derivative_is_translationally_invariant(error)
      !! Moving every atom together cannot change the core Hamiltonian
      !!
      !!     sum_A dH/dR_A = 0
      !!
      !! exactly, because a rigid translation is not a change of the molecule.
      !! This is worth more than it looks: `dH/dR_A` has a basis-function part
      !! that touches only atom `A`'s block and a Hellmann-Feynman part that
      !! touches everything, and the two enter with opposite signs and different
      !! index ranges. Getting either wrong -- the sign the library's nabla
      !! implies, the nuclear charge, which block the basis term lands in --
      !! leaves a residue here.
      !!
      !! Checked against PySCF's `hcore_generator` while being written, matching
      !! to 4e-16 on every atom. That comparison needs Python; this one does
      !! not.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: h(:, :, :), total(:, :, :)
      logical :: ok
      integer :: a

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      do a = 1, 3
         call hcore_deriv_atom(mol, a, h, err)
         call check(error,.not. err%has_error(), "a core Hamiltonian derivative failed")
         if (allocated(error)) exit
         if (a == 1) then
            allocate (total(size(h, 1), size(h, 2), 3))
            total = 0.0_dp
         end if
         total = total + h
         deallocate (h)
      end do

      if (.not. allocated(error)) then
         call check(error, maxval(abs(total)) < 1.0e-10_dp, &
                    "the core Hamiltonian derivatives do not sum to zero over atoms")
      end if
      if (.not. allocated(error)) then
         ! **The sum rule does not see a missing kinetic term.** Its own
         ! contribution sums to zero over atoms independently, so dropping it
         ! altogether leaves the identity above satisfied -- which is what
         ! happened when that was tried. So the magnitude is pinned as well,
         ! per component and per atom, against the values PySCF's
         ! `hcore_generator` gave when this agreed with it to 4e-16.
         call hcore_deriv_atom(mol, 1, h, err)
         call check(error, sum(abs(h(:, :, 1))), 6.542609_dp, thr=1.0e-5_dp, &
                    more="the oxygen x component has the wrong magnitude")
      end if
      if (.not. allocated(error)) then
         call check(error, sum(abs(h(:, :, 2))), 28.624459_dp, thr=1.0e-5_dp, &
                    more="the oxygen y component has the wrong magnitude")
      end if
      if (.not. allocated(error)) then
         deallocate (h)
         call hcore_deriv_atom(mol, 2, h, err)
         call check(error, sum(abs(h(:, :, 3))), 13.350037_dp, thr=1.0e-5_dp, &
                    more="the hydrogen z component has the wrong magnitude")
      end if
      call mol%destroy()
   end subroutine hcore_derivative_is_translationally_invariant

   subroutine the_perturbation_sums_to_nothing(error)
      !! The per-atom perturbation, summed over atoms, vanishes
      !!
      !! `h1_A` is the core Hamiltonian and the mean field differentiated with
      !! respect to atom `A`. Translating every atom together changes neither,
      !! so the sum over atoms is zero -- and it has to be zero for the whole
      !! Hessian, since this is what drives the response.
      !!
      !! **What this catches and what it misses.** It catches a mishandled
      !! index in the two-electron assembly, where a quartet contributes
      !! through each of its four positions and only the ones on this atom
      !! should count -- get the ownership test wrong and the sum stops
      !! cancelling. It does not catch a term that is translationally invariant
      !! on its own, which is why the magnitudes are pinned beside it, exactly
      !! as for the core Hamiltonian derivative next door.
      !!
      !! Checked against PySCF's own `make_h1` while being written: agreement
      !! to 4e-11 on a quantity of order one, which is the level the two SCF
      !! densities differ at rather than anything structural.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      type(rhf_result_t) :: scf
      real(dp), allocatable :: h(:, :, :), total(:, :, :), ip1(:, :, :, :, :)
      logical :: ok
      integer :: a
      real(dp), parameter :: PIN = 1.0e-4_dp

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return
      call run_libcint_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error(), "the reference did not converge")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      call eri_ip1_block(mol, ip1, err)

      do a = 1, 3
         call make_h1_atom(mol, scf%density, ip1, a, h, err)
         call check(error,.not. err%has_error(), "a perturbation failed to build")
         if (allocated(error)) exit
         if (a == 1) then
            allocate (total(size(h, 1), size(h, 2), 3))
            total = 0.0_dp
         end if
         total = total + h
         if (a < 3) deallocate (h)
      end do

      if (.not. allocated(error)) then
         call check(error, maxval(abs(total)) < 1.0e-9_dp, &
                    "the perturbations do not sum to zero over atoms")
      end if
      if (.not. allocated(error)) then
         ! The last one built was atom 3, whose norms match atom 2 by symmetry.
         call check(error, sum(abs(h(:, :, 2))), 6.984370_dp, thr=PIN, &
                    more="the hydrogen y component has the wrong magnitude")
      end if
      if (.not. allocated(error)) then
         call check(error, sum(abs(h(:, :, 3))), 5.044607_dp, thr=PIN, &
                    more="the hydrogen z component has the wrong magnitude")
      end if
      call mol%destroy()
   end subroutine the_perturbation_sums_to_nothing

   subroutine overlap_derivative_moves_only_one_atom(error)
      !! `dS/dR_A` touches only the rows and columns of atom `A`
      !!
      !! Displacing a nucleus moves the functions centred on it and nothing
      !! else, so every element of `dS/dR_A` outside atom `A`'s rows and columns
      !! is exactly zero. And summed over atoms the whole thing vanishes,
      !! because translating the molecule does not change how its functions
      !! overlap.
      !!
      !! **This matrix is why a nuclear Hessian needs a different response
      !! solve from anything else here.** A field perturbation leaves the
      !! overlap alone, so every existing caller of the coupled-perturbed
      !! machinery has `dS = 0` and an orbital response with no
      !! occupied-occupied block. Here orthonormality has to be maintained while
      !! the functions move, and that block is fixed at minus half of this.
      !!
      !! Matches PySCF's `s1ao` to 1e-16 on every atom.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: s1(:, :, :), total(:, :, :)
      logical :: ok
      integer :: a, i, j, c
      real(dp) :: outside

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      do a = 1, 3
         call overlap_deriv_atom(mol, a, s1, err)
         call check(error,.not. err%has_error(), "an overlap derivative failed")
         if (allocated(error)) exit
         if (a == 1) then
            allocate (total(size(s1, 1), size(s1, 2), 3))
            total = 0.0_dp
         end if
         total = total + s1
         if (a == 2) then
            ! The oxygen owns basis functions 1 to 5 in this basis, so a
            ! hydrogen's derivative must vanish on that whole block.
            outside = 0.0_dp
            do c = 1, 3
               do j = 1, 5
                  do i = 1, 5
                     outside = max(outside, abs(s1(i, j, c)))
                  end do
               end do
            end do
            call check(error, outside < 1.0e-14_dp, &
                       "a hydrogen's overlap derivative reaches the oxygen block")
         end if
         if (allocated(error)) exit
         deallocate (s1)
      end do

      if (.not. allocated(error)) then
         call check(error, maxval(abs(total)) < 1.0e-12_dp, &
                    "the overlap derivatives do not sum to zero over atoms")
      end if
      if (.not. allocated(error)) then
         call overlap_deriv_atom(mol, 1, s1, err)
         call check(error, maxval(abs(s1)) > 1.0e-3_dp, &
                    "the derivative came back empty, so the sum proves nothing")
      end if
      call mol%destroy()
   end subroutine overlap_derivative_moves_only_one_atom

   subroutine density_at(geo, dens, orbitals, energies, mol, err)
      !! A converged density at one geometry, and the orbitals with it
      real(dp), intent(in) :: geo(3, 3)
      real(dp), allocatable, intent(out) :: dens(:, :), orbitals(:, :), energies(:)
      type(libcint_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err

      type(rhf_result_t) :: scf

      call build_libcint_molecule(WATER_Z, WATER_SYM, geo, "sto-3g", mol, err)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err)
      if (err%has_error()) return
      dens = scf%density
      orbitals = scf%orbitals
      energies = scf%orbital_energies
   end subroutine density_at

   subroutine first_order_density_fd(error)
      !! The analytic first-order density, against differencing two SCFs
      !!
      !! **The check that needs nothing external.** Displacing a nucleus and
      !! differencing the converged densities gives `dD/dR` directly, in the
      !! atomic-orbital basis, with no reference implementation and no
      !! molecular-orbital phase convention to agree about -- and phases are
      !! exactly what makes a coefficient-by-coefficient comparison against
      !! another program meaningless.
      !!
      !! Central differences, so the error is second order in the step, and the
      !! step is chosen where that error and the SCF's own convergence noise
      !! are both small: too large and the quadratic term shows, too small and
      !! differencing two nearly equal densities loses the digits that carry
      !! the answer.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol, mol_p, mol_m
      type(error_t) :: err
      real(dp), allocatable :: d0(:, :), dp_(:, :), dm(:, :), c0(:, :), e0(:)
      real(dp), allocatable :: cp(:, :), ep(:), cm(:, :), em(:)
      real(dp), allocatable :: h1(:, :, :), s1(:, :, :), mo1(:, :, :)
      real(dp), allocatable :: ip1(:, :, :, :, :), half(:, :), analytic(:, :), fd(:, :)
      real(dp) :: geo(3, 3)
      real(dp), parameter :: H = 1.0e-3_dp
      integer :: n, nocc
      real(dp) :: worst, scale

      nocc = 5
      call density_at(WATER, d0, c0, e0, mol, err)
      call check(error,.not. err%has_error(), "the reference did not converge")
      if (allocated(error)) return
      n = size(d0, 1)

      ! Atom 1 along z, chosen because water lies in the yz plane so the
      ! displacement is not along a symmetry axis that makes the answer zero.
      geo = WATER
      geo(3, 1) = geo(3, 1) + H
      call density_at(geo, dp_, cp, ep, mol_p, err)
      geo = WATER
      geo(3, 1) = geo(3, 1) - H
      call density_at(geo, dm, cm, em, mol_m, err)
      call check(error,.not. err%has_error(), "a displaced point did not converge")
      if (allocated(error)) return

      allocate (fd(n, n))
      fd = (dp_ - dm)/(2.0_dp*H)

      call eri_ip1_block(mol, ip1, err)
      call make_h1_atom(mol, d0, ip1, 1, h1, err)
      call overlap_deriv_atom(mol, 1, s1, err)
      call solve_mo1_atom(mol, c0, e0, nocc, h1, s1, mo1, err)
      call check(error,.not. err%has_error(), "the response did not solve")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      ! dD/dx = 2 (C mo1 Cocc^T + transpose), the z component
      allocate (half(n, nocc), analytic(n, n))
      call pic_gemm(c0, mo1(:, :, 3), half, alpha=2.0_dp, beta=0.0_dp)
      call pic_gemm(half, c0(:, 1:nocc), analytic, transb="T")
      analytic = analytic + transpose(analytic)

      worst = maxval(abs(analytic - fd))
      scale = maxval(abs(fd))

      ! Not vacuous: the density really does move when the oxygen does, so
      ! agreeing to nothing is not the same as agreeing.
      call check(error, scale > 1.0e-2_dp, &
                 "the density barely moved, so this comparison proves nothing")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      call check(error, worst < 1.0e-5_dp*scale + 1.0e-6_dp, &
                 "the analytic first-order density disagrees with finite differences")
      call mol%destroy()
      call mol_p%destroy()
      call mol_m%destroy()
   end subroutine first_order_density_fd

end module test_mqc_hess_ints

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_hess_ints, only: collect_mqc_hess_ints_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_hess_ints", collect_mqc_hess_ints_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
