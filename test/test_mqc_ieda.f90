!! Regrouping the kinetic energy onto atoms and atom pairs
module test_mqc_ieda
   !! The decomposition is a regrouping of `sum_pq gamma_pq T_pq`, so its
   !! defining property is that nothing is lost: the atomic and atom-pair pieces
   !! must add back up to the sum they came from. `kinetic_decomposition`
   !! enforces that itself and refuses to return otherwise, so every test here
   !! that gets an answer at all has already exercised the guard.
   !!
   !! What the guard cannot check is whether the *split* is right -- a
   !! consistent decomposition that puts energy on the wrong atom still sums
   !! correctly. So the numbers below are worked out by hand from a four-orbital
   !! matrix small enough to add up by eye, and compared element by element.
   !!
   !! **The factor of two is the thing being pinned.** An atom pair collects
   !! both `p in A, q in B` and `p in B, q in A`, and `inter` reports their sum,
   !! not one of them. Getting that wrong halves every interatomic number while
   !! leaving a table that still looks entirely reasonable; it is how the first
   !! version of this module was written, and the sum rule is what caught it.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_quao, only: quao_result_t
   use mqc_libcint_ieda, only: kinetic_decomposition, kinetic_total, &
                               nuclear_attraction_per_atom, nuclear_decomposition
   implicit none
   private

   public :: collect_mqc_ieda_tests

   !> Two atoms, two orbitals each, symmetric as a real interference matrix must
   !> be. Chosen so every block sums to a round number:
   !>
   !>     intra(1) = 1.0 + 0.2 + 0.2 + 2.0            =  3.4
   !>     intra(2) = 3.0 + 0.6 + 0.6 + 4.0            =  8.2
   !>     one-way  = -0.5 - 0.1 - 0.3 - 0.4           = -1.3
   !>     inter    = 2 * one-way                      = -2.6
   !>     total    = 3.4 + 8.2 - 1.3 - 1.3            =  9.0
   integer, parameter :: N_ORB = 4
   integer, parameter :: ATOM_OF(N_ORB) = [1, 1, 2, 2]
   real(dp), parameter :: T_INTF(N_ORB, N_ORB) = reshape([ &
                                                         1.0_dp, 0.2_dp, -0.5_dp, -0.1_dp, &
                                                         0.2_dp, 2.0_dp, -0.3_dp, -0.4_dp, &
                                                         -0.5_dp, -0.3_dp, 3.0_dp, 0.6_dp, &
                                                         -0.1_dp, -0.4_dp, 0.6_dp, 4.0_dp], [N_ORB, N_ORB])
   real(dp), parameter :: TOL = 1.0e-12_dp

   !> Hydrogen, 1.4 Bohr apart. Two identical nuclei, so which nucleus a
   !> per-nucleus integral belongs to can be read straight off the numbers
   !> without dividing out a charge.
   integer, parameter :: H2_Z(2) = [1, 1]
   character(len=2), parameter :: H2_SYM(2) = ["H ", "H "]
   real(dp), parameter :: H2(3, 2) = reshape([ &
                                             0.0_dp, 0.0_dp, 0.0_dp, &
                                             0.0_dp, 0.0_dp, 1.4_dp], [3, 2])

   !> A flat density over the same two atoms, so every block sums to four and
   !> the grouping can be worked out in the margin:
   !>
   !>     intra(1) = 4 * 1                                =  4   (nucleus 1)
   !>     intra(2) = 4 * 2                                =  8   (nucleus 2)
   !>     inter    = 4*2 + 4*1  (own density, foreign nucleus)
   !>              + 2 * (4*1 + 4*2)  (interference, both orders, both nuclei)
   !>              = 12 + 24                              = 36
   !>     total    = 4 + 8 + 36                           = 48 = 3 * sum(D)
   real(dp), parameter :: FLAT(N_ORB, N_ORB) = 1.0_dp

contains

   subroutine collect_mqc_ieda_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("kinetic_pieces_add_back_up", kinetic_pieces_add_back_up), &
                  new_unittest("atom_terms_are_the_right_blocks", atom_terms_are_the_right_blocks), &
                  new_unittest("pair_energy_counts_both_orderings", pair_energy_counts_both_orderings), &
                  new_unittest("inter_is_symmetric_with_no_diagonal", inter_is_symmetric_with_no_diagonal), &
                  new_unittest("rejects_an_orbital_off_the_molecule", rejects_an_orbital_off_the_molecule), &
                  new_unittest("rejects_a_matrix_of_the_wrong_size", rejects_a_matrix_of_the_wrong_size), &
                  new_unittest("nuclear_terms_land_in_the_right_bins", nuclear_terms_land_in_the_right_bins), &
                  new_unittest("nuclear_pieces_add_back_up", nuclear_pieces_add_back_up), &
                  new_unittest("rejects_a_partial_nuclear_field", rejects_a_partial_nuclear_field), &
                  new_unittest("attraction_is_on_the_right_nucleus", attraction_is_on_the_right_nucleus) &
                  ]
   end subroutine collect_mqc_ieda_tests

   subroutine two_atom_case(quao)
      !! The fixture: `kinetic_decomposition` reads nothing else from the type
      type(quao_result_t), intent(out) :: quao

      quao%n_quao = N_ORB
      allocate (quao%atom_of(N_ORB))
      quao%atom_of = ATOM_OF
   end subroutine two_atom_case

   subroutine kinetic_pieces_add_back_up(error)
      !! The sum rule, checked from outside as well as in
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp), allocatable :: intra(:), inter(:, :)

      call two_atom_case(quao)
      call kinetic_decomposition(quao, T_INTF, 2, intra, inter, err)
      call check(error,.not. err%has_error(), "the decomposition refused a valid case")
      if (allocated(error)) return

      call check(error, kinetic_total(intra, inter), sum(T_INTF), thr=TOL, &
                 more="the pieces do not add up to the kinetic energy")
      if (allocated(error)) return
      call check(error, kinetic_total(intra, inter), 9.0_dp, thr=TOL, &
                 more="the total is not the hand-computed 9.0")
   end subroutine kinetic_pieces_add_back_up

   subroutine atom_terms_are_the_right_blocks(error)
      !! Summing correctly is not enough; the energy must land on the right atom
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp), allocatable :: intra(:), inter(:, :)

      call two_atom_case(quao)
      call kinetic_decomposition(quao, T_INTF, 2, intra, inter, err)
      call check(error,.not. err%has_error(), "the decomposition refused a valid case")
      if (allocated(error)) return

      call check(error, intra(1), 3.4_dp, thr=TOL, more="atom 1 has the wrong block")
      if (allocated(error)) return
      call check(error, intra(2), 8.2_dp, thr=TOL, more="atom 2 has the wrong block")
   end subroutine atom_terms_are_the_right_blocks

   subroutine pair_energy_counts_both_orderings(error)
      !! `inter` is the whole pair, not half of it
      !!
      !! The one-way sum is -1.3 and the pair energy is -2.6. Asserting the
      !! latter is what makes a silent halving of every interatomic number fail
      !! here rather than in a table someone reads six months later.
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp), allocatable :: intra(:), inter(:, :)
      real(dp) :: one_way

      call two_atom_case(quao)
      call kinetic_decomposition(quao, T_INTF, 2, intra, inter, err)
      call check(error,.not. err%has_error(), "the decomposition refused a valid case")
      if (allocated(error)) return

      one_way = sum(T_INTF(1:2, 3:4))
      call check(error, one_way, -1.3_dp, thr=TOL, more="the fixture changed")
      if (allocated(error)) return
      call check(error, inter(1, 2), 2.0_dp*one_way, thr=TOL, &
                 more="the pair energy is not both orderings summed")
   end subroutine pair_energy_counts_both_orderings

   subroutine inter_is_symmetric_with_no_diagonal(error)
      !! Intra-atomic energy belongs in `intra` and must not also sit in `inter`
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp), allocatable :: intra(:), inter(:, :)

      call two_atom_case(quao)
      call kinetic_decomposition(quao, T_INTF, 2, intra, inter, err)
      call check(error,.not. err%has_error(), "the decomposition refused a valid case")
      if (allocated(error)) return

      call check(error, inter(1, 1), 0.0_dp, thr=TOL, more="inter has a diagonal")
      if (allocated(error)) return
      call check(error, inter(2, 2), 0.0_dp, thr=TOL, more="inter has a diagonal")
      if (allocated(error)) return
      call check(error, inter(1, 2), inter(2, 1), thr=TOL, more="inter is not symmetric")
   end subroutine inter_is_symmetric_with_no_diagonal

   subroutine rejects_an_orbital_off_the_molecule(error)
      !! An orbital on atom 2 with only one atom declared would be dropped
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp), allocatable :: intra(:), inter(:, :)

      call two_atom_case(quao)
      call kinetic_decomposition(quao, T_INTF, 1, intra, inter, err)
      call check(error, err%has_error(), &
                 "an orbital assigned outside the molecule was accepted")
   end subroutine rejects_an_orbital_off_the_molecule

   subroutine rejects_a_matrix_of_the_wrong_size(error)
      !! A matrix that did not come from these orbitals
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp), allocatable :: intra(:), inter(:, :)
      real(dp) :: small(2, 2)

      small = 0.0_dp
      call two_atom_case(quao)
      call kinetic_decomposition(quao, small, 2, intra, inter, err)
      call check(error, err%has_error(), "a mismatched matrix was accepted")
   end subroutine rejects_a_matrix_of_the_wrong_size

   subroutine flat_field(v)
      !! Nucleus 1 contributes one everywhere, nucleus 2 contributes two
      real(dp), allocatable, intent(out) :: v(:, :, :)

      allocate (v(N_ORB, N_ORB, 2))
      v(:, :, 1) = 1.0_dp
      v(:, :, 2) = 2.0_dp
   end subroutine flat_field

   subroutine nuclear_terms_land_in_the_right_bins(error)
      !! The three-way assignment, against numbers worked out by hand
      !!
      !! A nuclear attraction term carries three atomic labels, so there are
      !! three places it can go and two ways to get it wrong that still sum
      !! correctly: charging an atom's density in a foreign field to that atom
      !! rather than to the pair, and charging a three-centre interference to
      !! the nucleus rather than to the two atoms sharing the density.
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp), allocatable :: v(:, :, :), intra(:), inter(:, :)

      call two_atom_case(quao)
      call flat_field(v)
      call nuclear_decomposition(quao, FLAT, v, 2, intra, inter, err)
      call check(error,.not. err%has_error(), "the decomposition refused a valid case")
      if (allocated(error)) return

      call check(error, intra(1), 4.0_dp, thr=TOL, more="atom 1's own-nucleus term")
      if (allocated(error)) return
      call check(error, intra(2), 8.0_dp, thr=TOL, more="atom 2's own-nucleus term")
      if (allocated(error)) return
      call check(error, inter(1, 2), 36.0_dp, thr=TOL, more="the pair term")
      if (allocated(error)) return
      call check(error, inter(1, 2), inter(2, 1), thr=TOL, more="inter is not symmetric")
   end subroutine nuclear_terms_land_in_the_right_bins

   subroutine nuclear_pieces_add_back_up(error)
      !! Every term reaches exactly one bin
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp), allocatable :: v(:, :, :), intra(:), inter(:, :)

      call two_atom_case(quao)
      call flat_field(v)
      call nuclear_decomposition(quao, FLAT, v, 2, intra, inter, err)
      call check(error,.not. err%has_error(), "the decomposition refused a valid case")
      if (allocated(error)) return

      call check(error, kinetic_total(intra, inter), 3.0_dp*sum(FLAT), thr=TOL, &
                 more="the pieces do not add up to the nuclear attraction energy")
   end subroutine nuclear_pieces_add_back_up

   subroutine rejects_a_partial_nuclear_field(error)
      !! A field covering fewer atoms than the molecule has would be silently short
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp), allocatable :: v(:, :, :), intra(:), inter(:, :)

      call two_atom_case(quao)
      call flat_field(v)
      call nuclear_decomposition(quao, FLAT, v, 3, intra, inter, err)
      call check(error, err%has_error(), "a nuclear field missing an atom was accepted")
   end subroutine rejects_a_partial_nuclear_field

   subroutine attraction_is_on_the_right_nucleus(error)
      !! Which nucleus each per-nucleus integral actually belongs to
      !!
      !! `nuclear_attraction_per_atom` checks its own pieces against the
      !! ordinary nuclear attraction, but that sum is **invariant to permuting
      !! the nuclei** -- hand the array back in the wrong order and the check
      !! still passes. This is the part it cannot see, so it is checked here.
      !!
      !! Two identical nuclei make it readable directly: each basis function is
      !! attracted most strongly by the nucleus it sits on, and with equal
      !! charges there is no `Z` to divide out first. In a heteronuclear
      !! molecule the heavy nucleus dominates every diagonal element regardless
      !! of where the function sits, which is why this is not water.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: v_atom(:, :, :)

      call build_libcint_molecule(H2_Z, H2_SYM, H2, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "could not build H2")
      if (allocated(error)) return

      call nuclear_attraction_per_atom(mol, H2_Z, H2, v_atom, err)
      call check(error,.not. err%has_error(), &
                 "the per-nucleus attraction did not sum to the ordinary one")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call check(error, v_atom(1, 1, 1) < v_atom(2, 2, 1), &
                 "nucleus 1 does not attract its own basis function most strongly")
      if (.not. allocated(error)) then
         call check(error, v_atom(2, 2, 2) < v_atom(1, 1, 2), &
                    "nucleus 2 does not attract its own basis function most strongly")
      end if
      if (.not. allocated(error)) then
         call check(error, v_atom(1, 1, 1), v_atom(2, 2, 2), thr=1.0e-10_dp, &
                    more="the two identical nuclei are not equivalent")
      end if
      call mol%destroy()
   end subroutine attraction_is_on_the_right_nucleus

end module test_mqc_ieda

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_ieda, only: collect_mqc_ieda_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_ieda", collect_mqc_ieda_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
