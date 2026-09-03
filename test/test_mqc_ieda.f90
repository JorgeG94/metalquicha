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
   use mqc_czt_atomic_guess, only: free_atom_energies
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_quao, only: quao_result_t
   use mqc_czt_ieda, only: kinetic_decomposition, kinetic_total, &
                           nuclear_attraction_per_atom, nuclear_decomposition, &
                           combine_quao_sets, quao_eris, two_electron_decomposition, &
                           nuclear_repulsion_pairs, active_cumulant, &
                           quao_projection, transform_cumulant, project_no_sharing, &
                           screened_nucleus_split
   implicit none
   private

   public :: collect_mqc_ieda_tests

   !! Two atoms, two orbitals each, symmetric as a real interference matrix must
   !! be. Chosen so every block sums to a round number:
   !!
   !!     intra(1) = 1.0 + 0.2 + 0.2 + 2.0            =  3.4
   !!     intra(2) = 3.0 + 0.6 + 0.6 + 4.0            =  8.2
   !!     one-way  = -0.5 - 0.1 - 0.3 - 0.4           = -1.3
   !!     inter    = 2 * one-way                      = -2.6
   !!     total    = 3.4 + 8.2 - 1.3 - 1.3            =  9.0
   integer, parameter :: N_ORB = 4
   integer, parameter :: ATOM_OF(N_ORB) = [1, 1, 2, 2]
   real(dp), parameter :: T_INTF(N_ORB, N_ORB) = reshape([ &
                                                         1.0_dp, 0.2_dp, -0.5_dp, -0.1_dp, &
                                                         0.2_dp, 2.0_dp, -0.3_dp, -0.4_dp, &
                                                         -0.5_dp, -0.3_dp, 3.0_dp, 0.6_dp, &
                                                         -0.1_dp, -0.4_dp, 0.6_dp, 4.0_dp], [N_ORB, N_ORB])
   real(dp), parameter :: TOL = 1.0e-12_dp

   !! Hydrogen, 1.4 Bohr apart. Two identical nuclei, so which nucleus a
   !! per-nucleus integral belongs to can be read straight off the numbers
   !! without dividing out a charge.
   integer, parameter :: H2_Z(2) = [1, 1]
   character(len=2), parameter :: H2_SYM(2) = ["H ", "H "]
   real(dp), parameter :: H2(3, 2) = reshape([ &
                                             0.0_dp, 0.0_dp, 0.0_dp, &
                                             0.0_dp, 0.0_dp, 1.4_dp], [3, 2])

   !! A flat density over the same two atoms, so every block sums to four and
   !! the grouping can be worked out in the margin:
   !!
   !!     intra(1) = 4 * 1                                =  4   (nucleus 1)
   !!     intra(2) = 4 * 2                                =  8   (nucleus 2)
   !!     inter    = 4*2 + 4*1  (own density, foreign nucleus)
   !!              + 2 * (4*1 + 4*2)  (interference, both orders, both nuclei)
   !!              = 12 + 24                              = 36
   !!     total    = 4 + 8 + 36                           = 48 = 3 * sum(D)
   real(dp), parameter :: FLAT(N_ORB, N_ORB) = 1.0_dp

   !! A second symmetric matrix, so a density and an operator can be rotated
   !! independently of each other rather than being the same array twice.
   real(dp), parameter :: DENS(N_ORB, N_ORB) = reshape([ &
                                                       2.0_dp, 0.3_dp, 0.7_dp, 0.2_dp, &
                                                       0.3_dp, 1.5_dp, 0.1_dp, 0.9_dp, &
                                                       0.7_dp, 0.1_dp, 1.8_dp, 0.4_dp, &
                                                       0.2_dp, 0.9_dp, 0.4_dp, 1.1_dp], [N_ORB, N_ORB])

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
                  new_unittest("attraction_is_on_the_right_nucleus", attraction_is_on_the_right_nucleus), &
                  new_unittest("survives_rotation_within_an_atom", survives_rotation_within_an_atom), &
                  new_unittest("core_joins_with_a_density_of_two", core_joins_with_a_density_of_two), &
                  new_unittest("two_electron_bins_by_hand", two_electron_bins_by_hand), &
                  new_unittest("two_electron_energy_has_a_closed_form", two_electron_energy_has_a_closed_form), &
                  new_unittest("eris_keep_their_index_order", eris_keep_their_index_order), &
                  new_unittest("nuclear_repulsion_is_already_pairwise", nuclear_repulsion_is_already_pairwise), &
                  new_unittest("cumulant_vanishes_for_a_determinant", cumulant_vanishes_for_a_determinant), &
                  new_unittest("cumulant_survives_an_isometry", cumulant_survives_an_isometry), &
                  new_unittest("projection_measures_what_it_misses", projection_measures_what_it_misses), &
                  new_unittest("cumulant_reaches_the_energy", cumulant_reaches_the_energy), &
                  new_unittest("free_atoms_are_equivalent_and_bound", free_atoms_are_equivalent_and_bound), &
                  new_unittest("classical_share_is_separated_out", classical_share_is_separated_out), &
                  new_unittest("no_sharing_keeps_the_neutral_determinants", no_sharing_keeps_the_neutral_determinants), &
                  new_unittest("no_sharing_is_a_projection_not_a_solve", no_sharing_is_a_projection_not_a_solve), &
                  new_unittest("screened_nucleus_only_relabels", screened_nucleus_only_relabels), &
                  new_unittest("screened_nucleus_refuses_a_soft_core", screened_nucleus_refuses_a_soft_core) &
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

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: v_atom(:, :, :)

      call build_czt_molecule(H2_Z, H2_SYM, H2, "sto-3g", mol, err)
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

   subroutine survives_rotation_within_an_atom(error)
      !! Atom and pair totals do not depend on the intra-atomic orientation
      !!
      !! Rotating among one atom's orbitals transforms the density and the
      !! operator by the same rotation, and inside `sum_pq gamma_pq O_pq`
      !! restricted to a block it cancels. Individual orbital pairs move a lot;
      !! the atom and atom-pair totals must not move at all.
      !!
      !! This is what lets the core orbitals go in unoriented, which is worth
      !! having as an assertion rather than an argument -- the core set really
      !! is never oriented, so if the invariance failed the printed numbers
      !! would depend on how LAPACK happened to return a subspace.
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp), allocatable :: intra(:), inter(:, :), spun_intra(:), spun_inter(:, :)
      real(dp) :: u(N_ORB, N_ORB), d(N_ORB, N_ORB), o(N_ORB, N_ORB)
      real(dp) :: c, sn
      integer :: i

      ! Block diagonal: each atom's pair of orbitals rotated among themselves.
      c = cos(0.7_dp)
      sn = sin(0.7_dp)
      u = 0.0_dp
      u(1, 1) = c; u(1, 2) = -sn; u(2, 1) = sn; u(2, 2) = c
      u(3, 3) = c; u(3, 4) = -sn; u(4, 3) = sn; u(4, 4) = c

      call two_atom_case(quao)
      call kinetic_decomposition(quao, DENS*T_INTF, 2, intra, inter, err)
      call check(error,.not. err%has_error(), "the decomposition refused a valid case")
      if (allocated(error)) return

      d = matmul(transpose(u), matmul(DENS, u))
      o = matmul(transpose(u), matmul(T_INTF, u))
      call kinetic_decomposition(quao, d*o, 2, spun_intra, spun_inter, err)
      call check(error,.not. err%has_error(), "the decomposition refused the rotated case")
      if (allocated(error)) return

      ! The rotation has to have actually done something, or this passes for
      ! the wrong reason.
      call check(error, maxval(abs(d*o - DENS*T_INTF)) > 0.1_dp, &
                 "the rotation left the orbital-pair values alone, so this proves nothing")
      if (allocated(error)) return

      do i = 1, 2
         call check(error, spun_intra(i), intra(i), thr=1.0e-12_dp, &
                    more="an atom total moved under an intra-atomic rotation")
         if (allocated(error)) return
      end do
      call check(error, spun_inter(1, 2), inter(1, 2), thr=1.0e-12_dp, &
                 more="a pair total moved under an intra-atomic rotation")
   end subroutine survives_rotation_within_an_atom

   subroutine core_joins_with_a_density_of_two(error)
      !! Stacking the core in front of the valence
      !!
      !! The core density is exactly two on the diagonal and zero everywhere
      !! else, including in the core-valence blocks, because the core orbitals
      !! span a space the density projects onto entirely and which is orthogonal
      !! to the valence one.
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: core, valence, combined
      type(error_t) :: err

      core%n_quao = 2
      allocate (core%atom_of(2), core%orbitals(6, 2), core%population_bond_order(2, 2))
      core%atom_of = [1, 2]
      core%orbitals = 0.0_dp
      core%population_bond_order = 0.0_dp

      valence%n_quao = N_ORB
      allocate (valence%atom_of(N_ORB), valence%orbitals(6, N_ORB))
      allocate (valence%population_bond_order(N_ORB, N_ORB))
      valence%atom_of = ATOM_OF
      valence%orbitals = 1.0_dp
      valence%population_bond_order = DENS

      call combine_quao_sets(core, valence, combined, err)
      call check(error,.not. err%has_error(), "the sets would not combine")
      if (allocated(error)) return

      call check(error, combined%n_quao, 6, "the combined set is the wrong size")
      if (allocated(error)) return
      call check(error, all(combined%atom_of == [1, 2, 1, 1, 2, 2]), &
                 "the atom assignments were not carried across")
      if (allocated(error)) return
      call check(error, combined%population_bond_order(1, 1), 2.0_dp, thr=TOL, &
                 more="a core orbital does not hold two electrons")
      if (allocated(error)) return
      call check(error, maxval(abs(combined%population_bond_order(1:2, 3:6))), 0.0_dp, &
                 thr=TOL, more="the core-valence block is not zero")
      if (allocated(error)) return
      call check(error, maxval(abs(combined%population_bond_order(3:6, 3:6) - DENS)), &
                 0.0_dp, thr=TOL, more="the valence density was not carried across")
   end subroutine core_joins_with_a_density_of_two

   subroutine two_orbital_case(quao)
      !! One orbital on each of two atoms, so an atom label is an orbital label
      type(quao_result_t), intent(out) :: quao

      quao%n_quao = 2
      allocate (quao%atom_of(2))
      quao%atom_of = [1, 2]
   end subroutine two_orbital_case

   subroutine two_electron_bins_by_hand(error)
      !! Where each of the four labels sends a term, worked out on paper
      !!
      !! Two atoms, one orbital each, a density of two on each diagonal and
      !! every integral equal to one. With `gamma` diagonal only three families
      !! of term survive:
      !!
      !!   * all four indices on one atom -- Coulomb 4 less half of exchange 4,
      !!     halved, so 1 to that atom's own bin;
      !!   * `(1,1,2,2)` and `(2,2,1,1)` -- pure Coulomb, 2 each, and since each
      !!     distribution sits on one atom and the atoms differ they are the
      !!     Coulombic interaction of the pair;
      !!   * `(1,2,2,1)` and `(2,1,1,2)` -- pure exchange, -1 each, and both
      !!     distributions are spread across the two atoms, so each is split in
      !!     half between `{1,2}` and `{2,1}`, which is the same pair.
      !!
      !! Giving `intra = 1, 1` and a pair energy of `2 + 2 - 1 - 1 = 2`, with
      !! the four bins summing to `0.25 * (sum gamma)^2 = 4`.
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp) :: d(2, 2), eri(2, 2, 2, 2), energy
      real(dp) :: ones(N_ORB, N_ORB, N_ORB, N_ORB)
      real(dp), allocatable :: intra(:), inter(:, :)

      d = 0.0_dp
      d(1, 1) = 2.0_dp
      d(2, 2) = 2.0_dp
      eri = 1.0_dp

      call two_orbital_case(quao)
      call two_electron_decomposition(quao, d, eri, 2, intra, inter, energy, err)
      call check(error,.not. err%has_error(), "the decomposition refused a valid case")
      if (allocated(error)) return

      call check(error, energy, 4.0_dp, thr=TOL, more="the two-electron energy")
      if (allocated(error)) return
      call check(error, intra(1), 1.0_dp, thr=TOL, more="atom 1's own bin")
      if (allocated(error)) return
      call check(error, intra(2), 1.0_dp, thr=TOL, more="atom 2's own bin")
      if (allocated(error)) return
      call check(error, inter(1, 2), 2.0_dp, thr=TOL, more="the pair bin")
   end subroutine two_electron_bins_by_hand

   subroutine two_electron_energy_has_a_closed_form(error)
      !! With every integral equal to one the energy collapses analytically
      !!
      !!     E2 = (1/2) sum [gamma_pq gamma_rs - (1/2) gamma_ps gamma_rq]
      !!        = (1/2)(sum gamma)^2 - (1/4)(sum gamma)^2
      !!        = (1/4)(sum gamma)^2
      !!
      !! for **any** density, because both double sums factorize the same way.
      !! So this pins the Coulomb and exchange prefactors against a number that
      !! owes nothing to the implementation, using an off-diagonal density where
      !! the hand-worked case above could not.
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp) :: eri(N_ORB, N_ORB, N_ORB, N_ORB), energy
      real(dp), allocatable :: intra(:), inter(:, :)

      eri = 1.0_dp
      call two_atom_case(quao)
      call two_electron_decomposition(quao, DENS, eri, 2, intra, inter, energy, err)
      call check(error,.not. err%has_error(), "the decomposition refused a valid case")
      if (allocated(error)) return

      call check(error, energy, 0.25_dp*sum(DENS)**2, thr=1.0e-10_dp, &
                 more="the Coulomb and exchange prefactors do not give the closed form")
      if (allocated(error)) return
      call check(error, kinetic_total(intra, inter), energy, thr=1.0e-10_dp, &
                 more="the bins do not sum to the energy")
   end subroutine two_electron_energy_has_a_closed_form

   subroutine eris_keep_their_index_order(error)
      !! That four quarter-transformations leave the indices where they started
      !!
      !! The transformation cycles each contracted index to the back, so after
      !! four passes the order is restored -- provided every pass agrees about
      !! which index it is on. A permutation matrix makes that readable: with
      !! `C` swapping the two basis functions, the transformed array must be the
      !! original with all four indices swapped, and any pair of passes that
      !! disagree shows up as an element in the wrong place.
      !!
      !! An identity matrix would not catch it. A permutation does.
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp) :: eri(2, 2, 2, 2)
      real(dp), allocatable :: transformed(:, :, :, :)
      integer :: p, q, r, sx, swap(2)

      swap = [2, 1]
      do p = 1, 2
         do q = 1, 2
            do r = 1, 2
               do sx = 1, 2
                  eri(p, q, r, sx) = 1000.0_dp*p + 100.0_dp*q + 10.0_dp*r + real(sx, dp)
               end do
            end do
         end do
      end do

      call two_orbital_case(quao)
      allocate (quao%orbitals(2, 2), quao%population_bond_order(2, 2))
      quao%orbitals = 0.0_dp
      quao%orbitals(2, 1) = 1.0_dp
      quao%orbitals(1, 2) = 1.0_dp
      quao%population_bond_order = 0.0_dp

      call quao_eris(quao, eri, transformed, err)
      call check(error,.not. err%has_error(), "the transformation refused a valid case")
      if (allocated(error)) return

      do p = 1, 2
         do q = 1, 2
            do r = 1, 2
               do sx = 1, 2
                  call check(error, transformed(p, q, r, sx), &
                             eri(swap(p), swap(q), swap(r), swap(sx)), thr=1.0e-10_dp, &
                             more="an index came out of the transformation in the wrong place")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine eris_keep_their_index_order

   subroutine nuclear_repulsion_is_already_pairwise(error)
      !! The one term that arrives decomposed
      !!
      !! Two unit charges 1.4 Bohr apart repel by `1/1.4` hartree, and that has
      !! to survive the pair convention intact: the matrix holds the whole pair
      !! energy in both elements, so it is `1/1.4` in each and `kinetic_total`
      !! halves the matrix back to `1/1.4` for the molecule. Written the other
      !! way round it would be half the nuclear repulsion, which is large enough
      !! to notice and small enough to argue about.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: inter(:, :)
      real(dp), allocatable :: intra(:)

      call nuclear_repulsion_pairs(H2_Z, H2, inter, err)
      call check(error,.not. err%has_error(), "the nuclear repulsion was refused")
      if (allocated(error)) return

      call check(error, inter(1, 2), 1.0_dp/1.4_dp, thr=1.0e-12_dp, &
                 more="the pair holds the wrong repulsion")
      if (allocated(error)) return
      call check(error, inter(2, 1), inter(1, 2), thr=1.0e-12_dp, more="not symmetric")
      if (allocated(error)) return
      call check(error, inter(1, 1), 0.0_dp, thr=1.0e-12_dp, more="an atom repels itself")
      if (allocated(error)) return

      allocate (intra(2))
      intra = 0.0_dp
      call check(error, kinetic_total(intra, inter), 1.0_dp/1.4_dp, thr=1.0e-12_dp, &
                 more="the molecule's nuclear repulsion came out wrong")
   end subroutine nuclear_repulsion_is_already_pairwise

   subroutine cumulant_vanishes_for_a_determinant(error)
      !! Correlation adds nothing when there is none
      !!
      !! The cumulant is `d` less the two-particle density a determinant would
      !! have, so feeding it a determinant's own `d` must give exactly zero.
      !! That is the sharpest available check on the expression, because it
      !! fails for any wrong sign, transposed index or factor: none of those
      !! cancel against the same expression written the other way round.
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: d1(N_ORB, N_ORB), d2(N_ORB, N_ORB, N_ORB, N_ORB)
      real(dp), allocatable :: cumulant(:, :, :, :)
      integer :: p, q, r, sx

      ! Two electrons in each of the first two orbitals, none in the rest.
      d1 = 0.0_dp
      d1(1, 1) = 2.0_dp
      d1(2, 2) = 2.0_dp
      do p = 1, N_ORB
         do q = 1, N_ORB
            do r = 1, N_ORB
               do sx = 1, N_ORB
                  d2(p, q, r, sx) = d1(p, q)*d1(r, sx) - 0.5_dp*d1(p, sx)*d1(r, q)
               end do
            end do
         end do
      end do

      call active_cumulant(d1, d2, cumulant)
      call check(error, maxval(abs(cumulant)), 0.0_dp, thr=1.0e-14_dp, &
                 more="a determinant was given a non-zero cumulant")
   end subroutine cumulant_vanishes_for_a_determinant

   subroutine cumulant_survives_an_isometry(error)
      !! Changing basis through an isometry loses nothing
      !!
      !! `U` maps a small active space into a larger quasi-atomic one, so it is
      !! rectangular and cannot be inverted -- but `U^T U = 1`, and that is
      !! enough for the transformation to be undone by its own transpose. If it
      !! were not, the correlation energy would depend on which basis it was
      !! contracted in, which is the failure this guards.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp) :: u(N_ORB, 2), small(2, 2, 2, 2)
      real(dp), allocatable :: big(:, :, :, :), back(:, :, :, :)
      integer :: p, q, r, sx

      ! An isometry from two dimensions into four.
      u = 0.0_dp
      u(1, 1) = 1.0_dp
      u(3, 2) = 1.0_dp

      do p = 1, 2
         do q = 1, 2
            do r = 1, 2
               do sx = 1, 2
                  small(p, q, r, sx) = 1.0_dp*p + 0.1_dp*q + 0.01_dp*r + 0.001_dp*sx
               end do
            end do
         end do
      end do

      call transform_cumulant(u, small, big, err)
      call check(error,.not. err%has_error(), "the transformation was refused")
      if (allocated(error)) return
      call check(error, size(big, 1), N_ORB, "the transformed array is the wrong size")
      if (allocated(error)) return

      call transform_cumulant(transpose(u), big, back, err)
      call check(error,.not. err%has_error(), "the reverse transformation was refused")
      if (allocated(error)) return
      call check(error, maxval(abs(back - small)), 0.0_dp, thr=1.0e-12_dp, &
                 more="the round trip through the isometry did not return the cumulant")
   end subroutine cumulant_survives_an_isometry

   subroutine projection_measures_what_it_misses(error)
      !! The deficit is reported, not assumed away
      !!
      !! Orbitals inside the quasi-atomic span project to an isometry and the
      !! deficit is zero. Orbitals partly outside it do not, and the number that
      !! comes back is how much of them the analysis cannot see. Real active
      !! orbitals are the second case -- they were optimised to lower an energy,
      !! the valence-virtual ones to look atomic -- so a routine that treated
      !! this as an error would refuse every correlated wave function.
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp) :: overlap(3, 3), inside(3, 1), outside(3, 1)
      real(dp), allocatable :: u(:, :)
      real(dp) :: deficit

      ! Three atomic orbitals, orthonormal, of which the quasi-atomic set spans
      ! the first two.
      overlap = 0.0_dp
      overlap(1, 1) = 1.0_dp
      overlap(2, 2) = 1.0_dp
      overlap(3, 3) = 1.0_dp

      quao%n_quao = 2
      allocate (quao%atom_of(2), quao%orbitals(3, 2), quao%population_bond_order(2, 2))
      quao%atom_of = [1, 2]
      quao%orbitals = 0.0_dp
      quao%orbitals(1, 1) = 1.0_dp
      quao%orbitals(2, 2) = 1.0_dp
      quao%population_bond_order = 0.0_dp

      inside = 0.0_dp
      inside(1, 1) = 1.0_dp
      call quao_projection(quao, overlap, inside, u, deficit, err)
      call check(error,.not. err%has_error(), "the projection was refused")
      if (allocated(error)) return
      call check(error, deficit, 0.0_dp, thr=1.0e-14_dp, &
                 more="an orbital inside the span was reported as outside it")
      if (allocated(error)) return

      ! Half in the span, half out.
      outside = 0.0_dp
      outside(1, 1) = sqrt(0.5_dp)
      outside(3, 1) = sqrt(0.5_dp)
      call quao_projection(quao, overlap, outside, u, deficit, err)
      call check(error,.not. err%has_error(), "the projection was refused")
      if (allocated(error)) return
      call check(error, deficit, 0.5_dp, thr=1.0e-12_dp, &
                 more="the missing half of the orbital was not measured")
      if (allocated(error)) return

      ! **With an orthogonal metric the overlap could be dropped entirely and
      ! nothing above would notice**, since it is the identity. So the metric is
      ! made non-trivial here: atomic orbital 3 now overlaps orbital 1 by a half,
      ! which is the only reason a function along 3 projects onto a quasi-atomic
      ! set spanning only 1. Without `S` the projection is zero and the deficit
      ! comes back as one instead.
      overlap(1, 3) = 0.5_dp
      overlap(3, 1) = 0.5_dp
      deallocate (quao%atom_of, quao%orbitals, quao%population_bond_order)
      quao%n_quao = 1
      allocate (quao%atom_of(1), quao%orbitals(3, 1), quao%population_bond_order(1, 1))
      quao%atom_of = [1]
      quao%orbitals = 0.0_dp
      quao%orbitals(1, 1) = 1.0_dp
      quao%population_bond_order = 0.0_dp

      outside = 0.0_dp
      outside(3, 1) = 1.0_dp
      call quao_projection(quao, overlap, outside, u, deficit, err)
      call check(error,.not. err%has_error(), "the projection was refused")
      if (allocated(error)) return
      call check(error, u(1, 1), 0.5_dp, thr=1.0e-12_dp, &
                 more="the projection ignored the overlap metric")
      if (allocated(error)) return
      call check(error, deficit, 0.75_dp, thr=1.0e-12_dp, &
                 more="the deficit ignored the overlap metric")
   end subroutine projection_measures_what_it_misses

   subroutine cumulant_reaches_the_energy(error)
      !! That a cumulant handed in actually changes the answer
      !!
      !! With every integral one, the determinant part collapses to
      !! `(1/4)(sum gamma)^2` and a constant cumulant `c` adds `(1/2) c n^4` on
      !! top, both in closed form. Without this the optional argument could be
      !! ignored entirely -- every other test passes a wave function with no
      !! correlation in it, so a dropped cumulant would look exactly like a
      !! correct one.
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp) :: eri(N_ORB, N_ORB, N_ORB, N_ORB)
      real(dp) :: lambda(N_ORB, N_ORB, N_ORB, N_ORB)
      real(dp) :: energy, plain, expected
      real(dp), allocatable :: intra(:), inter(:, :)

      eri = 1.0_dp
      lambda = 0.1_dp
      call two_atom_case(quao)

      call two_electron_decomposition(quao, DENS, eri, 2, intra, inter, plain, err)
      call check(error,.not. err%has_error(), "the uncorrelated case was refused")
      if (allocated(error)) return
      deallocate (intra, inter)

      call two_electron_decomposition(quao, DENS, eri, 2, intra, inter, energy, err, &
                                      cumulant=lambda)
      call check(error,.not. err%has_error(), "the correlated case was refused")
      if (allocated(error)) return

      expected = 0.25_dp*sum(DENS)**2 + 0.5_dp*0.1_dp*real(N_ORB, dp)**4
      call check(error, energy, expected, thr=1.0e-10_dp, &
                 more="the cumulant did not reach the energy")
      if (allocated(error)) return
      call check(error, abs(energy - plain) > 1.0_dp, &
                 "the cumulant changed nothing, so this proves nothing")
      if (allocated(error)) return
      call check(error, kinetic_total(intra, inter), energy, thr=1.0e-10_dp, &
                 more="the correlated bins do not sum to the correlated energy")
   end subroutine cumulant_reaches_the_energy

   subroutine free_atoms_are_equivalent_and_bound(error)
      !! The reference an energy of formation is measured against
      !!
      !! Two hydrogens in the same basis must give the same free-atom energy --
      !! the reference depends on the element and its basis functions, not on
      !! where the atom sits -- and each must be bound. Getting the two mixed up
      !! is the failure that a sum rule cannot see, since swapping references
      !! between atoms of the same element changes nothing and swapping them
      !! between different elements still leaves the total intact.
      !!
      !! Hydrogen in STO-3G is a single function, so its unrestricted doublet is
      !! that function singly occupied and the energy is above the exact -0.5.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: energies(:)

      call build_czt_molecule(H2_Z, H2_SYM, H2, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "could not build H2")
      if (allocated(error)) return

      call free_atom_energies(mol, energies, err)
      call check(error,.not. err%has_error(), "the free atoms would not converge")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call check(error, size(energies), 2, "one energy per atom was expected")
      if (.not. allocated(error)) then
         call check(error, energies(1), energies(2), thr=1.0e-10_dp, &
                    more="two atoms of the same element got different references")
      end if
      if (.not. allocated(error)) then
         call check(error, energies(1) < 0.0_dp .and. energies(1) > -0.5_dp, &
                    "a hydrogen atom in a minimal basis should sit just above -0.5")
      end if
      call mol%destroy()
   end subroutine free_atoms_are_equivalent_and_bound

   subroutine classical_share_is_separated_out(error)
      !! Which part of a pair interaction an electrostatic model could produce
      !!
      !! Both decompositions carry the same numbers as before; what is new is
      !! that the classical part is reported separately. So the thing to pin is
      !! that the split is a partition -- the classical share plus what is left
      !! must be the pair energy that was already being reported, exactly -- and
      !! that each term lands on the correct side of it.
      !!
      !! For the nuclear fixture the classical share is one atom's density in
      !! the other's field: atom 1's four units in nucleus 2's field at strength
      !! two, and atom 2's four in nucleus 1's at strength one, giving twelve of
      !! the thirty-six. For the two-electron fixture it is the pure Coulomb
      !! terms `(1,1,2,2)` and `(2,2,1,1)`, two each, giving four of the two --
      !! the interference part is negative there, which is the exchange.
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp), allocatable :: v(:, :, :), intra(:), inter(:, :), coulomb(:, :)
      real(dp) :: d(2, 2), eri(2, 2, 2, 2), energy, expected
      real(dp) :: ones(N_ORB, N_ORB, N_ORB, N_ORB)

      ones = 1.0_dp
      call two_atom_case(quao)
      call flat_field(v)
      call nuclear_decomposition(quao, FLAT, v, 2, intra, inter, err, coulomb=coulomb)
      call check(error,.not. err%has_error(), "the decomposition refused a valid case")
      if (allocated(error)) return

      call check(error, coulomb(1, 2), 12.0_dp, thr=TOL, &
                 more="the nuclear classical share is wrong")
      if (allocated(error)) return
      call check(error, coulomb(1, 2), coulomb(2, 1), thr=TOL, &
                 more="the classical share is not symmetric")
      if (allocated(error)) return
      call check(error, abs(coulomb(1, 2)) < abs(inter(1, 2)), &
                 "the classical share is not a part of the pair energy")
      if (allocated(error)) return
      deallocate (intra, inter, coulomb)

      d = 0.0_dp
      d(1, 1) = 2.0_dp
      d(2, 2) = 2.0_dp
      eri = 1.0_dp
      call two_orbital_case(quao)
      call two_electron_decomposition(quao, d, eri, 2, intra, inter, energy, err, &
                                      coulomb=coulomb)
      call check(error,.not. err%has_error(), "the decomposition refused a valid case")
      if (allocated(error)) return
      call check(error, coulomb(1, 2), 4.0_dp, thr=TOL, &
                 more="the two-electron classical share is wrong")
      if (allocated(error)) return
      ! The remainder is exchange and is negative, so a split that quietly
      ! clamped it to zero would still look plausible.
      call check(error, inter(1, 2) - coulomb(1, 2), -2.0_dp, thr=TOL, &
                 more="the interference remainder is wrong")
      if (allocated(error)) return
      deallocate (intra, inter, coulomb)

      ! **A diagonal density never reaches the branch where the bra is spread
      ! and the ket is not**, so the fixture above cannot tell whether that
      ! branch is charged classically or as interference -- and getting it wrong
      ! there passed every check until this was added. An off-diagonal density
      ! reaches all four branches.
      !
      ! The reference is built here from the definition rather than from the
      ! routine: a term is classical when each of its two charge distributions
      ! sits wholly on one atom and the two atoms differ. Summing exactly those
      ! terms is an independent statement of the selection rule, using the same
      ! weight the closed-form test already pinned.
      call two_atom_case(quao)
      call classical_by_definition(DENS, expected)
      call two_electron_decomposition(quao, DENS, ones, 2, intra, inter, energy, err, &
                                      coulomb=coulomb)
      call check(error,.not. err%has_error(), "the decomposition refused a valid case")
      if (allocated(error)) return
      call check(error, coulomb(1, 2), expected, thr=1.0e-10_dp, &
                 more="a term was put on the wrong side of the classical split")
      if (allocated(error)) return
      call check(error, abs(inter(1, 2) - coulomb(1, 2)) > 1.0e-6_dp, &
                 "this density produced no interference, so it proves nothing")
   end subroutine classical_share_is_separated_out

   subroutine classical_by_definition(d, expected)
      !! The classical share of the pair, summed straight from the selection rule
      !!
      !! Both distributions wholly on one atom, on different atoms. Written from
      !! the definition rather than lifted from the routine under test, so that
      !! the two agreeing means something.
      real(dp), intent(in) :: d(N_ORB, N_ORB)
      real(dp), intent(out) :: expected

      integer :: p, q, r, sx

      expected = 0.0_dp
      do p = 1, N_ORB
         do q = 1, N_ORB
            if (ATOM_OF(p) /= ATOM_OF(q)) cycle
            do r = 1, N_ORB
               do sx = 1, N_ORB
                  if (ATOM_OF(r) /= ATOM_OF(sx)) cycle
                  if (ATOM_OF(p) == ATOM_OF(r)) cycle
                  expected = expected + 0.5_dp*(d(p, q)*d(r, sx) &
                                                - 0.5_dp*d(p, sx)*d(r, q))
               end do
            end do
         end do
      end do
   end subroutine classical_by_definition

   subroutine no_sharing_keeps_the_neutral_determinants(error)
      !! Which determinants survive, counted against an independent enumeration
      !!
      !! Water's valence shell is six orbitals -- oxygen's 2s and three 2p, and
      !! one on each hydrogen -- holding eight electrons, so the full valence
      !! space is `C(6,4)^2 = 225` determinants. Requiring six electrons on the
      !! oxygen and one on each hydrogen leaves 44 of them, which was counted by
      !! hand from the definition before this was written and not read off the
      !! routine.
      !!
      !! The count is the sharp test. A projector that kept too much would still
      !! renormalise to one and still produce a plausible density; only the
      !! number of survivors says whether the condition being applied is the
      !! right one.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: ci(:, :)
      real(dp) :: recovered
      integer :: n_kept
      integer, parameter :: WATER_ATOM_OF(6) = [1, 1, 1, 1, 2, 3]
      integer, parameter :: WATER_NEUTRAL(3) = [6, 1, 1]

      allocate (ci(15, 15))
      ci = 1.0_dp/15.0_dp     ! normalised, and flat so every determinant counts
      call project_no_sharing(WATER_ATOM_OF, 3, WATER_NEUTRAL, 4, 4, ci, &
                              recovered, n_kept, err)
      call check(error,.not. err%has_error(), "the projection was refused")
      if (allocated(error)) return

      call check(error, n_kept, 44, "the wrong number of determinants survived")
      if (allocated(error)) return
      ! A flat vector puts equal weight everywhere, so the surviving fraction is
      ! the determinant fraction exactly.
      call check(error, recovered, 44.0_dp/225.0_dp, thr=1.0e-12_dp, &
                 more="the recovered norm does not match the surviving fraction")
      if (allocated(error)) return
      call check(error, sum(ci**2), 1.0_dp, thr=1.0e-12_dp, &
                 more="the projection was not renormalised")
   end subroutine no_sharing_keeps_the_neutral_determinants

   subroutine no_sharing_is_a_projection_not_a_solve(error)
      !! The surviving amplitudes are the parent's, only rescaled
      !!
      !! This is the distinction the papers are emphatic about and the one an
      !! implementation is most likely to lose: optimising a wave function
      !! inside the neutral space gives a lower energy and a different state.
      !! So every kept coefficient must stay in the same ratio to every other
      !! kept coefficient as it was before, which a re-solve would not respect.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: ci(:, :), before(:, :)
      real(dp) :: recovered, scale
      integer :: n_kept, ia, ib, first_a, first_b
      integer, parameter :: WATER_ATOM_OF(6) = [1, 1, 1, 1, 2, 3]
      integer, parameter :: WATER_NEUTRAL(3) = [6, 1, 1]

      allocate (ci(15, 15), before(15, 15))
      do ib = 1, 15
         do ia = 1, 15
            ci(ia, ib) = sin(real(3*ia + 7*ib, dp))
         end do
      end do
      ci = ci/sqrt(sum(ci**2))
      before = ci

      call project_no_sharing(WATER_ATOM_OF, 3, WATER_NEUTRAL, 4, 4, ci, &
                              recovered, n_kept, err)
      call check(error,.not. err%has_error(), "the projection was refused")
      if (allocated(error)) return

      ! Find a survivor to fix the scale, then check every other survivor sits
      ! at the same multiple of what it was.
      first_a = 0
      first_b = 0
      do ib = 1, 15
         do ia = 1, 15
            if (abs(ci(ia, ib)) > 1.0e-12_dp) then
               first_a = ia
               first_b = ib
               exit
            end if
         end do
         if (first_a > 0) exit
      end do
      call check(error, first_a > 0, "every coefficient was struck out")
      if (allocated(error)) return

      scale = ci(first_a, first_b)/before(first_a, first_b)
      do ib = 1, 15
         do ia = 1, 15
            if (abs(ci(ia, ib)) < 1.0e-12_dp) cycle
            call check(error, ci(ia, ib), scale*before(ia, ib), thr=1.0e-12_dp, &
                       more="a surviving amplitude is not the parent's, rescaled")
            if (allocated(error)) return
         end do
      end do
      call check(error, abs(scale - 1.0_dp) > 1.0e-6_dp, &
                 "nothing was struck out, so this proves nothing")
   end subroutine no_sharing_is_a_projection_not_a_solve

   subroutine oxygen_and_hydrogen(quao, density)
      !! One core and one valence orbital on an oxygen, one valence on a hydrogen
      !!
      !! Core first across the whole set, which is what `combine_quao_sets`
      !! produces and what lets a count identify them.
      type(quao_result_t), intent(out) :: quao
      real(dp), intent(out) :: density(3, 3)

      quao%n_quao = 3
      allocate (quao%atom_of(3), quao%orbitals(3, 3), quao%population_bond_order(3, 3))
      quao%atom_of = [1, 1, 2]
      quao%orbitals = 0.0_dp
      quao%population_bond_order = 0.0_dp
      density = 0.0_dp
      density(1, 1) = 2.0_dp     ! the oxygen core
      density(2, 2) = 1.5_dp     ! the oxygen valence
      density(3, 3) = 0.8_dp     ! the hydrogen
   end subroutine oxygen_and_hydrogen

   subroutine screened_nucleus_only_relabels(error)
      !! The two shares add back to the term they divide
      !!
      !! Oxygen has `Z = 8` and a two-electron core, so its valence density's
      !! attraction to its own nucleus is divided one part in four to the core
      !! and three to the valence, while the core density keeps the whole charge.
      !! Hydrogen has no core, so nothing is taken from it -- which is the case
      !! that catches a division that forgot the element and used a fixed
      !! fraction.
      !!
      !! With every integral one the arithmetic is visible: the oxygen core term
      !! is 2, its valence term 1.5, so the columns are `2 + 0.25*1.5 = 2.375`
      !! and `0.75*1.5 = 1.125`, summing to 3.5.
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp) :: density(3, 3), v(3, 3, 2)
      real(dp), allocatable :: core_share(:), valence_share(:)

      call oxygen_and_hydrogen(quao, density)
      v = 1.0_dp
      call screened_nucleus_split(quao, density, v, 1, [8, 1], 2, &
                                  core_share, valence_share, err)
      call check(error,.not. err%has_error(), "the split was refused")
      if (allocated(error)) return

      call check(error, core_share(1), 2.375_dp, thr=1.0e-12_dp, &
                 more="the oxygen core share is wrong")
      if (allocated(error)) return
      call check(error, valence_share(1), 1.125_dp, thr=1.0e-12_dp, &
                 more="the oxygen valence share is wrong")
      if (allocated(error)) return
      call check(error, core_share(1) + valence_share(1), 3.5_dp, thr=1.0e-12_dp, &
                 more="the split did not conserve the oxygen term")
      if (allocated(error)) return

      call check(error, core_share(2), 0.0_dp, thr=1.0e-12_dp, &
                 more="hydrogen was given a core to screen with")
      if (allocated(error)) return
      call check(error, valence_share(2), 0.8_dp, thr=1.0e-12_dp, &
                 more="the hydrogen term did not survive intact")
   end subroutine screened_nucleus_only_relabels

   subroutine screened_nucleus_refuses_a_soft_core(error)
      !! A density connecting core and valence is refused rather than absorbed
      !!
      !! The split assumes the core is frozen, which makes the core-valence
      !! block of the density exactly zero. A correlated core breaks that, and
      !! the terms would otherwise be quietly swept into the core column where
      !! nothing would look wrong.
      type(error_type), allocatable, intent(out) :: error

      type(quao_result_t) :: quao
      type(error_t) :: err
      real(dp) :: density(3, 3), v(3, 3, 2)
      real(dp), allocatable :: core_share(:), valence_share(:)

      call oxygen_and_hydrogen(quao, density)
      density(1, 2) = 0.05_dp
      density(2, 1) = 0.05_dp
      v = 1.0_dp
      call screened_nucleus_split(quao, density, v, 1, [8, 1], 2, &
                                  core_share, valence_share, err)
      call check(error, err%has_error(), "a core-valence density was accepted")
   end subroutine screened_nucleus_refuses_a_soft_core

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
