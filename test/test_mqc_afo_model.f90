!! Unit tests for the model system a frozen orbital is derived from
module test_mqc_afo_model
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: system_geometry_t, to_bohr, to_angstrom
   use mqc_bond_perception, only: find_severed_bonds, severed_bond_t
   use mqc_czt_afo, only: afo_model_t, build_afo_model, cuts_outside_group, &
                          group_electron_shift, peptide_bond_advice
   implicit none

   !! Geometries are quoted to four decimals in Angstrom, so a distance
   !! reconstructed from them agrees to about there and not further.
   real(dp), parameter :: TOL = 1.0e-4_dp

   private
   public :: collect_mqc_afo_model

contains

   subroutine collect_mqc_afo_model(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("model_of_a_whole_molecule_needs_no_caps", test_ethane_whole), &
                  new_unittest("model_caps_every_bond_it_leaves", test_propane_capped), &
                  new_unittest("model_cap_sits_at_a_standard_bond_length", test_cap_length), &
                  new_unittest("model_keeps_the_hydrogens_of_what_it_took", test_hydrogens), &
                  new_unittest("model_keeps_a_carbonyl_oxygen_rather_than_capping_it", &
                               test_terminal_oxygen), &
                  new_unittest("a_backbone_cut_gives_a_closed_shell_model", test_backbone_model), &
                  new_unittest("a_charged_terminus_is_counted_into_the_model", &
                               test_charged_termini), &
                  new_unittest("a_peptide_bond_names_the_c_alpha_cut_instead", test_peptide_advice), &
                  new_unittest("model_moves_rigidly_with_the_system", test_rotation), &
                  new_unittest("a_monomer_sees_only_its_own_boundaries", test_group_monomer), &
                  new_unittest("a_dimer_restores_the_bond_between_its_members", test_group_dimer), &
                  new_unittest("the_whole_system_has_no_boundaries_left", test_group_whole), &
                  new_unittest("electron_shifts_cancel_over_the_fragments", test_group_electrons) &
                  ]
   end subroutine collect_mqc_afo_model

   subroutine test_ethane_whole(error)
      !! A radius covering the molecule leaves nothing to cap
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(system_geometry_t) :: sys
      type(severed_bond_t), allocatable :: cuts(:)
      type(afo_model_t) :: model
      integer :: n_cuts

      call ethane(sys)
      call find_severed_bonds(sys, [1, 1, 1, 1, 2, 2, 2, 2], cuts, n_cuts)
      call check(error, n_cuts, 1, "ethane split at the C-C should cut one bond")
      if (allocated(error)) return

      call build_afo_model(sys%element_numbers, sys%coordinates, cuts(1), model, err, &
                           radius=5.0_dp)
      call check(error,.not. err%has_error(), "building the model failed")
      if (allocated(error)) return
      call check(error, model%n_caps, 0, "a model holding the whole molecule was capped")
      if (allocated(error)) return
      call check(error, model%n_atoms, 8, "the model is not the whole of ethane")
      if (allocated(error)) return
      call check(error, model%z(model%bda_local), 6, "the bond's first atom is not the carbon")
      if (allocated(error)) return
      call check(error, model%z(model%baa_local), 6, "the bond's second atom is not the carbon")
      if (allocated(error)) return
      call check(error, mod(model%nelec, 2), 0, "the model has an odd electron count")
   end subroutine test_ethane_whole

   subroutine test_propane_capped(error)
      !! A radius that stops short of the third carbon caps it away
      !!
      !! 1.4 Angstrom reaches the hydrogens at 1.09 and not the carbon at 1.53,
      !! so the model is the cut bond, its five hydrogens, and one cap standing
      !! in for the methyl that was left out -- ethane, in effect, which is what
      !! a model system for a C-C bond should be.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(system_geometry_t) :: sys
      type(severed_bond_t), allocatable :: cuts(:)
      type(afo_model_t) :: model
      integer :: n_cuts

      call propane(sys)
      call find_severed_bonds(sys, [1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2], cuts, n_cuts)
      call check(error, n_cuts, 1, "the partition should cut only the C1-C2 bond")
      if (allocated(error)) return

      call build_afo_model(sys%element_numbers, sys%coordinates, cuts(1), model, err, &
                           radius=1.4_dp)
      call check(error,.not. err%has_error(), "building the model failed")
      if (allocated(error)) return
      call check(error, model%n_caps, 1, "the methyl left outside the model was not capped")
      if (allocated(error)) return
      call check(error, model%n_atoms, 8, "the model is not the seven kept atoms plus a cap")
      if (allocated(error)) return
      call check(error, model%z(model%n_atoms), 1, "the cap is not a hydrogen")
      if (allocated(error)) return
      call check(error, mod(model%nelec, 2), 0, &
                 "the capping did not close every valence it opened")
   end subroutine test_propane_capped

   subroutine test_cap_length(error)
      !! The cap hydrogen sits at the standard C-H length, along the real bond
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(system_geometry_t) :: sys
      type(severed_bond_t), allocatable :: cuts(:)
      type(afo_model_t) :: model
      integer :: n_cuts
      real(dp) :: from_c(3), along(3), length, projected

      call propane(sys)
      call find_severed_bonds(sys, [1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2], cuts, n_cuts)
      call build_afo_model(sys%element_numbers, sys%coordinates, cuts(1), model, err, &
                           radius=1.4_dp)

      ! The cap replaces C3, which hung off C2 -- local slot 2.
      from_c = model%xyz(:, model%n_atoms) - model%xyz(:, 2)
      length = to_angstrom(sqrt(sum(from_c**2)))
      call check(error, abs(length - 1.07_dp) < 1.0e-3_dp, &
                 "the cap is not at the summed covalent radii of carbon and hydrogen")
      if (allocated(error)) return

      ! And it lies on the C2->C3 line rather than merely at the right distance.
      along = sys%coordinates(:, 3) - sys%coordinates(:, 2)
      along = along/sqrt(sum(along**2))
      projected = dot_product(from_c, along)
      call check(error, abs(projected - sqrt(sum(from_c**2))) < TOL, &
                 "the cap is off the axis of the bond it stands in for")
   end subroutine test_cap_length

   subroutine test_hydrogens(error)
      !! Hydrogens come in with their heavy atom, however the radius falls
      !!
      !! Taken wholesale rather than by distance on purpose: a hydrogen just
      !! outside the radius would otherwise be swapped for a cap hydrogen a
      !! few hundredths of an Angstrom away, which is a different molecule for
      !! no reason and a discontinuity in anything that moves the geometry.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(system_geometry_t) :: sys
      type(severed_bond_t), allocatable :: cuts(:)
      type(afo_model_t) :: model
      integer :: n_cuts

      call propane(sys)
      call find_severed_bonds(sys, [1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2], cuts, n_cuts)

      ! Well inside every C-H bond, so only the wholesale rule can bring them in.
      call build_afo_model(sys%element_numbers, sys%coordinates, cuts(1), model, err, &
                           radius=0.5_dp)
      call check(error,.not. err%has_error(), "building the model failed")
      if (allocated(error)) return
      call check(error, count(model%z(:model%n_atoms - model%n_caps) == 1), 5, &
                 "the hydrogens on the two carbons were not taken with them")
   end subroutine test_hydrogens

   subroutine test_rotation(error)
      !! Rotate the system and the model is the same molecule, moved
      !!
      !! Nothing about which atoms are chosen or where the cap goes may depend
      !! on the frame, so the model's internal distances have to come back
      !! unchanged. This is the cheap half of the covariance the hybrid orbital
      !! itself will have to satisfy.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(system_geometry_t) :: sys, spun
      type(severed_bond_t), allocatable :: cuts(:), spun_cuts(:)
      type(afo_model_t) :: model, spun_model
      integer :: n_cuts, n_spun, i, j
      real(dp) :: rot(3, 3), c, s, worst, d1, d2

      call propane(sys)
      call propane(spun)

      ! A rotation about no axis in particular, so no coordinate is spared.
      c = cos(0.7_dp)
      s = sin(0.7_dp)
      rot = reshape([c, -s, 0.0_dp, s, c, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
      spun%coordinates = matmul(rot, spun%coordinates)
      c = cos(0.4_dp)
      s = sin(0.4_dp)
      rot = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, c, -s, 0.0_dp, s, c], [3, 3])
      spun%coordinates = matmul(rot, spun%coordinates)

      call find_severed_bonds(sys, [1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2], cuts, n_cuts)
      call find_severed_bonds(spun, [1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2], spun_cuts, n_spun)
      call build_afo_model(sys%element_numbers, sys%coordinates, cuts(1), model, err, &
                           radius=1.4_dp)
      call build_afo_model(spun%element_numbers, spun%coordinates, spun_cuts(1), &
                           spun_model, err, radius=1.4_dp)

      call check(error, spun_model%n_atoms, model%n_atoms, &
                 "rotating the system changed which atoms the model holds")
      if (allocated(error)) return
      call check(error, spun_model%n_caps, model%n_caps, &
                 "rotating the system changed how many caps the model needs")
      if (allocated(error)) return

      worst = 0.0_dp
      do i = 1, model%n_atoms
         do j = i + 1, model%n_atoms
            d1 = sqrt(sum((model%xyz(:, i) - model%xyz(:, j))**2))
            d2 = sqrt(sum((spun_model%xyz(:, i) - spun_model%xyz(:, j))**2))
            worst = max(worst, abs(d1 - d2))
         end do
      end do
      call check(error, worst < 1.0e-10_dp, &
                 "the model is not the same molecule after the system was rotated")
   end subroutine test_rotation

   subroutine chain(cuts)
      !! Three fragments in a row, cut twice: 1-2 and 2-3
      !!
      !! Written down rather than perceived from a geometry, because what is
      !! under test is the bookkeeping over a partition and not the chemistry
      !! that produced it.
      type(severed_bond_t), allocatable, intent(out) :: cuts(:)

      allocate (cuts(2))
      cuts(1)%atom_a = 1
      cuts(1)%atom_b = 2
      cuts(1)%frag_a = 1
      cuts(1)%frag_b = 2
      cuts(2)%atom_a = 3
      cuts(2)%atom_b = 4
      cuts(2)%frag_a = 2
      cuts(2)%frag_b = 3
   end subroutine chain

   subroutine test_group_monomer(error)
      !! The middle fragment is cut on both sides, the ends on one
      type(error_type), allocatable, intent(out) :: error
      type(severed_bond_t), allocatable :: cuts(:)
      integer, allocatable :: outside(:)
      integer :: n

      call chain(cuts)

      call cuts_outside_group(cuts, 2, [1], outside, n)
      call check(error, n, 1, "the first fragment should be cut on one side only")
      if (allocated(error)) return

      call cuts_outside_group(cuts, 2, [2], outside, n)
      call check(error, n, 2, "the middle fragment should be cut on both sides")
      if (allocated(error)) return

      call cuts_outside_group(cuts, 2, [3], outside, n)
      call check(error, n, 1, "the last fragment should be cut on one side only")
   end subroutine test_group_monomer

   subroutine test_group_dimer(error)
      !! The rule the hydrogen caps got wrong
      !!
      !! A dimer holding both ends of a bond has restored it and must carry
      !! nothing standing in for it, while still being cut against whatever is
      !! outside. Decide this per fragment and inherit it, and the dimer arrives
      !! holding a stand-in for a bond it contains.
      type(error_type), allocatable, intent(out) :: error
      type(severed_bond_t), allocatable :: cuts(:)
      integer, allocatable :: outside(:)
      integer :: n

      call chain(cuts)

      call cuts_outside_group(cuts, 2, [1, 2], outside, n)
      call check(error, n, 1, &
                 "the dimer 1-2 should have restored its own bond and kept the other")
      if (allocated(error)) return
      call check(error, outside(1), 2, "the dimer kept the wrong one of the two cuts")
      if (allocated(error)) return

      ! Two fragments that are not neighbours restore nothing.
      call cuts_outside_group(cuts, 2, [1, 3], outside, n)
      call check(error, n, 2, &
                 "a dimer of two non-adjacent fragments restored a bond it does not hold")
   end subroutine test_group_dimer

   subroutine test_group_whole(error)
      !! Every fragment together is the molecule, with no boundary anywhere
      !!
      !! The property the two-fragment exactness test rests on: an n-mer holding
      !! everything has nothing frozen and no electrons moved, so it is an
      !! ordinary calculation on the whole system.
      type(error_type), allocatable, intent(out) :: error
      type(severed_bond_t), allocatable :: cuts(:)
      integer, allocatable :: outside(:)
      integer :: n

      call chain(cuts)
      call cuts_outside_group(cuts, 2, [1, 2, 3], outside, n)
      call check(error, n, 0, "the whole system still has a boundary")
      if (allocated(error)) return
      call check(error, group_electron_shift(cuts, 2, [1, 2, 3]), 0, &
                 "the whole system gained or lost electrons")
   end subroutine test_group_whole

   subroutine test_group_electrons(error)
      !! Electrons are moved between fragments, never created
      !!
      !! Of the pair, the detached end gets nothing of the bond and is an
      !! electron short; the attached end gets all of it and is one over. Summed
      !! over every fragment those cancel, and that is the check -- an
      !! assignment that did not cancel would be losing electrons somewhere.
      type(error_type), allocatable, intent(out) :: error
      type(severed_bond_t), allocatable :: cuts(:)
      integer :: total

      call chain(cuts)

      call check(error, group_electron_shift(cuts, 2, [1]), -1, &
                 "the detached end of the first bond did not give up its electron")
      if (allocated(error)) return
      call check(error, group_electron_shift(cuts, 2, [2]), 0, &
                 "the middle fragment takes one bond and gives up another, so nets zero")
      if (allocated(error)) return
      call check(error, group_electron_shift(cuts, 2, [3]), 1, &
                 "the attached end of the second bond did not take its electron")
      if (allocated(error)) return

      total = group_electron_shift(cuts, 2, [1]) + group_electron_shift(cuts, 2, [2]) &
              + group_electron_shift(cuts, 2, [3])
      call check(error, total, 0, "the electron shifts do not cancel over the fragments")
   end subroutine test_group_electrons

   subroutine ethane(sys)
      type(system_geometry_t), intent(out) :: sys
      real(dp) :: xyz(3, 8)

      xyz = reshape([0.000_dp, 0.000_dp, 0.768_dp, &
                     -1.019_dp, 0.000_dp, 1.157_dp, &
                     0.510_dp, 0.883_dp, 1.157_dp, &
                     0.510_dp, -0.883_dp, 1.157_dp, &
                     0.000_dp, 0.000_dp, -0.768_dp, &
                     1.019_dp, 0.000_dp, -1.157_dp, &
                     -0.510_dp, -0.883_dp, -1.157_dp, &
                     -0.510_dp, 0.883_dp, -1.157_dp], [3, 8])

      sys%total_atoms = 8
      sys%n_monomers = 0
      allocate (sys%element_numbers(8))
      sys%element_numbers = [6, 1, 1, 1, 6, 1, 1, 1]
      allocate (sys%coordinates(3, 8))
      sys%coordinates = to_bohr(xyz)
   end subroutine ethane

   subroutine test_terminal_oxygen(error)
      !! A doubly bonded oxygen comes in whole rather than as a cap hydrogen
      !!
      !! One cap hydrogen closes one electron pair. A carbonyl oxygen is held
      !! by two, so swapping it for a cap opens a valence nothing closes and
      !! the model comes back a radical -- which is what the parity check
      !! reports, one step too late to be useful. Bond perception here is
      !! distance-based and cannot tell a double bond from a short single one,
      !! so the order is not guessed: a singly bonded heavy neighbour is taken
      !! wholesale, exactly as a hydrogen is.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(system_geometry_t) :: sys
      type(severed_bond_t), allocatable :: cuts(:)
      type(afo_model_t) :: model
      integer :: n_cuts, oxygens

      call glycine_tripeptide(sys)
      ! The C-alpha--C(=O) bond of residue 2, atoms 10 and 11. Its sphere
      ! reaches residue 1's carbonyl carbon but not that carbon's oxygen.
      call find_severed_bonds(sys, owner_c_alpha(), cuts, n_cuts)
      call check(error, n_cuts, 2, "the partition should cut two C-alpha--C bonds")
      if (allocated(error)) return

      call build_afo_model(sys%element_numbers, sys%coordinates, cuts(2), model, err)
      call check(error,.not. err%has_error(), "building the model failed")
      if (allocated(error)) return

      oxygens = count(model%z(:model%n_atoms - model%n_caps) == 8)
      call check(error, oxygens, 2, &
                 "the model did not take both carbonyl oxygens whole")
   end subroutine test_terminal_oxygen

   subroutine test_backbone_model(error)
      !! Every C-alpha--C(=O) cut of a peptide gives a model that can be solved
      !!
      !! An odd electron count is refused rather than paired up, so this is
      !! the check that stands between a protein and a frozen orbital at all.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(system_geometry_t) :: sys
      type(severed_bond_t), allocatable :: cuts(:)
      type(afo_model_t) :: model
      integer :: n_cuts, c

      call glycine_tripeptide(sys)
      call find_severed_bonds(sys, owner_c_alpha(), cuts, n_cuts)
      do c = 1, n_cuts
         call build_afo_model(sys%element_numbers, sys%coordinates, cuts(c), model, err)
         call check(error,.not. err%has_error(), "a backbone model could not be built")
         if (allocated(error)) return
         call check(error, mod(model%nelec, 2), 0, &
                    "a backbone model came out with an odd electron count")
         if (allocated(error)) return
         call check(error, model%charge, 0, "a neutral backbone model was given a charge")
         if (allocated(error)) return
      end do
   end subroutine test_backbone_model

   subroutine test_charged_termini(error)
      !! The ammonium and carboxylate a zwitterionic peptide ends in
      !!
      !! The sphere around the first C-alpha--C(=O) bond takes the N-terminal
      !! nitrogen whole, hydrogens and all, and a sphere around the last
      !! C-alpha--C bond takes the C-terminal carboxylate whole. Charged, each
      !! has one proton more or less than a neutral model, so without its
      !! charge the model is a radical and the cut is refused -- which is what
      !! every charged protein did.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(system_geometry_t) :: neutral, sys
      type(severed_bond_t), allocatable :: cuts(:)
      type(severed_bond_t) :: last
      type(afo_model_t) :: model
      integer :: n_cuts, c, owner(24)

      call glycine_tripeptide(neutral)

      ! Zwitterion: a third hydrogen on N1, and none on the terminal O23.
      sys%total_atoms = 24
      sys%n_monomers = 0
      sys%element_numbers = [neutral%element_numbers(1:23), 1]
      allocate (sys%coordinates(3, 24))
      sys%coordinates(:, 1:23) = neutral%coordinates(:, 1:23)
      sys%coordinates(:, 24) = to_bohr([0.1428_dp, -1.4730_dp, 0.1766_dp])

      ! The new hydrogen is on the first fragment's nitrogen.
      owner = owner_c_alpha()
      owner(24) = 1
      call find_severed_bonds(sys, owner, cuts, n_cuts)
      call check(error, n_cuts, 2, "the zwitterion has two C-alpha--C(=O) cuts")
      if (allocated(error)) return
      do c = 1, n_cuts
         call build_afo_model(sys%element_numbers, sys%coordinates, cuts(c), model, err)
         call check(error,.not. err%has_error(), "a zwitterion backbone model was refused")
         if (allocated(error)) then
            write (*, *) "   message: ", trim(err%get_message())
            return
         end if
         call check(error, mod(model%nelec, 2), 0, "a zwitterion backbone model is "// &
                    "a radical")
         if (allocated(error)) return
         ! Only the first cut's sphere reaches the ammonium.
         if (min(cuts(c)%atom_a, cuts(c)%atom_b) == 2) then
            call check(error, model%charge, 1, "the N-terminal model is not +1")
         else
            call check(error, model%charge, 0, "a model away from both termini is charged")
         end if
         if (allocated(error)) return
      end do

      ! The last residue's C-alpha--C bond, which this partition leaves whole.
      last%atom_a = 17
      last%atom_b = 18
      last%frag_a = 3
      last%frag_b = 4
      call build_afo_model(sys%element_numbers, sys%coordinates, last, model, err)
      call check(error,.not. err%has_error(), "the C-terminal model was refused")
      if (allocated(error)) then
         write (*, *) "   message: ", trim(err%get_message())
         return
      end if
      call check(error, mod(model%nelec, 2), 0, "the C-terminal model is a radical")
      if (allocated(error)) return
      call check(error, model%charge, -1, "the C-terminal model is not -1")
   end subroutine test_charged_termini

   subroutine test_peptide_advice(error)
      !! A refusal on an amide says which bond to cut instead
      !!
      !! The one sentence a chemist needs: the backbone C(=O)-N carries more
      !! than one localized orbital and is refused, and the convention is to
      !! cut the C-alpha--C(=O) bond one place along. Atoms are numbered from
      !! one, as everywhere in the backend.
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      character(len=:), allocatable :: advice

      call glycine_tripeptide(sys)

      ! The peptide bond between residues 1 and 2.
      advice = peptide_bond_advice(sys%element_numbers, sys%coordinates, 3, 9)
      call check(error, len_trim(advice) > 0, "an amide C-N drew no advice")
      if (allocated(error)) return
      call check(error, index(advice, "C-alpha") > 0, &
                 "the advice does not name the C-alpha--C bond")
      if (allocated(error)) return
      call check(error, index(advice, "atoms 2 and 3") > 0, &
                 "the advice does not name the bond to cut instead")
      if (allocated(error)) return

      ! The C-alpha--C(=O) bond it points at is not itself an amide.
      advice = peptide_bond_advice(sys%element_numbers, sys%coordinates, 2, 3)
      call check(error, len_trim(advice), 0, "a C-C bond drew peptide advice")
      if (allocated(error)) return

      ! Neither is the amine nitrogen at the N terminus, which has no carbonyl.
      advice = peptide_bond_advice(sys%element_numbers, sys%coordinates, 1, 2)
      call check(error, len_trim(advice), 0, "an amine N-C drew peptide advice")
   end subroutine test_peptide_advice

   pure function owner_c_alpha() result(owner)
      !! Glycine tripeptide split at both C-alpha--C(=O) bonds
      !!
      !! The FMO convention for a protein: the peptide bond stays whole and
      !! each fragment is one residue's carbonyl with the next residue's amine
      !! and C-alpha.
      integer :: owner(24)

      owner = [1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 3, 3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3]
   end function owner_c_alpha

   subroutine glycine_tripeptide(sys)
      !! Glycine tripeptide, the geometry shipped as `sample_inputs/gly3.xyz`
      !!
      !! Atoms in the file's order: for each residue N, C-alpha, C(=O), O, then
      !! its hydrogens, and a terminal OH on the last carbonyl.
      type(system_geometry_t), intent(out) :: sys
      real(dp) :: xyz(3, 24)

      xyz = reshape([ &
                    0.0172_dp, -0.4777_dp, -0.0078_dp, &    ! 1  N
                    1.3251_dp, 0.1638_dp, 0.0713_dp, &      ! 2  C alpha
                    1.8818_dp, 0.1765_dp, 1.4668_dp, &      ! 3  C carbonyl
                    1.1564_dp, 0.4759_dp, 2.4031_dp, &      ! 4  O
                    2.0041_dp, -0.3893_dp, -0.6156_dp, &    ! 5  H
                    1.2934_dp, 1.2141_dp, -0.2903_dp, &     ! 6  H
                    -0.6558_dp, -0.0682_dp, 0.6786_dp, &    ! 7  H
                    -0.3827_dp, -0.2692_dp, -0.9506_dp, &   ! 8  H
                    3.2094_dp, -0.0781_dp, 1.6702_dp, &     ! 9  N
                    3.8490_dp, -0.0589_dp, 2.9843_dp, &     ! 10 C alpha
                    5.3502_dp, -0.0789_dp, 2.9477_dp, &     ! 11 C carbonyl
                    5.9543_dp, -0.1657_dp, 1.8893_dp, &     ! 12 O
                    3.5421_dp, 0.8561_dp, 3.5394_dp, &      ! 13 H
                    3.4987_dp, -0.9403_dp, 3.5644_dp, &     ! 14 H
                    3.7846_dp, -0.3120_dp, 0.8286_dp, &     ! 15 H
                    6.0352_dp, 0.0004_dp, 4.1282_dp, &      ! 16 N
                    7.4955_dp, -0.0139_dp, 4.2014_dp, &     ! 17 C alpha
                    8.0730_dp, 0.0278_dp, 5.5910_dp, &      ! 18 C carbonyl
                    7.3557_dp, 0.0642_dp, 6.5759_dp, &      ! 19 O
                    7.8695_dp, -0.9354_dp, 3.7022_dp, &     ! 20 H
                    7.8868_dp, 0.8596_dp, 3.6345_dp, &      ! 21 H
                    5.4671_dp, 0.0787_dp, 5.0035_dp, &      ! 22 H
                    9.3769_dp, 0.0222_dp, 5.7818_dp, &      ! 23 O
                    9.9377_dp, -0.0106_dp, 4.9381_dp], [3, 24])   ! 24 H

      sys%total_atoms = 24
      sys%n_monomers = 0
      allocate (sys%element_numbers(24))
      sys%element_numbers = [7, 6, 6, 8, 1, 1, 1, 1, 7, 6, 6, 8, 1, 1, 1, &
                             7, 6, 6, 8, 1, 1, 1, 8, 1]
      allocate (sys%coordinates(3, 24))
      sys%coordinates = to_bohr(xyz)
   end subroutine glycine_tripeptide

   subroutine propane(sys)
      !! Idealised propane: C-C 1.526, C-H 1.090, C-C-C 112 degrees
      !!
      !! Atom order is the two carbons of the cut bond, the third carbon, then
      !! hydrogens, so that a partition can be written down by eye.
      type(system_geometry_t), intent(out) :: sys
      real(dp) :: xyz(3, 11)

      xyz = reshape([1.5260_dp, 0.0000_dp, 0.0000_dp, &     ! C1
                     0.0000_dp, 0.0000_dp, 0.0000_dp, &     ! C2
                     -0.5716_dp, 1.4149_dp, 0.0000_dp, &    ! C3
                     2.1553_dp, -0.8900_dp, 0.0000_dp, &    ! H on C1
                     2.1553_dp, 0.4450_dp, -0.7707_dp, &    ! H on C1
                     2.1553_dp, 0.4450_dp, 0.7707_dp, &     ! H on C1
                     -0.3519_dp, -0.5217_dp, 0.8900_dp, &   ! H on C2
                     -0.3519_dp, -0.5217_dp, -0.8900_dp, &  ! H on C2
                     0.0178_dp, 2.3318_dp, 0.0000_dp, &     ! H on C3
                     -1.2200_dp, 1.8317_dp, -0.7707_dp, &   ! H on C3
                     -1.2200_dp, 1.8317_dp, 0.7707_dp], [3, 11])

      sys%total_atoms = 11
      sys%n_monomers = 0
      allocate (sys%element_numbers(11))
      sys%element_numbers = [6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1]
      allocate (sys%coordinates(3, 11))
      sys%coordinates = to_bohr(xyz)
   end subroutine propane

end module test_mqc_afo_model

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_afo_model, only: collect_mqc_afo_model
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_afo_model", collect_mqc_afo_model) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
