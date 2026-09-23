!! The frozen-orbital expansion, end to end through the FMO path
module test_mqc_afo_fmo
   !! One test carries this file, and it is the one the hydrogen caps would have
   !! failed on their first run.
   !!
   !! With **two** fragments an FMO2 expansion is not an approximation. The
   !! monomer terms cancel against the pair correction exactly --
   !! `E = E_1 + E_2 + (E_12 - E_1 - E_2) = E_12` -- and with nothing outside
   !! the pair there is no field and no response term either. So whatever the
   !! partition did, however the bond was detached and whatever was frozen where,
   !! the answer has to come back as an ordinary calculation on the whole
   !! molecule.
   !!
   !! That makes it a test of the bookkeeping alone, with the physics held out:
   !! the dimer holds both ends of the cut bond, so it must carry no ghost, no
   !! frozen orbital and no electron shift, while each monomer carries all
   !! three. Get the boundary set from a fragment instead of from the group and
   !! the dimer arrives holding a stand-in for a bond it contains, which is
   !! exactly how the capped version came back 11 Hartree adrift.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: to_bohr
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: run_czt_rhf, rhf_result_t
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   implicit none

   !! Both sides are the same SCF on the same molecule in the same basis, so
   !! what separates them is the assembly and not the convergence.
   real(dp), parameter :: TOL = 1.0e-9_dp

   private
   public :: collect_mqc_afo_fmo

contains

   subroutine collect_mqc_afo_fmo(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("two_fragments_across_a_cut_bond_are_exact", test_exact), &
                  new_unittest("afo_is_refused_with_an_exact_embedding", test_refuse_esp), &
                  new_unittest("a_ring_cut_is_refused_by_name", test_refuse_ring), &
                  new_unittest("three_fragments_at_full_order_are_exact", test_three_exact), &
                  new_unittest("truncating_at_pairs_costs_the_three_body_term", test_three_pairs), &
                  new_unittest("one_atom_detached_from_two_bonds_is_ghosted_once", test_shared_bda), &
                  new_unittest("two_fragments_keep_the_identity_under_point_charges", test_ptc_two), &
                  new_unittest("three_fragments_keep_it_with_a_field_on_top", test_ptc_three), &
                  new_unittest("embedded_truncation_at_pairs_costs_the_three_body_term", test_ptc_pairs), &
                  new_unittest("a_ring_kept_whole_takes_a_field_and_a_frozen_orbital", test_ptc_ring_kept) &
                  ]
   end subroutine collect_mqc_afo_fmo

   subroutine test_exact(error)
      !! Ethane as two methyls must equal ethane
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: whole
      integer :: z(8)
      character(len=2) :: sym(8)
      real(dp) :: xyz(3, 8)

      call ethane(z, sym, xyz)

      opts%basis = "sto-3g"
      opts%esp = "none"
      opts%expansion = "mbe"
      opts%bond_breaking = "afo"
      opts%level = 2
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp

      call run_fmo2(z, sym, xyz, [1, 1, 1, 1, 2, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the frozen-orbital expansion failed")
      if (allocated(error)) then
         if (err%has_error()) write (*, *) "   message: ", trim(err%get_message())
         return
      end if
      call check(error, res%converged, "the expansion did not converge")
      if (allocated(error)) return

      call build_czt_molecule(z, sym, xyz, "sto-3g", mol, err)
      call run_czt_rhf(mol, 18, 200, 1.0e-11_dp, 1.0e-9_dp, .false., whole, err)
      call check(error,.not. err%has_error(), "the reference calculation failed")
      if (allocated(error)) return

      call check(error, abs(res%energy - whole%energy) < TOL, &
                 "two fragments across a detached bond did not reproduce the whole "// &
                 "molecule -- the boundary bookkeeping is wrong, not the physics")
      if (allocated(error)) then
         write (*, *) "   fmo         =", res%energy
         write (*, *) "   unfragmented=", whole%energy
         write (*, *) "   difference  =", res%energy - whole%energy
      end if
   end subroutine test_exact

   subroutine test_three_exact(error)
      !! Three fragments, two detached bonds, expanded to full order
      !!
      !! The strongest statement this machinery can make about itself. An
      !! expansion carried to the fragment count is exact by inclusion and
      !! exclusion whatever the partition did, so the answer has to be the whole
      !! molecule -- and getting there exercises far more than the two-fragment
      !! case does. Three monomers each carry boundaries; three dimers do, one
      !! of them the pair of end fragments that are not bonded to each other and
      !! whose group holds a ghost of a carbon belonging to neither; and the
      !! trimer carries none. Ghosts, electron shifts, frozen occupied orbitals
      !! and frozen virtual ones all have to compose across those seven groups
      !! for the sum to land.
      type(error_type), allocatable, intent(out) :: error

      call three_fragment_error(3, error, 1.0e-9_dp, "none", &
                                "a full-order expansion over detached bonds did not "// &
                                "reproduce the whole molecule")
   end subroutine test_three_exact

   subroutine test_three_pairs(error)
      !! What truncating at pairs costs, measured rather than assumed
      !!
      !! Unlike the cases above this is a real approximation and the number is
      !! not small: about 0.18 Hartree on propane in STO-3G. That is the
      !! three-body term, and over covalent bonds it is large -- which is worth
      !! knowing rather than hiding, because the same quantity is a rounding
      !! error for a water cluster.
      !!
      !! It is an approximation and not a mistake, and the test above is what
      !! says so: the identical machinery at full order lands on the exact
      !! answer to 1e-13. Without that, a number this size would be
      !! indistinguishable from the 11 Hartree the capped attempt produced.
      type(error_type), allocatable, intent(out) :: error

      call three_fragment_error(2, error, 0.5_dp, "none", &
                                "pairwise truncation is further off than the three-body "// &
                                "term explains")
   end subroutine test_three_pairs

   subroutine three_fragment_error(level, error, tol, esp, what)
      !! Propane as three fragments at the given order, against the whole molecule
      integer, intent(in) :: level
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: tol
      character(len=*), intent(in) :: esp
         !! The embedding the fragments are solved in. `"none"` and `"ptc"` are
         !! the same expansion over differently polarised fragments, and at full
         !! order both have to land on the same number.
      character(len=*), intent(in) :: what

      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: whole
      integer :: z(11)
      character(len=2) :: sym(11)
      real(dp) :: xyz(3, 11)

      call propane(z, sym, xyz)
      opts%basis = "sto-3g"
      opts%esp = esp
      opts%expansion = merge("mbe", "fmo", esp == "none")
      opts%bond_breaking = "afo"
      opts%level = level
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp

      call run_fmo2(z, sym, xyz, [1, 2, 3, 1, 1, 1, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "the expansion failed")
      if (allocated(error)) then
         if (err%has_error()) write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call build_czt_molecule(z, sym, xyz, "sto-3g", mol, err)
      call run_czt_rhf(mol, 26, 200, 1.0e-11_dp, 1.0e-9_dp, .false., whole, err)
      call check(error,.not. err%has_error(), "the reference calculation failed")
      if (allocated(error)) return

      write (*, *) "   esp ", esp, " level", level, " error (hartree) =", &
         res%energy - whole%energy
      call check(error, abs(res%energy - whole%energy) < tol, what)
      if (allocated(error)) return
      if (esp /= "none") call check_charge_sum(res, 0.0_dp, error)
   end subroutine three_fragment_error

   subroutine test_shared_bda(error)
      !! The same atom detached from two bonds at once
      !!
      !! Number propane with its middle carbon first and that carbon is the
      !! detached end of *both* C-C bonds. The dimer of the two end methyls then
      !! holds the attached end of both, and has to bring the middle carbon in
      !! as a ghost twice over -- or rather must not: two copies of one atom's
      !! functions make the overlap exactly singular, and while the canonical
      !! orthogonalisation absorbs that and the answer survives, it survives by
      !! leaning on a linear-dependence threshold rather than by being right.
      !! One copy carries both hybrids, which are different orbitals of one atom
      !! and independent for that reason.
      !!
      !! Whether the situation arises at all depends on the order atoms were
      !! written down in, which is why it is worth a test of its own: the same
      !! molecule numbered the other way never reaches this path.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: whole
      integer :: z(11)
      character(len=2) :: sym(11)
      real(dp) :: xyz(3, 11)

      call propane_middle_first(z, sym, xyz)
      opts%basis = "sto-3g"
      opts%esp = "none"
      opts%expansion = "mbe"
      opts%bond_breaking = "afo"
      opts%level = 3
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp

      call run_fmo2(z, sym, xyz, [1, 2, 3, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "the expansion failed")
      if (allocated(error)) then
         if (err%has_error()) write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call build_czt_molecule(z, sym, xyz, "sto-3g", mol, err)
      call run_czt_rhf(mol, 26, 200, 1.0e-11_dp, 1.0e-9_dp, .false., whole, err)
      call check(error,.not. err%has_error(), "the reference calculation failed")
      if (allocated(error)) return

      write (*, *) "   shared-detached-atom error (hartree) =", res%energy - whole%energy
      call check(error, abs(res%energy - whole%energy) < 1.0e-9_dp, &
                 "a group ghosting one atom for two boundaries did not reproduce the "// &
                 "whole molecule at full order")
   end subroutine test_shared_bda

   subroutine test_ptc_two(error)
      !! Ethane as two methyls, with point charges and a frozen orbital at once
      !!
      !! The cheapest case in which an embedding field and a detached bond are
      !! both present, and the one that says whether they can coexist at all.
      !! Two fragments make the pair correction cancel the monomers exactly, so
      !! the answer is again the whole molecule however the monomers were
      !! polarised on the way -- the dimer holds both ends of the cut and sees
      !! nothing outside itself, so it carries neither a boundary nor a field.
      !!
      !! What the field adds over `esp = "none"` is everything the monomer loop
      !! does: charges read off a density that lives partly on an atom the
      !! fragment does not own, a field built from them, and an outer SCF that
      !! has to settle. None of that may move the total.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: whole
      integer :: z(8)
      character(len=2) :: sym(8)
      real(dp) :: xyz(3, 8)

      call ethane(z, sym, xyz)

      opts%basis = "sto-3g"
      opts%esp = "ptc"
      opts%expansion = "fmo"
      opts%bond_breaking = "afo"
      opts%level = 2
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp

      call run_fmo2(z, sym, xyz, [1, 1, 1, 1, 2, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the embedded frozen-orbital expansion failed")
      if (allocated(error)) then
         if (err%has_error()) write (*, *) "   message: ", trim(err%get_message())
         return
      end if
      call check(error, res%converged, "the embedded expansion did not converge")
      if (allocated(error)) return

      call build_czt_molecule(z, sym, xyz, "sto-3g", mol, err)
      call run_czt_rhf(mol, 18, 200, 1.0e-11_dp, 1.0e-9_dp, .false., whole, err)
      call check(error,.not. err%has_error(), "the reference calculation failed")
      if (allocated(error)) return

      write (*, *) "   ptc/afo two-fragment error (hartree) =", res%energy - whole%energy
      call check(error, abs(res%energy - whole%energy) < TOL, &
                 "a detached bond inside a point-charge field did not reproduce the "// &
                 "whole molecule")
      if (allocated(error)) then
         write (*, *) "   fmo         =", res%energy
         write (*, *) "   unfragmented=", whole%energy
         return
      end if

      call check_charge_sum(res, 0.0_dp, error)
   end subroutine test_ptc_two

   subroutine test_ptc_three(error)
      !! Propane at full order, embedded in point charges
      !!
      !! Three fragments and two detached bonds, so every group shape the
      !! machinery has appears: monomers carrying a boundary, a dimer of the two
      !! end fragments whose ghost carbon belongs to neither of them, and a
      !! trimer carrying none. With a field on top, each of those groups also
      !! has to be told which charges are somebody else's -- the atom a group
      !! holds as a ghost is one it already describes, and counting it in the
      !! field as well would describe it twice.
      type(error_type), allocatable, intent(out) :: error

      call three_fragment_error(3, error, 1.0e-9_dp, "ptc", &
                                "a full-order embedded expansion over detached bonds "// &
                                "did not reproduce the whole molecule")
   end subroutine test_ptc_three

   subroutine test_ptc_pairs(error)
      !! What truncating at pairs costs once the fragments are polarised
      !!
      !! The embedded counterpart of `test_three_pairs`, and the same
      !! reasoning: this one is a real approximation, so the number is printed
      !! and only loosely bounded. It is here to be read next to the full-order
      !! case above, which is what says a number of this size is the three-body
      !! term and not a bookkeeping fault.
      !!
      !! It comes out at 0.489 Hartree, against 0.180 unembedded, so the point
      !! charges make the *truncated* expansion worse. That is expected rather
      !! than wrong and the reason is in `fmo.rst`; the bound here is an order
      !! of magnitude and not a pin, because what it has to separate is a
      !! three-body term from an eleven-Hartree assembly fault.
      type(error_type), allocatable, intent(out) :: error

      call three_fragment_error(2, error, 1.0_dp, "ptc", &
                                "embedded pairwise truncation is further off than the "// &
                                "three-body term explains")
   end subroutine test_ptc_pairs

   subroutine test_ptc_ring_kept(error)
      !! Methylcyclopropane, cut at the bond that is not in the ring
      !!
      !! Cyclopropane itself cannot stand in here. Any partition of a three-ring
      !! severs two bonds between the same pair of fragments, which is the ring
      !! cut `test_a_ring_cut_is_refused_by_name` asserts is refused -- a frozen
      !! orbital stands in for one detached bond and not for two. So the ring is
      !! kept whole and the exocyclic bond is the one detached, which is the
      !! case a strained ring actually contributes: a fragment whose geometry is
      !! nothing like the model system the hybrid was taken off, carrying a
      !! field from its neighbour at the same time.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: whole
      integer :: z(12)
      character(len=2) :: sym(12)
      real(dp) :: xyz(3, 12)

      call methylcyclopropane(z, sym, xyz)

      opts%basis = "sto-3g"
      opts%esp = "ptc"
      opts%expansion = "fmo"
      opts%bond_breaking = "afo"
      opts%level = 2
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the embedded expansion failed")
      if (allocated(error)) then
         if (err%has_error()) write (*, *) "   message: ", trim(err%get_message())
         return
      end if

      call build_czt_molecule(z, sym, xyz, "sto-3g", mol, err)
      call run_czt_rhf(mol, 32, 200, 1.0e-11_dp, 1.0e-9_dp, .false., whole, err)
      call check(error,.not. err%has_error(), "the reference calculation failed")
      if (allocated(error)) return

      write (*, *) "   ptc/afo ring-kept error (hartree) =", res%energy - whole%energy
      call check(error, abs(res%energy - whole%energy) < TOL, &
                 "an exocyclic bond detached inside a point-charge field did not "// &
                 "reproduce the whole molecule")
      if (allocated(error)) return

      call check_charge_sum(res, 0.0_dp, error)
   end subroutine test_ptc_ring_kept

   subroutine check_charge_sum(res, total, error)
      !! Every electron is somebody's, and the atomic charges have to say so
      !!
      !! The population of a detached atom arrives from two fragments: the one
      !! that owns it, whose hybrid there is frozen empty, and the one that
      !! holds it as a ghost with the bond pair in that same hybrid. Those two
      !! shares are added. Drop either and the charges no longer sum to the
      !! molecular charge -- and since they are what the field is built from,
      !! every fragment would then be embedded in a system carrying a charge the
      !! molecule does not have.
      type(fmo_result_t), intent(in) :: res
      real(dp), intent(in) :: total
      type(error_type), allocatable, intent(out) :: error

      call check(error, allocated(res%charges), "an embedded run reported no charges")
      if (allocated(error)) return
      write (*, *) "   sum of atomic charges =", sum(res%charges)
      call check(error, abs(sum(res%charges) - total) < 1.0e-8_dp, &
                 "the fragment charges do not sum to the molecular charge, so the "// &
                 "detached atom's population was counted once instead of twice")
   end subroutine check_charge_sum

   subroutine methylcyclopropane(z, sym, xyz)
      !! Methylcyclopropane, carbons first: ring C1-C3, then the methyl carbon
      !!
      !! Idealised rather than optimised -- a 1.51 Angstrom equilateral ring in
      !! the xy plane, 1.09 Angstrom C-H, and the exocyclic bond raised out of
      !! the ring plane. The fragmented and unfragmented answers are compared at
      !! the same geometry, so what matters is that the connectivity is
      !! unambiguous and the ring is a ring.
      integer, intent(out) :: z(12)
      character(len=2), intent(out) :: sym(12)
      real(dp), intent(out) :: xyz(3, 12)
      real(dp) :: ang(3, 12)
      integer :: i

      z = [6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1]
      do i = 1, 12
         if (z(i) == 6) then
            sym(i) = "C "
         else
            sym(i) = "H "
         end if
      end do
      ang = reshape([0.8718_dp, 0.0000_dp, 0.0000_dp, &    ! C1, carries the methyl
                     -0.4359_dp, 0.7550_dp, 0.0000_dp, &   ! C2
                     -0.4359_dp, -0.7550_dp, 0.0000_dp, &  ! C3
                     1.7379_dp, 0.0000_dp, 1.2370_dp, &    ! C4, the methyl carbon
                     1.4970_dp, 0.0000_dp, -0.8929_dp, &   ! H on C1
                     -0.7485_dp, 1.2964_dp, 0.8929_dp, &   ! H on C2
                     -0.7485_dp, 1.2964_dp, -0.8929_dp, &
                     -0.7485_dp, -1.2964_dp, 0.8929_dp, &  ! H on C3
                     -0.7485_dp, -1.2964_dp, -0.8929_dp, &
                     1.9466_dp, 1.0274_dp, 1.5351_dp, &    ! H on C4
                     1.2178_dp, -0.5137_dp, 2.0454_dp, &
                     2.6755_dp, -0.5137_dp, 1.0248_dp], [3, 12])
      xyz = to_bohr(ang)
   end subroutine methylcyclopropane

   subroutine propane_middle_first(z, sym, xyz)
      !! The same propane, with the middle carbon numbered first
      integer, intent(out) :: z(11)
      character(len=2), intent(out) :: sym(11)
      real(dp), intent(out) :: xyz(3, 11)
      real(dp) :: ang(3, 11)
      integer :: i

      z = [6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1]
      do i = 1, 11
         if (z(i) == 6) then
            sym(i) = "C "
         else
            sym(i) = "H "
         end if
      end do
      ang = reshape([0.0000_dp, 0.0000_dp, 0.0000_dp, &     ! C middle
                     1.5260_dp, 0.0000_dp, 0.0000_dp, &     ! C end
                     -0.5716_dp, 1.4149_dp, 0.0000_dp, &    ! C end
                     -0.3519_dp, -0.5217_dp, 0.8900_dp, &   ! H on the middle
                     -0.3519_dp, -0.5217_dp, -0.8900_dp, &  ! H on the middle
                     2.1553_dp, -0.8900_dp, 0.0000_dp, &
                     2.1553_dp, 0.4450_dp, -0.7707_dp, &
                     2.1553_dp, 0.4450_dp, 0.7707_dp, &
                     0.0178_dp, 2.3318_dp, 0.0000_dp, &
                     -1.2200_dp, 1.8317_dp, -0.7707_dp, &
                     -1.2200_dp, 1.8317_dp, 0.7707_dp], [3, 11])
      xyz = to_bohr(ang)
   end subroutine propane_middle_first

   subroutine test_refuse_esp(error)
      !! An exact embedding stays refused, and for a reason that has not moved
      !!
      !! Point charges now run alongside a detached bond because the detached
      !! atom's share of the field is one number per atom there and comes back
      !! out exactly. The exact embedding builds a Coulomb contraction over the
      !! neighbour's whole density matrix, which has no per-atom part to remove,
      !! so there is nothing to subtract that would not be the point-charge
      !! approximation under another name.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(8)
      character(len=2) :: sym(8)
      real(dp) :: xyz(3, 8)

      call ethane(z, sym, xyz)
      opts%basis = "sto-3g"
      opts%esp = "exact"
      opts%bond_breaking = "afo"

      call run_fmo2(z, sym, xyz, [1, 1, 1, 1, 2, 2, 2, 2], opts, res, err)
      call check(error, err%has_error(), &
                 "a detached bond was accepted alongside an embedding field")
   end subroutine test_refuse_esp

   subroutine test_refuse_ring(error)
      !! Cyclopropane into three CH2 joins each pair of fragments twice
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      integer :: i

      z = [6, 1, 1, 6, 1, 1, 6, 1, 1]
      do i = 1, 9
         if (z(i) == 6) then
            sym(i) = "C "
         else
            sym(i) = "H "
         end if
      end do
      xyz = reshape([0.8718_dp, 0.0000_dp, 0.0000_dp, &
                     1.3818_dp, 0.0000_dp, 0.9500_dp, &
                     1.3818_dp, 0.0000_dp, -0.9500_dp, &
                     -0.4359_dp, 0.7550_dp, 0.0000_dp, &
                     -0.6909_dp, 1.1967_dp, 0.9500_dp, &
                     -0.6909_dp, 1.1967_dp, -0.9500_dp, &
                     -0.4359_dp, -0.7550_dp, 0.0000_dp, &
                     -0.6909_dp, -1.1967_dp, 0.9500_dp, &
                     -0.6909_dp, -1.1967_dp, -0.9500_dp], [3, 9])
      xyz = to_bohr(xyz)

      opts%basis = "sto-3g"
      opts%esp = "none"
      opts%expansion = "mbe"
      opts%bond_breaking = "afo"

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error, err%has_error(), "a ring cut was accepted")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "ring") > 0, &
                 "a ring cut was refused for some other reason than being a ring")
   end subroutine test_refuse_ring

   subroutine propane(z, sym, xyz)
      !! Idealised propane, carbons first so a partition reads by eye
      integer, intent(out) :: z(11)
      character(len=2), intent(out) :: sym(11)
      real(dp), intent(out) :: xyz(3, 11)
      real(dp) :: ang(3, 11)
      integer :: i

      z = [6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1]
      do i = 1, 11
         if (z(i) == 6) then
            sym(i) = "C "
         else
            sym(i) = "H "
         end if
      end do
      ang = reshape([1.5260_dp, 0.0000_dp, 0.0000_dp, &
                     0.0000_dp, 0.0000_dp, 0.0000_dp, &
                     -0.5716_dp, 1.4149_dp, 0.0000_dp, &
                     2.1553_dp, -0.8900_dp, 0.0000_dp, &
                     2.1553_dp, 0.4450_dp, -0.7707_dp, &
                     2.1553_dp, 0.4450_dp, 0.7707_dp, &
                     -0.3519_dp, -0.5217_dp, 0.8900_dp, &
                     -0.3519_dp, -0.5217_dp, -0.8900_dp, &
                     0.0178_dp, 2.3318_dp, 0.0000_dp, &
                     -1.2200_dp, 1.8317_dp, -0.7707_dp, &
                     -1.2200_dp, 1.8317_dp, 0.7707_dp], [3, 11])
      xyz = to_bohr(ang)
   end subroutine propane

   subroutine ethane(z, sym, xyz)
      integer, intent(out) :: z(8)
      character(len=2), intent(out) :: sym(8)
      real(dp), intent(out) :: xyz(3, 8)
      real(dp) :: ang(3, 8)
      integer :: i

      z = [6, 1, 1, 1, 6, 1, 1, 1]
      do i = 1, 8
         if (z(i) == 6) then
            sym(i) = "C "
         else
            sym(i) = "H "
         end if
      end do
      ang = reshape([0.000_dp, 0.000_dp, 0.768_dp, &
                     -1.019_dp, 0.000_dp, 1.157_dp, &
                     0.510_dp, 0.883_dp, 1.157_dp, &
                     0.510_dp, -0.883_dp, 1.157_dp, &
                     0.000_dp, 0.000_dp, -0.768_dp, &
                     1.019_dp, 0.000_dp, -1.157_dp, &
                     -0.510_dp, -0.883_dp, -1.157_dp, &
                     -0.510_dp, 0.883_dp, -1.157_dp], [3, 8])
      xyz = to_bohr(ang)
   end subroutine ethane

end module test_mqc_afo_fmo

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_afo_fmo, only: collect_mqc_afo_fmo
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_afo_fmo", collect_mqc_afo_fmo)]
   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
