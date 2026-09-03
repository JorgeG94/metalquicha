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
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: run_libcint_rhf, rhf_result_t
   use mqc_libcint_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
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
                  new_unittest("afo_is_refused_with_an_embedding_field", test_refuse_esp), &
                  new_unittest("a_ring_cut_is_refused_by_name", test_refuse_ring), &
                  new_unittest("three_fragments_at_full_order_are_exact", test_three_exact), &
                  new_unittest("truncating_at_pairs_costs_the_three_body_term", test_three_pairs), &
                  new_unittest("one_atom_detached_from_two_bonds_is_ghosted_once", test_shared_bda) &
                  ]
   end subroutine collect_mqc_afo_fmo

   subroutine test_exact(error)
      !! Ethane as two methyls must equal ethane
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(libcint_molecule_t) :: mol
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

      call build_libcint_molecule(z, sym, xyz, "sto-3g", mol, err)
      call run_libcint_rhf(mol, 18, 200, 1.0e-11_dp, 1.0e-9_dp, .false., whole, err)
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

      call three_fragment_error(3, error, 1.0e-9_dp, &
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

      call three_fragment_error(2, error, 0.5_dp, &
                                "pairwise truncation is further off than the three-body "// &
                                "term explains")
   end subroutine test_three_pairs

   subroutine three_fragment_error(level, error, tol, what)
      !! Propane as three fragments at the given order, against the whole molecule
      integer, intent(in) :: level
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: tol
      character(len=*), intent(in) :: what

      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: whole
      integer :: z(11)
      character(len=2) :: sym(11)
      real(dp) :: xyz(3, 11)

      call propane(z, sym, xyz)
      opts%basis = "sto-3g"
      opts%esp = "none"
      opts%expansion = "mbe"
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

      call build_libcint_molecule(z, sym, xyz, "sto-3g", mol, err)
      call run_libcint_rhf(mol, 26, 200, 1.0e-11_dp, 1.0e-9_dp, .false., whole, err)
      call check(error,.not. err%has_error(), "the reference calculation failed")
      if (allocated(error)) return

      write (*, *) "   level", level, " error (hartree) =", res%energy - whole%energy
      call check(error, abs(res%energy - whole%energy) < tol, what)
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
      type(libcint_molecule_t) :: mol
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

      call build_libcint_molecule(z, sym, xyz, "sto-3g", mol, err)
      call run_libcint_rhf(mol, 26, 200, 1.0e-11_dp, 1.0e-9_dp, .false., whole, err)
      call check(error,.not. err%has_error(), "the reference calculation failed")
      if (allocated(error)) return

      write (*, *) "   shared-detached-atom error (hartree) =", res%energy - whole%energy
      call check(error, abs(res%energy - whole%energy) < 1.0e-9_dp, &
                 "a group ghosting one atom for two boundaries did not reproduce the "// &
                 "whole molecule at full order")
   end subroutine test_shared_bda

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
      !! A frozen orbital and an embedding field both describe the bond
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
