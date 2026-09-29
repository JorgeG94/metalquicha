!! Empirical dispersion under FMO and EE-MBE, added per fragment and per n-mer
module test_mqc_fmo_dispersion
   !! GAMESS applies Grimme's correction to every monomer and every n-mer on
   !! that group's own atoms (`DFTDSM` in `dftdis.src`) and never distributes
   !! the whole system's. So a pair term carries `E_D(IJ) - E_D(I) - E_D(J)`,
   !! and what these tests check is that, at every level of the expansion:
   !!
   !! * at full level the corrections telescope to `E_D` of the whole system,
   !!   so the total is the supermolecular Kohn-Sham energy plus that;
   !! * below full level the change the correction makes to the total is
   !!   exactly `sum_I E_D(I) + sum_IJ [E_D(IJ) - E_D(I) - E_D(J)]`, computed
   !!   here directly on atom subsets, separated pairs included;
   !! * PIEDA reports the same interaction as `Edi`, inside the pair energy.
   !!
   !! Every case skips when the build lacks libxc or the dispersion library.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_fragment, only: to_bohr
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: run_czt_rhf, rhf_result_t
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   use mqc_dispersion_apply, only: dispersion_apply, dispersion_kind_available
   implicit none
   private

   public :: collect_mqc_fmo_dispersion

   real(dp), parameter :: TOL_SUPER = 1.0e-8_dp
      !! A fragment expansion against the supermolecule: two SCFs on the same
      !! atoms and grid, as in `test_mqc_fmo_dft`.
   real(dp), parameter :: TOL_INCREMENT = 1.0e-10_dp
      !! The same SCFs with and without the correction: the density does not
      !! see it, so only rounding separates them.
   integer, parameter :: GRID_LEVEL = 3
   real(dp), parameter :: RESDIM_ONE_FAR = 1.5_dp
      !! On the stacked waters: the 1-2 and 2-3 pairs are inside it and the 1-3
      !! pair, twice as far, is beyond it. Checked by the tests that use it.

contains

   subroutine collect_mqc_fmo_dispersion(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fmo2_water_dimer_pbe_d3bj_is_the_supermolecule", test_dimer), &
                  new_unittest("fmo3_water_trimer_pbe_d3bj_is_the_supermolecule", test_trimer), &
                  new_unittest("propane_cut_pbe_d3bj_is_the_supermolecule", test_propane_cut), &
                  new_unittest("fmo2_cyclic_trimer_changes_by_the_per_group_sum", test_fmo2_sum), &
                  new_unittest("separated_pair_carries_its_dispersion_increment", &
                               test_separated_pair), &
                  new_unittest("fmo3_with_a_separated_pair_still_telescopes", &
                               test_fmo3_separated), &
                  new_unittest("eembe_trimer_changes_by_the_per_group_sum", test_eembe_sum), &
                  new_unittest("pieda_reports_the_dispersion_as_edi_inside_the_pair", &
                               test_pieda), &
                  new_unittest("dispersion_with_pieda_dispersion_is_refused", &
                               test_pieda_dispersion_refused), &
                  new_unittest("d4_takes_each_fragments_declared_charge", test_d4_charges), &
                  new_unittest("print_per_group_dispersion_for_gamess", test_print_gamess) &
                  ]
   end subroutine collect_mqc_fmo_dispersion

   function ready(kind) result(ok)
      !! Whether this build can run a Kohn-Sham fragment with `kind`
      character(len=*), intent(in) :: kind
      logical :: ok

      ok = xc_available() .and. dispersion_kind_available(kind)
   end function ready

   subroutine base_options(opts, kind)
      !! PBE, sto-3g, tight fragment SCFs, dispersion `kind` ("none" for off)
      type(fmo_options_t), intent(out) :: opts
      character(len=*), intent(in) :: kind

      opts%basis = "sto-3g"
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-11_dp
      opts%scf_density_tol = 1.0e-9_dp
      opts%outer_tol = 1.0e-10_dp
      opts%method%functional = "pbe"
      opts%method%grid_level = GRID_LEVEL
      opts%dispersion = kind
   end subroutine base_options

   function disp(kind, charge, z, xyz, atoms) result(energy)
      !! The correction of the atoms `atoms` taken as a molecule of their own
      !!
      !! `huge` when the library refuses, so that any comparison fails loudly.
      character(len=*), intent(in) :: kind
      real(dp), intent(in) :: charge
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: atoms(:)
      real(dp) :: energy

      type(error_t) :: derr

      call dispersion_apply(kind, "pbe", charge, z(atoms), xyz(:, atoms), energy, error=derr)
      if (derr%has_error()) energy = huge(1.0_dp)
   end function disp

   function members(owner, ids) result(atoms)
      !! The atoms whose fragment is one of `ids`, in system order
      integer, intent(in) :: owner(:), ids(:)
      integer, allocatable :: atoms(:)

      integer :: i

      atoms = pack([(i, i=1, size(owner))], [(any(owner(i) == ids), i=1, size(owner))])
   end function members

   function per_group_sum(kind, z, xyz, owner, n_frag, level) result(total)
      !! `sum_I D(I) + sum_IJ [D(IJ) - D(I) - D(J)]` at level two, plus the
      !! three-body remainder at level three, each on the group's own atoms
      character(len=*), intent(in) :: kind
      integer, intent(in) :: z(:), owner(:), n_frag, level
      real(dp), intent(in) :: xyz(:, :)
      real(dp) :: total

      real(dp) :: d1(n_frag), d2(n_frag, n_frag), d3
      integer :: i, j, k

      total = 0.0_dp
      do i = 1, n_frag
         d1(i) = disp(kind, 0.0_dp, z, xyz, members(owner, [i]))
      end do
      total = sum(d1)
      d2 = 0.0_dp
      do i = 1, n_frag
         do j = i + 1, n_frag
            d2(i, j) = disp(kind, 0.0_dp, z, xyz, members(owner, [i, j]))
            total = total + d2(i, j) - d1(i) - d1(j)
         end do
      end do
      if (level < 3) return
      do i = 1, n_frag
         do j = i + 1, n_frag
            do k = j + 1, n_frag
               d3 = disp(kind, 0.0_dp, z, xyz, members(owner, [i, j, k]))
               total = total + d3 - d2(i, j) - d2(i, k) - d2(j, k) + d1(i) + d1(j) + d1(k)
            end do
         end do
      end do
   end function per_group_sum

   subroutine test_dimer(error)
      !! FMO2 of two waters is exact, so PBE-D3(BJ) is the supermolecule's PBE
      !! plus D3(BJ) of the whole dimer
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6), whole, d_whole

      if (.not. ready("d3bj")) return
      call water_stack(2, z, sym, xyz)
      call base_options(opts, "d3bj")
      opts%level = 2

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO2 run failed: "//err%get_message())
      if (allocated(error)) return
      call supermolecule_energy(z, sym, xyz, 20, whole, err)
      d_whole = disp("d3bj", 0.0_dp, z, xyz, [1, 2, 3, 4, 5, 6])
      call check(error,.not. err%has_error(), "the reference failed: "//err%get_message())
      if (allocated(error)) return

      call check(error, abs(res%energy - (whole + d_whole)) < TOL_SUPER, &
                 "FMO2 PBE-D3(BJ) is not the supermolecule plus the whole dimer's D3(BJ)")
      write (*, *) "   fmo2 - (super + D) =", res%energy - (whole + d_whole), " D =", d_whole
   end subroutine test_dimer

   subroutine test_trimer(error)
      !! The same at three fragments and level three
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9), whole, d_whole

      if (.not. ready("d3bj")) return
      call water_stack(3, z, sym, xyz)
      call base_options(opts, "d3bj")
      opts%level = 3

      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2, 3, 3, 3], opts, res, err)
      call check(error,.not. err%has_error(), "the FMO3 run failed: "//err%get_message())
      if (allocated(error)) return
      call supermolecule_energy(z, sym, xyz, 30, whole, err)
      d_whole = disp("d3bj", 0.0_dp, z, xyz, [1, 2, 3, 4, 5, 6, 7, 8, 9])
      call check(error,.not. err%has_error(), "the reference failed: "//err%get_message())
      if (allocated(error)) return

      call check(error, abs(res%energy - (whole + d_whole)) < TOL_SUPER, &
                 "FMO3 PBE-D3(BJ) is not the supermolecule plus the whole trimer's D3(BJ)")
      write (*, *) "   fmo3 - (super + D) =", res%energy - (whole + d_whole), " D =", d_whole
   end subroutine test_trimer

   subroutine test_propane_cut(error)
      !! Propane cut across one C-C bond, in the frozen-orbital treatment: the
      !! dimer holds both ends, so the telescoping is exact and the dispersion,
      !! evaluated on real atoms only, gives the whole molecule's
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(11)
      character(len=2) :: sym(11)
      real(dp) :: xyz(3, 11), whole, d_whole
      integer :: i

      if (.not. ready("d3bj")) return
      call propane(z, sym, xyz)
      call base_options(opts, "d3bj")
      opts%bond_breaking = "afo"
      opts%level = 2

      call run_fmo2(z, sym, xyz, [1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2], opts, res, err)
      call check(error,.not. err%has_error(), "the cut run failed: "//err%get_message())
      if (allocated(error)) return
      call supermolecule_energy(z, sym, xyz, 26, whole, err)
      d_whole = disp("d3bj", 0.0_dp, z, xyz, [(i, i=1, 11)])
      call check(error,.not. err%has_error(), "the reference failed: "//err%get_message())
      if (allocated(error)) return

      call check(error, abs(res%energy - (whole + d_whole)) < TOL_SUPER, &
                 "cut propane PBE-D3(BJ) is not the supermolecule plus its D3(BJ)")
      write (*, *) "   cut - (super + D) =", res%energy - (whole + d_whole), " D =", d_whole
   end subroutine test_propane_cut

   subroutine test_fmo2_sum(error)
      !! Cyclic water trimer, FMO2, every pair solved: below full level, so the
      !! total is not the supermolecule's and the check is the per-group sum
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: with, without
      integer :: z(9), owner(9), k, i, j
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9), expect, incr

      if (.not. ready("d3bj")) return
      call water_cyclic(z, sym, xyz)
      owner = [1, 1, 1, 2, 2, 2, 3, 3, 3]

      call base_options(opts, "none")
      opts%level = 2
      opts%resppc = -1.0_dp
      call run_fmo2(z, sym, xyz, owner, opts, without, err)
      call check(error,.not. err%has_error(), "the run without D failed: "//err%get_message())
      if (allocated(error)) return
      opts%dispersion = "d3bj"
      call run_fmo2(z, sym, xyz, owner, opts, with, err)
      call check(error,.not. err%has_error(), "the run with D failed: "//err%get_message())
      if (allocated(error)) return

      expect = per_group_sum("d3bj", z, xyz, owner, 3, 2)
      call check(error,.not. err%has_error(), "the reference failed: "//err%get_message())
      if (allocated(error)) return
      call check(error, abs((with%energy - without%energy) - expect) < TOL_INCREMENT, &
                 "the change in the total is not sum D(I) + sum [D(IJ) - D(I) - D(J)]")
      write (*, *) "   change - per-group sum =", (with%energy - without%energy) - expect, &
         " sum =", expect
      if (allocated(error)) return

      ! Each monomer and each pair carries its own share.
      do k = 1, 3
         call check(error, abs((with%monomer_energy(k) - without%monomer_energy(k)) &
                               - disp("d3bj", 0.0_dp, z, xyz, members(owner, [k]))) &
                    < TOL_INCREMENT, "a monomer's energy does not carry its own D3(BJ)")
         if (allocated(error)) return
      end do
      do k = 1, size(with%pairs)
         i = with%pairs(k)%i
         j = with%pairs(k)%j
         incr = disp("d3bj", 0.0_dp, z, xyz, members(owner, [i, j])) &
                - disp("d3bj", 0.0_dp, z, xyz, members(owner, [i])) &
                - disp("d3bj", 0.0_dp, z, xyz, members(owner, [j]))
         call check(error, abs((with%pairs(k)%energy - without%pairs(k)%energy) - incr) &
                    < TOL_INCREMENT, "a pair's term does not carry D(IJ) - D(I) - D(J)")
         if (allocated(error)) return
      end do
   end subroutine test_fmo2_sum

   subroutine test_separated_pair(error)
      !! Stacked trimer with `resdim` between the near pairs and the far one:
      !! the far pair runs no SCF, and its term must still carry its increment
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: with, without
      integer :: z(9), owner(9), k
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9), expect

      if (.not. ready("d3bj")) return
      call water_stack(3, z, sym, xyz)
      owner = [1, 1, 1, 2, 2, 2, 3, 3, 3]

      call base_options(opts, "none")
      opts%level = 2
      opts%resdim = RESDIM_ONE_FAR
      call run_fmo2(z, sym, xyz, owner, opts, without, err)
      call check(error,.not. err%has_error(), "the run without D failed: "//err%get_message())
      if (allocated(error)) return
      opts%dispersion = "d3bj"
      call run_fmo2(z, sym, xyz, owner, opts, with, err)
      call check(error,.not. err%has_error(), "the run with D failed: "//err%get_message())
      if (allocated(error)) return

      call check(error, count(with%pairs%separated) == 1, &
                 "the setup should leave exactly one separated pair")
      if (allocated(error)) return
      expect = per_group_sum("d3bj", z, xyz, owner, 3, 2)
      call check(error,.not. err%has_error(), "the reference failed: "//err%get_message())
      if (allocated(error)) return
      call check(error, abs((with%energy - without%energy) - expect) < TOL_INCREMENT, &
                 "with a separated pair the change is not the per-group sum")
      write (*, *) "   change - per-group sum =", (with%energy - without%energy) - expect
      if (allocated(error)) return
      do k = 1, size(with%pairs)
         if (.not. with%pairs(k)%separated) cycle
         call check(error, abs(with%pairs(k)%energy - without%pairs(k)%energy) > 1.0e-6_dp, &
                    "the separated pair's term carries no dispersion")
      end do
   end subroutine test_separated_pair

   subroutine test_fmo3_separated(error)
      !! FMO3 on three fragments with one pair separated is not the
      !! supermolecule -- the separated pair is electrostatics -- but the
      !! dispersion still telescopes to the whole trimer's, which is what
      !! including the separated pair's increment makes true
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: with, without
      integer :: z(9), owner(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9), d_whole

      if (.not. ready("d3bj")) return
      call water_stack(3, z, sym, xyz)
      owner = [1, 1, 1, 2, 2, 2, 3, 3, 3]

      call base_options(opts, "none")
      opts%level = 3
      opts%resdim = RESDIM_ONE_FAR
      call run_fmo2(z, sym, xyz, owner, opts, without, err)
      call check(error,.not. err%has_error(), "the run without D failed: "//err%get_message())
      if (allocated(error)) return
      opts%dispersion = "d3bj"
      call run_fmo2(z, sym, xyz, owner, opts, with, err)
      call check(error,.not. err%has_error(), "the run with D failed: "//err%get_message())
      if (allocated(error)) return

      d_whole = disp("d3bj", 0.0_dp, z, xyz, [1, 2, 3, 4, 5, 6, 7, 8, 9])
      call check(error, abs((with%energy - without%energy) - d_whole) < TOL_INCREMENT, &
                 "FMO3 with a separated pair: the dispersion is not the whole system's")
      write (*, *) "   change - D(whole) =", (with%energy - without%energy) - d_whole
   end subroutine test_fmo3_separated

   subroutine test_eembe_sum(error)
      !! EE-MBE (point-charge field, total energies): the same per-group sum
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: with, without
      integer :: z(9), owner(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9), expect

      if (.not. ready("d3bj")) return
      call water_stack(3, z, sym, xyz)
      owner = [1, 1, 1, 2, 2, 2, 3, 3, 3]

      call base_options(opts, "none")
      opts%level = 2
      opts%esp = "ptc"
      opts%expansion = "mbe"
      call run_fmo2(z, sym, xyz, owner, opts, without, err)
      call check(error,.not. err%has_error(), "the run without D failed: "//err%get_message())
      if (allocated(error)) return
      opts%dispersion = "d3bj"
      call run_fmo2(z, sym, xyz, owner, opts, with, err)
      call check(error,.not. err%has_error(), "the run with D failed: "//err%get_message())
      if (allocated(error)) return

      expect = per_group_sum("d3bj", z, xyz, owner, 3, 2)
      call check(error,.not. err%has_error(), "the reference failed: "//err%get_message())
      if (allocated(error)) return
      call check(error, abs((with%energy - without%energy) - expect) < TOL_INCREMENT, &
                 "EE-MBE: the change in the total is not the per-group sum")
      write (*, *) "   change - per-group sum =", (with%energy - without%energy) - expect
   end subroutine test_eembe_sum

   subroutine test_pieda(error)
      !! PIEDA with dispersion on: `Edi` is the dispersion interaction, inside
      !! the pair energy; the other three terms are what they were without it
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: with, without
      integer :: z(9), owner(9), k, i, j
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9), edi

      if (.not. ready("d3bj")) return
      call water_stack(3, z, sym, xyz)
      owner = [1, 1, 1, 2, 2, 2, 3, 3, 3]

      call base_options(opts, "none")
      opts%level = 2
      opts%pieda = .true.
      opts%resdim = RESDIM_ONE_FAR
      call run_fmo2(z, sym, xyz, owner, opts, without, err)
      call check(error,.not. err%has_error(), "the run without D failed: "//err%get_message())
      if (allocated(error)) return
      call check(error,.not. without%edi_in_energy, "Edi is in the energy without any")
      if (allocated(error)) return
      opts%dispersion = "d3bj"
      call run_fmo2(z, sym, xyz, owner, opts, with, err)
      call check(error,.not. err%has_error(), "the run with D failed: "//err%get_message())
      if (allocated(error)) return
      call check(error, with%edi_in_energy, "the dispersion Edi should be inside the energy")
      if (allocated(error)) return

      do k = 1, size(with%pairs)
         i = with%pairs(k)%i
         j = with%pairs(k)%j
         edi = disp("d3bj", 0.0_dp, z, xyz, members(owner, [i, j])) &
               - disp("d3bj", 0.0_dp, z, xyz, members(owner, [i])) &
               - disp("d3bj", 0.0_dp, z, xyz, members(owner, [j]))
         call check(error, with%pairs(k)%pieda, "a pair was not decomposed")
         if (allocated(error)) return
         call check(error, abs(with%pairs(k)%edi - edi) < 1.0e-12_dp, &
                    "Edi is not D(IJ) - D(I) - D(J)")
         if (allocated(error)) return
         call check(error, abs(with%pairs(k)%ees + with%pairs(k)%eex + with%pairs(k)%ect_mix &
                               + with%pairs(k)%edi - with%pairs(k)%energy) < 1.0e-12_dp, &
                    "the PIEDA terms do not sum to the pair energy")
         if (allocated(error)) return
         call check(error, abs((with%pairs(k)%energy - without%pairs(k)%energy) - edi) &
                    < TOL_INCREMENT, "the pair energy did not move by Edi")
         if (allocated(error)) return
         call check(error, abs(with%pairs(k)%ees - without%pairs(k)%ees) < TOL_INCREMENT, &
                    "dispersion moved Ees")
         if (allocated(error)) return
         call check(error, abs(with%pairs(k)%eex - without%pairs(k)%eex) < TOL_INCREMENT, &
                    "dispersion moved Eex")
         if (allocated(error)) return
         call check(error, abs(with%pairs(k)%ect_mix - without%pairs(k)%ect_mix) < TOL_INCREMENT, &
                    "dispersion moved Ect+mix")
         if (allocated(error)) return
         write (*, "(a,i0,a,i0,a,es12.4,a,l1)") "   pair ", i, "-", j, "  Edi =", &
            with%pairs(k)%edi, "  separated ", with%pairs(k)%separated
      end do
   end subroutine test_pieda

   subroutine test_pieda_dispersion_refused(error)
      !! The pair's dispersion is already in its energy, so PIEDA's own would
      !! count it twice
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: xyz(3, 6)

      if (.not. ready("d3bj")) return
      call water_stack(2, z, sym, xyz)
      call base_options(opts, "d3bj")
      opts%pieda = .true.
      opts%pieda_dispersion = "d3bj"
      call run_fmo2(z, sym, xyz, [1, 1, 1, 2, 2, 2], opts, res, err)
      call check(error, err%has_error(), "both dispersions were accepted")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "pieda_dispersion") > 0, &
                 "the refusal does not name pieda_dispersion")
   end subroutine test_pieda_dispersion_refused

   subroutine test_d4_charges(error)
      !! An oxonium ion and a water under D4: the monomers carry their own
      !! charges, and at full level the total is the supermolecule's plus D4 of
      !! the whole ion
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(fmo_options_t) :: opts, plain
      type(fmo_result_t) :: res, res_plain
      integer :: z(7)
      character(len=2) :: sym(7)
      real(dp) :: xyz(3, 7), whole, d_whole, d_ion, d_ion_neutral
      integer, parameter :: OWNER(7) = [1, 1, 1, 1, 2, 2, 2]

      if (.not. ready("d4")) return
      call oxonium_water(z, sym, xyz)

      call base_options(opts, "d4")
      opts%level = 2
      opts%net_charge = [1, 0]
      call run_fmo2(z, sym, xyz, OWNER, opts, res, err)
      call check(error,.not. err%has_error(), "the D4 run failed: "//err%get_message())
      if (allocated(error)) return
      plain = opts
      plain%dispersion = "none"
      call run_fmo2(z, sym, xyz, OWNER, plain, res_plain, err)
      call check(error,.not. err%has_error(), "the run without D failed: "//err%get_message())
      if (allocated(error)) return

      call supermolecule_energy(z, sym, xyz, 20, whole, err)
      d_whole = disp("d4", 1.0_dp, z, xyz, [1, 2, 3, 4, 5, 6, 7])
      d_ion = disp("d4", 1.0_dp, z, xyz, [1, 2, 3, 4])
      d_ion_neutral = disp("d4", 0.0_dp, z, xyz, [1, 2, 3, 4])
      call check(error,.not. err%has_error(), "the reference failed: "//err%get_message())
      if (allocated(error)) return

      call check(error, abs(res%energy - (whole + d_whole)) < TOL_SUPER, &
                 "FMO2 PBE-D4 is not the supermolecule plus D4 of the whole ion")
      if (allocated(error)) return
      call check(error, abs(d_ion - d_ion_neutral) > 1.0e-7_dp, &
                 "the ion's charge does not move D4: the test cannot tell")
      if (allocated(error)) return
      call check(error, abs((res%monomer_energy(1) - res_plain%monomer_energy(1)) - d_ion) &
                 < TOL_INCREMENT, "the ion was given a D4 charge other than +1")
      write (*, *) "   D4 ion (+1) =", d_ion, " (neutral) =", d_ion_neutral
   end subroutine test_d4_charges

   subroutine test_print_gamess(error)
      !! Our per-group dispersion on `w3_pbe_d3.inp`'s geometry, beside what
      !! GAMESS's log gave. Printed, not asserted: `idcver=3` is not shown to be
      !! the same variant as our "d3bj".
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer :: z(9)
      character(len=2) :: sym(9)
      real(dp) :: xyz(3, 9)
      integer :: owner(9)

      if (.not. dispersion_kind_available("d3bj")) return
      call water_cyclic(z, sym, xyz)
      owner = [1, 1, 1, 2, 2, 2, 3, 3, 3]

      write (*, *) "   D3(BJ), PBE, Hartree              mqc      GAMESS (log)"
      write (*, "(a,es14.5,a)") "    monomer 1                  ", &
         disp("d3bj", 0.0_dp, z, xyz, members(owner, [1])), "   -8.9e-06 (a lone water)"
      write (*, "(a,es14.5)") "    monomer 2                  ", &
         disp("d3bj", 0.0_dp, z, xyz, members(owner, [2]))
      write (*, "(a,es14.5)") "    monomer 3                  ", &
         disp("d3bj", 0.0_dp, z, xyz, members(owner, [3]))
      write (*, "(a,es14.5,a)") "    dimer 1-2                  ", &
         disp("d3bj", 0.0_dp, z, xyz, members(owner, [1, 2])), "   -7.45e-04 (a dimer)"
      write (*, "(a,es14.5)") "    dimer 1-3                  ", &
         disp("d3bj", 0.0_dp, z, xyz, members(owner, [1, 3]))
      write (*, "(a,es14.5)") "    dimer 2-3                  ", &
         disp("d3bj", 0.0_dp, z, xyz, members(owner, [2, 3]))
      write (*, "(a,es14.5,a)") "    trimer (supermolecule)     ", &
         disp("d3bj", 0.0_dp, z, xyz, members(owner, [1, 2, 3])), "   -2.225e-03"
      call check(error,.not. err%has_error(), "the dispersion failed: "//err%get_message())
   end subroutine test_print_gamess

   subroutine supermolecule_energy(z, sym, xyz, nelec, energy, error)
      !! An ordinary restricted PBE energy on the whole system
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: nelec
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: mol
      type(xc_context_t) :: xc
      type(rhf_result_t) :: scf

      energy = 0.0_dp
      call build_czt_molecule(z, sym, xyz, "sto-3g", mol, error)
      if (error%has_error()) return
      call xc_context_create(mol, "pbe", xc, error, level=GRID_LEVEL, polarized=.false.)
      if (error%has_error()) return
      call run_czt_rhf(mol, nelec, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf, error, xc=xc)
      call xc%destroy()
      if (error%has_error()) return
      energy = scf%energy
   end subroutine supermolecule_energy

   subroutine water_stack(n, z, sym, xyz)
      !! `n` waters stacked 2.9 A apart along z, in Bohr
      integer, intent(in) :: n
      integer, intent(out) :: z(3*n)
      character(len=2), intent(out) :: sym(3*n)
      real(dp), intent(out) :: xyz(3, 3*n)

      real(dp) :: ang(3, 3)
      integer :: w, k

      ang = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                     0.0_dp, -0.7572_dp, 0.5865_dp, &
                     0.0_dp, 0.7572_dp, 0.5865_dp], [3, 3])
      do w = 1, n
         do k = 1, 3
            z(3*(w - 1) + k) = merge(8, 1, k == 1)
            sym(3*(w - 1) + k) = merge("O ", "H ", k == 1)
            xyz(:, 3*(w - 1) + k) = ang(:, k)
            xyz(3, 3*(w - 1) + k) = xyz(3, 3*(w - 1) + k) + real(w - 1, dp)*2.9_dp
         end do
      end do
      xyz = to_bohr(xyz)
   end subroutine water_stack

   subroutine water_cyclic(z, sym, xyz)
      !! `sample_inputs/water3_cyclic.xyz`, the geometry of
      !! `tools/fmo_validation/gamess/w3_pbe_d3.inp`
      integer, intent(out) :: z(9)
      character(len=2), intent(out) :: sym(9)
      real(dp), intent(out) :: xyz(3, 9)

      z = [8, 1, 1, 8, 1, 1, 8, 1, 1]
      sym = ["O ", "H ", "H ", "O ", "H ", "H ", "O ", "H ", "H "]
      xyz = to_bohr(reshape([-0.8309834169_dp, 1.3972845867_dp, -0.2245777899_dp, &
                             -1.7252647239_dp, 1.0653378323_dp, -0.1532773281_dp, &
                             -0.7392089607_dp, 2.0451809087_dp, 0.4613985156_dp, &
                             -0.3320345076_dp, -1.3821619786_dp, 0.2567741911_dp, &
                             -0.1833783620_dp, -0.4480121188_dp, 0.1145770742_dp, &
                             0.1095285385_dp, -1.8294391791_dp, -0.4524662378_dp, &
                             -3.0234372871_dp, -0.3756747342_dp, 0.2555351867_dp, &
                             -2.2955882864_dp, -0.9880884300_dp, 0.3537585219_dp, &
                             -3.5920719939_dp, -0.7444718872_dp, -0.4067791338_dp], [3, 9]))
   end subroutine water_cyclic

   subroutine oxonium_water(z, sym, xyz)
      !! H3O+ (atoms 1-4) hydrogen-bonded to a water (5-7), in Bohr
      integer, intent(out) :: z(7)
      character(len=2), intent(out) :: sym(7)
      real(dp), intent(out) :: xyz(3, 7)

      z = [8, 1, 1, 1, 8, 1, 1]
      sym = ["O ", "H ", "H ", "H ", "O ", "H ", "H "]
      xyz = to_bohr(reshape([0.0000_dp, 0.0000_dp, 0.1200_dp, &
                             0.9400_dp, 0.0000_dp, -0.2300_dp, &
                             -0.4700_dp, 0.8140_dp, -0.2300_dp, &
                             -0.4700_dp, -0.8140_dp, -0.2300_dp, &
                             -2.4000_dp, 0.0000_dp, -0.9000_dp, &
                             -2.9000_dp, 0.7600_dp, -1.2000_dp, &
                             -2.9000_dp, -0.7600_dp, -1.2000_dp], [3, 7]))
   end subroutine oxonium_water

   subroutine propane(z, sym, xyz)
      !! Idealised propane, carbons first -- as in `test_mqc_fmo_dft`
      integer, intent(out) :: z(11)
      character(len=2), intent(out) :: sym(11)
      real(dp), intent(out) :: xyz(3, 11)

      integer :: i

      z = [6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1]
      do i = 1, 11
         sym(i) = merge("C ", "H ", z(i) == 6)
      end do
      xyz = to_bohr(reshape([1.5260_dp, 0.0000_dp, 0.0000_dp, &
                             0.0000_dp, 0.0000_dp, 0.0000_dp, &
                             -0.5716_dp, 1.4149_dp, 0.0000_dp, &
                             2.1553_dp, -0.8900_dp, 0.0000_dp, &
                             2.1553_dp, 0.4450_dp, -0.7707_dp, &
                             2.1553_dp, 0.4450_dp, 0.7707_dp, &
                             -0.3519_dp, -0.5217_dp, 0.8900_dp, &
                             -0.3519_dp, -0.5217_dp, -0.8900_dp, &
                             0.0178_dp, 2.3318_dp, 0.0000_dp, &
                             -1.2200_dp, 1.8317_dp, -0.7707_dp, &
                             -1.2200_dp, 1.8317_dp, 0.7707_dp], [3, 11]))
   end subroutine propane

end module test_mqc_fmo_dispersion

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fmo_dispersion, only: collect_mqc_fmo_dispersion
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fmo_dispersion", collect_mqc_fmo_dispersion)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
