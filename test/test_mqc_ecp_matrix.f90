module test_mqc_ecp_matrix
   !! The ECP matrix for a lanthanide, where the local channel is l = 5.
   !!
   !! Ce-Lu are the only elements def2-ECP tabulates with a local channel above
   !! l = 3, and they are a separate path through the angular code: five
   !! projected channels running l = 0 to 4 rather than three running 0 to 2.
   !! `validation/check_ecp.f90` covers iodine, which is the l = 3 case and the
   !! one thirty-six of the fifty elements use; this is the other fourteen.
   !!
   !! Tested at the matrix rather than through an SCF for cost. Ytterbium in
   !! def2-SVP is 81 basis functions with a g shell in it, and a Hartree-Fock
   !! energy on it runs for half a minute against the tenth of a second the
   !! rest of the validation suite averages -- too slow to sit in CI for
   !! coverage that a one-electron integral settles outright. What an SCF would
   !! add over this is the *rest* of the program, which every other ECP case
   !! already exercises.
   !!
   !! The reference is PySCF's `mol.intor('ECPscalar_sph')` on the same system,
   !! fed **our** `def2-ecp.json` through `bse_ecp_to_pyscf` in
   !! `tools/cpu_validation`. It has to be ours: PySCF ships no def2 ECP for the
   !! lanthanides at all, so there is no table of theirs to compare against.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_basis_utils, only: find_basis_file
   use mqc_cgto, only: molecular_basis_type
   use mqc_ecp, only: molecular_ecp_type
   use mqc_error, only: error_t
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use mqc_json_ecp_reader, only: build_molecular_ecp_json
   use mqc_czt_ecp, only: ecp_matrix, ecp_deriv_matrix
   use mqc_czt_integrals, only: czt_molecule_t, atom_ao_blocks
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_ecp_matrix_tests

   !! PySCF 2.14, `mol.intor('ECPscalar_sph')` on Yb plus a hydrogen 25 Bohr
   !! away -- far enough that the overlap screen has to reject as well as accept.
   real(dp), parameter :: REF_MAX = 8.43440930940146e+02_dp
   real(dp), parameter :: REF_TRACE = 2.78522522914258e+03_dp

   !! Measured agreement is 6e-12 on the trace and 1e-13 on the largest element,
   !! against references written to fifteen significant figures -- so the bound
   !! below carries two orders over what the comparison can resolve.
   !!
   !! It is calibrated against losing a channel, which is the fault this exists
   !! to catch and which is not equally visible in all of them: dropping
   !! ytterbium's l = 0 projector moves the trace by 2.5e+03 and the l = 3 one by
   !! 1.3e+01, but the l = 4 projector -- the one no lighter element has, and so
   !! the whole reason this case is here -- is worth only 3.0e-06. That is the
   !! number the tolerance has to sit under, and 1e-09 clears it by three orders
   !! where the 1e-06 this started at would have had a factor of three.
   real(dp), parameter :: TOL = 1.0e-9_dp

contains

   subroutine collect_mqc_ecp_matrix_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("lanthanide_local_channel_is_l5", test_channels), &
                  new_unittest("lanthanide_ecp_matrix_matches_pyscf", test_matrix), &
                  new_unittest("lanthanide_ecp_matrix_is_symmetric", test_symmetric), &
                  new_unittest("ecp_derivative_matches_differences", test_derivative), &
                  new_unittest("ecp_derivative_splits_by_atom", test_derivative_by_atom) &
                  ]
   end subroutine collect_mqc_ecp_matrix_tests

   subroutine build_ytterbium(basis, ecp, mol, error)
      !! Yb with a distant hydrogen, read from the shipped def2 files
      type(molecular_basis_type), intent(out) :: basis
      type(molecular_ecp_type), intent(out) :: ecp
      type(czt_molecule_t), intent(out) :: mol
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      character(len=:), allocatable :: basis_file, ecp_file
      character(len=8) :: symbols(2)
      real(dp) :: coordinates(3, 2)
      integer :: numbers(2)

      symbols = ["Yb", "H "]
      numbers = [70, 1]
      coordinates = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                             0.0_dp, 0.0_dp, 25.0_dp], [3, 2])

      call find_basis_file("def2-svp", basis_file, err)
      call check(error,.not. err%has_error(), "def2-svp not found")
      if (allocated(error)) return
      call find_basis_file("def2-ecp", ecp_file, err)
      call check(error,.not. err%has_error(), "def2-ecp not found")
      if (allocated(error)) return

      call build_molecular_basis_json(basis_file, symbols, basis, err)
      call check(error,.not. err%has_error(), "could not read def2-svp")
      if (allocated(error)) return
      call build_molecular_ecp_json(ecp_file, symbols, ecp, err)
      call check(error,.not. err%has_error(), "could not read def2-ecp")
      if (allocated(error)) return

      call mol%build(numbers, coordinates, basis, err, ecp=ecp)
      call check(error,.not. err%has_error(), "could not build the molecule")
   end subroutine build_ytterbium

   subroutine test_channels(error)
      !! The layout, before any integral. A channel that came back claiming the
      !! wrong l still integrates to something, and the something looks
      !! plausible -- asserting the layout separates a reader fault from an
      !! integral one.
      type(error_type), allocatable, intent(out) :: error
      type(molecular_basis_type) :: basis
      type(molecular_ecp_type) :: ecp
      type(czt_molecule_t) :: mol
      integer :: i

      call build_ytterbium(basis, ecp, mol, error)
      if (allocated(error)) return

      call check(error, ecp%atoms(1)%local%ang_mom, 5, &
                 "ytterbium's local channel should be l = 5")
      if (allocated(error)) return
      call check(error, ecp%atoms(1)%n_projected, 5, &
                 "ytterbium should carry 5 projected channels")
      if (allocated(error)) return
      do i = 1, 5
         call check(error, ecp%atoms(1)%projected(i)%ang_mom, i - 1, &
                    "the projected channels should run l = 0 to 4")
         if (allocated(error)) return
      end do

      call check(error, ecp%atoms(1)%core_electrons, 28, &
                 "ytterbium's def2-ECP replaces 28 electrons")
      if (allocated(error)) return
      call check(error,.not. ecp%atoms(2)%has_ecp, &
                 "hydrogen should carry no ECP")
   end subroutine test_channels

   subroutine test_matrix(error)
      type(error_type), allocatable, intent(out) :: error
      type(molecular_basis_type) :: basis
      type(molecular_ecp_type) :: ecp
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: m(:, :)
      real(dp) :: got_max, got_trace
      integer :: i

      call build_ytterbium(basis, ecp, mol, error)
      if (allocated(error)) return

      call check(error, mol%necpbas > 0, "no ECP rows were laid out")
      if (allocated(error)) return

      call ecp_matrix(mol%nao, mol%nbas, mol%natm, mol%cartesian, mol%atm, &
                      mol%bas_with_ecp, mol%env, mol%shell_offset, mol%necpbas, m)

      got_max = maxval(abs(m))
      got_trace = sum([(m(i, i), i=1, mol%nao)])

      call check(error, abs(got_max - REF_MAX) < TOL, &
                 "max |ECP| disagrees with PySCF")
      if (allocated(error)) return
      call check(error, abs(got_trace - REF_TRACE) < TOL, &
                 "the ECP trace disagrees with PySCF")
   end subroutine test_matrix

   subroutine test_symmetric(error)
      !! Both scalars above can match while a shell-pair block sits in the
      !! wrong place, which symmetry catches and they do not.
      type(error_type), allocatable, intent(out) :: error
      type(molecular_basis_type) :: basis
      type(molecular_ecp_type) :: ecp
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: m(:, :)

      call build_ytterbium(basis, ecp, mol, error)
      if (allocated(error)) return

      call ecp_matrix(mol%nao, mol%nbas, mol%natm, mol%cartesian, mol%atm, &
                      mol%bas_with_ecp, mol%env, mol%shell_offset, mol%necpbas, m)

      call check(error, maxval(abs(m - transpose(m))) < 1.0e-12_dp, &
                 "the ECP operator is symmetric and this matrix is not")
   end subroutine test_symmetric

   subroutine build_hydrogen_iodide(coordinates, basis, ecp, mol, error)
      !! HI at the given coordinates (Bohr), read from the shipped def2 files
      real(dp), intent(in) :: coordinates(3, 2)
      type(molecular_basis_type), intent(out) :: basis
      type(molecular_ecp_type), intent(out) :: ecp
      type(czt_molecule_t), intent(out) :: mol
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      character(len=:), allocatable :: basis_file, ecp_file
      character(len=8) :: symbols(2)

      symbols = ["I ", "H "]
      call find_basis_file("def2-svp", basis_file, err)
      call check(error,.not. err%has_error(), "def2-svp not found")
      if (allocated(error)) return
      call find_basis_file("def2-ecp", ecp_file, err)
      call check(error,.not. err%has_error(), "def2-ecp not found")
      if (allocated(error)) return

      call build_molecular_basis_json(basis_file, symbols, basis, err)
      call check(error,.not. err%has_error(), "could not read def2-svp")
      if (allocated(error)) return
      call build_molecular_ecp_json(ecp_file, symbols, ecp, err)
      call check(error,.not. err%has_error(), "could not read def2-ecp")
      if (allocated(error)) return

      call mol%build([53, 1], coordinates, basis, err, ecp=ecp)
      call check(error,.not. err%has_error(), "could not build the molecule")
   end subroutine build_hydrogen_iodide

   subroutine ecp_at(coordinates, m, error)
      !! The ECP matrix of HI at one geometry
      real(dp), intent(in) :: coordinates(3, 2)
      real(dp), allocatable, intent(out) :: m(:, :)
      type(error_type), allocatable, intent(out) :: error
      type(molecular_basis_type) :: basis
      type(molecular_ecp_type) :: ecp
      type(czt_molecule_t) :: mol

      call build_hydrogen_iodide(coordinates, basis, ecp, mol, error)
      if (allocated(error)) return
      call ecp_matrix(mol%nao, mol%nbas, mol%natm, mol%cartesian, mol%atm, &
                      mol%bas_with_ecp, mol%env, mol%shell_offset, mol%necpbas, m)
   end subroutine ecp_at

   subroutine test_derivative(error)
      !! The derivative integrals against the ECP matrix differenced in each
      !! atom's coordinates.
      !!
      !! Moving atom A moves the basis functions on it and, if it has one, its
      !! potential. With `D = <nabla mu|U|nu>` over every centre and `D_A` over
      !! atom A's alone,
      !!
      !!     dU/dR_A (mu, nu) = -D(mu, nu) [mu on A] - D(nu, mu) [nu on A]
      !!                        + D_A(mu, nu) + D_A(nu, mu)
      !!
      !! which is how every gradient in the program assembles it, so this is the
      !! check that the assembly and the integrals mean the same thing. The
      !! geometry is off the axes so that no component vanishes by symmetry.
      !! A four-point difference at 5e-3 Bohr leaves an error near 1e-9.
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: STEP = 5.0e-3_dp
      real(dp), parameter :: TOL_DERIV = 1.0e-7_dp
      type(molecular_basis_type) :: basis
      type(molecular_ecp_type) :: ecp
      type(czt_molecule_t) :: mol
      real(dp) :: at_rest(3, 2), moved(3, 2)
      real(dp) :: worst, scale
      real(dp), allocatable :: whole(:, :, :), one(:, :, :), analytic(:, :)
      real(dp), allocatable :: plus1(:, :), minus1(:, :), plus2(:, :), minus2(:, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: iatom, comp, p0, p1

      at_rest = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.6_dp, -0.4_dp, 2.9_dp], [3, 2])
      call build_hydrogen_iodide(at_rest, basis, ecp, mol, error)
      if (allocated(error)) return

      allocate (whole(mol%nao, mol%nao, 3), one(mol%nao, mol%nao, 3))
      allocate (analytic(mol%nao, mol%nao))
      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)
      call ecp_deriv_matrix(mol%nao, mol%nbas, mol%natm, mol%cartesian, mol%atm, &
                            mol%bas_with_ecp, mol%env, mol%shell_offset, mol%necpbas, &
                            0, whole)
      scale = maxval(abs(whole))
      call check(error, scale > 1.0_dp, "the derivative integrals are suspiciously small")
      if (allocated(error)) return

      worst = 0.0_dp
      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)
         call ecp_deriv_matrix(mol%nao, mol%nbas, mol%natm, mol%cartesian, mol%atm, &
                               mol%bas_with_ecp, mol%env, mol%shell_offset, mol%necpbas, &
                               iatom, one)
         do comp = 1, 3
            analytic = one(:, :, comp) + transpose(one(:, :, comp))
            analytic(p0:p1, :) = analytic(p0:p1, :) - whole(p0:p1, :, comp)
            analytic(:, p0:p1) = analytic(:, p0:p1) - transpose(whole(p0:p1, :, comp))

            moved = at_rest
            moved(comp, iatom) = at_rest(comp, iatom) + STEP
            call ecp_at(moved, plus1, error)
            if (allocated(error)) return
            moved = at_rest
            moved(comp, iatom) = at_rest(comp, iatom) - STEP
            call ecp_at(moved, minus1, error)
            if (allocated(error)) return
            moved = at_rest
            moved(comp, iatom) = at_rest(comp, iatom) + 2.0_dp*STEP
            call ecp_at(moved, plus2, error)
            if (allocated(error)) return
            moved = at_rest
            moved(comp, iatom) = at_rest(comp, iatom) - 2.0_dp*STEP
            call ecp_at(moved, minus2, error)
            if (allocated(error)) return

            worst = max(worst, maxval(abs(analytic - (8.0_dp*(plus1 - minus1) &
                                                      - (plus2 - minus2))/(12.0_dp*STEP))))
         end do
      end do

      call check(error, worst < TOL_DERIV*scale, &
                 "the ECP derivative integrals disagree with the differenced matrix")
   end subroutine test_derivative

   subroutine test_derivative_by_atom(error)
      !! An atom's own potential: iodine's is the whole of it, hydrogen's is
      !! nothing.
      type(error_type), allocatable, intent(out) :: error
      type(molecular_basis_type) :: basis
      type(molecular_ecp_type) :: ecp
      type(czt_molecule_t) :: mol
      real(dp) :: at_rest(3, 2)
      real(dp), allocatable :: whole(:, :, :), one(:, :, :)

      at_rest = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.6_dp, -0.4_dp, 2.9_dp], [3, 2])
      call build_hydrogen_iodide(at_rest, basis, ecp, mol, error)
      if (allocated(error)) return

      allocate (whole(mol%nao, mol%nao, 3), one(mol%nao, mol%nao, 3))
      call ecp_deriv_matrix(mol%nao, mol%nbas, mol%natm, mol%cartesian, mol%atm, &
                            mol%bas_with_ecp, mol%env, mol%shell_offset, mol%necpbas, &
                            0, whole)
      call ecp_deriv_matrix(mol%nao, mol%nbas, mol%natm, mol%cartesian, mol%atm, &
                            mol%bas_with_ecp, mol%env, mol%shell_offset, mol%necpbas, &
                            1, one)
      call check(error, maxval(abs(one - whole)) < 1.0e-12_dp*maxval(abs(whole)), &
                 "iodine carries the only potential, and its share is not the whole")
      if (allocated(error)) return

      call ecp_deriv_matrix(mol%nao, mol%nbas, mol%natm, mol%cartesian, mol%atm, &
                            mol%bas_with_ecp, mol%env, mol%shell_offset, mol%necpbas, &
                            2, one)
      call check(error, all(one == 0.0_dp), &
                 "hydrogen has no potential and was given a derivative of one")
   end subroutine test_derivative_by_atom

end module test_mqc_ecp_matrix

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_ecp_matrix, only: collect_mqc_ecp_matrix_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_ecp_matrix", collect_mqc_ecp_matrix_tests) &
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
