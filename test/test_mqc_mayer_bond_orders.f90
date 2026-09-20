!! Mayer bond orders, against numbers computed somewhere else
module test_mqc_mayer_bond_orders
   !! Unlike a Mulliken charge, a Mayer bond order has a reference: it is a
   !! closed-form function of the converged density and the overlap, so another
   !! program handed the same wave function must produce the same matrix. Every
   !! number pinned here came out of PySCF and numpy, through
   !! `tools/cpu_validation/mayer_bond_orders_reference.py`, which feeds PySCF
   !! *this repository's* basis JSON -- PySCF's internal Pople tables differ in
   !! the eighth decimal of the exponents and fake a disagreement that looks
   !! exactly like a bug here.
   !!
   !! **The open-shell case is the one that earns its place.** The open-shell
   !! formula is not the closed-shell expression applied to the total density,
   !! and the difference vanishes on a closed shell -- so a wrong
   !! implementation passes every restricted test there is. Triplet O2 is where
   !! it shows: the open-shell formula gives exactly 2.00 and the closed-shell
   !! one gives 1.50 on the same density. `open_shell_formula_is_not_the_closed`
   !! computes both and refuses to let them be confused.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: run_czt_rhf, run_czt_uhf, rhf_result_t
   use mqc_czt_charges, only: ao_to_atom
   use mqc_czt_bond_orders, only: mayer_bond_orders, mayer_bond_orders_open_shell, &
                                  mayer_valences
   use mqc_population_analysis, only: mayer_atomic_bond_orders, &
                                      mayer_atomic_bond_orders_open_shell, &
                                      mayer_atomic_valences
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_mayer_bond_orders_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

   ! How far our matrix may sit from PySCF's. This is not a convergence
   ! tolerance: both codes solve the same equations in the same basis, so the
   ! gap is whatever the two SCFs disagree by, and 1e-7 is loose against the
   ! 1e-10 density threshold used below while being tight enough that no wrong
   ! formula fits inside it -- the open-shell error is five parts in ten.
   real(dp), parameter :: TOL = 1.0e-7_dp

contains

   subroutine collect_mqc_mayer_bond_orders_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("water_matches_pyscf", test_water), &
                  new_unittest("ethane_matches_pyscf_and_chemistry", test_ethane), &
                  new_unittest("open_shell_matches_pyscf", test_triplet_oxygen), &
                  new_unittest("open_shell_formula_is_not_the_closed", test_formulas_differ), &
                  new_unittest("equal_spin_densities_give_the_closed_shell", test_equal_spins), &
                  new_unittest("doublet_matches_pyscf", test_doublet_cation), &
                  new_unittest("matrix_is_symmetric_with_no_diagonal", test_shape), &
                  new_unittest("valence_is_the_row_sum", test_valence), &
                  new_unittest("an_ao_reordering_changes_nothing", test_reordering), &
                  new_unittest("mismatched_shapes_are_refused", test_refusal) &
                  ]
   end subroutine collect_mqc_mayer_bond_orders_tests

   subroutine water(mol, err, basis)
      !! The geometry `test_mqc_czt_charges.f90` uses, and the reference script
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      character(len=*), intent(in) :: basis

      real(dp) :: c(3, 3)
      integer :: z(3)
      character(len=2) :: symbols(3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule(z, symbols, c, basis, mol, err)
   end subroutine water

   subroutine ethane(mol, err)
      !! Staggered ethane, C-C 1.526 A and C-H 1.088 A at 111.2 degrees
      !!
      !! The coordinates are the reference script's, printed to six decimals,
      !! so the two cannot drift apart.
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err

      real(dp) :: c(3, 8)
      integer :: z(8)
      character(len=2) :: symbols(8)

      z = [6, 6, 1, 1, 1, 1, 1, 1]
      symbols = ["C ", "C ", "H ", "H ", "H ", "H ", "H ", "H "]
      c = reshape([0.000000_dp, 0.000000_dp, 0.000000_dp, &
                   0.000000_dp, 0.000000_dp, 1.526000_dp, &
                   1.014368_dp, 0.000000_dp, -0.393448_dp, &
                   -0.507184_dp, 0.878469_dp, -0.393448_dp, &
                   -0.507184_dp, -0.878469_dp, -0.393448_dp, &
                   0.507184_dp, 0.878469_dp, 1.919448_dp, &
                   -1.014368_dp, 0.000000_dp, 1.919448_dp, &
                   0.507184_dp, -0.878469_dp, 1.919448_dp], [3, 8])*ANG
      call build_czt_molecule(z, symbols, c, "sto-3g", mol, err)
   end subroutine ethane

   subroutine oxygen(mol, err)
      !! O2 at 1.208 A, whose ground state is a triplet
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err

      real(dp) :: c(3, 2)
      integer :: z(2)
      character(len=2) :: symbols(2)

      z = [8, 8]
      symbols = ["O ", "O "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 1.208_dp*ANG], [3, 2])
      call build_czt_molecule(z, symbols, c, "sto-3g", mol, err)
   end subroutine oxygen

   subroutine closed_shell_orders(mol, nelec, orders, err)
      !! Converge a restricted SCF and take its Mayer orders
      type(czt_molecule_t), intent(in) :: mol
      integer, intent(in) :: nelec
      real(dp), allocatable, intent(out) :: orders(:, :)
      type(error_t), intent(inout) :: err

      type(rhf_result_t) :: scf
      real(dp), allocatable :: overlap(:, :)

      ! The commutator threshold is set rather than derived: a consumer of the
      ! density needs it, because the density error goes as the gradient norm
      ! and not as its square.
      call run_czt_rhf(mol, nelec, 80, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, &
                       grad_tol=1.0e-8_dp)
      if (err%has_error()) return
      call mol%overlap(overlap)
      call mayer_bond_orders(mol, scf%density, overlap, orders, err)
   end subroutine closed_shell_orders

   subroutine test_water(error)
      !! Water, STO-3G, against PySCF
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: orders(:, :)

      call water(mol, err, "sto-3g")
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return
      call closed_shell_orders(mol, 10, orders, err)
      call check(error,.not. err%has_error(), "the reference SCF and the bond orders")
      if (allocated(error)) return

      call check(error, abs(orders(1, 2) - 0.954168619763_dp) < TOL, "O-H")
      if (allocated(error)) return
      call check(error, abs(orders(1, 3) - 0.954014119750_dp) < TOL, "the other O-H")
      if (allocated(error)) return
      ! Two hydrogens 1.5 A apart with no bond between them, and a number that
      ! is small rather than zero: Mayer's definition has no distance in it and
      ! says so honestly.
      call check(error, abs(orders(2, 3) - 0.012473116986_dp) < TOL, "H...H")
   end subroutine test_water

   subroutine test_ethane(error)
      !! Ethane: one C-C, six C-H, and nothing across the molecule
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: orders(:, :), valences(:)
      integer :: iatom

      call ethane(mol, err)
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return
      call closed_shell_orders(mol, 18, orders, err)
      call check(error,.not. err%has_error(), "the reference SCF and the bond orders")
      if (allocated(error)) return

      call check(error, abs(orders(1, 2) - 1.010966289912_dp) < TOL, "C-C")
      if (allocated(error)) return
      do iatom = 3, 5
         call check(error, abs(orders(1, iatom) - 0.983590579729_dp) < TOL, &
                    "C1-H")
         if (allocated(error)) return
         ! The hydrogens on the *other* carbon. Near zero, which is the whole
         ! claim: a bond order that fell out of a distance criterion would put
         ! something here, and one that fell out of the density does not.
         call check(error, abs(orders(1, iatom + 3) - 0.002048125249_dp) < TOL, &
                    "C1 to a hydrogen on C2")
         if (allocated(error)) return
      end do

      call mayer_valences(orders, valences)
      call check(error, abs(valences(1) - 3.967882404846_dp) < TOL, &
                 "a carbon in ethane is tetravalent")
      if (allocated(error)) return
      call check(error, abs(valences(3) - 0.996635494930_dp) < TOL, &
                 "a hydrogen in ethane is monovalent")
   end subroutine test_ethane

   subroutine triplet_oxygen_densities(mol, scf, err)
      !! Converge the triplet, which is O2's ground state
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      call oxygen(mol, err)
      if (err%has_error()) return
      call run_czt_uhf(mol, 16, 3, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, &
                       grad_tol=1.0e-8_dp)
   end subroutine triplet_oxygen_densities

   subroutine test_triplet_oxygen(error)
      !! Triplet O2: a double bond, and PySCF says 2.000000 exactly
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :), orders(:, :)

      call triplet_oxygen_densities(mol, scf, err)
      call check(error,.not. err%has_error(), "the reference UHF")
      if (allocated(error)) return
      call check(error, scf%converged, "the reference UHF did not converge")
      if (allocated(error)) return

      call mol%overlap(overlap)
      call mayer_bond_orders_open_shell(mol, scf%density, scf%density_beta, overlap, &
                                        orders, err)
      call check(error,.not. err%has_error(), "the open-shell bond orders")
      if (allocated(error)) return
      call check(error, abs(orders(1, 2) - 2.000000000000_dp) < TOL, &
                 "triplet O2 is a double bond")
   end subroutine test_triplet_oxygen

   subroutine test_formulas_differ(error)
      !! The wrong formula, computed on purpose, and shown to be wrong
      !!
      !! This is the test that makes the open-shell entry worth having. It
      !! feeds the *total* density to the closed-shell routine -- exactly the
      !! mistake being guarded against -- and checks that what comes back is
      !! far from the reference. If someone ever implements the open-shell
      !! case as the closed-shell expression, this fails loudly where the
      !! restricted tests above would all still pass.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :), total(:, :)
      real(dp), allocatable :: orders(:, :), naive(:, :)

      call triplet_oxygen_densities(mol, scf, err)
      call check(error,.not. err%has_error(), "the reference UHF")
      if (allocated(error)) return

      call mol%overlap(overlap)
      call mayer_bond_orders_open_shell(mol, scf%density, scf%density_beta, overlap, &
                                        orders, err)
      call check(error,.not. err%has_error(), "the open-shell bond orders")
      if (allocated(error)) return

      total = scf%density + scf%density_beta
      call mayer_bond_orders(mol, total, overlap, naive, err)
      call check(error,.not. err%has_error(), "the closed-shell expression")
      if (allocated(error)) return

      ! PySCF, same density, closed-shell expression: 1.498294548252.
      call check(error, abs(naive(1, 2) - 1.498294548252_dp) < TOL, &
                 "the wrong formula does not give the number it is known to give")
      if (allocated(error)) return
      call check(error, abs(orders(1, 2) - naive(1, 2)) > 0.5_dp, &
                 "the two formulas must disagree on an open shell; if they do "// &
                 "not, one of them is not what it claims to be")
   end subroutine test_formulas_differ

   subroutine test_equal_spins(error)
      !! Da = Db = D/2 sends the open-shell formula back to the closed-shell one
      !!
      !! Machine precision, not a tolerance: it is the same arithmetic
      !! rearranged, over the same converged density, so anything above
      !! rounding is an algebra error in one of the two entries.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :), half(:, :)
      real(dp), allocatable :: closed(:, :), open_form(:, :)

      call water(mol, err, "sto-3g")
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return
      call run_czt_rhf(mol, 10, 80, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, &
                       grad_tol=1.0e-8_dp)
      call check(error,.not. err%has_error(), "the reference SCF")
      if (allocated(error)) return

      call mol%overlap(overlap)
      call mayer_bond_orders(mol, scf%density, overlap, closed, err)
      half = 0.5_dp*scf%density
      call mayer_bond_orders_open_shell(mol, half, half, overlap, open_form, err)
      call check(error,.not. err%has_error(), "both bond order entries")
      if (allocated(error)) return
      call check(error, maxval(abs(closed - open_form)) < 1.0e-12_dp, &
                 "the open-shell formula must reduce to the closed-shell one")
   end subroutine test_equal_spins

   subroutine test_doublet_cation(error)
      !! The water cation, where the two formulas differ by only 9e-4
      !!
      !! The opposite end of the range from triplet O2, and the reason the
      !! tolerance here is 1e-7 rather than something comfortable: a doublet
      !! with one unpaired electron out of nine is where a wrong formula is
      !! nearly right, and a loose test would let it through.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :), orders(:, :)

      call water(mol, err, "sto-3g")
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return
      call run_czt_uhf(mol, 9, 2, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, &
                       grad_tol=1.0e-8_dp)
      call check(error,.not. err%has_error(), "the reference UHF")
      if (allocated(error)) return

      call mol%overlap(overlap)
      call mayer_bond_orders_open_shell(mol, scf%density, scf%density_beta, overlap, &
                                        orders, err)
      call check(error,.not. err%has_error(), "the open-shell bond orders")
      if (allocated(error)) return
      call check(error, abs(orders(1, 2) - 0.782819137661_dp) < TOL, "O-H")
      if (allocated(error)) return
      call check(error, abs(orders(1, 3) - 0.782810453965_dp) < TOL, "the other O-H")
   end subroutine test_doublet_cation

   subroutine test_shape(error)
      !! Symmetric, and nothing on the diagonal
      !!
      !! The diagonal is not a bond and is deliberately zeroed. Leaving what
      !! the block sum puts there -- roughly twice the atom's own population --
      !! would silently double every valence.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: orders(:, :)
      integer :: iatom, jatom

      call water(mol, err, "6-31g")
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return
      call closed_shell_orders(mol, 10, orders, err)
      call check(error,.not. err%has_error(), "the reference SCF and the bond orders")
      if (allocated(error)) return

      call check(error, size(orders, 1) == mol%natm .and. size(orders, 2) == mol%natm, &
                 "one entry per atom pair")
      if (allocated(error)) return
      do iatom = 1, mol%natm
         call check(error, abs(orders(iatom, iatom)) < 1.0e-14_dp, &
                    "an atom is not bonded to itself")
         if (allocated(error)) return
         do jatom = 1, mol%natm
            call check(error, abs(orders(iatom, jatom) - orders(jatom, iatom)) < 1.0e-12_dp, &
                       "the matrix must be symmetric")
            if (allocated(error)) return
         end do
      end do

      ! And the 6-31G numbers, since this molecule is converged anyway: the
      ! same bond in a bigger basis is 0.80 rather than 0.95, which is the
      ! basis-set dependence worth having on record.
      call check(error, abs(orders(1, 2) - 0.803341305998_dp) < TOL, "O-H in 6-31G")
   end subroutine test_shape

   subroutine test_valence(error)
      !! The valence is the row sum, on a matrix nobody had to converge
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: orders(3, 3)
      real(dp), allocatable :: valences(:)

      orders = reshape([0.0_dp, 1.0_dp, 0.25_dp, &
                        1.0_dp, 0.0_dp, 0.5_dp, &
                        0.25_dp, 0.5_dp, 0.0_dp], [3, 3])
      call mayer_atomic_valences(orders, valences)
      call check(error, size(valences) == 3, "one valence per atom")
      if (allocated(error)) return
      call check(error, abs(valences(1) - 1.25_dp) < 1.0e-14_dp, "first row")
      if (allocated(error)) return
      call check(error, abs(valences(2) - 1.5_dp) < 1.0e-14_dp, "second row")
      if (allocated(error)) return
      call check(error, abs(valences(3) - 0.75_dp) < 1.0e-14_dp, "third row")
   end subroutine test_valence

   subroutine test_reordering(error)
      !! Permuting the basis functions permutes nothing about the answer
      !!
      !! The bond order is a sum over the AO indices of an atom, so it cannot
      !! depend on the order they are stored in -- as long as the owner array
      !! is permuted with them. Reversing the whole AO index is the harshest
      !! permutation available and costs one more pass over a 7x7 matrix.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: overlap(:, :), orders(:, :), shuffled(:, :)
      real(dp), allocatable :: d_perm(:, :), s_perm(:, :)
      integer, allocatable :: owner(:), owner_perm(:), perm(:)
      integer :: nao, mu, nu

      call water(mol, err, "sto-3g")
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return
      call run_czt_rhf(mol, 10, 80, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, &
                       grad_tol=1.0e-8_dp)
      call check(error,.not. err%has_error(), "the reference SCF")
      if (allocated(error)) return
      call mol%overlap(overlap)
      call ao_to_atom(mol, owner)
      nao = size(owner)

      call mayer_atomic_bond_orders(owner, mol%natm, scf%density, overlap, orders, err)
      call check(error,.not. err%has_error(), "the bond orders in the natural order")
      if (allocated(error)) return

      allocate (perm(nao), owner_perm(nao), d_perm(nao, nao), s_perm(nao, nao))
      do mu = 1, nao
         perm(mu) = nao - mu + 1
      end do
      do mu = 1, nao
         owner_perm(mu) = owner(perm(mu))
         do nu = 1, nao
            d_perm(mu, nu) = scf%density(perm(mu), perm(nu))
            s_perm(mu, nu) = overlap(perm(mu), perm(nu))
         end do
      end do

      call mayer_atomic_bond_orders(owner_perm, mol%natm, d_perm, s_perm, shuffled, err)
      call check(error,.not. err%has_error(), "the bond orders in the shuffled order")
      if (allocated(error)) return
      call check(error, maxval(abs(orders - shuffled)) < 1.0e-13_dp, &
                 "the matrix must not depend on the basis function order")
   end subroutine test_reordering

   subroutine test_refusal(error)
      !! A density that is not the size of the basis is refused, not indexed
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp) :: density(4, 4), overlap(4, 4)
      real(dp), allocatable :: orders(:, :)
      integer :: owner(5)

      density = 0.0_dp
      overlap = 0.0_dp
      owner = [1, 1, 2, 2, 2]

      call err%clear()
      call mayer_atomic_bond_orders(owner, 2, density, overlap, orders, err)
      call check(error, err%has_error(), &
                 "a density smaller than the owner array must be refused")
      if (allocated(error)) return
      call check(error,.not. allocated(orders), &
                 "a refused call must not leave a matrix behind")
      if (allocated(error)) return

      call err%clear()
      call mayer_atomic_bond_orders_open_shell(owner, 2, density, density, overlap, &
                                               orders, err)
      call check(error, err%has_error(), &
                 "the open-shell entry must refuse the same mismatch")
   end subroutine test_refusal

end module test_mqc_mayer_bond_orders

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mayer_bond_orders, only: collect_mqc_mayer_bond_orders_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_mayer_bond_orders", &
                               collect_mqc_mayer_bond_orders_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
