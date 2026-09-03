!! RI-MP2 over a density-fitted reference, against its own energy
program check_df_ref_gradient
   !! **What this differentiates, and why it is a third thing.** An `ri-mp2`
   !! deck fits the correlation and leaves the SCF exact;
   !! `keywords.scf.density_fitting` on its own fits the SCF and leaves the
   !! correlation alone. Both together is neither, and it has its own gradient:
   !! the response operator, the two potentials built from the relaxed density
   !! and the whole two-electron derivative term all move to the fitted
   !! integrals, and the last of those stops being a four-centre contraction at
   !! all.
   !!
   !! **The reference is our own energy, extrapolated.** PySCF does implement a
   !! DF-MP2 gradient over a DF-HF reference, and comparing against it is worth
   !! doing -- but not first. It would conflate "is this the derivative of the
   !! energy we computed" with every question about whether two programs fit the
   !! same thing, and only the first of those has a clean answer. So this
   !! differentiates the energy this program produces, and a plain central
   !! difference is not enough to do it: at 2.5e-3 the truncation error is 5e-6,
   !! which is the size of a wrong factor in one of the two exchange `Gamma`
   !! halves. Two steps combined as `(4 D(h/2) - D(h))/3` remove the `h^2` term
   !! and leave a check that can tell those apart.
   !!
   !! The pieces underneath are validated separately and beforehand:
   !! `check_fitted_reference_gradient` takes the derivative of the fitted
   !! two-electron operator against a finite difference with two indefinite
   !! matrices, and `check_df_cphf` takes the fitted response operator. What is
   !! left for this program is the assembly.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mp2, only: mp2_result_t, run_czt_ri_mp2
   use mqc_czt_ri_mp2_gradient, only: czt_ri_mp2_gradient
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none

   integer, parameter :: N_DIM = 3
   real(dp), parameter :: STEP = 2.5e-3_dp
      !! The coarse step. The fine one is half of it, and the two are combined.
   real(dp), parameter :: TOL = 1.0e-8_dp
   real(dp), parameter :: DENSITY_TOL = 1.0e-12_dp
      !! There is no quadrature anywhere in this, unlike the double hybrid next
      !! door, so the SCF can be driven as far as the arithmetic allows and the
      !! finite difference is limited by its own truncation rather than by noise.
   integer :: n_bad

   n_bad = 0
   write (*, "(a)") "== RI-MP2 over a density-fitted reference"

   call check_case("H2 / sto-3g", [1, 1], ["H", "H"], &
                   reshape([0.0_dp, 0.0_dp, -0.7_dp, &
                            0.0_dp, 0.0_dp, 0.7_dp], [N_DIM, 2]), &
                   "sto-3g", "cc-pvdz-rifit", 2, n_bad)

   call check_case("H2O / sto-3g", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", "cc-pvdz-rifit", 10, n_bad)

   ! Asymmetric, with p functions on two centres: the exchange term carries one
   ! density on each side of the metric-weighted integrals, and swapping them
   ! survives a symmetric molecule.
   call check_case("HCN / sto-3g", [1, 6, 7], ["H", "C", "N"], &
                   reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                   "sto-3g", "cc-pvdz-rifit", 14, n_bad)

   write (*, "(a)") ""
   if (n_bad == 0) then
      write (*, "(a)") "all density-fitted reference gradient checks passed"
   else
      write (*, "(a,i0,a)") "FAILED: ", n_bad, " case(s)"
      stop 1
   end if

contains

   subroutine check_case(label, numbers, symbols, coords, basis, aux_basis, nelec, n_bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: nelec
      integer, intent(inout) :: n_bad

      type(error_t) :: error
      real(dp), allocatable :: analytic(:, :), numeric(:, :)
      real(dp) :: coarse, fine, e0, worst
      integer :: natm, ia, comp

      natm = size(numbers)
      write (*, "(a)") ""
      write (*, "(a,a)") "== ", label
      flush (output_unit)

      call analytic_gradient(numbers, symbols, coords, basis, aux_basis, nelec, &
                             analytic, error)
      if (fail(error, n_bad)) return
      call energy_at(numbers, symbols, coords, basis, aux_basis, nelec, e0, error)
      if (fail(error, n_bad)) return

      allocate (numeric(3, natm))
      do ia = 1, natm
         do comp = 1, 3
            call central(numbers, symbols, coords, basis, aux_basis, nelec, &
                         ia, comp, STEP, coarse, error)
            if (fail(error, n_bad)) return
            call central(numbers, symbols, coords, basis, aux_basis, nelec, &
                         ia, comp, 0.5_dp*STEP, fine, error)
            if (fail(error, n_bad)) return
            ! Richardson: the h^2 error cancels between the two, leaving h^4.
            numeric(comp, ia) = (4.0_dp*fine - coarse)/3.0_dp
         end do
      end do

      write (*, "(a,f20.12)") "   E(RI-MP2 / DF-RHF) = ", e0
      write (*, "(a)") "   atom      analytic (x, y, z)                        "// &
         "finite difference"
      do ia = 1, natm
         write (*, "(i6,3f14.9,a,3f14.9)") ia, analytic(:, ia), "   ", numeric(:, ia)
      end do
      worst = maxval(abs(analytic - numeric))
      write (*, "(a,es14.4)") "   largest deviation: ", worst
      write (*, "(a,es14.4)") "   |sum over atoms|:  ", &
         maxval(abs(sum(analytic, dim=2)))
      flush (output_unit)

      if (worst > TOL) then
         write (*, "(a)") "  FAIL: analytic and finite difference disagree"
         n_bad = n_bad + 1
      end if
      ! Blind to most of what can go wrong here -- it held at 1e-14 through a
      ! version of this gradient that was out by a factor of thirty -- but a
      ! term landing on the wrong atom fails it and nothing else does.
      if (maxval(abs(sum(analytic, dim=2))) > 1.0e-9_dp) then
         write (*, "(a)") "  FAIL: does not sum to zero over atoms"
         n_bad = n_bad + 1
      end if
   end subroutine check_case

   subroutine analytic_gradient(numbers, symbols, coords, basis, aux_basis, nelec, &
                                gradient, error)
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: nelec
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf

      call converge(numbers, symbols, coords, basis, aux_basis, nelec, mol, aux, &
                    scf, error)
      if (error%has_error()) return
      call czt_ri_mp2_gradient(mol, aux, scf%orbitals, scf%orbital_energies, &
                               nelec/2, gradient, error, fitted_reference=.true.)
      call mol%destroy()
      call aux%destroy()
   end subroutine analytic_gradient

   subroutine central(numbers, symbols, coords, basis, aux_basis, nelec, &
                      ia, comp, h, deriv, error)
      !! One central difference of the same energy the gradient claims
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: nelec, ia, comp
      real(dp), intent(in) :: h
      real(dp), intent(out) :: deriv
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: shifted(:, :)
      real(dp) :: plus, minus

      deriv = 0.0_dp
      allocate (shifted(size(coords, 1), size(coords, 2)))
      shifted = coords
      shifted(comp, ia) = coords(comp, ia) + h
      call energy_at(numbers, symbols, shifted, basis, aux_basis, nelec, plus, error)
      if (error%has_error()) return
      shifted(comp, ia) = coords(comp, ia) - h
      call energy_at(numbers, symbols, shifted, basis, aux_basis, nelec, minus, error)
      if (error%has_error()) return
      deriv = (plus - minus)/(2.0_dp*h)
   end subroutine central

   subroutine energy_at(numbers, symbols, coords, basis, aux_basis, nelec, energy, error)
      !! The density-fitted SCF plus the fitted correlation, at one geometry
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: nelec
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(mp2_result_t) :: mp2

      energy = 0.0_dp
      call converge(numbers, symbols, coords, basis, aux_basis, nelec, mol, aux, &
                    scf, error)
      if (error%has_error()) return
      call run_czt_ri_mp2(mol, aux, scf%orbitals, scf%orbital_energies, nelec/2, &
                          scf%energy, mp2, error)
      if (error%has_error()) return
      energy = mp2%total
      call mol%destroy()
      call aux%destroy()
   end subroutine energy_at

   subroutine converge(numbers, symbols, coords, basis, aux_basis, nelec, mol, aux, &
                       scf, error)
      !! The density-fitted SCF, with one auxiliary basis serving both halves
      !!
      !! `aux` is handed to `run_czt_rhf`, which is what makes the reference
      !! fitted; the same molecule then names the correlation fit. That is not a
      !! simplification for the test's sake -- `model.aux_basis` is the only
      !! place a deck can name a fitting set, so it is the only combination that
      !! can be asked for.
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: nelec
      type(czt_molecule_t), intent(out) :: mol, aux
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: error

      call build_czt_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) return
      call build_czt_molecule(numbers, symbols, coords, aux_basis, aux, error)
      if (error%has_error()) return
      call run_czt_rhf(mol, nelec, 300, 1.0e-14_dp, DENSITY_TOL, .false., scf, &
                       error, aux=aux)
   end subroutine converge

   function fail(error, n_bad) result(bad)
      type(error_t), intent(inout) :: error
      integer, intent(inout) :: n_bad
      logical :: bad

      bad = error%has_error()
      if (bad) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
      end if
   end function fail

end program check_df_ref_gradient
