!! RIJCOSX through the SCF: a fitted Coulomb matrix beside a chain-of-spheres exchange
module test_mqc_czt_cosx_scf
   !! `test_mqc_czt_cosx` settles that the exchange matrix is right. This
   !! settles that the SCF uses it: that it replaces the exchange and only the
   !! exchange, that the energy it converges to is the one the kernel's
   !! convergence promises, and that what it does not implement it refuses.
   !!
   !! **Two errors, kept apart.** An RIJCOSX energy differs from Hartree-Fock by
   !! the fit of the Coulomb matrix and by the quadrature of the exchange. The
   !! grid is varied with the fit held fixed, so the first group of checks sees
   !! the quadrature alone; the comparison with the unfitted energy then bounds
   !! the two together.
   !!
   !! **One identity is exact.** A pure functional has no exact exchange, so
   !! asking for chain-of-spheres must change nothing at all. That is the check
   !! that the Coulomb half is untouched, and it holds to rounding.
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available
   implicit none
   private

   public :: collect_mqc_czt_cosx_scf_tests

   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [3, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   character(len=*), parameter :: BASIS = "def2-svp"
   character(len=*), parameter :: AUX_BASIS = "def2-universal-jkfit"
   integer, parameter :: N_ELEC = 10

contains

   subroutine collect_mqc_czt_cosx_scf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("rijcosx_energy_converges_with_the_grid", test_grid), &
                  new_unittest("rijcosx_final_grid_energy_is_second_order", test_final_grid), &
                  new_unittest("rijcosx_is_close_to_hartree_fock", test_against_exact), &
                  new_unittest("rijcosx_changes_nothing_for_a_pure_functional", test_pure), &
                  new_unittest("rijcosx_scales_with_a_hybrid", test_hybrid), &
                  new_unittest("rijcosx_refuses_what_it_does_not_build", test_refusals) &
                  ]
   end subroutine collect_mqc_czt_cosx_scf_tests

   subroutine molecules(mol, aux, error)
      !! Water in the orbital basis and in the auxiliary one
      type(czt_molecule_t), intent(out) :: mol, aux
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, BASIS, mol, err)
      if (.not. err%has_error()) then
         call build_czt_molecule(WATER_Z, WATER_SYM, WATER, AUX_BASIS, aux, err)
      end if
      call check(error,.not. err%has_error(), "the molecules should build: "//err%get_full_trace())
   end subroutine molecules

   subroutine converge(mol, energy, error, aux, cosx, level, xc, final_grid)
      !! One closed-shell SCF, tightly converged, and its total energy
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(out) :: energy
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t), intent(in), optional :: aux
      logical, intent(in), optional :: cosx
      integer, intent(in), optional :: level
      type(xc_context_t), intent(inout), optional :: xc
      logical, intent(in), optional :: final_grid

      type(rhf_result_t) :: scf
      type(error_t) :: err

      energy = 0.0_dp
      call run_czt_rhf(mol, N_ELEC, 60, 1.0e-11_dp, 1.0e-9_dp, .false., scf, err, &
                       aux=aux, xc=xc, cosx=cosx, cosx_grid_level=level, &
                       cosx_final_grid=final_grid)
      call check(error,.not. err%has_error(), "the SCF should run: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, scf%converged, "the SCF should converge")
      if (allocated(error)) return
      energy = scf%energy
   end subroutine converge

   subroutine test_grid(error)
      !! With the fit fixed, the SCF energy settles as the exchange grid level rises
      !!
      !! The three levels, each with its final grid. Measured when this was
      !! written: level 1 is 8.6e-6 hartree from level 3 and level 2 is 2.9e-6.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      real(dp) :: e1, e2, e3

      call molecules(mol, aux, error)
      if (allocated(error)) return
      call converge(mol, e1, error, aux=aux, cosx=.true., level=1)
      if (allocated(error)) return
      call converge(mol, e2, error, aux=aux, cosx=.true., level=2)
      if (allocated(error)) return
      call converge(mol, e3, error, aux=aux, cosx=.true., level=3)
      if (allocated(error)) return

      call check(error, abs(e2 - e3) < abs(e1 - e3), &
                 "the SCF energy does not settle as the exchange grid level rises")
      if (allocated(error)) return
      call check(error, abs(e2 - e3) < 2.0e-5_dp, &
                 "the default exchange grid is too far from the largest")
      if (allocated(error)) return
      ! Guards the guard: a grid level that never reached the exchange would
      ! give three identical energies and pass both checks above.
      call check(error, abs(e1 - e3) > 1.0e-7_dp, &
                 "the grid level did not change the energy, so it was not used")
      call mol%destroy()
      call aux%destroy()
   end subroutine test_grid

   subroutine test_final_grid(error)
      !! The energy from a density converged on the small grid is second order from self-consistent
      !!
      !! With the final step on, the iterations run on a grid too coarse
      !! radially to take an energy from, and the energy is one build on the
      !! larger grid. With it off the iterations run on the larger grid
      !! throughout. The two differ only by how the density was converged,
      !! which the energy sees at second order: measured when this was
      !! written, 5.4e-8 hartree, where the small grid's own exchange
      !! energy is wrong by thousandths.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      real(dp) :: e_with, e_without

      call molecules(mol, aux, error)
      if (allocated(error)) return
      call converge(mol, e_with, error, aux=aux, cosx=.true., level=2)
      if (allocated(error)) return
      call converge(mol, e_without, error, aux=aux, cosx=.true., level=2, final_grid=.false.)
      if (allocated(error)) return

      call check(error, abs(e_with - e_without) < 2.0e-6_dp, &
                 "the energy from the small-grid density is too far from the "// &
                 "self-consistent one on the larger grid")
      if (allocated(error)) return
      ! Guards the guard: identical would mean the switch selected nothing.
      call check(error, abs(e_with - e_without) > 1.0e-12_dp, &
                 "switching the final grid off did not change how the SCF ran")
      call mol%destroy()
      call aux%destroy()
   end subroutine test_final_grid

   subroutine test_against_exact(error)
      !! RIJCOSX at the default grid, against Hartree-Fock with no approximation
      !!
      !! With a JK fitting set. Measured when this was written: the fitted
      !! energy is 5.6e-5 hartree above the exact one and RIJCOSX 2.6e-5 below
      !! it. The bound is an order of magnitude below what a missing or doubled
      !! exchange would leave, which is hartrees.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      real(dp) :: e_exact, e_fit, e_cosx

      call molecules(mol, aux, error)
      if (allocated(error)) return
      call converge(mol, e_exact, error)
      if (allocated(error)) return
      call converge(mol, e_fit, error, aux=aux)
      if (allocated(error)) return
      call converge(mol, e_cosx, error, aux=aux, cosx=.true.)
      if (allocated(error)) return

      call check(error, abs(e_cosx - e_exact) < 2.0e-4_dp, &
                 "RIJCOSX is too far from the unapproximated Hartree-Fock energy")
      if (allocated(error)) return
      call check(error, abs(e_cosx - e_fit) < 2.0e-4_dp, &
                 "RIJCOSX is too far from the fully fitted energy")
      if (allocated(error)) return
      ! Guards the guard: identical to the fitted energy would mean the
      ! exchange still came from the fit.
      call check(error, abs(e_cosx - e_fit) > 1.0e-7_dp, &
                 "RIJCOSX reproduced the fitted energy exactly, so the exchange "// &
                 "was not the seminumerical one")
      call mol%destroy()
      call aux%destroy()
   end subroutine test_against_exact

   subroutine test_pure(error)
      !! PBE has no exact exchange, so the request must be a no-op to rounding
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(xc_context_t) :: ctx
      type(error_t) :: err
      real(dp) :: e_fit, e_cosx

      if (.not. xc_available()) then
         call skip_test(error, "this build has no libxc")
         return
      end if
      call molecules(mol, aux, error)
      if (allocated(error)) return

      call xc_context_create(mol, "pbe", ctx, err, level=1)
      call check(error,.not. err%has_error(), "the functional should build: "//err%get_full_trace())
      if (allocated(error)) return
      call converge(mol, e_fit, error, aux=aux, xc=ctx)
      if (allocated(error)) return
      call converge(mol, e_cosx, error, aux=aux, cosx=.true., xc=ctx)
      if (allocated(error)) return

      call check(error, abs(e_cosx - e_fit) < 1.0e-10_dp, &
                 "chain-of-spheres changed the energy of a functional with no exact exchange")
      call ctx%destroy()
      call mol%destroy()
      call aux%destroy()
   end subroutine test_pure

   subroutine test_hybrid(error)
      !! PBE0 takes a quarter of the exchange, and the quadrature error scales with it
      !!
      !! Same grids as the Hartree-Fock case, and no final grid so each energy
      !! is on one quadrature. The difference between the coarse and the fine
      !! exchange grid should then be about a quarter of what it is there. A fraction applied twice, or not at all, is off by a factor
      !! of four in one direction or the other. Measured 0.29: the two
      !! differences are a few microhartree, so the ratio is not sharper than that.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(xc_context_t) :: ctx
      type(error_t) :: err
      real(dp) :: hf1, hf3, ks1, ks3, e_fit, ratio

      if (.not. xc_available()) then
         call skip_test(error, "this build has no libxc")
         return
      end if
      call molecules(mol, aux, error)
      if (allocated(error)) return

      call converge(mol, hf1, error, aux=aux, cosx=.true., level=1, final_grid=.false.)
      if (allocated(error)) return
      call converge(mol, hf3, error, aux=aux, cosx=.true., level=3, final_grid=.false.)
      if (allocated(error)) return

      call xc_context_create(mol, "pbe0", ctx, err, level=1)
      call check(error,.not. err%has_error(), "the functional should build: "//err%get_full_trace())
      if (allocated(error)) return
      call converge(mol, e_fit, error, aux=aux, xc=ctx)
      if (allocated(error)) return
      call converge(mol, ks1, error, aux=aux, cosx=.true., level=1, xc=ctx, final_grid=.false.)
      if (allocated(error)) return
      call converge(mol, ks3, error, aux=aux, cosx=.true., level=3, final_grid=.false., xc=ctx)
      if (allocated(error)) return

      call check(error, abs(ks3 - e_fit) < 1.0e-4_dp, &
                 "PBE0 with chain-of-spheres is too far from the fitted energy")
      if (allocated(error)) return
      ratio = (ks1 - ks3)/(hf1 - hf3)
      call check(error, ratio > 0.15_dp .and. ratio < 0.40_dp, &
                 "the quadrature error of a 25% hybrid is not about a quarter of "// &
                 "Hartree-Fock's, so the exchange fraction is misapplied")
      call ctx%destroy()
      call mol%destroy()
      call aux%destroy()
   end subroutine test_hybrid

   subroutine test_refusals(error)
      !! No fit, a second-order SCF and a range-separated functional are refused
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(xc_context_t) :: ctx
      type(rhf_result_t) :: scf
      type(error_t) :: err

      call molecules(mol, aux, error)
      if (allocated(error)) return

      call run_czt_rhf(mol, N_ELEC, 60, 1.0e-8_dp, 1.0e-6_dp, .false., scf, err, cosx=.true.)
      call check(error, err%has_error(), "chain-of-spheres without a fit was accepted")
      if (allocated(error)) return
      call err%clear()

      call run_czt_rhf(mol, N_ELEC, 60, 1.0e-8_dp, 1.0e-6_dp, .false., scf, err, &
                       aux=aux, cosx=.true., second_order=.true.)
      call check(error, err%has_error(), "chain-of-spheres with a second-order SCF was accepted")
      if (allocated(error)) return
      call err%clear()

      if (.not. xc_available()) return
      call xc_context_create(mol, "cam-b3lyp", ctx, err, level=1)
      call check(error,.not. err%has_error(), "the functional should build: "//err%get_full_trace())
      if (allocated(error)) return
      call run_czt_rhf(mol, N_ELEC, 60, 1.0e-8_dp, 1.0e-6_dp, .false., scf, err, &
                       aux=aux, cosx=.true., xc=ctx)
      call check(error, err%has_error(), &
                 "chain-of-spheres with a range-separated functional was accepted")
      call ctx%destroy()
      call mol%destroy()
      call aux%destroy()
   end subroutine test_refusals

end module test_mqc_czt_cosx_scf

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_cosx_scf, only: collect_mqc_czt_cosx_scf_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_cosx_scf", collect_mqc_czt_cosx_scf_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
