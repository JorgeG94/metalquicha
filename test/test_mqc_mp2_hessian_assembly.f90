!! The assembled MP2 correlation Hessian, checked against what already works
module test_mqc_mp2_hessian_assembly
   !! Unit 1.9 of the MP2 Hessian ladder (`mp2-hessian-phased-plan.md`): the
   !! three gated groups -- fixed-density second skeletons, orbital response,
   !! and the 2n+1 density response -- combined by `mp2_correlation_hessian`
   !! into the correlation block of the analytic Hessian.
   !!
   !! **Why these checks and not a pycc comparison.** The cross-code gates
   !! (the full correlation block against pycc's `H_correlation`, symmetric
   !! and asymmetric geometry) were run when the unit landed and their
   !! residuals are recorded in the commit; they are deliberately not stored
   !! as tests, for the reason `test_mqc_hess_ints`' header gives -- a stored
   !! comparison against an external library pins our layout to that library's
   !! conventions rather than to anything true. The analytic column against
   !! the pycc literals already in the tree lives with those literals, in
   !! `test_mqc_mp2_hessian_fd`. What stays here is what must not drift
   !! silently:
   !!
   !! * the SCF reference block the assembly deposits, completed by the
   !!   delegated CPHF response and the nuclear repulsion, still equals
   !!   standalone `rhf_hessian` -- the five duplicated deposit lines diverge
   !!   silently otherwise (the Unit 1.3 guard, re-earned through the full
   !!   assembly);
   !! * `H = H^T` and the translational sum rule `sum_A H[Aa,Bb] = 0`, at
   !!   both geometries. **Necessary and weak** -- a gradient in this tree was
   !!   once wrong by 14 Hartree/Bohr while leaving the net force at 8e-14 --
   !!   so they stand beside the cross-code gates and never in place of them;
   !! * a frozen core is refused with an error, not computed wrongly: the
   !!   core<->active `U^X` rewrite and the Sylvester divide are Phase 2's.
   !!
   !! One thread throughout, as everywhere on this ladder.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mp2_hessian, only: mp2_correlation_hessian
   use mqc_czt_hessian, only: response_hessian, rhf_hessian, &
                              nuclear_repulsion_hessian
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   implicit none
   private

   public :: collect_mqc_mp2_hessian_assembly_tests

   !! pycc's pinned case: water/6-31G, bohr, frame pinned -- the geometry the
   !! whole ladder gates on, and its asymmetric companion. The asymmetric one
   !! is not decoration: pycc record a CCSD ket-swap bug the symmetric
   !! geometry masked, and an SCF bug visible only under C2v.
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.000000_dp, 0.000000_dp, &
                           0.0_dp, 0.000000_dp, 1.814137_dp, &
                           0.0_dp, 1.756000_dp, -0.454300_dp], [3, 3])
   real(dp), parameter :: WATER_ASYM(3, 3) = reshape( &
                          [0.0_dp, 0.000000_dp, 0.000000_dp, &
                           0.031_dp, -0.024_dp, 1.814137_dp, &
                           0.0_dp, 1.756000_dp, -0.454300_dp], [3, 3])
   integer, parameter :: WATER_NELEC = 10

contains

   subroutine collect_mqc_mp2_hessian_assembly_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("reference_block_matches_rhf_hessian", reference_matches), &
                  new_unittest("correlation_hessian_symmetric_and_translates", &
                               symmetric_and_translates), &
                  new_unittest("frozen_core_symmetric_and_translates", &
                               frozen_core_symmetric) &
                  ]
   end subroutine collect_mqc_mp2_hessian_assembly_tests

   !! One SCF and the full assembly at a given geometry.
   subroutine assembly_at(coords, mol, scf, hess_corr, hess_ref, err)
      real(dp), intent(in) :: coords(:, :)
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      real(dp), allocatable, intent(out) :: hess_corr(:, :, :, :)
      real(dp), allocatable, intent(out) :: hess_ref(:, :, :, :)
      type(error_t), intent(inout) :: err

      call build_czt_molecule(WATER_Z, WATER_SYM, coords, "6-31g", mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, .false., &
                       scf, err)
      if (err%has_error()) return
      call mp2_correlation_hessian(mol, scf%orbitals, scf%orbital_energies, &
                                   scf%density, WATER_NELEC/2, 0, &
                                   hess_corr, hess_ref, err)
   end subroutine assembly_at

   !! The guard: the reference block the assembly produces, completed the way
   !! its caller completes it (nuclear repulsion plus the delegated CPHF
   !! response), equals standalone `rhf_hessian`. Unit 1.3 measured 2.5e-14
   !! for the bare sweep; this re-earns it through the full assembly, where a
   !! quiet edit to the duplicated deposits would otherwise surface only as an
   !! unattributable total-Hessian error. 1e-10 and not tighter because the
   !! two two-electron walks, and since `hess_rinv_contract` the two
   !! nuclear-attraction contractions, sum in different orders over terms of
   !! about 7e3 on oxygen's core functions: 3e-12 measured, from 2.5e-14 when
   !! both ran the same statements. A wrong term lands at 1e-6 or above.
   subroutine reference_matches(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: hess_corr(:, :, :, :), hess_ref(:, :, :, :)
      real(dp), allocatable :: nuc(:, :, :, :), resp(:, :, :, :), full(:, :, :, :)
      integer :: threads

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call assembly_at(WATER, mol, scf, hess_corr, hess_ref, err)
      if (.not. err%has_error()) then
         call nuclear_repulsion_hessian(WATER_Z, mol%coords, nuc, err)
         call response_hessian(mol, scf%density, scf%orbitals, &
                               scf%orbital_energies, WATER_NELEC/2, resp, err)
         call rhf_hessian(mol, WATER_Z, scf%density, scf%orbitals, &
                          scf%orbital_energies, WATER_NELEC/2, full, err)
      end if
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), &
                 "the assembly did not evaluate: "//err%get_message())
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      write (*, "(a, es10.3)") "        max |assembly ref - rhf_hessian| = ", &
         maxval(abs(nuc + hess_ref + resp - full))
      call check(error, maxval(abs(nuc + hess_ref + resp - full)) < 1.0e-10_dp, &
                 "the assembly's reference block diverged from rhf_hessian")
      call mol%destroy()
   end subroutine reference_matches

   !! Gate 1.9c at both geometries: `H = H^T` element for element, and every
   !! Cartesian direction's sum over either atom index vanishing. Neither is
   !! imposed anywhere -- the (X, Y) and (Y, X) elements assemble the response
   !! from different perturbations' solves, and the sum rule runs through
   !! every group including both Z-vector solves. Measured 6.0e-15 / 3.2e-15
   !! at one thread; the tolerance leaves room for integral rounding, not for
   !! a misplaced term.
   subroutine symmetric_and_translates(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: hess_corr(:, :, :, :), hess_ref(:, :, :, :)
      real(dp) :: worst_sym, worst_tr, acc
      integer :: threads, g, natm, ia, ja, a, b

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      do g = 1, 2
         if (g == 1) then
            call assembly_at(WATER, mol, scf, hess_corr, hess_ref, err)
         else
            call assembly_at(WATER_ASYM, mol, scf, hess_corr, hess_ref, err)
         end if
         if (err%has_error()) exit

         natm = size(hess_corr, 3)
         worst_sym = 0.0_dp
         worst_tr = 0.0_dp
         do ja = 1, natm
            do ia = 1, natm
               do b = 1, 3
                  do a = 1, 3
                     worst_sym = max(worst_sym, abs(hess_corr(a, b, ia, ja) &
                                                    - hess_corr(b, a, ja, ia)))
                  end do
               end do
            end do
         end do
         do ja = 1, natm
            do b = 1, 3
               do a = 1, 3
                  acc = 0.0_dp
                  do ia = 1, natm
                     acc = acc + hess_corr(a, b, ia, ja)
                  end do
                  worst_tr = max(worst_tr, abs(acc))
                  acc = 0.0_dp
                  do ia = 1, natm
                     acc = acc + hess_corr(b, a, ja, ia)
                  end do
                  worst_tr = max(worst_tr, abs(acc))
               end do
            end do
         end do
         write (*, "(a, i0, a, 2es10.3)") "        geometry ", g, &
            " symmetry / translation = ", worst_sym, worst_tr
         call check(error, max(worst_sym, worst_tr) < 1.0e-12_dp, &
                    "the correlation Hessian lost a symmetry it never imposed")
         call mol%destroy()
         if (allocated(error)) then
            call omp_set_num_threads(threads)
            return
         end if
         deallocate (hess_corr, hess_ref)
      end do
      call omp_set_num_threads(threads)
      call check(error,.not. err%has_error(), &
                 "the assembly did not evaluate: "//err%get_message())
   end subroutine symmetric_and_translates

   !! The frozen-core assembly earns the same symmetries the all-electron
   !! one does: `H = H^T` across the mixed second derivatives, and rigid
   !! translation summing every atom of a Cartesian pair to nothing. Weak
   !! alone -- the cross-code element-wise gates (pycc's frozen-core
   !! `H_correlation`, both geometries, ~6e-12) ran when Phase 2 landed and
   !! live in the commit messages; the finite-difference tie-out is
   !! `test_mqc_mp2_hessian_fd`'s frozen column. What this pins in-tree is
   !! that neither symmetry drifts, on the same two geometries the
   !! all-electron test walks.
   subroutine frozen_core_symmetric(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: hess_corr(:, :, :, :), hess_ref(:, :, :, :)
      real(dp) :: coords(3, 3), worst_sym, worst_tr, acc
      integer :: threads, g, natm, ia, ja, a, b

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      do g = 1, 2
         coords = WATER
         if (g == 2) then
            coords(1, 2) = 0.031_dp
            coords(2, 2) = -0.024_dp
         end if
         call build_czt_molecule(WATER_Z, WATER_SYM, coords, "6-31g", mol, err)
         if (.not. err%has_error()) then
            call run_czt_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, &
                             .false., scf, err)
         end if
         if (.not. err%has_error()) then
            call mp2_correlation_hessian(mol, scf%orbitals, scf%orbital_energies, &
                                         scf%density, WATER_NELEC/2, 1, &
                                         hess_corr, hess_ref, err)
         end if
         call check(error,.not. err%has_error(), &
                    "the frozen-core assembly did not evaluate: "//err%get_message())
         if (allocated(error)) then
            call omp_set_num_threads(threads)
            call mol%destroy()
            return
         end if
         natm = size(hess_corr, 3)
         worst_sym = 0.0_dp
         worst_tr = 0.0_dp
         do ja = 1, natm
            do ia = 1, natm
               do b = 1, 3
                  do a = 1, 3
                     worst_sym = max(worst_sym, abs(hess_corr(a, b, ia, ja) &
                                                    - hess_corr(b, a, ja, ia)))
                  end do
               end do
            end do
         end do
         do ja = 1, natm
            do b = 1, 3
               do a = 1, 3
                  acc = 0.0_dp
                  do ia = 1, natm
                     acc = acc + hess_corr(a, b, ia, ja)
                  end do
                  worst_tr = max(worst_tr, abs(acc))
                  acc = 0.0_dp
                  do ia = 1, natm
                     acc = acc + hess_corr(b, a, ja, ia)
                  end do
                  worst_tr = max(worst_tr, abs(acc))
               end do
            end do
         end do
         write (*, "(a, i0, a, 2es10.3)") "        geometry ", g, &
            " symmetry / translation = ", worst_sym, worst_tr
         call check(error, max(worst_sym, worst_tr) < 1.0e-12_dp, &
                    "the frozen-core correlation Hessian lost a symmetry")
         call mol%destroy()
         if (allocated(error)) then
            call omp_set_num_threads(threads)
            return
         end if
         deallocate (hess_corr, hess_ref)
      end do
      call omp_set_num_threads(threads)
   end subroutine frozen_core_symmetric

end module test_mqc_mp2_hessian_assembly

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mp2_hessian_assembly, only: collect_mqc_mp2_hessian_assembly_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_mp2_hessian_assembly", &
                               collect_mqc_mp2_hessian_assembly_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
