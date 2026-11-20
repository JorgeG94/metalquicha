!! Nuclear magnetic shielding tensors, against PySCF and against themselves
module test_mqc_czt_nmr
   !! The pinned tensors are `pyscf.prop.nmr.RHF` on water in 6-31G at this
   !! geometry, fed this repository's basis: all nine elements of each
   !! nucleus's tensor, which fixes every sign and index convention at once.
   !! The tensor is not symmetric, so its two off-diagonal elements in the
   !! molecular plane are separate checks.
   !!
   !! Beside them, what needs no other program:
   !!
   !!   * moving the molecule leaves every tensor unchanged, which is the
   !!     property gauge-including orbitals exist to give and which a
   !!     common-origin calculation in this basis does not have;
   !!   * the diamagnetic and paramagnetic parts add up to the total;
   !!   * the summary's isotropic value is a third of the trace and its
   !!     anisotropy that of a tensor whose principal values are known.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_nmr, only: czt_nmr_shielding
   use mqc_nmr_summary, only: shielding_summary
   implicit none
   private

   public :: collect_mqc_czt_nmr_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   ! PySCF, water/6-31G, ppm. Rows are the field's direction.
   real(dp), parameter :: O_REF(3, 3) = reshape( &
                          [331.646944_dp, 0.0_dp, -17.748779_dp, &
                           0.0_dp, 307.327544_dp, 0.0_dp, &
                           -17.701335_dp, 0.0_dp, 341.403792_dp], [3, 3])
   real(dp), parameter :: H1_REF(3, 3) = reshape( &
                          [27.797327_dp, 0.0_dp, -2.421815_dp, &
                           0.0_dp, 25.909842_dp, 0.0_dp, &
                           -0.918204_dp, 0.0_dp, 44.357334_dp], [3, 3])
   real(dp), parameter :: H2_REF(3, 3) = reshape( &
                          [42.551405_dp, 0.0_dp, -4.740960_dp, &
                           0.0_dp, 25.938788_dp, 0.0_dp, &
                           -6.247660_dp, 0.0_dp, 29.677511_dp], [3, 3])
   real(dp), parameter :: PPM_TOL = 5.0e-5_dp
      !! The references are printed to six decimals.

contains

   subroutine collect_mqc_czt_nmr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("water_tensors_match_pyscf", test_water), &
                  new_unittest("translation_leaves_the_tensors_unchanged", test_translation), &
                  new_unittest("summary_of_a_known_tensor", test_summary) &
                  ]
   end subroutine collect_mqc_czt_nmr_tests

   subroutine shield_water(shift, sigma, dia, para, err)
      !! Water in 6-31G, moved by `shift`, and its shielding tensors
      real(dp), intent(in) :: shift(3)
      real(dp), allocatable, intent(out) :: sigma(:, :, :), dia(:, :, :), para(:, :, :)
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp) :: c(3, 3)
      integer :: i

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      do i = 1, 3
         c(:, i) = c(:, i) + shift
      end do
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "6-31g", mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, 10, 200, 1.0e-13_dp, 1.0e-11_dp, .false., scf, err, &
                       grad_tol=1.0e-10_dp)
      if (err%has_error()) return
      call czt_nmr_shielding(mol, scf%orbitals, scf%orbital_energies, scf%n_occupied, &
                             sigma, err, diamagnetic=dia, paramagnetic=para)
      call mol%destroy()
   end subroutine shield_water

   subroutine test_water(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: sigma(:, :, :), dia(:, :, :), para(:, :, :)
      integer :: iatom, x

      call shield_water([0.0_dp, 0.0_dp, 0.0_dp], sigma, dia, para, err)
      call check(error,.not. err%has_error(), "shielding: "//err%get_full_trace())
      if (allocated(error)) return
      do iatom = 1, 3
         write (*, "(4x,a,i0)") "atom ", iatom
         do x = 1, 3
            write (*, "(6x,3f14.6)") sigma(x, :, iatom)
         end do
      end do
      call check(error, maxval(abs(sigma(:, :, 1) - O_REF)) < PPM_TOL, &
                 "oxygen's shielding tensor does not match PySCF's")
      if (allocated(error)) return
      call check(error, maxval(abs(sigma(:, :, 2) - H1_REF)) < PPM_TOL, &
                 "the first hydrogen's shielding tensor does not match PySCF's")
      if (allocated(error)) return
      call check(error, maxval(abs(sigma(:, :, 3) - H2_REF)) < PPM_TOL, &
                 "the second hydrogen's shielding tensor does not match PySCF's")
      if (allocated(error)) return
      call check(error, maxval(abs(sigma - dia - para)) < 1.0e-10_dp, &
                 "the two parts do not add up to the total")
      if (allocated(error)) return
      ! The tensor is not symmetric, and the reference says by how much; a
      ! transposed result would pass a symmetric check and fail this one.
      call check(error, abs(sigma(1, 3, 2) - sigma(3, 1, 2)) > 1.0_dp, &
                 "the hydrogen tensor should be visibly asymmetric")
   end subroutine test_water

   subroutine test_translation(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: here(:, :, :), there(:, :, :), dia(:, :, :), para(:, :, :)
      real(dp), allocatable :: dia_moved(:, :, :)

      call shield_water([0.0_dp, 0.0_dp, 0.0_dp], here, dia, para, err)
      call shield_water([3.0_dp, -5.0_dp, 7.0_dp], there, dia_moved, para, err)
      call check(error,.not. err%has_error(), "shielding: "//err%get_full_trace())
      if (allocated(error)) return
      write (*, "(4x,a,es10.2)") "largest change on moving the molecule ", &
         maxval(abs(there - here))
      call check(error, maxval(abs(there - here)) < 1.0e-6_dp, &
                 "the shielding changed when the molecule was moved")
   end subroutine test_translation

   subroutine test_summary(error)
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: tensor(3, 3), iso, aniso

      ! Principal values 10, 20 and 60: the mean is 30 and 60 is furthest
      ! from it, so the anisotropy is 60 - 15. The antisymmetric part must
      ! not enter either number.
      tensor = 0.0_dp
      tensor(1, 1) = 10.0_dp
      tensor(2, 2) = 20.0_dp
      tensor(3, 3) = 60.0_dp
      tensor(1, 2) = 4.0_dp
      tensor(2, 1) = -4.0_dp
      call shielding_summary(tensor, iso, aniso)
      call check(error, abs(iso - 30.0_dp) < 1.0e-12_dp, "isotropic value")
      if (allocated(error)) return
      call check(error, abs(aniso - 45.0_dp) < 1.0e-10_dp, "anisotropy")
      if (allocated(error)) return

      ! And when the outlier is the smallest principal value: -30, 20, 40
      ! average 10, -30 is furthest, so the anisotropy is -30 - 30.
      tensor = 0.0_dp
      tensor(1, 1) = 40.0_dp
      tensor(2, 2) = -30.0_dp
      tensor(3, 3) = 20.0_dp
      call shielding_summary(tensor, iso, aniso)
      call check(error, abs(aniso + 60.0_dp) < 1.0e-10_dp, &
                 "anisotropy when the outlying principal value is the lowest")
   end subroutine test_summary

end module test_mqc_czt_nmr

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_nmr, only: collect_mqc_czt_nmr_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_nmr", collect_mqc_czt_nmr_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
