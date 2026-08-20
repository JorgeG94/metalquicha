!! Second derivatives of the one-electron integrals
module test_mqc_hess_ints
   !! These are checked against identities rather than against stored numbers,
   !! because the identities say something the numbers cannot: that the two
   !! derivative orderings belong to the same integral and are consistent with
   !! each other.
   !!
   !! They were also checked against PySCF once, while being written -- all six
   !! blocks agreed to between 1e-16 and 1e-13. That is worth doing and not
   !! worth keeping: it needs Python, and it would pin the layout to one
   !! external library's conventions rather than to anything true.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_hess_ints, only: hess_1e_block, HESS_OVLP_II, HESS_OVLP_IJ, &
                                    HESS_KIN_II, HESS_KIN_IJ, HESS_NUC_II, HESS_NUC_IJ
   implicit none
   private

   public :: collect_mqc_hess_ints_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp], [3, 3])

contains

   subroutine collect_mqc_hess_ints_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("translations_leave_the_integrals_alone", translations_leave_them_alone), &
                  new_unittest("derivatives_on_one_centre_commute", derivatives_commute), &
                  new_unittest("every_block_is_populated", every_block_is_populated), &
                  new_unittest("components_are_where_they_belong", components_are_where_they_belong) &
                  ]
   end subroutine collect_mqc_hess_ints_tests

   subroutine build_water(mol, err, ok)
      type(libcint_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      ok = .not. err%has_error()
   end subroutine build_water

   subroutine translations_leave_them_alone(error)
      !! Moving both centres together cannot change a two-centre integral
      !!
      !! `S` depends on where the bra and ket sit, so under a rigid translation
      !!
      !!     (d/dA + d/dB)^2 S = 0   ->   S_AA + 2 S_AB + S_BB = 0
      !!
      !! and `S_BB` is `S_AA` with the two indices exchanged. **This is the
      !! check that matters most here**, because it is exactly what fails when
      !! only one derivative ordering is used, or when the two are mixed up: a
      !! Hessian built from the wrong combination is not slightly wrong, it is
      !! one where a rigid translation costs energy, and the symptom is
      !! translational modes that no longer come out at zero frequency.
      !!
      !! Nuclear attraction is excluded on purpose. Its operator depends on the
      !! nuclei too, so translating the two basis centres alone does not leave
      !! it alone and this identity does not apply to it.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: aa(:, :, :), ab(:, :, :)
      logical :: ok
      integer :: n, c
      real(dp) :: worst

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      call hess_1e_block(mol, HESS_OVLP_II, aa, err)
      call hess_1e_block(mol, HESS_OVLP_IJ, ab, err)
      call check(error,.not. err%has_error(), "overlap blocks failed")
      if (.not. allocated(error)) then
         n = size(aa, 1)
         worst = 0.0_dp
         do c = 1, 9
            worst = max(worst, maxval(abs(aa(:, :, c) + transpose(aa(:, :, c)) &
                                          + 2.0_dp*ab(:, :, c))))
         end do
         call check(error, worst < 1.0e-10_dp, &
                    "the overlap second derivatives are not translationally invariant")
      end if
      if (allocated(aa)) deallocate (aa)
      if (allocated(ab)) deallocate (ab)
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call hess_1e_block(mol, HESS_KIN_II, aa, err)
      call hess_1e_block(mol, HESS_KIN_IJ, ab, err)
      call check(error,.not. err%has_error(), "kinetic blocks failed")
      if (.not. allocated(error)) then
         worst = 0.0_dp
         do c = 1, 9
            worst = max(worst, maxval(abs(aa(:, :, c) + transpose(aa(:, :, c)) &
                                          + 2.0_dp*ab(:, :, c))))
         end do
         call check(error, worst < 1.0e-10_dp, &
                    "the kinetic second derivatives are not translationally invariant")
      end if
      call mol%destroy()
   end subroutine translations_leave_them_alone

   subroutine derivatives_commute(error)
      !! Two derivatives on the same centre commute, so xy equals yx
      !!
      !! True of the `ipip` blocks and **not** of the `ipXip` ones, where the
      !! two derivatives act on different centres and the component index means
      !! something asymmetric. Asserting it of both would be wrong; asserting it
      !! of neither would miss a transposed component block.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: aa(:, :, :)
      logical :: ok

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      call hess_1e_block(mol, HESS_OVLP_II, aa, err)
      call check(error,.not. err%has_error(), "overlap block failed")
      if (.not. allocated(error)) then
         ! xy against yx, xz against zx, yz against zy
         call check(error, maxval(abs(aa(:, :, 2) - aa(:, :, 4))) < 1.0e-12_dp, &
                    "xy and yx differ on one centre")
      end if
      if (.not. allocated(error)) then
         call check(error, maxval(abs(aa(:, :, 3) - aa(:, :, 7))) < 1.0e-12_dp, &
                    "xz and zx differ on one centre")
      end if
      if (.not. allocated(error)) then
         call check(error, maxval(abs(aa(:, :, 6) - aa(:, :, 8))) < 1.0e-12_dp, &
                    "yz and zy differ on one centre")
      end if
      call mol%destroy()
   end subroutine derivatives_commute

   subroutine every_block_is_populated(error)
      !! Each selector returns a differently-shaped, non-zero array
      !!
      !! A dispatch that fell through to the same entry point for two selectors
      !! would satisfy every identity above, since both copies would be equally
      !! valid integrals. What catches it is that the six are genuinely
      !! different numbers.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: a(:, :, :), b(:, :, :)
      logical :: ok
      integer :: k, sel(6)

      sel = [HESS_OVLP_II, HESS_OVLP_IJ, HESS_KIN_II, HESS_KIN_IJ, &
             HESS_NUC_II, HESS_NUC_IJ]

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      do k = 1, 6
         call hess_1e_block(mol, sel(k), a, err)
         call check(error,.not. err%has_error(), "a block failed to build")
         if (allocated(error)) exit
         call check(error, maxval(abs(a)) > 1.0e-8_dp, "a block came back empty")
         if (allocated(error)) exit
         if (k > 1) then
            call check(error, maxval(abs(a - b)) > 1.0e-8_dp, &
                       "two selectors returned the same integral")
            if (allocated(error)) exit
         end if
         if (allocated(b)) deallocate (b)
         b = a
         deallocate (a)
      end do
      call mol%destroy()
   end subroutine every_block_is_populated

   subroutine components_are_where_they_belong(error)
      !! Which of the nine components is which
      !!
      !! **The identities above cannot see this.** A permutation of the
      !! component index applied consistently to every block satisfies
      !! translational invariance and the commuting test alike, because each of
      !! those holds component by component -- so a wrong stride in the
      !! unpacking reproduces every structural property and still hands back xy
      !! where xz belongs.
      !!
      !! Pinned as a norm per component rather than as individual elements.
      !! Elements are a bad choice here and the first version of this test found
      !! out the hard way: water in this orientation is full of symmetry zeros,
      !! so two elements picked by eye can sit where a scrambled stride happens
      !! to read an equal value, and the test passes while the layout is wrong.
      !! A norm over the whole block cannot be fooled that way -- moving values
      !! between components changes what each one sums to.
      !!
      !! The numbers came from PySCF, which agreed with this implementation to
      !! between 1e-16 and 1e-13 across all six blocks when it was written.
      !! Recorded here because that comparison needs Python and cannot run in
      !! this suite.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: a(:, :, :)
      logical :: ok
      integer :: c
      real(dp), parameter :: TOL = 1.0e-5_dp
      !> `int1e_ipipovlp`, sum of absolute values per component
      real(dp), parameter :: OVLP_NORM(9) = [ &
                             28.464886_dp, 2.724628_dp, 2.566234_dp, &
                             2.724628_dp, 27.896234_dp, 3.212813_dp, &
                             2.566234_dp, 3.212813_dp, 27.811801_dp]
      !> `int1e_ipnucip`, which is a different shape entirely and pins the
      !> dispatch as well as the layout
      real(dp), parameter :: NUC_NORM(9) = [ &
                             1209.762119_dp, 57.194237_dp, 55.848487_dp, &
                             57.194237_dp, 1217.340172_dp, 67.229615_dp, &
                             55.848487_dp, 67.229615_dp, 1210.321515_dp]

      call build_water(mol, err, ok)
      call check(error, ok, "could not build water")
      if (allocated(error)) return

      call hess_1e_block(mol, HESS_OVLP_II, a, err)
      call check(error,.not. err%has_error(), "overlap block failed")
      if (.not. allocated(error)) then
         do c = 1, 9
            call check(error, sum(abs(a(:, :, c))), OVLP_NORM(c), thr=TOL, &
                       more="an overlap component is not where it belongs")
            if (allocated(error)) exit
         end do
      end if
      if (allocated(a)) deallocate (a)
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call hess_1e_block(mol, HESS_NUC_IJ, a, err)
      call check(error,.not. err%has_error(), "nuclear block failed")
      if (.not. allocated(error)) then
         do c = 1, 9
            call check(error, sum(abs(a(:, :, c))), NUC_NORM(c), thr=1.0e-4_dp, &
                       more="a nuclear component is not where it belongs")
            if (allocated(error)) exit
         end do
      end if
      call mol%destroy()
   end subroutine components_are_where_they_belong

end module test_mqc_hess_ints

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_hess_ints, only: collect_mqc_hess_ints_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_hess_ints", collect_mqc_hess_ints_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
