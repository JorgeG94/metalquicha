!! Multipole integrals, against the identities that fix their convention
module test_mqc_czt_multipole
   !! These matrices are the input to every distributed multipole, every
   !! localization centroid and every dipole this program reports, and a
   !! transposed component or a mistaken origin convention produces numbers of
   !! the right magnitude throughout. Three things pin them without any external
   !! reference:
   !!
   !!   * **Shifting the origin shifts the dipole by the overlap.** For any
   !!     origin `R`, `<mu| r - R |nu> = <mu| r |nu> - R <mu|nu>`, exactly.
   !!     That is one identity per component and it fails loudly if the origin
   !!     is ignored, applied twice, or applied with the wrong sign -- the last
   !!     of which is invisible at the origin everything is usually computed at.
   !!   * **The matrices are symmetric.** `<mu|x|nu> = <nu|x|mu>` because the
   !!     operator is multiplication by a real function. An index swapped inside
   !!     the shell-pair loop breaks this while leaving the trace, and therefore
   !!     the reported dipole, untouched.
   !!   * **The quadrupole's full Cartesian tensor is symmetric in its two
   !!     indices**: the `xy` block equals the `yx` block. This module returns
   !!     all nine components rather than six on purpose, so the redundancy is
   !!     there to be checked.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_multipole, only: multipole_matrices
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_czt_multipole_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_czt_multipole_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("dipole_origin_shift_is_the_overlap", test_origin_shift), &
                  new_unittest("multipole_matrices_are_symmetric", test_symmetry), &
                  new_unittest("quadrupole_tensor_is_index_symmetric", test_quadrupole) &
                  ]
   end subroutine collect_mqc_czt_multipole_tests

   subroutine water(mol, err)
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err

      real(dp) :: c(3, 3)
      integer :: z(3)
      character(len=2) :: symbols(3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      ! A polarization set, so d functions are exercised: the Cartesian tensor
      ! components are where a stride mistake shows up, and an s/p-only basis
      ! cannot see one.
      call build_czt_molecule(z, symbols, c, "6-31g*", mol, err)
   end subroutine water

   subroutine test_origin_shift(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: at_zero(:, :, :), shifted(:, :, :), overlap(:, :)
      real(dp) :: r(3), worst, expected
      integer :: comp, i, j

      call water(mol, err)
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return

      call mol%overlap(overlap)
      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, at_zero, err)
      call check(error,.not. err%has_error(), "dipole at the origin")
      if (allocated(error)) return

      r = [0.37_dp, -1.10_dp, 2.05_dp]
      call multipole_matrices(mol, r, 1, shifted, err)
      call check(error,.not. err%has_error(), "dipole about a shifted origin")
      if (allocated(error)) return

      worst = 0.0_dp
      do comp = 1, 3
         do j = 1, mol%nao
            do i = 1, mol%nao
               expected = at_zero(i, j, comp) - r(comp)*overlap(i, j)
               worst = max(worst, abs(shifted(i, j, comp) - expected))
            end do
         end do
      end do
      call check(error, worst < 1.0e-12_dp, &
                 "the dipole does not shift with its origin as -R*S")
   end subroutine test_origin_shift

   subroutine test_symmetry(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: matrices(:, :, :)
      real(dp) :: worst
      integer :: order, comp, i, j

      call water(mol, err)
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return

      do order = 1, 3
         call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], order, matrices, err)
         call check(error,.not. err%has_error(), "multipole matrices")
         if (allocated(error)) return

         worst = 0.0_dp
         do comp = 1, size(matrices, 3)
            do j = 1, mol%nao
               do i = 1, j
                  worst = max(worst, abs(matrices(i, j, comp) - matrices(j, i, comp)))
               end do
            end do
         end do
         call check(error, worst < 1.0e-13_dp, &
                    "a multipole matrix is not symmetric in its basis indices")
         if (allocated(error)) return
      end do
   end subroutine test_symmetry

   subroutine test_quadrupole(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: quad(:, :, :)
      real(dp) :: worst
      integer :: a, b, i, j

      call water(mol, err)
      call check(error,.not. err%has_error(), "building the molecule")
      if (allocated(error)) return

      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 2, quad, err)
      call check(error,.not. err%has_error(), "quadrupole matrices")
      if (allocated(error)) return
      call check(error, size(quad, 3) == 9, "the full Cartesian tensor is nine components")
      if (allocated(error)) return

      ! Component (a,b) is stored with b fastest, so xy is 2 and yx is 4.
      worst = 0.0_dp
      do a = 1, 3
         do b = 1, 3
            do j = 1, mol%nao
               do i = 1, mol%nao
                  worst = max(worst, abs(quad(i, j, 3*(a - 1) + b) &
                                         - quad(i, j, 3*(b - 1) + a)))
               end do
            end do
         end do
      end do
      call check(error, worst < 1.0e-13_dp, &
                 "the quadrupole tensor is not symmetric under swapping its Cartesian indices")
   end subroutine test_quadrupole

end module test_mqc_czt_multipole

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_multipole, only: collect_mqc_czt_multipole_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_multipole", collect_mqc_czt_multipole_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
