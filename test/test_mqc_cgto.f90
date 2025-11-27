module test_mqc_cgto
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_cgto, only: cgto_type, atomic_basis_type, molecular_basis_type
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_cgto_tests

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_cgto_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("cgto_allocate_and_destroy", test_cgto_allocate_destroy), &
                  new_unittest("cgto_num_basis_s_shell", test_cgto_num_basis_s), &
                  new_unittest("cgto_num_basis_p_shell", test_cgto_num_basis_p), &
                  new_unittest("cgto_num_basis_d_shell", test_cgto_num_basis_d), &
                  new_unittest("cgto_num_basis_f_shell", test_cgto_num_basis_f), &
                  new_unittest("atomic_basis_allocate_destroy", test_atomic_basis_allocate_destroy), &
                  new_unittest("atomic_basis_num_functions", test_atomic_basis_num_functions), &
                  new_unittest("molecular_basis_allocate_destroy", test_molecular_basis_allocate_destroy), &
                  new_unittest("molecular_basis_num_functions", test_molecular_basis_num_functions) &
                  ]
   end subroutine collect_mqc_cgto_tests

   subroutine test_cgto_allocate_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(cgto_type) :: shell

      call shell%allocate_arrays(3)

      call check(error, shell%nfunc, 3, "Should allocate 3 functions")
      if (allocated(error)) return

      call check(error, allocated(shell%exponents), "Exponents should be allocated")
      if (allocated(error)) return

      call check(error, allocated(shell%coefficients), "Coefficients should be allocated")
      if (allocated(error)) return

      call check(error, size(shell%exponents), 3, "Exponents array should have size 3")
      if (allocated(error)) return

      call shell%destroy()

      call check(error,.not. allocated(shell%exponents), "Exponents should be deallocated")
      if (allocated(error)) return

      call check(error, shell%nfunc, 0, "nfunc should be reset to 0")
   end subroutine test_cgto_allocate_destroy

   subroutine test_cgto_num_basis_s(error)
      type(error_type), allocatable, intent(out) :: error
      type(cgto_type) :: shell
      integer :: nbf

      shell%ang_mom = 0  ! S shell
      nbf = shell%num_basis_functions()

      ! S shell: (0+1)*(0+2)/2 = 1
      call check(error, nbf, 1, "S shell should have 1 basis function")

      call shell%destroy()
   end subroutine test_cgto_num_basis_s

   subroutine test_cgto_num_basis_p(error)
      type(error_type), allocatable, intent(out) :: error
      type(cgto_type) :: shell
      integer :: nbf

      shell%ang_mom = 1  ! P shell
      nbf = shell%num_basis_functions()

      ! P shell: (1+1)*(1+2)/2 = 3 (px, py, pz)
      call check(error, nbf, 3, "P shell should have 3 basis functions")

      call shell%destroy()
   end subroutine test_cgto_num_basis_p

   subroutine test_cgto_num_basis_d(error)
      type(error_type), allocatable, intent(out) :: error
      type(cgto_type) :: shell
      integer :: nbf

      shell%ang_mom = 2  ! D shell
      nbf = shell%num_basis_functions()

      ! D shell: (2+1)*(2+2)/2 = 6 (Cartesian)
      call check(error, nbf, 6, "D shell should have 6 Cartesian basis functions")

      call shell%destroy()
   end subroutine test_cgto_num_basis_d

   subroutine test_cgto_num_basis_f(error)
      type(error_type), allocatable, intent(out) :: error
      type(cgto_type) :: shell
      integer :: nbf

      shell%ang_mom = 3  ! F shell
      nbf = shell%num_basis_functions()

      ! F shell: (3+1)*(3+2)/2 = 10 (Cartesian)
      call check(error, nbf, 10, "F shell should have 10 Cartesian basis functions")

      call shell%destroy()
   end subroutine test_cgto_num_basis_f

   subroutine test_atomic_basis_allocate_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(atomic_basis_type) :: basis

      call basis%allocate_shells(2)

      call check(error, basis%nshells, 2, "Should allocate 2 shells")
      if (allocated(error)) return

      call check(error, allocated(basis%shells), "Shells should be allocated")
      if (allocated(error)) return

      call basis%destroy()

      call check(error,.not. allocated(basis%shells), "Shells should be deallocated")
      if (allocated(error)) return

      call check(error, basis%nshells, 0, "nshells should be reset")
   end subroutine test_atomic_basis_allocate_destroy

   subroutine test_atomic_basis_num_functions(error)
      type(error_type), allocatable, intent(out) :: error
      type(atomic_basis_type) :: basis
      integer :: nbf

      ! Create basis with 2 shells: S and P
      call basis%allocate_shells(2)
      basis%shells(1)%ang_mom = 0  ! S shell (1 function)
      basis%shells(2)%ang_mom = 1  ! P shell (3 functions)

      nbf = basis%num_basis_functions()

      ! Total: 1 + 3 = 4
      call check(error, nbf, 4, "Basis with S+P should have 4 functions")

      call basis%destroy()
   end subroutine test_atomic_basis_num_functions

   subroutine test_molecular_basis_allocate_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(molecular_basis_type) :: mol_basis

      call mol_basis%allocate_elements(3)

      call check(error, mol_basis%nelements, 3, "Should allocate 3 elements")
      if (allocated(error)) return

      call check(error, allocated(mol_basis%elements), "Elements should be allocated")
      if (allocated(error)) return

      call mol_basis%destroy()

      call check(error,.not. allocated(mol_basis%elements), "Elements should be deallocated")
      if (allocated(error)) return

      call check(error, mol_basis%nelements, 0, "nelements should be reset")
   end subroutine test_molecular_basis_allocate_destroy

   subroutine test_molecular_basis_num_functions(error)
      type(error_type), allocatable, intent(out) :: error
      type(molecular_basis_type) :: mol_basis
      integer :: nbf

      ! Create molecular basis with 2 atoms
      ! Atom 1: S shell (1 function)
      ! Atom 2: S + P shells (4 functions)
      call mol_basis%allocate_elements(2)

      call mol_basis%elements(1)%allocate_shells(1)
      mol_basis%elements(1)%shells(1)%ang_mom = 0  ! S

      call mol_basis%elements(2)%allocate_shells(2)
      mol_basis%elements(2)%shells(1)%ang_mom = 0  ! S
      mol_basis%elements(2)%shells(2)%ang_mom = 1  ! P

      nbf = mol_basis%num_basis_functions()

      ! Total: 1 + (1 + 3) = 5
      call check(error, nbf, 5, "Molecular basis should have 5 total functions")

      call mol_basis%destroy()
   end subroutine test_molecular_basis_num_functions

end module test_mqc_cgto

program tester_mqc_cgto
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_cgto, only: collect_mqc_cgto_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_cgto", collect_mqc_cgto_tests) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_cgto
