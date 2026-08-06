module test_mqc_json_basis_reader
   !! Pins the Basis Set Exchange JSON reader against transcribed reference data.
   !!
   !! This reader is the only one left -- the Gaussian94 and GAMESS readers
   !! were removed -- so it has no sibling to cross-check against and every
   !! calculation the code does now depends on it alone. The def2-SVP oxygen
   !! shells below were transcribed from the Gaussian94 file that used to sit
   !! in basis_sets/ (Turbomole 7.3 data, via BSE), which makes them an
   !! independent witness: they came out of a different file, in a different
   !! format, written by different tooling. If the JSON reader ever starts
   !! mangling exponents, dropping shells or mis-assigning angular momenta,
   !! these literals will not move with it.
   !!
   !! It also pins the behaviour that motivated going JSON-only: a combined SP
   !! shell must split into separate S and P shells sharing exponents, which
   !! Gaussian94 could only express by being downloaded pre-split.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_json_basis_reader, only: read_json_basis_element
   use mqc_cgto, only: atomic_basis_type
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_json_basis_reader_tests

   character(len=*), parameter :: DEF2_SVP = "../basis_sets/def2-svp.json"
   character(len=*), parameter :: STO_3G = "../basis_sets/sto-3g.json"

   real(dp), parameter :: TIGHT = 1.0e-10_dp  !! Basis data must survive exactly

   ! def2-SVP oxygen: 6 shells, S(5) S(1) S(1) P(3) P(1) D(1)
   integer, parameter :: O_NSHELLS = 6
   integer, parameter :: O_ANG_MOM(O_NSHELLS) = [0, 0, 0, 1, 1, 2]
   integer, parameter :: O_NPRIM(O_NSHELLS) = [5, 1, 1, 3, 1, 1]

   real(dp), parameter :: O_S1_EXPONENTS(5) = [ &
                          2266.1767785_dp, 340.87010191_dp, 77.363135167_dp, &
                          21.479644940_dp, 6.6589433124_dp]
   real(dp), parameter :: O_S1_COEFFICIENTS(5) = [ &
                          -0.53431809926e-02_dp, -0.39890039230e-01_dp, -0.17853911985_dp, &
                          -0.46427684959_dp, -0.44309745172_dp]
   real(dp), parameter :: O_P1_EXPONENTS(3) = [ &
                          17.721504317_dp, 3.8635505440_dp, 1.0480920883_dp]
   real(dp), parameter :: O_P1_COEFFICIENTS(3) = [ &
                          0.43394573193e-01_dp, 0.23094120765_dp, 0.51375311064_dp]
   real(dp), parameter :: O_D_EXPONENT = 1.2000000_dp

contains

   subroutine collect_mqc_json_basis_reader_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("def2svp_oxygen_shell_structure", test_oxygen_structure), &
                  new_unittest("def2svp_oxygen_primitives", test_oxygen_primitives), &
                  new_unittest("sp_shell_is_split", test_sp_shell_split), &
                  new_unittest("missing_element_errors", test_missing_element) &
                  ]
   end subroutine collect_mqc_json_basis_reader_tests

   subroutine read_oxygen(basis, ok, error)
      !! Parse def2-SVP oxygen, failing the test if the file will not load
      type(atomic_basis_type), intent(out) :: basis
      logical, intent(out) :: ok
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: parse_error

      call read_json_basis_element(DEF2_SVP, "O", basis, parse_error)
      ok = .not. parse_error%has_error()
      call check(error, ok, "def2-svp.json must parse (run cmake to extract it): "// &
                 trim(parse_error%get_message()))
   end subroutine read_oxygen

   subroutine test_oxygen_structure(error)
      !! The shell count, angular momenta and contraction lengths must match
      type(error_type), allocatable, intent(out) :: error
      type(atomic_basis_type) :: basis
      logical :: ok
      integer :: ishell

      call read_oxygen(basis, ok, error)
      if (allocated(error) .or. .not. ok) return

      call check(error, basis%nshells, O_NSHELLS, "def2-SVP oxygen has six shells")
      if (allocated(error)) return

      do ishell = 1, O_NSHELLS
         call check(error, basis%shells(ishell)%ang_mom, O_ANG_MOM(ishell), &
                    "wrong angular momentum on an oxygen shell")
         if (allocated(error)) return
         call check(error, basis%shells(ishell)%nfunc, O_NPRIM(ishell), &
                    "wrong primitive count on an oxygen shell")
         if (allocated(error)) return
      end do
   end subroutine test_oxygen_structure

   subroutine test_oxygen_primitives(error)
      !! Exponents and contraction coefficients must survive to full precision
      type(error_type), allocatable, intent(out) :: error
      type(atomic_basis_type) :: basis
      logical :: ok
      integer :: iprim

      call read_oxygen(basis, ok, error)
      if (allocated(error) .or. .not. ok) return
      if (basis%nshells < O_NSHELLS) return

      do iprim = 1, size(O_S1_EXPONENTS)
         call check(error, basis%shells(1)%exponents(iprim), O_S1_EXPONENTS(iprim), &
                    "contracted S exponent", thr=TIGHT)
         if (allocated(error)) return
         call check(error, basis%shells(1)%coefficients(iprim), O_S1_COEFFICIENTS(iprim), &
                    "contracted S coefficient", thr=TIGHT)
         if (allocated(error)) return
      end do

      do iprim = 1, size(O_P1_EXPONENTS)
         call check(error, basis%shells(4)%exponents(iprim), O_P1_EXPONENTS(iprim), &
                    "contracted P exponent", thr=TIGHT)
         if (allocated(error)) return
         call check(error, basis%shells(4)%coefficients(iprim), O_P1_COEFFICIENTS(iprim), &
                    "contracted P coefficient", thr=TIGHT)
         if (allocated(error)) return
      end do

      ! The polarization D is the last shell, and the one a truncated parse
      ! would most easily lose.
      call check(error, basis%shells(6)%exponents(1), O_D_EXPONENT, &
                 "polarization D exponent", thr=TIGHT)
   end subroutine test_oxygen_primitives

   subroutine test_sp_shell_split(error)
      !! STO-3G oxygen is one S shell plus one SP shell in the JSON. The reader
      !! must present that as three shells -- S, S, P -- with the S and P of
      !! the combined entry sharing exponents.
      type(error_type), allocatable, intent(out) :: error
      type(atomic_basis_type) :: basis
      type(error_t) :: parse_error
      integer :: iprim

      call read_json_basis_element(STO_3G, "O", basis, parse_error)
      call check(error,.not. parse_error%has_error(), "sto-3g.json must parse: "// &
                 trim(parse_error%get_message()))
      if (allocated(error)) return

      call check(error, basis%nshells, 3, "STO-3G oxygen must split into three shells")
      if (allocated(error)) return
      call check(error, basis%shells(1)%ang_mom, 0, "shell 1 is S")
      if (allocated(error)) return
      call check(error, basis%shells(2)%ang_mom, 0, "shell 2 is the S of the SP pair")
      if (allocated(error)) return
      call check(error, basis%shells(3)%ang_mom, 1, "shell 3 is the P of the SP pair")
      if (allocated(error)) return

      do iprim = 1, basis%shells(2)%nfunc
         call check(error, basis%shells(2)%exponents(iprim), basis%shells(3)%exponents(iprim), &
                    "a split SP shell must share its exponents", thr=1.0e-12_dp)
         if (allocated(error)) return
      end do
   end subroutine test_sp_shell_split

   subroutine test_missing_element(error)
      !! An element absent from the file must be an error, not an empty basis
      type(error_type), allocatable, intent(out) :: error
      type(atomic_basis_type) :: basis
      type(error_t) :: parse_error

      ! Uranium, not iron: def2-svp.json carries the whole periodic table up to
      ! radon, so only an element def2-SVP never defined at all is a stable
      ! stand-in for "absent from the file".
      call read_json_basis_element(DEF2_SVP, "U", basis, parse_error)
      call check(error, parse_error%has_error(), "a missing element must report an error")
   end subroutine test_missing_element

end module test_mqc_json_basis_reader

program tester_mqc_json_basis_reader
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_json_basis_reader, only: collect_mqc_json_basis_reader_tests
   implicit none

   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [ &
                new_testsuite("mqc_json_basis_reader", collect_mqc_json_basis_reader_tests) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_json_basis_reader
