module test_mqc_czt_cartesian
   !! Pins that the angular form the reader decided actually reaches libcint.
   !!
   !! `test_mqc_basis_cartesian` checks the decision; this checks that the
   !! decision is acted on. The two are worth separating because the failure
   !! that hid for so long lived entirely in the gap between them: the basis
   !! file said Cartesian all along, and every integral call asked for
   !! spherical anyway.
   !!
   !! The observable is the basis function count. Six functions per d shell
   !! against five is the whole difference, and it is not subtle once you look
   !! at it -- water in 6-31G* is 19 functions Cartesian and 18 spherical.
   !! Nothing else in the program prints a number that distinguishes the two,
   !! which is why nothing caught it.
   !!
   !! Reference counts are arithmetic, not measurements. 6-31G* water is
   !! H(2s) x2 + O(3s 2p 1d): 2*2 + 3 + 2*3 + 6 = 19 Cartesian, and 18 with a
   !! 5d shell. cc-pVDZ water is 2*5 + (3*1 + 2*3 + 5) = 24 spherical.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use mqc_cgto, only: molecular_basis_type
   use mqc_czt_integrals, only: czt_molecule_t, build_df_tensor
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_czt_cartesian_tests

   character(len=*), parameter :: POPLE_STAR = "../basis_sets/6-31g_st_.json"
   character(len=*), parameter :: CC_PVDZ = "../basis_sets/cc-pvdz.json"
   character(len=*), parameter :: STO_3G = "../basis_sets/sto-3g.json"

   !! Water, in Bohr. The geometry the Cartesian reference energies were
   !! generated on, so the counts here name the same molecule as those.
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   real(dp), parameter :: WATER_COORDS(3, 3) = reshape([ &
                                                       0.0_dp, 0.0_dp, 0.225374_dp, &
                                                       0.0_dp, 1.442316_dp, -0.901497_dp, &
                                                       0.0_dp, -1.442316_dp, -0.901497_dp], [3, 3])

contains

   subroutine collect_mqc_czt_cartesian_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("pople_star_water_is_19_functions", test_cartesian_count), &
                  new_unittest("ccpvdz_water_is_24_functions", test_spherical_count), &
                  new_unittest("sto3g_water_is_7_functions", test_minimal_count), &
                  new_unittest("overlap_carries_the_cartesian_signature", test_cartesian_overlap), &
                  new_unittest("mixed_fitting_bases_refused", test_mixed_fit_refused), &
                  new_unittest("matched_fitting_bases_accepted", test_matched_fit_ok) &
                  ]
   end subroutine collect_mqc_czt_cartesian_tests

   subroutine build_water(path, mol, ok, error)
      !! Pack water in the given basis, failing the test if it will not load
      character(len=*), intent(in) :: path
      type(czt_molecule_t), intent(out) :: mol
      logical, intent(out) :: ok
      type(error_type), allocatable, intent(out) :: error

      type(molecular_basis_type) :: basis
      type(error_t) :: build_error
      character(len=2) :: symbols(3)

      symbols = ["O ", "H ", "H "]
      call build_molecular_basis_json(path, symbols, basis, build_error)
      ok = .not. build_error%has_error()
      call check(error, ok, trim(path)//" must parse: "//trim(build_error%get_message()))
      if (.not. ok) return

      call mol%build(WATER_Z, WATER_COORDS, basis, build_error)
      ok = .not. build_error%has_error()
      call check(error, ok, "libcint must pack water: "//trim(build_error%get_message()))
      call basis%destroy()
   end subroutine build_water

   subroutine test_cartesian_count(error)
      !! The number this whole change exists for
      !!
      !! 18 here means the routing did not happen and 6-31G* is being computed
      !! as a smaller basis with the same name -- which converges perfectly
      !! well, 1.4 mHartree away from anything published.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      logical :: ok

      call build_water(POPLE_STAR, mol, ok, error)
      if (allocated(error) .or. .not. ok) return

      call check(error, mol%cartesian, "6-31G* water must reach libcint as Cartesian")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      call check(error, mol%nao, 19, "6-31G* water is 19 Cartesian basis functions, not 18")
      call mol%destroy()
   end subroutine test_cartesian_count

   subroutine test_spherical_count(error)
      !! The spherical path has to be exactly what it was
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      logical :: ok

      call build_water(CC_PVDZ, mol, ok, error)
      if (allocated(error) .or. .not. ok) return

      call check(error,.not. mol%cartesian, "cc-pVDZ water is spherical")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      call check(error, mol%nao, 24, "cc-pVDZ water is 24 spherical basis functions")
      call mol%destroy()
   end subroutine test_spherical_count

   subroutine test_minimal_count(error)
      !! A set with nothing above p is unaffected either way
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      logical :: ok

      call build_water(STO_3G, mol, ok, error)
      if (allocated(error) .or. .not. ok) return

      call check(error, mol%nao, 7, "STO-3G water is 7 basis functions in either form")
      call mol%destroy()
   end subroutine test_minimal_count

   subroutine test_cartesian_overlap(error)
      !! S must come back 19 square, and carry the Cartesian d signature
      !!
      !! The count alone would also be satisfied by a Cartesian shell dimension
      !! handed to a spherical integral call, which is a buffer overrun rather
      !! than a wrong answer and would not necessarily show. So this looks at
      !! the values.
      !!
      !! The signature is the d block's diagonal. A Cartesian d shell is not
      !! orthonormal under libcint's convention -- xx, yy and zz come out three
      !! times xy, xz and yz, because the two kinds of function are normalised
      !! against different angular integrals -- while a spherical d shell has
      !! five identical diagonal entries. Thirteen entries at one and six
      !! splitting three to one is a shape only a Cartesian build produces.
      !!
      !! The unequal diagonal is not an error and does not need fixing: an SCF
      !! solves F C = S C e, so any nonsingular scaling of the functions leaves
      !! the energy alone. The two values are pinned against PySCF's
      !! `int1e_ovlp_cart` on this molecule, which produces the same pair.
      !!
      !! The s and p entries are one only to about 1e-9, and that is the basis
      !! file rather than the integrals: BSE writes contraction coefficients to
      !! ten significant figures, so a contracted function is normalised to ten
      !! figures and no further.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: s(:, :)
      real(dp) :: diagonal(19)
      logical :: ok
      integer :: i, unit_count

      ! PySCF, int1e_ovlp_cart on this geometry in 6-31G*.
      real(dp), parameter :: D_AXIAL = 2.5132741229_dp    !! xx, yy, zz
      real(dp), parameter :: D_OFF_AXIS = 0.8377580410_dp !! xy, xz, yz

      call build_water(POPLE_STAR, mol, ok, error)
      if (allocated(error) .or. .not. ok) return

      call mol%overlap(s)
      call check(error, size(s, 1), 19, "the overlap must be 19 by 19")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call check(error, maxval(abs(s - transpose(s))), 0.0_dp, thr=1.0e-12_dp, &
                 message="the overlap must be symmetric")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      do i = 1, 19
         diagonal(i) = s(i, i)
      end do
      call check(error, minval(diagonal) > 0.0_dp, &
                 "no basis function may have a non-positive self-overlap")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      ! s and p are normalised in either convention, so the ones count them.
      unit_count = count(abs(diagonal - 1.0_dp) < 1.0e-8_dp)
      call check(error, unit_count, 13, &
                 "13 s and p functions must be normalised, leaving the 6 Cartesian d")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call check(error, count(abs(diagonal - D_AXIAL) < 1.0e-8_dp), 3, &
                 "xx, yy and zz share one self-overlap")
      if (.not. allocated(error)) then
         call check(error, count(abs(diagonal - D_OFF_AXIS) < 1.0e-8_dp), 3, &
                    "xy, xz and yz share the other, three times smaller")
      end if
      call mol%destroy()
   end subroutine test_cartesian_overlap

   subroutine test_mixed_fit_refused(error)
      !! Cartesian orbitals fitted with a spherical auxiliary basis is refused
      !!
      !! libcint builds all three centres of a fitting integral in one form, so
      !! there is nothing to route -- one of the two bases would have to be
      !! read as something it is not, and that is the failure mode this whole
      !! change exists to remove rather than relocate.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: orb, aux
      real(dp), allocatable :: b(:, :)
      type(error_t) :: df_error
      logical :: ok

      call build_water(POPLE_STAR, orb, ok, error)
      if (allocated(error) .or. .not. ok) return
      call build_water(CC_PVDZ, aux, ok, error)
      if (allocated(error) .or. .not. ok) then
         call orb%destroy()
         return
      end if

      call build_df_tensor(orb, aux, b, df_error)
      call check(error, df_error%has_error(), &
                 "a Cartesian orbital basis and a spherical auxiliary one must be refused")
      if (.not. allocated(error)) then
         call check(error, index(df_error%get_message(), "Cartesian") > 0, &
                    "the refusal has to name the forms: "//trim(df_error%get_message()))
      end if

      call orb%destroy()
      call aux%destroy()
   end subroutine test_mixed_fit_refused

   subroutine test_matched_fit_ok(error)
      !! Two Cartesian bases fit, and the tensor comes back the right shape
      !!
      !! The point of the check is that the refusal above is about the mismatch
      !! and not about Cartesian fitting being unsupported: the three- and
      !! two-centre Cartesian entry points exist and are wired up.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: orb, aux
      real(dp), allocatable :: b(:, :)
      type(error_t) :: df_error
      logical :: ok

      call build_water(POPLE_STAR, orb, ok, error)
      if (allocated(error) .or. .not. ok) return
      call build_water(POPLE_STAR, aux, ok, error)
      if (allocated(error) .or. .not. ok) then
         call orb%destroy()
         return
      end if

      call build_df_tensor(orb, aux, b, df_error)
      call check(error,.not. df_error%has_error(), &
                 "fitting a Cartesian basis with a Cartesian auxiliary must work: "// &
                 trim(df_error%get_message()))
      if (.not. allocated(error)) then
         call check(error, size(b, 1), 19*19, "B is (nao*nao, naux)")
      end if
      if (.not. allocated(error)) then
         call check(error, size(b, 2), 19, "B is (nao*nao, naux)")
      end if

      call orb%destroy()
      call aux%destroy()
   end subroutine test_matched_fit_ok

end module test_mqc_czt_cartesian

program tester_mqc_czt_cartesian
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_cartesian, only: collect_mqc_czt_cartesian_tests
   implicit none

   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [ &
                new_testsuite("mqc_czt_cartesian", collect_mqc_czt_cartesian_tests) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_czt_cartesian
