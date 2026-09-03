!! Valence virtual orbitals
module test_mqc_quao_vvo
   !! The chemically meaningful part of an SCF's virtual space, recovered by
   !! projection onto free-atom minimal-basis orbitals.
   !!
   !! Two properties are worth testing and they are different claims. The first
   !! is that the projection separates cleanly -- the retained singular values
   !! sit near one and the rejected ones an order of magnitude below, so the cut
   !! falls in a gap rather than through a continuum. Paper I Table I reports
   !! 0.99999 against 0.105-0.272 across eight molecules. The second is that the
   !! result barely moves with the orbital basis, which is the property the whole
   !! construction exists for: an ordinary virtual space changes beyond
   !! recognition between cc-pVDZ and cc-pVTZ, and the valence-virtual space
   !! should not.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t
   use mqc_aambs, only: aambs_dimensions, aambs_dimensions_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, &
                                mixed_basis_overlap
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_quao, only: build_aambs_molecule, mo_aambs_overlap, &
                           valence_virtual_orbitals, vvo_result_t
   implicit none
   private

   public :: collect_mqc_quao_vvo_tests

   integer, parameter :: N_DIM = 3
   real(dp), parameter :: WATER(N_DIM, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]

contains

   subroutine collect_mqc_quao_vvo_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("singular_values_separate", test_gap), &
                  new_unittest("vvos_are_orthonormal_virtuals", test_orthonormal), &
                  new_unittest("count_is_basis_independent", test_basis_independence) &
                  ]
   end subroutine collect_mqc_quao_vvo_tests

   subroutine water_vvos(basis_name, vvo, orbitals, overlap, dims, err, ok)
      !! An SCF on water and the valence-virtual orbitals that follow
      character(len=*), intent(in) :: basis_name
      type(vvo_result_t), intent(out) :: vvo
      real(dp), allocatable, intent(out) :: orbitals(:, :)
      real(dp), allocatable, intent(out) :: overlap(:, :)
      type(aambs_dimensions_t), intent(out) :: dims
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(czt_molecule_t) :: mol, aambs
      type(rhf_result_t) :: scf
      real(dp), allocatable :: mixed(:, :), s_mbs(:, :), projection(:, :)

      ok = .false.
      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, basis_name, mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      if (err%has_error()) return
      if (.not. scf%converged) return

      call build_aambs_molecule(WATER_Z, WATER_SYM, WATER, aambs, err)
      if (err%has_error()) return

      call mol%overlap(overlap)
      call aambs%overlap(s_mbs)
      call mixed_basis_overlap(mol, aambs, mixed, err)
      if (err%has_error()) return

      call mo_aambs_overlap(scf%orbitals, mixed, s_mbs, projection, err)
      if (err%has_error()) return

      call aambs_dimensions(WATER_Z, 10, dims, err)
      if (err%has_error()) return

      call valence_virtual_orbitals(scf%orbitals, projection, dims, vvo, err)
      orbitals = scf%orbitals
      ok = .not. err%has_error()

      call mol%destroy()
      call aambs%destroy()
   end subroutine water_vvos

   subroutine test_gap(error)
      !! The cut falls in a gap, not through a continuum
      type(error_type), allocatable, intent(out) :: error
      type(vvo_result_t) :: vvo
      type(aambs_dimensions_t) :: dims
      type(error_t) :: err
      real(dp), allocatable :: orbitals(:, :), overlap(:, :)
      logical :: ok

      call water_vvos("cc-pvdz", vvo, orbitals, overlap, dims, err, ok)
      call check(error, ok, "the water valence-virtual extraction should succeed")
      if (allocated(error)) return

      ! Water in a minimal basis has seven orbitals and five are occupied, so
      ! two antibonding valence orbitals are missing from the occupied space.
      call check(error, vvo%n_vvo, 2, "water should have two valence-virtual orbitals")
      if (allocated(error)) return

      ! Paper I's numbers are 0.99999-ish retained against 0.105-0.272 rejected.
      ! The bounds here are loose enough not to encode one basis's arithmetic and
      ! tight enough that a projection onto the wrong space fails them.
      call check(error, vvo%smallest_retained > 0.95_dp, &
                 "the retained singular values should be close to one")
      if (allocated(error)) return
      call check(error, vvo%largest_rejected < 0.5_dp, &
                 "the rejected singular values should be well below one")
      if (allocated(error)) return
      call check(error, vvo%smallest_retained - vvo%largest_rejected > 0.4_dp, &
                 "there should be a clear gap where the cut falls -- without one, "// &
                 "the valence space is not separable and the count is arbitrary")
      if (allocated(error)) return

      ! Descending, which every downstream index assumes.
      call check(error, vvo%singular_values(1) >= vvo%singular_values(2), &
                 "singular values should come back descending")
   end subroutine test_gap

   subroutine test_orthonormal(error)
      !! They are orthonormal, and they live in the virtual space
      !!
      !! Orthonormality is inherited rather than imposed -- they are orthogonal
      !! combinations of orthonormal virtuals -- so a failure here means the
      !! rotation matrix is not orthogonal, which would mean the eigenvectors
      !! were mixed up with something else.
      type(error_type), allocatable, intent(out) :: error
      type(vvo_result_t) :: vvo
      type(aambs_dimensions_t) :: dims
      type(error_t) :: err
      real(dp), allocatable :: orbitals(:, :), overlap(:, :), work(:, :), gram(:, :)
      logical :: ok
      integer :: i, n

      call water_vvos("cc-pvdz", vvo, orbitals, overlap, dims, err, ok)
      call check(error, ok, "the extraction should succeed")
      if (allocated(error)) return

      n = vvo%n_vvo
      allocate (work(size(overlap, 1), n), gram(n, n))
      call pic_gemm(overlap, vvo%orbitals, work)
      call pic_gemm(vvo%orbitals, work, gram, transa="T")

      do i = 1, n
         call check(error, abs(gram(i, i) - 1.0_dp) < 1.0e-10_dp, &
                    "each valence-virtual orbital should be normalized")
         if (allocated(error)) return
      end do
      call check(error, maxval(abs(gram - transpose(gram))) < 1.0e-12_dp, &
                 "the Gram matrix should be symmetric")
      if (allocated(error)) return
      if (n > 1) then
         call check(error, abs(gram(1, 2)) < 1.0e-10_dp, &
                    "valence-virtual orbitals should be mutually orthogonal")
         if (allocated(error)) return
      end if

      ! Built from virtual columns only, so they are orthogonal to every
      ! occupied orbital by construction. Asserted because it is the property
      ! that makes them a *virtual* space rather than a rotation of everything.
      ! `work` already holds S C_vvo from the Gram matrix above.
      call check(error, maxval(abs(matmul(transpose(orbitals(:, 1:dims%n_occupied)), &
                                          work))) < 1.0e-10_dp, &
                 "valence-virtual orbitals must be orthogonal to the occupied space")
   end subroutine test_orthonormal

   subroutine test_basis_independence(error)
      !! The count does not depend on the orbital basis
      !!
      !! This is the claim the construction is built on. The virtual space of
      !! cc-pVTZ is far larger than that of cc-pVDZ and looks nothing like it,
      !! but the number of valence-virtual orbitals is a property of the atoms,
      !! and the separation should stay clean in both.
      type(error_type), allocatable, intent(out) :: error
      type(vvo_result_t) :: dz, tz
      type(aambs_dimensions_t) :: dims_dz, dims_tz
      type(error_t) :: err
      real(dp), allocatable :: c_dz(:, :), c_tz(:, :), s_dz(:, :), s_tz(:, :)
      logical :: ok

      call water_vvos("cc-pvdz", dz, c_dz, s_dz, dims_dz, err, ok)
      call check(error, ok, "the double-zeta extraction should succeed")
      if (allocated(error)) return
      call water_vvos("cc-pvtz", tz, c_tz, s_tz, dims_tz, err, ok)
      call check(error, ok, "the triple-zeta extraction should succeed")
      if (allocated(error)) return

      call check(error, dz%n_vvo, tz%n_vvo, &
                 "the valence-virtual count is a property of the atoms, not of "// &
                 "the orbital basis")
      if (allocated(error)) return

      ! The virtual space grew by a large factor and the separation should not
      ! have degraded.
      call check(error, size(tz%singular_values) > 2*size(dz%singular_values), &
                 "the triple-zeta virtual space should be much larger")
      if (allocated(error)) return
      call check(error, tz%smallest_retained > 0.95_dp, &
                 "the retained singular values should stay near one in the "// &
                 "larger basis")
      if (allocated(error)) return
      call check(error, tz%smallest_retained - tz%largest_rejected > 0.4_dp, &
                 "the gap should survive the basis-set increase")
   end subroutine test_basis_independence

end module test_mqc_quao_vvo

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_quao_vvo, only: collect_mqc_quao_vvo_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_quao_vvo", collect_mqc_quao_vvo_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
