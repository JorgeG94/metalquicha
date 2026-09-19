!! The wavefunction stability analysis
module test_mqc_czt_stability
   !! Does the eigensolver find the lowest curvature of the matrix we handed it?
   !!
   !! Two questions, and they are separable, so they are separated.
   !!
   !! **Is the operator the electronic Hessian?** `mqc_czt_ov_hessian` recovers
   !! `(A+B)` from the response operator by arithmetic -- one factor of two in
   !! the trial density, one sign in the energy denominators -- and getting any
   !! of it wrong produces a curvature that is plausible and meaningless.
   !! `test_matches_dense_hessian` settles it by building the same matrix twice
   !! for the same molecule: once column by column through the operator, which
   !! is a direct Fock build per column, and once through `build_hessian`, which
   !! transforms a stored integral tensor. Two routes with nothing in common but
   !! the answer.
   !!
   !! **Does the eigensolver agree with a dense diagonalisation?** Against a
   !! real molecule this is a weak question: a converged closed-shell SCF is a
   !! minimum, so every solver says "stable" and the number it reports has no
   !! independent reference. So the comparison is made against a synthetic
   !! operator with the same structure and a spectrum chosen to be negative --
   !! the matrix is written down in the test, diagonalised with `pic_syev`, and
   !! the lowest eigenvalue compared with what the Davidson returns. Synthetic
   !! because the *number* has to be known independently, not because a real one
   !! was unavailable; the real one is checked above, as a matrix.
   !!
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp, default_int
   use pic_lapack_interfaces, only: pic_syev
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_cphf, only: build_hessian
   use mqc_czt_hessian, only: nuclear_response_t
   use mqc_czt_response, only: response_operator_t
   use mqc_czt_ov_hessian, only: ov_hessian_t, stability_result_t, build_scf_ov_hessian
   use mqc_czt_native_stability, only: native_stability_of_hessian, &
                                       native_scf_stability
   implicit none
   private

   public :: collect_mqc_czt_stability_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

   !! The synthetic cases are three occupied and four virtual orbitals: twelve
   !! rotations, which is more than the eigensolver's default trial space and
   !! so exercises the subspace expansion rather than filling the whole space
   !! on the first pass, while still diagonalising in no time at all.
   integer, parameter :: TOY_OCC = 3
   integer, parameter :: TOY_VIR = 4
   integer, parameter :: TOY_MO = TOY_OCC + TOY_VIR
   integer, parameter :: TOY_OV = TOY_OCC*TOY_VIR

   type, extends(response_operator_t) :: toy_response_t
      !! A response operator whose two-electron part is a matrix we wrote down
      !!
      !! Mimics `nuclear_response_t` exactly where it matters: the vector is
      !! `n_mo` by `n_occ`, the occupied rows are zeroed in the image, and the
      !! virtual rows come back divided by `e_i - e_a`. So the Hessian
      !! `mqc_czt_ov_hessian` recovers from it is `diag(gaps) + coupling`, and
      !! that is known to the last bit without any integrals.
      real(dp) :: gaps(TOY_OV) = 1.0_dp
      real(dp) :: coupling(TOY_OV, TOY_OV) = 0.0_dp
   contains
      procedure :: apply => toy_apply
      procedure :: length => toy_length
   end type toy_response_t

contains

   subroutine collect_mqc_czt_stability_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("stability_operator_is_the_dense_hessian", &
                               test_matches_dense_hessian), &
                  new_unittest("stability_refuses_without_an_operator", &
                               test_refuses_without_operator), &
                  new_unittest("native_lowest_eigenvalue", test_native_lowest_eigenvalue), &
                  new_unittest("native_reports_a_stable_curvature", &
                               test_native_reports_stable_curvature), &
                  new_unittest("native_water_is_a_minimum", test_native_water) &
                  ]
   end subroutine collect_mqc_czt_stability_tests

   ! ---- the synthetic operator -------------------------------------------

   pure function toy_length(this) result(n)
      !! Fixed by the module parameters, so `this` is the interface and not a
      !! source of information.
      class(toy_response_t), intent(in) :: this
      integer :: n

      n = TOY_MO*TOY_OCC
   end function toy_length

   subroutine toy_apply(this, vector, image, error)
      !! `Delta^-1 K`, on the layout the real operator uses
      class(toy_response_t), intent(inout) :: this
      real(dp), intent(in) :: vector(:)
      real(dp), intent(out) :: image(:)
      type(error_t), intent(inout) :: error

      real(dp) :: x(TOY_OV), y(TOY_OV)
      integer :: i, a, row

      if (error%has_error()) return

      do i = 1, TOY_OCC
         do a = 1, TOY_VIR
            x(a + (i - 1)*TOY_VIR) = vector(TOY_OCC + a + (i - 1)*TOY_MO)
         end do
      end do
      y = matmul(this%coupling, x)
      image = 0.0_dp
      do i = 1, TOY_OCC
         do a = 1, TOY_VIR
            row = a + (i - 1)*TOY_VIR
            ! The real operator divides by `e_i - e_a`, which is `-gaps`. That
            ! sign is the one this test exists to pin.
            image(TOY_OCC + a + (i - 1)*TOY_MO) = -y(row)/this%gaps(row)
         end do
      end do
   end subroutine toy_apply

   subroutine toy_hessian(operator, hessian, diagonal_shift)
      !! A synthetic electronic Hessian, and the operator that reproduces it
      !!
      !! `diagonal_shift` is the whole difference between the stable case and
      !! the unstable one: added to the coupling's diagonal it moves the
      !! spectrum bodily, so one number selects a positive definite matrix or
      !! one with a clear negative eigenvalue, and the test never has to assert
      !! what that eigenvalue is -- `pic_syev` says.
      type(toy_response_t), intent(out), target :: operator
      type(ov_hessian_t), intent(out) :: hessian
      real(dp), intent(in) :: diagonal_shift

      integer :: i, j

      ! Deterministic rather than random: a test that fails has to fail the
      ! same way on the next run.
      do i = 1, TOY_OV
         operator%gaps(i) = 0.4_dp + 0.15_dp*real(i, dp)
      end do
      do j = 1, TOY_OV
         do i = 1, TOY_OV
            operator%coupling(i, j) = 0.21_dp/real(i + j, dp)
         end do
         operator%coupling(j, j) = operator%coupling(j, j) + diagonal_shift
      end do

      hessian%response => operator
      hessian%n_occ = TOY_OCC
      hessian%n_vir = TOY_VIR
      hessian%n_mo = TOY_MO
      hessian%gaps = operator%gaps
   end subroutine toy_hessian

   subroutine dense_from_operator(hessian, dense, ok)
      !! The operator's matrix, one column per unit vector
      !!
      !! This is the only honest way to interrogate a matrix-free operator: ask
      !! it what it does to every basis vector.
      type(ov_hessian_t), intent(inout) :: hessian
      real(dp), allocatable, intent(out) :: dense(:, :)
      logical, intent(out) :: ok

      real(dp), allocatable :: unit(:), column(:)
      integer :: n, j

      n = hessian%length()
      allocate (dense(n, n), unit(n), column(n))
      ok = .true.
      do j = 1, n
         unit = 0.0_dp
         unit(j) = 1.0_dp
         call hessian%apply(unit, column)
         if (hessian%error%has_error()) then
            ok = .false.
            return
         end if
         dense(:, j) = column
      end do
   end subroutine dense_from_operator

   function lowest_eigenvalue(matrix, ok) result(low)
      !! The smallest eigenvalue of a symmetric matrix, densely
      real(dp), intent(in) :: matrix(:, :)
      logical, intent(out) :: ok
      real(dp) :: low

      real(dp), allocatable :: work(:, :), values(:)
      integer(default_int) :: info

      allocate (work, source=matrix)
      allocate (values(size(matrix, 1)))
      call pic_syev(work, values, jobz="N", uplo="U", info=info)
      ok = info == 0
      low = values(1)
   end function lowest_eigenvalue

   ! ---- the cases ---------------------------------------------------------

   subroutine test_matches_dense_hessian(error)
      !! The operator is `(A+B)`, and `build_hessian` is the second opinion
      !!
      !! Nothing in this case needs the optional backend: the question is what
      !! the operator computes, not what diagonalises it, so it runs on every
      !! build with CPU integrals.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(nuclear_response_t), target :: response
      type(ov_hessian_t) :: hessian
      real(dp), allocatable :: mine(:, :), aplus(:, :), aminus(:, :)
      real(dp), allocatable :: eri(:, :, :, :), zero_h(:, :), bounds(:, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      integer :: n_ao, n_mo, n_occ, n_vir, i, a
      logical :: ok
      real(dp) :: worst

      call water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF failed")
      if (allocated(error)) return

      n_ao = mol%nao
      n_mo = size(scf%orbitals, 2)
      n_occ = scf%n_occupied
      n_vir = n_mo - n_occ

      call build_scf_ov_hessian(mol, scf%orbitals, scf%orbital_energies, n_occ, &
                                response, hessian, err)
      call check(error,.not. err%has_error(), "the Hessian would not build: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, hessian%length() == n_vir*n_occ, "the rotation space is the "// &
                 "wrong size, so the redundant rotations were not excluded")
      if (allocated(error)) return

      call dense_from_operator(hessian, mine, ok)
      call check(error, ok, "applying the Hessian failed: "//hessian%error%get_message())
      if (allocated(error)) return

      ! The independent route: `(A+B)` from the stored integral tensor.
      allocate (c_occ(n_ao, n_occ), c_vir(n_ao, n_vir), gaps(n_vir, n_occ))
      c_occ = scf%orbitals(:, 1:n_occ)
      c_vir = scf%orbitals(:, n_occ + 1:n_mo)
      do i = 1, n_occ
         do a = 1, n_vir
            gaps(a, i) = scf%orbital_energies(n_occ + a) - scf%orbital_energies(i)
         end do
      end do
      call mol%eris(eri)
      allocate (zero_h(n_ao, n_ao), bounds(0, 0))
      zero_h = 0.0_dp
      call build_hessian(mol, .false., eri, bounds, zero_h, c_occ, c_vir, gaps, &
                         aplus, aminus, 8, err)
      call check(error,.not. err%has_error(), "the dense build failed: "//err%get_message())
      if (allocated(error)) return

      ! One route is a direct Fock build per column and the other a contraction
      ! over a stored tensor, so the two differ by the direct build's screening
      ! and its threaded reduction order and not by anything structural.
      worst = maxval(abs(mine - aplus))
      call check(error, worst < 1.0e-9_dp, "the electronic Hessian recovered from the "// &
                 "response operator is not the one built densely; worst element "// &
                 "differs by "//real_to_text(worst))
      if (allocated(error)) return

      ! And it is symmetric, which the recovery could break on its own even if
      ! it agreed on the diagonal.
      worst = maxval(abs(mine - transpose(mine)))
      call check(error, worst < 1.0e-10_dp, "the recovered Hessian is not symmetric; "// &
                 "worst asymmetry "//real_to_text(worst))

      call mol%destroy()
   end subroutine test_matches_dense_hessian

   subroutine test_refuses_without_operator(error)
      !! An electronic Hessian with nothing behind it says so
      type(error_type), allocatable, intent(out) :: error

      type(ov_hessian_t) :: hessian
      real(dp) :: x(TOY_OV), hx(TOY_OV)

      hessian%n_occ = TOY_OCC
      hessian%n_vir = TOY_VIR
      hessian%n_mo = TOY_MO
      allocate (hessian%gaps(TOY_OV))
      hessian%gaps = 1.0_dp

      x = 1.0_dp
      call hessian%apply(x, hx)
      call check(error, hessian%error%has_error(), "an unattached Hessian should "// &
                 "refuse rather than return zeros")
      if (allocated(error)) return
      call check(error, all(hx == 0.0_dp), "a refused application should not leave "// &
                 "a half-computed image behind")
   end subroutine test_refuses_without_operator

   ! ---- the native eigensolver -------------------------------------------

   subroutine test_native_lowest_eigenvalue(error)
      !! The native analysis against a dense diagonalisation of the same matrix
      !!
      !! The check that depends on neither eigensolver being right: the matrix
      !! is interrogated column by column through the operator and handed to
      !! `pic_syev`, and the number that comes back is the reference. Needs no
      !! optional dependency, so it runs on every build with CPU integrals.
      type(error_type), allocatable, intent(out) :: error

      type(toy_response_t), target :: operator
      type(ov_hessian_t) :: hessian
      type(stability_result_t) :: result
      type(error_t) :: err
      real(dp), allocatable :: dense(:, :), image(:)
      real(dp) :: dense_low, residual
      logical :: ok

      call toy_hessian(operator, hessian, -2.0_dp)
      call dense_from_operator(hessian, dense, ok)
      call check(error, ok, "applying the synthetic Hessian failed")
      if (allocated(error)) return
      dense_low = lowest_eigenvalue(dense, ok)
      call check(error, ok, "the dense diagonalisation failed")
      if (allocated(error)) return
      call check(error, dense_low < -1.0e-2_dp, "the synthetic Hessian was meant to "// &
                 "have a clearly negative eigenvalue and does not, so this case is "// &
                 "not testing what it says")
      if (allocated(error)) return

      call native_stability_of_hessian(hessian, result, err, conv_tol=1.0e-10_dp)
      call check(error,.not. err%has_error(), "the native stability analysis failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, result%ran, "the native analysis reported that it did not run")
      if (allocated(error)) return
      call check(error,.not. result%stable, "a Hessian with eigenvalue "// &
                 real_to_text(dense_low)//" was called stable")
      if (allocated(error)) return
      call check(error, abs(result%lowest_curvature - dense_low) < 1.0e-8_dp, &
                 "the native lowest curvature "// &
                 real_to_text(result%lowest_curvature)//" disagrees with the dense "// &
                 "diagonalisation's "//real_to_text(dense_low))
      if (allocated(error)) return
      call check(error, result%n_products > 0, "the native analysis reported no "// &
                 "Hessian-vector products")
      if (allocated(error)) return

      call check(error, allocated(result%rotation), "an unstable reference should "// &
                 "come back with a rotation that lowers the energy")
      if (allocated(error)) return
      allocate (image(size(result%rotation)))
      call hessian%apply(result%rotation, image)
      residual = sqrt(sum((image - dense_low*result%rotation)**2))
      call check(error, residual < 1.0e-7_dp, "the returned rotation is not an "// &
                 "eigenvector of the Hessian; residual norm "//real_to_text(residual))
   end subroutine test_native_lowest_eigenvalue

   subroutine test_native_reports_stable_curvature(error)
      !! A minimum comes back with its eigenvalue, not just a verdict
      !!
      !! This is the whole difference from the borrowed path, which can only
      !! recover the eigenvalue for an unstable reference because the library
      !! returns the eigenvector and not the eigenvalue. How stiff a minimum is
      !! is as much of an answer as which way a saddle falls, so the native
      !! path reports it always -- and a test says so, because a field that is
      !! merely usually filled is a field nobody can rely on.
      type(error_type), allocatable, intent(out) :: error

      type(toy_response_t), target :: operator
      type(ov_hessian_t) :: hessian
      type(stability_result_t) :: result
      type(error_t) :: err
      real(dp), allocatable :: dense(:, :)
      real(dp) :: dense_low
      logical :: ok

      call toy_hessian(operator, hessian, 0.1_dp)
      call dense_from_operator(hessian, dense, ok)
      call check(error, ok, "applying the synthetic Hessian failed")
      if (allocated(error)) return
      dense_low = lowest_eigenvalue(dense, ok)
      call check(error, ok, "the dense diagonalisation failed")
      if (allocated(error)) return
      call check(error, dense_low > 0.0_dp, "the synthetic Hessian was meant to be "// &
                 "positive definite and is not")
      if (allocated(error)) return

      call native_stability_of_hessian(hessian, result, err, conv_tol=1.0e-10_dp)
      call check(error,.not. err%has_error(), "the native stability analysis failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, result%stable, "a positive definite Hessian was called unstable")
      if (allocated(error)) return
      call check(error, result%has_curvature, "the native analysis must report the "// &
                 "eigenvalue for a stable reference too; that gap is why it exists")
      if (allocated(error)) return
      call check(error, abs(result%lowest_curvature - dense_low) < 1.0e-8_dp, &
                 "the native lowest curvature "// &
                 real_to_text(result%lowest_curvature)//" disagrees with the dense "// &
                 "diagonalisation's "//real_to_text(dense_low))
      if (allocated(error)) return
      call check(error,.not. allocated(result%rotation), "there is no downhill "// &
                 "direction from a minimum, so none should be returned")
      if (allocated(error)) return
      call check(error, result%n_parameters == TOY_OV, "the wrong number of rotations "// &
                 "was searched")
   end subroutine test_native_reports_stable_curvature

   subroutine test_native_water(error)
      !! A well-behaved closed shell, end to end, with no optional dependency
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(nuclear_response_t), target :: response
      type(ov_hessian_t) :: hessian
      type(stability_result_t) :: result
      real(dp), allocatable :: dense(:, :)
      real(dp) :: dense_low
      logical :: ok

      call water(mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "the reference SCF failed")
      if (allocated(error)) return

      call native_scf_stability(mol, scf%orbitals, scf%orbital_energies, &
                                scf%n_occupied, result, err, conv_tol=1.0e-9_dp)
      call check(error,.not. err%has_error(), "the native stability analysis failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, result%stable, "water in STO-3G should be a minimum")
      if (allocated(error)) return
      call check(error, result%has_curvature, "the eigenvalue should come back for a "// &
                 "stable reference")
      if (allocated(error)) return

      ! The verdict, and then the reason it is the right verdict.
      call build_scf_ov_hessian(mol, scf%orbitals, scf%orbital_energies, &
                                scf%n_occupied, response, hessian, err)
      call check(error,.not. err%has_error(), "the Hessian would not build")
      if (allocated(error)) return
      call dense_from_operator(hessian, dense, ok)
      call check(error, ok, "applying the Hessian failed")
      if (allocated(error)) return
      dense_low = lowest_eigenvalue(dense, ok)
      call check(error, ok, "the dense diagonalisation failed")
      if (allocated(error)) return
      call check(error, dense_low > 0.0_dp, "the verdict agreed with the eigensolver "// &
                 "but the matrix has a negative eigenvalue "//real_to_text(dense_low))
      if (allocated(error)) return
      ! Printed, not only compared. The tolerance says the two agree; the
      ! numbers say by how much, which is what anyone asking "how well does the
      ! matrix-free solver do on a real molecule" actually wants, and it costs
      ! one line to not have to instrument the test to find out.
      call logger%info("  water/STO-3G lowest curvature: native "// &
                       real_to_text(result%lowest_curvature)//", dense "// &
                       real_to_text(dense_low)//", difference "// &
                       real_to_text(result%lowest_curvature - dense_low))
      call check(error, abs(result%lowest_curvature - dense_low) < 1.0e-7_dp, &
                 "the native curvature "//real_to_text(result%lowest_curvature)// &
                 " disagrees with the dense diagonalisation's "//real_to_text(dense_low))
   end subroutine test_native_water

   ! ---- shared setup ------------------------------------------------------

   subroutine water(mol, scf, err)
      !! Water in STO-3G, the same geometry the coupled-perturbed tests use
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      real(dp) :: c(3, 3)

      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], c, "sto-3g", mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, 10, 200, 1.0e-14_dp, 1.0e-10_dp, .false., scf, err, &
                       in_core=.true.)
   end subroutine water

   function real_to_text(value) result(text)
      !! A number in a failure message, so the message says how badly
      real(dp), intent(in) :: value
      character(len=:), allocatable :: text

      character(len=32) :: buffer

      write (buffer, "(es22.14)") value
      text = trim(adjustl(buffer))
   end function real_to_text

end module test_mqc_czt_stability

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_stability, only: collect_mqc_czt_stability_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_stability", collect_mqc_czt_stability_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
