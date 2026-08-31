!! The SCF dipole derivative, against a difference of the SCF dipole
module test_mqc_scf_dipole_deriv
   !! Unit 2 of the infrared ladder: the whole dipole derivative
   !!
   !! `test_mqc_dipole_deriv` pins the integrals one rung down. This is the
   !! quantity an intensity is actually built from, and it differs from those
   !! integrals by the term that is easy to leave out and impossible to notice:
   !! the electrons relax when a nucleus moves, so
   !!
   !!     d mu / dR = (nuclear) + (basis functions ride the nucleus)
   !!                           + (density responds)
   !!
   !! and the third term comes from the coupled-perturbed solve. Without it the
   !! derivative is smooth, has the right sign, obeys the translational sum
   !! rule, and is wrong by tens of percent -- so it is checked here against a
   !! difference of the converged dipole rather than against a symmetry.
   !!
   !! Hartree-Fock and Kohn-Sham both, because the response operator differs
   !! between them and a `k_scale` or a kernel dropped on one side of it would
   !! pass the first and fail the second.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_hessian, only: response_hessian
   use mqc_libcint_multipole, only: multipole_matrices
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_available
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   implicit none
   private

   public :: collect_mqc_scf_dipole_deriv_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.000000_dp, 0.000000_dp, &
                           0.0_dp, 0.000000_dp, 1.814137_dp, &
                           0.0_dp, 1.756000_dp, -0.454300_dp], [3, 3])
   integer, parameter :: WATER_NELEC = 10
   character(len=*), parameter :: BASIS = "6-31g"
   !> The 7-point O(h^6) stencil the rest of this repository's derivative
   !> checks use, and for the reason they use it. A plain central difference
   !> was tried here first and has a floor near 1e-7 on this quantity: 2.22e-07
   !> at h = 1e-3, 1.19e-07 at 5e-4, then *back up* to 7.67e-07 at 2.5e-4 as
   !> round-off overtakes truncation. That floor is thick enough to hide a real
   !> error, so it is not a tolerance worth setting.
   !> One column, as `test_mqc_mp2_hessian_fd` takes one: six converged
   !> gradients prove the derivation as well as fifty-four would, and the full
   !> matrix is covered against PySCF without displacing anything.
   integer, parameter :: PERT_ATOM = 1
   integer, parameter :: PERT_CART = 3

   real(dp), parameter :: STEP = 2.0e-3_dp

   !> PySCF's own analytic `d mu_a / dR_(X,b)`, `(3, 9)` with the Cartesian
   !> index fastest, water/6-31G about the nuclear centroid.
   !>
   !> **Assembled from PySCF primitives rather than taken from a PySCF
   !> routine**, because there is not one: `pyscf.prop.infrared` is a separate
   !> package and is not installed here. What is used is PySCF's `int1e_irp`
   !> through libcint and its coupled-perturbed solve in `pyscf.hessian`, put
   !> together the way this repository puts its own together.
   !>
   !> So this is independent in *code* and not in *derivation* -- a wrong
   !> formula would be wrong identically on both sides. That half is what
   !> `column_differences_the_dipole` is for: it knows nothing about the
   !> formula, only about the converged dipole. The two together are the pair
   !> worth having, one saying the physics is right and the other saying the
   !> implementation is, at a precision no finite difference reaches.
   real(dp), parameter :: PYSCF_HF(3, 9) = transpose(reshape([ &
      -9.322692995854e-01_dp, -1.897276338758e-16_dp, -8.157149066743e-17_dp, 4.661036045585e-01_dp, 2.711586012525e-17_dp, 6.469395321244e-17_dp, 4.661656950269e-01_dp, 1.626117737506e-16_dp, 1.687753745499e-17_dp, &
      2.164412485395e-16_dp, -4.252399592844e-01_dp, 5.401963535989e-02_dp, -5.366410207491e-17_dp, 3.422515533543e-01_dp, -2.520011370934e-02_dp, -1.627771464646e-16_dp, 8.298840593013e-02_dp, -2.881952165054e-02_dp, &
      -6.507053263866e-16_dp, 5.413253152013e-02_dp, -4.530989731287e-01_dp, 5.937490087651e-16_dp, -9.596935101097e-02_dp, 9.690371397993e-02_dp, 5.695631762138e-17_dp, 4.183681949083e-02_dp, 3.561952591488e-01_dp], [9, 3]))

   !> The same at B3LYP, **on a level-7 grid**, and the level is the point.
   !> Against PySCF at level 5 this lands 1.49e-04 apart and at level 7 it is
   !> 3.73e-08 -- a factor of four thousand for two refinements, which is the
   !> two codes' quadratures differing rather than the derivative. A test at
   !> level 5 would need a 1e-3 bound and would prove nothing.
   real(dp), parameter :: PYSCF_KS(3, 9) = transpose(reshape([ &
      -8.714634688830e-01_dp, -5.772849780286e-17_dp, -2.158916310143e-16_dp, 4.356968462644e-01_dp, 4.135566907787e-17_dp, 1.655671191805e-16_dp, 4.357666103563e-01_dp, -1.288656863473e-17_dp, 5.561266005859e-17_dp, &
      -2.523098133985e-16_dp, -2.909840396750e-01_dp, 6.903060749963e-02_dp, 9.841889239141e-19_dp, 3.029711715825e-01_dp, -4.775041067811e-02_dp, 2.061008687331e-16_dp, -1.198713361976e-02_dp, -2.128019242858e-02_dp, &
      -2.388844393797e-17_dp, 6.915066873745e-02_dp, -3.266422555407e-01_dp, -1.369474381003e-16_dp, -1.028342714492e-01_dp, 5.828927547819e-03_dp, -7.969224953476e-17_dp, 3.368360691565e-02_dp, 3.208133225513e-01_dp], [9, 3]))
   integer, parameter :: KS_GRID_LEVEL = 7

   !> The same at Hartree-Fock on **6-31G***, which is where `l = 2` enters.
   !> Cartesian in our convention, so 6d. Nothing else in this file reaches a d
   !> function, and both the component layout and the assembly underneath were
   !> wrong at s and p before they were measured, so neither is taken on trust
   !> here.
   real(dp), parameter :: PYSCF_HF_D(3, 9) = transpose(reshape([ &
      -7.891501532599e-01_dp, 2.099238493154e-16_dp, 2.040924995416e-16_dp, 3.945559076098e-01_dp, -9.280084835285e-18_dp, -2.070202645103e-16_dp, 3.945942456500e-01_dp, -2.006437644801e-16_dp, 2.927764968698e-18_dp, &
      -5.755266533450e-16_dp, -4.777531502625e-01_dp, 3.207177071903e-02_dp, 5.182133519400e-17_dp, 3.195688419636e-01_dp, -1.576429940788e-02_dp, 5.237053181510e-16_dp, 1.581843082990e-01_dp, -1.630747131115e-02_dp, &
      1.483046759390e-17_dp, 3.217310415450e-02_dp, -4.942485656217e-01_dp, 1.870192737331e-16_dp, -5.811966754130e-02_dp, 1.664107959184e-01_dp, -2.018497413269e-16_dp, 2.594656338679e-02_dp, 3.278377697033e-01_dp], [9, 3]))
   character(len=*), parameter :: BASIS_D = "6-31g*"

contains

   subroutine collect_mqc_scf_dipole_deriv_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("one_column_differences_the_dipole", &
                               test_hf), &
                  new_unittest("hartree_fock_matrix_against_pyscf", &
                               test_against_pyscf_hf), &
                  new_unittest("kohn_sham_matrix_against_pyscf", &
                               test_against_pyscf_ks), &
                  new_unittest("d_functions_against_pyscf", &
                               test_against_pyscf_d) &
                  ]
   end subroutine collect_mqc_scf_dipole_deriv_tests

   !> The converged dipole at one geometry, about a fixed origin
   subroutine dipole_at(coords, functional, origin, mu, err, ok)
      real(dp), intent(in) :: coords(3, 3)
      character(len=*), intent(in) :: functional
      real(dp), intent(in) :: origin(3)
      real(dp), intent(out) :: mu(3)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: ctx
      real(dp), allocatable :: dip(:, :, :)
      integer :: ia, a

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, BASIS, mol, err)
      if (err%has_error()) return
      if (len_trim(functional) > 0) then
         call xc_context_create(mol, functional, ctx, err, level=KS_GRID_LEVEL)
         if (err%has_error()) then
            call mol%destroy()
            return
         end if
         call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, &
                              .false., scf, err, xc=ctx)
      else
         call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, &
                              .false., scf, err)
      end if
      if (.not. err%has_error()) call multipole_matrices(mol, origin, 1, dip, err)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if

      do a = 1, 3
         mu(a) = -sum(scf%density*dip(:, :, a))
      end do
      ! The nuclear half, on the ECP-reduced charge for the reason the assembly
      ! gives: the core electrons are absent from both halves or from neither.
      do ia = 1, 3
         mu = mu + mol%charges(ia)*(mol%coords(:, ia) - origin)
      end do
      call mol%destroy()
      ok = .true.
   end subroutine dipole_at

   subroutine compare(functional, worst, err, ok)
      character(len=*), intent(in) :: functional
      real(dp), intent(out) :: worst
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: ctx
      real(dp), allocatable :: hess(:, :, :, :), ddip(:, :)
      real(dp), parameter :: WEIGHT(6) = [-1.0_dp, 9.0_dp, -45.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]
      integer, parameter :: OFFSET(6) = [-3, -2, -1, 1, 2, 3]
      real(dp) :: origin(3), plus(3), fd(3), coords(3, 3)
      integer :: ia, b, a, threads, point

      ok = .false.
      worst = 0.0_dp
      threads = omp_get_max_threads()
      call omp_set_num_threads(1)

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, BASIS, mol, err)
      if (err%has_error()) then
         call omp_set_num_threads(threads)
         return
      end if
      ! The same origin the assembly picks, so the finite difference below is a
      ! difference of the same function.
      origin = 0.0_dp
      do ia = 1, 3
         origin = origin + mol%coords(:, ia)
      end do
      origin = origin/3.0_dp

      if (len_trim(functional) > 0) then
         call xc_context_create(mol, functional, ctx, err, level=KS_GRID_LEVEL)
         if (.not. err%has_error()) &
            call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, &
                                 .false., scf, err, xc=ctx)
         if (.not. err%has_error()) &
            call response_hessian(mol, scf%density, scf%orbitals, scf%orbital_energies, &
                                  WATER_NELEC/2, hess, err, xc=ctx, &
                                  reference=scf%density, k_scale=ctx%exx_fraction, &
                                  dipole_derivatives=ddip)
      else
         call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, &
                              .false., scf, err)
         if (.not. err%has_error()) &
            call response_hessian(mol, scf%density, scf%orbitals, scf%orbital_energies, &
                                  WATER_NELEC/2, hess, err, dipole_derivatives=ddip)
      end if
      call mol%destroy()
      if (err%has_error()) then
         call omp_set_num_threads(threads)
         return
      end if

      do ia = PERT_ATOM, PERT_ATOM
         do b = PERT_CART, PERT_CART
            fd = 0.0_dp
            do point = 1, 6
               coords = WATER
               coords(b, ia) = coords(b, ia) + real(OFFSET(point), dp)*STEP
               call dipole_at(coords, functional, origin, plus, err, ok)
               if (.not. ok) exit
               fd = fd + WEIGHT(point)*plus
            end do
            if (.not. ok) exit
            fd = fd/(60.0_dp*STEP)
            do a = 1, 3
               worst = max(worst, abs(ddip(a, 3*(ia - 1) + b) - fd(a)))
            end do
         end do
         if (.not. ok) exit
      end do
      call omp_set_num_threads(threads)
   end subroutine compare

   subroutine test_hf(error)
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: worst
      type(error_t) :: err
      logical :: ok

      call compare("", worst, err, ok)
      call check(error, ok, "the dipole derivative did not evaluate: "//err%get_message())
      if (allocated(error)) return
      write (*, "(a, es10.3)") "        HF    max |analytic - fd|  = ", worst
      ! Ten times the measured worst, 7.37e-08. That is the stencil's own
      ! floor and not the derivative's: it divides by 60h = 0.12 and so
      ! multiplies the converged dipole's last digits by about eight. The
      ! PySCF comparison below reaches 5.8e-08 on the same quantity without
      ! displacing anything, which is why the tight check is that one and this
      ! one is here for the derivation rather than for the precision.
      call check(error, worst < 8.0e-7_dp, &
                 "the Hartree-Fock dipole derivative does not difference the dipole")
   end subroutine test_hf

   subroutine test_against_pyscf_hf(error)
      !! The whole matrix, against PySCF's own analytic derivative
      type(error_type), allocatable, intent(out) :: error

      call against_pyscf("", BASIS, PYSCF_HF, 5.0e-7_dp, "HF   ", error)
   end subroutine test_against_pyscf_hf

   subroutine test_against_pyscf_ks(error)
      !! The same over a Kohn-Sham reference, where the response operator differs
      type(error_type), allocatable, intent(out) :: error

      if (.not. xc_available()) return
      call against_pyscf("b3lyp", BASIS, PYSCF_KS, 5.0e-7_dp, "B3LYP", error)
   end subroutine test_against_pyscf_ks

   subroutine test_against_pyscf_d(error)
      !! The whole matrix on a basis with d functions
      type(error_type), allocatable, intent(out) :: error

      call against_pyscf("", BASIS_D, PYSCF_HF_D, 5.0e-7_dp, "HF/d ", error)
   end subroutine test_against_pyscf_d

   subroutine against_pyscf(functional, basis, reference, tol, label, error)
      character(len=*), intent(in) :: functional, basis, label
      real(dp), intent(in) :: reference(3, 9)
      real(dp), intent(in) :: tol
         !! Ten times the measured worst -- 5.82e-08 for Hartree-Fock and
         !! 3.73e-08 for B3LYP. Not tighter, because the reference cannot
         !! support it: `pyscf.hessian.rhf.solve_mo1` calls `cphf.solve`
         !! without passing `tol`, so its coupled-perturbed threshold is the
         !! library default and no attribute reaches it. That floor is this
         !! repository's own finding, already recorded against the Hessian
         !! tolerances, and it is what both numbers are sitting on.
      type(error_type), allocatable, intent(out) :: error

      real(dp), allocatable :: ddip(:, :)
      real(dp) :: net(3, 3)
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, b

      call derivative(functional, basis, ddip, err, ok)
      call check(error, ok, "the dipole derivative did not evaluate: "//err%get_message())
      if (allocated(error)) return

      write (*, "(a, a, a, es10.3)") "        ", label, " max |ours - pyscf|  = ", &
         maxval(abs(ddip - reference))
      call check(error, maxval(abs(ddip - reference)) < tol, &
                 "the dipole derivative disagrees with PySCF's own")
      if (allocated(error)) return

      ! Translating a neutral molecule cannot change its dipole, so the columns
      ! sum to the total charge -- zero here. Cheap, exact, and it catches a
      ! term deposited on the wrong atom, which is the one mistake the
      ! comparison above would also catch but only if PySCF made it differently.
      net = 0.0_dp
      do ia = 1, 3
         do b = 1, 3
            do a = 1, 3
               net(a, b) = net(a, b) + ddip(a, 3*(ia - 1) + b)
            end do
         end do
      end do
      write (*, "(a, a, a, es10.3)") "        ", label, " translational sum   = ", &
         maxval(abs(net))
      call check(error, maxval(abs(net)) < sum_rule_tol(functional), &
                 "the dipole derivative does not sum to the molecular charge")
   end subroutine against_pyscf

   !> How exactly the columns are expected to sum, which differs between the two
   !>
   !> Hartree-Fock has no quadrature in it, so the sum rule is an identity among
   !> integrals: measured 8.99e-15, held at machine precision, and a term left
   !> on the wrong atom cannot hide anywhere near that.
   !>
   !> Kohn-Sham is limited by the grid instead -- 3.68e-08 at level 7, the same
   !> order as this repository's entire disagreement with PySCF on the same
   !> matrix, which is what says the two share a cause. Holding it at the
   !> Hartree-Fock bound would make it a check on the quadrature rather than on
   !> the assembly.
   pure function sum_rule_tol(functional) result(tol)
      character(len=*), intent(in) :: functional
      real(dp) :: tol

      if (len_trim(functional) > 0) then
         tol = 5.0e-7_dp
      else
         tol = 1.0e-12_dp
      end if
   end function sum_rule_tol

   !> `d mu / dR` at the reference geometry, nothing displaced
   subroutine derivative(functional, basis, ddip, err, ok)
      character(len=*), intent(in) :: functional, basis
      real(dp), allocatable, intent(out) :: ddip(:, :)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: ctx
      real(dp), allocatable :: hess(:, :, :, :)

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, basis, mol, err)
      if (err%has_error()) return
      if (len_trim(functional) > 0) then
         call xc_context_create(mol, functional, ctx, err, level=KS_GRID_LEVEL)
         if (.not. err%has_error()) &
            call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, &
                                 .false., scf, err, xc=ctx)
         if (.not. err%has_error()) &
            call response_hessian(mol, scf%density, scf%orbitals, scf%orbital_energies, &
                                  WATER_NELEC/2, hess, err, xc=ctx, &
                                  reference=scf%density, k_scale=ctx%exx_fraction, &
                                  dipole_derivatives=ddip)
      else
         call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, &
                              .false., scf, err)
         if (.not. err%has_error()) &
            call response_hessian(mol, scf%density, scf%orbitals, scf%orbital_energies, &
                                  WATER_NELEC/2, hess, err, dipole_derivatives=ddip)
      end if
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine derivative

end module test_mqc_scf_dipole_deriv

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_scf_dipole_deriv, only: collect_mqc_scf_dipole_deriv_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_scf_dipole_deriv", collect_mqc_scf_dipole_deriv_tests)]
   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
