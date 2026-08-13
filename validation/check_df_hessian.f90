!! How closely a fitted response Hessian reproduces the exact one
!!
!!     ./build/check_df_hessian
!!
!! `build_hessian` costs `n_ov` Fock builds and is the whole expense of a fragment
!! potential. `build_hessian_df` replaces that with two matrix products over an
!! auxiliary index. The question this asks is not whether it is faster -- it plainly
!! is -- but whether what it produces is the same matrix.
!!
!! **Measured against our own exact build, not against GAMESS.** Fitting introduces
!! an error of its own, and the potential already has to agree with GAMESS to about
!! 1e-9 on every energy term. Comparing a fitted potential straight to GAMESS mixes
!! the fitting error into a number that has several other contributions and gives no
!! way to attribute it. Against the exact build on the same geometry and the same
!! orbital basis, the difference is the fitting error and nothing else.
!!
!! What is reported per case:
!!
!!   * the largest element-wise difference in `(A+B)` and `(A-B)`, absolute and
!!     relative to the largest element of the exact matrix
!!   * the difference in the static polarizability trace that follows from them,
!!     which is the quantity a potential actually carries, and where cancellation
!!     between elements can leave the tensor far better than the matrices are
!!   * seconds for each build, which is the point of the exercise
!!
!! An auxiliary basis must match the orbital basis in angular form -- libcint builds
!! all three centres of a fitting integral in one form. That rules out the Cartesian
!! Pople sets here, since the fitting sets on hand are spherical, so these cases use
!! Dunning orbital bases with `cc-pVTZ-JKFIT`. A JK-fit set rather than an RI-fit one
!! because both an exchange-like `(ab|ij)` and a Coulomb-like `(ai|bj)` are fitted.
program check_df_hessian
   use pic_types, only: dp, int64
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_cphf, only: build_hessian, build_hessian_df
   use mqc_libcint_multipole, only: multipole_matrices
   use pic_blas_interfaces, only: pic_gemm
   use omp_lib, only: omp_get_max_threads
   implicit none

   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp
   integer :: failures

   failures = 0
   write (*, "(A,I0,A)") "  water, ", omp_get_max_threads(), " threads"
   write (*, "(A)") ""
   call one_case("cc-pvdz", "cc-pvtz-jkfit")
   call one_case("cc-pvtz", "cc-pvtz-jkfit")
   ! The same orbital basis fitted with the correlation set instead. The two are
   ! built for different integrals -- a JK set to fit the `(ab|ij)` a Fock build
   ! needs, an RI set to fit the `(ai|bj)` a correlation treatment needs -- and the
   ! response Hessian carries both, so which one serves it better is a measurement
   ! rather than a matter of which label reads correctly.
   call one_case("cc-pvtz", "cc-pvtz-rifit")
   ! The size where the exact build is the whole cost of a potential, fitted with a
   ! set built to serve any orbital basis rather than matched to this one -- there is
   ! no cc-pVQZ-JKFIT here, and the universal set is what a matched one improves on.
   call one_case("cc-pvqz", "def2-universal-jkfit")
   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[df-hessian] every case built both ways"
   else
      write (*, "(A,I0,A)") "[df-hessian] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine one_case(basis, aux_basis)
      character(len=*), intent(in) :: basis, aux_basis

      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: bounds(:, :), zero_h(:, :), eri(:, :, :, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :), gaps(:, :)
      real(dp), allocatable :: ap(:, :), am(:, :), ap_df(:, :), am_df(:, :)
      real(dp) :: c(3, 3), t0, t1, exact_s, df_s
      integer :: z(3), n_occ, n_vir, n_ao, n_mo, a, i
      character(len=2) :: symbols(3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])

      call build_libcint_molecule(z, symbols, c, basis, mol, err)
      if (err%has_error()) then
         write (*, "(A,A,A,A)") "  ", basis, ": ", err%get_message()
         failures = failures + 1
         return
      end if
      call build_libcint_molecule(z, symbols, c, aux_basis, aux, err)
      if (err%has_error()) then
         write (*, "(A,A,A,A)") "  ", aux_basis, ": ", err%get_message()
         failures = failures + 1
         return
      end if

      call run_libcint_rhf(mol, 10, 200, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
      if (err%has_error() .or. .not. scf%converged) then
         write (*, "(A,A,A)") "  ", basis, ": SCF failed"
         failures = failures + 1
         return
      end if

      n_ao = mol%nao
      n_mo = size(scf%orbitals, 2)
      n_occ = 5
      n_vir = n_mo - n_occ
      allocate (c_occ(n_ao, n_occ), c_vir(n_ao, n_vir), gaps(n_vir, n_occ))
      c_occ = scf%orbitals(:, 1:n_occ)
      c_vir = scf%orbitals(:, n_occ + 1:n_mo)
      do i = 1, n_occ
         do a = 1, n_vir
            gaps(a, i) = scf%orbital_energies(n_occ + a) - scf%orbital_energies(i)
         end do
      end do
      allocate (zero_h(n_ao, n_ao))
      zero_h = 0.0_dp

      write (*, "(A,A,A,A,A,I0,A,I0,A,I0)") "  ", basis, " with ", aux_basis, &
         ":  orbitals ", n_ao, ", pairs ", n_vir*n_occ, ", auxiliary ", aux%nao

      ! In core, so the exact build is as fast as it can be here and the comparison
      ! is not flattered by making the reference take the slow path.
      call mol%eris(eri)
      allocate (bounds(0, 0))
      call seconds(t0)
      call build_hessian(mol, .false., eri, bounds, zero_h, c_occ, c_vir, gaps, &
                         ap, am, 64, err)
      call seconds(t1)
      exact_s = t1 - t0
      if (err%has_error()) then
         write (*, "(A,A)") "    exact build failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call seconds(t0)
      call build_hessian_df(mol, aux, c_occ, c_vir, gaps, ap_df, am_df, err)
      call seconds(t1)
      df_s = t1 - t0
      if (err%has_error()) then
         write (*, "(A,A)") "    fitted build failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      call compare("(A+B)", ap, ap_df)
      call compare("(A-B)", am, am_df)
      call trace_check(mol, ap, am, ap_df, am_df, c_occ, c_vir)
      write (*, "(A,F9.3,A,F9.3,A,F7.1,A)") "    exact ", exact_s, " s   fitted ", &
         df_s, " s   speedup ", exact_s/max(df_s, 1.0e-9_dp), "x"
      write (*, "(A)") ""

      deallocate (c_occ, c_vir, gaps, zero_h, bounds, eri, ap, am, ap_df, am_df)
      call mol%destroy()
      call aux%destroy()
   end subroutine one_case

   subroutine compare(what, exact, fitted)
      character(len=*), intent(in) :: what
      real(dp), intent(in) :: exact(:, :), fitted(:, :)

      real(dp) :: worst, scale

      worst = maxval(abs(exact - fitted))
      scale = maxval(abs(exact))
      write (*, "(A,A,A,ES10.3,A,ES10.3,A,ES10.3)") "    ", what, &
         "  largest element ", scale, "   worst difference ", worst, &
         "   relative ", worst/max(scale, 1.0e-30_dp)
   end subroutine compare

   subroutine trace_check(mol, ap, am, ap_df, am_df, c_occ, c_vir)
      !! The static polarizability trace both ways
      !!
      !! At `nu = 0` the response reduces to `(A+B) S = -2h` with `h` the dipole in
      !! the occupied-virtual basis, and `alpha_kl = -2 sum h^k S^l`.
      !!
      !! **The perturbation has to be the physical one.** An earlier version of this
      !! used the orbital energy gaps as a stand-in source, on the reasoning that a
      !! wrong operator shows up against any right-hand side. It does, but the
      !! reported number then failed to reproduce run to run at the first decimal
      !! place, while every matrix comparison beside it was bit-identical. The gaps
      !! are not invariant under a rotation among degenerate virtual orbitals, and
      !! LAPACK is free to return a different basis for a degenerate eigenspace on
      !! different runs. The dipole is invariant, being a physical operator projected
      !! into whatever basis the diagonalization chose, so the trace it gives is a
      !! property of the molecule rather than of the run.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: ap(:, :), am(:, :), ap_df(:, :), am_df(:, :)
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :)

      real(dp), allocatable :: dip(:, :, :), h(:, :, :), work(:, :), sol(:)
      real(dp) :: exact_trace, df_trace
      type(error_t) :: err
      integer :: n_ao, n_occ, n_vir, n_ov, k

      n_ao = size(c_occ, 1)
      n_occ = size(c_occ, 2)
      n_vir = size(c_vir, 2)
      n_ov = n_vir*n_occ

      call multipole_matrices(mol, [0.0_dp, 0.0_dp, 0.0_dp], 1, dip, error=err)
      if (err%has_error()) then
         write (*, "(A,A)") "    dipole integrals failed: ", err%get_message()
         failures = failures + 1
         return
      end if

      allocate (h(n_vir, n_occ, 3), work(n_ao, n_occ))
      do k = 1, 3
         call pic_gemm(dip(:, :, k), c_occ, work)
         call pic_gemm(c_vir, work, h(:, :, k), transa="T")
      end do

      exact_trace = 0.0_dp
      df_trace = 0.0_dp
      do k = 1, 3
         call solve(ap, -2.0_dp*reshape(h(:, :, k), [n_ov]), sol)
         exact_trace = exact_trace - 2.0_dp*sum(reshape(h(:, :, k), [n_ov])*sol)
         call solve(ap_df, -2.0_dp*reshape(h(:, :, k), [n_ov]), sol)
         df_trace = df_trace - 2.0_dp*sum(reshape(h(:, :, k), [n_ov])*sol)
      end do
      write (*, "(A,F12.6,A,F12.6,A,ES10.3)") "    static trace alpha:  exact ", &
         exact_trace/3.0_dp, "   fitted ", df_trace/3.0_dp, "   relative ", &
         abs(exact_trace - df_trace)/abs(exact_trace)

      ! `(A-B)` enters only at finite frequency, so it gets the same treatment
      ! against the same physical source.
      exact_trace = 0.0_dp
      df_trace = 0.0_dp
      do k = 1, 3
         call solve(am, -2.0_dp*reshape(h(:, :, k), [n_ov]), sol)
         exact_trace = exact_trace - 2.0_dp*sum(reshape(h(:, :, k), [n_ov])*sol)
         call solve(am_df, -2.0_dp*reshape(h(:, :, k), [n_ov]), sol)
         df_trace = df_trace - 2.0_dp*sum(reshape(h(:, :, k), [n_ov])*sol)
      end do
      write (*, "(A,F12.6,A,F12.6,A,ES10.3)") "    the same through (A-B):     ", &
         exact_trace/3.0_dp, "   fitted ", df_trace/3.0_dp, "   relative ", &
         abs(exact_trace - df_trace)/abs(exact_trace)

      deallocate (dip, h, work)
   end subroutine trace_check

   subroutine solve(matrix, rhs, solution)
      use pic_lapack_interfaces, only: pic_getrf, pic_getrs
      real(dp), intent(in) :: matrix(:, :), rhs(:)
      real(dp), allocatable, intent(out) :: solution(:)

      real(dp), allocatable :: lu(:, :), work(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n, info

      n = size(matrix, 1)
      allocate (lu(n, n), ipiv(n), work(n, 1))
      lu = matrix
      work(:, 1) = rhs
      call pic_getrf(lu, ipiv, info)
      call pic_getrs(lu, ipiv, work, info=info)
      allocate (solution(n))
      solution = work(:, 1)
      deallocate (lu, ipiv, work)
   end subroutine solve

   subroutine seconds(t)
      real(dp), intent(out) :: t
      integer(int64) :: count, rate
      call system_clock(count, rate)
      t = real(count, dp)/real(rate, dp)
   end subroutine seconds

end program check_df_hessian
