!! The fitted coupled-perturbed operator, against the exact one on the same integrals
program check_df_cphf
   !! **What is being asserted, and what is not.**
   !!
   !! `response_operator_df` never assembles the response density's exchange
   !! from a density. It substitutes the factorization the trial rotation
   !! already carries, `Dt = X C_occ^T + C_occ X^T`, and contracts through the
   !! auxiliary index instead. That rearrangement is either right or wrong as
   !! algebra, and how good density fitting is has nothing to do with which.
   !!
   !! So the check here is exact rather than tolerated: form the fitted
   !! integrals in four-index form, `(uv|ls)_RI = sum_P B(uv,P) B(ls,P)`, hand
   !! *those* to the ordinary exact-ERI solver, and require the two responses to
   !! agree to rounding. Any factor, transpose or missing permutation in the
   !! factored build shows up immediately, on a case small enough to form the
   !! four-index object at all.
   !!
   !! The distance from the *exact* response is printed alongside, and is a
   !! different quantity entirely -- that one is the fitting error, and it is
   !! meant to be small rather than zero. Confusing the two is how a wrong
   !! operator gets excused as "within the fitting error".
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, &
                                build_df_tensor
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_cphf, only: cphf_solve
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none

   integer, parameter :: N_DIM = 3

   integer :: n_bad

   n_bad = 0
   write (*, "(a)") "== the fitted CPHF operator"

   call check_case("H2O / sto-3g", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "sto-3g", "def2-universal-jkfit", 10, n_bad)

   ! Asymmetric, so a transposed factor cannot hide behind symmetry.
   call check_case("HCN / sto-3g", [1, 6, 7], ["H", "C", "N"], &
                   reshape([0.0_dp, 0.0_dp, -2.0_dp, &
                            0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 2.2_dp], [N_DIM, 3]), &
                   "sto-3g", "def2-universal-jkfit", 14, n_bad)

   call check_case("H2O / cc-pvdz", [8, 1, 1], ["O", "H", "H"], &
                   reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, -1.4308_dp, 1.1078_dp, &
                            0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                   "cc-pvdz", "def2-universal-jkfit", 10, n_bad)

   write (*, "(a)") ""
   if (n_bad == 0) then
      write (*, "(a)") "all fitted CPHF checks passed"
   else
      write (*, "(a,i0,a)") "FAILED: ", n_bad, " case(s)"
      stop 1
   end if

contains

   subroutine check_case(label, numbers, symbols, coords, basis, aux_basis, nelec, n_bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      integer, intent(in) :: nelec
      integer, intent(inout) :: n_bad

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: error
      real(dp), allocatable :: bmat(:, :), eri_ri(:, :, :, :), eri(:, :, :, :)
      real(dp), allocatable :: rhs(:, :, :)
      real(dp), allocatable :: u_fitted(:, :, :), u_refitted(:, :, :), u_exact(:, :, :)
      integer :: n_ao, n_occ, n_vir, naux, u, v, l, s, p, a, i
      real(dp) :: algebra, fitting

      write (*, "(a)") ""
      write (*, "(a,a)") "== ", label
      flush (output_unit)

      call build_czt_molecule(numbers, symbols, coords, basis, mol, error)
      if (fail(error, n_bad)) return
      call build_czt_molecule(numbers, symbols, coords, aux_basis, aux, error)
      if (fail(error, n_bad)) return
      call run_czt_rhf(mol, nelec, 300, 1.0e-14_dp, 1.0e-12_dp, .false., scf, error)
      if (fail(error, n_bad)) return

      n_ao = mol%nao
      n_occ = nelec/2

      call build_df_tensor(mol, aux, bmat, error)
      if (fail(error, n_bad)) return
      naux = size(bmat, 2)

      ! One right-hand side is enough -- the operator is what is being tested,
      ! not the perturbation -- and a fixed pattern in the occupied-virtual
      ! block is both reproducible and free of any symmetry that could let a
      ! transposed term cancel.
      n_vir = size(scf%orbitals, 2) - n_occ
      allocate (rhs(n_vir, n_occ, 1))
      do i = 1, n_occ
         do a = 1, n_vir
            rhs(a, i, 1) = 1.0_dp/real(a + 2*i, dp)
         end do
      end do

      ! The fitted integrals, in four-index form. Only ever affordable because
      ! these cases are small -- which is the point of the case list.
      allocate (eri_ri(n_ao, n_ao, n_ao, n_ao))
      eri_ri = 0.0_dp
      do p = 1, naux
         do s = 1, n_ao
            do l = 1, n_ao
               do v = 1, n_ao
                  do u = 1, n_ao
                     eri_ri(u, v, l, s) = eri_ri(u, v, l, s) &
                                          + bmat(u + (v - 1)*n_ao, p) &
                                          *bmat(l + (s - 1)*n_ao, p)
                  end do
               end do
            end do
         end do
      end do

      call cphf_solve(mol, scf%orbitals, scf%orbital_energies, n_occ, response=u_fitted, &
                      error=error, mo_rhs=rhs, tol=1.0e-12_dp, max_iter=200, bmat=bmat)
      if (fail(error, n_bad)) return

      call cphf_solve(mol, scf%orbitals, scf%orbital_energies, n_occ, &
                      response=u_refitted, error=error, mo_rhs=rhs, tol=1.0e-12_dp, &
                      max_iter=200, eri_in=eri_ri)
      if (fail(error, n_bad)) return

      call mol%eris(eri)
      call cphf_solve(mol, scf%orbitals, scf%orbital_energies, n_occ, response=u_exact, &
                      error=error, mo_rhs=rhs, tol=1.0e-12_dp, max_iter=200, eri_in=eri)
      if (fail(error, n_bad)) return

      algebra = maxval(abs(u_fitted - u_refitted))
      fitting = maxval(abs(u_fitted - u_exact))
      write (*, "(a,i0,a,i0)") "  n_ao = ", n_ao, "   naux = ", naux
      write (*, "(a,es14.4)") "  factored against four-index, same integrals: ", algebra
      write (*, "(a,es14.4)") "  fitted against exact (the fitting error):    ", fitting
      flush (output_unit)

      ! The first is an identity and is held to rounding. The second is an
      ! approximation and is only sanity-checked: a JKFIT set that moved the
      ! response by more than a percent would mean the tensor is not what it
      ! claims to be, which is a different failure from a wrong operator.
      if (algebra > 1.0e-9_dp) then
         write (*, "(a)") "  FAIL: the factored build is not the fitted operator"
         n_bad = n_bad + 1
      end if
      if (fitting > 1.0e-2_dp*maxval(abs(u_exact))) then
         write (*, "(a)") "  FAIL: the fitting error is too large to be fitting error"
         n_bad = n_bad + 1
      end if

      call mol%destroy()
      call aux%destroy()
   end subroutine check_case

   function fail(error, n_bad) result(bad)
      type(error_t), intent(inout) :: error
      integer, intent(inout) :: n_bad
      logical :: bad

      bad = error%has_error()
      if (bad) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
      end if
   end function fail

end program check_df_cphf
