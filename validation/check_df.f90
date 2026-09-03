!! Manual check that density fitting reproduces PySCF's
!!
!!     cmake -B build -DMQC_ENABLE_CZT=ON && ./build/check_df
!!
!! Density fitting is the approximation cuEST always makes -- it has no
!! conventional four-index path -- so a CPU reference that computes exact
!! ERIs can only ever agree with the GPU to the fitting error, which is
!! around 1e-6 and swamps any bug worth finding. Fitting on this side too is
!! what turns "roughly agrees" into "any disagreement is a bug".
!!
!! Two things are asserted, and the second is the sharper one. The fitted
!! energy matches PySCF's density_fit on the same auxiliary basis, and so
!! does the fitting *error* -- the gap between fitted and exact. A fitted
!! energy could come out right with a wrong metric if the errors happened to
!! cancel; the error itself agreeing to four digits could not.
program check_df
   use pic_types, only: dp
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_error, only: error_t
   implicit none
   type(czt_molecule_t) :: mol, aux
   type(rhf_result_t) :: exact, fitted
   type(error_t) :: error
   real(dp) :: coords(3, 3)
   integer :: z(3)
   character(len=8) :: names(3)
   integer :: failures
   z = [8, 1, 1]
   names = ["O", "H", "H"]
   coords(:, 1) = [0.0_dp, 0.0_dp, -0.1364652_dp]
   coords(:, 2) = [0.0_dp, 1.4304924_dp, 1.0826636_dp]
   coords(:, 3) = [0.0_dp, -1.4304924_dp, 1.0826636_dp]

   call build_czt_molecule(z, names, coords, "cc-pvdz", mol, error)
   if (error%has_error()) then
      print *, "orbital: ", error%get_message()
      stop
   end if
   call build_czt_molecule(z, names, coords, "cc-pvtz-jkfit", aux, error)
   if (error%has_error()) then
      print *, "aux: ", error%get_message()
      stop
   end if
   write (*, "(a,i0,a,i0)") "cc-pVDZ nao=", mol%nao, "   cc-pVTZ-JKFIT naux=", aux%nao

   call run_czt_rhf(mol, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., exact, error)
   call run_czt_rhf(mol, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., fitted, error, aux=aux)
   write (*, "(a,f18.10)") "  exact ERI   E = ", exact%energy
   write (*, "(a,f18.10)") "  density fit E = ", fitted%energy
   write (*, "(a,es12.4)") "  fitting error = ", fitted%energy - exact%energy

   ! PySCF 2.14, same geometry, same BSE basis data, scf.RHF().density_fit()
   failures = 0
   call expect(aux%nao == 139, "cc-pVTZ-JKFIT gives 139 auxiliary functions")
   call expect(abs(exact%energy - (-76.0220988827_dp)) < 1.0e-9_dp, &
               "exact ERI energy matches PySCF")
   call expect(abs(fitted%energy - (-76.0220956698_dp)) < 1.0e-9_dp, &
               "density-fitted energy matches PySCF")
   call expect(abs((fitted%energy - exact%energy) - 3.213e-06_dp) < 1.0e-9_dp, &
               "the fitting error itself matches PySCF, not just the total")

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[df] all ok -- fitted J and K agree with PySCF"
   else
      write (*, "(A,I0,A)") "[df] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine expect(condition, what)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: what

      if (.not. condition) then
         write (*, "(A)") "  FAILED: "//what
         failures = failures + 1
      end if
   end subroutine expect

end program check_df
