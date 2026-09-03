!! How long an analytic SCF gradient takes, against system size
program bench_scf_gradient
   !! Reports a warm timing. The first call in a process pays for lazy symbol
   !! resolution and first-touch page faults, which on a small molecule is most
   !! of what a single measurement would show -- 0.11 s of it, against 0.008 s
   !! of actual work for water in a minimal basis.
   use pic_types, only: dp, int64
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_gradient, only: czt_scf_gradient
   implicit none

   ! Wall clock, not cpu_time. The two-electron contraction is threaded, and
   ! cpu_time sums over threads -- it reports a *rise* when threading helps,
   ! which is the opposite of the thing being measured.

   call one("sto-3g", 1)
   call one("6-31g", 1)
   call one("cc-pvdz", 1)
   call one("cc-pvdz", 2)
   call one("cc-pvdz", 3)

contains

   subroutine one(basis, nwater)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nwater

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: error
      real(dp), allocatable :: g(:, :), coords(:, :)
      integer, allocatable :: numbers(:)
      character(len=2), allocatable :: symbols(:)
      integer(int64) :: tick0, tick1, rate
      integer :: iw, base

      allocate (numbers(3*nwater), symbols(3*nwater), coords(3, 3*nwater))
      do iw = 1, nwater
         base = 3*(iw - 1)
         numbers(base + 1:base + 3) = [8, 1, 1]
         symbols(base + 1) = "O "
         symbols(base + 2) = "H "
         symbols(base + 3) = "H "
         coords(:, base + 1) = [0.0_dp, 0.0_dp, 0.0_dp] + [6.0_dp*(iw - 1), 0.0_dp, 0.0_dp]
         coords(:, base + 2) = [0.0_dp, -1.4308_dp, 1.1078_dp] + [6.0_dp*(iw - 1), 0.0_dp, 0.0_dp]
         coords(:, base + 3) = [0.0_dp, 1.4308_dp, 1.1078_dp] + [6.0_dp*(iw - 1), 0.0_dp, 0.0_dp]
      end do

      call build_czt_molecule(numbers, symbols, coords, basis, mol, error)
      if (error%has_error()) then
         write (*, "(a,a)") "build failed: ", error%get_message()
         return
      end if
      call run_czt_rhf(mol, 10*nwater, 200, 1.0e-10_dp, 1.0e-8_dp, .false., scf, error)
      if (error%has_error()) then
         write (*, "(a,a)") "scf failed: ", error%get_message()
         return
      end if

      call czt_scf_gradient(mol, scf%density, orbitals=scf%orbitals, &
                            orbital_energies=scf%orbital_energies, &
                            n_occupied=scf%n_occupied, gradient=g, error=error)
      deallocate (g)
      call system_clock(tick0, rate)
      call czt_scf_gradient(mol, scf%density, orbitals=scf%orbitals, &
                            orbital_energies=scf%orbital_energies, &
                            n_occupied=scf%n_occupied, gradient=g, error=error)
      call system_clock(tick1)
      write (*, "(a10,a,i0,a,i4,a,i4,a,f10.4,a)") basis, "  (H2O)x", nwater, &
         "  nao=", mol%nao, "  nbas=", mol%nbas, "   gradient ", &
         real(tick1 - tick0, dp)/real(rate, dp), " s"
      deallocate (g)
   end subroutine one

end program bench_scf_gradient
