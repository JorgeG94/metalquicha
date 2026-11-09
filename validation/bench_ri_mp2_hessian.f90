!! Wall clock and peak memory of the RI-MP2 correlation Hessian against the conventional one
program bench_ri_mp2_hessian
   !! One method on one case per process, so that the peak resident set a run
   !! reports is that method's and not the larger of the two's:
   !!
   !!     bench_ri_mp2_hessian <case> <ri|conventional|rhf|partial>
   !!
   !! The SCF is converged first and its time and resident set are printed, so the
   !! Hessian's own figures are the differences. The Hessian is the correlation
   !! block only, in both, so the two answer the same question; `|H|` is printed
   !! to show they agree to the fitting error. Frozen core is off: the
   !! conventional Hessian runs the same case either way and the all-electron
   !! one is the stricter.
   use pic_types, only: dp, int64
   use pic_logger, only: logger => global_logger, performance_level
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_hessian, only: rhf_hessian, partial_hessian
   use mqc_czt_mp2_hessian, only: mp2_correlation_hessian
   use mqc_czt_ri_mp2_hessian, only: ri_mp2_correlation_hessian
   use omp_lib, only: omp_get_max_threads
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none

   character(len=64) :: which, method
   type(czt_molecule_t) :: mol, aux
   type(rhf_result_t) :: scf
   type(error_t) :: error
   real(dp), allocatable :: hess(:, :, :, :), hess_ref(:, :, :, :)
   real(dp), allocatable :: coords(:, :)
   integer, allocatable :: numbers(:)
   character(len=2), allocatable :: symbols(:)
   character(len=16) :: basis, aux_basis
   integer(int64) :: t0, t1, t2, rate
   real(dp) :: rss0, rss1
   integer :: nmol, i, nelec

   ! The stage timers of the fitted Hessian are logged at this level.
   call logger%configure(level=performance_level)
   call get_command_argument(1, which)
   call get_command_argument(2, method)
   select case (trim(which))
   case ("h2o-dz")
      nmol = 1
      basis = "cc-pvdz"
      aux_basis = "cc-pvdz-rifit"
   case ("h2o2-dz")
      nmol = 2
      basis = "cc-pvdz"
      aux_basis = "cc-pvdz-rifit"
   case ("h2o3-dz")
      nmol = 3
      basis = "cc-pvdz"
      aux_basis = "cc-pvdz-rifit"
   case ("h2o-tz")
      nmol = 1
      basis = "cc-pvtz"
      aux_basis = "cc-pvtz-rifit"
   case default
      write (*, "(a)") "usage: bench_ri_mp2_hessian h2o-dz|h2o2-dz|h2o3-dz|h2o-tz ri|conventional|rhf|partial"
      stop 1
   end select

   allocate (numbers(3*nmol), symbols(3*nmol), coords(3, 3*nmol))
   do i = 1, nmol
      numbers(3*i - 2:3*i) = [8, 1, 1]
      symbols(3*i - 2:3*i) = ["O ", "H ", "H "]
      coords(:, 3*i - 2) = [0.0_dp, 0.0_dp, 5.6_dp*(i - 1)]
      coords(:, 3*i - 1) = [0.0_dp, -1.4308_dp, 1.1078_dp + 5.6_dp*(i - 1)]
      coords(:, 3*i) = [0.0_dp, 1.4308_dp, 1.1078_dp + 5.6_dp*(i - 1)]
      ! Tilt each further molecule so no symmetry hides a wrong block.
      if (i > 1) coords(1, 3*i - 2:3*i) = coords(1, 3*i - 2:3*i) + 0.3_dp*(i - 1)
   end do
   nelec = 10*nmol

   call system_clock(t0, rate)
   call build_czt_molecule(numbers, symbols, coords, trim(basis), mol, error)
   call build_czt_molecule(numbers, symbols, coords, trim(aux_basis), aux, error)
   call run_czt_rhf(mol, nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, error)
   if (error%has_error()) then
      write (*, "(a,a)") "FAIL: ", error%get_message()
      stop 1
   end if
   call system_clock(t1)
   rss0 = peak_rss_gb()
   write (*, "(a,i0,a,i0,a,i0,a,i0,a)") "== ", nmol, " x H2O  nao=", mol%nao, "  naux=", aux%nao, &
      "  3N=", 9*nmol, "  threads="//trim(thread_text())
   write (*, "(a,f9.2,a,f8.3,a)") "   SCF and setup  ", real(t1 - t0, dp)/real(rate, dp), &
      " s   peak RSS ", rss0, " GB"

   select case (trim(method))
   case ("ri")
      call ri_mp2_correlation_hessian(mol, aux, scf%orbitals, scf%orbital_energies, nelec/2, &
                                      hess, error)
   case ("rhf")
      ! The reference Hessian alone: what the two-electron second-derivative
      ! sweep costs, for comparing builds of that kernel.
      call rhf_hessian(mol, numbers, scf%density, scf%orbitals, scf%orbital_energies, &
                       nelec/2, hess, error)
   case ("partial")
      ! The fixed-density second derivatives only, no response: the sweep over
      ! shell quartets that `hess_2e_contract` is, for timing that kernel alone.
      call partial_hessian(mol, scf%density, scf%density, hess, error)
   case ("conventional")
      call mp2_correlation_hessian(mol, scf%orbitals, scf%orbital_energies, scf%density, &
                                   nelec/2, 0, hess, hess_ref, error)
   case default
      write (*, "(a)") "the method is ri, conventional, rhf or partial"
      stop 1
   end select
   call system_clock(t2)
   if (error%has_error()) then
      write (*, "(a,a)") "FAIL: ", error%get_message()
      stop 1
   end if
   rss1 = peak_rss_gb()
   write (*, "(a,a,f9.2,a,f8.3,a,es14.6)") "   ", trim(method)//" Hessian  ", &
      real(t2 - t1, dp)/real(rate, dp), " s   peak RSS ", rss1, " GB   |H_corr| = ", &
      sqrt(sum(hess**2))
   flush (output_unit)

contains

   function peak_rss_gb() result(gb)
      !! VmHWM from /proc/self/status: the high-water mark of the resident set
      real(dp) :: gb
      character(len=256) :: line
      integer :: unit, stat
      real(dp) :: kb

      gb = -1.0_dp
      open (newunit=unit, file="/proc/self/status", status="old", action="read", iostat=stat)
      if (stat /= 0) return
      do
         read (unit, "(a)", iostat=stat) line
         if (stat /= 0) exit
         if (line(1:6) == "VmHWM:") then
            read (line(7:), *, iostat=stat) kb
            if (stat == 0) gb = kb/1024.0_dp/1024.0_dp
            exit
         end if
      end do
      close (unit)
   end function peak_rss_gb

   function thread_text() result(text)
      character(len=16) :: text

      write (text, "(i0)") omp_get_max_threads()
   end function thread_text

end program bench_ri_mp2_hessian
