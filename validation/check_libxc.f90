!! Manual check that libxc is linked and answering
!!
!!     cmake -B build -DMQC_ENABLE_LIBXC=ON && ./build/check_libxc
!!
!! Evaluates the uniform electron gas, which is the one case where an
!! exchange-correlation functional has a closed form to be checked against
!! rather than a tabulated number to be trusted.
!!
!! Slater exchange of a uniform density is
!!
!!     e_x = -(3/4) (3/pi)^(1/3) rho^(4/3)
!!
!! exactly, by construction -- LDA exchange *is* that expression. So if
!! LDA_X returns it, the library is linked, the module is the right one, the
!! spin convention is unpolarised as asked, and the units are Hartree per
!! unit volume. Four things, one number, no table.
!!
!! The correlation functional is checked only for being negative and small
!! against exchange, which is all that can be said without a table -- but
!! that much would still catch a functional id pointing somewhere else.
program check_libxc
   use pic_types, only: dp
   use xc_f03_lib_m, only: xc_f03_func_t, xc_f03_func_info_t, xc_f03_version, &
                           xc_f03_func_init, xc_f03_func_end, xc_f03_lda_exc, &
                           xc_f03_func_get_info, xc_f03_func_info_get_name, &
                           XC_UNPOLARIZED
   ! The functional ids live in their own module -- xc_f03_lib_m carries the
   ! interface and the spin constants, xc_f03_funcs_m the several thousand
   ! XC_* identifiers, which are generated at build time.
   use xc_f03_funcs_m, only: XC_LDA_X, XC_LDA_C_VWN
   implicit none

   integer, parameter :: NPOINT = 4
   real(dp), parameter :: PI = 3.14159265358979323846_dp

   type(xc_f03_func_t) :: func
   type(xc_f03_func_info_t) :: info
   real(dp) :: rho(NPOINT), exc(NPOINT)
   real(dp) :: expected
   integer :: i, failures

   failures = 0
   rho = [0.01_dp, 0.1_dp, 1.0_dp, 10.0_dp]

   write (*, "(A,A)") "libxc version: ", trim(version_text())

   ! ---- Slater exchange, against its closed form ---------------------------
   call xc_f03_func_init(func, XC_LDA_X, XC_UNPOLARIZED)
   info = xc_f03_func_get_info(func)
   write (*, "(A,A)") "functional:    ", trim(xc_f03_func_info_get_name(info))
   call xc_f03_lda_exc(func, int(NPOINT, 8), rho, exc)
   call xc_f03_func_end(func)

   write (*, "(A)") "  rho          e_x (libxc)       -(3/4)(3/pi)^(1/3) rho^(1/3)"
   do i = 1, NPOINT
      ! libxc returns the energy *per particle*, so the rho^(4/3) of the
      ! energy density becomes rho^(1/3) here. Getting that wrong is the
      ! obvious way to misread this call and it would be a factor of rho.
      expected = -0.75_dp*(3.0_dp/PI)**(1.0_dp/3.0_dp)*rho(i)**(1.0_dp/3.0_dp)
      write (*, "(F8.3,2F22.14)") rho(i), exc(i), expected
      call expect(abs(exc(i) - expected) < 1.0e-12_dp, "LDA exchange matches its closed form")
   end do

   ! ---- VWN correlation, for sign and scale --------------------------------
   call xc_f03_func_init(func, XC_LDA_C_VWN, XC_UNPOLARIZED)
   info = xc_f03_func_get_info(func)
   write (*, "(A,A)") "functional:    ", trim(xc_f03_func_info_get_name(info))
   call xc_f03_lda_exc(func, int(NPOINT, 8), rho, exc)
   call xc_f03_func_end(func)
   do i = 1, NPOINT
      call expect(exc(i) < 0.0_dp, "correlation energy is negative")
      call expect(abs(exc(i)) < abs(-0.75_dp*(3.0_dp/PI)**(1.0_dp/3.0_dp) &
                                    *rho(i)**(1.0_dp/3.0_dp)), &
                  "correlation is smaller than exchange")
   end do

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[libxc] all ok -- linked, callable, and LDA exchange is exact"
   else
      write (*, "(A,I0,A)") "[libxc] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   function version_text() result(text)
      character(len=32) :: text
      integer :: major, minor, micro

      call xc_f03_version(major, minor, micro)
      write (text, "(i0,'.',i0,'.',i0)") major, minor, micro
   end function version_text

   subroutine expect(condition, what)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: what

      if (.not. condition) then
         write (*, "(A)") "  FAILED: "//what
         failures = failures + 1
      end if
   end subroutine expect

end program check_libxc
