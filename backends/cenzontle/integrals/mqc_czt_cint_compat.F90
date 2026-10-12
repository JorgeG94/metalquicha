!! Integrals libfint's `libcint_fortran` has and the libcint fork's does not
module mqc_czt_cint_compat
   !! The multipole, point-charge-field and second-derivative integrals under their
   !! `libcint_fortran` names, whichever integral library was linked.
   !!
   !! libfint's `libcint_fortran` exports these as Fortran procedures. The libcint
   !! fork's does not, so a libcint build declares them here against libcint's own
   !! CINT2 C entry points (`cint1e_<X>_sph`, `cint2e_<X>_sph`, ...) under the same
   !! Fortran names. A caller `use`s this module and cannot tell which it got.
   !!
   !! Signatures are libfint's: `ret = f(buf, shls, atm, natm, bas, nbas, env)`, and
   !! with a trailing `opt` for the two-electron ones, which callers always pass
   !! (`c_null_ptr` is a valid optimizer). `buf` is `(di*dj*ncomp)` with the
   !! component slowest; `ret` is zero when the block was screened away.
   !! `libcint_1e_grids_*` takes `shls(4) = [ish, jsh, first_grid, end_grid]`, all
   !! 0-based, and reads the points from `env` at the 0-based offset stored in
   !! `env(LIBCINT_PTR_GRIDS + 1)`.
   !!
   !! Every `LIBCINT_PTR_*` slot is a 0-based offset into `env`, so an assignment to
   !! one carries a `+ 1`.
   !!
   !! **The assumed-size dummies in the libcint branch are the point, not an
   !! oversight.** A `bind(C)` dummy has to match what C passes, which is a bare
   !! pointer; an assumed-shape argument would pass a descriptor and corrupt the
   !! call. A wrong declaration there is a silent stack corruption at run time, not
   !! a compile error: nothing on the far side of a `bind(C)` boundary checks it.
   !
   ! When libcint support is dropped, delete everything from `#else` to `#endif`
   ! and the `#ifdef`; callers can then `use libcint_fortran` directly and this
   ! module goes away.
#ifdef MQC_WITH_LIBFINT
   use libcint_fortran, only: &
      LIBCINT_PTR_GRIDS, &
      libcint_1e_r_cart, libcint_1e_r_sph, &
      libcint_1e_rr_cart, libcint_1e_rr_sph, &
      libcint_1e_rrr_cart, libcint_1e_rrr_sph, &
      libcint_1e_irp_cart, libcint_1e_irp_sph, &
      libcint_1e_drinv_cart, libcint_1e_drinv_sph, &
      libcint_1e_ipipovlp_cart, libcint_1e_ipipovlp_sph, &
      libcint_1e_ipovlpip_cart, libcint_1e_ipovlpip_sph, &
      libcint_1e_ipipkin_cart, libcint_1e_ipipkin_sph, &
      libcint_1e_ipkinip_cart, libcint_1e_ipkinip_sph, &
      libcint_1e_ipipnuc_cart, libcint_1e_ipipnuc_sph, &
      libcint_1e_ipnucip_cart, libcint_1e_ipnucip_sph, &
      libcint_1e_ipiprinv_cart, libcint_1e_ipiprinv_sph, &
      libcint_1e_iprinvip_cart, libcint_1e_iprinvip_sph, &
      libcint_1e_grids_cart, libcint_1e_grids_sph, &
      libcint_2e_ipip1_cart, libcint_2e_ipip1_sph, &
      libcint_2e_ipvip1_cart, libcint_2e_ipvip1_sph, &
      libcint_2e_ip1ip2_cart, libcint_2e_ip1ip2_sph
#else
   use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_null_ptr
#endif
   implicit none
   private

   public :: LIBCINT_PTR_GRIDS
   public :: libcint_1e_r_cart, libcint_1e_r_sph
   public :: libcint_1e_rr_cart, libcint_1e_rr_sph
   public :: libcint_1e_rrr_cart, libcint_1e_rrr_sph
   public :: libcint_1e_irp_cart, libcint_1e_irp_sph
   public :: libcint_1e_drinv_cart, libcint_1e_drinv_sph
   public :: libcint_1e_ipipovlp_cart, libcint_1e_ipipovlp_sph
   public :: libcint_1e_ipovlpip_cart, libcint_1e_ipovlpip_sph
   public :: libcint_1e_ipipkin_cart, libcint_1e_ipipkin_sph
   public :: libcint_1e_ipkinip_cart, libcint_1e_ipkinip_sph
   public :: libcint_1e_ipipnuc_cart, libcint_1e_ipipnuc_sph
   public :: libcint_1e_ipnucip_cart, libcint_1e_ipnucip_sph
   public :: libcint_1e_ipiprinv_cart, libcint_1e_ipiprinv_sph
   public :: libcint_1e_iprinvip_cart, libcint_1e_iprinvip_sph
   public :: libcint_1e_grids_cart, libcint_1e_grids_sph
   public :: libcint_2e_ipip1_cart, libcint_2e_ipip1_sph
   public :: libcint_2e_ipvip1_cart, libcint_2e_ipvip1_sph
   public :: libcint_2e_ip1ip2_cart, libcint_2e_ip1ip2_sph

#ifndef MQC_WITH_LIBFINT
   integer(c_int), parameter :: LIBCINT_PTR_GRIDS = 12
      !! `PTR_GRIDS` from libcint's `cint.h`, which the fork's Fortran interface does
      !! not export. 0-based like every other `PTR_*`.

   interface

      function libcint_1e_r_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_r_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_r_cart

      function libcint_1e_r_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_r_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_r_sph

      function libcint_1e_rr_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rr_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_rr_cart

      function libcint_1e_rr_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rr_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_rr_sph

      function libcint_1e_rrr_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rrr_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_rrr_cart

      function libcint_1e_rrr_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rrr_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_rrr_sph

      function libcint_1e_irp_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_irp_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_irp_cart

      function libcint_1e_irp_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_irp_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_irp_sph

      function libcint_1e_drinv_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_drinv_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_drinv_cart

      function libcint_1e_drinv_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_drinv_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_drinv_sph

      function libcint_1e_ipipovlp_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipipovlp_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipipovlp_cart

      function libcint_1e_ipipovlp_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipipovlp_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipipovlp_sph

      function libcint_1e_ipovlpip_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipovlpip_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipovlpip_cart

      function libcint_1e_ipovlpip_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipovlpip_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipovlpip_sph

      function libcint_1e_ipipkin_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipipkin_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipipkin_cart

      function libcint_1e_ipipkin_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipipkin_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipipkin_sph

      function libcint_1e_ipkinip_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipkinip_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipkinip_cart

      function libcint_1e_ipkinip_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipkinip_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipkinip_sph

      function libcint_1e_ipipnuc_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipipnuc_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipipnuc_cart

      function libcint_1e_ipipnuc_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipipnuc_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipipnuc_sph

      function libcint_1e_ipnucip_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipnucip_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipnucip_cart

      function libcint_1e_ipnucip_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipnucip_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipnucip_sph

      function libcint_1e_ipiprinv_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipiprinv_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipiprinv_cart

      function libcint_1e_ipiprinv_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipiprinv_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_ipiprinv_sph

      function libcint_1e_iprinvip_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_iprinvip_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_iprinvip_cart

      function libcint_1e_iprinvip_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_iprinvip_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function libcint_1e_iprinvip_sph

      function libcint_2e_ipip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ipip1_cart")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function libcint_2e_ipip1_cart

      function libcint_2e_ipip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ipip1_sph")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function libcint_2e_ipip1_sph

      function libcint_2e_ipvip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ipvip1_cart")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function libcint_2e_ipvip1_cart

      function libcint_2e_ipvip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ipvip1_sph")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function libcint_2e_ipvip1_sph

      function libcint_2e_ip1ip2_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ip1ip2_cart")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function libcint_2e_ip1ip2_cart

      function libcint_2e_ip1ip2_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ip1ip2_sph")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function libcint_2e_ip1ip2_sph

      ! libcint's grid integral takes an optimizer as an eighth argument, which
      ! libfint's does not, so it is private and wrapped below.
      function cint1e_grids_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint1e_grids_sph")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function cint1e_grids_sph

      function cint1e_grids_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint1e_grids_cart")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function cint1e_grids_cart

   end interface

contains

   function libcint_1e_grids_sph(buf, shls, atm, natm, bas, nbas, env) result(ret)
      !! The grid integral in libfint's seven-argument form; libcint's optimizer is null.
      ! allow(assumed-size)
      real(c_double), intent(out) :: buf(*)
      ! allow(assumed-size)
      integer(c_int), intent(in) :: shls(4), atm(*), bas(*)
      integer(c_int), intent(in) :: natm, nbas
      ! allow(assumed-size)
      real(c_double), intent(in) :: env(*)
      integer(c_int) :: ret

      ret = cint1e_grids_sph(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
   end function libcint_1e_grids_sph

   function libcint_1e_grids_cart(buf, shls, atm, natm, bas, nbas, env) result(ret)
      !! The grid integral in libfint's seven-argument form; libcint's optimizer is null.
      ! allow(assumed-size)
      real(c_double), intent(out) :: buf(*)
      ! allow(assumed-size)
      integer(c_int), intent(in) :: shls(4), atm(*), bas(*)
      integer(c_int), intent(in) :: natm, nbas
      ! allow(assumed-size)
      real(c_double), intent(in) :: env(*)
      integer(c_int) :: ret

      ret = cint1e_grids_cart(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
   end function libcint_1e_grids_cart
#endif

end module mqc_czt_cint_compat
