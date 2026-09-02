!! The second-derivative integrals, declared against whichever library is linked
module mqc_libcint_hess_abi
   !! `bind(C)` declarations for the twenty-four entry points an analytic
   !! Hessian is made of.
   !!
   !! **This is what lets one source file serve both backends.** libcint and
   !! libfint export these under the same names with the same signatures --
   !! libfint's `cint_c_abi` exists to make that true, and covers the second
   !! derivatives since v0.1.1 -- so a declaration here resolves against
   !! whichever library the build linked.
   !!
   !! **The CINT2 spelling, deliberately.** Both libraries offer these twice: as
   !! `int1e_ipipovlp_sph(out, dims, shls, ..., opt, cache)` and as
   !! `cint1e_ipipovlp_sph(buf, shls, ..., env)`. The second has no `dims`, no
   !! `cache`, and -- for one-electron integrals -- no optimizer, so neither
   !! library needs special handling for the scratch buffer libcint would
   !! otherwise want sized by a null-output call. The buffer comes back packed
   !! with the natural shell dimensions.
   !!
   !! Two-electron entry points take an optimizer and one-electron ones do not.
   !! That is libcint's convention; a null pointer is a valid optimizer and is
   !! what the callers pass.
   !!
   !! **Assumed-size arrays are the point, not an oversight.** A `bind(C)` dummy
   !! has to match what C passes, which is a bare pointer; an assumed-shape
   !! argument would pass a descriptor and corrupt the call.
   !!
   !! A wrong declaration here is a silent stack corruption at run time, not a
   !! compile error: there is nothing on the far side of a `bind(C)` boundary to
   !! check it against.
   use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr
   implicit none
   private

   public :: cint1e_ipipovlp_sph, cint1e_ipipovlp_cart
   public :: cint1e_ipovlpip_sph, cint1e_ipovlpip_cart
   public :: cint1e_ipipkin_sph, cint1e_ipipkin_cart
   public :: cint1e_ipkinip_sph, cint1e_ipkinip_cart
   public :: cint1e_ipipnuc_sph, cint1e_ipipnuc_cart
   public :: cint1e_ipnucip_sph, cint1e_ipnucip_cart
   public :: cint1e_ipiprinv_sph, cint1e_ipiprinv_cart
   public :: cint1e_iprinvip_sph, cint1e_iprinvip_cart
   public :: cint2e_ipip1_sph, cint2e_ipip1_cart
   public :: cint2e_ipvip1_sph, cint2e_ipvip1_cart
   public :: cint2e_ip1ip2_sph, cint2e_ip1ip2_cart
   public :: cint2e_ip1_sph, cint2e_ip1_cart

   interface

      function cint1e_ipipovlp_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipipovlp_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipipovlp_sph

      function cint1e_ipipovlp_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipipovlp_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipipovlp_cart

      function cint1e_ipovlpip_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipovlpip_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipovlpip_sph

      function cint1e_ipovlpip_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipovlpip_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipovlpip_cart

      function cint1e_ipipkin_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipipkin_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipipkin_sph

      function cint1e_ipipkin_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipipkin_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipipkin_cart

      function cint1e_ipkinip_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipkinip_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipkinip_sph

      function cint1e_ipkinip_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipkinip_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipkinip_cart

      function cint1e_ipipnuc_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipipnuc_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipipnuc_sph

      function cint1e_ipipnuc_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipipnuc_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipipnuc_cart

      function cint1e_ipnucip_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipnucip_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipnucip_sph

      function cint1e_ipnucip_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipnucip_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipnucip_cart

      function cint1e_ipiprinv_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipiprinv_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipiprinv_sph

      function cint1e_ipiprinv_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipiprinv_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_ipiprinv_cart

      function cint1e_iprinvip_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_iprinvip_sph")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_iprinvip_sph

      function cint1e_iprinvip_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_iprinvip_cart")
         import :: c_double, c_int
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: ret
      end function cint1e_iprinvip_cart

      function cint2e_ipip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ipip1_sph")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function cint2e_ipip1_sph

      function cint2e_ipip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ipip1_cart")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function cint2e_ipip1_cart

      function cint2e_ipvip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ipvip1_sph")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function cint2e_ipvip1_sph

      function cint2e_ipvip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ipvip1_cart")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function cint2e_ipvip1_cart

      function cint2e_ip1ip2_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ip1ip2_sph")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function cint2e_ip1ip2_sph

      function cint2e_ip1ip2_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ip1ip2_cart")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function cint2e_ip1ip2_cart

      function cint2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ip1_sph")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function cint2e_ip1_sph

      function cint2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ip1_cart")
         import :: c_double, c_int, c_ptr
         implicit none
         ! allow(assumed-size)
         real(c_double), intent(out) :: buf(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: shls(*)
         ! allow(assumed-size)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         ! allow(assumed-size)
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         ! allow(assumed-size)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: ret
      end function cint2e_ip1_cart

   end interface

end module mqc_libcint_hess_abi
