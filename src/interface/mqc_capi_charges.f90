!! Atomic partial charges over the C boundary
module mqc_capi_charges
   !! Mulliken and CHELPG charges on a system handle, from an RHF density.
   !!
   !! Same shape as [[mqc_capi_bond_orders]] and for the same reason: charges
   !! are an input to deciding what to calculate, and a caller weighing
   !! fragmentations asks for them once per trial partition.
   !!
   !! **Unlike bond orders, these cost a real SCF.** They come from
   !! `run_czt_rhf` in whatever basis the caller names, so the basis is a
   !! parameter rather than a default hidden inside, and it is the knob that
   !! decides what `compute` costs.
   !!
   !! **Which scheme to ask for.** CHELPG, unless the question is specifically
   !! about basis-function populations: a Mulliken charge moves by a factor of
   !! two between 6-31G and aug-cc-pVDZ where the CHELPG one barely moves, so
   !! anything treating a charge as a physical quantity wants the stable one.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_associated, c_f_pointer
   use pic_types, only: dp
   use mqc_capi_system, only: system_handle_t, last_message
   use mqc_capi_status, only: MQC_OK, MQC_FAIL, MQC_BAD_HANDLE
   use mqc_error, only: error_t
   use mqc_elements, only: element_number_to_symbol
   use mqc_czt_bridge, only: run_czt_charges
   implicit none
   private

   public :: mqc_system_compute_charges
   public :: mqc_system_get_charges
   public :: mqc_system_charge_on
   public :: mqc_system_has_charges
   public :: mqc_system_charge_scheme

contains

   function mqc_system_compute_charges(handle, scheme_len, scheme, basis_len, basis) &
      result(status) bind(C, name="mqc_system_compute_charges")
      !! Run one RHF over the whole system and keep its atomic charges
      !!
      !! `scheme` is "chelpg" or "mulliken"; an empty string takes chelpg, for
      !! the reason in the module header. `basis` is any basis the build
      !! carries; an empty string takes 6-31g.
      !!
      !! Closed shell only: an odd electron count is refused rather than
      !! silently paired up.
      type(c_ptr), value :: handle
      integer(c_int), value :: scheme_len
      character(kind=c_char), intent(in) :: scheme(scheme_len)
      integer(c_int), value :: basis_len
      character(kind=c_char), intent(in) :: basis(basis_len)
      integer(c_int) :: status

      type(system_handle_t), pointer :: h
      type(error_t) :: error
      real(dp), allocatable :: q(:)
      character(len=2), allocatable :: symbols(:)
      character(len=:), allocatable :: which, basis_name
      integer :: i, n

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null system handle"
         return
      end if
      call c_f_pointer(handle, h)

      if (h%geom%total_atoms <= 0) then
         last_message = "mqc_system_compute_charges: set the geometry first"
         status = MQC_FAIL
         return
      end if

      which = from_c(scheme, scheme_len)
      if (len_trim(which) == 0) which = "chelpg"
      basis_name = from_c(basis, basis_len)
      if (len_trim(basis_name) == 0) basis_name = "6-31g"

      if (which /= "chelpg" .and. which /= "mulliken") then
         last_message = "mqc_system_compute_charges: unknown scheme '"//which// &
                        "'; use 'chelpg' or 'mulliken'"
         status = MQC_FAIL
         return
      end if

      n = h%geom%total_atoms
      allocate (symbols(n))
      do i = 1, n
         symbols(i) = element_number_to_symbol(h%geom%element_numbers(i))
      end do

      ! Through the bridge, not into the backend: the odd-electron refusal and
      ! the SCF are both on the other side of it, so this reports rather than
      ! decides.
      call run_czt_charges(h%geom%element_numbers, symbols, h%geom%coordinates, &
                           basis_name, which, h%geom%charge, q, error)
      if (failed(error, status)) return

      if (allocated(h%charges)) deallocate (h%charges)
      allocate (h%charges(n), source=q)
      h%charge_scheme = which
      status = MQC_OK
   end function mqc_system_compute_charges

   function mqc_system_has_charges(handle) result(has) &
      bind(C, name="mqc_system_has_charges")
      !! Whether `compute` has run on this handle. Zero for a null handle too.
      type(c_ptr), value :: handle
      integer(c_int) :: has

      type(system_handle_t), pointer :: h

      has = 0
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, h)
      if (allocated(h%charges)) has = 1
   end function mqc_system_has_charges

   subroutine mqc_system_charge_scheme(handle, buffer_len, buffer) &
      bind(C, name="mqc_system_charge_scheme")
      !! Which scheme produced the charges currently on the handle
      !!
      !! Empty if none have been computed. Worth asking, because "the charge on
      !! atom 3" is not a well-defined number without it.
      use, intrinsic :: iso_c_binding, only: c_null_char
      type(c_ptr), value :: handle
      integer(c_int), value :: buffer_len
      character(kind=c_char), intent(inout) :: buffer(buffer_len)

      type(system_handle_t), pointer :: h
      character(len=:), allocatable :: text
      integer :: n, i

      if (buffer_len <= 0) return
      text = ""
      if (c_associated(handle)) then
         call c_f_pointer(handle, h)
         if (allocated(h%charges)) text = trim(h%charge_scheme)
      end if

      n = min(len(text), int(buffer_len) - 1)
      do i = 1, n
         buffer(i) = text(i:i)
      end do
      buffer(n + 1) = c_null_char
   end subroutine mqc_system_charge_scheme

   function mqc_system_get_charges(handle, n, out) result(status) &
      bind(C, name="mqc_system_get_charges")
      !! Copy the charges into a caller buffer of n doubles, one per atom
      type(c_ptr), value :: handle
      integer(c_int), value :: n
      real(c_double), intent(out) :: out(n)
      integer(c_int) :: status

      type(system_handle_t), pointer :: h

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null system handle"
         return
      end if
      call c_f_pointer(handle, h)

      if (.not. allocated(h%charges)) then
         last_message = "mqc_system_get_charges: compute them first"
         status = MQC_FAIL
         return
      end if
      if (n /= size(h%charges)) then
         last_message = "mqc_system_get_charges: buffer is the wrong size for this system"
         status = MQC_FAIL
         return
      end if

      out = h%charges
      status = MQC_OK
   end function mqc_system_get_charges

   function mqc_system_charge_on(handle, i) result(q) &
      bind(C, name="mqc_system_charge_on")
      !! One atom's charge, 0-based
      !!
      !! Unlike a bond order, a partial charge is legitimately negative, so
      !! there is no out-of-band value to signal failure with. A bad index or an
      !! uncomputed handle gives a quiet NaN, which propagates rather than
      !! silently reading as a neutral atom.
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      type(c_ptr), value :: handle
      integer(c_int), value :: i
      real(c_double) :: q

      type(system_handle_t), pointer :: h

      q = ieee_value(1.0_dp, ieee_quiet_nan)
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, h)
      if (.not. allocated(h%charges)) return
      if (i < 0 .or. i >= size(h%charges)) return
      q = h%charges(i + 1)
   end function mqc_system_charge_on

   function from_c(buf, n) result(s)
      !! A C character buffer as a Fortran string
      integer(c_int), intent(in) :: n
      character(kind=c_char), intent(in) :: buf(n)
      character(len=:), allocatable :: s

      integer :: i

      s = ""
      do i = 1, n
         s = s//buf(i)
      end do
      s = trim(adjustl(s))
   end function from_c

   function failed(error, status) result(bad)
      !! Turn an error_t into a C status, and say whether to stop
      type(error_t), intent(inout) :: error
      integer(c_int), intent(inout) :: status
      logical :: bad

      bad = error%has_error()
      if (bad) then
         last_message = "mqc_system_compute_charges: "//error%get_message()
         status = MQC_FAIL
      end if
   end function failed

end module mqc_capi_charges
