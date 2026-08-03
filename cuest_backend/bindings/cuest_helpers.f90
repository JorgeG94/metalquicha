! ============================================================================
!  cuest_helpers.f90 -- optional convenience layer over the raw `cuest` module.
!
!  The cuEST parameter/query API is generic: values are passed through
!  `void* + size_t` (see cuestParametersConfigure / cuestParametersQuery /
!  cuestQuery).  From Fortran that means building a C_LOC(...) of a TARGET
!  variable and passing C_SIZEOF(...).  These thin wrappers hide that idiom
!  for the common scalar types so callers can write, e.g.
!
!      ist = cuest_param_set_i64(CUEST_HANDLE_PARAMETERS, p, &
!                                CUEST_HANDLE_PARAMETERS_MAX_GAUSS_HERMITE, 20_c_int64_t)
!
!  Every wrapper simply forwards to the underlying cuEST call and returns the
!  cuestStatus_t it produced.  Nothing here allocates or owns memory.
! ============================================================================
module cuest_helpers
   use, intrinsic :: iso_c_binding
   use cuest
   implicit none
   private

   public :: cuest_param_set_i64, cuest_param_set_i32, cuest_param_set_f64
   public :: cuest_param_get_i64, cuest_param_get_i32, cuest_param_get_f64
   public :: cuest_query_i64, cuest_query_i32, cuest_query_f64
   public :: cuest_param_set_str, cuest_status_name

contains

   !> Human-readable name for a cuestStatus_t code (mirrors the C sample's
   !  cuestGetErrorEnum helper; handy in error messages).
   function cuest_status_name(st) result(name)
      integer(c_int), intent(in) :: st
      character(len=:), allocatable :: name
      select case (st)
      case (CUEST_STATUS_SUCCESS); name = "CUEST_STATUS_SUCCESS"
      case (CUEST_STATUS_EXCEPTION); name = "CUEST_STATUS_EXCEPTION"
      case (CUEST_STATUS_NULL_POINTER); name = "CUEST_STATUS_NULL_POINTER"
      case (CUEST_STATUS_INVALID_ARGUMENT); name = "CUEST_STATUS_INVALID_ARGUMENT"
      case (CUEST_STATUS_INVALID_SIZE); name = "CUEST_STATUS_INVALID_SIZE"
      case (CUEST_STATUS_INVALID_TYPE); name = "CUEST_STATUS_INVALID_TYPE"
      case (CUEST_STATUS_INVALID_PARAMETER); name = "CUEST_STATUS_INVALID_PARAMETER"
      case (CUEST_STATUS_INVALID_ATTRIBUTE); name = "CUEST_STATUS_INVALID_ATTRIBUTE"
      case (CUEST_STATUS_INVALID_HANDLE); name = "CUEST_STATUS_INVALID_HANDLE"
      case (CUEST_STATUS_UNKNOWN_ERROR); name = "CUEST_STATUS_UNKNOWN_ERROR"
      case (CUEST_STATUS_UNSUPPORTED_ARGUMENT); name = "CUEST_STATUS_UNSUPPORTED_ARGUMENT"
      case (CUEST_STATUS_UNSUPPORTED_ARCHITECTURE); name = "CUEST_STATUS_UNSUPPORTED_ARCHITECTURE"
      case (CUEST_STATUS_INVALID_PLAN); name = "CUEST_STATUS_INVALID_PLAN"
      case (CUEST_STATUS_HOME_NOT_FOUND); name = "CUEST_STATUS_HOME_NOT_FOUND"
      case default; name = "<unknown>"
      end select
   end function cuest_status_name

   !> Set a string-valued attribute (e.g. CUEST_HANDLE_PARAMETERS_JIT_CACHE_DIR).
   !  cuEST expects a char** whose target is a NUL-terminated string, of size
   !  sizeof(char*); cuEST copies the string during this call.
   integer(c_int) function cuest_param_set_str(ptype, params, attr, str) result(st)
      integer(c_int), value :: ptype, attr
      type(c_ptr), value :: params
      character(*), intent(in) :: str
      character(kind=c_char), target :: cbuf(len_trim(str) + 1)
      type(c_ptr), target :: sp
      integer :: i, n
      n = len_trim(str)
      do i = 1, n
         cbuf(i) = str(i:i)
      end do
      cbuf(n + 1) = c_null_char
      sp = c_loc(cbuf)
      st = cuestParametersConfigure(ptype, params, attr, c_loc(sp), c_sizeof(sp))
   end function cuest_param_set_str

   ! ---- configure a parameters object attribute -------------------------
   integer(c_int) function cuest_param_set_i64(ptype, params, attr, val) result(st)
      integer(c_int), value :: ptype, attr
      type(c_ptr), value :: params
      integer(c_int64_t), intent(in), target :: val
      st = cuestParametersConfigure(ptype, params, attr, c_loc(val), c_sizeof(val))
   end function cuest_param_set_i64

   integer(c_int) function cuest_param_set_i32(ptype, params, attr, val) result(st)
      integer(c_int), value :: ptype, attr
      type(c_ptr), value :: params
      integer(c_int32_t), intent(in), target :: val
      st = cuestParametersConfigure(ptype, params, attr, c_loc(val), c_sizeof(val))
   end function cuest_param_set_i32

   integer(c_int) function cuest_param_set_f64(ptype, params, attr, val) result(st)
      integer(c_int), value :: ptype, attr
      type(c_ptr), value :: params
      real(c_double), intent(in), target :: val
      st = cuestParametersConfigure(ptype, params, attr, c_loc(val), c_sizeof(val))
   end function cuest_param_set_f64

   ! ---- query a parameters object attribute -----------------------------
   integer(c_int) function cuest_param_get_i64(ptype, params, attr, val) result(st)
      integer(c_int), value :: ptype, attr
      type(c_ptr), value :: params
      integer(c_int64_t), intent(out), target :: val
      st = cuestParametersQuery(ptype, params, attr, c_loc(val), c_sizeof(val))
   end function cuest_param_get_i64

   integer(c_int) function cuest_param_get_i32(ptype, params, attr, val) result(st)
      integer(c_int), value :: ptype, attr
      type(c_ptr), value :: params
      integer(c_int32_t), intent(out), target :: val
      st = cuestParametersQuery(ptype, params, attr, c_loc(val), c_sizeof(val))
   end function cuest_param_get_i32

   integer(c_int) function cuest_param_get_f64(ptype, params, attr, val) result(st)
      integer(c_int), value :: ptype, attr
      type(c_ptr), value :: params
      real(c_double), intent(out), target :: val
      st = cuestParametersQuery(ptype, params, attr, c_loc(val), c_sizeof(val))
   end function cuest_param_get_f64

   ! ---- query an attribute from a live cuEST object (cuestQuery) ---------
   integer(c_int) function cuest_query_i64(handle, otype, object, attr, val) result(st)
      type(c_ptr), value :: handle, object
      integer(c_int), value :: otype, attr
      integer(c_int64_t), intent(out), target :: val
      st = cuestQuery(handle, otype, object, attr, c_loc(val), c_sizeof(val))
   end function cuest_query_i64

   integer(c_int) function cuest_query_i32(handle, otype, object, attr, val) result(st)
      type(c_ptr), value :: handle, object
      integer(c_int), value :: otype, attr
      integer(c_int32_t), intent(out), target :: val
      st = cuestQuery(handle, otype, object, attr, c_loc(val), c_sizeof(val))
   end function cuest_query_i32

   integer(c_int) function cuest_query_f64(handle, otype, object, attr, val) result(st)
      type(c_ptr), value :: handle, object
      integer(c_int), value :: otype, attr
      real(c_double), intent(out), target :: val
      st = cuestQuery(handle, otype, object, attr, c_loc(val), c_sizeof(val))
   end function cuest_query_f64

end module cuest_helpers
