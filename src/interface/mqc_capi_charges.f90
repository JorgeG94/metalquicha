!! Atomic partial charges over the C boundary
module mqc_capi_charges
   !! Mulliken and CHELPG charges on a system handle, from an RHF density.
   !!
   !! Same shape as [[mqc_capi_bond_orders]] and for the same reason: charges
   !! are an input to deciding what to calculate, not the output of a
   !! calculation someone asked for. A caller weighing fragmentations wants to
   !! know how much charge sits near each candidate cut, and wants to ask that
   !! many times over many trial partitions -- a Python loop, not a deck.
   !!
   !! **Unlike bond orders, these cost a real SCF.** Bond orders come from xTB,
   !! which is cheap enough to run on the whole system without thinking about
   !! it. Charges here come from `run_libcint_rhf` in whatever basis the caller
   !! names, so `compute` on a large system in a large basis is a large
   !! calculation. That is the caller's decision to make, which is why the basis
   !! is a parameter and not a default hidden inside. Measured: 32 atoms in
   !! 6-31G is about twenty seconds, nearly all of it the SCF -- so the basis is
   !! the knob that matters and the choice of scheme is not.
   !!
   !! **Which scheme to ask for.** CHELPG, unless the reason for asking is
   !! specifically about basis-function populations. Measured on water between
   !! 6-31G and aug-cc-pVDZ, the Mulliken charge on oxygen moves from -0.79 to
   !! -0.30 while CHELPG moves from -0.94 to -0.74 -- the same molecule, and one
   !! of the two answers has changed by more than half. For anything downstream
   !! that treats a charge as a physical quantity, embedding especially, the
   !! basis-stable one is the only defensible input.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_associated, c_f_pointer
   use pic_types, only: dp
   use mqc_capi_system, only: system_handle_t, last_message, MQC_OK, MQC_FAIL, MQC_BAD_HANDLE
   use mqc_error, only: error_t
   use mqc_elements, only: element_number_to_symbol
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_charges, only: mulliken_charges, chelpg_charges
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   implicit none
   private

   public :: mqc_system_compute_charges
   public :: mqc_system_get_charges
   public :: mqc_system_charge_on
   public :: mqc_system_has_charges
   public :: mqc_system_charge_scheme

   integer, parameter :: SCF_MAX_ITER = 100
   real(dp), parameter :: SCF_ENERGY_TOL = 1.0e-9_dp
   real(dp), parameter :: SCF_DENSITY_TOL = 1.0e-7_dp

contains

   function mqc_system_compute_charges(handle, scheme_len, scheme, basis_len, basis) &
      result(status) bind(C, name="mqc_system_compute_charges")
      !! Run one RHF over the whole system and keep its atomic charges
      !!
      !! `scheme` is "chelpg" or "mulliken"; an empty string takes chelpg, for
      !! the reason in the module header. `basis` is any basis the build carries;
      !! an empty string takes 6-31g, which is small enough to be affordable on
      !! the systems this gets pointed at and large enough that the charges mean
      !! something.
      !!
      !! Closed shell only. An odd electron count is refused rather than
      !! silently paired up, because a caller fragmenting a radical needs to
      !! know that this cannot answer for it.
      type(c_ptr), value :: handle
      integer(c_int), value :: scheme_len
      character(kind=c_char), intent(in) :: scheme(scheme_len)
      integer(c_int), value :: basis_len
      character(kind=c_char), intent(in) :: basis(basis_len)
      integer(c_int) :: status

      type(system_handle_t), pointer :: h
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: error
      real(dp), allocatable :: overlap(:, :), q(:)
      character(len=2), allocatable :: symbols(:)
      character(len=:), allocatable :: which, basis_name
      integer :: i, n, nelec

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
      nelec = sum(h%geom%element_numbers) - h%geom%charge
      if (mod(nelec, 2) /= 0) then
         last_message = "mqc_system_compute_charges: this is a closed-shell RHF and the "// &
                        "system has an odd number of electrons"
         status = MQC_FAIL
         return
      end if

      allocate (symbols(n))
      do i = 1, n
         symbols(i) = element_number_to_symbol(h%geom%element_numbers(i))
      end do

      call build_libcint_molecule(h%geom%element_numbers, symbols, h%geom%coordinates, &
                                  basis_name, mol, error)
      if (failed(error, status)) return

      call run_libcint_rhf(mol, nelec, SCF_MAX_ITER, SCF_ENERGY_TOL, SCF_DENSITY_TOL, &
                           .false., scf, error)
      if (failed(error, status)) return
      if (.not. scf%converged) then
         last_message = "mqc_system_compute_charges: the SCF did not converge, so there "// &
                        "is no density to partition"
         status = MQC_FAIL
         return
      end if

      if (which == "mulliken") then
         call mol%overlap(overlap)
         call mulliken_charges(mol, scf%density, overlap, q, error)
      else
         call chelpg_charges(mol, scf%density, q, error, &
                             total_charge=real(h%geom%charge, dp))
      end if
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
      !! atom 3" is not a well-defined number without it -- two schemes
      !! disagree by design, and a caller comparing systems or caching results
      !! has to be able to tell which question was answered.
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
