!! Bond orders over the C boundary, for callers building their own connectivity
module mqc_capi_bond_orders
   !! Wiberg-Mayer bond orders from xTB, on a system handle.
   !!
   !! Deliberately not driven from a deck. A bond order is not something a
   !! calculation is *for* -- it is an input to deciding what calculation to run,
   !! and that deciding is exploratory: which bonds to cut, how large a cap has
   !! to be before the order stops moving, whether a fragmentation reproduces
   !! the parent's connectivity. That is a Python loop over many trial
   !! partitions, not a JSON key.
   !!
   !! So the shape here is compute-then-read, like the rest of the system API:
   !! `compute` runs one xTB single point over the whole system and stores the
   !! matrix on the handle, and `get` copies it out. A caller that wants twenty
   !! trial fragmentations computes once and reads twenty times.
   !!
   !! **What the numbers are good for.** Measured on cases with known answers,
   !! GFN2 gives 1.03 / 2.03 / 3.00 for the C-C bonds of ethane, ethene and
   !! ethyne, and 0.019 for a water dimer's hydrogen bond. A real single bond
   !! and a hydrogen bond are separated by a factor of fifty, which is the
   !! discrimination a covalent-radius rule cannot make at all -- it sees only
   !! the distance, and calls both of them bonds or neither.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_associated, c_f_pointer
   use pic_types, only: dp
   use mqc_capi_system, only: system_handle_t, last_message
   use mqc_capi_status, only: MQC_OK, MQC_FAIL, MQC_BAD_HANDLE
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
#ifndef MQC_WITHOUT_TBLITE
   use mqc_method_xtb, only: xtb_method_t
#endif
   implicit none
   private

   public :: mqc_system_compute_bond_orders
   public :: mqc_system_get_bond_orders
   public :: mqc_system_bond_order
   public :: mqc_system_has_bond_orders

contains

   function mqc_system_compute_bond_orders(handle, variant_len, variant, accuracy) &
      result(status) bind(C, name="mqc_system_compute_bond_orders")
      !! Run one xTB single point over the whole system and keep its bond orders
      !!
      !! The whole system, not the monomers: the point of these is to decide
      !! where the monomers should be, so a partition cannot be an input. A
      !! caller wanting fragment-local orders builds a handle per fragment.
      !!
      !! `variant` is "gfn2" or "gfn1"; an empty string takes gfn2, which is the
      !! one whose bond orders were checked. `accuracy` <= 0 takes tblite's
      !! default.
      type(c_ptr), value :: handle
      integer(c_int), value :: variant_len
      character(kind=c_char), intent(in) :: variant(variant_len)
      real(c_double), value :: accuracy
      integer(c_int) :: status

      type(system_handle_t), pointer :: h
#ifndef MQC_WITHOUT_TBLITE
      type(physical_fragment_t) :: whole
      type(calculation_result_t) :: res
      type(xtb_method_t) :: xtb
      character(len=:), allocatable :: name
      integer :: i
#endif

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null system handle"
         return
      end if
      call c_f_pointer(handle, h)

#ifdef MQC_WITHOUT_TBLITE
      ! Bond orders come from an xTB single point and there is no other route to
      ! them here. Declining with the same signature, so a C or Python caller
      ! still links and gets told why rather than finding the symbol missing.
      last_message = "mqc_system_compute_bond_orders: bond orders need tblite; "// &
                     "build with -DMQC_ENABLE_TBLITE=ON"
      status = MQC_FAIL
      return
#else
      if (h%geom%total_atoms <= 0) then
         last_message = "mqc_system_compute_bond_orders: set the geometry first"
         status = MQC_FAIL
         return
      end if

      name = ""
      do i = 1, variant_len
         name = name//variant(i)
      end do
      if (len_trim(name) == 0) name = "gfn2"

      ! The system as one fragment. Charge and multiplicity come from the
      ! geometry rather than from any monomer, because this is the parent.
      whole%n_atoms = h%geom%total_atoms
      whole%charge = h%geom%charge
      whole%multiplicity = h%geom%multiplicity
      allocate (whole%element_numbers(whole%n_atoms), source=h%geom%element_numbers)
      allocate (whole%coordinates(3, whole%n_atoms), source=h%geom%coordinates)
      whole%nelec = sum(h%geom%element_numbers) - h%geom%charge

      xtb%variant = trim(adjustl(name))
      xtb%want_bond_orders = .true.
      if (accuracy > 0.0_dp) xtb%accuracy = accuracy
      call xtb%calc_energy(whole, res)

      if (res%has_error) then
         last_message = "mqc_system_compute_bond_orders: "//res%error%get_message()
         status = MQC_FAIL
         return
      end if
      if (.not. res%has_bond_orders) then
         last_message = "mqc_system_compute_bond_orders: the calculation returned no "// &
                        "bond orders"
         status = MQC_FAIL
         return
      end if

      if (allocated(h%bond_orders)) deallocate (h%bond_orders)
      allocate (h%bond_orders(whole%n_atoms, whole%n_atoms), source=res%bond_orders)
      status = MQC_OK
#endif
   end function mqc_system_compute_bond_orders

   function mqc_system_has_bond_orders(handle) result(has) &
      bind(C, name="mqc_system_has_bond_orders")
      !! Whether `compute` has run on this handle. Zero for a null handle too.
      type(c_ptr), value :: handle
      integer(c_int) :: has

      type(system_handle_t), pointer :: h

      has = 0
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, h)
      if (allocated(h%bond_orders)) has = 1
   end function mqc_system_has_bond_orders

   function mqc_system_get_bond_orders(handle, n, out) result(status) &
      bind(C, name="mqc_system_get_bond_orders")
      !! Copy the matrix into a caller buffer of n*n doubles
      !!
      !! Column-major, which is what the caller's `n_atoms` and Fortran agree on
      !! and what numpy reads with `order="F"`. Symmetric, so a caller that gets
      !! the order wrong sees the right answer anyway -- which is a reason to
      !! say it here rather than rely on it.
      type(c_ptr), value :: handle
      integer(c_int), value :: n
      real(c_double), intent(out) :: out(n*n)
      integer(c_int) :: status

      type(system_handle_t), pointer :: h
      integer :: i, j, k

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null system handle"
         return
      end if
      call c_f_pointer(handle, h)

      if (.not. allocated(h%bond_orders)) then
         last_message = "mqc_system_get_bond_orders: compute them first"
         status = MQC_FAIL
         return
      end if
      if (n /= size(h%bond_orders, 1)) then
         last_message = "mqc_system_get_bond_orders: buffer is the wrong size for "// &
                        "this system"
         status = MQC_FAIL
         return
      end if

      k = 0
      do j = 1, n
         do i = 1, n
            k = k + 1
            out(k) = h%bond_orders(i, j)
         end do
      end do
      status = MQC_OK
   end function mqc_system_get_bond_orders

   function mqc_system_bond_order(handle, i, j) result(order) &
      bind(C, name="mqc_system_bond_order")
      !! One element, 0-based, for a caller probing a single pair
      !!
      !! Negative for anything wrong -- no handle, no orders, indices out of
      !! range -- because a bond order is never negative and a caller asking one
      !! pair at a time does not want a status code per question.
      type(c_ptr), value :: handle
      integer(c_int), value :: i, j
      real(c_double) :: order

      type(system_handle_t), pointer :: h

      order = -1.0_dp
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, h)
      if (.not. allocated(h%bond_orders)) return
      if (i < 0 .or. j < 0) return
      if (i >= size(h%bond_orders, 1) .or. j >= size(h%bond_orders, 2)) return
      order = h%bond_orders(i + 1, j + 1)
   end function mqc_system_bond_order

end module mqc_capi_bond_orders
