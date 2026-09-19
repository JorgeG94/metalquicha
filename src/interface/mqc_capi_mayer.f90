!! Mayer bond orders over the C boundary
module mqc_capi_mayer
   !! The ab initio half of `compute_bond_orders`: one RHF in the basis the
   !! caller names, and Mayer's orders off its converged density.
   !!
   !! Split from [[mqc_capi_bond_orders]] rather than added to it because it
   !! needs the integrals backend and that module does not. A build without
   !! one keeps the xTB entry and simply lacks these symbols, which a Python
   !! caller sees as a named refusal rather than a link error -- the same
   !! arrangement [[mqc_capi_charges]] has, for the same reason.
   !!
   !! **What this costs.** An SCF, not an xTB single point. The two variants
   !! behind one Python call differ by orders of magnitude in price and by a
   !! whole Hamiltonian in meaning; the cheap one is for deciding where to cut
   !! a molecule, and this one is for finding out whether that decision was
   !! sound.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_associated, &
                                          c_f_pointer
   use pic_types, only: dp
   use mqc_capi_system, only: system_handle_t, last_message
   use mqc_capi_status, only: MQC_OK, MQC_FAIL, MQC_BAD_HANDLE
   use mqc_error, only: error_t
   use mqc_elements, only: element_number_to_symbol
   use mqc_czt_bridge, only: run_czt_mayer_bond_orders
   implicit none
   private

   public :: mqc_system_compute_mayer_bond_orders

contains

   function mqc_system_compute_mayer_bond_orders(handle, basis_len, basis) &
      result(status) bind(C, name="mqc_system_compute_mayer_bond_orders")
      !! Run one RHF over the whole system and keep Mayer's bond orders
      !!
      !! `basis` is any basis the build carries; an empty string takes 6-31g,
      !! which is a real default and not a dead initialiser -- a C caller has
      !! no JSON layer behind it to supply one.
      !!
      !! Closed shell only: an odd electron count is refused rather than
      !! quietly paired up. The open-shell formula exists and is reachable
      !! from a deck; it is the RHF behind this entry that has no place to put
      !! a multiplicity.
      type(c_ptr), value :: handle
      integer(c_int), value :: basis_len
      character(kind=c_char), intent(in) :: basis(basis_len)
      integer(c_int) :: status

      type(system_handle_t), pointer :: h
      type(error_t) :: error
      real(dp), allocatable :: orders(:, :), valences(:)
      character(len=2), allocatable :: symbols(:)
      character(len=:), allocatable :: basis_name
      integer :: i, n

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null system handle"
         return
      end if
      call c_f_pointer(handle, h)

      if (h%geom%total_atoms <= 0) then
         last_message = "mqc_system_compute_mayer_bond_orders: set the geometry first"
         status = MQC_FAIL
         return
      end if

      basis_name = ""
      do i = 1, basis_len
         basis_name = basis_name//basis(i)
      end do
      basis_name = trim(adjustl(basis_name))
      if (len_trim(basis_name) == 0) basis_name = "6-31g"

      n = h%geom%total_atoms
      allocate (symbols(n))
      do i = 1, n
         symbols(i) = element_number_to_symbol(h%geom%element_numbers(i))
      end do

      ! Through the bridge, not into the backend: the odd-electron refusal and
      ! the SCF both live on the other side of it.
      call run_czt_mayer_bond_orders(h%geom%element_numbers, symbols, &
                                     h%geom%coordinates, basis_name, h%geom%charge, &
                                     orders, valences, error)
      if (error%has_error()) then
         last_message = "mqc_system_compute_mayer_bond_orders: "//error%get_message()
         status = MQC_FAIL
         return
      end if

      if (allocated(h%bond_orders)) deallocate (h%bond_orders)
      if (allocated(h%bond_order_valences)) deallocate (h%bond_order_valences)
      call move_alloc(orders, h%bond_orders)
      call move_alloc(valences, h%bond_order_valences)
      h%bond_order_scheme = "mayer"
      status = MQC_OK
   end function mqc_system_compute_mayer_bond_orders

end module mqc_capi_mayer
