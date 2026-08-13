!! The C entry point that actually computes something
module mqc_capi_run
   !! Runs a calculation on a system handle, optionally with a supplied term
   !! list, across every rank of the session.
   !!
   !! This is where the three handles meet. A caller builds a system, builds or
   !! screens a term list, says what to compute in a JSON settings document,
   !! and hands all three over; the job's other ranks -- which have seen none
   !! of it -- are brought up to date and take part.
   !!
   !! **The settings travel as JSON text rather than as a struct.** That looks
   !! like a round trip and is a deliberate one. The alternative is
   !! broadcasting `method_config_t`, which is sixty-odd fields across seven
   !! nested sub-configs, and every one of them would have to be transcribed
   !! into a send and a receive that stay in step forever. A field missed there
   !! does not fail; it leaves the workers running a slightly different method
   !! than rank 0 asked for, and the run converges and reports a number. Text
   !! has no field list to forget, and it puts the workers through the same
   !! validator as any input deck rather than around it.
   !!
   !! The system does *not* travel that way, because it is the opposite case:
   !! it carries conventions -- the monomer padding, the 0-based atoms against
   !! 1-based monomers -- that a document would have to re-derive on arrival,
   !! which is where the silent wrongness would live for that one.
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_ptr, &
                                                                             c_null_char, c_f_pointer, c_associated
   use pic_types, only: int64, default_int
   use mqc_capi_system, only: system_handle_t
   use mqc_capi_session, only: current_session
   use mqc_fraglist, only: fraglist_t
   use mqc_session, only: mqc_session_t
   use mqc_result_types, only: calculation_result_t
   use mqc_error, only: error_t
   use mqc_config_types, only: bond_t
   use mqc_bond_perception, only: missing_broken_bonds
   implicit none
   private

   public :: mqc_run, mqc_run_last_error

   integer(c_int), parameter :: MQC_OK = 0
   integer(c_int), parameter :: MQC_FAIL = 1
   integer(c_int), parameter :: MQC_BAD_HANDLE = 2

   integer, parameter :: MESSAGE_LEN = 512
   character(len=MESSAGE_LEN), save :: last_message = ""

contains

   function mqc_run(system, terms, settings_len, settings, label_len, label, &
                    write_output, energy) result(status) bind(C, name="mqc_run")
      !! Compute, and give the caller back a number
      !!
      !! `terms` may be null, which means "generate the list yourself" -- the
      !! ordinary fragmented run, with distance screening and sorting. A list
      !! that *is* supplied is taken as given, in the order given, and must be
      !! closed under subsets; `run_calculation` refuses one that is not.
      !!
      !! `energy` is the total. The per-term breakdown is not returned, because
      !! it is what a screening pass wants and it is already written to the
      !! fragment breakdown file -- which is why `write_output` and a returned
      !! result are independent rather than a single switch.
      integer(c_int), value :: settings_len, label_len
      type(c_ptr), value :: system
      type(c_ptr), value :: terms  !! May be null: generate the list instead
      character(kind=c_char), intent(in) :: settings(settings_len)
      character(kind=c_char), intent(in) :: label(label_len)
      integer(c_int), value :: write_output
      real(c_double), intent(out) :: energy
      integer(c_int) :: status

      type(system_handle_t), pointer :: h
      type(fraglist_t), pointer :: list
      type(mqc_session_t), pointer :: session
      type(calculation_result_t) :: result
      type(error_t) :: error
      integer, allocatable :: supplied(:, :)
      integer(int64) :: n_supplied, iterm

      energy = 0.0_c_double
      status = MQC_BAD_HANDLE
      if (.not. c_associated(system)) then
         last_message = "mqc_run: null system handle"
         return
      end if
      call c_f_pointer(system, h)

      ! A caller who never mentioned bonds is not a caller with no bonds, and
      ! the difference is invisible further down: `validate_system` audits the
      ! declared list against the geometry, and an empty list has nothing to
      ! audit. So a covalent molecule cut into monomers with no bonds declared
      ! would fragment into radicals, cap nothing, and return a number.
      ! `bonds_declared` is the only place that distinction survives, so this
      ! is the only place the check can be made.
      if (.not. h%bonds_declared) then
         call refuse_undeclared_cuts(h, status)
         if (status /= MQC_OK) return
      end if

      session => current_session()
      if (.not. session%active) then
         last_message = "mqc_run: no session; call mqc_session_begin first"
         status = MQC_FAIL
         return
      end if

      ! An empty list still has to be an allocated array of some shape --
      ! `run_calculation` takes it by assumed shape, and n_supplied is what
      ! says how much of it means anything.
      if (c_associated(terms)) then
         call c_f_pointer(terms, list)
         n_supplied = list%n_terms
         allocate (supplied(max(list%n_terms, 1_int64), max(list%max_level, 1)))
         supplied = 0
         do iterm = 1, list%n_terms
            supplied(iterm, 1:list%max_level) = int(list%terms(iterm, 1:list%max_level))
         end do
      else
         n_supplied = 0
         allocate (supplied(1, 1))
         supplied = 0
      end if

      call session%run(chars_to_string(settings_len, settings), h%geom, &
                       supplied, n_supplied, chars_to_string(label_len, label), &
                       write_output /= 0_c_int, result, error)
      deallocate (supplied)

      if (error%has_error()) then
         last_message = error%get_message()
         status = MQC_FAIL
         return
      end if

      ! A failure inside a fragment lands on the result, not on the session
      ! error: the driver keeps going so it can report every bad term at the
      ! end rather than the first. Checking only `error` therefore returned
      ! MQC_OK with whatever total had accumulated -- for a gradient the CPU
      ! backend refuses, that was a clean 0.0 handed back as a success, which
      ! is the worst shape a failure can take.
      if (result%has_error) then
         last_message = result%error%get_message()
         status = MQC_FAIL
         return
      end if

      energy = real(result%energy%total(), c_double)
      status = MQC_OK
   end function mqc_run

   subroutine mqc_run_last_error(buffer_len, buffer) bind(C, name="mqc_run_last_error")
      !! Copy the most recent failure message out as a C string
      integer(c_int), value :: buffer_len
      character(kind=c_char), intent(inout) :: buffer(buffer_len)

      integer :: n, i

      if (buffer_len <= 0) return
      n = min(len_trim(last_message), int(buffer_len) - 1)
      do i = 1, n
         buffer(i) = last_message(i:i)
      end do
      buffer(n + 1) = c_null_char
   end subroutine mqc_run_last_error

   subroutine refuse_undeclared_cuts(h, status)
      !! Refuse a partition that cuts bonds nobody declared
      !!
      !! Perception against an empty declared list, so every cross-monomer bond
      !! the geometry implies comes back as missing. For a molecular cluster
      !! that is none and this costs one O(N^2) sweep; for a peptide chopped at
      !! the backbone it is the difference between capped fragments and
      !! radicals.
      type(system_handle_t), intent(in) :: h
      integer(c_int), intent(out) :: status

      type(bond_t), allocatable :: none(:)
      integer, allocatable :: missing_i(:), missing_j(:)
      integer :: n_missing

      status = MQC_OK
      if (h%geom%total_atoms <= 0 .or. h%geom%n_monomers <= 0) return

      allocate (none(0))
      call missing_broken_bonds(h%geom, none, 0, missing_i, missing_j, n_missing)
      deallocate (none)
      if (n_missing <= 0) return

      write (last_message, "(A,I0,A,I0,A,I0,A)") &
         "the monomers cut ", n_missing, " bond(s) that were never declared, "// &
         "starting with atoms ", missing_i(1), " and ", missing_j(1), &
         "; declare them with mqc_system_set_bonds so the fragments are capped, "// &
         "or set unchecked_input if the partition is deliberate"
      status = MQC_FAIL
   end subroutine refuse_undeclared_cuts

   pure function chars_to_string(n, chars) result(text)
      !! A C character array as a Fortran string
      !!
      !! Length-counted rather than null-terminated: the caller already knows
      !! how long the settings document is, and a JSON payload has no business
      !! being scanned for a terminator it might legitimately never reach.
      integer(c_int), intent(in) :: n
      character(kind=c_char), intent(in) :: chars(n)
      character(len=:), allocatable :: text

      integer :: i

      allocate (character(len=max(int(n), 0)) :: text)
      do i = 1, int(n)
         text(i:i) = chars(i)
      end do
   end function chars_to_string

end module mqc_capi_run
