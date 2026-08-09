!! Manual check that a session-driven run matches an ordinary one
!!
!!     mpirun -np 4 ./build/check_session_run validation/inputs/prism.json
!!     ./build/mqc validation/inputs/prism.json          # the number to match
!!
!! **Only rank 0 reads the deck.** That is the whole point, and it is what
!! separates this from `check_run`, where every rank parsed the input the way
!! `main` does. Here the workers are parked inside `mqc_session_begin` before
!! the filename is even looked at; everything they know about the calculation
!! arrived over a broadcast. If the geometry transport or the settings
!! transport is wrong, this is where it shows -- as a wrong energy or a hang,
!! rather than as anything that looks like a bug in the science.
!!
!! The settings are given as a JSON literal rather than re-serialized from the
!! deck, deliberately: that exercises the string parser a driving script will
!! actually use, and it demonstrates the shape such a script must produce --
!! the deck, exactly, less the molecules.
program check_session_run
   use pic_types, only: dp, int64, default_int
   use mqc_session, only: mqc_session_t
   use mqc_config_types, only: mqc_config_t
   use mqc_json_config_reader, only: read_json_config_file
   use mqc_config_adapter, only: config_to_system_geometry
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_result_types, only: calculation_result_t
   use mqc_fraglist, only: fraglist_t
   use mqc_error, only: error_t
   implicit none

   !> prism.json with the molecules removed. Kept on one logical line per
   !  section so a mismatch against the deck is easy to see by eye.
   character(len=*), parameter :: SETTINGS = &
      '{"schema":{"name":"mqc-frag","version":"1.0"},'// &
      '"model":{"method":"XTB-GFN1","basis":"cc-pVDZ","aux_basis":"cc-pVDZ-RIFIT"},'// &
      '"keywords":{"scf":{"maxiter":300,"tolerance":1e-06},'// &
      '"fragmentation":{"method":"MBE","allow_overlapping_fragments":false,'// &
      '"level":2,"embedding":"none","cutoff_method":"distance",'// &
      '"distance_metric":"min"}},'// &
      '"system":{"logger":{"level":"info"}},'// &
      '"driver":"Energy"}'

   type(mqc_session_t) :: session
   type(mqc_config_t) :: config
   type(system_geometry_t) :: sys_geom
   type(calculation_result_t) :: result
   type(fraglist_t) :: terms
   type(error_t) :: error
   character(len=512) :: deck
   integer :: stat
   integer(int64) :: iterm
   integer, allocatable :: supplied(:, :)

   ! Nothing above this line runs anywhere but here on rank 0 -- every other
   ! rank enters `begin` and does not come back.
   call session%begin(error)
   if (error%has_error()) then
      write (*, "(A)") "FAILED to start session: "//error%get_message()
      error stop 1
   end if

   call get_command_argument(1, deck, status=stat)
   if (stat /= 0) then
      write (*, "(A)") "usage: check_session_run <deck.json>"
      call finish(1)
   end if

   call read_json_config_file(trim(deck), config, error)
   if (error%has_error()) then
      write (*, "(A)") "FAILED to read deck: "//error%get_message()
      call finish(1)
   end if
   call config_to_system_geometry(config, sys_geom, error)
   if (error%has_error()) then
      write (*, "(A)") "FAILED to build geometry: "//error%get_message()
      call finish(1)
   end if

   write (*, "(A,I0,A,I0,A)") "driving ", session%n_ranks, " rank(s) over ", &
      sys_geom%n_monomers, " monomers; only this rank has read the deck"

   ! ---- pass 1: let the driver generate its own list ------------------------
   allocate (supplied(1, 1))
   supplied = 0
   call session%run(SETTINGS, sys_geom, supplied, 0_int64, trim(deck), .false., &
                    result, error)
   if (error%has_error()) then
      write (*, "(A)") "FAILED pass 1: "//error%get_message()
      call finish(1)
   end if
   write (*, "(A,F20.10)") "[pass 1] generated list ", result%energy%total()
   deallocate (supplied)

   ! ---- pass 2: hand over the same list -------------------------------------
   ! A full list closed under subsets is the list the driver would have built,
   ! so the two energies must agree. Once they do, a *screened* list is a
   ! change of science rather than a change of plumbing.
   call terms%create(int(sys_geom%n_monomers, default_int), 2_default_int, error)
   if (error%has_error()) then
      write (*, "(A)") "FAILED to build terms: "//error%get_message()
      call finish(1)
   end if
   ! `create` gives the n-mers of the top level only -- 15 dimers here, not the
   ! 21 terms an MBE(2) needs. Leaving this out is exactly the mistake a naive
   ! screen makes, and it was caught by the closure check rather than by
   ! producing a number, which is the point of having one.
   call terms%close_subsets(error)
   if (error%has_error()) then
      write (*, "(A)") "FAILED to close subsets: "//error%get_message()
      call finish(1)
   end if
   allocate (supplied(terms%n_terms, terms%max_level))
   do iterm = 1, terms%n_terms
      supplied(iterm, :) = int(terms%terms(iterm, :))
   end do

   call session%run(SETTINGS, sys_geom, supplied, terms%n_terms, trim(deck), .true., &
                    result, error)
   if (error%has_error()) then
      write (*, "(A)") "FAILED pass 2: "//error%get_message()
      call finish(1)
   end if
   write (*, "(A,I0,A,F20.10)") "[pass 2] supplied ", terms%n_terms, " terms ", &
      result%energy%total()

   write (*, "(A)") ""
   write (*, "(A)") "Both must equal ./build/mqc on the same deck. A session that"
   write (*, "(A)") "moved the geometry or the settings wrongly lands here, not in"
   write (*, "(A)") "anything that looks like a bug in the science."

   ! ---- pass 3: a bad list must come back, not take the process with it -----
   ! Dropping the monomers leaves the shape a naive energy screen produces.
   ! Reaching the driver, that aborts the job -- correct for `main`, fatal for
   ! a caller who has just spent an hour building state. The session has to
   ! refuse it before the workers are committed and still be usable after.
   call drop_monomers(terms)
   deallocate (supplied)
   allocate (supplied(terms%n_terms, terms%max_level))
   do iterm = 1, terms%n_terms
      supplied(iterm, :) = int(terms%terms(iterm, :))
   end do
   write (*, "(A)") ""
   call session%run(SETTINGS, sys_geom, supplied, terms%n_terms, "", .false., result, error)
   if (error%has_error()) then
      write (*, "(A)") "[pass 3] refused, as it must be: "//error%get_message()
      call error%clear()
   else
      write (*, "(A)") "[pass 3] NOT REFUSED -- a list missing its monomers was run"
      call finish(1)
   end if

   ! Still alive, and still able to compute: a refusal must cost the caller
   ! the calculation, not the session.
   deallocate (supplied)
   allocate (supplied(1, 1))
   supplied = 0
   call session%run(SETTINGS, sys_geom, supplied, 0_int64, "", .false., result, error)
   if (error%has_error()) then
      write (*, "(A)") "[pass 4] session unusable after a refusal: "//error%get_message()
      call finish(1)
   end if
   write (*, "(A,F20.10)") "[pass 4] session still works after refusing ", &
      result%energy%total()

   call terms%destroy()
   call finish(0)

contains

   subroutine drop_monomers(list)
      !! Keep every n-mer and throw the monomers away
      type(fraglist_t), intent(inout) :: list

      logical, allocatable :: mask(:)
      type(error_t) :: err
      integer(int64) :: i

      allocate (mask(list%n_terms))
      do i = 1, list%n_terms
         mask(i) = (list%level_of(i) /= 1)
      end do
      call list%keep(mask, err)
      deallocate (mask)
   end subroutine drop_monomers

   subroutine finish(code)
      integer, intent(in) :: code
      type(error_t) :: err
      call session%end_session(err)
      if (code /= 0) error stop 1
      stop
   end subroutine finish

end program check_session_run
