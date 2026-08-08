!! Manual check that a supplied fragment list runs and gives the same answer
!!
!!     mpirun -np 4 ./build/check_run validation/inputs/prism.json
!!
!! Reads a deck the ordinary way, then runs it twice: once letting the driver
!! generate its own fragment list, and once handing it a list built here. Both
!! must produce the same energy, because the list handed over is the same list
!! the driver would have made.
!!
!! That is the whole point of the check. Everything about supplying a list --
!! the fraglist object, the C interface over it, the screening -- is only worth
!! anything if a list that should reproduce the default actually does. Once
!! that holds, a *different* list is a change of science rather than a change
!! of plumbing.
!!
!! It runs before any Python exists deliberately: this proves the MPI half
!! while the FFI is not yet in the way of seeing what broke.
program check_run
   use pic_types, only: dp, int64, default_int
   use pic_mpi_lib, only: comm_world, pic_mpi_init, pic_mpi_finalize
   use mqc_resources, only: resources_t
   use mqc_config_types, only: mqc_config_t
   use mqc_json_config_reader, only: read_json_config_file
   use mqc_config_adapter, only: driver_config_t, config_to_driver, config_to_system_geometry
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_driver, only: run_calculation
   use mqc_fraglist, only: fraglist_t
   use mqc_error, only: error_t
   implicit none

   type(resources_t) :: resources
   type(mqc_config_t) :: config
   type(driver_config_t) :: driver
   type(system_geometry_t) :: sys_geom
   type(fraglist_t) :: terms
   type(error_t) :: error
   character(len=512) :: deck
   integer :: stat, rank
   integer(int64) :: iterm
   integer(default_int), allocatable :: supplied(:, :)

   call pic_mpi_init()
   resources%mpi_comms%world_comm = comm_world()
   resources%mpi_comms%node_comm = resources%mpi_comms%world_comm%split()
   rank = resources%mpi_comms%world_comm%rank()

   call get_command_argument(1, deck, status=stat)
   if (stat /= 0) then
      if (rank == 0) write (*, "(A)") "usage: check_run <deck.json>"
      call finish(1)
   end if

   ! Every rank reads the deck, as the ordinary program does. A session would
   ! broadcast instead; that is a separate piece and not what this checks.
   call read_json_config_file(trim(deck), config, error)
   if (error%has_error()) then
      if (rank == 0) write (*, "(A)") "FAILED to read deck: "//error%get_message()
      call finish(1)
   end if
   call config_to_driver(config, driver)
   if (config%nmol > 0) then
      ! Multi-molecule decks go through run_multi_molecule_calculations, which
      ! runs one calculation per molecule; a single supplied term list has no
      ! meaning across several. Not a limitation of supplying lists -- just not
      ! what this driver checks.
      if (rank == 0) write (*, "(A)") "check_run takes single-molecule decks; "// &
         trim(deck)//" defines several"
      call finish(1)
   end if

   call config_to_system_geometry(config, sys_geom, error)
   if (error%has_error()) then
      if (rank == 0) write (*, "(A)") "FAILED to build geometry: "//error%get_message()
      call finish(1)
   end if

   if (driver%nlevel < 2) then
      if (rank == 0) write (*, "(A)") "check_run needs a fragmented deck (nlevel >= 2)"
      call finish(1)
   end if

   ! ---- pass 1: let the driver generate its own list ------------------------
   if (rank == 0) write (*, "(A)") "[pass 1] driver-generated fragment list"
   call run_calculation(resources, driver, sys_geom, sys_geom%bonds)

   ! ---- pass 2: hand over the same list -------------------------------------
   ! Monomers first, then every n-mer, which is the order the driver builds --
   ! and then closed under subsets, which for a full list changes nothing but
   ! is what any screened list would need.
   if (rank == 0) then
      call terms%create(int(sys_geom%n_monomers, default_int), &
                        int(driver%nlevel, default_int), error)
      if (error%has_error()) then
         write (*, "(A)") "FAILED to build terms: "//error%get_message()
         call finish(1)
      end if
      call terms%close_subsets(error)
      if (error%has_error()) then
         write (*, "(A)") "FAILED to close subsets: "//error%get_message()
         call finish(1)
      end if

      allocate (supplied(terms%n_terms, terms%max_level))
      do iterm = 1, terms%n_terms
         supplied(iterm, :) = terms%terms(iterm, :)
      end do
      write (*, "(A,I0,A)") "[pass 2] supplied fragment list, ", terms%n_terms, " terms"
   else
      ! Only rank 0 needs the list: the coordinator hands out fragment
      ! definitions and workers build them from sys_geom.
      allocate (supplied(1, 1))
      supplied = 0
   end if

   call run_calculation(resources, driver, sys_geom, sys_geom%bonds, &
                        supplied_terms=supplied, n_supplied_terms=terms%n_terms)

   if (rank == 0) then
      write (*, "(A)") ""
      write (*, "(A)") "Compare the two total energies above: a supplied list that"
      write (*, "(A)") "reproduces the default must give the same number."
      write (*, "(A)") ""
      write (*, "(A)") "Pass a second argument to instead drop every monomer and keep the"
      write (*, "(A)") "higher terms -- the shape a naive energy screen produces. That"
      write (*, "(A)") "must be refused rather than run."
   end if

   ! ---- pass 3, on request: a screen that looks fine and is not -------------
   if (command_argument_count() > 1) then
      if (rank == 0) then
         ! Drop the monomers and keep the n-mers. Removing the top level
         ! instead would leave a list that is still closed -- the subsets
         ! of what remains are all present -- and rightly pass.
         call drop_all_of_level(terms, 1_default_int)
         deallocate (supplied)
         allocate (supplied(terms%n_terms, terms%max_level))
         do iterm = 1, terms%n_terms
            supplied(iterm, :) = terms%terms(iterm, :)
         end do
         write (*, "(A)") ""
         write (*, "(A,I0,A)") "[pass 3] monomers dropped, ", terms%n_terms, &
            " terms -- this should be refused"
      end if
      call run_calculation(resources, driver, sys_geom, sys_geom%bonds, &
                           supplied_terms=supplied, n_supplied_terms=terms%n_terms)
      if (rank == 0) write (*, "(A)") "REACHED HERE -- the bad list was not refused"
   end if

   call terms%destroy()
   call finish(0)

contains

   subroutine drop_all_of_level(list, level)
      !! Remove every term of one order, leaving the rest
      type(fraglist_t), intent(inout) :: list
      integer(default_int), intent(in) :: level

      logical, allocatable :: mask(:)
      type(error_t) :: err
      integer(int64) :: i

      allocate (mask(list%n_terms))
      do i = 1, list%n_terms
         mask(i) = (list%level_of(i) /= level)
      end do
      call list%keep(mask, err)
      deallocate (mask)
   end subroutine drop_all_of_level

   subroutine finish(code)
      integer, intent(in) :: code
      call resources%mpi_comms%node_comm%finalize()
      call resources%mpi_comms%world_comm%finalize()
      call pic_mpi_finalize()
      if (code /= 0) error stop 1
      stop
   end subroutine finish

end program check_run
