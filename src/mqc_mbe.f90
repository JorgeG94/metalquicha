!! Many-Body Expansion (MBE) calculation module
module mqc_mbe
   !! Implements hierarchical many-body expansion for fragment-based quantum chemistry
   !! calculations with MPI parallelization and energy/gradient computation.
   use pic_types, only: int32, int64, dp
   use pic_timer, only: timer_type
   use pic_mpi_lib, only: comm_t, send, recv, iprobe, MPI_Status, MPI_ANY_SOURCE, MPI_ANY_TAG, abort_comm
   use pic_logger, only: logger => global_logger, verbose_level, debug_level, info_level
   use pic_io, only: to_char
   use mqc_mbe_io, only: print_detailed_breakdown
   use mqc_json_output_types, only: json_output_data_t, OUTPUT_MODE_MBE
   use mqc_thermochemistry, only: thermochemistry_result_t, compute_thermochemistry
   use mqc_mpi_tags, only: TAG_WORKER_REQUEST, TAG_WORKER_FRAGMENT, TAG_WORKER_FINISH, &
                           TAG_WORKER_SCALAR_RESULT, &
                           TAG_NODE_REQUEST, TAG_NODE_FRAGMENT, TAG_NODE_FINISH, &
                           TAG_NODE_SCALAR_RESULT
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, build_fragment_from_indices, to_angstrom
   use mqc_frag_utils, only: get_next_combination, fragment_lookup_t
   use mqc_vibrational_analysis, only: compute_vibrational_analysis, print_vibrational_analysis
   use mqc_program_limits, only: MAX_MBE_LEVEL

   implicit none
   private

   ! Public interface
   public :: compute_mbe   !! MBE energy with optional gradient and hessian
   public :: compute_gmbe  !! GMBE energy with optional gradient and hessian

contains

   function compute_mbe_delta(fragment_idx, fragment, lookup, energies, delta_energies, n, world_comm) result(delta_E)
      !! Bottom-up computation of n-body correction (non-recursive, uses pre-computed subset deltas)
      !! deltaE(i1,i2,...,in) = E(i1,i2,...,in) - sum of all subset deltaE values
      !! All subsets must have been computed already (guaranteed by processing fragments in order)
      integer(int64), intent(in) :: fragment_idx  !! Index of this fragment (already known)
      integer, intent(in) :: fragment(:), n
      type(fragment_lookup_t), intent(in) :: lookup  !! Pre-built hash table for lookups
      real(dp), intent(in) :: energies(:), delta_energies(:)  !! Pre-computed delta values
      type(comm_t), intent(in), optional :: world_comm  !! MPI communicator for abort
      real(dp) :: delta_E

      integer :: subset_size, i
      integer :: indices(MAX_MBE_LEVEL), subset(MAX_MBE_LEVEL)  ! Stack arrays to avoid heap contention
      integer(int64) :: subset_idx
      logical :: has_next

      ! Start with the full n-mer energy
      delta_E = energies(fragment_idx)

      ! Subtract all proper subsets (size 1 to n-1)
      do subset_size = 1, n - 1
         ! Initialize first combination
         do i = 1, subset_size
            indices(i) = i
         end do

         ! Loop through all combinations
         do
            ! Build current subset
            do i = 1, subset_size
               subset(i) = fragment(indices(i))
            end do

            ! Look up subset index
            subset_idx = lookup%find(subset(1:subset_size), subset_size)
            if (subset_idx < 0) then
               block
                  use pic_io, only: to_char
                  character(len=512) :: error_msg
                  integer :: j
                  write (error_msg, '(a,i0,a,*(i0,1x))') "Subset not found! Fragment idx=", fragment_idx, &
                     " seeking subset: ", (subset(j), j=1, subset_size)
                  call logger%error(trim(error_msg))
                  write (error_msg, '(a,*(i0,1x))') "  Full fragment: ", (fragment(j), j=1, n)
                  call logger%error(trim(error_msg))
                  if (present(world_comm)) then
                     call abort_comm(world_comm, 1)
                  else
                     error stop "Subset not found in bottom-up MBE!"
                  end if
               end block
            end if

            ! Subtract pre-computed delta energy
            delta_E = delta_E - delta_energies(subset_idx)

            ! Get next combination
            call get_next_combination(indices, subset_size, n, has_next)
            if (.not. has_next) exit
         end do
      end do

   end function compute_mbe_delta

   subroutine map_fragment_to_system_gradient(frag_grad, monomers, sys_geom, sys_grad, bonds, world_comm)
      !! Map fragment gradient to system gradient coordinates with hydrogen cap redistribution
      !!
      !! This function rebuilds the fragment to get local→global mappings and cap information,
      !! then redistributes gradients including hydrogen caps to their original atoms.
      !!
      !! If bonds are not present, uses the old simple mapping (no caps possible).
      use mqc_physical_fragment, only: build_fragment_from_indices, redistribute_cap_gradients
      use mqc_config_parser, only: bond_t
      use mqc_error, only: error_t
      use pic_logger, only: verbose_level
      real(dp), intent(in) :: frag_grad(:, :)  !! (3, natoms_frag)
      integer, intent(in) :: monomers(:)  !! Monomer indices in fragment
      type(system_geometry_t), intent(in) :: sys_geom
      real(dp), intent(inout) :: sys_grad(:, :)  !! (3, total_atoms)
      type(bond_t), intent(in), optional :: bonds(:)  !! Bond information for caps
      type(comm_t), intent(in), optional :: world_comm  !! MPI communicator for abort

      type(physical_fragment_t) :: fragment
      type(error_t) :: error
      integer :: i_mon, i_atom, frag_atom_idx, sys_atom_idx
      integer :: current_log_level

      ! Explicitly zero out the entire sys_grad array
      sys_grad = 0.0_dp

      ! Debug output
      call logger%configuration(level=current_log_level)
      if (current_log_level >= debug_level) then
         block
            character(len=256) :: debug_msg
            write (debug_msg, '(a,i0,a,*(i0,1x))') "  Mapping fragment with ", size(monomers), " monomers: ", monomers
            call logger%debug(trim(debug_msg))
            write (debug_msg, '(a,i0,a)') "  Fragment has ", size(frag_grad, 2), " atoms"
            call logger%debug(trim(debug_msg))
         end block
      end if

      if (present(bonds)) then
         ! Rebuild fragment to get local→global mapping and cap information
         call build_fragment_from_indices(sys_geom, monomers, fragment, error, bonds)
         if (error%has_error()) then
            call logger%error(error%get_full_trace())
            if (present(world_comm)) then
               call abort_comm(world_comm, 1)
            else
               error stop "Failed to build fragment in gradient mapping"
            end if
         end if

         ! Use new gradient redistribution with cap handling
         call redistribute_cap_gradients(fragment, frag_grad, sys_grad)

         ! Clean up
         call fragment%destroy()
      else
         ! Old code path for fragments without hydrogen caps
         ! Map fragment gradient to system positions (fixed-size monomers only)
         frag_atom_idx = 0
         do i_mon = 1, size(monomers)
            do i_atom = 1, sys_geom%atoms_per_monomer
               frag_atom_idx = frag_atom_idx + 1
               sys_atom_idx = (monomers(i_mon) - 1)*sys_geom%atoms_per_monomer + i_atom

               if (current_log_level >= debug_level .and. i_atom == 1) then
                  block
                     character(len=256) :: debug_msg
                     write (debug_msg, '(a,i0,a,i0,a,i0)') &
                        "    Monomer ", monomers(i_mon), ": frag atoms ", &
                        frag_atom_idx, " -> sys atom ", sys_atom_idx
                     call logger%debug(trim(debug_msg))
                  end block
               end if

               sys_grad(:, sys_atom_idx) = frag_grad(:, frag_atom_idx)
            end do
         end do
      end if

   end subroutine map_fragment_to_system_gradient

 subroutine compute_mbe_gradient(fragment_idx, fragment, lookup, results, delta_gradients, n, sys_geom, bonds, world_comm)
      !! Bottom-up computation of n-body gradient correction
      !! Exactly mirrors the energy MBE logic: deltaG = G - sum(all subset deltaGs)
      !! All gradients are in system coordinates, so subtraction is simple
      use mqc_result_types, only: calculation_result_t
      use mqc_config_parser, only: bond_t
      integer(int64), intent(in) :: fragment_idx
      integer, intent(in) :: fragment(:), n
      type(fragment_lookup_t), intent(in) :: lookup
      type(calculation_result_t), intent(in) :: results(:)
      real(dp), intent(inout) :: delta_gradients(:, :, :)  !! (3, total_atoms, fragment_count)
      type(system_geometry_t), intent(in) :: sys_geom
      type(bond_t), intent(in), optional :: bonds(:)  !! Bond information for caps
      type(comm_t), intent(in), optional :: world_comm  !! MPI communicator for abort

      integer :: subset_size, i
      integer :: indices(MAX_MBE_LEVEL), subset(MAX_MBE_LEVEL)  ! Stack arrays to avoid heap contention
      integer(int64) :: subset_idx
      logical :: has_next

      ! Start with the full n-mer gradient mapped to system coordinates
      call map_fragment_to_system_gradient(results(fragment_idx)%gradient, fragment, &
                                           sys_geom, delta_gradients(:, :, fragment_idx), bonds, world_comm)

      ! Subtract all proper subsets (size 1 to n-1)
      ! This is EXACTLY like the energy calculation, but for each gradient component
      do subset_size = 1, n - 1
         ! Initialize first combination
         do i = 1, subset_size
            indices(i) = i
         end do

         ! Loop through all combinations
         do
            ! Build current subset
            do i = 1, subset_size
               subset(i) = fragment(indices(i))
            end do

            ! Look up subset index
            subset_idx = lookup%find(subset(1:subset_size), subset_size)
            if (subset_idx < 0) then
               call logger%error("Subset not found in MBE gradient computation")
               if (present(world_comm)) then
                  call abort_comm(world_comm, 1)
               else
                  error stop "Subset not found in MBE gradient!"
               end if
            end if

            ! Subtract pre-computed delta gradient (simple array subtraction in system coords)
            delta_gradients(:, :, fragment_idx) = delta_gradients(:, :, fragment_idx) - &
                                                  delta_gradients(:, :, subset_idx)

            ! Get next combination
            call get_next_combination(indices, subset_size, n, has_next)
            if (.not. has_next) exit
         end do
      end do

   end subroutine compute_mbe_gradient

   subroutine compute_mbe_dipole(fragment_idx, fragment, lookup, results, delta_dipoles, n, world_comm)
      !! Bottom-up computation of n-body dipole correction
      !! Exactly mirrors the energy MBE logic: deltaDipole = Dipole - sum(all subset deltaDipoles)
      !! Dipoles are additive vectors in the system frame, no coordinate mapping needed
      use mqc_result_types, only: calculation_result_t
      integer(int64), intent(in) :: fragment_idx
      integer, intent(in) :: fragment(:), n
      type(fragment_lookup_t), intent(in) :: lookup
      type(calculation_result_t), intent(in) :: results(:)
      real(dp), intent(inout) :: delta_dipoles(:, :)  !! (3, fragment_count)
      type(comm_t), intent(in), optional :: world_comm  !! MPI communicator for abort

      integer :: subset_size, i
      integer :: indices(MAX_MBE_LEVEL), subset(MAX_MBE_LEVEL)  ! Stack arrays to avoid heap contention
      integer(int64) :: subset_idx
      logical :: has_next

      ! Start with the full n-mer dipole
      delta_dipoles(:, fragment_idx) = results(fragment_idx)%dipole

      ! Subtract all proper subsets (size 1 to n-1)
      do subset_size = 1, n - 1
         ! Initialize first combination
         do i = 1, subset_size
            indices(i) = i
         end do

         ! Loop through all combinations
         do
            ! Build current subset
            do i = 1, subset_size
               subset(i) = fragment(indices(i))
            end do

            ! Look up subset index
            subset_idx = lookup%find(subset(1:subset_size), subset_size)
            if (subset_idx < 0) then
               call logger%error("Subset not found in MBE dipole computation")
               if (present(world_comm)) then
                  call abort_comm(world_comm, 1)
               else
                  error stop "Subset not found in MBE dipole!"
               end if
            end if

            ! Subtract pre-computed delta dipole
            delta_dipoles(:, fragment_idx) = delta_dipoles(:, fragment_idx) - &
                                             delta_dipoles(:, subset_idx)

            ! Get next combination
            call get_next_combination(indices, subset_size, n, has_next)
            if (.not. has_next) exit
         end do
      end do

   end subroutine compute_mbe_dipole

   subroutine map_fragment_to_system_hessian(frag_hess, monomers, sys_geom, sys_hess, bonds)
      !! Map fragment Hessian to system Hessian coordinates with hydrogen cap redistribution
      use mqc_physical_fragment, only: build_fragment_from_indices, redistribute_cap_hessian
      use mqc_config_parser, only: bond_t
      use mqc_error, only: error_t
      real(dp), intent(in) :: frag_hess(:, :)  !! (3*natoms_frag, 3*natoms_frag)
      integer, intent(in) :: monomers(:)
      type(system_geometry_t), intent(in) :: sys_geom
      real(dp), intent(inout) :: sys_hess(:, :)  !! (3*total_atoms, 3*total_atoms)
      type(bond_t), intent(in), optional :: bonds(:)

      type(physical_fragment_t) :: fragment
      type(error_t) :: error

      ! Zero out
      sys_hess = 0.0_dp

      if (present(bonds)) then
         ! Rebuild fragment to get local→global mapping and cap information
         call build_fragment_from_indices(sys_geom, monomers, fragment, error, bonds)
         call redistribute_cap_hessian(fragment, frag_hess, sys_hess)
         call fragment%destroy()
      else
         ! Old code path for fragments without hydrogen caps
         ! Map fragment Hessian to system positions (fixed-size monomers only)
         block
            integer :: i_mon, j_mon, i_atom, j_atom
            integer :: frag_atom_i, frag_atom_j, sys_atom_i, sys_atom_j
            integer :: frag_row_start, frag_col_start, sys_row_start, sys_col_start
            integer :: n_monomers

            n_monomers = size(monomers)
            frag_atom_i = 0

            ! Map each monomer's atoms
            do i_mon = 1, n_monomers
               do i_atom = 1, sys_geom%atoms_per_monomer
                  frag_atom_i = frag_atom_i + 1
                  sys_atom_i = (monomers(i_mon) - 1)*sys_geom%atoms_per_monomer + i_atom
                  frag_row_start = (frag_atom_i - 1)*3 + 1
                  sys_row_start = (sys_atom_i - 1)*3 + 1

                  ! Map this atom's Hessian blocks with all other atoms in fragment
                  frag_atom_j = 0
                  do j_mon = 1, n_monomers
                     do j_atom = 1, sys_geom%atoms_per_monomer
                        frag_atom_j = frag_atom_j + 1
                        sys_atom_j = (monomers(j_mon) - 1)*sys_geom%atoms_per_monomer + j_atom
                        frag_col_start = (frag_atom_j - 1)*3 + 1
                        sys_col_start = (sys_atom_j - 1)*3 + 1

                        ! Copy the 3×3 block for this atom pair
                        sys_hess(sys_row_start:sys_row_start + 2, sys_col_start:sys_col_start + 2) = &
                           frag_hess(frag_row_start:frag_row_start + 2, frag_col_start:frag_col_start + 2)
                     end do
                  end do
               end do
            end do
         end block
      end if

   end subroutine map_fragment_to_system_hessian

   subroutine map_fragment_to_system_dipole_derivatives(frag_dipole_derivs, monomers, sys_geom, sys_dipole_derivs, bonds)
      !! Map fragment dipole derivatives to system coordinates with hydrogen cap redistribution
      !!
      !! Dipole derivatives have shape (3, 3*N) where each column corresponds to
      !! the derivative of dipole w.r.t. one Cartesian coordinate of one atom.
      use mqc_physical_fragment, only: build_fragment_from_indices, redistribute_cap_dipole_derivatives
      use mqc_config_parser, only: bond_t
      use mqc_error, only: error_t
      real(dp), intent(in) :: frag_dipole_derivs(:, :)  !! (3, 3*natoms_frag)
      integer, intent(in) :: monomers(:)
      type(system_geometry_t), intent(in) :: sys_geom
      real(dp), intent(inout) :: sys_dipole_derivs(:, :)  !! (3, 3*total_atoms)
      type(bond_t), intent(in), optional :: bonds(:)

      type(physical_fragment_t) :: fragment
      type(error_t) :: error

      ! Zero out
      sys_dipole_derivs = 0.0_dp

      if (present(bonds)) then
         ! Rebuild fragment to get local→global mapping and cap information
         call build_fragment_from_indices(sys_geom, monomers, fragment, error, bonds)
         call redistribute_cap_dipole_derivatives(fragment, frag_dipole_derivs, sys_dipole_derivs)
         call fragment%destroy()
      else
         ! Code path for fragments without hydrogen caps
         ! Map fragment dipole derivative columns to system positions (fixed-size monomers only)
         block
            integer :: i_mon, i_atom
            integer :: frag_atom_idx, sys_atom_idx
            integer :: frag_col_start, sys_col_start
            integer :: n_monomers

            n_monomers = size(monomers)
            frag_atom_idx = 0

            ! Map each monomer's atoms
            do i_mon = 1, n_monomers
               do i_atom = 1, sys_geom%atoms_per_monomer
                  frag_atom_idx = frag_atom_idx + 1
                  sys_atom_idx = (monomers(i_mon) - 1)*sys_geom%atoms_per_monomer + i_atom
                  frag_col_start = (frag_atom_idx - 1)*3 + 1
                  sys_col_start = (sys_atom_idx - 1)*3 + 1

                  ! Copy the 3 columns (x, y, z derivatives) for this atom
                  sys_dipole_derivs(:, sys_col_start:sys_col_start + 2) = &
                     frag_dipole_derivs(:, frag_col_start:frag_col_start + 2)
               end do
            end do
         end block
      end if

   end subroutine map_fragment_to_system_dipole_derivatives

   subroutine compute_mbe_hessian(fragment_idx, fragment, lookup, results, delta_hessians, n, sys_geom, bonds)
      !! Bottom-up computation of n-body Hessian correction
      !! Mirrors MBE gradient logic but for second derivatives
      use mqc_result_types, only: calculation_result_t
      use mqc_config_parser, only: bond_t
      use mqc_error, only: error_t
      integer(int64), intent(in) :: fragment_idx
      integer, intent(in) :: fragment(:), n
      type(fragment_lookup_t), intent(in) :: lookup
      type(calculation_result_t), intent(in) :: results(:)
      real(dp), intent(inout) :: delta_hessians(:, :, :)  !! (3*total_atoms, 3*total_atoms, fragment_count)
      type(system_geometry_t), intent(in) :: sys_geom
      type(bond_t), intent(in), optional :: bonds(:)

      integer :: subset_size, i, hess_dim
      integer :: indices(MAX_MBE_LEVEL), subset(MAX_MBE_LEVEL)  ! Stack arrays to avoid heap contention
      integer(int64) :: subset_idx
      logical :: has_next

      hess_dim = 3*sys_geom%total_atoms

      ! Start with the full n-mer Hessian mapped to system coordinates
      call map_fragment_to_system_hessian(results(fragment_idx)%hessian, fragment, &
                                          sys_geom, delta_hessians(:, :, fragment_idx), bonds)

      ! Subtract all proper subsets (size 1 to n-1)
      do subset_size = 1, n - 1
         ! Initialize first combination
         do i = 1, subset_size
            indices(i) = i
         end do

         has_next = .true.
         do while (has_next)
            do i = 1, subset_size
               subset(i) = fragment(indices(i))
            end do
            subset_idx = lookup%find(subset(1:subset_size), subset_size)

            if (subset_idx > 0) then
               ! Subtract this subset's delta Hessian
               delta_hessians(:, :, fragment_idx) = delta_hessians(:, :, fragment_idx) - &
                                                    delta_hessians(:, :, subset_idx)
            end if

            call get_next_combination(indices, subset_size, n, has_next)
         end do
      end do

   end subroutine compute_mbe_hessian

   subroutine compute_mbe_dipole_derivatives(fragment_idx, fragment, lookup, results, delta_dipole_derivs, n, sys_geom, bonds)
      !! Bottom-up computation of n-body dipole derivative correction
      !! Mirrors MBE Hessian logic but for dipole derivatives
      use mqc_result_types, only: calculation_result_t
      use mqc_config_parser, only: bond_t
      use mqc_error, only: error_t
      integer(int64), intent(in) :: fragment_idx
      integer, intent(in) :: fragment(:), n
      type(fragment_lookup_t), intent(in) :: lookup
      type(calculation_result_t), intent(in) :: results(:)
      real(dp), intent(inout) :: delta_dipole_derivs(:, :, :)  !! (3, 3*total_atoms, fragment_count)
      type(system_geometry_t), intent(in) :: sys_geom
      type(bond_t), intent(in), optional :: bonds(:)

      integer :: subset_size, i
      integer :: indices(MAX_MBE_LEVEL), subset(MAX_MBE_LEVEL)  ! Stack arrays to avoid heap contention
      integer(int64) :: subset_idx
      logical :: has_next

      ! Start with the full n-mer dipole derivatives mapped to system coordinates
      call map_fragment_to_system_dipole_derivatives(results(fragment_idx)%dipole_derivatives, fragment, &
                                                     sys_geom, delta_dipole_derivs(:, :, fragment_idx), bonds)

      ! Subtract all proper subsets (size 1 to n-1)
      do subset_size = 1, n - 1
         ! Initialize first combination
         do i = 1, subset_size
            indices(i) = i
         end do

         has_next = .true.
         do while (has_next)
            do i = 1, subset_size
               subset(i) = fragment(indices(i))
            end do
            subset_idx = lookup%find(subset(1:subset_size), subset_size)

            if (subset_idx > 0) then
               ! Subtract this subset's delta dipole derivatives
               delta_dipole_derivs(:, :, fragment_idx) = delta_dipole_derivs(:, :, fragment_idx) - &
                                                         delta_dipole_derivs(:, :, subset_idx)
            end if

            call get_next_combination(indices, subset_size, n, has_next)
         end do
      end do

   end subroutine compute_mbe_dipole_derivatives

   !---------------------------------------------------------------------------
   ! Helper subroutines for reducing code duplication
   !---------------------------------------------------------------------------

   subroutine build_mbe_lookup_table(polymers, fragment_count, max_level, lookup, error)
      !! Build hash table for fast fragment lookups
      use mqc_error, only: error_t
      integer, intent(in) :: polymers(:, :)
      integer(int64), intent(in) :: fragment_count
      integer, intent(in) :: max_level
      type(fragment_lookup_t), intent(inout) :: lookup
      type(error_t), intent(out), optional :: error

      integer(int64) :: i
      integer :: fragment_size
      type(timer_type) :: lookup_timer
      type(error_t) :: insert_error

      call lookup_timer%start()
      call lookup%init(fragment_count)
      do i = 1_int64, fragment_count
         fragment_size = count(polymers(i, :) > 0)
         call lookup%insert(polymers(i, :), fragment_size, i, insert_error)
         if (insert_error%has_error()) then
            if (present(error)) then
               error = insert_error
               call error%add_context("build_mbe_lookup_table")
            end if
            return
         end if
      end do
      call lookup_timer%stop()
      call logger%debug("Time to build lookup table: "//to_char(lookup_timer%get_elapsed_time())//" s")
      call logger%debug("Hash table size: "//to_char(lookup%table_size)// &
                        ", entries: "//to_char(lookup%n_entries))
   end subroutine build_mbe_lookup_table

   subroutine print_mbe_energy_breakdown(sum_by_level, max_level, total_energy)
      !! Print MBE energy breakdown to logger
      real(dp), intent(in) :: sum_by_level(:)
      integer, intent(in) :: max_level
      real(dp), intent(in) :: total_energy

      integer :: nlevel
      character(len=256) :: energy_line

      call logger%info("MBE Energy breakdown:")
      do nlevel = 1, max_level
         if (abs(sum_by_level(nlevel)) > 1e-15_dp) then
            write (energy_line, '(a,i0,a,f20.10)') "  ", nlevel, "-body:  ", sum_by_level(nlevel)
            call logger%info(trim(energy_line))
         end if
      end do
      write (energy_line, '(a,f20.10)') "  Total:   ", total_energy
      call logger%info(trim(energy_line))
   end subroutine print_mbe_energy_breakdown

   subroutine print_mbe_gradient_info(total_gradient, sys_geom, current_log_level)
      !! Print MBE gradient information
      real(dp), intent(in) :: total_gradient(:, :)
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: current_log_level

      integer :: iatom
      character(len=256) :: grad_line

      call logger%info("MBE gradient computation completed")
      call logger%info("  Total gradient norm: "//to_char(sqrt(sum(total_gradient**2))))

      if (current_log_level >= info_level .and. sys_geom%total_atoms < 100) then
         call logger%info(" ")
         call logger%info("Total MBE Gradient (Hartree/Bohr):")
         do iatom = 1, sys_geom%total_atoms
            write (grad_line, '(a,i5,a,3f20.12)') "  Atom ", iatom, ": ", &
               total_gradient(1, iatom), total_gradient(2, iatom), total_gradient(3, iatom)
            call logger%info(trim(grad_line))
         end do
         call logger%info(" ")
      end if
   end subroutine print_mbe_gradient_info

   subroutine compute_mbe(polymers, fragment_count, max_level, results, &
                          mbe_result, sys_geom, bonds, world_comm, json_data)
      !! Compute many-body expansion (MBE) energy with optional gradient, hessian, and dipole
      !!
      !! This is the core routine that handles all MBE computations.
      !! The caller pre-allocates desired components in mbe_result before calling:
      !!   - mbe_result%gradient allocated: compute gradient (requires sys_geom)
      !!   - mbe_result%hessian allocated: compute hessian (requires sys_geom)
      !!   - mbe_result%dipole allocated: compute total dipole moment
      !! If json_data is present, populates it for centralized JSON output
      use mqc_result_types, only: calculation_result_t, mbe_result_t
      use mqc_config_parser, only: bond_t

      ! Required arguments
      integer(int64), intent(in) :: fragment_count
      integer, intent(in) :: polymers(:, :), max_level
      type(calculation_result_t), intent(in) :: results(:)
      type(mbe_result_t), intent(inout) :: mbe_result  !! Pre-allocated by caller

      ! Optional arguments
      type(system_geometry_t), intent(in), optional :: sys_geom  !! Required for gradient/hessian
      type(bond_t), intent(in), optional :: bonds(:)             !! Bond info for H-cap handling
      type(comm_t), intent(in), optional :: world_comm           !! MPI communicator for abort
      type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data

      ! Local variables
      integer(int64) :: i
      integer :: fragment_size, nlevel, current_log_level, hess_dim
      real(dp), allocatable :: sum_by_level(:), delta_energies(:), energies(:)
      real(dp), allocatable :: delta_gradients(:, :, :), delta_hessians(:, :, :)
      real(dp), allocatable :: delta_dipoles(:, :)  !! (3, fragment_count)
      real(dp), allocatable :: delta_dipole_derivs(:, :, :)  !! (3, 3*total_atoms, fragment_count)
      real(dp), allocatable :: ir_intensities(:)  !! IR intensities in km/mol
      real(dp) :: delta_E
      logical :: do_detailed_print, compute_grad, compute_hess, compute_dipole, compute_dipole_derivs
      type(fragment_lookup_t) :: lookup

      ! Determine what to compute based on allocated components in mbe_result
      compute_grad = allocated(mbe_result%gradient)
      compute_hess = allocated(mbe_result%hessian)
      compute_dipole = allocated(mbe_result%dipole)
      compute_dipole_derivs = .false.  ! Will be set true if all fragments have dipole_derivatives

      ! Validate inputs for gradient computation
      if (compute_grad) then
         do i = 1_int64, fragment_count
            if (.not. results(i)%has_gradient) then
               call logger%error("Fragment "//to_char(i)//" does not have gradient!")
               if (present(world_comm)) then
                  call abort_comm(world_comm, 1)
               else
                  error stop "Missing gradient in compute_mbe_internal"
               end if
            end if
         end do
      end if

      ! Validate inputs for hessian computation
      if (compute_hess) then
         do i = 1_int64, fragment_count
            if (.not. results(i)%has_hessian) then
               call logger%error("Fragment "//to_char(i)//" does not have Hessian!")
               if (present(world_comm)) then
                  call abort_comm(world_comm, 1)
               else
                  error stop "Missing Hessian in compute_mbe_internal"
               end if
            end if
         end do
         hess_dim = 3*sys_geom%total_atoms
      end if

      ! Validate inputs for dipole computation
      if (compute_dipole) then
         do i = 1_int64, fragment_count
            if (.not. results(i)%has_dipole) then
               call logger%warning("Fragment "//to_char(i)//" does not have dipole - skipping dipole MBE")
               compute_dipole = .false.
               exit
            end if
         end do
      end if

      ! Check if dipole derivatives are available (for IR intensities)
      if (compute_hess) then
         compute_dipole_derivs = .true.
         do i = 1_int64, fragment_count
            if (.not. results(i)%has_dipole_derivatives) then
               compute_dipole_derivs = .false.
               exit
            end if
         end do
      end if

      ! Get logger configuration
      call logger%configuration(level=current_log_level)
      do_detailed_print = (current_log_level >= verbose_level)

      ! Allocate energy arrays (always needed)
      allocate (sum_by_level(max_level))
      allocate (delta_energies(fragment_count))
      allocate (energies(fragment_count))
      sum_by_level = 0.0_dp
      delta_energies = 0.0_dp

      ! Extract total energies from results
      do i = 1_int64, fragment_count
         energies(i) = results(i)%energy%total()
      end do

      ! Allocate gradient delta arrays if needed
      if (compute_grad) then
         allocate (delta_gradients(3, sys_geom%total_atoms, fragment_count))
         delta_gradients = 0.0_dp
         mbe_result%gradient = 0.0_dp
      end if

      ! Allocate hessian delta arrays if needed
      if (compute_hess) then
         allocate (delta_hessians(hess_dim, hess_dim, fragment_count))
         delta_hessians = 0.0_dp
         mbe_result%hessian = 0.0_dp
      end if

      ! Allocate dipole delta arrays if needed
      if (compute_dipole) then
         allocate (delta_dipoles(3, fragment_count))
         delta_dipoles = 0.0_dp
         mbe_result%dipole = 0.0_dp
      end if

      ! Allocate dipole derivative delta arrays if needed (for IR intensities)
      if (compute_dipole_derivs) then
         allocate (delta_dipole_derivs(3, hess_dim, fragment_count))
         delta_dipole_derivs = 0.0_dp
         allocate (mbe_result%dipole_derivatives(3, hess_dim))
         mbe_result%dipole_derivatives = 0.0_dp
      end if

      ! Build hash table for fast fragment lookups
      block
         use mqc_error, only: error_t
         type(error_t) :: lookup_error
         call build_mbe_lookup_table(polymers, fragment_count, max_level, lookup, lookup_error)
         if (lookup_error%has_error()) then
            call logger%error("Failed to build lookup table: "//lookup_error%get_message())
            if (present(world_comm)) then
               call abort_comm(world_comm, 1)
            else
               error stop "Failed to build lookup table"
            end if
         end if
      end block

      ! Bottom-up computation: process fragments by size (1-body, then 2-body, etc.)
      ! This makes the algorithm independent of input fragment order
      ! We process by n-mer level to ensure all subsets are computed before they're needed
      do nlevel = 1, max_level
         do i = 1_int64, fragment_count
            fragment_size = count(polymers(i, :) > 0)

            ! Only process fragments of the current nlevel
            if (fragment_size /= nlevel) cycle

            if (fragment_size == 1) then
               ! 1-body: delta = value (no subsets to subtract)
               delta_energies(i) = energies(i)
               sum_by_level(1) = sum_by_level(1) + delta_energies(i)

               if (compute_grad) then
                  call map_fragment_to_system_gradient(results(i)%gradient, &
                                      polymers(i, 1:fragment_size), sys_geom, delta_gradients(:, :, i), bonds, world_comm)
               end if

               if (compute_hess) then
                  call map_fragment_to_system_hessian(results(i)%hessian, &
                                                   polymers(i, 1:fragment_size), sys_geom, delta_hessians(:, :, i), bonds)
               end if

               if (compute_dipole) then
                  ! For 1-body, delta dipole is just the fragment dipole
                  delta_dipoles(:, i) = results(i)%dipole
               end if

               if (compute_dipole_derivs) then
                  ! For 1-body, delta dipole derivatives are just the fragment values mapped to system
                  call map_fragment_to_system_dipole_derivatives(results(i)%dipole_derivatives, &
                                              polymers(i, 1:fragment_size), sys_geom, delta_dipole_derivs(:, :, i), bonds)
               end if

            else if (fragment_size >= 2 .and. fragment_size <= max_level) then
               ! n-body: delta = value - sum(all subset deltas)
               delta_E = compute_mbe_delta(i, polymers(i, 1:fragment_size), lookup, &
                                           energies, delta_energies, fragment_size, world_comm)
               delta_energies(i) = delta_E
               sum_by_level(fragment_size) = sum_by_level(fragment_size) + delta_E

               if (compute_grad) then
                  call compute_mbe_gradient(i, polymers(i, 1:fragment_size), lookup, &
                                            results, delta_gradients, fragment_size, sys_geom, bonds, world_comm)
               end if

               if (compute_hess) then
                  call compute_mbe_hessian(i, polymers(i, 1:fragment_size), lookup, &
                                           results, delta_hessians, fragment_size, sys_geom, bonds)
               end if

               if (compute_dipole) then
                  call compute_mbe_dipole(i, polymers(i, 1:fragment_size), lookup, &
                                          results, delta_dipoles, fragment_size, world_comm)
               end if

               if (compute_dipole_derivs) then
                  call compute_mbe_dipole_derivatives(i, polymers(i, 1:fragment_size), lookup, &
                                                      results, delta_dipole_derivs, fragment_size, sys_geom, bonds)
               end if
            end if
         end do
      end do

      ! Clean up lookup table
      call lookup%destroy()

      ! Compute totals and set status flags
      mbe_result%total_energy = sum(sum_by_level)
      mbe_result%has_energy = .true.

      if (compute_grad) then
         do i = 1_int64, fragment_count
            fragment_size = count(polymers(i, :) > 0)
            if (fragment_size <= max_level) then
               mbe_result%gradient = mbe_result%gradient + delta_gradients(:, :, i)
            end if
         end do
         mbe_result%has_gradient = .true.
      end if

      if (compute_hess) then
         do i = 1_int64, fragment_count
            fragment_size = count(polymers(i, :) > 0)
            if (fragment_size <= max_level) then
               mbe_result%hessian = mbe_result%hessian + delta_hessians(:, :, i)
            end if
         end do
         mbe_result%has_hessian = .true.
      end if

      if (compute_dipole) then
         do i = 1_int64, fragment_count
            fragment_size = count(polymers(i, :) > 0)
            if (fragment_size <= max_level) then
               mbe_result%dipole = mbe_result%dipole + delta_dipoles(:, i)
            end if
         end do
         mbe_result%has_dipole = .true.
      end if

      if (compute_dipole_derivs) then
         do i = 1_int64, fragment_count
            fragment_size = count(polymers(i, :) > 0)
            if (fragment_size <= max_level) then
               mbe_result%dipole_derivatives = mbe_result%dipole_derivatives + delta_dipole_derivs(:, :, i)
            end if
         end do
         mbe_result%has_dipole_derivatives = .true.
      end if

      ! Print energy breakdown (always)
      call print_mbe_energy_breakdown(sum_by_level, max_level, mbe_result%total_energy)

      ! Print gradient info if computed
      if (compute_grad) then
         call print_mbe_gradient_info(mbe_result%gradient, sys_geom, current_log_level)
      end if

      ! Print hessian info if computed
      if (compute_hess) then
         call logger%info("MBE Hessian computation completed")
         call logger%info("  Total Hessian Frobenius norm: "//to_char(sqrt(sum(mbe_result%hessian**2))))

         ! Compute and print full vibrational analysis with thermochemistry
         block
            real(dp), allocatable :: frequencies(:), reduced_masses(:), force_constants(:)
            real(dp), allocatable :: cart_disp(:, :), fc_mdyne(:)
            type(thermochemistry_result_t) :: thermo_result
            integer :: n_at, n_modes

            call logger%info("  Computing vibrational analysis (projecting trans/rot modes)...")
            if (compute_dipole_derivs) then
               call compute_vibrational_analysis(mbe_result%hessian, sys_geom%element_numbers, frequencies, &
                                                 reduced_masses, force_constants, cart_disp, &
                                                 coordinates=sys_geom%coordinates, &
                                                 project_trans_rot=.true., &
                                                 force_constants_mdyne=fc_mdyne, &
                                                 dipole_derivatives=mbe_result%dipole_derivatives, &
                                                 ir_intensities=ir_intensities)
            else
               call compute_vibrational_analysis(mbe_result%hessian, sys_geom%element_numbers, frequencies, &
                                                 reduced_masses, force_constants, cart_disp, &
                                                 coordinates=sys_geom%coordinates, &
                                                 project_trans_rot=.true., &
                                                 force_constants_mdyne=fc_mdyne)
            end if

            if (allocated(frequencies)) then
               ! Compute thermochemistry
               n_at = size(sys_geom%element_numbers)
               n_modes = size(frequencies)
               call compute_thermochemistry(sys_geom%coordinates, sys_geom%element_numbers, &
                                            frequencies, n_at, n_modes, thermo_result)

               ! Print vibrational analysis to log
               if (allocated(ir_intensities)) then
                  call print_vibrational_analysis(frequencies, reduced_masses, force_constants, &
                                                  cart_disp, sys_geom%element_numbers, &
                                                  force_constants_mdyne=fc_mdyne, &
                                                  ir_intensities=ir_intensities, &
                                                  coordinates=sys_geom%coordinates, &
                                                  electronic_energy=mbe_result%total_energy)
               else
                  call print_vibrational_analysis(frequencies, reduced_masses, force_constants, &
                                                  cart_disp, sys_geom%element_numbers, &
                                                  force_constants_mdyne=fc_mdyne, &
                                                  coordinates=sys_geom%coordinates, &
                                                  electronic_energy=mbe_result%total_energy)
               end if

               ! Populate json_data for vibrational output if present
               if (present(json_data)) then
                  json_data%output_mode = OUTPUT_MODE_MBE
                  json_data%total_energy = mbe_result%total_energy
                  json_data%has_energy = mbe_result%has_energy
                  json_data%has_vibrational = .true.

                  ! Copy vibrational data
                  allocate (json_data%frequencies(n_modes))
                  allocate (json_data%reduced_masses(n_modes))
                  allocate (json_data%force_constants(n_modes))
                  json_data%frequencies = frequencies
                  json_data%reduced_masses = reduced_masses
                  json_data%force_constants = fc_mdyne
                  json_data%thermo = thermo_result

                  if (allocated(ir_intensities)) then
                     allocate (json_data%ir_intensities(n_modes))
                     json_data%ir_intensities = ir_intensities
                     json_data%has_ir_intensities = .true.
                  end if

                  ! Copy dipole if available
                  if (mbe_result%has_dipole) then
                     allocate (json_data%dipole(3))
                     json_data%dipole = mbe_result%dipole
                     json_data%has_dipole = .true.
                  end if

                  ! Copy gradient if available
                  if (mbe_result%has_gradient) then
                     allocate (json_data%gradient, source=mbe_result%gradient)
                     json_data%has_gradient = .true.
                  end if

                  ! Copy hessian if available
                  if (mbe_result%has_hessian) then
                     allocate (json_data%hessian, source=mbe_result%hessian)
                     json_data%has_hessian = .true.
                  end if
               end if

               if (allocated(ir_intensities)) deallocate (ir_intensities)
               deallocate (frequencies, reduced_masses, force_constants, cart_disp, fc_mdyne)
            end if
         end block
      end if

      ! Print dipole info if computed
      if (compute_dipole) then
         block
            character(len=256) :: dipole_line
            real(dp) :: dipole_magnitude
            dipole_magnitude = norm2(mbe_result%dipole)*2.541746_dp  ! Convert e*Bohr to Debye
            call logger%info("MBE Dipole moment:")
            write (dipole_line, '(a,3f15.8)') "  Dipole (e*Bohr): ", mbe_result%dipole
            call logger%info(trim(dipole_line))
            write (dipole_line, '(a,f15.8)') "  Dipole magnitude (Debye): ", dipole_magnitude
            call logger%info(trim(dipole_line))
         end block
      end if

      ! Print detailed breakdown if requested
      if (do_detailed_print) then
         call print_detailed_breakdown(polymers, fragment_count, max_level, energies, delta_energies)
      end if

      ! Populate json_data for non-Hessian case if present
      ! (Hessian case already handled above in the vibrational block)
      if (present(json_data) .and. .not. compute_hess) then
         json_data%output_mode = OUTPUT_MODE_MBE
         json_data%total_energy = mbe_result%total_energy
         json_data%has_energy = mbe_result%has_energy
         json_data%max_level = max_level
         json_data%fragment_count = fragment_count

         ! Copy fragment breakdown data
         allocate (json_data%polymers(fragment_count, max_level))
         json_data%polymers = polymers(1:fragment_count, 1:max_level)

         allocate (json_data%fragment_energies(fragment_count))
         json_data%fragment_energies = energies

         allocate (json_data%delta_energies(fragment_count))
         json_data%delta_energies = delta_energies

         allocate (json_data%sum_by_level(max_level))
         json_data%sum_by_level = sum_by_level

         ! Copy fragment distances if available
         allocate (json_data%fragment_distances(fragment_count))
         do i = 1_int64, fragment_count
            json_data%fragment_distances(i) = results(i)%distance
         end do

         ! Copy dipole if available
         if (mbe_result%has_dipole) then
            allocate (json_data%dipole(3))
            json_data%dipole = mbe_result%dipole
            json_data%has_dipole = .true.
         end if

         ! Copy gradient if available
         if (mbe_result%has_gradient) then
            allocate (json_data%gradient, source=mbe_result%gradient)
            json_data%has_gradient = .true.
         end if
      end if

      ! Cleanup
      deallocate (sum_by_level, delta_energies, energies)
      if (allocated(delta_gradients)) deallocate (delta_gradients)
      if (allocated(delta_hessians)) deallocate (delta_hessians)
      if (allocated(delta_dipoles)) deallocate (delta_dipoles)
      if (allocated(delta_dipole_derivs)) deallocate (delta_dipole_derivs)

   end subroutine compute_mbe

   !---------------------------------------------------------------------------
   ! GMBE Helper subroutines for reducing code duplication
   !---------------------------------------------------------------------------

   subroutine print_gmbe_energy_breakdown(monomer_energy, n_monomers, level_energies, level_counts, &
                                          max_level, total_energy, has_intersections)
      !! Print GMBE energy breakdown using inclusion-exclusion principle
      real(dp), intent(in) :: monomer_energy
      integer, intent(in) :: n_monomers
      real(dp), intent(in), optional :: level_energies(:)
      integer, intent(in), optional :: level_counts(:)
      integer, intent(in) :: max_level
      real(dp), intent(in) :: total_energy
      logical, intent(in) :: has_intersections

      integer :: k
      real(dp) :: sign_factor
      character(len=256) :: line

      if (has_intersections) then
         call logger%info("GMBE Energy breakdown (Inclusion-Exclusion Principle):")
         write (line, '(a,i0,a,f20.10)') "  Monomers (", n_monomers, "):  ", monomer_energy
         call logger%info(trim(line))

         do k = 2, max_level
            if (level_counts(k) > 0) then
               sign_factor = real((-1)**(k + 1), dp)
               if (sign_factor > 0.0_dp) then
                  write (line, '(a,i0,a,i0,a,f20.10)') "  ", k, "-way ∩ (", level_counts(k), "):  +", level_energies(k)
               else
                  write (line, '(a,i0,a,i0,a,f20.10)') "  ", k, "-way ∩ (", level_counts(k), "):  ", level_energies(k)
               end if
               call logger%info(trim(line))
            end if
         end do

         write (line, '(a,f20.10)') "  Total GMBE:      ", total_energy
         call logger%info(trim(line))
      else
         call logger%info("GMBE Energy breakdown:")
         write (line, '(a,i0,a,f20.10)') "  Monomers (", n_monomers, "): ", monomer_energy
         call logger%info(trim(line))
         write (line, '(a,f20.10)') "  Total GMBE:  ", total_energy
         call logger%info(trim(line))
      end if
   end subroutine print_gmbe_energy_breakdown

   subroutine print_gmbe_gradient_info(total_gradient, sys_geom, current_log_level)
      !! Print GMBE gradient information
      real(dp), intent(in) :: total_gradient(:, :)
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: current_log_level

      integer :: iatom
      character(len=256) :: grad_line

      call logger%info("GMBE gradient computation completed")
      call logger%info("  Total gradient norm: "//to_char(sqrt(sum(total_gradient**2))))

      if (current_log_level >= info_level .and. sys_geom%total_atoms < 100) then
         call logger%info(" ")
         call logger%info("Total GMBE Gradient (Hartree/Bohr):")
         do iatom = 1, sys_geom%total_atoms
            write (grad_line, '(a,i5,a,3f20.12)') "  Atom ", iatom, ": ", &
               total_gradient(1, iatom), total_gradient(2, iatom), total_gradient(3, iatom)
            call logger%info(trim(grad_line))
         end do
         call logger%info(" ")
      end if
   end subroutine print_gmbe_gradient_info

   subroutine compute_gmbe(monomers, n_monomers, monomer_results, &
                           n_intersections, intersection_results, &
                           intersection_sets, intersection_levels, &
                           total_energy, &
                           sys_geom, total_gradient, total_hessian, bonds, world_comm)
      !! Compute generalized many-body expansion (GMBE) energy with optional gradient and/or hessian
      !!
      !! This is the core routine that handles all GMBE computations using
      !! the inclusion-exclusion principle for overlapping fragments.
      !! Optional arguments control what derivatives are computed:
      !!   - If sys_geom and total_gradient are present: compute gradient
      !!   - If total_hessian is also present: compute hessian
      use mqc_result_types, only: calculation_result_t
      use mqc_physical_fragment, only: build_fragment_from_indices, build_fragment_from_atom_list, &
                                       redistribute_cap_gradients, redistribute_cap_hessian
      use mqc_gmbe_utils, only: find_fragment_intersection
      use mqc_config_parser, only: bond_t
      use mqc_error, only: error_t

      ! Required arguments
      integer, intent(in) :: monomers(:)
      integer, intent(in) :: n_monomers
      type(calculation_result_t), intent(in) :: monomer_results(:)
      integer, intent(in) :: n_intersections
      type(calculation_result_t), intent(in), optional :: intersection_results(:)
      integer, intent(in), optional :: intersection_sets(:, :)
      integer, intent(in), optional :: intersection_levels(:)
      real(dp), intent(out) :: total_energy

      ! Optional arguments for gradient computation
      type(system_geometry_t), intent(in), optional :: sys_geom
      real(dp), intent(out), optional :: total_gradient(:, :)

      ! Optional arguments for hessian computation
      real(dp), intent(out), optional :: total_hessian(:, :)

      ! Optional bond information for hydrogen cap handling
      type(bond_t), intent(in), optional :: bonds(:)

      ! Optional MPI communicator for abort
      type(comm_t), intent(in), optional :: world_comm

      ! Local variables
      integer :: i, k, max_level, current_log_level, hess_dim
      real(dp) :: monomer_energy, sign_factor
      real(dp), allocatable :: level_energies(:)
      integer, allocatable :: level_counts(:)
      type(physical_fragment_t) :: fragment
      type(error_t) :: error
      integer, allocatable :: monomer_idx(:)
      logical :: compute_grad, compute_hess, has_intersections

      ! Determine what to compute based on optional arguments
      compute_grad = present(sys_geom) .and. present(total_gradient)
      compute_hess = compute_grad .and. present(total_hessian)
      has_intersections = n_intersections > 0 .and. present(intersection_results) .and. &
                          present(intersection_sets) .and. present(intersection_levels)

      ! Initialize outputs
      if (compute_grad) then
         total_gradient = 0.0_dp
      end if
      if (compute_hess) then
         total_hessian = 0.0_dp
         hess_dim = 3*sys_geom%total_atoms
      end if

      ! Sum monomer energies (and gradients/hessians if requested)
      monomer_energy = 0.0_dp
      do i = 1, n_monomers
         monomer_energy = monomer_energy + monomer_results(i)%energy%total()

         ! Accumulate monomer derivatives if requested
         if (compute_grad .and. monomer_results(i)%has_gradient) then
            allocate (monomer_idx(1))
            monomer_idx(1) = monomers(i)
            call build_fragment_from_indices(sys_geom, monomer_idx, fragment, error, bonds)
            if (error%has_error()) then
               call logger%error(error%get_full_trace())
               if (present(world_comm)) then
                  call abort_comm(world_comm, 1)
               else
                  error stop "Failed to build monomer fragment in GMBE"
               end if
            end if
            call redistribute_cap_gradients(fragment, monomer_results(i)%gradient, total_gradient)
            if (compute_hess .and. monomer_results(i)%has_hessian) then
               call redistribute_cap_hessian(fragment, monomer_results(i)%hessian, total_hessian)
            end if
            call fragment%destroy()
            deallocate (monomer_idx)
         end if
      end do

      ! Start with monomer contribution
      total_energy = monomer_energy

      ! Handle intersections with inclusion-exclusion
      if (has_intersections) then
         max_level = maxval(intersection_levels)

         allocate (level_energies(2:max_level))
         allocate (level_counts(2:max_level))
         level_energies = 0.0_dp
         level_counts = 0

         ! Process intersection energies and derivatives
         do i = 1, n_intersections
            k = intersection_levels(i)
            sign_factor = real((-1)**(k + 1), dp)
            level_energies(k) = level_energies(k) + intersection_results(i)%energy%total()
            level_counts(k) = level_counts(k) + 1

            ! Handle intersection derivatives if requested
            if (compute_grad .and. (intersection_results(i)%has_gradient .or. &
                                    (compute_hess .and. intersection_results(i)%has_hessian))) then
               call process_intersection_derivatives(i, k, sign_factor, intersection_results, &
                                                     intersection_sets, sys_geom, total_gradient, total_hessian, bonds, &
                                                     compute_grad, compute_hess, hess_dim, world_comm)
            end if
         end do

         ! Apply inclusion-exclusion to energy
         do k = 2, max_level
            if (level_counts(k) > 0) then
               sign_factor = real((-1)**(k + 1), dp)
               total_energy = total_energy + sign_factor*level_energies(k)
            end if
         end do

         ! Print energy breakdown
         call print_gmbe_energy_breakdown(monomer_energy, n_monomers, level_energies, level_counts, &
                                          max_level, total_energy, .true.)

         ! Print debug info for intersections
         call print_gmbe_intersection_debug(n_intersections, n_monomers, intersection_sets, &
                                            intersection_levels, intersection_results)

         deallocate (level_energies, level_counts)
      else
         ! No intersections - just print monomer sum
         call print_gmbe_energy_breakdown(monomer_energy, n_monomers, level_energies, level_counts, &
                                          0, total_energy, .false.)
      end if

      ! Print gradient/hessian info
      if (compute_grad) then
         call logger%configuration(level=current_log_level)
         call print_gmbe_gradient_info(total_gradient, sys_geom, current_log_level)
      end if

      if (compute_hess) then
         call logger%info("GMBE Hessian computation completed")
         call logger%info("  Total Hessian Frobenius norm: "//to_char(sqrt(sum(total_hessian**2))))

         ! Compute and print full vibrational analysis
         block
            real(dp), allocatable :: frequencies(:), reduced_masses(:), force_constants(:)
            real(dp), allocatable :: cart_disp(:, :), fc_mdyne(:)

            call logger%info("  Computing vibrational analysis (projecting trans/rot modes)...")
            call compute_vibrational_analysis(total_hessian, sys_geom%element_numbers, frequencies, &
                                              reduced_masses, force_constants, cart_disp, &
                                              coordinates=sys_geom%coordinates, &
                                              project_trans_rot=.true., &
                                              force_constants_mdyne=fc_mdyne)

            if (allocated(frequencies)) then
               call print_vibrational_analysis(frequencies, reduced_masses, force_constants, &
                                               cart_disp, sys_geom%element_numbers, &
                                               force_constants_mdyne=fc_mdyne, &
                                               coordinates=sys_geom%coordinates, &
                                               electronic_energy=total_energy)
               deallocate (frequencies, reduced_masses, force_constants, cart_disp, fc_mdyne)
            end if
         end block
      end if

   end subroutine compute_gmbe

   subroutine process_intersection_derivatives(inter_idx, k, sign_factor, intersection_results, &
                                               intersection_sets, sys_geom, total_gradient, &
                                               total_hessian, bonds, compute_grad, compute_hess, hess_dim, world_comm)
      !! Process derivatives for a single intersection fragment
      use mqc_result_types, only: calculation_result_t
      use mqc_physical_fragment, only: build_fragment_from_atom_list, &
                                       redistribute_cap_gradients, redistribute_cap_hessian
      use mqc_gmbe_utils, only: find_fragment_intersection
      use mqc_config_parser, only: bond_t
      use mqc_error, only: error_t

      integer, intent(in) :: inter_idx, k
      real(dp), intent(in) :: sign_factor
      type(calculation_result_t), intent(in) :: intersection_results(:)
      integer, intent(in) :: intersection_sets(:, :)
      type(system_geometry_t), intent(in) :: sys_geom
      real(dp), intent(inout) :: total_gradient(:, :)
      real(dp), intent(inout), optional :: total_hessian(:, :)
      type(bond_t), intent(in), optional :: bonds(:)
      logical, intent(in) :: compute_grad, compute_hess
      integer, intent(in) :: hess_dim
      type(comm_t), intent(in), optional :: world_comm

      integer :: j, current_n, next_n, n_intersect
      integer, allocatable :: current_atoms(:), next_atoms(:), intersect_atoms(:)
      real(dp), allocatable :: term_gradient(:, :), term_hessian(:, :)
      logical :: has_intersection
      type(physical_fragment_t) :: inter_fragment
      type(error_t) :: inter_error

      ! Build the intersection atom list
      call get_monomer_atom_list(sys_geom, intersection_sets(1, inter_idx), current_atoms, current_n)

      if (current_n > 0) then
         do j = 2, k
            call get_monomer_atom_list(sys_geom, intersection_sets(j, inter_idx), next_atoms, next_n)
            has_intersection = find_fragment_intersection(current_atoms, current_n, &
                                                          next_atoms, next_n, &
                                                          intersect_atoms, n_intersect)
            deallocate (current_atoms, next_atoms)

            if (.not. has_intersection) then
               current_n = 0
               exit
            end if

            current_n = n_intersect
            allocate (current_atoms(current_n))
            current_atoms = intersect_atoms
            deallocate (intersect_atoms)
         end do
      end if

      if (current_n > 0) then
         call build_fragment_from_atom_list(sys_geom, current_atoms, current_n, &
                                            inter_fragment, inter_error, bonds)
         if (inter_error%has_error()) then
            call logger%error(inter_error%get_full_trace())
            if (present(world_comm)) then
               call abort_comm(world_comm, 1)
            else
               error stop "Failed to build intersection fragment in GMBE"
            end if
         end if

         if (compute_grad .and. intersection_results(inter_idx)%has_gradient) then
            allocate (term_gradient(3, sys_geom%total_atoms))
            term_gradient = 0.0_dp
            call redistribute_cap_gradients(inter_fragment, intersection_results(inter_idx)%gradient, &
                                            term_gradient)
            total_gradient = total_gradient + sign_factor*term_gradient
            deallocate (term_gradient)
         end if

         if (compute_hess .and. intersection_results(inter_idx)%has_hessian) then
            allocate (term_hessian(hess_dim, hess_dim))
            term_hessian = 0.0_dp
            call redistribute_cap_hessian(inter_fragment, intersection_results(inter_idx)%hessian, &
                                          term_hessian)
            total_hessian = total_hessian + sign_factor*term_hessian
            deallocate (term_hessian)
         end if

         call inter_fragment%destroy()
      else
         call logger%warning("GMBE intersection has no atoms; skipping derivatives")
      end if

      if (allocated(current_atoms)) deallocate (current_atoms)
   end subroutine process_intersection_derivatives

   subroutine print_gmbe_intersection_debug(n_intersections, n_monomers, intersection_sets, &
                                            intersection_levels, intersection_results)
      !! Print debug information about GMBE intersections
      use mqc_result_types, only: calculation_result_t

      integer, intent(in) :: n_intersections, n_monomers
      integer, intent(in) :: intersection_sets(:, :)
      integer, intent(in) :: intersection_levels(:)
      type(calculation_result_t), intent(in) :: intersection_results(:)

      integer :: i, j, set_size
      real(dp) :: sign_factor
      character(len=512) :: detail_line
      character(len=256) :: set_str

      if (n_intersections > 0) then
         call logger%debug("GMBE intersection details:")
         do i = 1, n_intersections
            set_str = "("
            set_size = 0
            do j = 1, n_monomers
               if (intersection_sets(j, i) > 0) then
                  if (set_size > 0) set_str = trim(set_str)//","
                  write (set_str, '(a,i0)') trim(set_str), intersection_sets(j, i)
                  set_size = set_size + 1
               end if
            end do
            set_str = trim(set_str)//")"

            sign_factor = real((-1)**(intersection_levels(i) + 1), dp)
            write (detail_line, '(a,i0,a,i0,a,a,a,f16.8)') &
               "  Intersection ", i, ": level=", intersection_levels(i), &
               " fragments=", trim(set_str), " energy=", intersection_results(i)%energy%total()
            call logger%debug(trim(detail_line))
         end do
      end if
   end subroutine print_gmbe_intersection_debug

   subroutine get_monomer_atom_list(sys_geom, monomer_idx, atom_list, n_atoms)
      !! Build 0-indexed atom list for a monomer, handling fixed or variable-sized fragments.
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: monomer_idx
      integer, allocatable, intent(out) :: atom_list(:)
      integer, intent(out) :: n_atoms

      integer :: i, base_idx

      if (allocated(sys_geom%fragment_atoms)) then
         n_atoms = sys_geom%fragment_sizes(monomer_idx)
         if (n_atoms > 0) then
            allocate (atom_list(n_atoms))
            atom_list = sys_geom%fragment_atoms(1:n_atoms, monomer_idx)
         else
            allocate (atom_list(0))
         end if
      else
         n_atoms = sys_geom%atoms_per_monomer
         if (n_atoms > 0) then
            allocate (atom_list(n_atoms))
            base_idx = (monomer_idx - 1)*sys_geom%atoms_per_monomer
            do i = 1, n_atoms
               atom_list(i) = base_idx + (i - 1)
            end do
         else
            allocate (atom_list(0))
         end if
      end if
   end subroutine get_monomer_atom_list

end module mqc_mbe
