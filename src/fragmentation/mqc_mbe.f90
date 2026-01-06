!! Many-Body Expansion (MBE) calculation module
module mqc_mbe
   !! Implements hierarchical many-body expansion for fragment-based quantum chemistry
   !! calculations with MPI parallelization and energy/gradient computation.
   use pic_types, only: int32, int64, dp
   use pic_timer, only: timer_type
   use pic_mpi_lib, only: comm_t, send, recv, iprobe, MPI_Status, MPI_ANY_SOURCE, MPI_ANY_TAG
   use pic_logger, only: logger => global_logger, verbose_level, debug_level, info_level
   use pic_io, only: to_char
   use mqc_mbe_io, only: print_detailed_breakdown, print_detailed_breakdown_json
   use mqc_mpi_tags, only: TAG_WORKER_REQUEST, TAG_WORKER_FRAGMENT, TAG_WORKER_FINISH, &
                           TAG_WORKER_SCALAR_RESULT, &
                           TAG_NODE_REQUEST, TAG_NODE_FRAGMENT, TAG_NODE_FINISH, &
                           TAG_NODE_SCALAR_RESULT
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, build_fragment_from_indices, to_angstrom
   use mqc_frag_utils, only: get_next_combination, fragment_lookup_t

   implicit none
   private

   ! Public interface
   public :: compute_mbe_energy, compute_mbe_energy_gradient, compute_mbe_energy_gradient_hessian
   public :: compute_gmbe_energy  !! GMBE energy with intersection correction
   public :: compute_gmbe_energy_gradient  !! GMBE energy and gradient with intersection correction
   public :: compute_gmbe_energy_gradient_hessian  !! GMBE energy, gradient, and Hessian with intersection correction

contains

   subroutine compute_mbe_energy(polymers, fragment_count, max_level, results, total_energy)
      !! Compute the many-body expansion (MBE) energy
      !! Total = sum(E(i)) + sum(deltaE(ij)) + sum(deltaE(ijk)) + ...
      !! General n-body correction:
      !! deltaE(i1,i2,...,in) = E(i1,i2,...,in) - sum of all lower-order terms
      !! Uses int64 for fragment_count to handle large fragment counts that overflow int32.
      !! Detailed breakdown is printed only if logger level is verbose or higher.
      use mqc_result_types, only: calculation_result_t
      integer(int64), intent(in) :: fragment_count
      integer, intent(in) :: polymers(:, :), max_level
      type(calculation_result_t), intent(in) :: results(:)
      real(dp), intent(out) :: total_energy

      integer(int64) :: i
      integer :: fragment_size, nlevel, current_log_level
      real(dp), allocatable :: sum_by_level(:), delta_energies(:), energies(:)
      real(dp) :: delta_E
      logical :: do_detailed_print
      type(fragment_lookup_t) :: lookup
      type(timer_type) :: lookup_timer

      call logger%configuration(level=current_log_level)
      do_detailed_print = (current_log_level >= verbose_level)

      allocate (sum_by_level(max_level))
      allocate (delta_energies(fragment_count))
      allocate (energies(fragment_count))
      sum_by_level = 0.0_dp
      delta_energies = 0.0_dp

      ! Extract total energies from results
      do i = 1_int64, fragment_count
         energies(i) = results(i)%energy%total()
      end do

      ! Build hash table for fast fragment lookups
      call lookup_timer%start()
      call lookup%init(fragment_count)
      do i = 1_int64, fragment_count
         fragment_size = count(polymers(i, :) > 0)
         call lookup%insert(polymers(i, :), fragment_size, i)
      end do
      call lookup_timer%stop()
      call logger%debug("Time to build lookup table: "//to_char(lookup_timer%get_elapsed_time())//" s")
      call logger%debug("Hash table size: "//to_char(lookup%table_size)// &
                        ", entries: "//to_char(lookup%n_entries))

      ! Bottom-up computation: process fragments by size (1-body, then 2-body, then 3-body, etc.)
      ! This makes the algorithm independent of input fragment order
      ! We process by n-mer level to ensure all subsets are computed before they're needed
      do nlevel = 1, max_level
         do i = 1_int64, fragment_count
            fragment_size = count(polymers(i, :) > 0)

            ! Only process fragments of the current nlevel
            if (fragment_size /= nlevel) cycle

            if (fragment_size == 1) then
               ! 1-body: deltaE = E (no subsets to subtract)
               delta_energies(i) = energies(i)
               sum_by_level(1) = sum_by_level(1) + delta_energies(i)
            else if (fragment_size >= 2 .and. fragment_size <= max_level) then
               ! n-body: deltaE = E - sum(all subset deltaEs)
               ! All subsets have already been computed in previous nlevel iterations
               delta_E = compute_mbe(i, polymers(i, 1:fragment_size), lookup, &
                                     energies, delta_energies, fragment_size)
               delta_energies(i) = delta_E
               sum_by_level(fragment_size) = sum_by_level(fragment_size) + delta_E
            end if
         end do
      end do

      ! Clean up lookup table
      call lookup%destroy()

      total_energy = sum(sum_by_level)

      ! Print text summary to console
      call logger%info("MBE Energy breakdown:")
      do nlevel = 1, max_level
         if (abs(sum_by_level(nlevel)) > 1e-15_dp) then
            block
               character(len=256) :: energy_line
               write (energy_line, '(a,i0,a,f20.10)') "  ", nlevel, "-body:  ", sum_by_level(nlevel)
               call logger%info(trim(energy_line))
            end block
         end if
      end do
      block
         character(len=256) :: total_line
         write (total_line, '(a,f20.10)') "  Total:   ", total_energy
         call logger%info(trim(total_line))
      end block

      ! Print detailed breakdown if requested
      if (do_detailed_print) then
         call print_detailed_breakdown(polymers, fragment_count, max_level, energies, delta_energies)
      end if

      ! Always write JSON file for machine-readable output
      call print_detailed_breakdown_json(polymers, fragment_count, max_level, energies, delta_energies, &
                                         sum_by_level, total_energy, results=results)

      deallocate (sum_by_level, delta_energies, energies)

   end subroutine compute_mbe_energy

   function compute_mbe(fragment_idx, fragment, lookup, energies, delta_energies, n) result(delta_E)
      !! Bottom-up computation of n-body correction (non-recursive, uses pre-computed subset deltas)
      !! deltaE(i1,i2,...,in) = E(i1,i2,...,in) - sum of all subset deltaE values
      !! All subsets must have been computed already (guaranteed by processing fragments in order)
      integer(int64), intent(in) :: fragment_idx  !! Index of this fragment (already known)
      integer, intent(in) :: fragment(:), n
      type(fragment_lookup_t), intent(in) :: lookup  !! Pre-built hash table for lookups
      real(dp), intent(in) :: energies(:), delta_energies(:)  !! Pre-computed delta values
      real(dp) :: delta_E

      integer :: subset_size, i
      integer, allocatable :: indices(:), subset(:)
      integer(int64) :: subset_idx
      logical :: has_next

      ! Start with the full n-mer energy
      delta_E = energies(fragment_idx)

      ! Subtract all proper subsets (size 1 to n-1)
      do subset_size = 1, n - 1
         allocate (indices(subset_size))
         allocate (subset(subset_size))

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
            subset_idx = lookup%find(subset, subset_size)
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
                  error stop "Subset not found in bottom-up MBE!"
               end block
            end if

            ! Subtract pre-computed delta energy
            delta_E = delta_E - delta_energies(subset_idx)

            ! Get next combination
            call get_next_combination(indices, subset_size, n, has_next)
            if (.not. has_next) exit
         end do

         deallocate (indices, subset)
      end do

   end function compute_mbe

   subroutine map_fragment_to_system_gradient(frag_grad, monomers, sys_geom, sys_grad, bonds)
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
            error stop "Failed to build fragment in gradient mapping"
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

   subroutine compute_mbe_gradient(fragment_idx, fragment, lookup, results, delta_gradients, n, sys_geom, bonds)
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

      integer :: subset_size, i
      integer, allocatable :: indices(:), subset(:)
      integer(int64) :: subset_idx
      logical :: has_next

      ! Start with the full n-mer gradient mapped to system coordinates
      call map_fragment_to_system_gradient(results(fragment_idx)%gradient, fragment, &
                                           sys_geom, delta_gradients(:, :, fragment_idx), bonds)

      ! Subtract all proper subsets (size 1 to n-1)
      ! This is EXACTLY like the energy calculation, but for each gradient component
      do subset_size = 1, n - 1
         allocate (indices(subset_size))
         allocate (subset(subset_size))

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
            subset_idx = lookup%find(subset, subset_size)
            if (subset_idx < 0) error stop "Subset not found in MBE gradient!"

            ! Subtract pre-computed delta gradient (simple array subtraction in system coords)
            delta_gradients(:, :, fragment_idx) = delta_gradients(:, :, fragment_idx) - &
                                                  delta_gradients(:, :, subset_idx)

            ! Get next combination
            call get_next_combination(indices, subset_size, n, has_next)
            if (.not. has_next) exit
         end do

         deallocate (indices, subset)
      end do

   end subroutine compute_mbe_gradient

   subroutine compute_mbe_energy_gradient(polymers, fragment_count, max_level, results, sys_geom, &
                                          total_energy, total_gradient, bonds)
      !! Compute both MBE energy and gradient in a single pass
      !! This is more efficient than calling compute_mbe_energy and compute_mbe_gradient separately
      !! as it only builds the lookup table once and processes all fragments in one loop
      use mqc_result_types, only: calculation_result_t
      use mqc_config_parser, only: bond_t
      integer(int64), intent(in) :: fragment_count
      integer, intent(in) :: polymers(:, :), max_level
      type(calculation_result_t), intent(in) :: results(:)
      type(system_geometry_t), intent(in) :: sys_geom
      real(dp), intent(out) :: total_energy
      real(dp), intent(out) :: total_gradient(:, :)  !! (3, total_atoms)
      type(bond_t), intent(in), optional :: bonds(:)  !! Bond information for caps

      integer(int64) :: i
      integer :: fragment_size, nlevel, current_log_level, iatom
      real(dp), allocatable :: sum_by_level(:), delta_energies(:), energies(:)
      real(dp), allocatable :: delta_gradients(:, :, :)  !! (3, total_atoms, fragment_count)
      real(dp) :: delta_E
      logical :: do_detailed_print
      type(fragment_lookup_t) :: lookup
      type(timer_type) :: lookup_timer

      ! Validate that all fragments have gradients
      do i = 1_int64, fragment_count
         if (.not. results(i)%has_gradient) then
            call logger%error("Fragment "//to_char(i)//" does not have gradient!")
            error stop "Missing gradient in compute_mbe_energy_gradient"
         end if
      end do

      call logger%configuration(level=current_log_level)
      do_detailed_print = (current_log_level >= verbose_level)

      ! Allocate arrays for energy
      allocate (sum_by_level(max_level))
      allocate (delta_energies(fragment_count))
      allocate (energies(fragment_count))
      sum_by_level = 0.0_dp
      delta_energies = 0.0_dp

      ! Extract total energies from results
      do i = 1_int64, fragment_count
         energies(i) = results(i)%energy%total()
      end do

      ! Allocate arrays for gradient
      allocate (delta_gradients(3, sys_geom%total_atoms, fragment_count))
      delta_gradients = 0.0_dp
      total_gradient = 0.0_dp

      ! Build hash table for fast fragment lookups (shared for both energy and gradient)
      call lookup_timer%start()
      call lookup%init(fragment_count)
      do i = 1_int64, fragment_count
         fragment_size = count(polymers(i, :) > 0)
         call lookup%insert(polymers(i, :), fragment_size, i)
      end do
      call lookup_timer%stop()
      call logger%debug("Time to build lookup table: "//to_char(lookup_timer%get_elapsed_time())//" s")
      call logger%debug("Hash table size: "//to_char(lookup%table_size)// &
                        ", entries: "//to_char(lookup%n_entries))

      ! Bottom-up computation: process fragments by size (1-body, 2-body, 3-body, etc.)
      ! Process by n-mer level to ensure all subsets are computed before they're needed
      do nlevel = 1, max_level
         do i = 1_int64, fragment_count
            fragment_size = count(polymers(i, :) > 0)

            ! Only process fragments of the current nlevel
            if (fragment_size /= nlevel) cycle

            if (fragment_size == 1) then
               ! 1-body: deltaE = E, deltaG = G (no subsets to subtract)
               delta_energies(i) = energies(i)
               sum_by_level(1) = sum_by_level(1) + delta_energies(i)

               ! Map fragment gradient to system coordinates
               call map_fragment_to_system_gradient(results(i)%gradient, polymers(i, 1:fragment_size), &
                                                    sys_geom, delta_gradients(:, :, i), bonds)

            else if (fragment_size >= 2 .and. fragment_size <= max_level) then
               ! n-body: deltaE = E - sum(all subset deltaEs), deltaG = G - sum(all subset deltaGs)
               ! Energy delta
               delta_E = compute_mbe(i, polymers(i, 1:fragment_size), lookup, &
                                     energies, delta_energies, fragment_size)
               delta_energies(i) = delta_E
               sum_by_level(fragment_size) = sum_by_level(fragment_size) + delta_E

               ! Gradient delta
               call compute_mbe_gradient(i, polymers(i, 1:fragment_size), lookup, &
                                         results, delta_gradients, fragment_size, sys_geom, bonds)
            end if
         end do
      end do

      ! Clean up lookup table
      call lookup%destroy()

      ! Compute total energy
      total_energy = sum(sum_by_level)

      ! Compute total gradient
      do i = 1_int64, fragment_count
         fragment_size = count(polymers(i, :) > 0)
         if (fragment_size <= max_level) then
            total_gradient = total_gradient + delta_gradients(:, :, i)
         end if
      end do

      ! Print energy breakdown
      call logger%info("MBE Energy breakdown:")
      do nlevel = 1, max_level
         if (abs(sum_by_level(nlevel)) > 1e-15_dp) then
            block
               character(len=256) :: energy_line
               write (energy_line, '(a,i0,a,f20.10)') "  ", nlevel, "-body:  ", sum_by_level(nlevel)
               call logger%info(trim(energy_line))
            end block
         end if
      end do
      block
         character(len=256) :: total_line
         write (total_line, '(a,f20.10)') "  Total:   ", total_energy
         call logger%info(trim(total_line))
      end block

      ! Print gradient info
      call logger%info("MBE gradient computation completed")
      call logger%info("  Total gradient norm: "//to_char(sqrt(sum(total_gradient**2))))

      ! Print detailed gradient if verbose and small system
      if (current_log_level >= info_level .and. sys_geom%total_atoms < 100) then
         call logger%info(" ")
         call logger%info("Total MBE Gradient (Hartree/Bohr):")
         do iatom = 1, sys_geom%total_atoms
            block
               character(len=256) :: grad_line
               write (grad_line, '(a,i5,a,3f20.12)') "  Atom ", iatom, ": ", &
                  total_gradient(1, iatom), total_gradient(2, iatom), total_gradient(3, iatom)
               call logger%info(trim(grad_line))
            end block
         end do
         call logger%info(" ")
      end if

      ! Print detailed breakdown if requested
      if (do_detailed_print) then
         call print_detailed_breakdown(polymers, fragment_count, max_level, energies, delta_energies)
      end if

      ! Always write JSON file for machine-readable output (include gradient)
      call print_detailed_breakdown_json(polymers, fragment_count, max_level, energies, delta_energies, &
                                         sum_by_level, total_energy, total_gradient, results=results)

      deallocate (sum_by_level, delta_energies, energies, delta_gradients)

   end subroutine compute_mbe_energy_gradient

   subroutine compute_mbe_energy_gradient_hessian(polymers, fragment_count, max_level, results, sys_geom, &
                                                  total_energy, total_gradient, total_hessian, bonds)
      !! Compute MBE energy, gradient, and Hessian in a single pass
      !! Most efficient for simultaneous energy+gradient+Hessian calculations
      use mqc_result_types, only: calculation_result_t
      use mqc_config_parser, only: bond_t
      use mqc_physical_fragment, only: build_fragment_from_indices, redistribute_cap_gradients, &
                                       redistribute_cap_hessian
      integer(int64), intent(in) :: fragment_count
      integer, intent(in) :: polymers(:, :), max_level
      type(calculation_result_t), intent(in) :: results(:)
      type(system_geometry_t), intent(in) :: sys_geom
      real(dp), intent(out) :: total_energy
      real(dp), intent(out) :: total_gradient(:, :)  !! (3, total_atoms)
      real(dp), intent(out) :: total_hessian(:, :)   !! (3*total_atoms, 3*total_atoms)
      type(bond_t), intent(in), optional :: bonds(:)

      integer(int64) :: i
      integer :: fragment_size, nlevel, current_log_level, iatom
      real(dp), allocatable :: sum_by_level(:), delta_energies(:), energies(:)
      real(dp), allocatable :: delta_gradients(:, :, :), delta_hessians(:, :, :)
      type(fragment_lookup_t) :: lookup
      logical :: do_detailed_print
      integer :: hess_dim
      type(physical_fragment_t) :: fragment
      real(dp), allocatable :: temp_hess(:, :)
      real(dp) :: delta_E

      hess_dim = 3*sys_geom%total_atoms

      ! Validate all fragments have Hessians
      do i = 1, fragment_count
         if (.not. results(i)%has_hessian) then
            call logger%error("Fragment "//to_char(i)//" does not have Hessian!")
            error stop "Missing Hessian in compute_mbe_energy_gradient_hessian"
         end if
      end do

      call logger%configuration(level=current_log_level)
      do_detailed_print = (current_log_level >= verbose_level)

      allocate (sum_by_level(max_level))
      allocate (delta_energies(fragment_count))
      allocate (energies(fragment_count))
      allocate (delta_gradients(3, sys_geom%total_atoms, fragment_count))
      allocate (delta_hessians(hess_dim, hess_dim, fragment_count))
      allocate (temp_hess(hess_dim, hess_dim))

      sum_by_level = 0.0_dp
      delta_energies = 0.0_dp
      energies = 0.0_dp
      delta_gradients = 0.0_dp
      delta_hessians = 0.0_dp
      total_gradient = 0.0_dp
      total_hessian = 0.0_dp

      ! Build lookup table for fragment indices
      call lookup%init(fragment_count)
      do i = 1_int64, fragment_count
         fragment_size = count(polymers(i, :) > 0)
         if (fragment_size > 0 .and. fragment_size <= max_level) then
            call lookup%insert(polymers(i, 1:fragment_size), fragment_size, i)
         end if
      end do

      ! Extract energies first
      do i = 1_int64, fragment_count
         energies(i) = results(i)%energy%total()
      end do

      ! Compute delta energies, gradients, and Hessians for each fragment
      ! Process by n-mer level to ensure all subsets are computed before they're needed
      do nlevel = 1, max_level
         do i = 1_int64, fragment_count
            fragment_size = count(polymers(i, :) > 0)

            ! Only process fragments of the current nlevel
            if (fragment_size /= nlevel) cycle

            if (fragment_size == 1) then
               ! 1-body: delta = value (no subsets)
               delta_energies(i) = energies(i)
               sum_by_level(1) = sum_by_level(1) + delta_energies(i)

               ! Map fragment gradient and Hessian to system coordinates
               call map_fragment_to_system_gradient(results(i)%gradient, polymers(i, 1:fragment_size), &
                                                    sys_geom, delta_gradients(:, :, i), bonds)
               call map_fragment_to_system_hessian(results(i)%hessian, polymers(i, 1:fragment_size), &
                                                   sys_geom, delta_hessians(:, :, i), bonds)

            else if (fragment_size >= 2 .and. fragment_size <= max_level) then
               ! n-body: delta = value - sum(all subset deltas)
               delta_E = compute_mbe(i, polymers(i, 1:fragment_size), lookup, &
                                     energies, delta_energies, fragment_size)
               delta_energies(i) = delta_E
               sum_by_level(fragment_size) = sum_by_level(fragment_size) + delta_E

               ! Gradient delta
               call compute_mbe_gradient(i, polymers(i, 1:fragment_size), lookup, &
                                         results, delta_gradients, fragment_size, sys_geom, bonds)

               ! Hessian delta
               call compute_mbe_hessian(i, polymers(i, 1:fragment_size), lookup, &
                                        results, delta_hessians, fragment_size, sys_geom, bonds)
            end if
         end do
      end do

      ! Clean up lookup table
      call lookup%destroy()

      ! Compute total energy, gradient, and Hessian
      total_energy = sum(sum_by_level)
      do i = 1_int64, fragment_count
         fragment_size = count(polymers(i, :) > 0)
         if (fragment_size <= max_level) then
            total_gradient = total_gradient + delta_gradients(:, :, i)
            total_hessian = total_hessian + delta_hessians(:, :, i)
         end if
      end do

      ! Print energy breakdown
      call logger%info("MBE Energy breakdown:")
      do nlevel = 1, max_level
         if (abs(sum_by_level(nlevel)) > 1e-15_dp) then
            block
               character(len=256) :: energy_line
               write (energy_line, '(a,i0,a,f20.10)') "  ", nlevel, "-body:  ", sum_by_level(nlevel)
               call logger%info(trim(energy_line))
            end block
         end if
      end do
      block
         character(len=256) :: total_line
         write (total_line, '(a,f20.10)') "  Total:   ", total_energy
         call logger%info(trim(total_line))
      end block

      ! Print gradient and Hessian info
      call logger%info("MBE gradient computation completed")
      call logger%info("  Total gradient norm: "//to_char(sqrt(sum(total_gradient**2))))
      call logger%info("MBE Hessian computation completed")
      call logger%info("  Total Hessian Frobenius norm: "//to_char(sqrt(sum(total_hessian**2))))

      ! Print detailed gradient if verbose and small system
      if (current_log_level >= info_level .and. sys_geom%total_atoms < 100) then
         call logger%info(" ")
         call logger%info("Total MBE Gradient (Hartree/Bohr):")
         do iatom = 1, sys_geom%total_atoms
            block
               character(len=256) :: grad_line
               write (grad_line, '(a,i5,a,3f20.12)') "  Atom ", iatom, ": ", &
                  total_gradient(1, iatom), total_gradient(2, iatom), total_gradient(3, iatom)
               call logger%info(trim(grad_line))
            end block
         end do
         call logger%info(" ")
      end if

      ! Always write JSON file for machine-readable output (include gradient and Hessian)
      call print_detailed_breakdown_json(polymers, fragment_count, max_level, energies, delta_energies, &
                                         sum_by_level, total_energy, total_gradient, total_hessian, results)

      deallocate (sum_by_level, delta_energies, energies, delta_gradients, delta_hessians, temp_hess)

   end subroutine compute_mbe_energy_gradient_hessian

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
      integer, allocatable :: indices(:), subset(:)
      integer(int64) :: subset_idx
      logical :: has_next

      hess_dim = 3*sys_geom%total_atoms

      ! Start with the full n-mer Hessian mapped to system coordinates
      call map_fragment_to_system_hessian(results(fragment_idx)%hessian, fragment, &
                                          sys_geom, delta_hessians(:, :, fragment_idx), bonds)

      ! Subtract all proper subsets (size 1 to n-1)
      do subset_size = 1, n - 1
         allocate (indices(subset_size))
         indices = [(i, i=1, subset_size)]
         allocate (subset(subset_size))

         has_next = .true.
         do while (has_next)
            subset = fragment(indices)
            subset_idx = lookup%find(subset, subset_size)

            if (subset_idx > 0) then
               ! Subtract this subset's delta Hessian
               delta_hessians(:, :, fragment_idx) = delta_hessians(:, :, fragment_idx) - &
                                                    delta_hessians(:, :, subset_idx)
            end if

            call get_next_combination(indices, subset_size, n, has_next)
         end do

         deallocate (indices, subset)
      end do

   end subroutine compute_mbe_hessian

   subroutine compute_gmbe_energy(monomers, n_monomers, monomer_results, &
                                  n_intersections, intersection_results, &
                                  intersection_sets, intersection_levels, total_energy)
      !! Compute GMBE (Generalized Many-Body Expansion) energy with full inclusion-exclusion
      !!
      !! For overlapping fragments, the GMBE formula follows the inclusion-exclusion principle:
      !!   E_total = sum(E_monomers)
      !!             - sum(E_2-way_intersections)
      !!             + sum(E_3-way_intersections)
      !!             - sum(E_4-way_intersections)
      !!             + ...
      !!
      !! Example for three overlapping fragments (gly0, gly1, gly2):
      !!   E = E(gly0) + E(gly1) + E(gly2)
      !!       - E(gly0∩gly1) - E(gly0∩gly2) - E(gly1∩gly2)
      !!       + E(gly0∩gly1∩gly2)
      !!
      !! This correctly accounts for all overlapping regions following inclusion-exclusion.
      use mqc_result_types, only: calculation_result_t

      integer, intent(in) :: monomers(:)              !! Monomer indices (1-based)
      integer, intent(in) :: n_monomers               !! Number of monomers
      type(calculation_result_t), intent(in) :: monomer_results(:)     !! Monomer energies
      integer, intent(in) :: n_intersections          !! Number of intersection fragments
      type(calculation_result_t), intent(in), optional :: intersection_results(:)  !! Intersection energies
      integer, intent(in), optional :: intersection_sets(:, :)  !! k-tuples that created each intersection (n_monomers, n_intersections)
      integer, intent(in), optional :: intersection_levels(:)   !! Level k of each intersection
      real(dp), intent(out) :: total_energy           !! Total GMBE energy

      integer :: i, k, max_level
      real(dp) :: monomer_energy
      real(dp), allocatable :: level_energies(:)
      integer, allocatable :: level_counts(:)
      real(dp) :: sign_factor

      ! Sum monomer energies
      monomer_energy = 0.0_dp
      do i = 1, n_monomers
         monomer_energy = monomer_energy + monomer_results(i)%energy%total()
      end do

      ! Start with monomer contribution
      total_energy = monomer_energy

      if (n_intersections > 0 .and. present(intersection_results) .and. &
          present(intersection_sets) .and. present(intersection_levels)) then
         ! Find maximum intersection level
         max_level = maxval(intersection_levels)

         ! Allocate arrays to track contributions by level
         allocate (level_energies(2:max_level))
         allocate (level_counts(2:max_level))
         level_energies = 0.0_dp
         level_counts = 0

         ! Sum intersection energies by level with alternating signs
         do i = 1, n_intersections
            k = intersection_levels(i)
            level_energies(k) = level_energies(k) + intersection_results(i)%energy%total()
            level_counts(k) = level_counts(k) + 1
         end do

         ! Apply inclusion-exclusion: sign is (-1)^(k+1) for k-way intersections
         do k = 2, max_level
            if (level_counts(k) > 0) then
               sign_factor = real((-1)**(k + 1), dp)
               total_energy = total_energy + sign_factor*level_energies(k)
            end if
         end do

         ! Print breakdown
         call logger%info("GMBE Energy breakdown (Inclusion-Exclusion Principle):")
         block
            character(len=256) :: line
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
         end block

         ! Print detailed intersection info at debug level
         if (n_intersections > 0) then
            call logger%debug("GMBE intersection details:")
            do i = 1, n_intersections
               block
                  character(len=512) :: detail_line
                  character(len=256) :: set_str
                  integer :: j, set_size

                  ! Build string of fragment indices in this intersection
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
               end block
            end do
         end if

         deallocate (level_energies, level_counts)
      else
         ! No intersections - just report monomer sum
         call logger%info("GMBE Energy breakdown:")
         block
            character(len=256) :: line
            write (line, '(a,i0,a,f20.10)') "  Monomers (", n_monomers, "): ", monomer_energy
            call logger%info(trim(line))
            write (line, '(a,f20.10)') "  Total GMBE:  ", total_energy
            call logger%info(trim(line))
         end block
      end if

   end subroutine compute_gmbe_energy

   subroutine compute_gmbe_energy_gradient(monomers, n_monomers, monomer_results, &
                                           n_intersections, intersection_results, &
                                           intersection_sets, intersection_levels, &
                                           sys_geom, total_energy, total_gradient, bonds)
      !! Compute GMBE (Generalized Many-Body Expansion) energy and gradient with full inclusion-exclusion
      !!
      !! Extends compute_gmbe_energy to also compute gradients for overlapping fragments.
      !! Gradient redistribution handles hydrogen caps at broken bonds.
      !!
      !! For overlapping fragments, the GMBE formula follows the inclusion-exclusion principle:
      !!   E_total = sum(E_monomers) - sum(E_2-way) + sum(E_3-way) - ...
      !!   ∇E_total = sum(∇E_monomers) - sum(∇E_2-way) + sum(∇E_3-way) - ...
      use mqc_result_types, only: calculation_result_t
      use mqc_physical_fragment, only: build_fragment_from_indices, build_fragment_from_atom_list, &
                                       redistribute_cap_gradients
      use mqc_config_parser, only: bond_t
      use mqc_error, only: error_t
      use pic_logger, only: info_level

      integer, intent(in) :: monomers(:)              !! Monomer indices (1-based)
      integer, intent(in) :: n_monomers               !! Number of monomers
      type(calculation_result_t), intent(in) :: monomer_results(:)     !! Monomer energies and gradients
      integer, intent(in) :: n_intersections          !! Number of intersection fragments
      type(calculation_result_t), intent(in), optional :: intersection_results(:)  !! Intersection energies and gradients
      integer, intent(in), optional :: intersection_sets(:, :)  !! k-tuples that created each intersection (n_monomers, n_intersections)
      integer, intent(in), optional :: intersection_levels(:)   !! Level k of each intersection
      type(system_geometry_t), intent(in) :: sys_geom
      real(dp), intent(out) :: total_energy           !! Total GMBE energy
      real(dp), intent(out) :: total_gradient(:, :)   !! Total GMBE gradient (3, total_atoms)
      type(bond_t), intent(in), optional :: bonds(:)  !! Bond information for caps

      integer :: i, k, max_level, current_log_level
      real(dp) :: monomer_energy
      real(dp), allocatable :: level_energies(:)
      integer, allocatable :: level_counts(:)
      real(dp) :: sign_factor
      type(physical_fragment_t) :: fragment
      type(error_t) :: error
      integer, allocatable :: monomer_idx(:)

      ! Zero out total gradient
      total_gradient = 0.0_dp

      ! Sum monomer energies and gradients
      monomer_energy = 0.0_dp
      do i = 1, n_monomers
         monomer_energy = monomer_energy + monomer_results(i)%energy%total()

         ! Accumulate monomer gradients
         if (monomer_results(i)%has_gradient) then
            ! Rebuild fragment to get local→global mapping
            allocate (monomer_idx(1))
            monomer_idx(1) = monomers(i)
            call build_fragment_from_indices(sys_geom, monomer_idx, fragment, error, bonds)
            if (error%has_error()) then
               call logger%error(error%get_full_trace())
               error stop "Failed to build monomer fragment in GMBE gradient"
            end if
            call redistribute_cap_gradients(fragment, monomer_results(i)%gradient, total_gradient)
            call fragment%destroy()
            deallocate (monomer_idx)
         end if
      end do

      ! Start with monomer contribution
      total_energy = monomer_energy

      if (n_intersections > 0 .and. present(intersection_results) .and. &
          present(intersection_sets) .and. present(intersection_levels)) then
         ! Find maximum intersection level
         max_level = maxval(intersection_levels)

         ! Allocate arrays to track contributions by level
         allocate (level_energies(2:max_level))
         allocate (level_counts(2:max_level))
         level_energies = 0.0_dp
         level_counts = 0

         ! Sum intersection energies by level with alternating signs
         do i = 1, n_intersections
            k = intersection_levels(i)
            sign_factor = real((-1)**(k + 1), dp)
            level_energies(k) = level_energies(k) + intersection_results(i)%energy%total()
            level_counts(k) = level_counts(k) + 1

            ! Accumulate intersection gradients with sign
            if (intersection_results(i)%has_gradient) then
               ! Intersection gradients are accumulated with the same sign as energy
               ! Note: We don't need to rebuild intersection fragments here because
               ! the gradient redistribution was already done during fragment calculation
               ! via build_fragment_from_atom_list
               ! TODO: This needs to be implemented properly - we need the fragment geometry
               ! to get the local→global mapping. For now, log a warning.
               call logger%warning("GMBE gradient with intersections not fully implemented yet!")
               call logger%warning("Intersection gradient redistribution requires storing fragment geometry")
            end if
         end do

         ! Apply inclusion-exclusion to energy: sign is (-1)^(k+1) for k-way intersections
         do k = 2, max_level
            if (level_counts(k) > 0) then
               sign_factor = real((-1)**(k + 1), dp)
               total_energy = total_energy + sign_factor*level_energies(k)
            end if
         end do

         ! Print breakdown
         call logger%info("GMBE Energy breakdown (Inclusion-Exclusion Principle):")
         block
            character(len=256) :: line
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
         end block

         deallocate (level_energies, level_counts)
      else
         ! No intersections - just report monomer sum
         call logger%info("GMBE Energy breakdown:")
         block
            character(len=256) :: line
            write (line, '(a,i0,a,f20.10)') "  Monomers (", n_monomers, "): ", monomer_energy
            call logger%info(trim(line))
            write (line, '(a,f20.10)') "  Total GMBE:  ", total_energy
            call logger%info(trim(line))
         end block
      end if

      ! Print gradient info (same format as MBE)
      call logger%info("GMBE gradient computation completed")
      call logger%info("  Total gradient norm: "//to_char(sqrt(sum(total_gradient**2))))

      ! Print detailed gradient if info level and small system
      call logger%configuration(level=current_log_level)
      if (current_log_level >= info_level .and. sys_geom%total_atoms < 100) then
         call logger%info(" ")
         call logger%info("Total GMBE Gradient (Hartree/Bohr):")
         do i = 1, sys_geom%total_atoms
            block
               character(len=256) :: grad_line
               write (grad_line, '(a,i5,a,3f20.12)') "  Atom ", i, ": ", &
                  total_gradient(1, i), total_gradient(2, i), total_gradient(3, i)
               call logger%info(trim(grad_line))
            end block
         end do
         call logger%info(" ")
      end if

   end subroutine compute_gmbe_energy_gradient

   subroutine compute_gmbe_energy_gradient_hessian(monomers, n_monomers, monomer_results, &
                                                   n_intersections, intersection_results, &
                                                   intersection_sets, intersection_levels, &
                                                   sys_geom, total_energy, total_gradient, total_hessian, bonds)
      !! Compute GMBE energy, gradient, and Hessian with full inclusion-exclusion
      !! TODO: Full implementation with intersection Hessians pending
      use mqc_result_types, only: calculation_result_t
      use mqc_physical_fragment, only: build_fragment_from_indices, redistribute_cap_gradients, &
                                       redistribute_cap_hessian
      use mqc_config_parser, only: bond_t
      use mqc_error, only: error_t
      use pic_logger, only: info_level

      integer, intent(in) :: monomers(:)
      integer, intent(in) :: n_monomers
      type(calculation_result_t), intent(in) :: monomer_results(:)
      integer, intent(in) :: n_intersections
      type(calculation_result_t), intent(in), optional :: intersection_results(:)
      integer, intent(in), optional :: intersection_sets(:, :)
      integer, intent(in), optional :: intersection_levels(:)
      type(system_geometry_t), intent(in) :: sys_geom
      real(dp), intent(out) :: total_energy
      real(dp), intent(out) :: total_gradient(:, :)
      real(dp), intent(out) :: total_hessian(:, :)
      type(bond_t), intent(in), optional :: bonds(:)

      integer :: i, k, max_level, current_log_level, hess_dim
      real(dp) :: monomer_energy
      real(dp), allocatable :: level_energies(:)
      integer, allocatable :: level_counts(:)
      real(dp) :: sign_factor
      type(physical_fragment_t) :: fragment
      type(error_t) :: error
      integer, allocatable :: monomer_idx(:)

      ! Zero out outputs
      total_gradient = 0.0_dp
      total_hessian = 0.0_dp
      hess_dim = 3*sys_geom%total_atoms

      ! Sum monomer energies, gradients, and Hessians
      monomer_energy = 0.0_dp
      do i = 1, n_monomers
         monomer_energy = monomer_energy + monomer_results(i)%energy%total()

         ! Accumulate monomer gradients
         if (monomer_results(i)%has_gradient) then
            allocate (monomer_idx(1))
            monomer_idx(1) = monomers(i)
            call build_fragment_from_indices(sys_geom, monomer_idx, fragment, error, bonds)
            if (error%has_error()) then
               call logger%error(error%get_full_trace())
               error stop "Failed to build monomer fragment in GMBE gradient+Hessian"
            end if
            call redistribute_cap_gradients(fragment, monomer_results(i)%gradient, total_gradient)
            if (monomer_results(i)%has_hessian) then
               call redistribute_cap_hessian(fragment, monomer_results(i)%hessian, total_hessian)
            end if
            call fragment%destroy()
            deallocate (monomer_idx)
         end if
      end do

      total_energy = monomer_energy

      ! Intersection Hessians not yet implemented
      if (n_intersections > 0 .and. present(intersection_results)) then
         call logger%warning("GMBE Hessian with intersections not yet fully implemented!")
         call logger%warning("Only monomer Hessians will be included")
      end if

      ! Print info
      call logger%info("GMBE Energy breakdown:")
      call logger%info("  Monomers ("//to_char(n_monomers)//"): "//to_char(monomer_energy))
      call logger%info("  Total GMBE:  "//to_char(total_energy))
      call logger%info("GMBE gradient computation completed")
      call logger%info("  Total gradient norm: "//to_char(sqrt(sum(total_gradient**2))))
      call logger%info("GMBE Hessian computation completed")
      call logger%info("  Total Hessian Frobenius norm: "//to_char(sqrt(sum(total_hessian**2))))

   end subroutine compute_gmbe_energy_gradient_hessian

end module mqc_mbe
