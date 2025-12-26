!! Many-Body Expansion (MBE) calculation module
module mqc_mbe
   !! Implements hierarchical many-body expansion for fragment-based quantum chemistry
   !! calculations with MPI parallelization and energy/gradient computation.
   use pic_types, only: int32, int64, dp
   use pic_timer, only: timer_type
   use pic_blas_interfaces, only: pic_gemm, pic_dot
   use pic_mpi_lib, only: comm_t, send, recv, iprobe, MPI_Status, MPI_ANY_SOURCE, MPI_ANY_TAG
   use pic_logger, only: logger => global_logger
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
   public :: compute_mbe_energy

contains

   subroutine compute_mbe_energy(polymers, fragment_count, max_level, results, total_energy)
      !! Compute the many-body expansion (MBE) energy
      !! Total = sum(E(i)) + sum(deltaE(ij)) + sum(deltaE(ijk)) + ...
      !! General n-body correction:
      !! deltaE(i1,i2,...,in) = E(i1,i2,...,in) - sum of all lower-order terms
      !! Uses int64 for fragment_count to handle large fragment counts that overflow int32.
      !! Detailed breakdown is printed only if logger level is verbose or higher.
      use pic_logger, only: verbose_level
      use mqc_result_types, only: calculation_result_t
      integer(int64), intent(in) :: fragment_count
      integer, intent(in) :: polymers(:, :), max_level
      type(calculation_result_t), intent(in) :: results(:)
      real(dp), intent(out) :: total_energy

      integer(int64) :: i
      integer :: fragment_size, body_level, current_log_level
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
      ! This eliminates recursion and redundant subset computations
      do i = 1_int64, fragment_count
         fragment_size = count(polymers(i, :) > 0)

         if (fragment_size == 1) then
            ! 1-body: deltaE = E (no subsets to subtract)
            delta_energies(i) = energies(i)
            sum_by_level(1) = sum_by_level(1) + delta_energies(i)
         else if (fragment_size >= 2 .and. fragment_size <= max_level) then
            ! n-body: deltaE = E - sum(all subset deltaEs)
            ! All subsets have already been computed in previous iterations
            delta_E = compute_mbe(i, polymers(i, 1:fragment_size), lookup, &
                                  energies, delta_energies, fragment_size)
            delta_energies(i) = delta_E
            sum_by_level(fragment_size) = sum_by_level(fragment_size) + delta_E
         end if
      end do

      ! Clean up lookup table
      call lookup%destroy()

      total_energy = sum(sum_by_level)

      ! Print text summary to console
      call logger%info("MBE Energy breakdown:")
      do body_level = 1, max_level
         if (abs(sum_by_level(body_level)) > 1e-15_dp) then
            block
               character(len=256) :: energy_line
               write (energy_line, '(a,i0,a,f20.10)') "  ", body_level, "-body:  ", sum_by_level(body_level)
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
                                         sum_by_level, total_energy)

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
            if (subset_idx < 0) error stop "Subset not found in bottom-up MBE!"

            ! Subtract pre-computed delta energy
            delta_E = delta_E - delta_energies(subset_idx)

            ! Get next combination
            call get_next_combination(indices, subset_size, n, has_next)
            if (.not. has_next) exit
         end do

         deallocate (indices, subset)
      end do

   end function compute_mbe

end module mqc_mbe
