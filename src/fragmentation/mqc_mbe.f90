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
                           TAG_WORKER_SCALAR_RESULT, TAG_WORKER_MATRIX_RESULT, &
                           TAG_NODE_REQUEST, TAG_NODE_FRAGMENT, TAG_NODE_FINISH, &
                           TAG_NODE_SCALAR_RESULT, TAG_NODE_MATRIX_RESULT
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, build_fragment_from_indices, to_angstrom
   use mqc_frag_utils, only: next_combination, find_fragment_index

   implicit none
   private

   ! Public interface
   public :: compute_mbe_energy

contains

   subroutine compute_mbe_energy(polymers, fragment_count, max_level, energies, total_energy)
      !! Compute the many-body expansion (MBE) energy
      !! Total = sum(E(i)) + sum(deltaE(ij)) + sum(deltaE(ijk)) + ...
      !! General n-body correction:
      !! deltaE(i1,i2,...,in) = E(i1,i2,...,in) - sum of all lower-order terms
      !! Uses int64 for fragment_count to handle large fragment counts that overflow int32.
      !! Detailed breakdown is printed only if logger level is verbose or higher.
      use pic_logger, only: verbose_level
      integer(int64), intent(in) :: fragment_count
      integer, intent(in) :: polymers(:, :), max_level
      real(dp), intent(in) :: energies(:)
      real(dp), intent(out) :: total_energy

      integer(int64) :: i
      integer :: fragment_size, body_level, current_log_level
      real(dp), allocatable :: sum_by_level(:), delta_energies(:)
      real(dp) :: delta_E
      logical :: do_detailed_print

      ! Query logger to decide if we should print detailed breakdown
      call logger%configuration(level=current_log_level)
      do_detailed_print = (current_log_level >= verbose_level)

      allocate (sum_by_level(max_level))
      allocate (delta_energies(fragment_count))
      sum_by_level = 0.0_dp
      delta_energies = 0.0_dp

      ! Sum over all fragments by their size
      do i = 1_int64, fragment_count
         fragment_size = count(polymers(i, :) > 0)

         if (fragment_size == 1) then
            ! 1-body terms: just the monomer energies
            delta_E = energies(i)
            delta_energies(i) = delta_E
            sum_by_level(1) = sum_by_level(1) + delta_E
         else if (fragment_size >= 2 .and. fragment_size <= max_level) then
            ! n-body corrections for n >= 2
            delta_E = compute_delta_nbody(polymers(i, 1:fragment_size), polymers, energies, &
                                          fragment_count, fragment_size)
            delta_energies(i) = delta_E
            sum_by_level(fragment_size) = sum_by_level(fragment_size) + delta_E
         end if
      end do

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

      deallocate (sum_by_level, delta_energies)

   end subroutine compute_mbe_energy

   recursive function compute_delta_nbody(fragment, polymers, energies, fragment_count, n) result(delta_E)
      !! Compute general n-body correction using recursion
      !! deltaE(i1,i2,...,in) = E(i1,i2,...,in)
      !!                        - sum over all proper subsets of deltaE or E
      !! For n=2: deltaE(ij) = E(ij) - E(i) - E(j)
      !! For n=3: deltaE(ijk) = E(ijk) - deltaE(ij) - deltaE(ik) - deltaE(jk) - E(i) - E(j) - E(k)
      !! For n>=4: same pattern recursively
      !! Uses int64 for fragment_count to handle large fragment counts that overflow int32.
      integer, intent(in) :: fragment(:), polymers(:, :), n
      integer(int64), intent(in) :: fragment_count
      real(dp), intent(in) :: energies(:)
      real(dp) :: delta_E

      integer(int64) :: idx_n
      integer :: subset_size, i, j, num_subsets
      integer, allocatable :: subset(:), indices(:)
      real(dp) :: E_n, subset_contribution

      ! Find the energy of this n-mer
      idx_n = find_fragment_index(fragment, polymers, fragment_count, n)
      E_n = energies(idx_n)

      ! Start with the full n-mer energy
      delta_E = E_n

      ! Subtract contributions from all proper subsets (size 1 to n-1)
      do subset_size = 1, n - 1
         ! Generate all subsets of this size and subtract their contributions
         call generate_and_subtract_subsets(fragment, subset_size, n, polymers, energies, &
                                            fragment_count, delta_E)
      end do

   end function compute_delta_nbody

   subroutine generate_and_subtract_subsets(fragment, subset_size, n, polymers, energies, &
                                            fragment_count, delta_E)
      !! Generate all subsets of given size and subtract their contribution
      !! Uses int64 for fragment_count to handle large fragment counts that overflow int32.
      integer, intent(in) :: fragment(:), subset_size, n, polymers(:, :)
      integer(int64), intent(in) :: fragment_count
      real(dp), intent(in) :: energies(:)
      real(dp), intent(inout) :: delta_E

      integer, allocatable :: indices(:), subset(:)
      integer :: i

      allocate (indices(subset_size))
      allocate (subset(subset_size))

      ! Initialize indices for first combination
      do i = 1, subset_size
         indices(i) = i
      end do

      ! Loop through all combinations
      do
         ! Build the current subset
         do i = 1, subset_size
            subset(i) = fragment(indices(i))
         end do

         ! Subtract this subset's contribution
         if (subset_size == 1) then
            ! 1-body: just subtract the monomer energy
            delta_E = delta_E - energies(find_fragment_index(subset, polymers, fragment_count, 1))
         else
            ! n-body: recursively compute and subtract deltaE for this subset
            delta_E = delta_E - compute_delta_nbody(subset, polymers, energies, fragment_count, subset_size)
         end if

         ! Get next combination
         if (.not. next_combination(indices, subset_size, n)) exit
      end do

      deallocate (indices, subset)

   end subroutine generate_and_subtract_subsets

end module mqc_mbe
