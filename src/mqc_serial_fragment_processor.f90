submodule(mqc_mbe_fragment_distribution_scheme) mqc_serial_fragment_processor
   implicit none

contains

   module subroutine serial_fragment_processor(total_fragments, polymers, max_level, sys_geom, method, calc_type, bonds)
      !! Process all fragments serially in single-rank mode
      !! This is used when running with only 1 MPI rank
      use mqc_error, only: error_t
      integer(int64), intent(in) :: total_fragments
      integer, intent(in) :: polymers(:, :), max_level
      type(system_geometry_t), intent(in) :: sys_geom
      integer(int32), intent(in) :: method
      integer(int32), intent(in) :: calc_type
      type(bond_t), intent(in), optional :: bonds(:)

      integer(int64) :: frag_idx
      integer :: fragment_size, current_log_level, iatom
      integer, allocatable :: fragment_indices(:)
      type(calculation_result_t), allocatable :: results(:)
      real(dp) :: mbe_total_energy
      real(dp), allocatable :: mbe_total_gradient(:, :)
      real(dp), allocatable :: mbe_total_hessian(:, :)
      type(physical_fragment_t) :: phys_frag
      type(timer_type) :: coord_timer
      integer(int32) :: calc_type_local
      type(error_t) :: error

      calc_type_local = calc_type

      call logger%info("Processing "//to_char(total_fragments)//" fragments serially...")
      call logger%info("  Calculation type: "//calc_type_to_string(calc_type_local))

      allocate (results(total_fragments))

      call omp_set_num_threads(1)
      call coord_timer%start()
      do frag_idx = 1_int64, total_fragments
         fragment_size = count(polymers(frag_idx, :) > 0)
         allocate (fragment_indices(fragment_size))
         fragment_indices = polymers(frag_idx, 1:fragment_size)

         call build_fragment_from_indices(sys_geom, fragment_indices, phys_frag, error, bonds)
         if (error%has_error()) then
            call logger%error(error%get_full_trace())
            error stop "Failed to build fragment in serial processing"
         end if

         call do_fragment_work(frag_idx, results(frag_idx), method, phys_frag, calc_type=calc_type_local)

         ! Debug output for gradients
         if (calc_type_local == CALC_TYPE_GRADIENT .and. results(frag_idx)%has_gradient) then
            call logger%configuration(level=current_log_level)
            if (current_log_level >= verbose_level) then
               block
                  character(len=512) :: debug_line
                  integer :: iatom_local
                  write (debug_line, '(a,i0,a,*(i0,1x))') "Fragment ", frag_idx, " monomers: ", fragment_indices
                  call logger%verbose(trim(debug_line))
                  write (debug_line, '(a,f25.15)') "  Energy: ", results(frag_idx)%energy%total()
                  call logger%verbose(trim(debug_line))
                  write (debug_line, '(a,f25.15)') "  Gradient norm: ", sqrt(sum(results(frag_idx)%gradient**2))
                  call logger%verbose(trim(debug_line))
                  if (size(results(frag_idx)%gradient, 2) <= 20) then
                     call logger%verbose("  Fragment gradient:")
                     do iatom_local = 1, size(results(frag_idx)%gradient, 2)
                        write (debug_line, '(a,i3,a,3f20.12)') "    Atom ", iatom_local, ": ", &
                           results(frag_idx)%gradient(1, iatom_local), &
                           results(frag_idx)%gradient(2, iatom_local), &
                           results(frag_idx)%gradient(3, iatom_local)
                        call logger%verbose(trim(debug_line))
                     end do
                  end if
               end block
            end if
         end if

         call phys_frag%destroy()
         deallocate (fragment_indices)

         if (mod(frag_idx, max(1_int64, total_fragments/10)) == 0 .or. frag_idx == total_fragments) then
            call logger%info("  Processed "//to_char(frag_idx)//"/"//to_char(total_fragments)// &
                             " fragments ["//to_char(coord_timer%get_elapsed_time())//" s]")
         end if
      end do
      call coord_timer%stop()
      call logger%info("Time to evaluate all fragments "//to_char(coord_timer%get_elapsed_time())//" s")
      call omp_set_num_threads(omp_get_max_threads())

      call logger%info("All fragments processed")

      call logger%info(" ")
      call logger%info("Computing Many-Body Expansion (MBE)...")
      call coord_timer%start()

      ! Use combined function if computing gradients or Hessians (more efficient)
      if (calc_type_local == CALC_TYPE_HESSIAN) then
         allocate (mbe_total_gradient(3, sys_geom%total_atoms))
         allocate (mbe_total_hessian(3*sys_geom%total_atoms, 3*sys_geom%total_atoms))
         call compute_mbe_energy_gradient_hessian(polymers, total_fragments, max_level, results, sys_geom, &
                                                  mbe_total_energy, mbe_total_gradient, mbe_total_hessian, bonds)
         deallocate (mbe_total_gradient, mbe_total_hessian)
      else if (calc_type_local == CALC_TYPE_GRADIENT) then
         allocate (mbe_total_gradient(3, sys_geom%total_atoms))
         call compute_mbe_energy_gradient(polymers, total_fragments, max_level, results, sys_geom, &
                                          mbe_total_energy, mbe_total_gradient, bonds)
         deallocate (mbe_total_gradient)
      else
         call compute_mbe_energy(polymers, total_fragments, max_level, results, mbe_total_energy)
      end if

      call coord_timer%stop()
      call logger%info("Time to compute MBE "//to_char(coord_timer%get_elapsed_time())//" s")

      deallocate (results)

   end subroutine serial_fragment_processor

end submodule mqc_serial_fragment_processor
