!! The record an optimization leaves behind
module mqc_optimization_output
   !! Writes `output_<name>_optimization.json` and
   !! `output_<name>_trajectory.xyz`, both separate from the calculation's own
   !! `output_<name>.json`.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use json_module, only: json_core, json_value
   use mqc_physical_fragment, only: to_angstrom
   use mqc_elements, only: element_number_to_symbol
   use mqc_error, only: error_t, ERROR_IO
   implicit none
   private

   public :: optimization_record_t
   public :: write_optimization_json
   public :: write_trajectory_xyz

   type :: optimization_record_t
      !! Everything an optimization is worth reporting afterwards
      logical :: converged = .false.
      integer :: n_steps = 0            !! Geometries the optimizer accepted
      integer :: n_evaluations = 0      !! Energy and gradient calculations run
      real(dp) :: final_energy = 0.0_dp
      real(dp) :: final_gradient_max = 0.0_dp
      real(dp) :: gradient_tolerance = 0.0_dp
      character(len=:), allocatable :: coordinates  !! cartesian, hdlc, dlc
      character(len=:), allocatable :: algorithm

      integer, allocatable :: element_numbers(:)     !! (n_atoms)
      real(dp), allocatable :: final_coordinates(:, :)  !! (3, n_atoms), Bohr

      real(dp), allocatable :: trajectory(:, :, :)
         !! (3, n_atoms, n_steps) in Bohr. Unallocated when the deck turned the
         !! trajectory off.
      real(dp), allocatable :: trajectory_energies(:)
         !! (n_steps), Hartree, one per trajectory frame
   end type optimization_record_t

contains

   subroutine write_optimization_json(path, record, error)
      !! Write the optimization record as a JSON document
      character(len=*), intent(in) :: path
      type(optimization_record_t), intent(in) :: record
      type(error_t), intent(inout) :: error

      type(json_core) :: json
      type(json_value), pointer :: root, main_obj, geom_obj, traj_arr, frame_obj, symbols_arr
      integer :: iatom, istep, n_atoms
      logical :: status_ok

      n_atoms = 0
      if (allocated(record%final_coordinates)) n_atoms = size(record%final_coordinates, 2)

      call json%initialize()
      call json%create_object(root, "")
      call json%create_object(main_obj, "optimization")
      call json%add(root, main_obj)

      ! A run that stopped on the step cap still writes a geometry and an
      ! energy, so `converged` is the only field that distinguishes the two.
      call json%add(main_obj, "converged", record%converged)
      call json%add(main_obj, "steps", record%n_steps)
      call json%add(main_obj, "energy_evaluations", record%n_evaluations)
      call json%add(main_obj, "final_energy", record%final_energy)
      call json%add(main_obj, "final_gradient_max", record%final_gradient_max)
      call json%add(main_obj, "gradient_tolerance", record%gradient_tolerance)
      if (allocated(record%coordinates)) call json%add(main_obj, "coordinates", record%coordinates)
      if (allocated(record%algorithm)) call json%add(main_obj, "algorithm", record%algorithm)

      ! Angstrom here, unlike everything internal, and the units are stated in
      ! the document.
      call json%create_object(geom_obj, "final_geometry")
      call json%add(main_obj, geom_obj)
      call json%add(geom_obj, "units", "angstrom")
      call json%add(geom_obj, "n_atoms", n_atoms)
      if (allocated(record%element_numbers)) then
         ! One at a time rather than an array constructor: a character array
         ! constructor pads every element to the longest, so "H" comes back as
         ! "H " and a reader comparing it to "H" finds no match.
         call json%create_array(symbols_arr, "symbols")
         call json%add(geom_obj, symbols_arr)
         do iatom = 1, n_atoms
            call json%add(symbols_arr, "", &
                          trim(element_number_to_symbol(record%element_numbers(iatom))))
         end do
         call json%add(geom_obj, "atomic_numbers", record%element_numbers)
      end if
      if (allocated(record%final_coordinates)) then
         call json%add(geom_obj, "coordinates", &
                       [(to_angstrom(record%final_coordinates(mod(iatom - 1, 3) + 1, &
                                                              (iatom - 1)/3 + 1)), &
                         iatom=1, 3*n_atoms)])
      end if

      if (allocated(record%trajectory)) then
         call json%create_array(traj_arr, "trajectory")
         call json%add(main_obj, traj_arr)
         do istep = 1, size(record%trajectory, 3)
            call json%create_object(frame_obj, "")
            call json%add(traj_arr, frame_obj)
            call json%add(frame_obj, "step", istep)
            call json%add(frame_obj, "energy", record%trajectory_energies(istep))
            call json%add(frame_obj, "coordinates", &
                          [(to_angstrom(record%trajectory(mod(iatom - 1, 3) + 1, &
                                                          (iatom - 1)/3 + 1, istep)), &
                            iatom=1, 3*n_atoms)])
            nullify (frame_obj)
         end do
      end if

      call json%print(root, path)
      status_ok = .not. json%failed()
      call json%destroy(root)

      if (.not. status_ok) then
         call error%set(ERROR_IO, "Could not write the optimization record to "//path)
         return
      end if

      call logger%info("Optimization record written to "//path)

   end subroutine write_optimization_json

   subroutine write_trajectory_xyz(path, record, error)
      !! Write the path the optimization took, as a multi-frame .xyz
      !!
      !! One frame per accepted geometry, in the order they were accepted, with
      !! the energy on each comment line. Nothing is written when the record
      !! carries no trajectory.
      character(len=*), intent(in) :: path
      type(optimization_record_t), intent(in) :: record
      type(error_t), intent(inout) :: error

      integer :: unit, stat, istep, iatom, n_atoms

      if (.not. allocated(record%trajectory)) return
      if (.not. allocated(record%element_numbers)) return

      n_atoms = size(record%element_numbers)

      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call error%set(ERROR_IO, "Could not write the trajectory to "//path)
         return
      end if

      do istep = 1, size(record%trajectory, 3)
         write (unit, "(i0)") n_atoms
         write (unit, "(a,i0,a,f20.12,a)") "step ", istep, "  E = ", &
            record%trajectory_energies(istep), " Hartree"
         do iatom = 1, n_atoms
            write (unit, "(a4,3f20.12)") element_number_to_symbol(record%element_numbers(iatom)), &
               to_angstrom(record%trajectory(1, iatom, istep)), &
               to_angstrom(record%trajectory(2, iatom, istep)), &
               to_angstrom(record%trajectory(3, iatom, istep))
         end do
      end do
      close (unit)

      call logger%info("Trajectory written to "//path)

   end subroutine write_trajectory_xyz

end module mqc_optimization_output
