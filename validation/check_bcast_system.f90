!! Manual check that a whole system geometry reaches every rank intact
!!
!!     mpirun -np 4 ./build/check_bcast_system
!!
!! Every rank prints one line; a clean run is "ok" from all of them. This is
!! the check that matters most of the broadcast set, because a worker builds
!! every fragment it evaluates out of its own copy of this type. A component
!! that arrives wrong does not fail -- it produces fragments that are subtly
!! not the ones rank 0 asked for, and an energy that looks reasonable.
!!
!! The system below is deliberately awkward: monomers of unequal size, so
!! atoms_per_monomer must come across as 0 rather than a size; a non-zero
!! charge and a non-singlet multiplicity, so neither can be assumed; bonds
!! both broken and preserved, so the flag is not uniform.
program check_bcast_system
   use pic_types, only: dp, int32
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init, pic_mpi_finalize
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_bcast_system, only: bcast_system_geometry
   implicit none

   integer, parameter :: N_ATOMS = 7  !! Atoms in the reference system

   type(comm_t) :: comm
   type(system_geometry_t) :: sys
   integer(int32) :: rank, fails
   integer :: i

   call pic_mpi_init()
   comm = comm_world()
   rank = comm%rank()
   fails = 0

   if (rank == 0) call build_reference(sys)
   call bcast_system_geometry(comm, sys, 0_int32)

   ! Scalars
   call expect(sys%total_atoms == 7, "total_atoms", fails)
   call expect(sys%n_monomers == 3, "n_monomers", fails)
   call expect(sys%atoms_per_monomer == 0, "atoms_per_monomer (variable)", fails)
   call expect(sys%charge == -1, "charge", fails)
   call expect(sys%multiplicity == 2, "multiplicity", fails)

   ! Atoms
   call expect(allocated(sys%element_numbers), "element_numbers allocated", fails)
   if (allocated(sys%element_numbers)) then
      call expect(size(sys%element_numbers) == 7, "element_numbers size", fails)
      call expect(all(sys%element_numbers == [8, 1, 1, 8, 1, 1, 17]), &
                  "element_numbers values", fails)
   end if

   call expect(allocated(sys%coordinates), "coordinates allocated", fails)
   if (allocated(sys%coordinates)) then
      call expect(size(sys%coordinates, 1) == 3 .and. size(sys%coordinates, 2) == 7, &
                  "coordinates shape", fails)
      ! Exact: a broadcast copies bits.
      call expect(maxval(abs(sys%coordinates - reference_coords())) == 0.0_dp, &
                  "coordinates values", fails)
   end if

   ! Partition -- unequal monomers, so the padding must survive too
   call expect(allocated(sys%fragment_sizes), "fragment_sizes allocated", fails)
   if (allocated(sys%fragment_sizes)) then
      call expect(all(sys%fragment_sizes == [3, 3, 1]), "fragment_sizes values", fails)
   end if
   call expect(allocated(sys%fragment_atoms), "fragment_atoms allocated", fails)
   if (allocated(sys%fragment_atoms)) then
      call expect(size(sys%fragment_atoms, 1) == 3 .and. &
                  size(sys%fragment_atoms, 2) == 3, "fragment_atoms shape", fails)
      call expect(all(sys%fragment_atoms(1:3, 1) == [0, 1, 2]), "monomer 1 atoms", fails)
      call expect(sys%fragment_atoms(1, 3) == 6, "monomer 3 first atom", fails)
   end if
   if (allocated(sys%fragment_charges)) then
      call expect(all(sys%fragment_charges == [0, 0, -1]), "fragment_charges", fails)
   end if
   if (allocated(sys%fragment_multiplicities)) then
      call expect(all(sys%fragment_multiplicities == [1, 1, 2]), &
                  "fragment_multiplicities", fails)
   end if

   ! Bonds -- one broken, one not, so a uniform flag would pass wrongly
   call expect(allocated(sys%bonds), "bonds allocated", fails)
   if (allocated(sys%bonds)) then
      call expect(size(sys%bonds) == 2, "bond count", fails)
      call expect(sys%bonds(1)%atom_i == 2 .and. sys%bonds(1)%atom_j == 3, &
                  "bond 1 atoms", fails)
      call expect(sys%bonds(1)%order == 1, "bond 1 order", fails)
      call expect(sys%bonds(1)%is_broken, "bond 1 is broken", fails)
      call expect(.not. sys%bonds(2)%is_broken, "bond 2 is not broken", fails)
   end if

   if (fails == 0) then
      write (*, "(A,I0,A)") "[rank ", rank, "] ok"
   else
      write (*, "(A,I0,A,I0,A)") "[rank ", rank, "] FAILED ", fails, " checks"
   end if

   call sys%destroy()
   call comm%finalize()
   call pic_mpi_finalize()
   if (fails /= 0) error stop 1

contains

   function reference_coords() result(coords)
      real(dp) :: coords(3, N_ATOMS)
      integer :: iatom, ixyz
      do iatom = 1, N_ATOMS
         do ixyz = 1, 3
            coords(ixyz, iatom) = real(iatom, dp) + 0.125_dp*real(ixyz, dp)
         end do
      end do
   end function reference_coords

   subroutine build_reference(s)
      type(system_geometry_t), intent(out) :: s

      s%total_atoms = 7
      s%n_monomers = 3
      s%atoms_per_monomer = 0
      s%charge = -1
      s%multiplicity = 2

      allocate (s%element_numbers(7))
      s%element_numbers = [8, 1, 1, 8, 1, 1, 17]
      allocate (s%coordinates(3, 7))
      s%coordinates = reference_coords()

      allocate (s%fragment_sizes(3))
      s%fragment_sizes = [3, 3, 1]
      allocate (s%fragment_atoms(3, 3))
      s%fragment_atoms = 0
      s%fragment_atoms(1:3, 1) = [0, 1, 2]
      s%fragment_atoms(1:3, 2) = [3, 4, 5]
      s%fragment_atoms(1, 3) = 6
      allocate (s%fragment_charges(3))
      s%fragment_charges = [0, 0, -1]
      allocate (s%fragment_multiplicities(3))
      s%fragment_multiplicities = [1, 1, 2]

      allocate (s%bonds(2))
      s%bonds(1)%atom_i = 2
      s%bonds(1)%atom_j = 3
      s%bonds(1)%order = 1
      s%bonds(1)%is_broken = .true.
      s%bonds(2)%atom_i = 0
      s%bonds(2)%atom_j = 1
      s%bonds(2)%order = 1
      s%bonds(2)%is_broken = .false.
   end subroutine build_reference

   subroutine expect(condition, what, counter)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: what
      integer(int32), intent(inout) :: counter

      if (condition) return
      counter = counter + 1
      write (*, "(A,I0,A,A)") "[rank ", rank, "] wrong: ", what
   end subroutine expect

end program check_bcast_system
