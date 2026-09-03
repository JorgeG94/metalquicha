!! Broadcast of a whole system geometry from rank 0
module mqc_bcast_system
   !! Sends a `system_geometry_t` to every rank, component by component.
   !!
   !! A worker builds each fragment from its own copy of the system, so a list
   !! of monomer indices is not enough on its own. The ordinary program never
   !! needs this -- every rank parses the same input file -- but a run driven
   !! from outside gives its workers no filename to parse.
   !!
   !! Nothing is derived on arrival: anything a worker reads is something rank
   !! 0 sent it. Two-dimensional arrays travel as their flattening plus both
   !! extents, so the receiving side never assumes a shape.
   use pic_types, only: dp, int32
   use pic_mpi_lib, only: comm_t, bcast
   use mqc_physical_fragment, only: system_geometry_t, bond_t
   use mqc_bcast, only: bcast_integer_array, bcast_real_array
   implicit none
   private

   public :: bcast_system_geometry

   integer, parameter :: N_SCALARS = 5
      !! n_monomers, atoms_per_monomer, total_atoms, charge, multiplicity,
      !! sent as one message.

contains

   subroutine bcast_system_geometry(comm, sys_geom, root)
      !! Give every rank rank 0's system
      !!
      !! Collective: every rank calls it. On `root` the argument is read; on
      !! every other rank it is overwritten, whatever it held before.
      ! TODO(mqc): `cap_scale` and `fragment_potentials` are the two components
      ! of `system_geometry_t` that never travel. A worker therefore caps broken
      ! bonds at the default scale of 1 whatever the deck asked for, and an EFP
      ! system arrives with no potentials at all.
      type(comm_t), intent(in) :: comm
      type(system_geometry_t), intent(inout) :: sys_geom
      integer(int32), intent(in) :: root

      integer(int32) :: scalars(N_SCALARS)
      logical :: is_root

      is_root = (comm%rank() == root)

      ! The five scalars go together: one broadcast rather than five, and they
      ! are the sizes everything below is checked against.
      if (is_root) then
         scalars = [int(sys_geom%n_monomers, int32), &
                    int(sys_geom%atoms_per_monomer, int32), &
                    int(sys_geom%total_atoms, int32), &
                    int(sys_geom%charge, int32), &
                    int(sys_geom%multiplicity, int32)]
      else
         scalars = 0
      end if
      call bcast(comm, scalars(1), int(N_SCALARS, int32), root)

      if (.not. is_root) then
         call sys_geom%destroy()
         sys_geom%n_monomers = scalars(1)
         sys_geom%atoms_per_monomer = scalars(2)
         sys_geom%total_atoms = scalars(3)
         sys_geom%charge = scalars(4)
         sys_geom%multiplicity = scalars(5)
      end if

      call bcast_int_1d(comm, sys_geom%element_numbers, root, is_root)
      call bcast_real_2d(comm, sys_geom%coordinates, root, is_root)
      call bcast_int_1d(comm, sys_geom%fragment_sizes, root, is_root)
      call bcast_int_2d(comm, sys_geom%fragment_atoms, root, is_root)
      call bcast_int_1d(comm, sys_geom%fragment_charges, root, is_root)
      call bcast_int_1d(comm, sys_geom%fragment_multiplicities, root, is_root)
      call bcast_bonds(comm, sys_geom%bonds, root, is_root)
   end subroutine bcast_system_geometry

   subroutine bcast_int_1d(comm, array, root, is_root)
      !! One default-integer array, allocated on arrival if it was on root
      type(comm_t), intent(in) :: comm
      integer, allocatable, intent(inout) :: array(:)
      integer(int32), intent(in) :: root
      logical, intent(in) :: is_root

      integer(int32), allocatable :: buffer(:)
      integer(int32) :: present_flag

      ! An unallocated component and an empty one are different states, and the
      ! rest of the code distinguishes them, so the flag travels separately
      ! from the length.
      present_flag = 0
      if (is_root) then
         if (allocated(array)) present_flag = 1
      end if
      call bcast(comm, present_flag, 1_int32, root)
      if (present_flag == 0) then
         if (.not. is_root .and. allocated(array)) deallocate (array)
         return
      end if

      if (is_root) then
         allocate (buffer(size(array)))
         buffer = int(array, int32)
      else
         allocate (buffer(0))
      end if
      call bcast_integer_array(comm, buffer, root)

      if (.not. is_root) then
         if (allocated(array)) deallocate (array)
         allocate (array(size(buffer)))
         array = int(buffer)
      end if
      deallocate (buffer)
   end subroutine bcast_int_1d

   subroutine bcast_int_2d(comm, array, root, is_root)
      !! One default-integer matrix, with both extents
      type(comm_t), intent(in) :: comm
      integer, allocatable, intent(inout) :: array(:, :)
      integer(int32), intent(in) :: root
      logical, intent(in) :: is_root

      integer(int32), allocatable :: buffer(:)
      integer(int32) :: shape_info(3)

      shape_info = 0
      if (is_root) then
         if (allocated(array)) then
            shape_info = [1_int32, int(size(array, 1), int32), int(size(array, 2), int32)]
         end if
      end if
      call bcast(comm, shape_info(1), 3_int32, root)
      if (shape_info(1) == 0) then
         if (.not. is_root .and. allocated(array)) deallocate (array)
         return
      end if

      if (is_root) then
         allocate (buffer(size(array)))
         buffer = int(reshape(array, [size(array)]), int32)
      else
         allocate (buffer(0))
      end if
      call bcast_integer_array(comm, buffer, root)

      if (.not. is_root) then
         if (allocated(array)) deallocate (array)
         allocate (array(shape_info(2), shape_info(3)))
         array = int(reshape(buffer, [shape_info(2), shape_info(3)]))
      end if
      deallocate (buffer)
   end subroutine bcast_int_2d

   subroutine bcast_real_2d(comm, array, root, is_root)
      !! One double matrix, with both extents -- the coordinates, in Bohr
      type(comm_t), intent(in) :: comm
      real(dp), allocatable, intent(inout) :: array(:, :)
      integer(int32), intent(in) :: root
      logical, intent(in) :: is_root

      real(dp), allocatable :: buffer(:)
      integer(int32) :: shape_info(3)

      shape_info = 0
      if (is_root) then
         if (allocated(array)) then
            shape_info = [1_int32, int(size(array, 1), int32), int(size(array, 2), int32)]
         end if
      end if
      call bcast(comm, shape_info(1), 3_int32, root)
      if (shape_info(1) == 0) then
         if (.not. is_root .and. allocated(array)) deallocate (array)
         return
      end if

      if (is_root) then
         allocate (buffer(size(array)))
         buffer = reshape(array, [size(array)])
      else
         allocate (buffer(0))
      end if
      call bcast_real_array(comm, buffer, root)

      if (.not. is_root) then
         if (allocated(array)) deallocate (array)
         allocate (array(shape_info(2), shape_info(3)))
         array = reshape(buffer, [shape_info(2), shape_info(3)])
      end if
      deallocate (buffer)
   end subroutine bcast_real_2d

   subroutine bcast_bonds(comm, bonds, root, is_root)
      !! The bond list, packed four integers to a bond
      !!
      !! Component by component rather than as raw bytes: a derived type's
      !! layout is the compiler's business and need not match across a
      !! heterogeneous job.
      type(comm_t), intent(in) :: comm
      type(bond_t), allocatable, intent(inout) :: bonds(:)
      integer(int32), intent(in) :: root
      logical, intent(in) :: is_root

      integer(int32), allocatable :: packed(:)
      integer(int32) :: present_flag, n, ibond

      present_flag = 0
      if (is_root) then
         if (allocated(bonds)) present_flag = 1
      end if
      call bcast(comm, present_flag, 1_int32, root)
      if (present_flag == 0) then
         if (.not. is_root .and. allocated(bonds)) deallocate (bonds)
         return
      end if

      if (is_root) then
         n = int(size(bonds), int32)
         allocate (packed(4*n))
         do ibond = 1, n
            packed(4*(ibond - 1) + 1) = int(bonds(ibond)%atom_i, int32)
            packed(4*(ibond - 1) + 2) = int(bonds(ibond)%atom_j, int32)
            packed(4*(ibond - 1) + 3) = int(bonds(ibond)%order, int32)
            packed(4*(ibond - 1) + 4) = merge(1_int32, 0_int32, bonds(ibond)%is_broken)
         end do
      else
         allocate (packed(0))
      end if
      call bcast_integer_array(comm, packed, root)

      if (.not. is_root) then
         n = int(size(packed)/4, int32)
         if (allocated(bonds)) deallocate (bonds)
         allocate (bonds(n))
         do ibond = 1, n
            bonds(ibond)%atom_i = int(packed(4*(ibond - 1) + 1))
            bonds(ibond)%atom_j = int(packed(4*(ibond - 1) + 2))
            bonds(ibond)%order = int(packed(4*(ibond - 1) + 3))
            bonds(ibond)%is_broken = (packed(4*(ibond - 1) + 4) /= 0)
         end do
      end if
      deallocate (packed)
   end subroutine bcast_bonds

end module mqc_bcast_system
