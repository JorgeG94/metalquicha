!! Orbitals on a regular grid, and the Gaussian cube file they are written to
module mqc_czt_cube
   !! A box around a molecule, orbitals evaluated on it, and the cube format.
   !!
   !! Nothing here knows what kind of orbital it was handed: a column of
   !! coefficients over the atomic orbitals is all `orbitals_on_grid` needs, so
   !! a molecular orbital goes through the same path as a quasi-atomic one.
   !!
   !! Everything is in Bohr, the unit the cube format itself uses. Grid points
   !! are numbered with z fastest and x slowest, which is the order a cube file
   !! lists them in, so a value array can be written as it stands.
   use pic_types, only: dp, int64
   use pic_blas_interfaces, only: pic_gemm
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION, ERROR_IO
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_ao, only: eval_ao_block, AO_POINT_BLOCK
   implicit none
   private

   public :: cube_grid_t
   public :: build_cube_grid
   public :: orbitals_on_grid
   public :: write_cube
   public :: write_cube_orbitals
   public :: cube_square_norm
   public :: orbitals_per_pass
   public :: MAX_CUBE_POINTS
   public :: MAX_ORBITALS_PER_FILE_BYTES

   integer(int64), parameter :: MAX_CUBE_POINTS = 200000000_int64
      !! Largest grid built. One orbital on it is 1.6 GB, so a request past this
      !! is a spacing or margin that was meant to be coarser.

   integer(int64), parameter :: PASS_BYTES = 536870912_int64
      !! Memory `orbitals_on_grid` callers aim to hold at once: a pass of
      !! orbitals is `n_points` by `n_orbitals` reals.

   integer(int64), parameter :: MAX_ORBITALS_PER_FILE_BYTES = 2147483648_int64
      !! Most memory `write_cube_orbitals` callers hold for one file. A cube with
      !! several orbitals needs all of them in memory at once, where one file
      !! per orbital does not.

   integer, parameter :: INDICES_PER_LINE = 10
      !! Orbital numbers on one line of a multi-orbital header.

   type :: cube_grid_t
      !! A box of `n(1)` by `n(2)` by `n(3)` points
      real(dp) :: origin(3) = 0.0_dp   !! Bohr, the first point
      real(dp) :: spacing = 0.0_dp     !! Bohr between neighbours, on every axis
      integer :: n(3) = 0              !! Points along x, y and z
      integer :: n_points = 0          !! `product(n)`
   end type cube_grid_t

contains

   subroutine build_cube_grid(coordinates, margin, spacing, grid, error)
      !! The smallest axis-aligned box that holds every nucleus with `margin`
      !! to spare, at `spacing`
      real(dp), intent(in) :: coordinates(:, :)   !! (3, natm), Bohr
      real(dp), intent(in) :: margin              !! Bohr, not negative
      real(dp), intent(in) :: spacing             !! Bohr, positive
      type(cube_grid_t), intent(out) :: grid
      type(error_t), intent(inout) :: error

      real(dp) :: low(3), high(3)
      real(dp) :: extent
      integer(int64) :: total
      integer :: axis

      if (error%has_error()) return
      if (size(coordinates, 2) < 1) then
         call error%set(ERROR_VALIDATION, "a cube grid needs at least one atom to surround")
         return
      end if
      if (.not. (spacing > 0.0_dp)) then
         call error%set(ERROR_VALIDATION, "cube spacing must be positive, not "// &
                        to_char(spacing)//" Bohr")
         return
      end if
      if (.not. (margin >= 0.0_dp)) then
         call error%set(ERROR_VALIDATION, "cube margin cannot be negative, not "// &
                        to_char(margin)//" Bohr")
         return
      end if

      low = minval(coordinates, dim=2) - margin
      high = maxval(coordinates, dim=2) + margin
      total = 1_int64
      do axis = 1, 3
         extent = high(axis) - low(axis)
         ! Rounded up, so the box is never shorter than asked for. The grid is
         ! capped before the product is formed from `int` factors.
         if (extent/spacing > real(MAX_CUBE_POINTS, dp)) then
            call error%set(ERROR_VALIDATION, "the cube grid is far too large: "// &
                           to_char(extent)//" Bohr at a spacing of "//to_char(spacing)// &
                           " Bohr. Raise the spacing or lower the margin.")
            return
         end if
         grid%n(axis) = int(ceiling(extent/spacing - 1.0e-9_dp)) + 1
         total = total*int(grid%n(axis), int64)
         if (total > MAX_CUBE_POINTS) then
            call error%set(ERROR_VALIDATION, "the cube grid would exceed "// &
                           to_char(int(MAX_CUBE_POINTS/1000000_int64))// &
                           " million points at a spacing of "//to_char(spacing)// &
                           " Bohr and a margin of "//to_char(margin)// &
                           " Bohr. Raise the spacing or lower the margin.")
            return
         end if
      end do
      grid%origin = low
      grid%spacing = spacing
      grid%n_points = int(total)
   end subroutine build_cube_grid

   pure function orbitals_per_pass(grid) result(count)
      !! How many orbitals `orbitals_on_grid` should be given at a time to stay
      !! near half a gigabyte; at least one
      type(cube_grid_t), intent(in) :: grid
      integer :: count

      count = int(max(1_int64, PASS_BYTES/(8_int64*max(1_int64, int(grid%n_points, int64)))))
   end function orbitals_per_pass

   subroutine orbitals_on_grid(mol, grid, orbitals, values, error)
      !! `sum_mu C(mu, i) chi_mu(r)` at every grid point, for every column `i`
      !!
      !! The basis functions are evaluated by `eval_ao_block` in libcint's own
      !! normalisation and angular conventions, in blocks of points, so the
      !! full `n_points` by `n_ao` array is never held.
      type(czt_molecule_t), intent(in) :: mol
      type(cube_grid_t), intent(in) :: grid
      real(dp), intent(in) :: orbitals(:, :)           !! (n_ao, n_orb)
      real(dp), allocatable, intent(out) :: values(:, :)
         !! (n_points, n_orb), z fastest
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: points(:, :), ao(:, :), block(:, :)
      integer :: first, last, count, p, ix, iy, iz, plane

      if (error%has_error()) return
      if (size(orbitals, 1) /= mol%nao) then
         call error%set(ERROR_VALIDATION, "orbital coefficients run over "// &
                        to_char(size(orbitals, 1))//" basis functions but the molecule "// &
                        "has "//to_char(mol%nao))
         return
      end if

      allocate (values(grid%n_points, size(orbitals, 2)))
      allocate (points(3, AO_POINT_BLOCK))
      plane = grid%n(2)*grid%n(3)

      first = 1
      do while (first <= grid%n_points)
         last = min(first + AO_POINT_BLOCK - 1, grid%n_points)
         count = last - first + 1
         do p = 0, count - 1
            ix = (first + p - 1)/plane
            iy = mod((first + p - 1)/grid%n(3), grid%n(2))
            iz = mod(first + p - 1, grid%n(3))
            points(1, p + 1) = grid%origin(1) + grid%spacing*real(ix, dp)
            points(2, p + 1) = grid%origin(2) + grid%spacing*real(iy, dp)
            points(3, p + 1) = grid%origin(3) + grid%spacing*real(iz, dp)
         end do
         call eval_ao_block(mol, points(:, 1:count), ao, error)
         if (error%has_error()) return
         allocate (block(count, size(orbitals, 2)))
         call pic_gemm(ao, orbitals, block)
         values(first:last, :) = block
         deallocate (block)
         first = last + 1
      end do
   end subroutine orbitals_on_grid

   pure function cube_square_norm(grid, values) result(norm)
      !! The integral of the squared function over the box, as the sum of
      !! squares times the volume of one voxel
      type(cube_grid_t), intent(in) :: grid
      real(dp), intent(in) :: values(:)   !! (n_points)
      real(dp) :: norm

      norm = sum(values**2)*grid%spacing**3
   end function cube_square_norm

   subroutine write_cube(path, grid, atomic_numbers, coordinates, values, title, &
                         subtitle, error)
      !! One function on `grid` as a Gaussian cube file
      !!
      !! The header is two comment lines, the atom count with the grid origin,
      !! one line per axis with its point count and step vector, and one line
      !! per atom with its atomic number, a charge column that repeats it, and
      !! its position. Values follow six to a line, with a new line at the end
      !! of every run of z. The atom count is positive, which tells a reader
      !! there is no orbital-index record after the atoms.
      character(len=*), intent(in) :: path
      type(cube_grid_t), intent(in) :: grid
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)   !! (3, natm), Bohr
      real(dp), intent(in), target, contiguous :: values(:)
         !! (n_points), z fastest
      character(len=*), intent(in) :: title, subtitle
         !! The two comment lines. Anything past a line break is dropped.
      type(error_t), intent(inout) :: error

      real(dp), pointer :: as_column(:, :)

      as_column(1:size(values), 1:1) => values
      call write_cube_orbitals(path, grid, atomic_numbers, coordinates, as_column, &
                               [0], title, subtitle, error, plain=.true.)
   end subroutine write_cube

   subroutine write_cube_orbitals(path, grid, atomic_numbers, coordinates, values, &
                                  orbital_numbers, title, subtitle, error, plain)
      !! Several functions on `grid` as one Gaussian cube file
      !!
      !! The multi-orbital layout: the atom count is written *negative*, which
      !! is how a reader knows an orbital record follows the atoms; that record
      !! is the number of orbitals and then their numbers, ten to a line. The
      !! data block holds every orbital at one grid point before moving to the
      !! next, so the orbital index is fastest, then z, then y, then x, six
      !! values to a line with a new line at the end of every run of z.
      character(len=*), intent(in) :: path
      type(cube_grid_t), intent(in) :: grid
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)   !! (3, natm), Bohr
      real(dp), intent(in) :: values(:, :)        !! (n_points, n_orbitals)
      integer, intent(in) :: orbital_numbers(:)
         !! (n_orbitals) what the file calls each orbital; a viewer shows these
      character(len=*), intent(in) :: title, subtitle
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: plain
         !! One orbital written as a single-function cube: positive atom count
         !! and no orbital record. `write_cube` is this.

      integer :: unit, status, iatom, ix, iy, first, last, n_orb, k, row
      logical :: single
      character(len=256) :: message

      if (error%has_error()) return
      single = .false.
      if (present(plain)) single = plain
      n_orb = size(values, 2)
      if (size(values, 1) /= grid%n_points) then
         call error%set(ERROR_VALIDATION, "cube values hold "//to_char(size(values, 1))// &
                        " points but the grid has "//to_char(grid%n_points))
         return
      end if
      if (size(atomic_numbers) /= size(coordinates, 2)) then
         call error%set(ERROR_VALIDATION, "cube file given "//to_char(size(atomic_numbers))// &
                        " atomic numbers for "//to_char(size(coordinates, 2))//" positions")
         return
      end if
      if (n_orb < 1) then
         call error%set(ERROR_VALIDATION, "a cube file needs at least one orbital")
         return
      end if
      if (.not. single .and. size(orbital_numbers) /= n_orb) then
         call error%set(ERROR_VALIDATION, "cube file given "//to_char(size(orbital_numbers))// &
                        " orbital numbers for "//to_char(n_orb)//" orbitals")
         return
      end if

      open (newunit=unit, file=path, status="replace", action="write", &
            iostat=status, iomsg=message)
      if (status /= 0) then
         call error%set(ERROR_IO, "cannot write cube file '"//path//"': "//trim(message))
         return
      end if

      write (unit, "(a)") one_line(title)
      write (unit, "(a)") one_line(subtitle)
      if (single) then
         write (unit, "(i5,3f12.6)") size(atomic_numbers), grid%origin
      else
         write (unit, "(i5,3f12.6)") - size(atomic_numbers), grid%origin
      end if
      write (unit, "(i5,3f12.6)") grid%n(1), grid%spacing, 0.0_dp, 0.0_dp
      write (unit, "(i5,3f12.6)") grid%n(2), 0.0_dp, grid%spacing, 0.0_dp
      write (unit, "(i5,3f12.6)") grid%n(3), 0.0_dp, 0.0_dp, grid%spacing
      do iatom = 1, size(atomic_numbers)
         write (unit, "(i5,4f12.6)") atomic_numbers(iatom), real(atomic_numbers(iatom), dp), &
            coordinates(:, iatom)
      end do
      if (.not. single) then
         write (unit, "(i5,9(i5))") n_orb, orbital_numbers(1:min(n_orb, INDICES_PER_LINE - 1))
         do k = INDICES_PER_LINE, n_orb, INDICES_PER_LINE
            write (unit, "(10(i5))") orbital_numbers(k:min(k + INDICES_PER_LINE - 1, n_orb))
         end do
      end if

      status = 0
      do ix = 0, grid%n(1) - 1
         do iy = 0, grid%n(2) - 1
            first = (ix*grid%n(2) + iy)*grid%n(3) + 1
            last = first + grid%n(3) - 1
            ! Transposed so the orbital runs fastest along the record.
            write (unit, "(6(1p,e13.5))", iostat=status, iomsg=message) &
               ((values(row, k), k=1, n_orb), row=first, last)
            if (status /= 0) exit
         end do
         if (status /= 0) exit
      end do
      close (unit)
      if (status /= 0) then
         call error%set(ERROR_IO, "writing cube file '"//path//"' failed: "//trim(message))
      end if
   end subroutine write_cube_orbitals

   pure function one_line(text) result(clean)
      !! `text` up to its first line break
      character(len=*), intent(in) :: text
      character(len=:), allocatable :: clean

      integer :: cut

      cut = scan(text, achar(10)//achar(13))
      if (cut > 0) then
         clean = text(1:cut - 1)
      else
         clean = text
      end if
   end function one_line

end module mqc_czt_cube
