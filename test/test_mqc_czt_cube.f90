!! Orbitals on a grid, and the Gaussian cube files they are written to
module test_mqc_czt_cube
   !! The grid, the cube writer and the quasi-atomic cube output that sits on
   !! them.
   !!
   !! The check that matters physically is normalisation: an orthonormal orbital
   !! integrates to one over the box, so a sum of squares times the voxel volume
   !! that comes out wrong means the basis functions were evaluated in a
   !! different normalisation or angular convention from the one the orbital
   !! coefficients were built in. Core orbitals are left out of it on purpose --
   !! an oxygen 1s is narrower than a 0.2 bohr voxel, so its quadrature error
   !! says nothing about the evaluation.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_cube_config, only: cube_config_t, CUBE_LISTED, CUBE_OCCUPIED, CUBE_ALL
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_bonding, only: run_quao_analysis
   use mqc_czt_cube, only: cube_grid_t, build_cube_grid, orbitals_on_grid, &
                           write_cube, write_cube_orbitals, cube_square_norm, &
                           orbitals_per_pass
   implicit none
   private

   public :: collect_mqc_czt_cube_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [3, 3])

   real(dp), parameter :: NORM_TOL = 2.0e-2_dp
      !! The quadrature error of a 0.2 bohr grid on a valence orbital

contains

   subroutine collect_mqc_czt_cube_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("grid_covers_the_molecule", test_grid_extent), &
                  new_unittest("grid_refuses_nonsense", test_grid_refusals), &
                  new_unittest("valence_orbitals_normalise", test_normalisation), &
                  new_unittest("cube_round_trips", test_round_trip), &
                  new_unittest("many_orbitals_in_one_cube", test_multi_orbital), &
                  new_unittest("writer_refuses_bad_input", test_writer_refusals), &
                  new_unittest("quao_cubes_are_written", test_quao_cubes) &
                  ]
   end subroutine collect_mqc_czt_cube_tests

   subroutine test_grid_extent(error)
      !! The box holds every atom with the margin to spare, on a lattice of the
      !! requested spacing
      type(error_type), allocatable, intent(out) :: error
      type(cube_grid_t) :: grid
      type(error_t) :: err
      real(dp) :: top(3)

      call build_cube_grid(WATER, 4.0_dp, 0.2_dp, grid, err)
      call check(error,.not. err%has_error(), "a water box should build")
      if (allocated(error)) return

      call check(error, all(abs(grid%origin - (minval(WATER, dim=2) - 4.0_dp)) < 1.0e-12_dp), &
                 "the origin should sit one margin below the lowest atom")
      if (allocated(error)) return
      top = grid%origin + grid%spacing*real(grid%n - 1, dp)
      call check(error, all(top >= maxval(WATER, dim=2) + 4.0_dp - 1.0e-9_dp), &
                 "the far corner should reach one margin past the highest atom")
      if (allocated(error)) return
      call check(error, all(top < maxval(WATER, dim=2) + 4.0_dp + 0.2_dp), &
                 "and overshoot by less than one voxel")
      if (allocated(error)) return
      call check(error, grid%n_points == product(grid%n), "n_points is the product")
      if (allocated(error)) return
      call check(error, orbitals_per_pass(grid) >= 1, "a pass holds at least one orbital")
   end subroutine test_grid_extent

   subroutine test_grid_refusals(error)
      !! Spacing and margin that make no grid are refused rather than rounded
      type(error_type), allocatable, intent(out) :: error
      type(cube_grid_t) :: grid
      type(error_t) :: err

      call build_cube_grid(WATER, 4.0_dp, 0.0_dp, grid, err)
      call check(error, err%has_error(), "zero spacing should be refused")
      if (allocated(error)) return
      call err%clear()

      call build_cube_grid(WATER, 4.0_dp, -0.2_dp, grid, err)
      call check(error, err%has_error(), "negative spacing should be refused")
      if (allocated(error)) return
      call err%clear()

      call build_cube_grid(WATER, -1.0_dp, 0.2_dp, grid, err)
      call check(error, err%has_error(), "negative margin should be refused")
      if (allocated(error)) return
      call err%clear()

      call build_cube_grid(WATER, 4.0_dp, 1.0e-3_dp, grid, err)
      call check(error, err%has_error(), "a grid of billions of points should be refused")
   end subroutine test_grid_refusals

   subroutine water_scf(basis, mol, scf, err)
      character(len=*), intent(in) :: basis
      type(czt_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, basis, mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, 10, 200, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err)
   end subroutine water_scf

   subroutine test_normalisation(error)
      !! Occupied valence molecular orbitals integrate to one and to nothing
      !! against each other
      !!
      !! In cc-pVDZ, so the d shells go through the spherical transform too.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(cube_grid_t) :: grid
      type(error_t) :: err
      real(dp), allocatable :: values(:, :)
      real(dp) :: cross
      integer :: imo

      call water_scf("cc-pvdz", mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "water should converge")
      if (allocated(error)) return

      call build_cube_grid(WATER, 4.0_dp, 0.2_dp, grid, err)
      call orbitals_on_grid(mol, grid, scf%orbitals(:, 2:5), values, err)
      call check(error,.not. err%has_error(), "the orbitals should evaluate")
      if (allocated(error)) return
      call check(error, size(values, 1) == grid%n_points .and. size(values, 2) == 4, &
                 "values should be (n_points, n_orbitals)")
      if (allocated(error)) return

      do imo = 1, 4
         call check(error, abs(cube_square_norm(grid, values(:, imo)) - 1.0_dp) < NORM_TOL, &
                    "a valence molecular orbital should integrate to one over the box")
         if (allocated(error)) return
      end do
      cross = sum(values(:, 1)*values(:, 4))*grid%spacing**3
      call check(error, abs(cross) < NORM_TOL, "different orbitals should be orthogonal")
      call mol%destroy()
   end subroutine test_normalisation

   subroutine test_round_trip(error)
      !! What is written is what is read: the header, every atom, every value
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: PATH = "test_mqc_czt_cube_round_trip.cube"
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(cube_grid_t) :: grid
      type(error_t) :: err
      real(dp), allocatable :: values(:, :), back(:)
      character(len=200) :: title, subtitle
      integer :: n(3)
      integer :: unit, natm, z, iatom, ix, iy, lines
      real(dp) :: origin(3), step(3, 3), position(3)
      real(dp) :: charge, worst

      call water_scf("6-31g", mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "water should converge")
      if (allocated(error)) return
      ! Not a multiple of six points along z, so the last line of every run is short.
      call build_cube_grid(WATER, 3.0_dp, 0.37_dp, grid, err)
      call check(error, mod(grid%n(3), 6) /= 0, "the test needs a ragged last line")
      if (allocated(error)) return
      call orbitals_on_grid(mol, grid, scf%orbitals(:, 5:5), values, err)
      call write_cube(PATH, grid, WATER_Z, WATER, values(:, 1), "first line", &
                      "second line", err)
      call check(error,.not. err%has_error(), "the file should write")
      if (allocated(error)) return

      open (newunit=unit, file=PATH, status="old", action="read")
      read (unit, "(a)") title
      read (unit, "(a)") subtitle
      call check(error, trim(title) == "first line" .and. trim(subtitle) == "second line", &
                 "the two comment lines should come back")
      if (allocated(error)) return
      read (unit, *) natm, origin
      call check(error, natm == 3, "the header should carry the atom count")
      if (allocated(error)) return
      call check(error, all(abs(origin - grid%origin) < 1.0e-5_dp), "origin should match")
      if (allocated(error)) return
      do ix = 1, 3
         read (unit, *) n(ix), step(:, ix)
      end do
      call check(error, all(n == grid%n), "the point counts should match")
      if (allocated(error)) return
      call check(error, abs(step(1, 1) - grid%spacing) < 1.0e-5_dp .and. &
                 abs(step(2, 2) - grid%spacing) < 1.0e-5_dp .and. &
                 abs(step(3, 3) - grid%spacing) < 1.0e-5_dp .and. &
                 abs(step(2, 1)) + abs(step(3, 1)) + abs(step(1, 2)) + abs(step(3, 2)) &
                 + abs(step(1, 3)) + abs(step(2, 3)) < 1.0e-12_dp, &
                 "the axes should be orthogonal with the grid spacing")
      if (allocated(error)) return
      do iatom = 1, natm
         read (unit, *) z, charge, position
         call check(error, z == WATER_Z(iatom) .and. &
                    all(abs(position - WATER(:, iatom)) < 1.0e-5_dp), &
                    "atoms should come back in order with their positions")
         if (allocated(error)) return
      end do

      allocate (back(grid%n_points))
      read (unit, *) back
      close (unit)
      worst = maxval(abs(back - values(:, 1)))/maxval(abs(values(:, 1)))
      call check(error, worst < 2.0e-5_dp, &
                 "values should survive five decimals in the mantissa")
      if (allocated(error)) return
      call check(error, abs(cube_square_norm(grid, back) - cube_square_norm(grid, values(:, 1))) &
                 < 1.0e-4_dp, "and so should the integral of their squares")
      if (allocated(error)) return

      ! z is the fastest index: the first run of z at the corner is the first
      ! n(3) values, which the read above already put in order. Count the lines
      ! to see that every run of z starts a fresh one.
      lines = 0
      open (newunit=unit, file=PATH, status="old", action="read")
      do
         read (unit, "(a)", end=10) title
         lines = lines + 1
      end do
10    close (unit, status="delete")
      iy = (grid%n(3) + 5)/6
      call check(error, lines == 6 + natm + grid%n(1)*grid%n(2)*iy, &
                 "every run of z should start a new line, six values to a line")
      call mol%destroy()
   end subroutine test_round_trip

   subroutine read_cube(path, natm, n_grid, indices, data)
      !! Parse a cube file, single- or multi-orbital, into `data(orbital, point)`
      !!
      !! Reads the data block as a stream, which is all a viewer's parser does,
      !! so it cannot depend on where the writer put its line breaks.
      character(len=*), intent(in) :: path
      integer, intent(out) :: natm
         !! As written: negative when the file holds an orbital record
      integer, intent(out) :: n_grid(3)
      integer, allocatable, intent(out) :: indices(:)
      real(dp), allocatable, intent(out) :: data(:, :)

      character(len=200) :: text
      real(dp) :: origin(3), step(3), position(3)
      real(dp) :: charge
      integer :: unit, ix, nmo, z, iatom

      open (newunit=unit, file=path, status="old", action="read")
      read (unit, "(a)") text
      read (unit, "(a)") text
      read (unit, *) natm, origin
      do ix = 1, 3
         read (unit, *) n_grid(ix), step
      end do
      do iatom = 1, abs(natm)
         read (unit, *) z, charge, position
      end do
      nmo = 1
      if (natm < 0) read (unit, *) nmo
      allocate (indices(nmo))
      indices = 0
      if (natm < 0) then
         ! The count was consumed above, so the numbers are read from a rewound copy.
         rewind (unit)
         do ix = 1, 2 + 1 + 3 + abs(natm)
            read (unit, "(a)") text
         end do
         read (unit, *) nmo, indices
      end if
      allocate (data(nmo, product(n_grid)))
      read (unit, *) data
      close (unit)
   end subroutine read_cube

   subroutine test_multi_orbital(error)
      !! The multi-orbital layout: a negative atom count, the orbital record,
      !! and the orbital index running fastest through the values
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: MULTI = "test_mqc_czt_cube_multi.cube"
      character(len=*), parameter :: SINGLE = "test_mqc_czt_cube_single.cube"
      integer, parameter :: N_ORB = 12
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(cube_grid_t) :: grid
      type(error_t) :: err
      real(dp), allocatable :: values(:, :), data(:, :), one(:, :)
      integer, allocatable :: indices(:), one_indices(:)
      integer :: n_grid(3), numbers(N_ORB)
      integer :: natm, k

      call water_scf("6-31g", mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "water should converge")
      if (allocated(error)) return
      call build_cube_grid(WATER, 2.5_dp, 0.5_dp, grid, err)
      call orbitals_on_grid(mol, grid, scf%orbitals(:, 1:N_ORB), values, err)
      call check(error,.not. err%has_error(), "the orbitals should evaluate")
      if (allocated(error)) return

      ! Twelve orbitals, so the record of numbers has to continue on a second line.
      numbers = [(2*k + 1, k=1, N_ORB)]
      call write_cube_orbitals(MULTI, grid, WATER_Z, WATER, values, numbers, "multi", &
                               "orbital numbers", err)
      call check(error,.not. err%has_error(), "the multi-orbital file should write")
      if (allocated(error)) return
      call read_cube(MULTI, natm, n_grid, indices, data)
      call execute_command_line("rm -f "//MULTI)
      call check(error, natm == -3, "a multi-orbital cube carries a negative atom count")
      if (allocated(error)) return
      call check(error, all(n_grid == grid%n), "the point counts should match")
      if (allocated(error)) return
      call check(error, size(data, 1) == N_ORB .and. all(indices == numbers), &
                 "the orbital record should come back, continuation line included")
      if (allocated(error)) return
      call check(error, maxval(abs(data - transpose(values))) &
                 < 2.0e-5_dp*maxval(abs(values)), &
                 "orbital k at point p should be the value written for it")
      if (allocated(error)) return

      ! Against the same orbital written alone.
      call write_cube(SINGLE, grid, WATER_Z, WATER, values(:, 5), "single", "one", err)
      call read_cube(SINGLE, natm, n_grid, one_indices, one)
      call execute_command_line("rm -f "//SINGLE)
      call check(error, natm == 3 .and. size(one, 1) == 1, &
                 "a single-orbital cube keeps a positive atom count")
      if (allocated(error)) return
      call check(error, maxval(abs(one(1, :) - data(5, :))) < 1.0e-12_dp, &
                 "an orbital should read back the same from both layouts")
      call mol%destroy()
   end subroutine test_multi_orbital

   subroutine test_writer_refusals(error)
      !! A path that cannot be opened or a value array of the wrong length is an
      !! error, not a crash or a short file
      type(error_type), allocatable, intent(out) :: error
      type(cube_grid_t) :: grid
      type(error_t) :: err
      real(dp), allocatable :: values(:)

      call build_cube_grid(WATER, 1.0_dp, 0.5_dp, grid, err)
      allocate (values(grid%n_points), source=0.0_dp)

      call write_cube("/no/such/directory/x.cube", grid, WATER_Z, WATER, values, "a", "b", err)
      call check(error, err%has_error(), "an unwritable path should be an error")
      if (allocated(error)) return
      call err%clear()

      call write_cube("test_mqc_czt_cube_unused.cube", grid, WATER_Z, WATER, values(2:), &
                      "a", "b", err)
      call check(error, err%has_error(), "values that miss a point should be an error")
   end subroutine test_writer_refusals

   subroutine count_files(pattern, names, count)
      !! The files matching `pattern`, found by the shell
      character(len=*), intent(in) :: pattern
      character(len=200), allocatable, intent(out) :: names(:)
      integer, intent(out) :: count

      integer :: unit, status
      character(len=200) :: line

      call execute_command_line("ls "//pattern//" > test_mqc_czt_cube_list.txt 2> /dev/null")
      count = 0
      open (newunit=unit, file="test_mqc_czt_cube_list.txt", status="old", action="read")
      do
         read (unit, "(a)", iostat=status) line
         if (status /= 0) exit
         count = count + 1
      end do
      allocate (names(count))
      rewind (unit)
      do status = 1, count
         read (unit, "(a)") names(status)
      end do
      close (unit, status="delete")
   end subroutine count_files

   subroutine test_quao_cubes(error)
      !! The analysis writes the orbitals it was asked for, named by atom and
      !! type, each one normalised on its grid
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: STEM = "test_mqc_czt_cube_quao"
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      type(cube_config_t) :: cube
      character(len=200), allocatable :: names(:)
      character(len=200) :: title, subtitle
      integer :: n(3)
      integer :: count, unit, natm, at_line
      real(dp) :: origin(3)
      real(dp) :: norm

      call water_scf("6-31g", mol, scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "water should converge")
      if (allocated(error)) return

      ! Two named orbitals: exactly two files.
      cube%selection = CUBE_LISTED
      cube%indices = [2, 6]
      call run_quao_analysis(mol, WATER_Z, WATER_SYM, WATER, scf%orbitals, 10, err, &
                             verbose=.false., cube=cube, cube_stem=STEM)
      call check(error,.not. err%has_error(), "the analysis should run")
      if (allocated(error)) return
      call count_files(STEM//"_*.cube", names, count)
      call check(error, count == 2, "two listed orbitals should give two files")
      if (allocated(error)) then
         call execute_command_line("rm -f "//STEM//"_*.cube")
         return
      end if
      call check(error, index(names(1), STEM//"_002_O0_") == 1, &
                 "the first file should be orbital 2, the first after the oxygen core, on atom 0")
      if (allocated(error)) then
         call execute_command_line("rm -f "//STEM//"_*.cube")
         return
      end if
      call check(error, index(names(2), STEM//"_006_H") == 1, &
                 "the second should be orbital 6, on a hydrogen")
      if (allocated(error)) then
         call execute_command_line("rm -f "//STEM//"_*.cube")
         return
      end if

      ! A valence orbital's comment line carries a norm near one.
      open (newunit=unit, file=trim(names(2)), status="old", action="read")
      read (unit, "(a)") title
      read (unit, "(a)") subtitle
      read (unit, *) natm, origin
      read (unit, *) n(1)
      close (unit)
      at_line = index(subtitle, "grid") + 5
      read (subtitle(at_line:), *) norm
      call check(error, natm == 3 .and. n(1) > 10, "the file should be a cube of water")
      if (.not. allocated(error)) then
         call check(error, abs(norm - 1.0_dp) < NORM_TOL, &
                    "the orbital should integrate to one on the grid")
      end if
      call execute_command_line("rm -f "//STEM//"_*.cube")
      if (allocated(error)) return

      ! The same two orbitals in one file: the orbital record, and each orbital
      ! agreeing with the file it would otherwise have had. (Phase is fixed per
      ! orbital, so the two layouts must match in sign too.)
      call read_cube_pair(STEM, mol, scf, error)
      if (allocated(error)) return

      ! All of them, then only the occupied ones.
      cube%selection = CUBE_ALL
      call run_quao_analysis(mol, WATER_Z, WATER_SYM, WATER, scf%orbitals, 10, err, &
                             verbose=.false., cube=cube, cube_stem=STEM)
      call count_files(STEM//"_*.cube", names, count)
      call execute_command_line("rm -f "//STEM//"_*.cube")
      call check(error, count == 6, "water has six quasi-atomic orbitals in 6-31G")
      if (allocated(error)) return

      cube%selection = CUBE_OCCUPIED
      call run_quao_analysis(mol, WATER_Z, WATER_SYM, WATER, scf%orbitals, 10, err, &
                             verbose=.false., cube=cube, cube_stem=STEM)
      call count_files(STEM//"_*.cube", names, count)
      call execute_command_line("rm -f "//STEM//"_*.cube")
      call check(error, count >= 4 .and. count <= 6, &
                 "the occupied selection should leave out at most the empty ones")
      if (allocated(error)) return

      ! An orbital that does not exist writes nothing.
      cube%selection = CUBE_LISTED
      cube%indices = [2, 99]
      call run_quao_analysis(mol, WATER_Z, WATER_SYM, WATER, scf%orbitals, 10, err, &
                             verbose=.false., cube=cube, cube_stem=STEM)
      call count_files(STEM//"_*.cube", names, count)
      call execute_command_line("rm -f "//STEM//"_*.cube")
      call check(error, count == 0, "a list naming a missing orbital should write nothing")
      call mol%destroy()
   end subroutine test_quao_cubes

   subroutine read_cube_pair(stem, mol, scf, error)
      !! Write orbitals 2 and 6 both ways and compare the values
      character(len=*), intent(in) :: stem
      type(czt_molecule_t), intent(in) :: mol
      type(rhf_result_t), intent(in) :: scf
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(cube_config_t) :: cube
      character(len=200), allocatable :: names(:)
      real(dp), allocatable :: separate(:, :), together(:, :), first(:, :)
      integer, allocatable :: indices(:)
      integer :: n_grid(3)
      integer :: count, natm, k
      real(dp) :: norm

      cube%selection = CUBE_LISTED
      cube%indices = [2, 6]
      call run_quao_analysis(mol, WATER_Z, WATER_SYM, WATER, scf%orbitals, 10, err, &
                             verbose=.false., cube=cube, cube_stem=stem)
      call count_files(stem//"_*.cube", names, count)
      call check(error, count == 2, "listed orbitals should give one file each")
      if (allocated(error)) then
         call execute_command_line("rm -f "//stem//"_*.cube")
         return
      end if
      call read_cube(trim(names(1)), natm, n_grid, indices, first)
      call read_cube(trim(names(2)), natm, n_grid, indices, separate)
      call execute_command_line("rm -f "//stem//"_*.cube")

      cube%single_file = .true.
      call run_quao_analysis(mol, WATER_Z, WATER_SYM, WATER, scf%orbitals, 10, err, &
                             verbose=.false., cube=cube, cube_stem=stem)
      call count_files(stem//"_*.cube", names, count)
      call check(error, count == 1, "single_file should give exactly one file")
      if (allocated(error)) then
         call execute_command_line("rm -f "//stem//"_*.cube")
         return
      end if
      call read_cube(trim(names(1)), natm, n_grid, indices, together)
      call execute_command_line("rm -f "//stem//"_*.cube")

      call check(error, natm < 0 .and. size(together, 1) == 2 .and. all(indices == [2, 6]), &
                 "the file should hold orbitals 2 and 6 under a negative atom count")
      if (allocated(error)) return
      call check(error, maxval(abs(together(1, :) - first(1, :))) &
                 < 2.0e-5_dp*maxval(abs(first)) .and. &
                 maxval(abs(together(2, :) - separate(1, :))) &
                 < 2.0e-5_dp*maxval(abs(separate)), &
                 "each orbital should have the same values in both layouts")
      if (allocated(error)) return
      do k = 1, 2
         norm = sum(together(k, :)**2)*(0.2_dp**3)
         call check(error, abs(norm - 1.0_dp) < NORM_TOL, &
                    "each orbital in the shared file should still integrate to one")
         if (allocated(error)) return
      end do
   end subroutine read_cube_pair

end module test_mqc_czt_cube

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_cube, only: collect_mqc_czt_cube_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_cube", collect_mqc_czt_cube_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
