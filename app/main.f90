program gms_basis_app
  use gms_cli_parser
  use gms_xyz_reader
  use gms_basis_reader, only: build_molecular_basis, ang_mom_int_to_char
  use gms_cgto
  use basis_file_reader
  use iso_fortran_env, only: real64
  implicit none

  type(cli_args_type) :: args
  type(geometry_type) :: geom
  type(basis_file_t) :: basis_file
  type(molecular_basis_type) :: mol_basis
  integer :: stat, iatom, ishell
  character(len=:), allocatable :: errmsg
  character(len=:), allocatable :: basis_filename
  character(len=10), allocatable :: element_names(:)

  ! Parse command line arguments
  call parse_command_line(args, stat, errmsg)

  if (stat == -1) then
    ! Help was requested
    stop 0
  else if (stat /= 0) then
    print *, "ERROR: ", errmsg
    stop 1
  end if

  print *
  print *, "=========================================="
  print *, "GAMESS Basis Set Generator"
  print *, "=========================================="
  print *
  print *, "Geometry file:  ", trim(args%xyz_file)
  print *, "Basis set:      ", trim(args%basis_name)
  print *

  ! Step 1: Read geometry
  print *, "Reading geometry..."
  call read_xyz_file(args%xyz_file, geom, stat, errmsg)

  if (stat /= 0) then
    print *, "ERROR reading geometry: ", errmsg
    call args%destroy()
    stop 1
  end if

  print *, "  Number of atoms: ", geom%natoms
  print *

  ! Step 2: Find basis file (tries multiple name variants)
  print *, "Searching for basis set file..."
  call find_basis_file(args%basis_name, basis_filename, stat, errmsg)

  if (stat /= 0) then
    print *, "ERROR: ", errmsg
    call geom%destroy()
    call args%destroy()
    stop 1
  end if

  print *, "  Found: ", trim(basis_filename)
  print *

  print *, "Loading basis set..."
  call open_basis_file(basis_file, basis_filename)

  print *, "  Basis set loaded successfully"
  print *

  ! Step 3: Convert element symbols to full names
  allocate(element_names(geom%natoms))
  do iatom = 1, geom%natoms
    element_names(iatom) = symbol_to_name(geom%elements(iatom))
  end do

  ! Step 4: Build molecular basis
  print *, "Building molecular basis set..."
  call build_molecular_basis(basis_file%data_section, element_names, mol_basis, stat, errmsg)

  if (stat /= 0) then
    print *, "ERROR building basis: ", errmsg
    call geom%destroy()
    call args%destroy()
    stop 1
  end if

  print *, "  Basis set built successfully"
  print *

  ! Step 5: Print summary
  print *, "=========================================="
  print *, "Summary"
  print *, "=========================================="
  print *, "Total atoms:            ", mol_basis%nelements
  print *, "Total basis functions:  ", mol_basis%num_basis_functions()
  print *

  ! Print first few atoms in detail
  print *, "Atom details:"
  do iatom = 1, min(5, mol_basis%nelements)
    print '(2X, I3, 2X, A4, ":", I4, " basis functions from", I3, " shells")', &
      iatom, &
      trim(mol_basis%elements(iatom)%element), &
      mol_basis%elements(iatom)%num_basis_functions(), &
      mol_basis%elements(iatom)%nshells
  end do

  if (mol_basis%nelements > 5) then
    print *, "  ... (", mol_basis%nelements - 5, " more atoms)"
  end if
  print *

  ! Cleanup
  call args%destroy()
  call geom%destroy()
  call mol_basis%destroy()

  print *, "=========================================="
  print *, "SUCCESS"
  print *, "=========================================="

contains

  !> Convert element symbol to full name for basis file lookup
  function symbol_to_name(symbol) result(name)
    character(len=*), intent(in) :: symbol
    character(len=10) :: name
    character(len=10) :: sym_upper

    ! Convert to uppercase
    sym_upper = symbol
    call to_upper(sym_upper)

    select case (trim(sym_upper))
    case ("H")
      name = "HYDROGEN"
    case ("C")
      name = "CARBON"
    case ("N")
      name = "NITROGEN"
    case ("O")
      name = "OXYGEN"
    case ("F")
      name = "FLUORINE"
    case ("P")
      name = "PHOSPHORUS"
    case ("S")
      name = "SULFUR"
    case ("CL")
      name = "CHLORINE"
    case ("BR")
      name = "BROMINE"
    case ("I")
      name = "IODINE"
    case default
      ! If not recognized, return the symbol as-is
      name = sym_upper
    end select
  end function symbol_to_name

  !> Convert string to uppercase
  subroutine to_upper(str)
    character(len=*), intent(inout) :: str
    integer :: i, ic

    do i = 1, len(str)
      ic = iachar(str(i:i))
      if (ic >= iachar('a') .and. ic <= iachar('z')) then
        str(i:i) = achar(ic - 32)
      end if
    end do
  end subroutine to_upper

end program gms_basis_app
