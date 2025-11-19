program main
  use gms_xyz_reader
  use gms_basis_reader, only: classify_line, parse_element_basis, &
                              build_molecular_basis, ang_mom_int_to_char
  use gms_cgto
  use basis_file_reader
  use iso_fortran_env, only: real64
  implicit none
  character(len=*), parameter :: test_basis = &
                                 "$DATA"//new_line('a')// &
                                 ""//new_line('a')// &
                                 "HYDROGEN"//new_line('a')// &
                                 "S 2"//new_line('a')// &
                                 "1 1.0 2.0"//new_line('a')// &
                                 "2 1.0 2.0"//new_line('a')// &
                                 ""//new_line('a')// &
                                 "L 2"//new_line('a')// &
                                 "1 1.0 2.0 3.0"//new_line('a')// &
                                 "2 4.0 5.0 6.0"//new_line('a')// &
                                 ""//new_line('a')// &
                                 "$END"

  print *, test_basis

  call test_classify_lines()

  call test_parse_element()

  call test_h2_molecule()

  call test_basis_file_reader()

  call test_xyz_reader()

  call test_prism_with_basis()

contains

  subroutine test_classify_lines()
    integer :: line_start, line_end, line_num
    character(len=256) :: line
    integer :: line_type

    print *, "Testing line classification:"
    print *, "============================"

    line_start = 1
    line_num = 0

    do while (line_start <= len(test_basis))
      ! Find next newline
      line_end = index(test_basis(line_start:), new_line('a'))
      if (line_end == 0) then
        line = test_basis(line_start:)
        line_start = len(test_basis) + 1
      else
        line = test_basis(line_start:line_start + line_end - 2)
        line_start = line_start + line_end
      end if

      line_num = line_num + 1
      line = adjustl(line)

      if (len_trim(line) == 0) cycle

      line_type = classify_line(line)

      print '(i3,a,i2,a,a)', line_num, ': Type=', line_type, ' | ', trim(line)
    end do

  end subroutine test_classify_lines

  subroutine test_parse_element()
    type(atomic_basis_type) :: h_basis
    integer :: stat, i, j
    character(len=:), allocatable :: errmsg

    print *, "Testing parse_element_basis:"
    print *, "============================"
    print *

    call parse_element_basis(test_basis, "HYDROGEN", h_basis, stat, errmsg)

    if (stat /= 0) then
      print *, "ERROR: ", errmsg
      return
    end if

    print *, "Successfully parsed basis for: ", h_basis%element
    print *, "Number of shells: ", h_basis%nshells
    print *

    ! Print each shell
    do i = 1, h_basis%nshells
      print '(a,i0,a,a,a,i0,a)', "Shell ", i, " (", ang_mom_int_to_char(h_basis%shells(i)%ang_mom), &
        "): ", h_basis%shells(i)%nfunc, " primitives"
      do j = 1, h_basis%shells(i)%nfunc
        print '(2x,i2,2x,f12.6,2x,f12.6)', j, h_basis%shells(i)%exponents(j), &
          h_basis%shells(i)%coefficients(j)
      end do
      print *
    end do

    call h_basis%destroy()

  end subroutine test_parse_element

  subroutine test_h2_molecule()
    type(molecular_basis_type) :: h2_basis
    integer :: stat, iatom, ishell, ifunc
    character(len=:), allocatable :: errmsg
    character(len=*), dimension(2), parameter :: h2_atoms = ["HYDROGEN", "HYDROGEN"]

    print *, "Testing H2 molecular basis:"
    print *, "============================"
    print *

    call build_molecular_basis(test_basis, h2_atoms, h2_basis, stat, errmsg)

    if (stat /= 0) then
      print *, "ERROR: ", errmsg
      return
    end if

    print *, "Successfully built molecular basis for H2"
    print *, "Number of atoms: ", h2_basis%nelements
    print *, "Total basis functions for H2: ", h2_basis%num_basis_functions()
    print *

    ! Print basis for each atom
    do iatom = 1, h2_basis%nelements
      print '(a,i0,a,a)', "Atom ", iatom, ": ", trim(h2_basis%elements(iatom)%element)
      print '(a,i0,a)', "  Number of shells: ", h2_basis%elements(iatom)%nshells

      do ishell = 1, h2_basis%elements(iatom)%nshells
        print '(a,i0,a,a,a,i0,a)', "  Shell ", ishell, " (", &
          ang_mom_int_to_char(h2_basis%elements(iatom)%shells(ishell)%ang_mom), &
          "): ", h2_basis%elements(iatom)%shells(ishell)%nfunc, " primitives"

        do ifunc = 1, h2_basis%elements(iatom)%shells(ishell)%nfunc
          print '(4x,i2,2x,f12.6,2x,f12.6)', ifunc, &
            h2_basis%elements(iatom)%shells(ishell)%exponents(ifunc), &
            h2_basis%elements(iatom)%shells(ishell)%coefficients(ifunc)
        end do
      end do
      print *
    end do

    call h2_basis%destroy()

  end subroutine test_h2_molecule

  subroutine test_basis_file_reader()
  type(basis_file_t) :: basis_file
  type(atomic_basis_type) :: h_basis, c_basis
  character(len=:), allocatable :: h_content, c_content, errmsg
  character(len=*), parameter :: path_to_basis = "basis_sets/6-31G.txt"
  integer :: stat, i, j
  
  print *, "Testing basis_file_reader:"
  print *, "=========================="
  print *
  
  ! Open the basis set file once
  print *, "Opening basis set file: ", path_to_basis
  call open_basis_file(basis_file, path_to_basis)
  print *, "File opened successfully"
  print *
  
  ! Test 1: Extract and parse Hydrogen
  print *, "--- Extracting Hydrogen ---"
  h_content = extract_element(basis_file, "HYDROGEN")
  print *, "Extracted content:"
  print *, h_content
  print *
  
  call parse_element_basis(h_content, "HYDROGEN", h_basis, stat, errmsg)
  if (stat /= 0) then
    print *, "ERROR parsing Hydrogen: ", errmsg
    return
  end if
  
  print *, "Successfully parsed basis for: ", h_basis%element
  print *, "Number of shells: ", h_basis%nshells
  print *
  
  print *, h_basis%shells(1)%ang_mom
  do i = 1, h_basis%nshells
    print '(a,i0,a,a,a,i0,a)', "Shell ", i, " (", ang_mom_int_to_char(h_basis%shells(i)%ang_mom), &
      "): ", h_basis%shells(i)%nfunc, " primitives"
    do j = 1, h_basis%shells(i)%nfunc
      print '(2x,i2,2x,f12.6,2x,f12.6)', j, h_basis%shells(i)%exponents(j), &
        h_basis%shells(i)%coefficients(j)
    end do
    print *
  end do
  
  ! Test 2: Extract and parse Carbon
  print *, "--- Extracting Carbon ---"
  c_content = extract_element(basis_file, "CARBON")
  print *, "Extracted content:"
  print *, c_content
  print *
  
  call parse_element_basis(c_content, "CARBON", c_basis, stat, errmsg)
  if (stat /= 0) then
    print *, "ERROR parsing Carbon: ", errmsg
    call h_basis%destroy()
    return
  end if
  
  print *, "Successfully parsed basis for: ", c_basis%element
  print *, "Number of shells: ", c_basis%nshells
  print *

  do i = 1, c_basis%nshells
    print '(a,i0,a,a,a,i0,a)', "Shell ", i, " (", ang_mom_int_to_char(c_basis%shells(i)%ang_mom), &
      "): ", c_basis%shells(i)%nfunc, " primitives"
    do j = 1, c_basis%shells(i)%nfunc
      print '(2x,i2,2x,f12.6,2x,f12.6)', j, c_basis%shells(i)%exponents(j), &
        c_basis%shells(i)%coefficients(j)
    end do
    print *
  end do
  
  ! Cleanup
  call h_basis%destroy()
  call c_basis%destroy()
  
  print *, "All tests passed!"
  
end subroutine test_basis_file_reader

subroutine test_xyz_reader()
  
  type(geometry_type) :: geom
  integer :: stat, i
  character(len=:), allocatable :: errmsg
  character(len=*), parameter :: test_xyz = &
    "3" // new_line('a') // &
    "Water molecule" // new_line('a') // &
    "O    0.000000    0.000000    0.119262" // new_line('a') // &
    "H    0.000000    0.763239   -0.477047" // new_line('a') // &
    "H    0.000000   -0.763239   -0.477047"
  
  print *, "=========================================="
  print *, "Testing XYZ Reader"
  print *, "=========================================="
  print *
  
  ! Test 1: Read water molecule
  print *, "Test 1: Reading water molecule"
  print *, "------------------------------------------"
  call read_xyz_string(test_xyz, geom, stat, errmsg)
  
  if (stat /= 0) then
    print *, "FAILED: ", errmsg
    return
  end if
  
  print *, "Number of atoms:", geom%natoms
  print *, "Comment:", trim(geom%comment)
  print *
  print *, "Geometry:"
  do i = 1, geom%natoms
    print '(2X, A4, 3F12.6)', geom%elements(i), geom%coords(:, i)
  end do
  print *
  
  ! Verify results
  if (geom%natoms /= 3) then
    print *, "FAILED: Expected 3 atoms, got", geom%natoms
    return
  end if
  
  if (geom%elements(1) /= "O") then
    print *, "FAILED: Expected O, got", geom%elements(1)
    return
  end if
  
  if (abs(geom%coords(3, 1) - 0.119262_real64) > 1.0e-6_real64) then
    print *, "FAILED: Coordinate mismatch"
    return
  end if
  
  print *, "PASSED"
  print *
  call geom%destroy()
  
  ! Test 2: Single atom
  print *, "Test 2: Single atom (carbon)"
  print *, "------------------------------------------"
  call read_xyz_string("1" // new_line('a') // &
                       "Single carbon atom" // new_line('a') // &
                       "C  1.0  2.0  3.0", geom, stat, errmsg)
  
  if (stat /= 0) then
    print *, "FAILED: ", errmsg
    return
  end if
  
  if (geom%natoms /= 1 .or. geom%elements(1) /= "C") then
    print *, "FAILED: Single atom test"
    return
  end if
  
  print *, "Element:", trim(geom%elements(1))
  print '(2X, A, 3F12.6)', "Coords:", geom%coords(:, 1)
  print *, "PASSED"
  print *
  call geom%destroy()
  
  ! Test 3: Empty comment line
  print *, "Test 3: Empty comment line"
  print *, "------------------------------------------"
  call read_xyz_string("2" // new_line('a') // &
                       new_line('a') // &
                       "He  0.0  0.0  0.0" // new_line('a') // &
                       "Ne  5.0  0.0  0.0", geom, stat, errmsg)
  
  if (stat /= 0) then
    print *, "FAILED: ", errmsg
    return
  end if
  
  print *, "Comment: '", trim(geom%comment), "'"
  print *, "Atoms:", geom%natoms
  print *, "PASSED"
  print *
  call geom%destroy()
  
  ! Test 4: Error handling - insufficient lines
  print *, "Test 4: Error handling (insufficient lines)"
  print *, "------------------------------------------"
  call read_xyz_string("3" // new_line('a') // &
                       "Should fail" // new_line('a') // &
                       "H  0.0  0.0  0.0", geom, stat, errmsg)
  
  if (stat == 0) then
    print *, "FAILED: Should have detected insufficient lines"
    call geom%destroy()
    return
  end if
  
  print *, "Correctly caught error:", trim(errmsg)
  print *, "PASSED"
  print *
  
  ! Test 5: Error handling - invalid atom count
  print *, "Test 5: Error handling (invalid atom count)"
  print *, "------------------------------------------"
  call read_xyz_string("not_a_number" // new_line('a') // &
                       "Comment", geom, stat, errmsg)
  
  if (stat == 0) then
    print *, "FAILED: Should have detected invalid atom count"
    call geom%destroy()
    return
  end if
  
  print *, "Correctly caught error:", trim(errmsg)
  print *, "PASSED"
  print *
  
  ! Test 6: Error handling - malformed coordinate line
  print *, "Test 6: Error handling (malformed coordinates)"
  print *, "------------------------------------------"
  call read_xyz_string("1" // new_line('a') // &
                       "Test" // new_line('a') // &
                       "C  1.0  invalid  3.0", geom, stat, errmsg)
  
  if (stat == 0) then
    print *, "FAILED: Should have detected malformed coordinates"
    call geom%destroy()
    return
  end if
  
  print *, "Correctly caught error:", trim(errmsg)
  print *, "PASSED"
  print *

    ! Test 7: Read from file (prism.xyz)
  print *, "Test 7: Reading from prism.xyz file"
  print *, "------------------------------------------"
  call read_xyz_file("prism.xyz", geom, stat, errmsg)
  
  if (stat /= 0) then
    print *, "WARNING: Could not read prism.xyz (file may not exist)"
    print *, "Error: ", errmsg
    print *, "SKIPPED"
    print *
  else
    print *, "Number of atoms:", geom%natoms
    print *, "Comment:", trim(geom%comment)
    print *
    print *, "Geometry:"
    do i = 1, min(geom%natoms, 10)  ! Print first 10 atoms max
      print '(2X, A4, 3F12.6)', geom%elements(i), geom%coords(:, i)
    end do
    if (geom%natoms > 10) then
      print *, "  ... (", geom%natoms - 10, " more atoms)"
    end if
    print *
    print *, "PASSED"
    print *
    call geom%destroy()
  end if

  print *, "=========================================="
  print *, "All XYZ Reader Tests PASSED!"
  print *, "=========================================="

end subroutine test_xyz_reader

subroutine test_prism_with_basis()
  type(geometry_type) :: geom
  type(basis_file_t) :: basis_file
  type(molecular_basis_type) :: mol_basis
  integer :: stat, iatom, ishell
  character(len=:), allocatable :: errmsg, element_content
  character(len=*), parameter :: xyz_file = "prism.xyz"
  character(len=*), parameter :: basis_file_name = "basis_sets/6-31G.txt"
  character(len=10), allocatable :: element_names(:)

  print *
  print *, "=========================================="
  print *, "Testing Prism Molecule with 6-31G Basis"
  print *, "=========================================="
  print *

  ! Read geometry from XYZ file
  print *, "Step 1: Reading geometry from ", xyz_file
  print *, "------------------------------------------"
  call read_xyz_file(xyz_file, geom, stat, errmsg)

  if (stat /= 0) then
    print *, "WARNING: Could not read ", xyz_file
    print *, "Error: ", errmsg
    print *, "SKIPPED"
    return
  end if

  print *, "Successfully read geometry"
  print *, "  Number of atoms: ", geom%natoms
  print *, "  Comment: ", trim(geom%comment)
  print *

  ! Show first few atoms
  print *, "First 5 atoms:"
  do iatom = 1, min(5, geom%natoms)
    print '(2X, I3, 2X, A4, 3F12.6)', iatom, geom%elements(iatom), geom%coords(:, iatom)
  end do
  if (geom%natoms > 5) then
    print *, "  ... (", geom%natoms - 5, " more atoms)"
  end if
  print *

  ! Convert element symbols to full names for basis lookup
  allocate(element_names(geom%natoms))
  do iatom = 1, geom%natoms
    element_names(iatom) = symbol_to_name(geom%elements(iatom))
  end do

  ! Open basis set file
  print *, "Step 2: Loading basis set from ", basis_file_name
  print *, "------------------------------------------"
  call open_basis_file(basis_file, basis_file_name)
  print *, "Successfully loaded basis set file"
  print *

  ! Build molecular basis
  print *, "Step 3: Building molecular basis set"
  print *, "------------------------------------------"
  call build_molecular_basis(basis_file%data_section, element_names, mol_basis, stat, errmsg)

  if (stat /= 0) then
    print *, "FAILED: ", errmsg
    call geom%destroy()
    return
  end if

  print *, "Successfully built molecular basis set!"
  print *

  ! Print summary statistics
  print *, "=========================================="
  print *, "Basis Set Summary"
  print *, "=========================================="
  print *, "Total number of atoms:      ", mol_basis%nelements
  print *, "Total basis functions:      ", mol_basis%num_basis_functions()
  print *

  ! Count shells by element type
  print *, "Basis functions per atom type:"
  do iatom = 1, min(3, mol_basis%nelements)
    print '(2X, A4, ":", I4, " basis functions from", I3, " shells")', &
      trim(mol_basis%elements(iatom)%element), &
      mol_basis%elements(iatom)%num_basis_functions(), &
      mol_basis%elements(iatom)%nshells
  end do
  if (mol_basis%nelements > 3) then
    print *, "  ... (", mol_basis%nelements - 3, " more atoms)"
  end if
  print *

  ! Show detailed basis for first atom
  print *, "Detailed basis for atom 1 (", trim(mol_basis%elements(1)%element), "):"
  print *, "------------------------------------------"
  do ishell = 1, mol_basis%elements(1)%nshells
    print '(2X, "Shell ", I0, " (", A, "): ", I0, " primitives, ", I0, " basis functions")', &
      ishell, &
      ang_mom_int_to_char(mol_basis%elements(1)%shells(ishell)%ang_mom), &
      mol_basis%elements(1)%shells(ishell)%nfunc, &
      mol_basis%elements(1)%shells(ishell)%num_basis_functions()
  end do
  print *

  ! Cleanup
  call geom%destroy()
  call mol_basis%destroy()

  print *, "=========================================="
  print *, "Test PASSED!"
  print *, "=========================================="

end subroutine test_prism_with_basis

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

end program main
