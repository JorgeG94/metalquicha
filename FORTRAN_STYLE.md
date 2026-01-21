# Metalquicha Fortran Style Guide

This document defines the coding conventions for all Fortran code within the metalquicha project. External dependencies (pic, tblite, etc.) may follow their own conventions.

## Naming Conventions

### Files
- Module files: `mqc_<name>.f90` (lowercase with underscores)
- Preprocessed files: `mqc_<name>.F90` (uppercase extension)
- Test files: `test_mqc_<name>.f90`

### Modules
- Module names match filenames: `mqc_<name>`
- One module per file (except for submodules)

```fortran
! Good
module mqc_physical_fragment

! Bad
module PhysicalFragment
module physical_fragment  ! missing mqc_ prefix
```

### Derived Types
- All derived types use `_t` suffix
- Use lowercase with underscores

```fortran
! Good
type :: calculation_result_t
type :: physical_fragment_t
type :: energy_t

! Bad
type :: CalculationResult
type :: calculation_result   ! missing _t suffix
type :: TFragment
```

### Variables and Procedures
- Use lowercase with underscores (snake_case)
- Use descriptive names - avoid single letters except for loop indices or if it is a very local, easy to deduce variable

```fortran
! Good
integer :: n_atoms
real(dp) :: total_energy
subroutine calculate_gradient(fragment, result)

! Bad
integer :: nAtoms      ! camelCase
integer :: na          ! too short
real(dp) :: E          ! too short, unclear
```

### Constants
- Use UPPERCASE with underscores for parameters

```fortran
! Good
integer, parameter :: MAX_ITERATIONS = 100
real(dp), parameter :: BOHR_TO_ANGSTROM = 0.529177210903_dp

! Bad
integer, parameter :: maxIterations = 100
```

## Required Practices

### Use Statements
- **Always use `only` clause** - enforced by `fortitude` linter

```fortran
! Good
use pic_types, only: dp, int32
use mqc_result_types, only: calculation_result_t, energy_t

! Bad - will fail linter
use pic_types
use mqc_result_types
```

### Implicit None
- **Always include `implicit none`** in modules and programs

```fortran
module mqc_example
   use pic_types, only: dp
   implicit none
   private
   ! ...
end module mqc_example
```

### Intent Declarations
- **Always declare intent** for all procedure arguments

```fortran
! Good
subroutine compute_energy(fragment, result)
   type(physical_fragment_t), intent(in) :: fragment
   type(calculation_result_t), intent(out) :: result

! Bad
subroutine compute_energy(fragment, result)
   type(physical_fragment_t) :: fragment  ! missing intent
```

### Private by Default (linter enforced)
- Modules should be `private` by default
- Explicitly declare `public` entities

```fortran
module mqc_example
   implicit none
   private

   public :: my_public_type_t
   public :: my_public_subroutine

   type :: my_public_type_t
      ! ...
   end type

   type :: internal_helper_t  ! stays private
      ! ...
   end type
end module
```

### Limit Procedure Arguments
- Public subroutines/functions should have **6 or fewer arguments**
- If more are needed, group related arguments into a derived type
- This improves readability, maintainability, and makes API changes easier

```fortran
! Bad - too many arguments
subroutine run_calculation(coords, elements, charge, multiplicity, &
                           method, basis, max_iter, tolerance, &
                           energy, gradient, error)
   real(dp), intent(in) :: coords(:,:)
   integer, intent(in) :: elements(:)
   integer, intent(in) :: charge, multiplicity
   character(*), intent(in) :: method, basis
   integer, intent(in) :: max_iter
   real(dp), intent(in) :: tolerance
   real(dp), intent(out) :: energy
   real(dp), intent(out) :: gradient(:,:)
   type(error_t), intent(out) :: error
   ! 11 arguments - hard to use and maintain

! Good - grouped into logical types
subroutine run_calculation(system, config, result, error)
   type(system_t), intent(in) :: system           ! coords, elements, charge, mult
   type(calc_config_t), intent(in) :: config      ! method, basis, max_iter, tol
   type(calculation_result_t), intent(out) :: result  ! energy, gradient
   type(error_t), intent(out), optional :: error
   ! 4 arguments - clean and extensible
```

- **Exception**: Simple utility functions (e.g., `to_bohr(value)`) can have minimal arguments
- Internal/private procedures have more flexibility but should still aim for clarity

## Forbidden Practices

### No GOTO Statements
- Use structured control flow (`do`, `if`, `select case`, `exit`, `cycle`)

```fortran
! Bad
   if (error) goto 999
   ! ...
999 continue
   print *, "Error occurred"

! Good
   if (error) then
      print *, "Error occurred"
      return
   end if
```

### No Arithmetic IF
- Use `if-then-else` or `select case`

```fortran
! Bad (arithmetic IF - Fortran 77)
   if (x) 10, 20, 30

! Good
   if (x < 0) then
      ! negative case
   else if (x == 0) then
      ! zero case
   else
      ! positive case
   end if
```

### No Implicit SAVE
- Avoid module-level variables that retain state
- If state is needed, use explicit `save` attribute and document why

```fortran
! Bad - implicit save behavior
module mqc_bad_example
   integer :: counter = 0  ! implicitly saved, retains value between calls
end module

! Good - use derived types to manage state
module mqc_good_example
   type :: counter_t
      integer :: value = 0
   end type
end module
```

### No COMMON Blocks
- Use module variables or derived types instead

```fortran
! Bad
common /shared_data/ x, y, z

! Good
module mqc_shared_data
   use pic_types, only: dp
   implicit none
   type :: shared_data_t
      real(dp) :: x, y, z
   end type
end module
```

### No Naked Print Statements
- **NEVER** use `print *` or `write(*,*)` for output
- **ALWAYS** use `pic_logger` (`global_logger`)
- This is enforced - naked prints will be rejected in code review

```fortran
! Bad - forbidden
print *, "Starting calculation"
write(*,*) "Energy:", energy
write(6,*) "Done"

! Good - use logger
use pic_logger, only: logger => global_logger
call logger%info("Starting calculation")
call logger%info("Energy: " // to_char(energy))
call logger%info("Done")
```

### No Emojis in Fortran Code
- No emojis in Fortran source files (`.f90`, `.F90`)
- This includes comments, strings, and documentation
- Keep output professional and portable

```fortran
! Bad
call logger%info("Calculation complete! 🎉")
!! 🚀 Fast implementation of MBE

! Good
call logger%info("Calculation complete")
!! Fast implementation of MBE
```

- Python scripts (validation, tooling) may use emojis if desired

### No EQUIVALENCE
- Use proper type conversions or `transfer()` if absolutely necessary

### No Fixed-Form Source
- All code must be free-form (`.f90` / `.F90`)
- No `.f` or `.F` files

### No Assumed-Size Arrays
- Use assumed-shape arrays with explicit interface

```fortran
! Bad
subroutine process_array(arr, n)
   real(dp) :: arr(*)  ! assumed-size
   integer :: n

! Good
subroutine process_array(arr)
   real(dp), intent(in) :: arr(:)  ! assumed-shape
```

### No External Statements for Internal Procedures
- Use `contains` for internal procedures
- Use explicit interfaces via modules

```fortran
! Bad
external :: my_function

! Good - use modules
use mqc_my_module, only: my_function
```

## Recommended Practices

### Use pic Library Utilities
- **Always prefer pic functionality** over implementing your own
- Ensures consistency, reduces code duplication, and leverages tested implementations

| pic Module | Functionality | Use Instead Of |
|------------|---------------|----------------|
| `pic_types` | Kind parameters (`dp`, `int32`, etc.) | Literal kinds (`real(8)`) |
| `pic_logger` | Logging (`global_logger`) | `print *` / `write(*,*)` (**forbidden**) |
| `pic_timer` | Performance timing | Manual `cpu_time` calls |
| `pic_strings` | String utilities | Custom string manipulation |
| `pic_sorting` | Sorting algorithms | Hand-rolled sorts |
| `pic_math` | Math utilities | Reimplementing common math |
| `pic_test_helpers` | Test utilities (`is_equal`) | Custom comparison functions |

```fortran
! Good - use pic utilities
use pic_logger, only: logger => global_logger
use pic_timer, only: timer_t
use pic_sorting, only: argsort

call logger%info("Starting calculation")
call logger%debug("Processing fragment " // to_char(i))
call logger%warning("Large system detected, may be slow")
call logger%error("Invalid input: negative charge")

type(timer_t) :: timer
call timer%start()
! ... work ...
call timer%stop()
call logger%info("Elapsed: " // timer%elapsed_string())

indices = argsort(energies)  ! sorted indices

! Bad - NEVER use naked print statements
print *, "Starting calculation"  ! no log levels, no MPI rank awareness
write(*,*) "Debug info"          ! can't be silenced, clutters output
```

**Why logger over print?**
- Log levels (debug/info/warning/error) - control verbosity
- Consistent formatting across the codebase
- Can be redirected to files or silenced entirely
- Easier to add MPI-awareness later (TODO: only rank 0 prints)

- Check pic documentation before implementing common functionality
- If pic is missing something useful, consider contributing it upstream

### No Naked MPI or BLAS Calls

#### pic-mpi (Critical for Portability)
- **ALWAYS use pic_mpi_lib** for MPI operations - **NEVER** call MPI directly
- This is **critical**: pic-mpi abstracts the MPI backend, allowing seamless switching between `mpi_f08` and legacy `mpi` modules
- Some compilers/systems only support one or the other - pic-mpi handles this automatically
- Direct MPI calls break this portability and will cause build failures on some systems

```fortran
! Good - use pic-mpi wrappers
use pic_mpi_lib, only: comm_t, send, recv, bcast, allgather, &
                       isend, irecv, wait, iprobe, request_t, &
                       MPI_Status, MPI_ANY_SOURCE, MPI_ANY_TAG, abort_comm

type(comm_t) :: comm
call bcast(comm, data, root=0)
call send(comm, data, dest=1, tag=TAG_DATA)
call recv(comm, buffer, source=0, tag=TAG_DATA)

! Bad - naked MPI calls
use mpi_f08
call MPI_Bcast(data, size(data), MPI_DOUBLE_PRECISION, 0, comm, ierr)
call MPI_Send(data, size(data), MPI_DOUBLE_PRECISION, 1, tag, comm, ierr)
```

#### pic-blas
- **Always use pic-blas** for BLAS/LAPACK - never call BLAS/LAPACK directly
- Provides cleaner interface and handles different BLAS implementations

```fortran
! Good - use pic-blas wrappers
use pic_blas, only: pic_gemm, pic_dot

call pic_gemm(A, B, C)
result = pic_dot(x, y)

! Bad - naked BLAS calls
call dgemm('N', 'N', m, n, k, 1.0d0, A, lda, B, ldb, 0.0d0, C, ldc)
result = ddot(n, x, 1, y, 1)
```

- Exceptions: Only if pic-mpi/pic-blas lacks needed functionality (then consider contributing)

### Kind Parameters
- Use `pic_types` for portable kind definitions
- Never use literal kind numbers

```fortran
! Good
use pic_types, only: dp, int32
real(dp) :: energy
integer(int32) :: count

! Bad
real(8) :: energy      ! non-portable
real*8 :: energy       ! obsolete syntax
double precision :: e  ! obsolete
```

### Array Operations
- Prefer whole-array operations over explicit loops when clear

```fortran
! Good
gradient = 0.0_dp
total = sum(energies)

! Also fine for complex operations
do i = 1, n_atoms
   gradient(:, i) = gradient(:, i) + contribution(:, i)
end do
```

### Block Constructs for Limiting Scope
- Use `block` to limit variable scope and improve readability
- Useful for temporary variables needed only in a small section
- Helps prevent accidental reuse of variables

```fortran
! Good - temporary variables scoped to where they're needed
subroutine process_fragments(fragments, total_energy)
   type(fragment_t), intent(in) :: fragments(:)
   real(dp), intent(out) :: total_energy

   integer :: i

   total_energy = 0.0_dp
   do i = 1, size(fragments)
      block
         real(dp) :: frag_energy
         real(dp) :: correction

         call compute_energy(fragments(i), frag_energy)
         call compute_correction(fragments(i), correction)
         total_energy = total_energy + frag_energy + correction
      end block
   end do
end subroutine

! Bad - temporaries pollute the entire procedure scope
subroutine process_fragments(fragments, total_energy)
   type(fragment_t), intent(in) :: fragments(:)
   real(dp), intent(out) :: total_energy

   integer :: i
   real(dp) :: frag_energy    ! visible everywhere, might be misused later
   real(dp) :: correction     ! visible everywhere

   total_energy = 0.0_dp
   do i = 1, size(fragments)
      call compute_energy(fragments(i), frag_energy)
      call compute_correction(fragments(i), correction)
      total_energy = total_energy + frag_energy + correction
   end do
end subroutine
```

### Associate Construct for Readability
- Use `associate` to create short aliases for long expressions
- Improves readability without runtime cost

```fortran
! Good
associate(coords => fragment%coordinates, &
          n => fragment%n_atoms)
   do i = 1, n
      distance = norm2(coords(:, i) - origin)
   end do
end associate

! Bad - repeated long expressions
do i = 1, fragment%n_atoms
   distance = norm2(fragment%coordinates(:, i) - origin)
end do
```

### Error Handling
- Use the `error_t` type from `mqc_error` module
- Check and propagate errors

```fortran
use mqc_error, only: error_t, create_error, has_error

subroutine my_subroutine(input, output, error)
   type(input_t), intent(in) :: input
   type(output_t), intent(out) :: output
   type(error_t), intent(out), optional :: error

   if (invalid_input) then
      if (present(error)) then
         call create_error(error, "Invalid input: reason")
      end if
      return
   end if
end subroutine
```

### Documentation
- Use `!!` for FORD documentation comments
- Document public interfaces

```fortran
type :: calculation_result_t
   !! Container for quantum chemistry calculation results
   !!
   !! Holds energy, gradient, Hessian, and associated metadata
   !! from a single-point or property calculation.

   type(energy_t) :: energy
      !! Total and component energies (SCF, correlation, etc.)
   real(dp), allocatable :: gradient(:,:)
      !! Nuclear gradient (3, n_atoms) in Hartree/Bohr
end type
```

### Memory Management
- Provide `destroy` procedures for types with allocatable components
- Clean up allocatable arrays when no longer needed

```fortran
type :: my_type_t
   real(dp), allocatable :: data(:)
contains
   procedure :: destroy => my_type_destroy
end type

subroutine my_type_destroy(this)
   class(my_type_t), intent(inout) :: this
   if (allocated(this%data)) deallocate(this%data)
end subroutine
```

### Prefer Allocatable Over Pointer
- Use `allocatable` instead of `pointer` when possible
- Allocatable arrays are automatically deallocated when out of scope
- Less risk of memory leaks and dangling pointers
- Compiler can optimize allocatable better

```fortran
! Good - automatic cleanup, no leak risk
type :: data_container_t
   real(dp), allocatable :: values(:)
   real(dp), allocatable :: matrix(:,:)
end type

! Bad - manual cleanup required, leak risk
type :: data_container_t
   real(dp), pointer :: values(:) => null()
   real(dp), pointer :: matrix(:,:) => null()
end type
```

- **When to use pointer**: Only when you need pointer semantics (aliasing, linked structures, polymorphic returns)

### Pure and Elemental Procedures
- Mark functions/subroutines as `pure` when they have no side effects
- Use `elemental` for scalar operations that should work on arrays
- Enables compiler optimizations and documents intent

```fortran
! Good - pure function, no side effects
pure function kinetic_energy(mass, velocity) result(energy)
   real(dp), intent(in) :: mass, velocity
   real(dp) :: energy
   energy = 0.5_dp * mass * velocity**2
end function

! Good - elemental, works on scalars and arrays
elemental function to_bohr(angstrom) result(bohr)
   real(dp), intent(in) :: angstrom
   real(dp) :: bohr
   bohr = angstrom * ANGSTROM_TO_BOHR
end function

! Usage: works on scalar or array
r_bohr = to_bohr(1.5_dp)           ! scalar
coords_bohr = to_bohr(coords_ang)   ! array
```

### Avoid Deep Nesting
- Maximum 3-4 levels of indentation
- Use early returns, `cycle`, `exit` to reduce nesting
- Extract deeply nested code to separate subroutines

```fortran
! Bad - deeply nested, hard to follow
do i = 1, n
   if (condition1) then
      if (condition2) then
         do j = 1, m
            if (condition3) then
               ! actual work buried here
            end if
         end do
      end if
   end if
end do

! Good - early cycle, flat structure
do i = 1, n
   if (.not. condition1) cycle
   if (.not. condition2) cycle

   do j = 1, m
      if (.not. condition3) cycle
      ! actual work at reasonable depth
   end do
end do

! Also good - extract to subroutine
do i = 1, n
   call process_item(items(i), result)
end do
```

### No Magic Numbers
- Use named constants for any non-obvious literal values
- Makes code self-documenting and easier to maintain

```fortran
! Bad - what do these numbers mean?
if (n_atoms > 50) then
   cutoff = 4.5
end if
tolerance = 1.0e-8

! Good - self-documenting
integer, parameter :: LARGE_SYSTEM_THRESHOLD = 50
real(dp), parameter :: DEFAULT_DISTANCE_CUTOFF = 4.5_dp  ! Angstrom
real(dp), parameter :: SCF_CONVERGENCE_TOLERANCE = 1.0e-8_dp

if (n_atoms > LARGE_SYSTEM_THRESHOLD) then
   cutoff = DEFAULT_DISTANCE_CUTOFF
end if
tolerance = SCF_CONVERGENCE_TOLERANCE
```

### Allocatable Character Strings
- Use `character(len=:), allocatable` for dynamic strings
- Avoid fixed-length strings that waste memory or truncate

```fortran
! Good - allocatable string
character(len=:), allocatable :: method_name
character(len=:), allocatable :: error_message

method_name = "GFN2-xTB"  ! automatically sized
error_message = "Error in " // trim(filename) // ": " // trim(reason)

! Bad - fixed length, may truncate or waste space
character(len=32) :: method_name   ! what if name is longer?
character(len=256) :: error_message  ! wastes memory for short messages
```

### Do Concurrent (Use with Caution)
- `do concurrent` hints to compiler that iterations are independent
- Can enable auto-parallelization and vectorization
- **Caution**: Compiler support varies; behavior may differ between compilers
- Test thoroughly when using; consider guarding with preprocessor if needed

```fortran
! Good candidate - independent iterations, no dependencies
do concurrent (i = 1:n)
   result(i) = input(i) ** 2
end do

! Bad candidate - has dependencies
do concurrent (i = 2:n)
   result(i) = result(i-1) + input(i)  ! depends on previous iteration!
end do

! Safe alternative if do concurrent causes issues
do i = 1, n
   result(i) = input(i) ** 2
end do
```

- **Restrictions**: No `exit`, `cycle`, `return`, or `goto` inside `do concurrent`
- **When in doubt**: Use regular `do` loop - correctness over optimization

### Submodules for Large Modules
- Use submodules to separate interface from implementation
- Reduces recompilation when only implementation changes
- Useful for large modules with many procedures

```fortran
! In mqc_heavy_module.f90 - interface only
module mqc_heavy_module
   implicit none
   private
   public :: compute_expensive_thing

   interface
      module subroutine compute_expensive_thing(input, output)
         real(dp), intent(in) :: input(:)
         real(dp), intent(out) :: output(:)
      end subroutine
   end interface
end module

! In mqc_heavy_module_impl.f90 - implementation
submodule (mqc_heavy_module) mqc_heavy_module_impl
contains
   module subroutine compute_expensive_thing(input, output)
      real(dp), intent(in) :: input(:)
      real(dp), intent(out) :: output(:)
      ! Long implementation here
      ! Changes here don't trigger recompilation of modules that use mqc_heavy_module
   end subroutine
end submodule
```

## Units

- **Internal units**: Bohr (length), Hartree (energy)
- **Conversion**: Use `to_bohr()` / `to_angstrom()` from `mqc_physical_fragment`
- **Document units** in comments when not obvious

```fortran
real(dp) :: bond_length      ! in Bohr
real(dp) :: energy           ! in Hartree
real(dp) :: frequency_cm1    ! in cm^-1 (document non-atomic units)
```

## File Structure Template

```fortran
!! Brief module description
module mqc_example
   !! Extended module documentation
   !!
   !! More details about the module purpose and usage.
   use pic_types, only: dp, int32
   use mqc_other_module, only: needed_type_t
   implicit none
   private

   public :: example_type_t
   public :: example_subroutine

   type :: example_type_t
      !! Type documentation
      integer :: n_items
         !! Number of items
      real(dp), allocatable :: values(:)
         !! Array of values in Hartree
   contains
      procedure :: compute => example_compute
      procedure :: destroy => example_destroy
   end type example_type_t

contains

   subroutine example_compute(this, input, output)
      !! Subroutine documentation
      class(example_type_t), intent(inout) :: this
      real(dp), intent(in) :: input
      real(dp), intent(out) :: output

      ! Implementation
   end subroutine example_compute

   subroutine example_destroy(this)
      !! Clean up allocated memory
      class(example_type_t), intent(inout) :: this

      if (allocated(this%values)) deallocate(this%values)
   end subroutine example_destroy

end module mqc_example
```

## Enforcement

- **Linter**: Run `fortitude check` before committing
- **CI**: GitHub Actions will fail if linter errors exist
- **Code Review**: Verify conventions are followed in PRs
