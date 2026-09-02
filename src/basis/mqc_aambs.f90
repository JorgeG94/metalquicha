!! The accurate atomic minimal basis, and the orbital-space dimensions it fixes
module mqc_aambs
   !! The quasi-atomic orbital analysis projects molecular orbitals onto free-atom
   !! minimal-basis orbitals. This module owns that basis and the counting that
   !! comes with it: how many minimal-basis orbitals a molecule has, how many of
   !! them are chemical core, and therefore how many valence-virtual orbitals have
   !! to be recovered from the virtual space.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_IO, ERROR_VALIDATION
   use mqc_basis_utils, only: basis_search_path
   use json_module, only: json_file
   use pic_io, only: to_char
   implicit none
   private

   public :: aambs_file            !! Locate aambs.json
   public :: aambs_element_counts  !! Chemical core and valence orbital counts for one Z
   public :: aambs_dimensions      !! The minimal-basis dimensions of a whole molecule
   public :: aambs_dimensions_t
   public :: aambs_shell_labels  !! Principal and angular numbers, per minimal-basis function

   integer, parameter :: MAX_Z = 54
      !! The non-relativistic table covers hydrogen to xenon. The relativistic
      !! set beyond it is not shipped.

   type :: aambs_dimensions_t
      !! Every orbital-space dimension the QUAO construction needs
      integer :: n_mbs = 0        !! Minimal-basis orbitals, summed over atoms
      integer :: n_core = 0       !! Chemical-core orbitals, summed over atoms
      integer :: n_valence = 0    !! `n_mbs - n_core`
      integer :: n_occupied = 0   !! Doubly occupied molecular orbitals
      integer :: n_valocc = 0     !! Occupied that are not core, `n_occupied - n_core`
      integer :: n_vvo = 0        !! Valence-virtual orbitals to recover, `n_mbs - n_occupied`
   end type aambs_dimensions_t

contains

   subroutine aambs_shell_labels(atomic_numbers, atom_of, principal, angular, error)
      !! Which atom, shell and angular momentum each minimal-basis function is
      !!
      !! In the order `build_aambs_molecule` lays them out -- atoms outermost,
      !! then the element's shells in the order `aambs.json` lists them, then
      !! `2l+1` functions per shell. That ordering is an invariant of the basis
      !! file, not of this routine.
      !!
      !! The principal quantum number is in the file and in nothing else:
      !! libcint has no notion of a `2p` as opposed to a `3p`.
      integer, intent(in) :: atomic_numbers(:)
      integer, allocatable, intent(out) :: atom_of(:)     !! 1-based atom index
      integer, allocatable, intent(out) :: principal(:)   !! n
      integer, allocatable, intent(out) :: angular(:)     !! l
      type(error_t), intent(inout) :: error

      type(json_file) :: json
      character(len=:), allocatable :: path, key
      integer, allocatable :: shell_n(:), shell_l(:)
      integer :: iatom, ishell, icomp, total, cursor, n_shells
      logical :: found
      ! TODO(mqc): `key`, `n_shells` and `found` are declared here and never
      ! used.

      if (error%has_error()) return
      call aambs_file(path, error)
      if (error%has_error()) return
      call json%initialize()
      call json%load_file(filename=path)
      if (json%failed()) then
         call error%set(ERROR_IO, "could not read the atomic minimal basis at "//path)
         call json%destroy()
         return
      end if

      ! Two passes: the first counts so the arrays can be sized exactly, the
      ! second fills them.
      total = 0
      do iatom = 1, size(atomic_numbers)
         call element_shells(json, atomic_numbers(iatom), shell_n, shell_l, error)
         if (error%has_error()) then
            call json%destroy()
            return
         end if
         do ishell = 1, size(shell_l)
            total = total + 2*shell_l(ishell) + 1
         end do
      end do

      allocate (atom_of(total), principal(total), angular(total))
      cursor = 0
      do iatom = 1, size(atomic_numbers)
         call element_shells(json, atomic_numbers(iatom), shell_n, shell_l, error)
         if (error%has_error()) then
            call json%destroy()
            return
         end if
         do ishell = 1, size(shell_l)
            do icomp = 1, 2*shell_l(ishell) + 1
               cursor = cursor + 1
               atom_of(cursor) = iatom
               principal(cursor) = shell_n(ishell)
               angular(cursor) = shell_l(ishell)
            end do
         end do
      end do

      call json%destroy()
   end subroutine aambs_shell_labels

   subroutine element_shells(json, atomic_number, shell_n, shell_l, error)
      !! The principal and angular numbers of one element's shells, in file order
      type(json_file), intent(inout) :: json
      integer, intent(in) :: atomic_number
      integer, allocatable, intent(out) :: shell_n(:), shell_l(:)
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: key
      integer, allocatable :: momenta(:)
      integer :: ishell, value
      logical :: found

      if (error%has_error()) return
      if (allocated(shell_n)) deallocate (shell_n)
      if (allocated(shell_l)) deallocate (shell_l)
      allocate (shell_n(0), shell_l(0))

      ishell = 0
      do
         ishell = ishell + 1
         key = "elements."//to_char(atomic_number)//".electron_shells("// &
               to_char(ishell)//")"
         call json%get(key//".n", value, found)
         if (.not. found) exit
         shell_n = [shell_n, value]
         call json%get(key//".angular_momentum", momenta, found)
         if (.not. found .or. size(momenta) /= 1) then
            call error%set(ERROR_VALIDATION, "shell "//to_char(ishell)//" of Z = "// &
                           to_char(atomic_number)//" in the atomic minimal basis "// &
                           "does not carry exactly one angular momentum. Every "// &
                           "free-atom shell should.")
            return
         end if
         shell_l = [shell_l, momenta(1)]
      end do

      if (size(shell_l) == 0) then
         call error%set(ERROR_VALIDATION, "Z = "//to_char(atomic_number)// &
                        " has no shells in the atomic minimal basis.")
      end if
   end subroutine element_shells

   subroutine aambs_file(filename, error)
      !! Locate `aambs/aambs.json` on the basis search path
      ! In a subdirectory rather than beside the other basis sets:
      ! `mqc_extract_basis_sets` deletes any `*.json` under `basis_sets/` that
      ! `MQC_BASIS_SETS` does not name, and this file is tracked rather than
      ! extracted. Written `*.json` under `basis_sets/` rather than the other
      ! way round because CMake's Ninja generator preprocesses every Fortran
      ! source to scan module dependencies, and `cpp` reads a literal
      ! slash-star as the start of a C comment that never closes.
      character(len=:), allocatable, intent(out) :: filename
      type(error_t), intent(out) :: error

      character(len=:), allocatable :: directories(:)
      character(len=:), allocatable :: candidate, tried
      logical :: exists
      integer :: idir

      directories = basis_search_path()
      tried = ""
      do idir = 1, size(directories)
         candidate = trim(directories(idir))//"/aambs/aambs.json"
         inquire (file=candidate, exist=exists)
         if (exists) then
            filename = candidate
            return
         end if
         if (len(tried) > 0) tried = tried//", "
         tried = tried//candidate
      end do

      call error%set(ERROR_IO, "the accurate atomic minimal basis was not found. "// &
                     "Tried "//tried//". It ships as basis_sets/aambs/aambs.json and "// &
                     "is tracked rather than extracted, so a missing file means the "// &
                     "checkout is incomplete rather than that a basis was not requested.")
   end subroutine aambs_file

   subroutine aambs_element_counts(atomic_number, n_core, n_valence, error)
      !! Chemical-core and valence orbital counts for one element
      integer, intent(in) :: atomic_number
      integer, intent(out) :: n_core
      integer, intent(out) :: n_valence
      type(error_t), intent(inout) :: error

      type(json_file) :: json
      character(len=:), allocatable :: path, key
      logical :: found_core, found_valence

      n_core = 0
      n_valence = 0
      if (error%has_error()) return

      if (atomic_number < 1 .or. atomic_number > MAX_Z) then
         call error%set(ERROR_VALIDATION, "the accurate atomic minimal basis covers "// &
                        "hydrogen to xenon (Z = 1 to "//to_char(MAX_Z)//"); Z = "// &
                        to_char(atomic_number)//" is outside it. GAMESS switches to a "// &
                        "relativistic table here, which is not shipped -- refused rather "// &
                        "than run against a basis that does not describe the atom.")
         return
      end if

      call aambs_file(path, error)
      if (error%has_error()) return

      call json%initialize()
      call json%load_file(filename=path)
      if (json%failed()) then
         call error%set(ERROR_IO, "could not read the atomic minimal basis at "//path)
         call json%destroy()
         return
      end if

      key = "elements."//to_char(atomic_number)
      call json%get(key//".n_core_orbitals", n_core, found_core)
      call json%get(key//".n_valence_orbitals", n_valence, found_valence)
      call json%destroy()

      if (.not. (found_core .and. found_valence)) then
         call error%set(ERROR_IO, "the atomic minimal basis has no orbital counts for "// &
                        "Z = "//to_char(atomic_number)//". These are not Basis Set "// &
                        "Exchange fields; they are added by the extraction and the "// &
                        "construction cannot proceed without them.")
         n_core = 0
         n_valence = 0
      end if
   end subroutine aambs_element_counts

   subroutine aambs_dimensions(atomic_numbers, n_electrons, dims, error)
      !! The orbital-space dimensions of a molecule in the minimal basis
      !!
      !! A molecule with more occupied orbitals than its minimal basis has room
      !! for gives a negative `n_vvo`; that is refused here rather than
      !! propagated into an allocate.
      integer, intent(in) :: atomic_numbers(:)
      integer, intent(in) :: n_electrons
      type(aambs_dimensions_t), intent(out) :: dims
      type(error_t), intent(inout) :: error

      integer :: iatom, core, valence

      if (error%has_error()) return

      if (mod(n_electrons, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "the quasi-atomic analysis is closed shell: "// &
                        to_char(n_electrons)//" electrons is odd. Paper I derives the "// &
                        "construction for a single-determinant restricted reference.")
         return
      end if

      do iatom = 1, size(atomic_numbers)
         call aambs_element_counts(atomic_numbers(iatom), core, valence, error)
         if (error%has_error()) return
         dims%n_core = dims%n_core + core
         dims%n_mbs = dims%n_mbs + core + valence
      end do

      dims%n_valence = dims%n_mbs - dims%n_core
      dims%n_occupied = n_electrons/2
      dims%n_valocc = dims%n_occupied - dims%n_core
      dims%n_vvo = dims%n_mbs - dims%n_occupied

      if (dims%n_valocc < 0) then
         call error%set(ERROR_VALIDATION, "this molecule has fewer occupied orbitals ("// &
                        to_char(dims%n_occupied)//") than chemical core orbitals ("// &
                        to_char(dims%n_core)//"), which cannot happen for a neutral or "// &
                        "modestly charged species and means the electron count and the "// &
                        "atoms disagree.")
         return
      end if

      if (dims%n_vvo < 0) then
         call error%set(ERROR_VALIDATION, "the minimal basis has "//to_char(dims%n_mbs)// &
                        " orbitals but the molecule has "//to_char(dims%n_occupied)// &
                        " occupied ones, so there are no valence-virtual orbitals to "// &
                        "recover. The construction assumes the occupied space fits "// &
                        "inside the minimal basis.")
         return
      end if
   end subroutine aambs_dimensions

end module mqc_aambs
