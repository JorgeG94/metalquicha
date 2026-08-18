!! The accurate atomic minimal basis, and the orbital-space dimensions it fixes
module mqc_aambs
   !! The quasi-atomic orbital analysis projects molecular orbitals onto free-atom
   !! minimal-basis orbitals. This module owns that basis and the counting that
   !! comes with it: how many minimal-basis orbitals a molecule has, how many of
   !! them are chemical core, and therefore how many valence-virtual orbitals have
   !! to be recovered from the virtual space.
   !!
   !! The counting is not a detail. Every space in the construction is defined by
   !! its dimension rather than by a threshold -- `N_VVO = N_mbs - N_occ` exactly,
   !! not "however many singular values look large" -- so a wrong count is not a
   !! small error in the answer, it is a different calculation. That is also what
   !! makes it testable before any integral exists: Paper I tabulates these
   !! numbers for eight molecules and they can be checked against nothing but the
   !! basis file.
   !!
   !! `basis_sets/aambs/aambs.json` is Basis Set Exchange schema plus two fields
   !! per element that BSE has no slot for and this module needs:
   !! `n_core_orbitals` and `n_valence_orbitals`. The core is the *chemical* core,
   !! which counts semicore as core -- so gallium has 14, not 10, because its 3d
   !! shell is chemically inert even though a frozen-core count would leave it
   !! out. See `basis_sets/aambs/PROVENANCE.md`.
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

   integer, parameter :: MAX_Z = 54
      !! The non-relativistic table covers hydrogen to xenon. Beyond that GAMESS
      !! switches to a relativistic set which is not shipped here, and which needs
      !! a relativistic Hamiltonian this code does not have.

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

   subroutine aambs_file(filename, error)
      !! Locate `aambs/aambs.json` on the basis search path
      !!
      !! In a subdirectory rather than beside the other basis sets, and
      !! deliberately: `mqc_extract_basis_sets` deletes any `basis_sets/*.json`
      !! that `MQC_BASIS_SETS` does not name, and this file is tracked rather
      !! than extracted from the Basis Set Exchange bundle. The subdirectory is
      !! what keeps a configure from removing it.
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
      !! `n_vvo` can come out negative, and that is a real condition rather than
      !! an arithmetic accident: it means the molecule has more occupied orbitals
      !! than its minimal basis has room for. Paper I's construction has nothing
      !! to say about that case, so it is refused here rather than propagated as
      !! a negative dimension into an allocate.
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
