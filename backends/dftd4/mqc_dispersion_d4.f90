!! Charge-dependent empirical dispersion, through dftd4's C API
module mqc_dispersion_d4
   !! -D4 on a Kohn-Sham energy and its nuclear gradient.
   !!
   !! The twin of `mqc_dispersion`, one directory over, and the differences from
   !! it are the whole reason this is a second wrapper rather than a second
   !! branch inside the first one.
   !!
   !! **D4 is charge dependent, and that is the difference that matters.**
   !! `dftd4_new_structure` takes a total molecular charge, and the model
   !! equilibrates atomic partial charges from it before any dispersion
   !! coefficient is interpolated. D3 ignores charge entirely, so its wrapper's
   !! signature has nowhere to carry one. A D4 energy computed at the wrong
   !! total charge is wrong by a few tenths of a millihartree on a small cation
   !! and looks perfectly plausible, which is why `charge` is a required
   !! argument here rather than an optional one with a zero default.
   !!
   !! **The three-body term is on.** `dftd4_load_rational_damping`'s `atm` flag
   !! is not "add a term to a fixed parametrisation": `.true.` selects
   !! `get_d4eeq_bjatm_parameter` and `.false.` selects `get_d4eeq_bj_parameter`,
   !! two separately fitted tables. `.true.` is what the library's own default
   !! is -- `get_rational_damping` sets `mbd = .true.` when no s9 is given, which
   !! is what the `dftd4` command line does -- and what its Python bindings do,
   !! `load_param(method, atm=True)`. It is also what "-D4" means in the
   !! literature: D4(EEQ)-ATM. So this is not D3's `.false.` carried over; the
   !! two libraries default the same flag differently and each is followed.
   !!
   !! **Units.** mqc is Bohr and Hartree throughout, and so is this API:
   !! `dftd4.h` says "quantities in Bohr" on the structure constructor, and
   !! `new_structure_api` declares `positions(3, natoms)` in the same storage
   !! order mqc uses. Nothing is scaled on the way in or out. The gradient comes
   !! back as dE/dR in Hartree/Bohr with the same (3, natoms) shape.
   !!
   !! **The C surface, not the Fortran one.** dftd4 is LGPL-3-or-later against
   !! this project's MIT, and is linked as a shared library so that the two stay
   !! separable. Using its .mod files would fix the compiler as well, so the
   !! interface below is declared here rather than imported. It mirrors
   !! `include/dftd4.h`; nothing in this build reads that header.
   !!
   !! Compiled only when MQC_ENABLE_DFTD4=ON. The twin that stands in otherwise
   !! is `src/methods/dft/mqc_dispersion_d4_stub.f90`, same module name, same two
   !! public procedures.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_bool, c_null_ptr, c_null_char, c_loc, &
                                                                             c_associated
   use pic_types, only: dp, default_int
   use mqc_error, only: error_t, ERROR_GENERIC, ERROR_VALIDATION
   use mqc_dispersion_names, only: d4_functional_alias, dispersion_kind_is_d4, DISPERSION_KINDS
   implicit none
   private

   public :: dispersion_d4_available
   public :: dispersion_d4_correction

   integer(c_int), parameter :: ERROR_BUFFER = 512_c_int

   logical(c_bool), parameter :: ATM = .true._c_bool
      !! The Axilrod-Teller-Muto parametrisation, which is what "-D4" names.
      !! See the module note: this selects a fitted table, not a term.

   ! The C API of dftd4, transcribed from include/dftd4.h.
   !
   ! Every handle is an opaque `type(c_ptr)`. The `delete_*` entry points take a
   ! pointer *to* the handle and null it, which is why those dummies are
   ! intent(inout) and not `value`.
   !
   ! `charge`, `lattice`, `periodic`, `gradient` and `sigma` are optional on the
   ! C side, which for a bind(C) dummy means "a null pointer is an absent
   ! argument". They are declared `type(c_ptr), value` here so that this side
   ! decides explicitly, by passing `c_null_ptr` or `c_loc(...)`, rather than
   ! relying on how a particular compiler passes an absent optional. `charge` is
   ! always passed.
   interface

      function dftd4_new_error() bind(c, name="dftd4_new_error") result(handle)
         import :: c_ptr
         implicit none
         type(c_ptr) :: handle
      end function dftd4_new_error

      function dftd4_check_error(handle) bind(c, name="dftd4_check_error") result(status)
         import :: c_ptr, c_int
         implicit none
         type(c_ptr), value :: handle
         integer(c_int) :: status
      end function dftd4_check_error

      subroutine dftd4_get_error(handle, buffer, buffersize) bind(c, name="dftd4_get_error")
         import :: c_ptr, c_char, c_int
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(in) :: buffersize
         character(kind=c_char), intent(inout) :: buffer(buffersize)
      end subroutine dftd4_get_error

      subroutine dftd4_delete_error(handle) bind(c, name="dftd4_delete_error")
         import :: c_ptr
         implicit none
         type(c_ptr), intent(inout) :: handle
      end subroutine dftd4_delete_error

      function dftd4_new_structure(handle, natoms, numbers, positions, charge, lattice, periodic) &
         bind(c, name="dftd4_new_structure") result(mol)
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: natoms
         integer(c_int), intent(in) :: numbers(natoms)
         real(c_double), intent(in) :: positions(3, natoms)
         type(c_ptr), value :: charge
            !! The total molecular charge, as a pointer to one double. Null
            !! means zero to the library; this wrapper never passes null.
         type(c_ptr), value :: lattice
         type(c_ptr), value :: periodic
         type(c_ptr) :: mol
      end function dftd4_new_structure

      subroutine dftd4_delete_structure(mol) bind(c, name="dftd4_delete_structure")
         import :: c_ptr
         implicit none
         type(c_ptr), intent(inout) :: mol
      end subroutine dftd4_delete_structure

      function dftd4_new_d4_model(handle, mol) bind(c, name="dftd4_new_d4_model") result(model)
         import :: c_ptr
         implicit none
         type(c_ptr), value :: handle
         type(c_ptr), value :: mol
         type(c_ptr) :: model
      end function dftd4_new_d4_model

      subroutine dftd4_delete_model(model) bind(c, name="dftd4_delete_model")
         import :: c_ptr
         implicit none
         type(c_ptr), intent(inout) :: model
      end subroutine dftd4_delete_model

      function dftd4_load_rational_damping(handle, method, atm) &
         bind(c, name="dftd4_load_rational_damping") result(param)
         import :: c_ptr, c_bool
         implicit none
         type(c_ptr), value :: handle
         type(c_ptr), value :: method
            !! A null-terminated C string. `type(c_ptr)` rather than an
            !! assumed-size `character(kind=c_char)` array so that nothing here
            !! declares an array whose extent it does not know.
         logical(c_bool), value :: atm
         type(c_ptr) :: param
      end function dftd4_load_rational_damping

      subroutine dftd4_delete_param(param) bind(c, name="dftd4_delete_param")
         import :: c_ptr
         implicit none
         type(c_ptr), intent(inout) :: param
      end subroutine dftd4_delete_param

      subroutine dftd4_get_dispersion(handle, mol, model, param, energy, gradient, sigma) &
         bind(c, name="dftd4_get_dispersion")
         import :: c_ptr, c_double
         implicit none
         type(c_ptr), value :: handle
         type(c_ptr), value :: mol
         type(c_ptr), value :: model
         type(c_ptr), value :: param
         real(c_double), intent(out) :: energy
         type(c_ptr), value :: gradient
         type(c_ptr), value :: sigma
      end subroutine dftd4_get_dispersion

   end interface

contains

   pure function dispersion_d4_available() result(available)
      !! .true. -- this build linked dftd4
      !!
      !! Read through `mqc_dispersion_apply`, which refuses a deck naming a
      !! correction this build cannot run, and by the version banner.
      logical :: available

      available = .true.
   end function dispersion_d4_available

   subroutine dispersion_d4_correction(kind, functional, charge, atomic_numbers, coordinates, &
                                       energy, gradient, error)
      !! The D4 dispersion energy, and its gradient when one is asked for
      !!
      !! `coordinates` is (3, natoms) in Bohr and `gradient` comes back in the
      !! same shape in Hartree/Bohr, both being what dftd4 already uses. The
      !! energy is negative for any geometry with two atoms in it.
      character(len=*), intent(in) :: kind
         !! Which correction: only "d4" reaches here. See `DISPERSION_KINDS`.
      character(len=*), intent(in) :: functional
         !! This program's spelling, from `model.functional`. Translated by
         !! `d4_functional_alias`, which refuses rather than guesses.
      real(dp), intent(in) :: charge
         !! The total molecular charge, in units of the elementary charge. Not
         !! optional: see the module note. A wrong one is a wrong energy that
         !! looks right.
      integer(default_int), intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)
      real(dp), intent(out) :: energy
      real(dp), intent(out), optional :: gradient(:, :)
      type(error_t), intent(out) :: error

      type(c_ptr) :: handle, mol, model, param
      integer(default_int) :: natoms
      integer(c_int), allocatable, target :: numbers(:)
      real(c_double), allocatable, target :: positions(:, :)
      real(c_double), allocatable, target :: forces(:, :)
      real(c_double), target :: virial(3, 3)
      real(c_double), target :: total_charge
      real(c_double) :: e_disp
      character(kind=c_char), allocatable, target :: method(:)
      character(len=32) :: alias
      logical :: want_gradient

      energy = 0.0_dp
      natoms = size(atomic_numbers, kind=default_int)
      want_gradient = present(gradient)
      if (want_gradient) gradient = 0.0_dp

      if (.not. dispersion_kind_is_d4(kind)) then
         call error%set(ERROR_VALIDATION, "unknown dispersion correction '"//trim(adjustl(kind))// &
                        "'. Known: "//DISPERSION_KINDS//".")
         return
      end if
      if (natoms <= 0) then
         call error%set(ERROR_VALIDATION, "a dispersion correction was asked for on a "// &
                        "fragment with no atoms")
         return
      end if
      if (size(coordinates, 2, kind=default_int) /= natoms) then
         call error%set(ERROR_GENERIC, "dispersion: coordinates and atomic numbers "// &
                        "disagree on how many atoms there are")
         return
      end if

      call d4_functional_alias(functional, alias, error)
      if (error%has_error()) return

      ! Fixed-kind copies: `c_int` and `c_double` need not be the default
      ! integer and `dp`, and a bind(C) dummy is not a place to find that out.
      ! Copies rather than views also mean `coordinates` stays intent(in) in
      ! fact as well as in declaration.
      allocate (numbers(natoms))
      numbers = int(atomic_numbers, c_int)
      allocate (positions(3, natoms))
      positions = real(coordinates, c_double)
      total_charge = real(charge, c_double)
      method = c_string(trim(alias))

      handle = dftd4_new_error()
      ! `c_loc(total_charge)` and never `c_null_ptr`. Null would be read as a
      ! neutral molecule, which is the one mistake this wrapper exists to make
      ! impossible: a cation's D4 energy computed at charge zero differs by
      ! enough to matter and by nothing that would show.
      mol = dftd4_new_structure(handle, int(natoms, c_int), numbers, positions, &
                                c_loc(total_charge), c_null_ptr, c_null_ptr)
      if (failed(handle, "building the structure", error)) then
         call dftd4_delete_error(handle)
         return
      end if

      ! The plain D4 model, not D4S. It is what `dftd4_load_rational_damping`'s
      ! parameters were fitted against and what the Python bindings build by
      ! default (`model="d4"`); D4S is a 4.0 addition with its own fits.
      model = dftd4_new_d4_model(handle, mol)
      if (failed(handle, "building the D4 model", error)) then
         call dftd4_delete_structure(mol)
         call dftd4_delete_error(handle)
         return
      end if

      param = dftd4_load_rational_damping(handle, c_loc(method), ATM)
      if (failed(handle, "loading D4 damping parameters for '"//trim(alias)//"'", error)) then
         call dftd4_delete_model(model)
         call dftd4_delete_structure(mol)
         call dftd4_delete_error(handle)
         return
      end if

      ! `get_dispersion_api` reads its sigma argument under the same `present`
      ! test as its gradient one, so the two travel together: either both are
      ! passed or neither is. The virial is computed and dropped -- mqc has no
      ! periodic path to spend it on.
      if (want_gradient) then
         allocate (forces(3, natoms))
         forces = 0.0_c_double
         virial = 0.0_c_double
         call dftd4_get_dispersion(handle, mol, model, param, e_disp, &
                                   c_loc(forces), c_loc(virial))
      else
         call dftd4_get_dispersion(handle, mol, model, param, e_disp, &
                                   c_null_ptr, c_null_ptr)
      end if

      if (.not. failed(handle, "evaluating the dispersion energy", error)) then
         energy = real(e_disp, dp)
         if (want_gradient) gradient = real(forces, dp)
      end if

      call dftd4_delete_param(param)
      call dftd4_delete_model(model)
      call dftd4_delete_structure(mol)
      call dftd4_delete_error(handle)
   end subroutine dispersion_d4_correction

   pure function c_string(text) result(buffer)
      !! A Fortran string as the null-terminated array of characters C expects
      character(len=*), intent(in) :: text
      character(kind=c_char), allocatable :: buffer(:)

      integer :: i

      allocate (buffer(len(text) + 1))
      do i = 1, len(text)
         buffer(i) = text(i:i)
      end do
      buffer(len(text) + 1) = c_null_char
   end function c_string

   function failed(handle, doing, error) result(bad)
      !! Whether dftd4 set its error handle, and if so what it said
      !!
      !! Every entry point is a no-op once the handle carries an error, so this
      !! is checked after each one rather than at the end: a message about
      !! damping parameters is a different bug report from one about a
      !! structure.
      type(c_ptr), intent(in) :: handle
      character(len=*), intent(in) :: doing
      type(error_t), intent(inout) :: error
      logical :: bad

      character(kind=c_char) :: buffer(ERROR_BUFFER)
      character(len=ERROR_BUFFER) :: message
      integer :: i

      bad = .false.
      if (.not. c_associated(handle)) then
         call error%set(ERROR_GENERIC, "dftd4 would not create an error handle")
         bad = .true.
         return
      end if
      if (dftd4_check_error(handle) == 0_c_int) return

      buffer = c_null_char
      call dftd4_get_error(handle, buffer, ERROR_BUFFER)
      message = ""
      do i = 1, int(ERROR_BUFFER)
         if (buffer(i) == c_null_char) exit
         message(i:i) = buffer(i)
      end do

      call error%set(ERROR_VALIDATION, "dftd4 refused "//trim(doing)//": "//trim(message))
      bad = .true.
   end function failed

end module mqc_dispersion_d4
