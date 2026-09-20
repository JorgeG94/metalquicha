!! Empirical dispersion, through s-dftd3's C API
module mqc_dispersion
   !! -D3(BJ) on a Kohn-Sham energy and its nuclear gradient.
   !!
   !! The correction itself is not computed here and deliberately is not: the
   !! reference C6 coefficients D3 interpolates are sixty-six thousand lines of
   !! tabulated data that exist only inside existing distributions. What is here
   !! is the wiring, which is where this project's mistakes would be -- units,
   !! atom order, the functional name, and the sign of what is added.
   !!
   !! **Units.** mqc is Bohr and Hartree throughout, and so is this API:
   !! `s-dftd3.h` says "quantities in Bohr" on the structure constructor, and
   !! `new_structure_api` declares `positions(3, natoms)` in the same storage
   !! order mqc uses. Nothing is scaled on the way in or out. The gradient comes
   !! back as dE/dR in Hartree/Bohr with the same (3, natoms) shape.
   !!
   !! **The C surface, not the Fortran one.** s-dftd3 is LGPL-3-or-later against
   !! this project's MIT, and is linked as a shared library so that the two stay
   !! separable. Using its .mod files would fix the compiler as well, so the
   !! interface below is declared here rather than imported. It mirrors
   !! `include/s-dftd3.h`; nothing in this build reads that header.
   !!
   !! Compiled only when MQC_ENABLE_DFTD3=ON. The twin that stands in otherwise
   !! is `src/methods/dft/mqc_dispersion_stub.f90`, same module name, same two
   !! public procedures.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_bool, c_null_ptr, c_null_char, c_loc, &
                                                                             c_associated
   use pic_types, only: dp, default_int
   use mqc_error, only: error_t, ERROR_GENERIC, ERROR_VALIDATION
   use mqc_dispersion_names, only: d3_functional_alias, dispersion_kind_is_known, DISPERSION_KINDS
   implicit none
   private

   public :: dispersion_available
   public :: dispersion_correction

   integer(c_int), parameter :: ERROR_BUFFER = 512_c_int

   ! The C API of s-dftd3, transcribed from include/s-dftd3.h.
   !
   ! Every handle is an opaque `type(c_ptr)`. The `delete_*` entry points take a
   ! pointer *to* the handle and null it, which is why those dummies are
   ! intent(inout) and not `value`.
   !
   ! `lattice`, `periodic`, `gradient` and `sigma` are optional on the C side,
   ! which for a bind(C) dummy means "a null pointer is an absent argument".
   ! They are declared `type(c_ptr), value` here so that this side decides
   ! explicitly, by passing `c_null_ptr` or `c_loc(...)`, rather than relying on
   ! how a particular compiler passes an absent optional.
   interface

      function dftd3_new_error() bind(c, name="dftd3_new_error") result(handle)
         import :: c_ptr
         implicit none
         type(c_ptr) :: handle
      end function dftd3_new_error

      function dftd3_check_error(handle) bind(c, name="dftd3_check_error") result(status)
         import :: c_ptr, c_int
         implicit none
         type(c_ptr), value :: handle
         integer(c_int) :: status
      end function dftd3_check_error

      subroutine dftd3_get_error(handle, buffer, buffersize) bind(c, name="dftd3_get_error")
         import :: c_ptr, c_char, c_int
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(in) :: buffersize
         character(kind=c_char), intent(inout) :: buffer(buffersize)
      end subroutine dftd3_get_error

      subroutine dftd3_delete_error(handle) bind(c, name="dftd3_delete_error")
         import :: c_ptr
         implicit none
         type(c_ptr), intent(inout) :: handle
      end subroutine dftd3_delete_error

      function dftd3_new_structure(handle, natoms, numbers, positions, lattice, periodic) &
         bind(c, name="dftd3_new_structure") result(mol)
         import :: c_ptr, c_int, c_double
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: natoms
         integer(c_int), intent(in) :: numbers(natoms)
         real(c_double), intent(in) :: positions(3, natoms)
         type(c_ptr), value :: lattice
         type(c_ptr), value :: periodic
         type(c_ptr) :: mol
      end function dftd3_new_structure

      subroutine dftd3_delete_structure(mol) bind(c, name="dftd3_delete_structure")
         import :: c_ptr
         implicit none
         type(c_ptr), intent(inout) :: mol
      end subroutine dftd3_delete_structure

      function dftd3_new_d3_model(handle, mol) bind(c, name="dftd3_new_d3_model") result(model)
         import :: c_ptr
         implicit none
         type(c_ptr), value :: handle
         type(c_ptr), value :: mol
         type(c_ptr) :: model
      end function dftd3_new_d3_model

      subroutine dftd3_delete_model(model) bind(c, name="dftd3_delete_model")
         import :: c_ptr
         implicit none
         type(c_ptr), intent(inout) :: model
      end subroutine dftd3_delete_model

      function dftd3_load_rational_damping(handle, method, atm) &
         bind(c, name="dftd3_load_rational_damping") result(param)
         import :: c_ptr, c_bool
         implicit none
         type(c_ptr), value :: handle
         type(c_ptr), value :: method
            !! A null-terminated C string. `type(c_ptr)` rather than an
            !! assumed-size `character(kind=c_char)` array so that nothing here
            !! declares an array whose extent it does not know.
         logical(c_bool), value :: atm
         type(c_ptr) :: param
      end function dftd3_load_rational_damping

      subroutine dftd3_delete_param(param) bind(c, name="dftd3_delete_param")
         import :: c_ptr
         implicit none
         type(c_ptr), intent(inout) :: param
      end subroutine dftd3_delete_param

      subroutine dftd3_get_dispersion(handle, mol, model, param, energy, gradient, sigma) &
         bind(c, name="dftd3_get_dispersion")
         import :: c_ptr, c_double
         implicit none
         type(c_ptr), value :: handle
         type(c_ptr), value :: mol
         type(c_ptr), value :: model
         type(c_ptr), value :: param
         real(c_double), intent(out) :: energy
         type(c_ptr), value :: gradient
         type(c_ptr), value :: sigma
      end subroutine dftd3_get_dispersion

   end interface

contains

   pure function dispersion_available() result(available)
      !! .true. -- this build linked s-dftd3
      !!
      !! Read by `read_dft_dispersion`, which refuses a deck naming the keyword
      !! on a build without the library, and by the version banner.
      logical :: available

      available = .true.
   end function dispersion_available

   subroutine dispersion_correction(kind, functional, atomic_numbers, coordinates, &
                                    energy, gradient, error)
      !! The dispersion energy, and its gradient when one is asked for
      !!
      !! `coordinates` is (3, natoms) in Bohr and `gradient` comes back in the
      !! same shape in Hartree/Bohr, both being what s-dftd3 already uses. The
      !! energy is negative for any geometry with two atoms in it.
      character(len=*), intent(in) :: kind
         !! Which correction: only "d3bj" today. See `DISPERSION_KINDS`.
      character(len=*), intent(in) :: functional
         !! This program's spelling, from `model.functional`. Translated by
         !! `d3_functional_alias`, which refuses rather than guesses.
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
      real(c_double) :: e_disp
      character(kind=c_char), allocatable, target :: method(:)
      character(len=32) :: alias
      logical :: want_gradient

      energy = 0.0_dp
      natoms = size(atomic_numbers, kind=default_int)
      want_gradient = present(gradient)
      if (want_gradient) gradient = 0.0_dp

      if (.not. dispersion_kind_is_known(kind)) then
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

      call d3_functional_alias(functional, alias, error)
      if (error%has_error()) return

      ! Fixed-kind copies: `c_int` and `c_double` need not be the default
      ! integer and `dp`, and a bind(C) dummy is not a place to find that out.
      ! Copies rather than views also mean `coordinates` stays intent(in) in
      ! fact as well as in declaration.
      allocate (numbers(natoms))
      numbers = int(atomic_numbers, c_int)
      allocate (positions(3, natoms))
      positions = real(coordinates, c_double)
      method = c_string(trim(alias))

      handle = dftd3_new_error()
      mol = dftd3_new_structure(handle, int(natoms, c_int), numbers, positions, &
                                c_null_ptr, c_null_ptr)
      if (failed(handle, "building the structure", error)) then
         call dftd3_delete_error(handle)
         return
      end if

      model = dftd3_new_d3_model(handle, mol)
      if (failed(handle, "building the D3 model", error)) then
         call dftd3_delete_structure(mol)
         call dftd3_delete_error(handle)
         return
      end if

      ! atm = .false., which is s-dftd3's own default and the one its Python
      ! bindings use: "-D3(BJ)" names the two-body correction, and the
      ! three-body Axilrod-Teller-Muto term is a separate choice with its own
      ! published parameters. Turning it on silently would make this program's
      ! "d3bj" a different number from everyone else's.
      param = dftd3_load_rational_damping(handle, c_loc(method), .false._c_bool)
      if (failed(handle, "loading D3(BJ) damping parameters for '"//trim(alias)//"'", error)) then
         call dftd3_delete_model(model)
         call dftd3_delete_structure(mol)
         call dftd3_delete_error(handle)
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
         call dftd3_get_dispersion(handle, mol, model, param, e_disp, &
                                   c_loc(forces), c_loc(virial))
      else
         call dftd3_get_dispersion(handle, mol, model, param, e_disp, &
                                   c_null_ptr, c_null_ptr)
      end if

      if (.not. failed(handle, "evaluating the dispersion energy", error)) then
         energy = real(e_disp, dp)
         if (want_gradient) gradient = real(forces, dp)
      end if

      call dftd3_delete_param(param)
      call dftd3_delete_model(model)
      call dftd3_delete_structure(mol)
      call dftd3_delete_error(handle)
   end subroutine dispersion_correction

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
      !! Whether s-dftd3 set its error handle, and if so what it said
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
         call error%set(ERROR_GENERIC, "s-dftd3 would not create an error handle")
         bad = .true.
         return
      end if
      if (dftd3_check_error(handle) == 0_c_int) return

      buffer = c_null_char
      call dftd3_get_error(handle, buffer, ERROR_BUFFER)
      message = ""
      do i = 1, int(ERROR_BUFFER)
         if (buffer(i) == c_null_char) exit
         message(i:i) = buffer(i)
      end do

      call error%set(ERROR_VALIDATION, "s-dftd3 refused "//trim(doing)//": "//trim(message))
      bad = .true.
   end function failed

end module mqc_dispersion
