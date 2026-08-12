!! Cartesian multipole moment integrals over the AO basis
module mqc_libcint_multipole
   !! Dipole, quadrupole and octopole integrals about a chosen origin.
   !!
   !! `<chi_u| (x-O)^a (y-O)^b (z-O)^c |chi_v>`, for a+b+c = 1, 2 and 3. libcint
   !! computes these through `int1e_r`, `int1e_rr` and `int1e_rrr`, which are the
   !! same shell-pair machinery as the overlap with one extra factor -- so this
   !! module is bindings and bookkeeping, not new integral code.
   !!
   !! **The entry points are declared here rather than imported.** The Fortran
   !! interface shipped with the libcint fork wraps overlap, kinetic, nuclear
   !! attraction and a handful of derivatives, but not the multipoles. Declaring
   !! the three `bind(C)` interfaces locally keeps the addition inside this
   !! repository: no fork change, no version of libcint we have to wait for. The
   !! symbols are in the library either way -- `cint1e_r_sph` and friends -- and
   !! the signature is the same seven-argument form the fork's own wrappers use --
   !! except that the array arguments are declared `type(c_ptr), value` and passed
   !! through `c_loc`, rather than as assumed-size dummies. That is what a C
   !! pointer actually is, it matches the cuEST bindings, and it keeps this file
   !! inside the prohibition on assumed-size arrays in FORTRAN_STYLE.
   !!
   !! **The origin is passed through `env`, and that is the trap.** libcint reads
   !! it from `env(PTR_COMMON_ORIG)`, and the slot constants are libcint's own
   !! 0-based offsets which the Fortran interface does *not* convert -- unlike the
   !! `atm`/`bas` column constants, which it does. So the index carries a `+ 1`.
   !! Getting it wrong is silent: the moments come back expanded about whatever
   !! happened to be in the neighbouring slot, usually the origin, which for a
   !! neutral molecule's dipole is the right answer anyway and for every
   !! quadrupole is not. The same mistake cost a day on the range-separated
   !! exchange, where omega went into `PTR_RINV_ZETA` and was ignored.
   use pic_types, only: dp
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_loc
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim
   use libcint_fortran, only: LIBCINT_PTR_COMMON_ORIG
   implicit none
   private

   public :: multipole_matrices
   public :: DIPOLE_COMPONENTS, QUADRUPOLE_COMPONENTS, OCTOPOLE_COMPONENTS

   !> Components libcint returns, as full Cartesian tensors.
   !>
   !> Full, not packed: 9 quadrupole components rather than the 6 unique ones and
   !> 27 octopole rather than 10. Anything writing a GAMESS-style potential has to
   !> pack them, and that packing is a place a tensor transposes without
   !> complaint -- so it belongs to whoever writes the file, with a test, and not
   !> here where it would be invisible.
   integer, parameter :: DIPOLE_COMPONENTS = 3
   integer, parameter :: QUADRUPOLE_COMPONENTS = 9
   integer, parameter :: OCTOPOLE_COMPONENTS = 27

   interface
      function cint1e_r_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_r_sph")
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value, intent(in) :: buf
         type(c_ptr), value, intent(in) :: shls
         type(c_ptr), value, intent(in) :: atm
         integer(c_int), value, intent(in) :: natm
         type(c_ptr), value, intent(in) :: bas
         integer(c_int), value, intent(in) :: nbas
         type(c_ptr), value, intent(in) :: env
         integer(c_int) :: ret
      end function cint1e_r_sph

      function cint1e_r_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_r_cart")
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value, intent(in) :: buf
         type(c_ptr), value, intent(in) :: shls
         type(c_ptr), value, intent(in) :: atm
         integer(c_int), value, intent(in) :: natm
         type(c_ptr), value, intent(in) :: bas
         integer(c_int), value, intent(in) :: nbas
         type(c_ptr), value, intent(in) :: env
         integer(c_int) :: ret
      end function cint1e_r_cart

      function cint1e_rr_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rr_sph")
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value, intent(in) :: buf
         type(c_ptr), value, intent(in) :: shls
         type(c_ptr), value, intent(in) :: atm
         integer(c_int), value, intent(in) :: natm
         type(c_ptr), value, intent(in) :: bas
         integer(c_int), value, intent(in) :: nbas
         type(c_ptr), value, intent(in) :: env
         integer(c_int) :: ret
      end function cint1e_rr_sph

      function cint1e_rr_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rr_cart")
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value, intent(in) :: buf
         type(c_ptr), value, intent(in) :: shls
         type(c_ptr), value, intent(in) :: atm
         integer(c_int), value, intent(in) :: natm
         type(c_ptr), value, intent(in) :: bas
         integer(c_int), value, intent(in) :: nbas
         type(c_ptr), value, intent(in) :: env
         integer(c_int) :: ret
      end function cint1e_rr_cart

      function cint1e_rrr_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rrr_sph")
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value, intent(in) :: buf
         type(c_ptr), value, intent(in) :: shls
         type(c_ptr), value, intent(in) :: atm
         integer(c_int), value, intent(in) :: natm
         type(c_ptr), value, intent(in) :: bas
         integer(c_int), value, intent(in) :: nbas
         type(c_ptr), value, intent(in) :: env
         integer(c_int) :: ret
      end function cint1e_rrr_sph

      function cint1e_rrr_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rrr_cart")
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value, intent(in) :: buf
         type(c_ptr), value, intent(in) :: shls
         type(c_ptr), value, intent(in) :: atm
         integer(c_int), value, intent(in) :: natm
         type(c_ptr), value, intent(in) :: bas
         integer(c_int), value, intent(in) :: nbas
         type(c_ptr), value, intent(in) :: env
         integer(c_int) :: ret
      end function cint1e_rrr_cart
   end interface

contains

   subroutine multipole_matrices(mol, origin, order, matrices, error)
      !! Multipole integrals of one order, over every shell pair
      !!
      !! Returns `(n_ao, n_ao, n_components)`, the components in the order
      !! libcint produces them: x,y,z for the dipole; xx,xy,xz,yx,...,zz for the
      !! quadrupole; and the 27 of the octopole with z fastest. That is the *full*
      !! Cartesian tensor, deliberately -- see the component parameters above.
      type(libcint_molecule_t), intent(in), target :: mol
         !! `target` so its atm/bas/env can be addressed for the C call.
      real(dp), intent(in) :: origin(3)
         !! Expansion origin, Bohr. A dipole is origin-independent only for a
         !! neutral charge distribution; every higher moment depends on it, so
         !! there is no sensible default and the caller states one.
      integer, intent(in) :: order                     !! 1, 2 or 3
      real(dp), allocatable, intent(out) :: matrices(:, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable, target :: buf(:), env_local(:)
      integer, target :: shls(2)
      integer :: ish, jsh, di, dj, i, j, k, ret, io, jo, block_max, n_comp
      character(len=8) :: order_text

      if (order < 1 .or. order > 3) then
         write (order_text, "(i0)") order
         call error%set(ERROR_VALIDATION, "multipole order "//trim(order_text)// &
                        " is not available; libcint is called here for the dipole, "// &
                        "quadrupole and octopole only (1, 2 and 3).")
         return
      end if

      select case (order)
      case (1)
         n_comp = DIPOLE_COMPONENTS
      case (2)
         n_comp = QUADRUPOLE_COMPONENTS
      case default
         n_comp = OCTOPOLE_COMPONENTS
      end select

      ! A copy, because the origin lives in `env` and the molecule is read-only
      ! here -- the same reason the range-separated exchange copies it. `+ 1`
      ! because the slot constant is libcint's own 0-based offset and the Fortran
      ! interface does not convert it.
      env_local = mol%env
      env_local(LIBCINT_PTR_COMMON_ORIG + 1:LIBCINT_PTR_COMMON_ORIG + 3) = origin

      block_max = 1
      do ish = 1, mol%nbas
         block_max = max(block_max, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do

      allocate (matrices(mol%nao, mol%nao, n_comp))
      matrices = 0.0_dp
      ! One shell-pair block holds di*dj elements per component.
      allocate (buf(block_max*block_max*n_comp))

      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, mol%nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]

            if (mol%cartesian) then
               select case (order)
               case (1)
                  ret = cint1e_r_cart(c_loc(buf), c_loc(shls), c_loc(mol%atm), mol%natm, &
                                      c_loc(mol%bas), mol%nbas, c_loc(env_local))
               case (2)
                  ret = cint1e_rr_cart(c_loc(buf), c_loc(shls), c_loc(mol%atm), mol%natm, &
                                       c_loc(mol%bas), mol%nbas, c_loc(env_local))
               case default
                  ret = cint1e_rrr_cart(c_loc(buf), c_loc(shls), c_loc(mol%atm), mol%natm, &
                                        c_loc(mol%bas), mol%nbas, c_loc(env_local))
               end select
            else
               select case (order)
               case (1)
                  ret = cint1e_r_sph(c_loc(buf), c_loc(shls), c_loc(mol%atm), mol%natm, &
                                     c_loc(mol%bas), mol%nbas, c_loc(env_local))
               case (2)
                  ret = cint1e_rr_sph(c_loc(buf), c_loc(shls), c_loc(mol%atm), mol%natm, &
                                      c_loc(mol%bas), mol%nbas, c_loc(env_local))
               case default
                  ret = cint1e_rrr_sph(c_loc(buf), c_loc(shls), c_loc(mol%atm), mol%natm, &
                                       c_loc(mol%bas), mol%nbas, c_loc(env_local))
               end select
            end if
            if (ret == 0) cycle   ! screened away, leave the block zero

            ! libcint fills the block with the bra index fastest, then the ket,
            ! then the component.
            do k = 1, n_comp
               do j = 1, dj
                  do i = 1, di
                     matrices(io + i, jo + j, k) = &
                        buf(i + (j - 1)*di + (k - 1)*di*dj)
                  end do
               end do
            end do
         end do
      end do

      deallocate (buf, env_local)
   end subroutine multipole_matrices

end module mqc_libcint_multipole
