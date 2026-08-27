!! Effective core potential integrals, over libfint's ECP entry points
module mqc_libcint_ecp
   !! An ECP replaces an atom's core electrons with a potential, which buys two
   !! things at once: the core is gone from the SCF, and the potential can be
   !! fitted to reproduce a relativistic calculation. Below about the fourth
   !! row the second matters more than the first.
   !!
   !! **The integrals are not libcint's.** libcint has no ECP code at all;
   !! these come from libfint, which translated them from PySCF's
   !! `pyscf/lib/gto/nr_ecp.c`. So this is the one integral family where the
   !! CPU backend cannot fall back to libcint, and a build against a libfint
   !! without `ECPscalar_sph` fails to link rather than misbehaving.
   !!
   !! **How the potential reaches the library.** Not as an argument. libcint's
   !! entry points take `atm`, `bas`, `nbas` and `env` and nothing else, so
   !! PySCF's convention -- which libfint follows -- puts the ECP shells in
   !! `bas` after the orbital shells and records where in two otherwise unused
   !! `env` slots:
   !!
   !!     env(AS_ECPBAS_OFFSET) = the first ECP row, 0-based
   !!     env(AS_NECPBAS)       = how many there are
   !!
   !! `nbas` stays the orbital count. `mqc_libcint_integrals` lays that out;
   !! this module only reads it.
   !!
   !! **What an ECP changes elsewhere, and this module does not do.** The
   !! nuclear charge an atom presents to everything else drops by the number
   !! of core electrons replaced, which changes the nuclear attraction, the
   !! nuclear repulsion and the electron count. Those belong to the molecule
   !! and the fragment respectively, and doing them here would be doing them
   !! twice.
   use pic_types, only: dp
   use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_loc, c_null_ptr
   implicit none
   private

   public :: ecp_matrix

   interface
      function cECPscalar_sph(buf, dims, shls, atm, natm, bas, nbas, env, opt, &
                              cache) result(ret) bind(C, name="ECPscalar_sph")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value, intent(in) :: buf, dims, shls, atm
         integer(c_int), value, intent(in) :: natm
         type(c_ptr), value, intent(in) :: bas
         integer(c_int), value, intent(in) :: nbas
         type(c_ptr), value, intent(in) :: env, opt, cache
         integer(c_int) :: ret
      end function cECPscalar_sph

      function cECPscalar_cart(buf, dims, shls, atm, natm, bas, nbas, env, opt, &
                               cache) result(ret) bind(C, name="ECPscalar_cart")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value, intent(in) :: buf, dims, shls, atm
         integer(c_int), value, intent(in) :: natm
         type(c_ptr), value, intent(in) :: bas
         integer(c_int), value, intent(in) :: nbas
         type(c_ptr), value, intent(in) :: env, opt, cache
         integer(c_int) :: ret
      end function cECPscalar_cart
   end interface

contains

   subroutine ecp_matrix(nao, nbas, natm, cartesian, atm, bas, env, &
                         shell_offset, necpbas, matrix)
      !! The ECP contribution to the core Hamiltonian, (nao, nao)
      !!
      !! Comes back zero when `necpbas` is zero, rather than refusing: the
      !! caller adds this unconditionally and a zero matrix is the right answer
      !! for a molecule no potential touches.
      !!
      !! **Plain arrays rather than `libcint_molecule_t`**, which is not a
      !! style preference. `molecule_core_hamiltonian` has to add this term,
      !! and taking the type here would make the two modules mutually
      !! dependent -- so the molecule passes its own arrays and this module
      !! depends on nothing of mqc's.
      !!
      !! `bas` is the table *with* the ECP rows appended, while `nbas` is the
      !! orbital count alone. Those disagreeing is the whole convention.
      !!
      !! The loop mirrors `one_electron` in `mqc_libcint_integrals` -- same
      !! shell-pair walk, same column-major block, same "return zero means
      !! screened" reading. Shell extents come from `shell_offset` differences
      !! rather than from `shell_dim`, which keeps this free of that module.
      integer, intent(in) :: nao, nbas, natm
      logical, intent(in) :: cartesian
      integer, intent(in) :: atm(:, :)
      integer, intent(in) :: bas(:, :)     !! (BAS_SLOTS, nbas + necpbas)
      real(dp), intent(in) :: env(:)
      integer, intent(in) :: shell_offset(:)  !! (nbas + 1), first AO per shell
      integer, intent(in) :: necpbas
      real(dp), allocatable, intent(out) :: matrix(:, :)

      real(dp), allocatable, target :: buf(:)
      integer(c_int), allocatable, target :: shls(:)
      integer(c_int), allocatable, target :: atm_c(:), bas_c(:)
      real(dp), allocatable, target :: env_c(:)
      integer :: ish, jsh, di, dj, i, j, io, jo, ret, nmax

      allocate (matrix(nao, nao))
      matrix = 0.0_dp
      if (necpbas <= 0) return

      ! A block is di*dj, so the bound is the square of the largest shell and
      ! not the largest shell -- the distinction only shows up once a shell
      ! beyond s is present, by which time it is heap corruption.
      nmax = 0
      do ish = 1, nbas
         nmax = max(nmax, shell_offset(ish + 1) - shell_offset(ish))
      end do
      allocate (buf(nmax*nmax))
      allocate (shls(2))

      atm_c = reshape(int(atm, c_int), [size(atm)])
      bas_c = reshape(int(bas, c_int), [size(bas)])
      env_c = env

      do ish = 1, nbas
         di = shell_offset(ish + 1) - shell_offset(ish)
         io = shell_offset(ish)
         do jsh = 1, nbas
            dj = shell_offset(jsh + 1) - shell_offset(jsh)
            jo = shell_offset(jsh)
            shls = int([ish - 1, jsh - 1], c_int)
            buf(1:di*dj) = 0.0_dp

            if (cartesian) then
               ret = cECPscalar_cart(c_loc(buf), c_null_ptr, c_loc(shls), &
                                     c_loc(atm_c), int(natm, c_int), &
                                     c_loc(bas_c), int(nbas, c_int), &
                                     c_loc(env_c), c_null_ptr, c_null_ptr)
            else
               ret = cECPscalar_sph(c_loc(buf), c_null_ptr, c_loc(shls), &
                                    c_loc(atm_c), int(natm, c_int), &
                                    c_loc(bas_c), int(nbas, c_int), &
                                    c_loc(env_c), c_null_ptr, c_null_ptr)
            end if
            ! Zero means the overlap screen rejected the pair and the library
            ! zeroed the block. Skipping the copy is an optimisation; the
            ! block is correct either way.
            if (ret == 0) cycle

            do j = 1, dj
               do i = 1, di
                  matrix(io + i, jo + j) = buf(i + (j - 1)*di)
               end do
            end do
         end do
      end do
   end subroutine ecp_matrix

end module mqc_libcint_ecp
