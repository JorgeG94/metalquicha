!! Second derivatives of the one- and two-electron integrals with respect to nuclei
module mqc_libcint_hess_ints
   !! What an analytic Hessian needs from the integral library, and nothing else.
   !!
   !! **These come from libfint's native modules rather than through the
   !! `libcint_fortran` compatibility layer.** That layer exists so the code
   !! written against libcint's C ABI kept working when the backend swapped, and
   !! it deliberately exposes only the surface that code used -- first
   !! derivatives and below. Nothing here has that history to preserve, so it
   !! calls the Fortran directly and gets the better interface for it: a
   !! `cint_ws` workspace instead of a hand-sized cache, a logical result
   !! instead of a return code, and no C optimizer lifecycle to leak.
   !!
   !! The consequence is that analytic Hessians are available only in a libfint
   !! build, which is the default. A libcint build has no way to reach these
   !! entry points without a `bind(C)` declaration per integral, and keeps the
   !! finite-difference path instead.
   !!
   !! **Both orderings of every derivative are needed, and they are different
   !! integrals.** `ipipovlp` puts both derivatives on the bra centre;
   !! `ipovlpip` puts one on each. A Hessian built from only the first is not a
   !! slightly worse Hessian, it is one that violates translational invariance,
   !! which shows up as translations that are no longer zero-frequency modes.
   use, intrinsic :: iso_fortran_env, only: real64
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim
   use cint_workspace, only: cint_ws
   use cint_bas, only: cint_cgto_spheric, cint_cgto_cart
   use cint_gen_hess, only: int1e_ipipovlp_sph, int1e_ipipovlp_cart, &
                            int1e_ipovlpip_sph, int1e_ipovlpip_cart, &
                            int1e_ipipkin_sph, int1e_ipipkin_cart, &
                            int1e_ipkinip_sph, int1e_ipkinip_cart, &
                            int1e_ipipnuc_sph, int1e_ipipnuc_cart, &
                            int1e_ipnucip_sph, int1e_ipnucip_cart
   use cint_gen_hess, only: int2e_ipip1_sph, int2e_ipip1_cart, &
                            int2e_ipvip1_sph, int2e_ipvip1_cart, &
                            int2e_ip1ip2_sph, int2e_ip1ip2_cart
   implicit none
   private

   public :: HESS_OVLP_II, HESS_OVLP_IJ, HESS_KIN_II, HESS_KIN_IJ
   public :: HESS_NUC_II, HESS_NUC_IJ
   public :: hess_1e_block
   public :: HESS_ERI_II, HESS_ERI_IJ, HESS_ERI_IK
   public :: hess_2e_block

   integer, parameter :: HESS_OVLP_II = 1   !! `int1e_ipipovlp`, both derivatives on the bra
   integer, parameter :: HESS_OVLP_IJ = 2   !! `int1e_ipovlpip`, one on each centre
   integer, parameter :: HESS_KIN_II = 3
   integer, parameter :: HESS_KIN_IJ = 4
   integer, parameter :: HESS_NUC_II = 5
   integer, parameter :: HESS_NUC_IJ = 6

   integer, parameter :: HESS_ERI_II = 1   !! `int2e_ipip1`, both derivatives on centre 1
   integer, parameter :: HESS_ERI_IJ = 2   !! `int2e_ipvip1`, one on centre 1 and one on centre 2
   integer, parameter :: HESS_ERI_IK = 3   !! `int2e_ip1ip2`, one on centre 1 and one on centre 3

   integer, parameter :: N_COMPONENTS = 9
      !! Every one of these is a 3x3 Cartesian block per shell pair, laid out by
      !! libcint as `xx xy xz yx yy yz zx zy zz` in the slowest index.

   !> One workspace per thread. **`threadprivate` is not optional** -- libfint's
   !> own C layer says so and says how it was found: a shared workspace means
   !> every thread writing over every other one's scratch, and it does not fail
   !> cleanly but segfaults somewhere unrelated.
   type(cint_ws), save :: ws
   !$omp threadprivate(ws)

contains

   subroutine hess_1e_block(mol, which, matrix, error)
      !! One second-derivative one-electron integral over the whole basis
      !!
      !! Returns `(n_ao, n_ao, 9)`. The nine components are the Cartesian pairs
      !! in libcint's order, so component `3*(a-1)+b` is the derivative with
      !! respect to `a` on the first centre and `b` on whichever centre this
      !! particular integral differentiates second.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: which          !! One of the `HESS_*` selectors
      real(dp), allocatable, intent(out) :: matrix(:, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: buf(:)
      ! The native entry points take `atm` and `bas` flat, where the
      ! compatibility layer took them shaped. Flattened once here rather than
      ! per shell pair: it is a few hundred integers against a loop that runs
      ! `nbas^2` times.
      integer, allocatable :: atm_flat(:), bas_flat(:)
      integer :: dims(0:3), shls(0:1)
      integer :: ish, jsh, di, dj, io, jo, i, j, comp, mx
      logical :: have

      if (error%has_error()) return

      mx = 0
      do ish = 1, mol%nbas
         mx = max(mx, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do

      ! Flat and indexed by hand, for the reason the gradient driver gives: the
      ! library packs a block with the *actual* shell dimensions, so a rank-3
      ! buffer declared at the largest shell has the wrong strides for every
      ! smaller one -- invisible in a basis of s functions and wrong everywhere
      ! else.
      allocate (buf(mx*mx*N_COMPONENTS))
      atm_flat = reshape(mol%atm, [size(mol%atm)])
      bas_flat = reshape(mol%bas, [size(mol%bas)])
      allocate (matrix(mol%nao, mol%nao, N_COMPONENTS))
      matrix = 0.0_dp

      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, mol%nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]
            dims = [di, dj, 1, 1]

            have = drive(mol, which, buf, dims, shls, atm_flat, bas_flat)
            if (.not. have) cycle

            do comp = 1, N_COMPONENTS
               do j = 1, dj
                  do i = 1, di
                     matrix(io + i, jo + j, comp) = &
                        buf(i + di*(j - 1 + dj*(comp - 1)))
                  end do
               end do
            end do
         end do
      end do

      deallocate (buf, atm_flat, bas_flat)
   end subroutine hess_1e_block

   function drive(mol, which, buf, dims, shls, atm, bas) result(have)
      !! Dispatch one shell pair to the right entry point
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: which
      real(dp), intent(inout) :: buf(0:)
      integer, intent(in) :: dims(0:), shls(0:)
      integer, intent(in), target :: atm(0:), bas(0:)
      logical :: have

      have = .false.
      if (mol%cartesian) then
         select case (which)
         case (HESS_OVLP_II)
            have = int1e_ipipovlp_cart(buf, dims, shls, atm, mol%natm, &
                                       bas, mol%nbas, mol%env, ws)
         case (HESS_OVLP_IJ)
            have = int1e_ipovlpip_cart(buf, dims, shls, atm, mol%natm, &
                                       bas, mol%nbas, mol%env, ws)
         case (HESS_KIN_II)
            have = int1e_ipipkin_cart(buf, dims, shls, atm, mol%natm, &
                                      bas, mol%nbas, mol%env, ws)
         case (HESS_KIN_IJ)
            have = int1e_ipkinip_cart(buf, dims, shls, atm, mol%natm, &
                                      bas, mol%nbas, mol%env, ws)
         case (HESS_NUC_II)
            have = int1e_ipipnuc_cart(buf, dims, shls, atm, mol%natm, &
                                      bas, mol%nbas, mol%env, ws)
         case default
            have = int1e_ipnucip_cart(buf, dims, shls, atm, mol%natm, &
                                      bas, mol%nbas, mol%env, ws)
         end select
      else
         select case (which)
         case (HESS_OVLP_II)
            have = int1e_ipipovlp_sph(buf, dims, shls, atm, mol%natm, &
                                      bas, mol%nbas, mol%env, ws)
         case (HESS_OVLP_IJ)
            have = int1e_ipovlpip_sph(buf, dims, shls, atm, mol%natm, &
                                      bas, mol%nbas, mol%env, ws)
         case (HESS_KIN_II)
            have = int1e_ipipkin_sph(buf, dims, shls, atm, mol%natm, &
                                     bas, mol%nbas, mol%env, ws)
         case (HESS_KIN_IJ)
            have = int1e_ipkinip_sph(buf, dims, shls, atm, mol%natm, &
                                     bas, mol%nbas, mol%env, ws)
         case (HESS_NUC_II)
            have = int1e_ipipnuc_sph(buf, dims, shls, atm, mol%natm, &
                                     bas, mol%nbas, mol%env, ws)
         case default
            have = int1e_ipnucip_sph(buf, dims, shls, atm, mol%natm, &
                                     bas, mol%nbas, mol%env, ws)
         end select
      end if
   end function drive

   subroutine hess_2e_block(mol, which, eri, error)
      !! One second-derivative two-electron integral over the whole basis
      !!
      !! Returns `(n_ao, n_ao, n_ao, n_ao, 9)`, which is `n_ao^4` nine times
      !! over and is why this exists as a first correct implementation rather
      !! than a usable one. The Hessian contracts these against densities the
      !! moment they are built, so nothing needs the whole array at once; a
      !! shell-driven contraction is the obvious next step and would change
      !! nothing above it.
      !!
      !! **The three selectors are three genuinely different integrals**, not
      !! one integral indexed three ways. `ipip1` puts both derivatives on the
      !! first centre, `ipvip1` puts one each on the bra pair, and `ip1ip2`
      !! puts one on the bra and one on the ket. Translational invariance is
      !! what relates them, and it relates them only if all three are present.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: which
      real(dp), allocatable, intent(out) :: eri(:, :, :, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: buf(:)
      integer, allocatable :: atm_flat(:), bas_flat(:)
      integer :: dims(0:3), shls(0:3)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: io, jo, ko, lo, i, j, k, l, comp, mx, idx
      logical :: have

      if (error%has_error()) return

      mx = 0
      do ish = 1, mol%nbas
         mx = max(mx, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do

      allocate (buf(mx**4*N_COMPONENTS))
      allocate (eri(mol%nao, mol%nao, mol%nao, mol%nao, N_COMPONENTS))
      eri = 0.0_dp
      atm_flat = reshape(mol%atm, [size(mol%atm)])
      bas_flat = reshape(mol%bas, [size(mol%bas)])

      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, mol%nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            do ksh = 1, mol%nbas
               dk = shell_dim(mol%cartesian, ksh - 1, mol%bas)
               ko = mol%shell_offset(ksh)
               do lsh = 1, mol%nbas
                  dl = shell_dim(mol%cartesian, lsh - 1, mol%bas)
                  lo = mol%shell_offset(lsh)
                  shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]
                  dims = [di, dj, dk, dl]

                  have = drive_2e(mol, which, buf, dims, shls, atm_flat, bas_flat)
                  if (.not. have) cycle

                  do comp = 1, N_COMPONENTS
                     do l = 1, dl
                        do k = 1, dk
                           do j = 1, dj
                              do i = 1, di
                                 idx = i + di*(j - 1 + dj*(k - 1 + dk*(l - 1 &
                                                                       + dl*(comp - 1))))
                                 eri(io + i, jo + j, ko + k, lo + l, comp) = buf(idx)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      deallocate (buf, atm_flat, bas_flat)
   end subroutine hess_2e_block

   function drive_2e(mol, which, buf, dims, shls, atm, bas) result(have)
      !! Dispatch one shell quartet to the right entry point
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: which
      real(dp), intent(inout) :: buf(0:)
      integer, intent(in) :: dims(0:), shls(0:)
      integer, intent(in), target :: atm(0:), bas(0:)
      logical :: have

      have = .false.
      if (mol%cartesian) then
         select case (which)
         case (HESS_ERI_II)
            have = int2e_ipip1_cart(buf, dims, shls, atm, mol%natm, &
                                    bas, mol%nbas, mol%env, ws)
         case (HESS_ERI_IJ)
            have = int2e_ipvip1_cart(buf, dims, shls, atm, mol%natm, &
                                     bas, mol%nbas, mol%env, ws)
         case default
            have = int2e_ip1ip2_cart(buf, dims, shls, atm, mol%natm, &
                                     bas, mol%nbas, mol%env, ws)
         end select
      else
         select case (which)
         case (HESS_ERI_II)
            have = int2e_ipip1_sph(buf, dims, shls, atm, mol%natm, &
                                   bas, mol%nbas, mol%env, ws)
         case (HESS_ERI_IJ)
            have = int2e_ipvip1_sph(buf, dims, shls, atm, mol%natm, &
                                    bas, mol%nbas, mol%env, ws)
         case default
            have = int2e_ip1ip2_sph(buf, dims, shls, atm, mol%natm, &
                                    bas, mol%nbas, mol%env, ws)
         end select
      end if
   end function drive_2e

end module mqc_libcint_hess_ints
