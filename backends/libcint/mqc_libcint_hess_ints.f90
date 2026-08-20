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
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim, atom_ao_blocks
   use cint_workspace, only: cint_ws
   use cint_bas, only: cint_cgto_spheric, cint_cgto_cart
   use cint_gen_hess, only: int1e_ipipovlp_sph, int1e_ipipovlp_cart, &
                            int1e_ipovlpip_sph, int1e_ipovlpip_cart, &
                            int1e_ipipkin_sph, int1e_ipipkin_cart, &
                            int1e_ipkinip_sph, int1e_ipkinip_cart, &
                            int1e_ipipnuc_sph, int1e_ipipnuc_cart, &
                            int1e_ipnucip_sph, int1e_ipnucip_cart
   use cint_gen_hess, only: int1e_ipiprinv_sph, int1e_ipiprinv_cart, &
                            int1e_iprinvip_sph, int1e_iprinvip_cart
   use cint_envs, only: PTR_RINV_ORIG
   use cint_gen_grad2, only: int2e_ip1_sph, int2e_ip1_cart
   use cint_gen_hess, only: int2e_ipip1_sph, int2e_ipip1_cart, &
                            int2e_ipvip1_sph, int2e_ipvip1_cart, &
                            int2e_ip1ip2_sph, int2e_ip1ip2_cart
   implicit none
   private

   public :: HESS_OVLP_II, HESS_OVLP_IJ, HESS_KIN_II, HESS_KIN_IJ
   public :: HESS_NUC_II, HESS_NUC_IJ
   public :: hess_1e_block
   public :: HESS_RINV_II, HESS_RINV_IJ
   public :: hess_rinv_block
   public :: HESS_ERI_II, HESS_ERI_IJ, HESS_ERI_IK
   public :: hess_2e_block
   public :: hess_2e_contract
   public :: eri_ip1_block
   public :: h1_contract

   integer, parameter :: HESS_OVLP_II = 1   !! `int1e_ipipovlp`, both derivatives on the bra
   integer, parameter :: HESS_OVLP_IJ = 2   !! `int1e_ipovlpip`, one on each centre
   integer, parameter :: HESS_KIN_II = 3
   integer, parameter :: HESS_KIN_IJ = 4
   integer, parameter :: HESS_NUC_II = 5
   integer, parameter :: HESS_NUC_IJ = 6

   integer, parameter :: HESS_RINV_II = 1   !! `int1e_ipiprinv`, both derivatives on the bra
   integer, parameter :: HESS_RINV_IJ = 2   !! `int1e_iprinvip`, one on each centre

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

      ! **The dispatch below ends in a `case default`**, so an unrecognised
      ! selector would silently return whichever integral that branch names --
      ! a real integral, of the right shape, and the wrong one. Checked here
      ! instead of letting the fallthrough decide.
      if (which < HESS_OVLP_II .or. which > HESS_NUC_IJ) then
         call error%set(ERROR_VALIDATION, "unknown one-electron second-derivative "// &
                        "selector; expected one of the HESS_* constants.")
         return
      end if

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

   subroutine hess_rinv_block(mol, iatom, which, matrix, error)
      !! A second derivative of `1/|r-R|` with the operator origin on one atom
      !!
      !! The nuclear attraction second derivative splits three ways, and only
      !! one of the three is `int1e_ipipnuc`. Moving atom `A` moves the basis
      !! functions centred on it -- that is the `ipipnuc`/`ipnucip` part, which
      !! carries the sum over **all** nuclei and their charges. It also moves
      !! the nucleus itself, and the electrons feel that through the origin of
      !! `1/|r-R_A|` alone. Those cross and diagonal terms need the operator
      !! pinned to one atom at a time, which is what `rinv` is for.
      !!
      !! Returned **unscaled**, exactly as `iprinv_deriv_at` returns the first
      !! derivative: the charge and the sign of the electron-nucleus attraction
      !! are the caller's, because the caller is also the one that knows which
      !! of the three terms it is assembling.
      !!
      !! `env` is copied so the origin can be moved without touching a molecule
      !! that other threads share.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: iatom
      integer, intent(in) :: which          !! `HESS_RINV_II` or `HESS_RINV_IJ`
      real(dp), allocatable, intent(out) :: matrix(:, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: buf(:), env(:)
      integer, allocatable :: atm_flat(:), bas_flat(:)
      integer :: dims(0:3), shls(0:1)
      integer :: ish, jsh, di, dj, io, jo, i, j, comp, mx
      logical :: have

      if (error%has_error()) return

      if (which /= HESS_RINV_II .and. which /= HESS_RINV_IJ) then
         call error%set(ERROR_VALIDATION, "unknown rinv second-derivative "// &
                        "selector; expected HESS_RINV_II or HESS_RINV_IJ.")
         return
      end if
      if (iatom < 1 .or. iatom > mol%natm) then
         call error%set(ERROR_VALIDATION, "rinv origin atom is outside the molecule.")
         return
      end if

      mx = 0
      do ish = 1, mol%nbas
         mx = max(mx, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do

      allocate (buf(mx*mx*N_COMPONENTS))
      allocate (matrix(mol%nao, mol%nao, N_COMPONENTS))
      matrix = 0.0_dp
      atm_flat = reshape(mol%atm, [size(mol%atm)])
      bas_flat = reshape(mol%bas, [size(mol%bas)])
      allocate (env(size(mol%env)), source=mol%env)
      env(PTR_RINV_ORIG + 1:PTR_RINV_ORIG + 3) = mol%coords(:, iatom)

      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, mol%nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]
            dims = [di, dj, 1, 1]

            if (mol%cartesian) then
               if (which == HESS_RINV_II) then
                  have = int1e_ipiprinv_cart(buf, dims, shls, atm_flat, mol%natm, &
                                             bas_flat, mol%nbas, env, ws)
               else
                  have = int1e_iprinvip_cart(buf, dims, shls, atm_flat, mol%natm, &
                                             bas_flat, mol%nbas, env, ws)
               end if
            else
               if (which == HESS_RINV_II) then
                  have = int1e_ipiprinv_sph(buf, dims, shls, atm_flat, mol%natm, &
                                            bas_flat, mol%nbas, env, ws)
               else
                  have = int1e_iprinvip_sph(buf, dims, shls, atm_flat, mol%natm, &
                                            bas_flat, mol%nbas, env, ws)
               end if
            end if
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

      deallocate (buf, env, atm_flat, bas_flat)
   end subroutine hess_rinv_block

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

      ! Same reason as the one-electron driver: the dispatch has a
      ! `case default` and an unrecognised selector would come back as a
      ! plausible wrong integral rather than as an error.
      if (which < HESS_ERI_II .or. which > HESS_ERI_IK) then
         call error%set(ERROR_VALIDATION, "unknown two-electron second-derivative "// &
                        "selector; expected one of the HESS_ERI_* constants.")
         return
      end if

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

   subroutine eri_ip1_block(mol, eri, error)
      !! `int2e_ip1` over the whole basis, as `(n_ao, n_ao, n_ao, n_ao, 3)`
      !!
      !! The first derivative with the nabla on the bra's first index. The
      !! gradient contracts this on the fly and never materialises it, which is
      !! the right thing there; the Hessian's per-atom perturbation needs the
      !! same integrals summed over four different index positions per atom, so
      !! having the array is what makes that assembly readable.
      !!
      !! Same `n_ao^4` caveat as the second derivatives next door: correct
      !! first, contracted on the fly once it matters.
      type(libcint_molecule_t), intent(in) :: mol
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

      allocate (buf(mx**4*3))
      allocate (eri(mol%nao, mol%nao, mol%nao, mol%nao, 3))
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

                  if (mol%cartesian) then
                     have = int2e_ip1_cart(buf, dims, shls, atm_flat, mol%natm, &
                                           bas_flat, mol%nbas, mol%env, ws)
                  else
                     have = int2e_ip1_sph(buf, dims, shls, atm_flat, mol%natm, &
                                          bas_flat, mol%nbas, mol%env, ws)
                  end if
                  if (.not. have) cycle

                  do comp = 1, 3
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
   end subroutine eri_ip1_block

   subroutine hess_2e_contract(mol, density, hess, error)
      !! The two-electron second derivatives, contracted as they are computed
      !!
      !! Same numbers `hess_2e_block` produces and never the same array. Each
      !! shell quartet is evaluated into a buffer, deposited into the Hessian,
      !! and forgotten, so the memory is one quartet rather than `nao^4` times
      !! nine times three. That is the difference between a routine that works
      !! on water and one that works: at 24 basis functions the arrays are
      !! 72 MB, at 50 they are 1.3 GB, and the growth is what it looks like.
      !!
      !! **The deposit rules are the same three the assembly derives**, with
      !! the same 4, 4 and 8 from folding the sixteen index pairs against a
      !! weight symmetric under all eight permutations of `(mn|ls)`. The loop
      !! runs over every quartet in every order, which is what makes that
      !! folding valid -- restricting it to unique quartets would need the
      !! *derivative* integrals to carry symmetries they do not have.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(inout) :: hess(:, :, :, :)   !! (3, 3, natm, natm), accumulated
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: buf_ii(:), buf_ij(:), buf_ik(:)
      real(dp), allocatable :: hloc(:, :, :, :)
      integer, allocatable :: atm_flat(:), bas_flat(:), owner(:)
      integer, allocatable :: offsets(:), counts(:), sh_dim(:), sh_off(:)
      integer :: dims(0:3), shls(0:3)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: io, jo, ko, lo, i, j, k, l, comp, mx, idx
      integer :: nao, natm, ia, ja, ka, a, b, c, pair, n_pairs
      real(dp) :: gam
      logical :: have_ii, have_ij, have_ik

      if (error%has_error()) return

      nao = mol%nao
      natm = mol%natm

      allocate (offsets(natm), counts(natm), owner(nao))
      call atom_ao_blocks(mol, offsets, counts)
      owner = 0
      do c = 1, natm
         do i = offsets(c) + 1, offsets(c) + counts(c)
            owner(i) = c
         end do
      end do

      ! Hoisted out of the quartet loop: `shell_dim` reads `bas` and the
      ! offsets walk it, and both are the same answer every one of `nbas^4`
      ! times round.
      allocate (sh_dim(mol%nbas), sh_off(mol%nbas))
      mx = 0
      do ish = 1, mol%nbas
         sh_dim(ish) = shell_dim(mol%cartesian, ish - 1, mol%bas)
         sh_off(ish) = mol%shell_offset(ish)
         mx = max(mx, sh_dim(ish))
      end do

      atm_flat = reshape(mol%atm, [size(mol%atm)])
      bas_flat = reshape(mol%bas, [size(mol%bas)])
      n_pairs = mol%nbas*mol%nbas

      ! One accumulator per thread, merged once at the end. The Hessian is
      ! `9*natm^2` doubles -- small enough that a private copy costs nothing and
      ! an atomic update per quartet would cost everything.
      !$omp parallel default(none) &
      !$omp shared(mol, density, hess, owner, sh_dim, sh_off, atm_flat, bas_flat, &
      !$omp        mx, natm, n_pairs, error) &
      !$omp private(buf_ii, buf_ij, buf_ik, hloc, dims, shls, ish, jsh, ksh, lsh, &
      !$omp         di, dj, dk, dl, io, jo, ko, lo, i, j, k, l, comp, idx, &
      !$omp         ia, ja, ka, a, b, gam, pair, have_ii, have_ij, have_ik)
      allocate (buf_ii(mx**4*N_COMPONENTS), buf_ij(mx**4*N_COMPONENTS), &
                buf_ik(mx**4*N_COMPONENTS))
      allocate (hloc(3, 3, natm, natm))
      hloc = 0.0_dp

      !$omp do schedule(dynamic)
      do pair = 1, n_pairs
         ish = (pair - 1)/mol%nbas + 1
         jsh = pair - mol%nbas*(ish - 1)
         di = sh_dim(ish)
         io = sh_off(ish)
         dj = sh_dim(jsh)
         jo = sh_off(jsh)
         do ksh = 1, mol%nbas
            dk = sh_dim(ksh)
            ko = sh_off(ksh)
            do lsh = 1, mol%nbas
               dl = sh_dim(lsh)
               lo = sh_off(lsh)
               shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]
               dims = [di, dj, dk, dl]

               have_ii = drive_2e(mol, HESS_ERI_II, buf_ii, dims, shls, atm_flat, bas_flat)
               have_ij = drive_2e(mol, HESS_ERI_IJ, buf_ij, dims, shls, atm_flat, bas_flat)
               have_ik = drive_2e(mol, HESS_ERI_IK, buf_ik, dims, shls, atm_flat, bas_flat)
               if (.not. (have_ii .or. have_ij .or. have_ik)) cycle

               do l = 1, dl
                  do k = 1, dk
                     ka = owner(ko + k)
                     do j = 1, dj
                        ja = owner(jo + j)
                        do i = 1, di
                           ia = owner(io + i)
                           gam = 0.5_dp*density(io + i, jo + j)*density(ko + k, lo + l) &
                                 - 0.125_dp*(density(io + i, lo + l)*density(ko + k, jo + j) &
                                          + density(io + i, ko + k)*density(jo + j, lo + l))
                           if (gam == 0.0_dp) cycle
                           do comp = 1, N_COMPONENTS
                              a = (comp - 1)/3 + 1
                              b = comp - 3*(a - 1)
                              idx = i + di*(j - 1 + dj*(k - 1 + dk*(l - 1 + dl*(comp - 1))))
                              hloc(a, b, ia, ia) = hloc(a, b, ia, ia) + 4.0_dp*gam*buf_ii(idx)
                              hloc(a, b, ia, ja) = hloc(a, b, ia, ja) + 4.0_dp*gam*buf_ij(idx)
                              hloc(a, b, ia, ka) = hloc(a, b, ia, ka) + 8.0_dp*gam*buf_ik(idx)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end do

      !$omp critical
      hess = hess + hloc
      !$omp end critical
      deallocate (buf_ii, buf_ij, buf_ik, hloc)
      !$omp end parallel

      deallocate (atm_flat, bas_flat, owner, offsets, counts, sh_dim, sh_off)
   end subroutine hess_2e_contract

   subroutine h1_contract(mol, density, h1, error)
      !! The skeleton derivative Fock for **every** atom, in one pass
      !!
      !! `make_h1_atom` builds one atom's `dF/dR` from a stored `int2e_ip1`
      !! array, which costs `nao^4` of memory and another `nao^4` of work per
      !! atom. This walks the shell quartets once and deposits into whichever
      !! atom each quartet's differentiated index belongs to, so the `natm`
      !! disappears from the work and the array disappears entirely.
      !!
      !! **The eight terms become seven deposits, and that is not a lost one.**
      !! `make_h1_atom` names each term by permuting an index into first place
      !! and indexing the stored array accordingly. Here the permuted orderings
      !! are *other quartets*, which the loop visits anyway, so each term is
      !! rewritten in the one ordering the buffer has. Two of the Coulomb terms
      !! land on the same element with weights `d(i,j)` and `d(j,i)`, which are
      !! equal, so they are written once with a two.
      !!
      !! Checked against `make_h1_atom` element by element in the test suite --
      !! that routine stays as the readable statement of what this computes.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: h1(:, :, :, :)   !! (nao, nao, 3, natm)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: buf(:), hloc(:, :, :, :)
      integer, allocatable :: atm_flat(:), bas_flat(:), owner(:)
      integer, allocatable :: offsets(:), counts(:), sh_dim(:), sh_off(:)
      integer :: dims(0:3), shls(0:3)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: io, jo, ko, lo, i, j, k, l, comp, mx, idx
      integer :: nao, natm, ii, jj, kk, ll, at, c, pair, n_pairs
      real(dp) :: b
      logical :: have

      if (error%has_error()) return

      nao = mol%nao
      natm = mol%natm

      allocate (offsets(natm), counts(natm), owner(nao))
      call atom_ao_blocks(mol, offsets, counts)
      owner = 0
      do c = 1, natm
         do i = offsets(c) + 1, offsets(c) + counts(c)
            owner(i) = c
         end do
      end do

      allocate (sh_dim(mol%nbas), sh_off(mol%nbas))
      mx = 0
      do ish = 1, mol%nbas
         sh_dim(ish) = shell_dim(mol%cartesian, ish - 1, mol%bas)
         sh_off(ish) = mol%shell_offset(ish)
         mx = max(mx, sh_dim(ish))
      end do

      atm_flat = reshape(mol%atm, [size(mol%atm)])
      bas_flat = reshape(mol%bas, [size(mol%bas)])
      n_pairs = mol%nbas*mol%nbas

      allocate (h1(nao, nao, 3, natm))
      h1 = 0.0_dp

      ! A private accumulator per thread is `nao^2 * 3 * natm` doubles, which
      ! is the same size as the answer. That is affordable at the sizes this
      ! runs at and is the thing to revisit first if it stops being: the
      ! alternative is an atomic per deposit, and there are seven per buffer
      ! element.
      !$omp parallel default(none) &
      !$omp shared(mol, density, h1, owner, sh_dim, sh_off, atm_flat, bas_flat, &
      !$omp        mx, nao, natm, n_pairs) &
      !$omp private(buf, hloc, dims, shls, ish, jsh, ksh, lsh, di, dj, dk, dl, &
      !$omp         io, jo, ko, lo, i, j, k, l, comp, idx, ii, jj, kk, ll, at, b, &
      !$omp         pair, have)
      allocate (buf(mx**4*3))
      allocate (hloc(nao, nao, 3, natm))
      hloc = 0.0_dp

      !$omp do schedule(dynamic)
      do pair = 1, n_pairs
         ish = (pair - 1)/mol%nbas + 1
         jsh = pair - mol%nbas*(ish - 1)
         di = sh_dim(ish)
         io = sh_off(ish)
         dj = sh_dim(jsh)
         jo = sh_off(jsh)
         do ksh = 1, mol%nbas
            dk = sh_dim(ksh)
            ko = sh_off(ksh)
            do lsh = 1, mol%nbas
               dl = sh_dim(lsh)
               lo = sh_off(lsh)
               shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]
               dims = [di, dj, dk, dl]

               if (mol%cartesian) then
                  have = int2e_ip1_cart(buf, dims, shls, atm_flat, mol%natm, &
                                        bas_flat, mol%nbas, mol%env, ws)
               else
                  have = int2e_ip1_sph(buf, dims, shls, atm_flat, mol%natm, &
                                       bas_flat, mol%nbas, mol%env, ws)
               end if
               if (.not. have) cycle

               do comp = 1, 3
                  do l = 1, dl
                     ll = lo + l
                     do k = 1, dk
                        kk = ko + k
                        do j = 1, dj
                           jj = jo + j
                           do i = 1, di
                              ii = io + i
                              idx = i + di*(j - 1 + dj*(k - 1 + dk*(l - 1 + dl*(comp - 1))))
                              b = buf(idx)
                              if (b == 0.0_dp) cycle
                              ! The nabla sits on the first index, so the whole
                              ! quartet belongs to that index's atom.
                              at = owner(ii)
                              hloc(ii, jj, comp, at) = hloc(ii, jj, comp, at) - density(kk, ll)*b
                              hloc(jj, ii, comp, at) = hloc(jj, ii, comp, at) - density(kk, ll)*b
                              hloc(kk, ll, comp, at) = hloc(kk, ll, comp, at) - 2.0_dp*density(ii, jj)*b
                              hloc(ii, ll, comp, at) = hloc(ii, ll, comp, at) + 0.5_dp*density(jj, kk)*b
                              hloc(jj, ll, comp, at) = hloc(jj, ll, comp, at) + 0.5_dp*density(ii, kk)*b
                              hloc(kk, jj, comp, at) = hloc(kk, jj, comp, at) + 0.5_dp*density(ll, ii)*b
                              hloc(kk, ii, comp, at) = hloc(kk, ii, comp, at) + 0.5_dp*density(ll, jj)*b
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end do

      !$omp critical
      h1 = h1 + hloc
      !$omp end critical
      deallocate (buf, hloc)
      !$omp end parallel

      deallocate (atm_flat, bas_flat, owner, offsets, counts, sh_dim, sh_off)
   end subroutine h1_contract

end module mqc_libcint_hess_ints
