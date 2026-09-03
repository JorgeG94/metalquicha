!! Second derivatives of the one- and two-electron integrals with respect to nuclei
module mqc_libcint_hess_ints
   !! What an analytic Hessian needs from the integral library, and nothing else.
   !!
   !! **These come from libfint's native modules rather than through the
   !! `libcint_fortran` compatibility layer**, which exposes first derivatives
   !! and below only. So analytic Hessians are available in a libfint build --
   !! the default -- and a libcint build keeps the finite-difference path.
   !!
   !! **Both orderings of every derivative are needed, and they are different
   !! integrals.** `ipipovlp` puts both derivatives on the bra centre;
   !! `ipovlpip` puts one on each. A Hessian built from only the first is not a
   !! slightly worse Hessian, it is one that violates translational invariance,
   !! which shows up as translations that are no longer zero-frequency modes.
   use, intrinsic :: iso_fortran_env, only: real64
   use, intrinsic :: iso_c_binding, only: c_null_ptr
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim, atom_ao_blocks, &
                                    eri_shell_table_t, eri_shell_table, eri_schwarz_collapse
   use mqc_libcint_direct, only: schwarz_bounds
   use libcint_fortran, only: LIBCINT_NPRIM_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_RINV_ORIG, &
                              LIBCINT_PTR_RANGE_OMEGA
   use mqc_libcint_hess_abi, only: &
      cint1e_ipipovlp_sph, cint1e_ipipovlp_cart, cint1e_ipovlpip_sph, cint1e_ipovlpip_cart, &
      cint1e_ipipkin_sph, cint1e_ipipkin_cart, cint1e_ipkinip_sph, cint1e_ipkinip_cart, &
      cint1e_ipipnuc_sph, cint1e_ipipnuc_cart, cint1e_ipnucip_sph, cint1e_ipnucip_cart, &
      cint1e_ipiprinv_sph, cint1e_ipiprinv_cart, cint1e_iprinvip_sph, cint1e_iprinvip_cart, &
      cint2e_ipip1_sph, cint2e_ipip1_cart, cint2e_ipvip1_sph, cint2e_ipvip1_cart, &
      cint2e_ip1ip2_sph, cint2e_ip1ip2_cart, cint2e_ip1_sph, cint2e_ip1_cart
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
   public :: hess_2e_skeleton_contract
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

   real(dp), parameter :: HESS_SCREEN_TOL = 1.0e-10_dp
      !! Estimated contribution below which a differentiated quartet is dropped.
      !!
      !! The same number the SCF gradient screens its differentiated quartets
      !! at. The bound is an *estimate*, not a rigorous Cauchy-Schwarz one:
      !! differentiating brings a factor of `2 alpha r` down from the
      !! exponential, and what stands in for it is `2 sqrt(alpha_max)` per
      !! differentiation, twice over for a second derivative.

   integer, parameter :: N_COMPONENTS = 9
      !! Every one of these is a 3x3 Cartesian block per shell pair, laid out by
      !! libcint as `xx xy xz yx yy yz zx zy zz` in the slowest index.

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
      ! compatibility layer took them shaped. Flattened once rather than per
      ! shell pair.
      integer, allocatable :: atm_flat(:), bas_flat(:)
      integer :: dims(0:3), shls(0:1)
      integer :: ish, jsh, di, dj, io, jo, i, j, comp, mx
      logical :: have

      if (error%has_error()) return

      ! **The dispatch below ends in a `case default`**, so an unrecognised
      ! selector would silently return a real integral of the right shape and
      ! the wrong kind. Checked here instead.
      if (which < HESS_OVLP_II .or. which > HESS_NUC_IJ) then
         call error%set(ERROR_VALIDATION, "unknown one-electron second-derivative "// &
                        "selector; expected one of the HESS_* constants.")
         return
      end if

      mx = 0
      do ish = 1, mol%nbas
         mx = max(mx, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do

      ! Flat and indexed by hand: the library packs a block with the *actual*
      ! shell dimensions, so a rank-3 buffer declared at the largest shell has
      ! the wrong strides for every smaller one.
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
            have = cint1e_ipipovlp_cart(buf, shls, atm, mol%natm, &
                                        bas, mol%nbas, mol%env) /= 0
         case (HESS_OVLP_IJ)
            have = cint1e_ipovlpip_cart(buf, shls, atm, mol%natm, &
                                        bas, mol%nbas, mol%env) /= 0
         case (HESS_KIN_II)
            have = cint1e_ipipkin_cart(buf, shls, atm, mol%natm, &
                                       bas, mol%nbas, mol%env) /= 0
         case (HESS_KIN_IJ)
            have = cint1e_ipkinip_cart(buf, shls, atm, mol%natm, &
                                       bas, mol%nbas, mol%env) /= 0
         case (HESS_NUC_II)
            have = cint1e_ipipnuc_cart(buf, shls, atm, mol%natm, &
                                       bas, mol%nbas, mol%env) /= 0
         case default
            have = cint1e_ipnucip_cart(buf, shls, atm, mol%natm, &
                                       bas, mol%nbas, mol%env) /= 0
         end select
      else
         select case (which)
         case (HESS_OVLP_II)
            have = cint1e_ipipovlp_sph(buf, shls, atm, mol%natm, &
                                       bas, mol%nbas, mol%env) /= 0
         case (HESS_OVLP_IJ)
            have = cint1e_ipovlpip_sph(buf, shls, atm, mol%natm, &
                                       bas, mol%nbas, mol%env) /= 0
         case (HESS_KIN_II)
            have = cint1e_ipipkin_sph(buf, shls, atm, mol%natm, &
                                      bas, mol%nbas, mol%env) /= 0
         case (HESS_KIN_IJ)
            have = cint1e_ipkinip_sph(buf, shls, atm, mol%natm, &
                                      bas, mol%nbas, mol%env) /= 0
         case (HESS_NUC_II)
            have = cint1e_ipipnuc_sph(buf, shls, atm, mol%natm, &
                                      bas, mol%nbas, mol%env) /= 0
         case default
            have = cint1e_ipnucip_sph(buf, shls, atm, mol%natm, &
                                      bas, mol%nbas, mol%env) /= 0
         end select
      end if
   end function drive

   subroutine hess_rinv_block(mol, iatom, which, matrix, error)
      !! A second derivative of `1/|r-R|` with the operator origin on one atom
      !!
      !! The nuclear attraction second derivative splits three ways, and only
      !! one of the three is `int1e_ipipnuc`. Moving atom `A` moves the basis
      !! functions centred on it -- the `ipipnuc`/`ipnucip` part, which carries
      !! the sum over **all** nuclei and their charges -- and it moves the
      !! nucleus itself, which the electrons feel through the origin of
      !! `1/|r-R_A|` alone. Those cross and diagonal terms need the operator
      !! pinned to one atom at a time, which is what `rinv` is for.
      !!
      !! Returned **unscaled**, as `iprinv_deriv_at` returns the first
      !! derivative: the charge and the sign of the electron-nucleus attraction
      !! are the caller's.
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
      env(LIBCINT_PTR_RINV_ORIG + 1:LIBCINT_PTR_RINV_ORIG + 3) = mol%coords(:, iatom)

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
                  have = cint1e_ipiprinv_cart(buf, shls, atm_flat, mol%natm, &
                                              bas_flat, mol%nbas, env) /= 0
               else
                  have = cint1e_iprinvip_cart(buf, shls, atm_flat, mol%natm, &
                                              bas_flat, mol%nbas, env) /= 0
               end if
            else
               if (which == HESS_RINV_II) then
                  have = cint1e_ipiprinv_sph(buf, shls, atm_flat, mol%natm, &
                                             bas_flat, mol%nbas, env) /= 0
               else
                  have = cint1e_iprinvip_sph(buf, shls, atm_flat, mol%natm, &
                                             bas_flat, mol%nbas, env) /= 0
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
      !! over. `hess_2e_contract` is the same integrals contracted as they are
      !! computed, and is what a Hessian of any size uses.
      !!
      !! **The three selectors are three genuinely different integrals**, not
      !! one integral indexed three ways. `ipip1` puts both derivatives on the
      !! first centre, `ipvip1` puts one each on the bra pair, and `ip1ip2`
      !! puts one on the bra and one on the ket. Translational invariance
      !! relates them only if all three are present.
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

      ! As in the one-electron driver: the dispatch has a `case default`, so an
      ! unrecognised selector would come back as a plausible wrong integral.
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

                  have = drive_2e(mol, which, buf, dims, shls, atm_flat, bas_flat, mol%env)
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

   function drive_2e(mol, which, buf, dims, shls, atm, bas, env) result(have)
      !! Dispatch one shell quartet to the right entry point
      !!
      !! `env` is passed rather than taken from `mol` so a caller can hand in a
      !! copy with `PTR_RANGE_OMEGA` set. Range separation in libcint is an
      !! environment slot, not a different procedure.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: which
      real(dp), intent(inout) :: buf(0:)
      integer, intent(in) :: dims(0:), shls(0:)
      integer, intent(in), target :: atm(0:), bas(0:)
      real(dp), intent(in), target :: env(0:)
      logical :: have

      have = .false.
      if (mol%cartesian) then
         select case (which)
         case (HESS_ERI_II)
            have = cint2e_ipip1_cart(buf, shls, atm, mol%natm, &
                                     bas, mol%nbas, env, c_null_ptr) /= 0
         case (HESS_ERI_IJ)
            have = cint2e_ipvip1_cart(buf, shls, atm, mol%natm, &
                                      bas, mol%nbas, env, c_null_ptr) /= 0
         case default
            have = cint2e_ip1ip2_cart(buf, shls, atm, mol%natm, &
                                      bas, mol%nbas, env, c_null_ptr) /= 0
         end select
      else
         select case (which)
         case (HESS_ERI_II)
            have = cint2e_ipip1_sph(buf, shls, atm, mol%natm, &
                                    bas, mol%nbas, env, c_null_ptr) /= 0
         case (HESS_ERI_IJ)
            have = cint2e_ipvip1_sph(buf, shls, atm, mol%natm, &
                                     bas, mol%nbas, env, c_null_ptr) /= 0
         case default
            have = cint2e_ip1ip2_sph(buf, shls, atm, mol%natm, &
                                     bas, mol%nbas, env, c_null_ptr) /= 0
         end select
      end if
   end function drive_2e

   subroutine eri_ip1_block(mol, eri, error)
      !! `int2e_ip1` over the whole basis, as `(n_ao, n_ao, n_ao, n_ao, 3)`
      !!
      !! The first derivative with the nabla on the bra's first index. The
      !! gradient contracts this on the fly and never materialises it; the
      !! Hessian's per-atom perturbation needs the same integrals summed over
      !! four index positions per atom, which the array makes readable at the
      !! cost of `n_ao^4`. `h1_contract` is the contracted form.
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
                     have = cint2e_ip1_cart(buf, shls, atm_flat, mol%natm, &
                                            bas_flat, mol%nbas, mol%env, c_null_ptr) /= 0
                  else
                     have = cint2e_ip1_sph(buf, shls, atm_flat, mol%natm, &
                                           bas_flat, mol%nbas, mol%env, c_null_ptr) /= 0
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

   subroutine hess_2e_contract(mol, density, hess, error, screen_tol, k_scale, &
                               j_scale, omega)
      !! The two-electron second derivatives, contracted as they are computed
      !!
      !! Same numbers `hess_2e_block` produces and never the same array. Each
      !! shell quartet is evaluated into a buffer, deposited into the Hessian,
      !! and forgotten, so the memory is one quartet rather than `nao^4` times
      !! nine times three.
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
      real(dp), intent(in), optional :: screen_tol
      real(dp), intent(in), optional :: k_scale
         !! Fraction of exact exchange. One is Hartree-Fock and the default;
         !! zero is a pure density functional, which has none; a hybrid is its
         !! mixing fraction.
      real(dp), intent(in), optional :: j_scale
         !! Fraction of Coulomb, one by default. Zero is what a long-range pass
         !! wants: a range-separated functional's second exchange term has no
         !! Coulomb of its own to add, the full-range pass having already
         !! supplied it.
      real(dp), intent(in), optional :: omega
         !! Range-separation parameter. Reaches libcint through
         !! `env(PTR_RANGE_OMEGA)`, so the attenuated integrals are the same
         !! entry points over the same quartets -- no new procedures, which is
         !! what makes this a threading job rather than an integral one.
      real(dp), allocatable :: buf_ii(:), buf_ij(:), buf_ik(:), buf_ik2(:)
      real(dp), allocatable :: hloc(:, :, :, :)
      integer, allocatable :: atm_flat(:), bas_flat(:), owner(:)
      integer, allocatable :: offsets(:), counts(:), sh_dim(:), sh_off(:)
      real(dp), allocatable :: bounds(:, :), dsh(:, :), sa(:)
      real(dp) :: qq, wbound, est, tol, kx, jx
      real(dp), allocatable, target :: env_use(:)
      integer :: dims(0:3), shls(0:3), dims2(0:3), shls2(0:3)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: io, jo, ko, lo, i, j, k, l, comp, mx, idx
      integer :: nao, natm, ia, ja, ka, la, a, b, c, pair, n_pairs
      real(dp) :: gam, ket_w
      integer :: idx2
      logical :: have_ii, have_ij, have_ik, have_ik2

      if (error%has_error()) return

      kx = 1.0_dp
      if (present(k_scale)) kx = k_scale
      jx = 1.0_dp
      if (present(j_scale)) jx = j_scale
      ! A local copy so an omega pass can set the range-separation slot
      ! without disturbing the molecule every other caller shares.
      allocate (env_use(0:size(mol%env) - 1), source=mol%env)
      if (present(omega)) env_use(LIBCINT_PTR_RANGE_OMEGA) = omega

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

      ! Hoisted out of the quartet loop: both are the same answer every one of
      ! `nbas^4` times round.
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

      ! ---- what a quartet can contribute, before computing it ---------------
      !
      ! Two factors, and only one of them is rigorous. `bounds` is the exact
      ! Schwarz bound on the undifferentiated quartet. `sa` stands in for what
      ! differentiation does -- a factor of `2 alpha r` from the exponential,
      ! about `2 sqrt(alpha_max)` per derivative and taken twice here -- which
      ! is why `HESS_SCREEN_TOL` carries margin.
      !
      ! `dsh` bounds the contraction weight: a quartet whose integrals are
      ! large contributes nothing if every density element it multiplies is
      ! negligible.
      tol = HESS_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol
      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) return

      allocate (dsh(mol%nbas, mol%nbas), sa(mol%nbas))
      do ish = 1, mol%nbas
         do jsh = 1, mol%nbas
            dsh(ish, jsh) = maxval(abs(density(sh_off(ish) + 1:sh_off(ish) + sh_dim(ish), &
                                               sh_off(jsh) + 1:sh_off(jsh) + sh_dim(jsh))))
         end do
         sa(ish) = sqrt(maxval(mol%env(mol%bas(LIBCINT_PTR_EXP, ish) + 1: &
                                       mol%bas(LIBCINT_PTR_EXP, ish) &
                                       + mol%bas(LIBCINT_NPRIM_OF, ish))))
      end do

      ! One accumulator per thread, merged once at the end. The Hessian is
      ! `9*natm^2` doubles, so a private copy costs nothing where an atomic
      ! update per quartet would cost everything.
      !$omp parallel default(none) &
      !$omp shared(kx, jx, env_use, mol, density, hess, owner, sh_dim, sh_off, atm_flat, bas_flat, &
      !$omp        mx, natm, n_pairs, error, bounds, dsh, sa, tol) &
      !$omp private(buf_ii, buf_ij, buf_ik, buf_ik2, hloc, dims, shls, dims2, shls2, &
      !$omp         ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo, i, j, k, l, &
      !$omp         comp, idx, idx2, ia, ja, ka, la, a, b, gam, ket_w, pair, &
      !$omp         have_ii, have_ij, have_ik, have_ik2, qq, wbound, est)
      allocate (buf_ii(mx**4*N_COMPONENTS), buf_ij(mx**4*N_COMPONENTS), &
                buf_ik(mx**4*N_COMPONENTS), buf_ik2(mx**4*N_COMPONENTS))
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
            ! **Only the ket pair is restricted, and only two of the three
            ! integrals get to use it.** `ipip1` and `ipvip1` leave the ket
            ! untouched, so `(kl)` and `(lk)` are the same integral, carry the
            ! same weight and deposit into the same atom pair; one call covers
            ! both orderings.
            !
            ! `ip1ip2` does not get that: its second derivative sits on the
            ! *first* ket index, so swapping `k` and `l` moves the derivative
            ! to a different atom. Both orderings are computed.
            !
            ! No permutation of the bra is available at all, `ipip1`
            ! distinguishing its two bra indices and `ip1ip2` distinguishing
            ! bra from ket.
            do lsh = 1, ksh
               dl = sh_dim(lsh)
               lo = sh_off(lsh)
               shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]
               dims = [di, dj, dk, dl]
               ket_w = 2.0_dp
               if (lsh == ksh) ket_w = 1.0_dp

               ! The largest weight any element of this quartet can carry,
               ! written the way the contraction weight is: a Coulomb term and
               ! the two exchange terms the eightfold symmetrisation produces.
               wbound = 0.5_dp*dsh(ish, jsh)*dsh(ksh, lsh) &
                        + abs(kx)*0.125_dp*(dsh(ish, lsh)*dsh(ksh, jsh) &
                                            + dsh(ish, ksh)*dsh(jsh, lsh))
               qq = bounds(ish, jsh)*bounds(ksh, lsh)
               ! Four times the deposit weights (4, 4 and 8) against the
               ! per-integral derivative factors, which differ because the two
               ! derivatives sit on different centres in each.
               est = 16.0_dp*wbound*qq*(sa(ish)*sa(ish) &
                                        + sa(ish)*sa(jsh) &
                                        + 2.0_dp*sa(ish)*sa(ksh))
               if (est < tol) cycle

               have_ii = drive_2e(mol, HESS_ERI_II, buf_ii, dims, shls, atm_flat, bas_flat, env_use)
               have_ij = drive_2e(mol, HESS_ERI_IJ, buf_ij, dims, shls, atm_flat, bas_flat, env_use)
               have_ik = drive_2e(mol, HESS_ERI_IK, buf_ik, dims, shls, atm_flat, bas_flat, env_use)
               have_ik2 = .false.
               if (lsh /= ksh) then
                  shls2 = [ish - 1, jsh - 1, lsh - 1, ksh - 1]
                  dims2 = [di, dj, dl, dk]
                  have_ik2 = drive_2e(mol, HESS_ERI_IK, buf_ik2, dims2, shls2, &
                                      atm_flat, bas_flat, env_use)
               end if
               if (.not. (have_ii .or. have_ij .or. have_ik .or. have_ik2)) cycle

               do l = 1, dl
                  la = owner(lo + l)
                  do k = 1, dk
                     ka = owner(ko + k)
                     do j = 1, dj
                        ja = owner(jo + j)
                        do i = 1, di
                           ia = owner(io + i)
                           ! Symmetric under `k <-> l`, which is what lets the
                           ! swapped ordering reuse it.
                           gam = jx*0.5_dp*density(io + i, jo + j)*density(ko + k, lo + l) &
                                 - kx*0.125_dp &
                                 *(density(io + i, lo + l)*density(ko + k, jo + j) &
                                   + density(io + i, ko + k)*density(jo + j, lo + l))
                           if (gam == 0.0_dp) cycle
                           do comp = 1, N_COMPONENTS
                              a = (comp - 1)/3 + 1
                              b = comp - 3*(a - 1)
                              idx = i + di*(j - 1 + dj*(k - 1 + dk*(l - 1 + dl*(comp - 1))))
                              hloc(a, b, ia, ia) = hloc(a, b, ia, ia) &
                                                   + 4.0_dp*ket_w*gam*buf_ii(idx)
                              hloc(a, b, ia, ja) = hloc(a, b, ia, ja) &
                                                   + 4.0_dp*ket_w*gam*buf_ij(idx)
                              hloc(a, b, ia, ka) = hloc(a, b, ia, ka) + 8.0_dp*gam*buf_ik(idx)
                              if (have_ik2) then
                                 idx2 = i + di*(j - 1 + dj*(l - 1 + dl*(k - 1 &
                                                                        + dk*(comp - 1))))
                                 hloc(a, b, ia, la) = hloc(a, b, ia, la) &
                                                      + 8.0_dp*gam*buf_ik2(idx2)
                              end if
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
      deallocate (buf_ii, buf_ij, buf_ik, buf_ik2, hloc)
      !$omp end parallel

      deallocate (atm_flat, bas_flat, owner, offsets, counts, sh_dim, sh_off)
      deallocate (bounds, dsh, sa)
   end subroutine hess_2e_contract

   subroutine hess_2e_skeleton_contract(mol, gamma, dref, owner, hess_corr, &
                                        hess_ref, error, screen_tol, k_scale)
      !! The two-electron second derivatives against a **general** four-index
      !! density, contracted as they are computed
      !!
      !! What `hess_2e_contract` is to a separable density, this is to one that
      !! is not: the MP2 Hessian's effective two-particle density has bra-ket
      !! symmetry and nothing else, so the eightfold folding to weights of 4, 4
      !! and 8 does not apply. The general fold reindexes each of the sixteen
      !! ordered placements so its derivative sits on the slot the integral
      !! differentiates, and the loop over every quartet in every order turns
      !! each into a weight built from the density's permutation orbit --
      !!
      !!     w4 = G(m,n,l,s) + G(n,m,l,s) + G(l,s,m,n) + G(l,s,n,m)
      !!     w8 = w4 + G(m,n,s,l) + G(n,m,s,l) + G(s,l,m,n) + G(s,l,n,m)
      !!
      !! with `w4` multiplying `ipip1` and `ipvip1` and `w8` multiplying
      !! `ip1ip2`. A fully symmetric `G` makes every orbit term equal and the
      !! weights collapse to `4 G` and `8 G`, which is the check that this is
      !! `hess_2e_contract`'s rule and not a fourth one.
      !!
      !! **The SCF reference rides the same sweep.** Both densities deposit
      !! from each buffer -- the reference with the closed-form 4/4/8 weights
      !! -- so the integrals are generated once. `hess_corr` and `hess_ref` are
      !! accumulated into, matching `hess_2e_contract`.
      !!
      !! No ket restriction: `hess_2e_contract` halves its ket loop because a
      !! symmetric weight lets one `ipvip1` call stand for two orderings, and a
      !! general density grants no such thing.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: gamma(:, :, :, :)
         !! Chemist ordered `(mn|ls)` slot layout, as `build_effective_2pdm_ao`
         !! hands it over. Nothing below assumes its bra-ket symmetry.
      real(dp), intent(in) :: dref(:, :)
         !! The reference's separable AO density, factor of two included.
      integer, intent(in) :: owner(:)   !! Owning atom per AO, 1-based
      real(dp), intent(inout) :: hess_corr(:, :, :, :)   !! (3, 3, natm, natm)
      real(dp), intent(inout) :: hess_ref(:, :, :, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: screen_tol
      real(dp), intent(in), optional :: k_scale
         !! Exact-exchange fraction of the reference deposit's separable
         !! weight. Absent is one, Hartree-Fock; a Kohn-Sham reference passes
         !! its functional's fraction so the exchange half of `gref` matches
         !! the operator that produced `dref`. The general-density weights
         !! `w4`/`w8` are untouched: the correlation's two-particle density
         !! carries its own exchange structure, folded where it was built.

      real(dp), allocatable :: buf_ii(:), buf_ij(:), buf_ik(:)
      real(dp), allocatable :: hloc_c(:, :, :, :), hloc_r(:, :, :, :)
      integer, allocatable :: atm_flat(:), bas_flat(:)
      integer, allocatable :: sh_dim(:), sh_off(:)
      real(dp), allocatable :: bounds(:, :), dsh(:, :), gsh(:, :), sa(:)
      real(dp) :: qq, wbound, est, tol, w4, w8, gref, kx
      integer :: dims(0:3), shls(0:3)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: io, jo, ko, lo, i, j, k, l, comp, mx, idx
      integer :: natm, ia, ja, ka, a, b, pair, n_pairs
      integer :: ig, jg, kg, lg
      logical :: have_ii, have_ij, have_ik

      if (error%has_error()) return

      natm = mol%natm

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

      ! Screened as `hess_2e_contract` screens, with the weight bound adapted:
      ! the general density's orbit weights are bounded by per-shell-pair
      ! maxima of `|gamma|` over its leading pair, four block maxima covering
      ! `w4` and twice that covering `w8`.
      tol = HESS_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol
      kx = 1.0_dp
      if (present(k_scale)) kx = k_scale
      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) return

      allocate (dsh(mol%nbas, mol%nbas), gsh(mol%nbas, mol%nbas), sa(mol%nbas))
      do ish = 1, mol%nbas
         do jsh = 1, mol%nbas
            dsh(ish, jsh) = maxval(abs(dref(sh_off(ish) + 1:sh_off(ish) + sh_dim(ish), &
                                            sh_off(jsh) + 1:sh_off(jsh) + sh_dim(jsh))))
            gsh(ish, jsh) = maxval(abs(gamma(sh_off(ish) + 1:sh_off(ish) + sh_dim(ish), &
                                             sh_off(jsh) + 1:sh_off(jsh) + sh_dim(jsh), :, :)))
         end do
         sa(ish) = sqrt(maxval(mol%env(mol%bas(LIBCINT_PTR_EXP, ish) + 1: &
                                       mol%bas(LIBCINT_PTR_EXP, ish) &
                                       + mol%bas(LIBCINT_NPRIM_OF, ish))))
      end do

      !$omp parallel default(none) &
      !$omp shared(mol, gamma, dref, owner, hess_corr, hess_ref, sh_dim, sh_off, &
      !$omp        atm_flat, bas_flat, mx, natm, n_pairs, bounds, dsh, gsh, sa, tol, kx) &
      !$omp private(buf_ii, buf_ij, buf_ik, hloc_c, hloc_r, dims, shls, &
      !$omp         ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo, i, j, k, l, &
      !$omp         comp, idx, ia, ja, ka, a, b, pair, ig, jg, kg, lg, &
      !$omp         w4, w8, gref, have_ii, have_ij, have_ik, qq, wbound, est)
      allocate (buf_ii(mx**4*N_COMPONENTS), buf_ij(mx**4*N_COMPONENTS), &
                buf_ik(mx**4*N_COMPONENTS))
      allocate (hloc_c(3, 3, natm, natm), hloc_r(3, 3, natm, natm))
      hloc_c = 0.0_dp
      hloc_r = 0.0_dp

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

               wbound = 0.5_dp*dsh(ish, jsh)*dsh(ksh, lsh) &
                        + 0.125_dp*(dsh(ish, lsh)*dsh(ksh, jsh) &
                                    + dsh(ish, ksh)*dsh(jsh, lsh)) &
                        + 2.0_dp*(gsh(ish, jsh) + gsh(jsh, ish) &
                                  + gsh(ksh, lsh) + gsh(lsh, ksh))
               qq = bounds(ish, jsh)*bounds(ksh, lsh)
               est = 16.0_dp*wbound*qq*(sa(ish)*sa(ish) &
                                        + sa(ish)*sa(jsh) &
                                        + 2.0_dp*sa(ish)*sa(ksh))
               if (est < tol) cycle

               have_ii = drive_2e(mol, HESS_ERI_II, buf_ii, dims, shls, &
                                  atm_flat, bas_flat, mol%env)
               have_ij = drive_2e(mol, HESS_ERI_IJ, buf_ij, dims, shls, &
                                  atm_flat, bas_flat, mol%env)
               have_ik = drive_2e(mol, HESS_ERI_IK, buf_ik, dims, shls, &
                                  atm_flat, bas_flat, mol%env)
               if (.not. (have_ii .or. have_ij .or. have_ik)) cycle

               do l = 1, dl
                  lg = lo + l
                  do k = 1, dk
                     kg = ko + k
                     ka = owner(kg)
                     do j = 1, dj
                        jg = jo + j
                        ja = owner(jg)
                        do i = 1, di
                           ig = io + i
                           ia = owner(ig)
                           w4 = gamma(ig, jg, kg, lg) + gamma(jg, ig, kg, lg) &
                                + gamma(kg, lg, ig, jg) + gamma(kg, lg, jg, ig)
                           w8 = w4 + gamma(ig, jg, lg, kg) + gamma(jg, ig, lg, kg) &
                                + gamma(lg, kg, ig, jg) + gamma(lg, kg, jg, ig)
                           gref = 0.5_dp*dref(ig, jg)*dref(kg, lg) &
                                  - kx*0.125_dp*(dref(ig, lg)*dref(kg, jg) &
                                                 + dref(ig, kg)*dref(jg, lg))
                           do comp = 1, N_COMPONENTS
                              a = (comp - 1)/3 + 1
                              b = comp - 3*(a - 1)
                              idx = i + di*(j - 1 + dj*(k - 1 + dk*(l - 1 + dl*(comp - 1))))
                              hloc_c(a, b, ia, ia) = hloc_c(a, b, ia, ia) + w4*buf_ii(idx)
                              hloc_c(a, b, ia, ja) = hloc_c(a, b, ia, ja) + w4*buf_ij(idx)
                              hloc_c(a, b, ia, ka) = hloc_c(a, b, ia, ka) + w8*buf_ik(idx)
                              hloc_r(a, b, ia, ia) = hloc_r(a, b, ia, ia) &
                                                     + 4.0_dp*gref*buf_ii(idx)
                              hloc_r(a, b, ia, ja) = hloc_r(a, b, ia, ja) &
                                                     + 4.0_dp*gref*buf_ij(idx)
                              hloc_r(a, b, ia, ka) = hloc_r(a, b, ia, ka) &
                                                     + 8.0_dp*gref*buf_ik(idx)
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
      hess_corr = hess_corr + hloc_c
      hess_ref = hess_ref + hloc_r
      !$omp end critical
      deallocate (buf_ii, buf_ij, buf_ik, hloc_c, hloc_r)
      !$omp end parallel

      deallocate (atm_flat, bas_flat, sh_dim, sh_off, bounds, dsh, gsh, sa)
   end subroutine hess_2e_skeleton_contract

   subroutine h1_contract(mol, density, h1, error, screen_tol, k_scale, &
                          j_scale, omega)
      !! The skeleton derivative Fock for **every** atom, in one pass
      !!
      !! `make_h1_atom` builds one atom's `dF/dR` from a stored `int2e_ip1`
      !! array, at `nao^4` of memory and another `nao^4` of work per atom. This
      !! walks the shell quartets once and deposits into whichever atom each
      !! quartet's differentiated index belongs to, so the `natm` leaves the
      !! work and the array is not needed.
      !!
      !! **The eight terms become seven deposits, and none is lost.** The
      !! permuted orderings `make_h1_atom` indexes for are *other quartets*,
      !! which this loop visits anyway, so each term is rewritten in the one
      !! ordering the buffer has. Two of the Coulomb terms land on the same
      !! element with the equal weights `d(i,j)` and `d(j,i)`, so they are
      !! written once with a two.
      !!
      !! Checked against `make_h1_atom` element by element in the test suite,
      !! that routine staying as the readable statement of what this computes.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: h1(:, :, :, :)   !! (nao, nao, 3, natm)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: screen_tol
      real(dp), intent(in), optional :: k_scale
         !! Fraction of exact exchange. One is Hartree-Fock and the default;
         !! zero is a pure density functional, which has none; a hybrid is its
         !! mixing fraction.
      real(dp), intent(in), optional :: j_scale
         !! Fraction of Coulomb, one by default. Zero is what a long-range pass
         !! wants: a range-separated functional's second exchange term has no
         !! Coulomb of its own to add, the full-range pass having already
         !! supplied it.
      real(dp), intent(in), optional :: omega
         !! Range-separation parameter. Reaches libcint through
         !! `env(PTR_RANGE_OMEGA)`, so the attenuated integrals are the same
         !! entry points over the same quartets -- no new procedures, which is
         !! what makes this a threading job rather than an integral one.
      real(dp), allocatable :: buf(:), hloc(:, :, :, :)
      real(dp), allocatable :: bounds(:, :), bq(:, :), dsh(:, :), sa(:)
      real(dp) :: wmax, est, tol, kx, jx
      real(dp), allocatable, target :: env_use(:)
      integer, allocatable :: atm_flat(:), bas_flat(:), owner(:)
      integer, allocatable :: offsets(:), counts(:)
      type(eri_shell_table_t) :: tab
      integer :: dims(0:3), shls(0:3)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: io, jo, ko, lo, i, j, k, l, comp, mx, idx
      integer :: nao, natm, ii, jj, kk, ll, at, c, pair, n_pairs
      real(dp) :: b
      logical :: have

      if (error%has_error()) return

      kx = 1.0_dp
      if (present(k_scale)) kx = k_scale
      jx = 1.0_dp
      if (present(j_scale)) jx = j_scale
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

      ! **The fused-sp view where the molecule has one.** `int2e_ip1` is one of
      ! the two drivers libfint carries an L shell through. `hess_2e_contract`
      ! above must *not* do this: `ipip1`, `ipvip1` and `ip1ip2` behave as
      ! though L shells did not exist and would read a fused entry as a p shell
      ! over the s coefficients -- not an error, a silently wrong integral.
      !
      ! A fused shell's AO order is its split shells' order, so `owner` and the
      ! density are indexed as before and the deposits need no map.
      call eri_shell_table(mol, tab)
      mx = tab%block_max

      ! A local copy so an omega pass can set the range-separation slot without
      ! disturbing the table every other caller shares. **From `tab%env` and
      ! not `mol%env`**: this loop walks the fused-sp view, whose `bas` points
      ! into its own `env`, so a copy of the molecule's is the wrong array to
      ! hand the integral.
      allocate (env_use(0:size(tab%env) - 1), source=tab%env)
      if (present(omega)) env_use(LIBCINT_PTR_RANGE_OMEGA) = omega

      atm_flat = reshape(mol%atm, [size(mol%atm)])
      bas_flat = reshape(tab%bas, [size(tab%bas)])
      n_pairs = tab%nbas*tab%nbas

      ! Screened the same way `hess_2e_contract` is, with one derivative rather
      ! than two, so the estimate carries one factor of `2 sqrt(alpha_max)`.
      ! The bounds come from the split shells and are re-blocked onto the view,
      ! because that is the table this loop walks.
      tol = HESS_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol
      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) return
      call eri_schwarz_collapse(mol, bounds, bq)

      allocate (dsh(tab%nbas, tab%nbas), sa(tab%nbas))
      do ish = 1, tab%nbas
         do jsh = 1, tab%nbas
            dsh(ish, jsh) = maxval(abs(density(tab%offs(ish) + 1:tab%offs(ish) + tab%dims(ish), &
                                               tab%offs(jsh) + 1:tab%offs(jsh) + tab%dims(jsh))))
         end do
         sa(ish) = sqrt(maxval(tab%env(tab%bas(LIBCINT_PTR_EXP, ish) + 1: &
                                       tab%bas(LIBCINT_PTR_EXP, ish) &
                                       + tab%bas(LIBCINT_NPRIM_OF, ish))))
      end do

      allocate (h1(nao, nao, 3, natm))
      h1 = 0.0_dp

      ! A private accumulator per thread is `nao^2 * 3 * natm` doubles, the
      ! same size as the answer. The first thing to revisit if the memory stops
      ! being affordable; the alternative is an atomic per deposit, and there
      ! are seven per buffer element.
      !$omp parallel default(none) &
      !$omp shared(kx, jx, env_use, mol, density, h1, owner, tab, atm_flat, bas_flat, &
      !$omp        mx, nao, natm, n_pairs, bq, dsh, sa, tol) &
      !$omp private(buf, hloc, dims, shls, ish, jsh, ksh, lsh, di, dj, dk, dl, &
      !$omp         io, jo, ko, lo, i, j, k, l, comp, idx, ii, jj, kk, ll, at, b, &
      !$omp         pair, have, wmax, est)
      allocate (buf(mx**4*3))
      allocate (hloc(nao, nao, 3, natm))
      hloc = 0.0_dp

      !$omp do schedule(dynamic)
      do pair = 1, n_pairs
         ish = (pair - 1)/tab%nbas + 1
         jsh = pair - tab%nbas*(ish - 1)
         di = tab%dims(ish)
         io = tab%offs(ish)
         dj = tab%dims(jsh)
         jo = tab%offs(jsh)
         do ksh = 1, tab%nbas
            dk = tab%dims(ksh)
            ko = tab%offs(ksh)
            do lsh = 1, tab%nbas
               dl = tab%dims(lsh)
               lo = tab%offs(lsh)
               shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]
               dims = [di, dj, dk, dl]

               ! The largest density element any of the seven deposits can
               ! pick up, against the bound on a once-differentiated quartet.
               wmax = max(dsh(ksh, lsh), dsh(ish, jsh), dsh(jsh, ksh), &
                          dsh(ish, ksh), dsh(lsh, ish), dsh(lsh, jsh))
               est = 4.0_dp*sa(ish)*bq(ish, jsh)*bq(ksh, lsh)*wmax
               if (est < tol) cycle

               if (mol%cartesian) then
                  have = cint2e_ip1_cart(buf, shls, atm_flat, mol%natm, &
                                         bas_flat, tab%nbas, env_use, c_null_ptr) /= 0
               else
                  have = cint2e_ip1_sph(buf, shls, atm_flat, mol%natm, &
                                        bas_flat, tab%nbas, env_use, c_null_ptr) /= 0
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
                              hloc(ii, jj, comp, at) = hloc(ii, jj, comp, at) - jx*density(kk, ll)*b
                              hloc(jj, ii, comp, at) = hloc(jj, ii, comp, at) - jx*density(kk, ll)*b
                              hloc(kk, ll, comp, at) = hloc(kk, ll, comp, at) &
                                                       - jx*2.0_dp*density(ii, jj)*b
                              hloc(ii, ll, comp, at) = hloc(ii, ll, comp, at) &
                                                       + kx*0.5_dp*density(jj, kk)*b
                              hloc(jj, ll, comp, at) = hloc(jj, ll, comp, at) &
                                                       + kx*0.5_dp*density(ii, kk)*b
                              hloc(kk, jj, comp, at) = hloc(kk, jj, comp, at) &
                                                       + kx*0.5_dp*density(ll, ii)*b
                              hloc(kk, ii, comp, at) = hloc(kk, ii, comp, at) &
                                                       + kx*0.5_dp*density(ll, jj)*b
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

      deallocate (atm_flat, bas_flat, owner, offsets, counts, bounds, bq, dsh, sa)
   end subroutine h1_contract

end module mqc_libcint_hess_ints
