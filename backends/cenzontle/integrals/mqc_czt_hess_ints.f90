!! Second derivatives of the one- and two-electron integrals with respect to nuclei
module mqc_czt_hess_ints
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
   use, intrinsic :: iso_fortran_env, only: real64, int64
   use, intrinsic :: iso_c_binding, only: c_null_ptr
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, shell_dim, atom_ao_blocks, &
                                eri_shell_table_t, eri_shell_table, eri_schwarz_collapse
   use mqc_czt_direct, only: schwarz_bounds, direct_stats_t, pair_degeneracy, pair_work_order, &
                             omp_threads
   use pic_logger, only: logger => global_logger
   use mqc_program_limits, only: MAX_LINE_LENGTH
   use libcint_fortran, only: LIBCINT_NPRIM_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_RINV_ORIG, &
                              LIBCINT_PTR_RANGE_OMEGA
   use mqc_czt_hess_abi, only: &
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
   public :: hess_rinv_contract
   public :: HESS_ERI_II, HESS_ERI_IJ, HESS_ERI_IK
   public :: hess_2e_block
   public :: hess_2e_contract
   public :: hess_2e_skeleton_contract
   public :: eri_ip1_block
   public :: h1_contract

   integer, parameter :: N_ERI_HESS_BLOCKS = 6
      !! The distinct second-derivative blocks of one shell quartet: ii, jj,
      !! kk on the diagonal and ij, ik, jk off it. Translational invariance
      !! supplies the fourth centre, so there is no ll and no l cross term.
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

   real(dp), parameter :: H1_ACCUMULATOR_BUDGET = 8.0e9_dp
      !! Bytes the per-thread copies of the skeleton Fock derivative may take
      !! together before `h1_contract` splits the atoms into groups. A copy is
      !! `nao^2 * 3 * natm` doubles per thread, which is 30 GB for a 74-atom
      !! molecule in 6-31G* on 64 threads.

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
      type(czt_molecule_t), intent(in) :: mol
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
      atm_flat = reshape(mol%atm, [size(mol%atm)])
      bas_flat = reshape(mol%bas, [size(mol%bas)])
      allocate (matrix(mol%nao, mol%nao, N_COMPONENTS))
      matrix = 0.0_dp

      ! Threaded over the bra shell: every block a thread writes sits in its
      ! own rows, so nothing is shared but the molecule.
      !$omp parallel default(none) shared(mol, which, matrix, atm_flat, bas_flat, mx) &
      !$omp    private(ish, jsh, di, dj, io, jo, i, j, comp, shls, dims, buf, have)
      allocate (buf(mx*mx*N_COMPONENTS))
      !$omp do schedule(dynamic)
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
      !$omp end do
      deallocate (buf)
      !$omp end parallel

      deallocate (atm_flat, bas_flat)
   end subroutine hess_1e_block

   function drive(mol, which, buf, dims, shls, atm, bas) result(have)
      !! Dispatch one shell pair to the right entry point
      type(czt_molecule_t), intent(in) :: mol
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
      type(czt_molecule_t), intent(in) :: mol
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

      allocate (matrix(mol%nao, mol%nao, N_COMPONENTS))
      matrix = 0.0_dp
      atm_flat = reshape(mol%atm, [size(mol%atm)])
      bas_flat = reshape(mol%bas, [size(mol%bas)])
      allocate (env(size(mol%env)), source=mol%env)
      env(LIBCINT_PTR_RINV_ORIG + 1:LIBCINT_PTR_RINV_ORIG + 3) = mol%coords(:, iatom)

      !$omp parallel default(none) shared(mol, which, matrix, atm_flat, bas_flat, env, mx) &
      !$omp    private(ish, jsh, di, dj, io, jo, i, j, comp, shls, dims, buf, have)
      allocate (buf(mx*mx*N_COMPONENTS))
      !$omp do schedule(dynamic)
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
      !$omp end do
      deallocate (buf)
      !$omp end parallel

      deallocate (env, atm_flat, bas_flat)
   end subroutine hess_rinv_block

   subroutine hess_rinv_contract(mol, density, cross, error)
      !! Both `1/|r-R_C|` second derivatives, contracted, for every basis atom
      !! and every nucleus in one pass
      !!
      !! `cross(:, :, A, C)` is `sum_ij D_ij [ipiprinv + iprinvip]_ij` with the
      !! operator origin on nucleus `C` and the bra function `i` centred on
      !! atom `A` -- what `partial_hessian` used to assemble from two
      !! `hess_rinv_block` calls per nucleus, each a serial pass forming an
      !! `(nao, nao, 9)` matrix. Here the shell pairs are walked once, threaded,
      !! with the nuclei inside, and each block is contracted on the spot.
      !! Nothing of size `nao^2` is stored.
      !!
      !! Returned **unscaled**, as `hess_rinv_block` is: the charge and the
      !! sign of the electron-nucleus attraction are the caller's. The nine
      !! components are libcint's, `a` on the bra centre and `b` on the
      !! nucleus for `iprinvip` and both on the bra for `ipiprinv`.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: cross(:, :, :, :)   !! (3, 3, natm, natm)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: buf_ii(:), buf_ij(:), env(:), cloc(:, :, :, :)
      integer, allocatable :: atm_flat(:), bas_flat(:), offsets(:), counts(:), owner(:)
      real(dp) :: acc(N_COMPONENTS)
      integer :: dims(0:3), shls(0:1)
      integer :: ish, jsh, di, dj, io, jo, i, j, comp, mx, c, a, b, natm, nao, ia
      logical :: have_ii, have_ij

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

      mx = 0
      do ish = 1, mol%nbas
         mx = max(mx, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do
      atm_flat = reshape(mol%atm, [size(mol%atm)])
      bas_flat = reshape(mol%bas, [size(mol%bas)])

      allocate (cross(3, 3, natm, natm))
      cross = 0.0_dp

      ! `env` is private because the operator origin lives in it and each
      ! thread moves it from nucleus to nucleus.
      !$omp parallel default(none) shared(mol, density, cross, owner, atm_flat, bas_flat, mx, natm) &
      !$omp    private(ish, jsh, di, dj, io, jo, i, j, comp, c, a, b, ia, shls, dims, &
      !$omp            buf_ii, buf_ij, env, cloc, acc, have_ii, have_ij)
      allocate (buf_ii(mx*mx*N_COMPONENTS), buf_ij(mx*mx*N_COMPONENTS))
      allocate (env(size(mol%env)), cloc(3, 3, natm, natm))
      env = mol%env
      cloc = 0.0_dp
      !$omp do schedule(dynamic)
      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         ! Every function of a shell sits on the shell's atom.
         ia = owner(io + 1)
         do jsh = 1, mol%nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]
            dims = [di, dj, 1, 1]

            do c = 1, natm
               env(LIBCINT_PTR_RINV_ORIG + 1:LIBCINT_PTR_RINV_ORIG + 3) = mol%coords(:, c)
               if (mol%cartesian) then
                  have_ii = cint1e_ipiprinv_cart(buf_ii, shls, atm_flat, mol%natm, &
                                                 bas_flat, mol%nbas, env) /= 0
                  have_ij = cint1e_iprinvip_cart(buf_ij, shls, atm_flat, mol%natm, &
                                                 bas_flat, mol%nbas, env) /= 0
               else
                  have_ii = cint1e_ipiprinv_sph(buf_ii, shls, atm_flat, mol%natm, &
                                                bas_flat, mol%nbas, env) /= 0
                  have_ij = cint1e_iprinvip_sph(buf_ij, shls, atm_flat, mol%natm, &
                                                bas_flat, mol%nbas, env) /= 0
               end if
               if (.not. (have_ii .or. have_ij)) cycle
               if (.not. have_ii) buf_ii(1:di*dj*N_COMPONENTS) = 0.0_dp
               if (.not. have_ij) buf_ij(1:di*dj*N_COMPONENTS) = 0.0_dp

               do comp = 1, N_COMPONENTS
                  acc(comp) = 0.0_dp
                  do j = 1, dj
                     do i = 1, di
                        acc(comp) = acc(comp) + density(io + i, jo + j) &
                                    *(buf_ii(i + di*(j - 1 + dj*(comp - 1))) &
                                      + buf_ij(i + di*(j - 1 + dj*(comp - 1))))
                     end do
                  end do
               end do
               do comp = 1, N_COMPONENTS
                  a = (comp - 1)/3 + 1
                  b = comp - 3*(a - 1)
                  cloc(a, b, ia, c) = cloc(a, b, ia, c) + acc(comp)
               end do
            end do
         end do
      end do
      !$omp end do
      !$omp critical
      cross = cross + cloc
      !$omp end critical
      deallocate (buf_ii, buf_ij, env, cloc)
      !$omp end parallel

      deallocate (atm_flat, bas_flat, offsets, counts, owner)
   end subroutine hess_rinv_contract

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
      type(czt_molecule_t), intent(in) :: mol
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
      type(czt_molecule_t), intent(in) :: mol
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
      type(czt_molecule_t), intent(in) :: mol
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
                               j_scale, omega, stats, bounds)
      !! The two-electron second derivatives, contracted as they are computed
      !!
      !! Same numbers `hess_2e_block` produces and never the same array. Each
      !! shell quartet is evaluated into a buffer, deposited into the Hessian,
      !! and forgotten, so the memory is one quartet rather than `nao^4` times
      !! nine times three.
      !!
      !! **Unique quartets, sixteen blocks each.** The second derivative of
      !! `(ij|kl)` has a 3x3 block for every ordered pair of its four centres.
      !! Summed over all ordered quartets against a weight with the full
      !! eightfold symmetry, that sum is `deg * gam` times the sixteen-block
      !! sum on each *unique* quartet, because the sixteen-block sum is itself
      !! invariant under relabelling the quartet. So the loop walks the Fock
      !! build's unique quartets and deposits all sixteen blocks from each.
      !!
      !! **Six kernel calls give the sixteen blocks.** `ipip1` on `(ij|kl)`,
      !! `(ji|kl)` and `(kl|ij)` are the diagonal blocks of the first three
      !! centres; `ipvip1` on `(ij|kl)` is the `ij` block and `ip1ip2` on
      !! `(ij|kl)` and `(ji|kl)` the `ik` and `jk` blocks; their transposes
      !! are the `ji`, `ki`, `kj` blocks. Every block involving the fourth
      !! centre follows from translational invariance -- the derivatives with
      !! respect to the four centres sum to zero -- so `l` costs no kernel.
      !! That is six calls per unique quartet where the ordered-bra loop it
      !! replaces made four per half-ordered one: 0.75 against 2 per `nbas^4`.
      type(czt_molecule_t), intent(in) :: mol
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
      type(direct_stats_t), intent(out), optional :: stats
         !! Unique shell quartets walked, skipped by the estimate, and handed
         !! to the derivative kernels. A computed quartet is six kernel calls.
      real(dp), intent(in), optional :: bounds(:, :)
         !! From `schwarz_bounds`, when the caller already has them.
      real(dp), allocatable :: buf_ii(:), buf_jj(:), buf_kk(:)
      real(dp), allocatable :: buf_ij(:), buf_ik(:), buf_jk(:)
      real(dp), allocatable :: hloc(:, :, :, :)
      integer, allocatable :: atm_flat(:), bas_flat(:), owner(:)
      integer, allocatable :: offsets(:), counts(:), sh_dim(:), sh_off(:)
      integer, allocatable :: pair_i(:), pair_j(:), order(:)
      real(dp), allocatable :: qb(:, :), dsh(:, :), sa(:)
      real(dp) :: qq, wbound, est, tol, kx, jx, deg, gam, w
      real(dp) :: t(3, 3, 4, 4)
      real(dp), allocatable, target :: env_use(:)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: io, jo, ko, lo, i, j, k, l, comp, mx, n4
      integer :: o_ijkl, o_jikl, o_klij
      integer :: nao, natm, a, b, x, y, c
      integer :: at(4)
      integer :: npair, ipair, itask, ij, kl
      integer(int64) :: n_total, n_computed, n_screened
      logical :: have(N_ERI_HESS_BLOCKS)

      if (error%has_error()) return

      kx = 1.0_dp
      if (present(k_scale)) kx = k_scale
      jx = 1.0_dp
      if (present(j_scale)) jx = j_scale
      n_total = 0_int64
      n_computed = 0_int64
      n_screened = 0_int64
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

      ! The unique quartets as the Fock build enumerates them: shell pairs
      ! `i >= j`, and `kl <= ij`. Over the split shells, not the fused-sp
      ! view: these kernels read a fused entry as a p shell over the s
      ! coefficients, silently.
      npair = mol%nbas*(mol%nbas + 1)/2
      allocate (pair_i(npair), pair_j(npair))
      ipair = 0
      do ish = 1, mol%nbas
         do jsh = 1, ish
            ipair = ipair + 1
            pair_i(ipair) = ish
            pair_j(ipair) = jsh
         end do
      end do
      call pair_work_order(pair_i, pair_j, sh_dim, order)

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
      if (present(bounds)) then
         qb = bounds
      else
         call schwarz_bounds(mol, qb, error)
         if (error%has_error()) return
      end if

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
      !$omp        mx, natm, npair, pair_i, pair_j, order, qb, dsh, sa, tol) &
      !$omp private(buf_ii, buf_jj, buf_kk, buf_ij, buf_ik, buf_jk, hloc, t, &
      !$omp         ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo, i, j, k, l, &
      !$omp         comp, n4, o_ijkl, o_jikl, o_klij, at, a, b, x, y, deg, gam, w, &
      !$omp         itask, ij, kl, have, qq, wbound, est) &
      !$omp reduction(+:n_total, n_computed, n_screened)
      allocate (buf_ii(mx**4*N_COMPONENTS), buf_jj(mx**4*N_COMPONENTS), &
                buf_kk(mx**4*N_COMPONENTS), buf_ij(mx**4*N_COMPONENTS), &
                buf_ik(mx**4*N_COMPONENTS), buf_jk(mx**4*N_COMPONENTS))
      allocate (hloc(3, 3, natm, natm))
      hloc = 0.0_dp

      !$omp do schedule(dynamic)
      do itask = 1, npair
         ij = order(itask)
         ish = pair_i(ij)
         jsh = pair_j(ij)
         di = sh_dim(ish)
         io = sh_off(ish)
         dj = sh_dim(jsh)
         jo = sh_off(jsh)
         do kl = 1, ij
            ksh = pair_i(kl)
            lsh = pair_j(kl)
            dk = sh_dim(ksh)
            ko = sh_off(ksh)
            dl = sh_dim(lsh)
            lo = sh_off(lsh)
            deg = pair_degeneracy(ish, jsh, ksh, lsh)

            ! The largest weight any element of this quartet can carry,
            ! written the way the contraction weight is: a Coulomb term and
            ! the two exchange terms the eightfold symmetrisation produces,
            ! times the degeneracy the fold puts on it.
            wbound = deg*(0.5_dp*dsh(ish, jsh)*dsh(ksh, lsh) &
                          + abs(kx)*0.125_dp*(dsh(ish, lsh)*dsh(ksh, jsh) &
                                              + dsh(ish, ksh)*dsh(jsh, lsh)))
            qq = qb(ish, jsh)*qb(ksh, lsh)
            ! Sixteen blocks, each two derivatives on some pair of centres.
            est = 2.0_dp*wbound*qq*(sa(ish) + sa(jsh) + sa(ksh) + sa(lsh))**2
            n_total = n_total + 1_int64
            if (est < tol) then
               n_screened = n_screened + 1_int64
               cycle
            end if
            n_computed = n_computed + 1_int64

            have(1) = drive_2e(mol, HESS_ERI_II, buf_ii, [di, dj, dk, dl], &
                               [ish - 1, jsh - 1, ksh - 1, lsh - 1], atm_flat, bas_flat, env_use)
            have(2) = drive_2e(mol, HESS_ERI_II, buf_jj, [dj, di, dk, dl], &
                               [jsh - 1, ish - 1, ksh - 1, lsh - 1], atm_flat, bas_flat, env_use)
            have(3) = drive_2e(mol, HESS_ERI_II, buf_kk, [dk, dl, di, dj], &
                               [ksh - 1, lsh - 1, ish - 1, jsh - 1], atm_flat, bas_flat, env_use)
            have(4) = drive_2e(mol, HESS_ERI_IJ, buf_ij, [di, dj, dk, dl], &
                               [ish - 1, jsh - 1, ksh - 1, lsh - 1], atm_flat, bas_flat, env_use)
            have(5) = drive_2e(mol, HESS_ERI_IK, buf_ik, [di, dj, dk, dl], &
                               [ish - 1, jsh - 1, ksh - 1, lsh - 1], atm_flat, bas_flat, env_use)
            have(6) = drive_2e(mol, HESS_ERI_IK, buf_jk, [dj, di, dk, dl], &
                               [jsh - 1, ish - 1, ksh - 1, lsh - 1], atm_flat, bas_flat, env_use)
            if (.not. any(have)) cycle
            ! A kernel that found nothing leaves its buffer as it was.
            n4 = di*dj*dk*dl
            if (.not. have(1)) buf_ii(1:n4*N_COMPONENTS) = 0.0_dp
            if (.not. have(2)) buf_jj(1:n4*N_COMPONENTS) = 0.0_dp
            if (.not. have(3)) buf_kk(1:n4*N_COMPONENTS) = 0.0_dp
            if (.not. have(4)) buf_ij(1:n4*N_COMPONENTS) = 0.0_dp
            if (.not. have(5)) buf_ik(1:n4*N_COMPONENTS) = 0.0_dp
            if (.not. have(6)) buf_jk(1:n4*N_COMPONENTS) = 0.0_dp

            do l = 1, dl
               at(4) = owner(lo + l)
               do k = 1, dk
                  at(3) = owner(ko + k)
                  do j = 1, dj
                     at(2) = owner(jo + j)
                     do i = 1, di
                        at(1) = owner(io + i)
                        gam = jx*0.5_dp*density(io + i, jo + j)*density(ko + k, lo + l) &
                              - kx*0.125_dp &
                              *(density(io + i, lo + l)*density(ko + k, jo + j) &
                                + density(io + i, ko + k)*density(jo + j, lo + l))
                        if (gam == 0.0_dp) cycle
                        w = deg*gam

                        ! Where this function quartet sits in each buffer,
                        ! before the component stride: the three shell orders
                        ! the kernels were called in.
                        o_ijkl = i + di*(j - 1 + dj*(k - 1 + dk*(l - 1)))
                        o_jikl = j + dj*(i - 1 + di*(k - 1 + dk*(l - 1)))
                        o_klij = k + dk*(l - 1 + dl*(i - 1 + di*(j - 1)))
                        do comp = 1, N_COMPONENTS
                           a = (comp - 1)/3 + 1
                           b = comp - 3*(a - 1)
                           t(a, b, 1, 1) = buf_ii(o_ijkl + n4*(comp - 1))
                           t(a, b, 2, 2) = buf_jj(o_jikl + n4*(comp - 1))
                           t(a, b, 3, 3) = buf_kk(o_klij + n4*(comp - 1))
                           t(a, b, 1, 2) = buf_ij(o_ijkl + n4*(comp - 1))
                           t(a, b, 1, 3) = buf_ik(o_ijkl + n4*(comp - 1))
                           t(a, b, 2, 3) = buf_jk(o_jikl + n4*(comp - 1))
                        end do
                        t(:, :, 2, 1) = transpose(t(:, :, 1, 2))
                        t(:, :, 3, 1) = transpose(t(:, :, 1, 3))
                        t(:, :, 3, 2) = transpose(t(:, :, 2, 3))
                        ! The fourth centre: minus the other three, in each
                        ! derivative index in turn.
                        do x = 1, 3
                           t(:, :, x, 4) = -(t(:, :, x, 1) + t(:, :, x, 2) + t(:, :, x, 3))
                           t(:, :, 4, x) = transpose(t(:, :, x, 4))
                        end do
                        t(:, :, 4, 4) = -(t(:, :, 4, 1) + t(:, :, 4, 2) + t(:, :, 4, 3))

                        do y = 1, 4
                           do x = 1, 4
                              hloc(:, :, at(x), at(y)) = hloc(:, :, at(x), at(y)) + w*t(:, :, x, y)
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
      deallocate (buf_ii, buf_jj, buf_kk, buf_ij, buf_ik, buf_jk, hloc)
      !$omp end parallel

      if (present(stats)) then
         stats%quartets_total = n_total
         stats%quartets_computed = n_computed
         stats%quartets_screened = n_screened
      end if

      deallocate (atm_flat, bas_flat, owner, offsets, counts, sh_dim, sh_off)
      deallocate (qb, dsh, sa, pair_i, pair_j, order)
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
      type(czt_molecule_t), intent(in) :: mol
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
                          j_scale, omega, stats, memory_budget, bounds)
      !! The skeleton derivative Fock for **every** atom, in one pass
      !!
      !! `make_h1_atom` builds one atom's `dF/dR` from a stored `int2e_ip1`
      !! array, at `nao^4` of memory and another `nao^4` of work per atom. This
      !! walks the shell quartets once and deposits into whichever atom each
      !! quartet's differentiated function belongs to, so the `natm` leaves the
      !! work and the array is not needed.
      !!
      !! **Unique quartets, the Fock build's six updates, one centre at a
      !! time.** The derivative of `F` with respect to an atom is the Fock
      !! functional applied to the derivative integrals, and for a fixed
      !! function `X` the number `d(ij|kl)/dX` is the same however the quartet
      !! is ordered. So each unique quartet is folded exactly as
      !! `build_fock_direct_many` folds an undifferentiated one -- `deg` for
      !! the orbit, six updates, symmetrised at the end -- once per centre,
      !! into that centre's atom. `int2e_ip1` on `(ij|kl)`, `(ji|kl)` and
      !! `(kl|ij)` give the first three centres; the fourth is minus their sum.
      !! Three calls per unique quartet where the ordered loop it replaces
      !! made one per ordered one: 3/8 against 1 per `nbas^4`.
      !!
      !! **The accumulator is bounded, not the thread count.** A private copy
      !! of `h1` per thread is `nao^2 * 3 * natm` doubles each, which is tens
      !! of gigabytes at a few hundred atoms on a full node. The atoms are
      !! taken in groups sized so the copies fit `memory_budget`, each group
      !! walking the quartets that touch it; a quartet with centres in several
      !! groups is computed once per group, which is the price of the bound
      !! and is paid only when the bound bites. The copies are one shared
      !! array so the reduction runs in parallel over atoms and components
      !! rather than through a critical section.
      !!
      !! Checked against `make_h1_atom` element by element in the test suite,
      !! that routine staying as the readable statement of what this computes.
      ! TODO(mqc): with the budget binding, a quartet is recomputed once per
      ! group it touches, up to four times; the cheaper alternative is to
      ! deposit those into a smaller structure and was not tried.
!$    use omp_lib, only: omp_get_thread_num
      type(czt_molecule_t), intent(in) :: mol
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
      type(direct_stats_t), intent(out), optional :: stats
         !! Unique shell quartets walked, skipped by the estimate, and handed
         !! to `int2e_ip1`, over the fused-sp view, summed over atom groups. A
         !! computed quartet is three kernel calls.
      real(dp), intent(in), optional :: memory_budget
         !! Bytes the per-thread accumulators may take together;
         !! `H1_ACCUMULATOR_BUDGET` when absent.
      real(dp), intent(in), optional :: bounds(:, :)
         !! From `schwarz_bounds`, over the split shells, when the caller
         !! already has them.

      real(dp), allocatable :: buf_i(:), buf_j(:), buf_k(:)
      real(dp), allocatable :: acc(:, :, :, :, :)
      real(dp), allocatable :: d_half(:, :)
      real(dp), allocatable :: qb(:, :), bq(:, :), dsh(:, :), sa(:)
      real(dp) :: wmax, est, tol, kx, jx, deg, budget, val, jsc, ksc
      real(dp) :: v(3, 4)
      real(dp), allocatable, target :: env_use(:)
      integer, allocatable :: atm_flat(:), bas_flat(:), owner(:), shell_atom(:)
      integer, allocatable :: offsets(:), counts(:)
      integer, allocatable :: pair_i(:), pair_j(:), order(:)
      type(eri_shell_table_t) :: tab
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: io, jo, ko, lo, i, j, k, l, comp, mx, n4
      integer :: o_ijkl, o_jikl, o_klij
      integer :: nao, natm, c, x, slot
      integer :: at(4)
      integer :: npair, ipair, itask, ij, kl, nthreads, tid
      integer :: n_groups, group, first, last, wide
      integer(int64) :: n_total, n_computed, n_screened, bytes_full
      logical :: have(3)
      character(len=MAX_LINE_LENGTH) :: line

      if (error%has_error()) return

      kx = 1.0_dp
      if (present(k_scale)) kx = k_scale
      jx = 1.0_dp
      if (present(j_scale)) jx = j_scale
      nao = mol%nao
      natm = mol%natm
      n_total = 0_int64
      n_computed = 0_int64
      n_screened = 0_int64

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
      allocate (shell_atom(tab%nbas))
      do ish = 1, tab%nbas
         shell_atom(ish) = owner(tab%offs(ish) + 1)
      end do

      ! A local copy so an omega pass can set the range-separation slot without
      ! disturbing the table every other caller shares. **From `tab%env` and
      ! not `mol%env`**: this loop walks the fused-sp view, whose `bas` points
      ! into its own `env`, so a copy of the molecule's is the wrong array to
      ! hand the integral.
      allocate (env_use(0:size(tab%env) - 1), source=tab%env)
      if (present(omega)) env_use(LIBCINT_PTR_RANGE_OMEGA) = omega

      atm_flat = reshape(mol%atm, [size(mol%atm)])
      bas_flat = reshape(tab%bas, [size(tab%bas)])

      ! The unique quartets, as the Fock build enumerates and orders them.
      npair = tab%nbas*(tab%nbas + 1)/2
      allocate (pair_i(npair), pair_j(npair))
      ipair = 0
      do ish = 1, tab%nbas
         do jsh = 1, ish
            ipair = ipair + 1
            pair_i(ipair) = ish
            pair_j(ipair) = jsh
         end do
      end do
      call pair_work_order(pair_i, pair_j, tab%dims, order)

      ! Screened the same way `hess_2e_contract` is, with one derivative rather
      ! than two, so the estimate carries one factor of `2 sqrt(alpha_max)`
      ! per centre. The bounds come from the split shells and are re-blocked
      ! onto the view, because that is the table this loop walks.
      tol = HESS_SCREEN_TOL
      if (present(screen_tol)) tol = screen_tol
      if (present(bounds)) then
         qb = bounds
      else
         call schwarz_bounds(mol, qb, error)
         if (error%has_error()) return
      end if
      call eri_schwarz_collapse(mol, qb, bq)

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

      ! The six-update form is written for the density without its two.
      allocate (d_half(nao, nao))
      d_half = 0.5_dp*density

      allocate (h1(nao, nao, 3, natm))
      h1 = 0.0_dp

      ! How many atoms the copies can hold at once.
      budget = H1_ACCUMULATOR_BUDGET
      if (present(memory_budget)) budget = memory_budget
      nthreads = omp_threads()
      bytes_full = int(nao, int64)*int(nao, int64)*3_int64*int(natm, int64)*8_int64
      n_groups = int(ceiling(real(nthreads, dp)*real(bytes_full, dp)/budget))
      n_groups = max(1, min(n_groups, natm))
      wide = (natm + n_groups - 1)/n_groups
      n_groups = (natm + wide - 1)/wide
      if (n_groups > 1) then
         write (line, "(a,i0,a,i0,a)") "    skeleton Fock derivative: accumulator held to ", &
            n_groups, " atom groups of ", wide, " to fit the memory budget"
         call logger%performance(trim(line))
      end if
      allocate (acc(nao, nao, 3, wide, nthreads))

      do group = 1, n_groups
         first = (group - 1)*wide + 1
         last = min(group*wide, natm)

         !$omp parallel default(none) &
         !$omp shared(kx, jx, env_use, mol, density, d_half, h1, owner, shell_atom, tab, &
         !$omp        atm_flat, bas_flat, mx, nao, natm, npair, pair_i, pair_j, order, &
         !$omp        bq, dsh, sa, tol, acc, first, last, wide, nthreads) &
         !$omp private(buf_i, buf_j, buf_k, tid, ish, jsh, ksh, lsh, di, dj, dk, dl, &
         !$omp         io, jo, ko, lo, i, j, k, l, comp, n4, o_ijkl, o_jikl, o_klij, &
         !$omp         itask, ij, kl, deg, wmax, est, have, at, x, slot, v, val, jsc, ksc) &
         !$omp reduction(+:n_total, n_computed, n_screened)
         tid = 1
!$       tid = omp_get_thread_num() + 1
         allocate (buf_i(mx**4*3), buf_j(mx**4*3), buf_k(mx**4*3))
         ! Each thread zeroes its own copy, so the pages are first touched by
         ! the thread that will write them.
         acc(:, :, :, :, tid) = 0.0_dp
         !$omp barrier

         !$omp do schedule(dynamic)
         do itask = 1, npair
            ij = order(itask)
            ish = pair_i(ij)
            jsh = pair_j(ij)
            di = tab%dims(ish)
            io = tab%offs(ish)
            dj = tab%dims(jsh)
            jo = tab%offs(jsh)
            do kl = 1, ij
               ksh = pair_i(kl)
               lsh = pair_j(kl)
               ! Nothing to deposit in this group from a quartet none of
               ! whose centres belong to it.
               if (.not. (in_group(shell_atom(ish)) .or. in_group(shell_atom(jsh)) &
                          .or. in_group(shell_atom(ksh)) .or. in_group(shell_atom(lsh)))) cycle
               dk = tab%dims(ksh)
               ko = tab%offs(ksh)
               dl = tab%dims(lsh)
               lo = tab%offs(lsh)
               deg = pair_degeneracy(ish, jsh, ksh, lsh)

               ! The largest density element any of the six updates can pick
               ! up, against the bound on a once-differentiated quartet,
               ! summed over the four centres the deposits come from.
               wmax = max(dsh(ksh, lsh), dsh(ish, jsh), dsh(jsh, ksh), &
                          dsh(ish, ksh), dsh(lsh, ish), dsh(lsh, jsh))
               est = 4.0_dp*deg*(sa(ish) + sa(jsh) + sa(ksh) + sa(lsh)) &
                     *bq(ish, jsh)*bq(ksh, lsh)*wmax
               n_total = n_total + 1_int64
               if (est < tol) then
                  n_screened = n_screened + 1_int64
                  cycle
               end if
               n_computed = n_computed + 1_int64

               have(1) = ip1_block(mol, buf_i, [ish - 1, jsh - 1, ksh - 1, lsh - 1], &
                                   atm_flat, bas_flat, tab%nbas, env_use)
               have(2) = ip1_block(mol, buf_j, [jsh - 1, ish - 1, ksh - 1, lsh - 1], &
                                   atm_flat, bas_flat, tab%nbas, env_use)
               have(3) = ip1_block(mol, buf_k, [ksh - 1, lsh - 1, ish - 1, jsh - 1], &
                                   atm_flat, bas_flat, tab%nbas, env_use)
               if (.not. any(have)) cycle
               n4 = di*dj*dk*dl
               if (.not. have(1)) buf_i(1:n4*3) = 0.0_dp
               if (.not. have(2)) buf_j(1:n4*3) = 0.0_dp
               if (.not. have(3)) buf_k(1:n4*3) = 0.0_dp

               do l = 1, dl
                  at(4) = owner(lo + l)
                  do k = 1, dk
                     at(3) = owner(ko + k)
                     do j = 1, dj
                        at(2) = owner(jo + j)
                        do i = 1, di
                           at(1) = owner(io + i)
                           o_ijkl = i + di*(j - 1 + dj*(k - 1 + dk*(l - 1)))
                           o_jikl = j + dj*(i - 1 + di*(k - 1 + dk*(l - 1)))
                           o_klij = k + dk*(l - 1 + dl*(i - 1 + di*(j - 1)))
                           do comp = 1, 3
                              v(comp, 1) = buf_i(o_ijkl + n4*(comp - 1))
                              v(comp, 2) = buf_j(o_jikl + n4*(comp - 1))
                              v(comp, 3) = buf_k(o_klij + n4*(comp - 1))
                           end do
                           v(:, 4) = -(v(:, 1) + v(:, 2) + v(:, 3))

                           ! The Fock build's six updates, once per centre,
                           ! into that centre's atom. The library's nabla is
                           ! minus the derivative with respect to the atom.
                           do x = 1, 4
                              if (.not. in_group(at(x))) cycle
                              slot = at(x) - first + 1
                              do comp = 1, 3
                                 val = -deg*v(comp, x)
                                 if (val == 0.0_dp) cycle
                                 jsc = jx*val
                                 ksc = kx*0.25_dp*val
                                 acc(io + i, jo + j, comp, slot, tid) = &
                                    acc(io + i, jo + j, comp, slot, tid) + jsc*d_half(ko + k, lo + l)
                                 acc(ko + k, lo + l, comp, slot, tid) = &
                                    acc(ko + k, lo + l, comp, slot, tid) + jsc*d_half(io + i, jo + j)
                                 acc(io + i, ko + k, comp, slot, tid) = &
                                    acc(io + i, ko + k, comp, slot, tid) - ksc*d_half(jo + j, lo + l)
                                 acc(jo + j, lo + l, comp, slot, tid) = &
                                    acc(jo + j, lo + l, comp, slot, tid) - ksc*d_half(io + i, ko + k)
                                 acc(io + i, lo + l, comp, slot, tid) = &
                                    acc(io + i, lo + l, comp, slot, tid) - ksc*d_half(jo + j, ko + k)
                                 acc(jo + j, ko + k, comp, slot, tid) = &
                                    acc(jo + j, ko + k, comp, slot, tid) - ksc*d_half(io + i, lo + l)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         !$omp end do

         ! The reduction, in parallel over the group's atoms and components:
         ! every copy summed, then symmetrised, which is what turns the six
         ! updates into the full sum over the eight permutations.
         !$omp do collapse(2)
         do slot = 1, last - first + 1
            do comp = 1, 3
               do j = 1, nao
                  do i = 1, nao
                     h1(i, j, comp, first + slot - 1) = sum(acc(i, j, comp, slot, 1:nthreads))
                  end do
               end do
               do j = 1, nao
                  do i = 1, j - 1
                     val = 0.5_dp*(h1(i, j, comp, first + slot - 1) + h1(j, i, comp, first + slot - 1))
                     h1(i, j, comp, first + slot - 1) = val
                     h1(j, i, comp, first + slot - 1) = val
                  end do
               end do
            end do
         end do
         !$omp end do
         deallocate (buf_i, buf_j, buf_k)
         !$omp end parallel
      end do

      if (present(stats)) then
         stats%quartets_total = n_total
         stats%quartets_computed = n_computed
         stats%quartets_screened = n_screened
      end if

      deallocate (acc, d_half, atm_flat, bas_flat, owner, shell_atom, offsets, counts)
      deallocate (qb, bq, dsh, sa, pair_i, pair_j, order)

   contains

      pure function in_group(atom) result(is_in)
         !! Whether `atom` is in the group of atoms being accumulated
         integer, intent(in) :: atom
         logical :: is_in

         is_in = atom >= first .and. atom <= last
      end function in_group
   end subroutine h1_contract

   function ip1_block(mol, buf, shls, atm, bas, nbas, env) result(have)
      !! `int2e_ip1` on one shell quartet, spherical or Cartesian as the
      !! molecule is
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(inout) :: buf(0:)
      integer, intent(in) :: shls(0:)
      integer, intent(in), target :: atm(0:), bas(0:)
      integer, intent(in) :: nbas
      real(dp), intent(in), target :: env(0:)
      logical :: have

      if (mol%cartesian) then
         have = cint2e_ip1_cart(buf, shls, atm, mol%natm, bas, nbas, env, c_null_ptr) /= 0
      else
         have = cint2e_ip1_sph(buf, shls, atm, mol%natm, bas, nbas, env, c_null_ptr) /= 0
      end if
   end function ip1_block

end module mqc_czt_hess_ints
