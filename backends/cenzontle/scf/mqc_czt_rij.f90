!! The fitted Coulomb matrix on its own, without the fitted three-index tensor
module mqc_czt_rij
   !! `J(mu, nu) = sum_P (mu nu|P) c_P`, with `c = V^-1 gamma` and
   !! `gamma_P = sum_(mu nu) (mu nu|P) D(mu, nu)`.
   !!
   !! This is the Coulomb half of a density-fitted Fock build and nothing else.
   !! The fitted exchange needs `B = (mu nu|P) V^-1/2` whole, n_ao^2 by n_aux,
   !! and a build that also wants the exchange keeps that tensor. A build that
   !! takes its exchange from somewhere else does not need it: the Coulomb
   !! matrix is two contractions of the raw three-centre integrals with a
   !! linear solve between them.
   !!
   !! The integrals are kept only for shell pairs that survive the Schwarz
   !! screen, and only the lower triangle of those, which is a fraction of the
   !! square tensor and shrinks relative to it as the molecule grows. When
   !! even that does not fit in memory they are not kept at all, and each
   !! build makes two passes over them.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_trsm, pic_gemv
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_memory, only: memory_budget
   use mqc_program_limits, only: DF_PAIR_SCREEN
   use mqc_czt_integrals, only: czt_molecule_t, shell_dim, two_centre, build_df_shell_table, &
                                shell_pair_bounds, aux_shell_bounds, fit_metric_factor, &
                                max_block, angular_form_name
   use libcint_fortran, only: libcint_3c2e_sph, libcint_3c2e_cart, &
                              libcint_3c2e_sph_optimizer, libcint_3c2e_cart_optimizer, &
                              libcint_del_optimizer
   implicit none
   private

   public :: rij_t
   public :: build_rij
   public :: rij_coulomb

   real(dp), parameter :: RIJ_MEMORY_SHARE = 0.5_dp
      !! Share of the machine's available memory the stored integrals may take
   real(dp), parameter :: RIJ_MEMORY_BLIND = 4.0e9_dp
      !! Bytes they may take on a machine that reports nothing

   type :: rij_t
      !! What a fitted Coulomb build needs, held for an SCF
      logical :: in_core = .false.
         !! Whether the three-centre integrals are stored. False, they are
         !! recomputed in every build.
      logical :: cartesian = .false.
      integer :: n_ao = 0
      integer :: n_aux = 0
      integer :: n_pair = 0                    !! Orbital shell pairs kept, lower triangle
      integer :: n_row = 0                     !! Function pairs those hold
      integer :: n_shell_orb = 0
      integer :: n_shell_aux = 0
      integer :: dummy = 0                     !! libcint's "no fourth centre" shell, 1-based
      integer, allocatable :: bas(:, :)        !! Orbital then auxiliary shells, one table
      real(dp), allocatable :: env(:)
      integer, allocatable :: pair_i(:)        !! (n_pair) shells of each kept pair, i >= j
      integer, allocatable :: pair_j(:)
      integer, allocatable :: pair_row(:)      !! (n_pair) 0-based first row of each pair
      real(dp), allocatable :: pair_bound(:)   !! (n_pair) Schwarz bound of each pair
      real(dp), allocatable :: aux_bound(:)    !! (n_shell_aux) the auxiliary half
      integer, allocatable :: aux_offset(:)    !! (n_shell_aux) 0-based first function
      real(dp), allocatable :: three(:, :)
         !! (n_row, n_aux) the kept integrals; allocated only when `in_core`.
         !! A pair's block runs the first shell's functions fastest.
      real(dp), allocatable :: factor(:, :)    !! The metric, factored
      logical :: cholesky = .false.
         !! True: `factor` is U with V = U^T U. False: it is V^-1/2 in full.
   contains
      procedure :: active => rij_active
      procedure :: destroy => rij_destroy
   end type rij_t

contains

   pure function rij_active(this) result(active)
      !! Whether this has been built
      class(rij_t), intent(in) :: this
      logical :: active

      active = allocated(this%factor)
   end function rij_active

   pure subroutine rij_destroy(this)
      !! Release everything
      class(rij_t), intent(inout) :: this

      if (allocated(this%bas)) deallocate (this%bas)
      if (allocated(this%env)) deallocate (this%env)
      if (allocated(this%pair_i)) deallocate (this%pair_i)
      if (allocated(this%pair_j)) deallocate (this%pair_j)
      if (allocated(this%pair_row)) deallocate (this%pair_row)
      if (allocated(this%pair_bound)) deallocate (this%pair_bound)
      if (allocated(this%aux_bound)) deallocate (this%aux_bound)
      if (allocated(this%aux_offset)) deallocate (this%aux_offset)
      if (allocated(this%three)) deallocate (this%three)
      if (allocated(this%factor)) deallocate (this%factor)
      this%in_core = .false.
      this%n_pair = 0
      this%n_row = 0
   end subroutine rij_destroy

   subroutine build_rij(orb, aux, rij, error, in_core)
      !! Set up the fitted Coulomb build for a molecule
      !!
      !! Factors the metric, screens the orbital shell pairs, and stores the
      !! surviving three-centre integrals if they fit.
      type(czt_molecule_t), intent(in) :: orb   !! Orbital basis
      type(czt_molecule_t), intent(in) :: aux   !! Auxiliary basis, same atoms
      type(rij_t), intent(out) :: rij
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: in_core
         !! Store the integrals (true) or recompute them every build (false).
         !! Absent, they are stored when they fit in the memory this run may
         !! plan on.

      real(dp), allocatable :: metric(:, :), q_pair(:, :)
      real(dp) :: q_aux_max, bytes
      integer :: ish, jsh, ksh, n_pair_max, row

      if (orb%cartesian .neqv. aux%cartesian) then
         call error%set(ERROR_VALIDATION, "density fitting: the orbital basis is "// &
                        angular_form_name(orb%cartesian)//" and the auxiliary basis is "// &
                        angular_form_name(aux%cartesian)//". The fitting integrals are "// &
                        "built in one form for all three centres, so the two must agree.")
         return
      end if

      rij%cartesian = orb%cartesian
      rij%n_ao = orb%nao
      rij%n_aux = aux%nao
      rij%n_shell_orb = orb%nbas
      rij%n_shell_aux = aux%nbas
      call build_df_shell_table(orb, aux, rij%bas, rij%env, rij%dummy)

      call two_centre(aux, metric)
      rij%aux_bound = aux_shell_bounds(aux, metric)
      call fit_metric_factor(metric, rij%factor, rij%cholesky, error)
      if (error%has_error()) return
      deallocate (metric)

      allocate (rij%aux_offset(aux%nbas))
      do ksh = 1, aux%nbas
         rij%aux_offset(ksh) = aux%shell_offset(ksh)
      end do

      ! (mu nu|P) is symmetric in its orbital pair, so the lower triangle of
      ! shell pairs carries everything, and a pair whose Schwarz bound cannot
      ! reach even the largest auxiliary shell carries nothing.
      call shell_pair_bounds(orb, q_pair)
      q_aux_max = maxval(rij%aux_bound)
      n_pair_max = orb%nbas*(orb%nbas + 1)/2
      allocate (rij%pair_i(n_pair_max), rij%pair_j(n_pair_max), rij%pair_row(n_pair_max), &
                rij%pair_bound(n_pair_max))
      rij%n_pair = 0
      row = 0
      do ish = 1, orb%nbas
         do jsh = 1, ish
            if (q_pair(ish, jsh)*q_aux_max < DF_PAIR_SCREEN) cycle
            rij%n_pair = rij%n_pair + 1
            rij%pair_i(rij%n_pair) = ish
            rij%pair_j(rij%n_pair) = jsh
            rij%pair_row(rij%n_pair) = row
            rij%pair_bound(rij%n_pair) = q_pair(ish, jsh)
            row = row + shell_dim(orb%cartesian, ish - 1, rij%bas) &
                  *shell_dim(orb%cartesian, jsh - 1, rij%bas)
         end do
      end do
      rij%n_row = row

      if (present(in_core)) then
         rij%in_core = in_core
      else
         bytes = 8.0_dp*real(rij%n_row, dp)*real(rij%n_aux, dp)
         rij%in_core = bytes <= memory_budget(RIJ_MEMORY_BLIND, RIJ_MEMORY_SHARE)
      end if
      if (rij%in_core) call store_integrals(orb, rij)
   end subroutine build_rij

   subroutine store_integrals(orb, rij)
      !! Fill `rij%three` with the kept three-centre integrals
      type(czt_molecule_t), intent(in) :: orb
      type(rij_t), intent(inout) :: rij

      real(dp), allocatable :: buf(:)
      type(c_ptr) :: opt
      integer :: ipair, ksh, di, dj, dk, ko, k, ij, n_ij, row0, kcol

      allocate (rij%three(rij%n_row, rij%n_aux))
      ! Zeroed in parallel: the first touch decides which socket's memory each
      ! page lives on.
      !$omp parallel do default(none) shared(rij) private(kcol) schedule(static)
      do kcol = 1, rij%n_aux
         rij%three(:, kcol) = 0.0_dp
      end do
      !$omp end parallel do

      call make_optimizer(orb, rij, opt)
      ! A shell pair owns its rows, so the threads write disjoint memory.
      !$omp parallel default(none) shared(orb, rij, opt) &
      !$omp    private(ipair, ksh, di, dj, dk, ko, k, ij, n_ij, row0, buf)
      allocate (buf(max_block(orb)**2*maxval(aux_dims(rij))))
      !$omp do schedule(dynamic)
      do ipair = 1, rij%n_pair
         di = shell_dim(rij%cartesian, rij%pair_i(ipair) - 1, rij%bas)
         dj = shell_dim(rij%cartesian, rij%pair_j(ipair) - 1, rij%bas)
         n_ij = di*dj
         row0 = rij%pair_row(ipair)
         do ksh = 1, rij%n_shell_aux
            if (rij%pair_bound(ipair)*rij%aux_bound(ksh) < DF_PAIR_SCREEN) cycle
            if (.not. triplet(orb, rij, opt, ipair, ksh, buf)) cycle
            dk = shell_dim(rij%cartesian, rij%n_shell_orb + ksh - 1, rij%bas)
            ko = rij%aux_offset(ksh)
            do k = 1, dk
               do ij = 1, n_ij
                  rij%three(row0 + ij, ko + k) = buf(ij + (k - 1)*n_ij)
               end do
            end do
         end do
      end do
      !$omp end do
      deallocate (buf)
      !$omp end parallel
      call libcint_del_optimizer(opt)
   end subroutine store_integrals

   function aux_dims(rij) result(dims)
      !! Functions in each auxiliary shell
      type(rij_t), intent(in) :: rij
      integer :: dims(rij%n_shell_aux)

      integer :: ksh

      do ksh = 1, rij%n_shell_aux
         dims(ksh) = shell_dim(rij%cartesian, rij%n_shell_orb + ksh - 1, rij%bas)
      end do
   end function aux_dims

   subroutine make_optimizer(orb, rij, opt)
      !! libcint's per-shell-pair setup, built once and shared by every thread
      type(czt_molecule_t), intent(in) :: orb
      type(rij_t), intent(in) :: rij
      type(c_ptr), intent(out) :: opt

      opt = c_null_ptr
      if (rij%cartesian) then
         call libcint_3c2e_cart_optimizer(opt, orb%atm, orb%natm, rij%bas, &
                                          rij%n_shell_orb + rij%n_shell_aux + 1, rij%env)
      else
         call libcint_3c2e_sph_optimizer(opt, orb%atm, orb%natm, rij%bas, &
                                         rij%n_shell_orb + rij%n_shell_aux + 1, rij%env)
      end if
   end subroutine make_optimizer

   function triplet(orb, rij, opt, ipair, ksh, buf) result(nonzero)
      !! `(ij|K)` for one kept shell pair and one auxiliary shell
      !!
      !! `buf` runs the first orbital shell fastest, then the second, then the
      !! auxiliary one. False means libcint screened the whole block away.
      type(czt_molecule_t), intent(in) :: orb
      type(rij_t), intent(in) :: rij
      type(c_ptr), intent(in) :: opt
      integer, intent(in) :: ipair, ksh
      real(dp), intent(inout) :: buf(:)
      logical :: nonzero

      integer :: shls(4)
      integer :: ret

      shls = [rij%pair_i(ipair) - 1, rij%pair_j(ipair) - 1, rij%n_shell_orb + ksh - 1, &
              rij%dummy - 1]
      if (rij%cartesian) then
         ret = libcint_3c2e_cart(buf, shls, orb%atm, orb%natm, rij%bas, &
                                 rij%n_shell_orb + rij%n_shell_aux + 1, rij%env, opt)
      else
         ret = libcint_3c2e_sph(buf, shls, orb%atm, orb%natm, rij%bas, &
                                rij%n_shell_orb + rij%n_shell_aux + 1, rij%env, opt)
      end if
      nonzero = ret /= 0
   end function triplet

   subroutine rij_coulomb(orb, rij, density, j, error)
      !! The fitted Coulomb matrix of a density
      !!
      !! Only the symmetric part of `density` contributes, as it does to any
      !! Coulomb matrix.
      type(czt_molecule_t), intent(in) :: orb
      type(rij_t), intent(in) :: rij
      real(dp), intent(in) :: density(:, :)    !! (n_ao, n_ao)
      real(dp), allocatable, intent(out) :: j(:, :)   !! (n_ao, n_ao), symmetric
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: d_row(:), j_row(:), gamma(:), fit(:)
      integer :: ipair, a, b, di, dj, io, jo, row
      real(dp) :: weight

      if (.not. rij%active()) then
         call error%set(ERROR_VALIDATION, "the fitted Coulomb build was not set up")
         return
      end if
      if (size(density, 1) /= rij%n_ao .or. size(density, 2) /= rij%n_ao) then
         call error%set(ERROR_VALIDATION, "fitted Coulomb: the density does not match the basis")
         return
      end if

      ! The density over the kept function pairs. An off-diagonal shell pair
      ! stands for itself and its mirror.
      allocate (d_row(rij%n_row), j_row(rij%n_row), gamma(rij%n_aux), fit(rij%n_aux))
      !$omp parallel do default(none) shared(orb, rij, density, d_row) &
      !$omp    private(ipair, a, b, di, dj, io, jo, row, weight) schedule(static)
      do ipair = 1, rij%n_pair
         di = shell_dim(rij%cartesian, rij%pair_i(ipair) - 1, rij%bas)
         dj = shell_dim(rij%cartesian, rij%pair_j(ipair) - 1, rij%bas)
         io = orb%shell_offset(rij%pair_i(ipair))
         jo = orb%shell_offset(rij%pair_j(ipair))
         weight = 1.0_dp
         if (rij%pair_i(ipair) /= rij%pair_j(ipair)) weight = 2.0_dp
         row = rij%pair_row(ipair)
         do b = 1, dj
            do a = 1, di
               row = row + 1
               d_row(row) = 0.5_dp*weight*(density(io + a, jo + b) + density(jo + b, io + a))
            end do
         end do
      end do
      !$omp end parallel do

      if (rij%in_core) then
         call contract_stored(rij, d_row, gamma)
      else
         call contract_direct(orb, rij, d_row, gamma, to_aux=.true.)
      end if
      call solve_metric(rij, gamma, fit)
      if (rij%in_core) then
         call expand_stored(rij, fit, j_row)
      else
         call contract_direct(orb, rij, j_row, fit, to_aux=.false.)
      end if

      allocate (j(rij%n_ao, rij%n_ao))
      j = 0.0_dp
      do ipair = 1, rij%n_pair
         di = shell_dim(rij%cartesian, rij%pair_i(ipair) - 1, rij%bas)
         dj = shell_dim(rij%cartesian, rij%pair_j(ipair) - 1, rij%bas)
         io = orb%shell_offset(rij%pair_i(ipair))
         jo = orb%shell_offset(rij%pair_j(ipair))
         row = rij%pair_row(ipair)
         do b = 1, dj
            do a = 1, di
               row = row + 1
               j(io + a, jo + b) = j_row(row)
               j(jo + b, io + a) = j_row(row)
            end do
         end do
      end do
   end subroutine rij_coulomb

   subroutine contract_stored(rij, d_row, gamma)
      !! `gamma_P = sum_rows (row|P) d(row)`, from the stored integrals
      type(rij_t), intent(in) :: rij
      real(dp), intent(in) :: d_row(:)
      real(dp), intent(out) :: gamma(:)

      integer :: kcol

      !$omp parallel do default(none) shared(rij, d_row, gamma) private(kcol) schedule(static)
      do kcol = 1, rij%n_aux
         gamma(kcol) = dot_product(rij%three(:, kcol), d_row)
      end do
      !$omp end parallel do
   end subroutine contract_stored

   subroutine expand_stored(rij, fit, j_row)
      !! `J(row) = sum_P (row|P) c_P`, from the stored integrals
      type(rij_t), intent(in) :: rij
      real(dp), intent(in) :: fit(:)
      real(dp), intent(out) :: j_row(:)

      real(dp), allocatable :: local(:)
      integer :: kcol

      j_row = 0.0_dp
      ! Each thread walks whole columns, which is the order they are stored
      ! in, and adds its share at the end.
      !$omp parallel default(none) shared(rij, fit, j_row) private(kcol, local)
      allocate (local(rij%n_row))
      local = 0.0_dp
      !$omp do schedule(static)
      do kcol = 1, rij%n_aux
         local = local + fit(kcol)*rij%three(:, kcol)
      end do
      !$omp end do
      !$omp critical(mqc_rij_expand)
      j_row = j_row + local
      !$omp end critical(mqc_rij_expand)
      deallocate (local)
      !$omp end parallel
   end subroutine expand_stored

   subroutine contract_direct(orb, rij, row_vector, aux_vector, to_aux)
      !! One pass over the three-centre integrals, contracting as they are produced
      !!
      !! `to_aux` true: `aux_vector(P) = sum_rows (row|P) row_vector(row)`.
      !! False: `row_vector(row) = sum_P (row|P) aux_vector(P)`.
      type(czt_molecule_t), intent(in) :: orb
      type(rij_t), intent(in) :: rij
      real(dp), intent(inout) :: row_vector(:)
      real(dp), intent(inout) :: aux_vector(:)
      logical, intent(in) :: to_aux

      real(dp), allocatable :: buf(:), local(:)
      type(c_ptr) :: opt
      integer :: ipair, ksh, di, dj, dk, ko, k, ij, n_ij, row0
      real(dp) :: total

      call make_optimizer(orb, rij, opt)
      if (to_aux) then
         aux_vector = 0.0_dp
      else
         row_vector = 0.0_dp
      end if

      !$omp parallel default(none) shared(orb, rij, opt, row_vector, aux_vector, to_aux) &
      !$omp    private(ipair, ksh, di, dj, dk, ko, k, ij, n_ij, row0, buf, local, total)
      allocate (buf(max_block(orb)**2*maxval(aux_dims(rij))))
      if (to_aux) then
         allocate (local(rij%n_aux))
         local = 0.0_dp
      end if
      !$omp do schedule(dynamic)
      do ipair = 1, rij%n_pair
         di = shell_dim(rij%cartesian, rij%pair_i(ipair) - 1, rij%bas)
         dj = shell_dim(rij%cartesian, rij%pair_j(ipair) - 1, rij%bas)
         n_ij = di*dj
         row0 = rij%pair_row(ipair)
         do ksh = 1, rij%n_shell_aux
            if (rij%pair_bound(ipair)*rij%aux_bound(ksh) < DF_PAIR_SCREEN) cycle
            if (.not. triplet(orb, rij, opt, ipair, ksh, buf)) cycle
            dk = shell_dim(rij%cartesian, rij%n_shell_orb + ksh - 1, rij%bas)
            ko = rij%aux_offset(ksh)
            if (to_aux) then
               do k = 1, dk
                  total = 0.0_dp
                  do ij = 1, n_ij
                     total = total + buf(ij + (k - 1)*n_ij)*row_vector(row0 + ij)
                  end do
                  local(ko + k) = local(ko + k) + total
               end do
            else
               ! A shell pair owns its rows, so this needs no accumulator.
               do k = 1, dk
                  do ij = 1, n_ij
                     row_vector(row0 + ij) = row_vector(row0 + ij) &
                                             + buf(ij + (k - 1)*n_ij)*aux_vector(ko + k)
                  end do
               end do
            end if
         end do
      end do
      !$omp end do
      if (to_aux) then
         !$omp critical(mqc_rij_direct)
         aux_vector = aux_vector + local
         !$omp end critical(mqc_rij_direct)
         deallocate (local)
      end if
      deallocate (buf)
      !$omp end parallel
      call libcint_del_optimizer(opt)
   end subroutine contract_direct

   subroutine solve_metric(rij, gamma, fit)
      !! `c = V^-1 gamma`, however the metric was factored
      type(rij_t), intent(in) :: rij
      real(dp), intent(in) :: gamma(:)
      real(dp), intent(out) :: fit(:)

      real(dp), allocatable :: rhs(:, :), half(:)

      if (rij%cholesky) then
         ! V = U^T U: one solve against each triangle.
         allocate (rhs(rij%n_aux, 1))
         rhs(:, 1) = gamma
         call pic_trsm(rij%factor, rhs, side="L", uplo="U", transa="T")
         call pic_trsm(rij%factor, rhs, side="L", uplo="U", transa="N")
         fit = rhs(:, 1)
      else
         ! V^-1 = V^-1/2 V^-1/2
         allocate (half(rij%n_aux))
         half = 0.0_dp
         fit = 0.0_dp
         call pic_gemv(rij%factor, gamma, half)
         call pic_gemv(rij%factor, half, fit)
      end if
   end subroutine solve_metric

end module mqc_czt_rij
