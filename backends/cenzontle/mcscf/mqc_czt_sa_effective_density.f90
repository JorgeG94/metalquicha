!! The SA-CASSCF gradient and coupling columns as an effective density
module mqc_czt_sa_effective_density
   !! `sa_density_t` is what `sa_gradients_on_state` hands
   !! `contract_effective_density`: every root's relaxed gradient column and
   !! every extra Lagrangian column, over one shared density pool.
   !!
   !! The pool is the stack the SA path has always swept: `d_core`, then the
   !! SA active density, then per root its own active density followed by the
   !! response's active and core parts, then per extra column the response's
   !! two. A root column's separable part is
   !!
   !!     (core, core, 1/2, 1/4) (core, own, 1, 1/2) (own, own, 1/2, 1/4)
   !!     (core, active_bar, 1, 1/2) (core, core_bar, 1, 1/2)
   !!     (active_sa, core_bar, 1, 1/2)
   !!
   !! and an extra column has the last three. Every column's active two-body
   !! density, with the SA one on each response leg, is the Gamma.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_effective_density, only: effective_density_t, separable_pair_t, column_t
   use mqc_czt_mcscf_gradient, only: ao_gamma_block => gamma_block
   implicit none
   private

   public :: sa_density_t
   public :: sa_pool_slots

   integer, parameter :: MAX_PAIRS = 6
      !! A root column's separable pairs; an extra column has three

   type, extends(effective_density_t) :: sa_density_t
      !! Root and extra columns of one SA-CASSCF state
      integer :: n_roots = 0
         !! Columns `1..n_roots` are root gradients, the rest extra columns
      type(column_t), allocatable :: columns(:)
         !! (n_columns)
      real(dp), allocatable :: p(:, :, :)
         !! (n_ao, n_ao, n_columns): one-particle density of each column
      real(dp), allocatable :: w(:, :, :)
         !! (n_ao, n_ao, n_columns): energy-weighted density of each column
      real(dp), allocatable :: pool(:, :, :)
         !! (n_ao, n_ao, n_pool), in the slot order `sa_pool_slots` gives
      integer, allocatable :: ao_map(:)
         !! (n_ao): compressed position of each AO carrying Gamma, 0 for none
      real(dp), allocatable :: c_sig(:, :)
         !! (n_sig, n_active): active orbitals over the AOs carrying Gamma
      real(dp), allocatable :: cbar_sig(:, :, :)
         !! (n_sig, n_active, n_columns): each column's response leg, likewise
      real(dp), allocatable :: dm2_column(:, :, :, :, :)
         !! (n_active, n_active, n_active, n_active, n_columns): each column's
         !! own active two-body density on the all-`C` legs
      real(dp), allocatable :: dm2_sa(:, :, :, :)
         !! (n_active, n_active, n_active, n_active): the SA two-body density on
         !! each response leg
   contains
      procedure :: n_columns => sa_n_columns
      procedure :: column => sa_column
      procedure :: energy_weighted => sa_energy_weighted
      procedure :: one_particle => sa_one_particle
      procedure :: density_pool => sa_density_pool
      procedure :: separable_pairs => sa_separable_pairs
      procedure :: gamma_ao_map => sa_gamma_ao_map
      procedure :: gamma_block => sa_gamma_block
   end type sa_density_t

contains

   pure subroutine sa_pool_slots(ic, n_roots, own, s_bar, s_core_bar)
      !! Column `ic`'s slots in the pool; `own` is 0 for an extra column
      integer, intent(in) :: ic, n_roots
      integer, intent(out) :: own, s_bar, s_core_bar

      if (ic <= n_roots) then
         own = 2 + 3*(ic - 1) + 1
         s_bar = own + 1
      else
         own = 0
         s_bar = 2 + 3*n_roots + 2*(ic - n_roots - 1) + 1
      end if
      s_core_bar = s_bar + 1
   end subroutine sa_pool_slots

   function sa_n_columns(self) result(n)
      !! Root columns, then extra columns
      class(sa_density_t), intent(in) :: self
      integer :: n
      n = 0
      if (allocated(self%columns)) n = size(self%columns)
   end function sa_n_columns

   function sa_column(self, ic) result(col)
      !! The root, or pair of roots, column `ic` belongs to
      class(sa_density_t), intent(in) :: self
      integer, intent(in) :: ic
         !! Column, 1-based
      type(column_t) :: col
      col = self%columns(ic)
   end function sa_column

   subroutine sa_energy_weighted(self, ic, matrix, error)
      !! W of column `ic`
      class(sa_density_t), intent(in) :: self
      integer, intent(in) :: ic
         !! Column, 1-based
      real(dp), allocatable, intent(out) :: matrix(:, :)
         !! (n_ao, n_ao)
      type(error_t), intent(inout) :: error
      if (.not. in_range(self, ic, error)) return
      matrix = self%w(:, :, ic)
   end subroutine sa_energy_weighted

   subroutine sa_one_particle(self, ic, matrix, error)
      !! P of column `ic`
      class(sa_density_t), intent(in) :: self
      integer, intent(in) :: ic
         !! Column, 1-based
      real(dp), allocatable, intent(out) :: matrix(:, :)
         !! (n_ao, n_ao)
      type(error_t), intent(inout) :: error
      if (.not. in_range(self, ic, error)) return
      matrix = self%p(:, :, ic)
   end subroutine sa_one_particle

   subroutine sa_density_pool(self, pool, error)
      !! The densities every column's pairs index into
      class(sa_density_t), intent(in) :: self
      real(dp), allocatable, intent(out) :: pool(:, :, :)
         !! (n_ao, n_ao, n_pool)
      type(error_t), intent(inout) :: error
      associate (unused_error => error)
      end associate
      pool = self%pool
   end subroutine sa_density_pool

   subroutine sa_separable_pairs(self, ic, pairs, error)
      !! Column `ic`'s separable pairs, as in the module header
      class(sa_density_t), intent(in) :: self
      integer, intent(in) :: ic
         !! Column, 1-based
      type(separable_pair_t), allocatable, intent(out) :: pairs(:)
      type(error_t), intent(inout) :: error

      integer, parameter :: CORE = 1, ACTIVE_SA = 2
      integer :: own, s_bar, s_core_bar, n

      if (.not. in_range(self, ic, error)) return
      call sa_pool_slots(ic, self%n_roots, own, s_bar, s_core_bar)
      allocate (pairs(MAX_PAIRS))
      n = 0
      if (own > 0) then
         call add(CORE, CORE, 0.5_dp, 0.25_dp)
         call add(CORE, own, 1.0_dp, 0.5_dp)
         call add(own, own, 0.5_dp, 0.25_dp)
      end if
      call add(CORE, s_bar, 1.0_dp, 0.5_dp)
      call add(CORE, s_core_bar, 1.0_dp, 0.5_dp)
      call add(ACTIVE_SA, s_core_bar, 1.0_dp, 0.5_dp)
      pairs = pairs(1:n)

   contains

      subroutine add(left, right, coulomb, exchange)
         integer, intent(in) :: left, right
         real(dp), intent(in) :: coulomb, exchange
         n = n + 1
         pairs(n) = separable_pair_t(left=left, right=right, coulomb=coulomb, exchange=exchange)
      end subroutine add

   end subroutine sa_separable_pairs

   subroutine sa_gamma_ao_map(self, ao_map)
      !! The AOs some active leg reaches, in compressed numbering
      class(sa_density_t), intent(in) :: self
      integer, allocatable, intent(out) :: ao_map(:)
         !! (n_ao); unallocated when no AO carries Gamma
      if (allocated(self%ao_map)) ao_map = self%ao_map
   end subroutine sa_gamma_ao_map

   subroutine sa_gamma_block(self, p_lo, p_hi, gamma, error)
      !! Every column's AO Gamma over compressed first-index positions
      !! `p_lo..p_hi`: one transform over the doubled leg `[C | Cbar]` gives
      !! the column's own two-body density in the all-`C` block and `dm2_sa`
      !! in each block with `Cbar` on exactly one leg
      class(sa_density_t), intent(in) :: self
      integer, intent(in) :: p_lo, p_hi
         !! Compressed first-index range, 1-based, inclusive
      real(dp), allocatable, intent(inout) :: gamma(:, :, :, :, :)
         !! (p_hi - p_lo + 1, n_sig, n_sig, n_sig, n_columns)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: c_both(:, :), g_both(:, :, :, :), tmp(:, :, :, :)
      integer :: jc, na, na2, n_sig, n_col

      associate (unused_error => error)
      end associate
      na = size(self%c_sig, 2)
      na2 = 2*na
      n_sig = size(self%c_sig, 1)
      n_col = self%n_columns()
      if (allocated(gamma)) deallocate (gamma)
      allocate (gamma(p_hi - p_lo + 1, n_sig, n_sig, n_sig, n_col))
      allocate (c_both(n_sig, na2), g_both(na2, na2, na2, na2))
      do jc = 1, n_col
         c_both(:, 1:na) = self%c_sig
         c_both(:, na + 1:na2) = self%cbar_sig(:, :, jc)
         g_both = 0.0_dp
         g_both(1:na, 1:na, 1:na, 1:na) = self%dm2_column(:, :, :, :, jc)
         g_both(na + 1:, 1:na, 1:na, 1:na) = self%dm2_sa
         g_both(1:na, na + 1:, 1:na, 1:na) = self%dm2_sa
         g_both(1:na, 1:na, na + 1:, 1:na) = self%dm2_sa
         g_both(1:na, 1:na, 1:na, na + 1:) = self%dm2_sa
         call ao_gamma_block(c_both, g_both, p_lo, p_hi, tmp)
         gamma(:, :, :, :, jc) = tmp
      end do
   end subroutine sa_gamma_block

   function in_range(self, ic, error) result(ok)
      !! Whether `ic` is one of this state's columns; sets `error` if not
      class(sa_density_t), intent(in) :: self
      integer, intent(in) :: ic
      type(error_t), intent(inout) :: error
      logical :: ok
      ok = ic >= 1 .and. ic <= self%n_columns()
      if (.not. ok) call error%set(ERROR_VALIDATION, "sa_density_t: no such column")
   end function in_range

end module mqc_czt_sa_effective_density
