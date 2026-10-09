!! The effective densities an analytic gradient or coupling is contracted from
module mqc_czt_effective_density
   !! An analytic nuclear derivative -- of one state's energy, or of the
   !! coupling between two states -- is the contraction of method-specific
   !! effective densities with derivative integrals that know nothing of the
   !! method:
   !!
   !!     dE/dx = tr(W S^x) + tr(P h^x) + sum_pairs D_l . (J^x, K^x)(D_r)
   !!           + tr(Gamma (mn|ls)^x) + ...
   !!
   !! `effective_density_t` is what a method hands that contraction: the
   !! objects on the left of each trace, for one or more *columns*. A column
   !! is one energy gradient or one coupling. A method extends this type and
   !! fills it from its own relaxed densities; the contraction sees only the
   !! deferred bindings below.
   !!
   !! With exact four-index integrals a column carries five objects:
   !!
   !! 1. `energy_weighted`: W, against the overlap derivative
   !! 2. `one_particle`: P, against the core-Hamiltonian derivative
   !! 3. `separable_pairs`: the separable two-particle part, as pairs of
   !!    densities from `density_pool` with Coulomb and exchange weights
   !! 4. `gamma_block`: the non-separable two-particle density in the AO
   !!    basis, one block of its first index at a time; absent for a method
   !!    whose two-particle part is all separable, such as Hartree-Fock
   !! 5. `overlap_antisymmetric`: the antisymmetric overlap term a coupling
   !!    carries; absent for an energy gradient
   !!
   !! `column` says which state a gradient column belongs to, or which pair of
   !! states a coupling column couples, so a caller can file each contracted
   !! column under the right root or pair.
   !!
   !! Only traces against these objects are described here. A term that is not
   !! one -- the nuclear repulsion, which a gradient column has and a coupling
   !! does not, exchange-correlation, dispersion, the explicit terms of a
   !! continuum solvent or an ECP -- is added by the code that assembles the
   !! gradient. A term the method adds to its Lagrangian, such as a solvent's
   !! response, is already in the densities it hands over.
   !!
   !! The density-fitted form, which replaces object 4 with a three-index and a
   !! two-index density, is not defined here yet.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_GENERIC
   implicit none
   private

   public :: effective_density_t
   public :: separable_pair_t
   public :: column_t
   public :: COLUMN_KIND_GRADIENT, COLUMN_KIND_COUPLING

   integer, parameter :: COLUMN_KIND_GRADIENT = 1
      !! The nuclear gradient of one state's energy
   integer, parameter :: COLUMN_KIND_COUPLING = 2
      !! The nonadiabatic coupling between two states

   type :: column_t
      !! What one column is the derivative of
      integer :: kind = COLUMN_KIND_GRADIENT
         !! `COLUMN_KIND_GRADIENT` or `COLUMN_KIND_COUPLING`
      integer :: bra = 1
         !! The state, 1-based; the first state of a coupling
      integer :: ket = 1
         !! Equal to `bra` for a gradient; the second state of a coupling,
         !! in the `[state_i, state_j]` order of `keywords.mcscf.nac_pairs`
   end type column_t

   type :: separable_pair_t
      !! One separable two-particle term of a column:
      !! `coulomb * D_left . J^x(D_right) - exchange * D_left . K^x(D_right)`,
      !! with both densities taken from the provider's `density_pool`.
      !!
      !! Spin enters through the pool: an unrestricted method pools the total,
      !! alpha and beta densities, pairs the total with itself for Coulomb and
      !! each spin density with itself for exchange. The pool is symmetric, so
      !! a method whose density is not (a coupled-cluster response density, a
      !! transition density) pools its symmetric part. That is exact whenever
      !! the other density of the pair is symmetric, because the antisymmetric
      !! part then drops out of both `J^x` and `K^x`.
      integer :: left = 0
         !! Index of the left density in `density_pool`, 1-based
      integer :: right = 0
         !! Index of the right density in `density_pool`, 1-based
      real(dp) :: coulomb = 0.0_dp
         !! Weight of the Coulomb derivative term
      real(dp) :: exchange = 0.0_dp
         !! Weight of the exchange derivative term
   end type separable_pair_t
   ! TODO(mqc): there is no attenuated-exchange weight, so the second K^x pass
   ! a range-separated hybrid makes at the screened omega cannot be expressed;
   ! a Kohn-Sham provider for such a functional would miss that term.

   type, abstract :: effective_density_t
      !! The effective densities of `n_columns()` gradient or coupling columns,
      !! all in the AO basis of one molecule
   contains
      procedure(column_count_i), deferred :: n_columns
         !! How many columns this provider carries
      procedure(column_describe_i), deferred :: column
         !! Which state, or pair of states, one column belongs to
      procedure(column_matrix_i), deferred :: energy_weighted
         !! Object 1: W of one column, (n_ao, n_ao)
      procedure(column_matrix_i), deferred :: one_particle
         !! Object 2: P of one column, (n_ao, n_ao)
      procedure(density_pool_i), deferred :: density_pool
         !! The densities the separable pairs of every column index into,
         !! (n_ao, n_ao, n_pool). Shared, so a density common to several
         !! columns is held and contracted once.
      procedure(separable_pairs_i), deferred :: separable_pairs
         !! Object 3: the separable pairs of one column
      procedure :: gamma_ao_map => no_gamma_ao_map
         !! Which AOs carry any column's Gamma, as the compressed index
         !! `gamma_block` is written in. Unallocated unless overridden.
      procedure :: gamma_block => no_gamma_block
         !! Object 4: every column's Gamma over one first-index block. An
         !! error unless overridden, since there is no block to ask for.
      procedure :: overlap_antisymmetric => no_antisymmetric_overlap
         !! Object 5: unallocated unless overridden
   end type effective_density_t

   abstract interface
      function column_count_i(self) result(n)
         !! The number of columns
         import :: effective_density_t
         implicit none
         class(effective_density_t), intent(in) :: self
         integer :: n
      end function column_count_i

      function column_describe_i(self, ic) result(col)
         !! The state, or pair of states, column `ic` belongs to
         import :: effective_density_t, column_t
         implicit none
         class(effective_density_t), intent(in) :: self
         integer, intent(in) :: ic
            !! Column, 1-based
         type(column_t) :: col
      end function column_describe_i

      subroutine column_matrix_i(self, ic, matrix, error)
         !! One AO matrix of column `ic`
         import :: effective_density_t, dp, error_t
         implicit none
         class(effective_density_t), intent(in) :: self
         integer, intent(in) :: ic
            !! Column, 1-based
         real(dp), allocatable, intent(out) :: matrix(:, :)
            !! (n_ao, n_ao)
         type(error_t), intent(inout) :: error
      end subroutine column_matrix_i

      subroutine density_pool_i(self, pool, error)
         !! The shared density pool
         import :: effective_density_t, dp, error_t
         implicit none
         class(effective_density_t), intent(in) :: self
         real(dp), allocatable, intent(out) :: pool(:, :, :)
            !! (n_ao, n_ao, n_pool), each symmetric
         type(error_t), intent(inout) :: error
      end subroutine density_pool_i

      subroutine separable_pairs_i(self, ic, pairs, error)
         !! The separable pairs of column `ic`
         import :: effective_density_t, separable_pair_t, error_t
         implicit none
         class(effective_density_t), intent(in) :: self
         integer, intent(in) :: ic
            !! Column, 1-based
         type(separable_pair_t), allocatable, intent(out) :: pairs(:)
            !! Empty when the column has no separable part
         type(error_t), intent(inout) :: error
      end subroutine separable_pairs_i

      subroutine gamma_ao_map_i(self, ao_map)
         !! The compressed AO numbering of `gamma_block`
         import :: effective_density_t
         implicit none
         class(effective_density_t), intent(in) :: self
         integer, allocatable, intent(out) :: ao_map(:)
            !! (n_ao); position of each AO in the compressed numbering, or 0
            !! for one that carries no Gamma. Unallocated or all zero: no
            !! column has Gamma, and `gamma_block` must not be called.
      end subroutine gamma_ao_map_i

      subroutine gamma_block_i(self, p_lo, p_hi, gamma, error)
         !! Every column's AO Gamma over compressed first-index positions
         !! `p_lo..p_hi`
         import :: effective_density_t, dp, error_t
         implicit none
         class(effective_density_t), intent(in) :: self
         integer, intent(in) :: p_lo, p_hi
            !! Compressed first-index range, 1-based, inclusive
         real(dp), allocatable, intent(inout) :: gamma(:, :, :, :, :)
            !! (p_hi - p_lo + 1, n_sig, n_sig, n_sig, n_columns)
         type(error_t), intent(inout) :: error
      end subroutine gamma_block_i
   end interface

contains

   subroutine no_antisymmetric_overlap(self, ic, matrix, error)
      !! No antisymmetric overlap term: `matrix` comes back unallocated
      class(effective_density_t), intent(in) :: self
      integer, intent(in) :: ic
         !! Column, 1-based
      real(dp), allocatable, intent(out) :: matrix(:, :)
         !! Unallocated; (n_ao, n_ao) where overridden
      type(error_t), intent(inout) :: error
      ! Arguments a default must accept and does not use.
      associate (unused_self => self, unused_ic => ic, unused_error => error)
      end associate
   end subroutine no_antisymmetric_overlap

   subroutine no_gamma_ao_map(self, ao_map)
      !! No non-separable two-particle density: `ao_map` comes back unallocated
      class(effective_density_t), intent(in) :: self
      integer, allocatable, intent(out) :: ao_map(:)
         !! Unallocated; (n_ao) where overridden
      associate (unused_self => self)
      end associate
   end subroutine no_gamma_ao_map

   subroutine no_gamma_block(self, p_lo, p_hi, gamma, error)
      !! Refuses: a provider with no Gamma has no block to fill, and the empty
      !! `gamma_ao_map` already told the caller not to ask
      class(effective_density_t), intent(in) :: self
      integer, intent(in) :: p_lo, p_hi
         !! Compressed first-index range, 1-based, inclusive
      real(dp), allocatable, intent(inout) :: gamma(:, :, :, :, :)
         !! Deallocated, so nothing stale reads as a filled block
      type(error_t), intent(inout) :: error
      associate (unused_self => self, unused_lo => p_lo, unused_hi => p_hi)
      end associate
      if (allocated(gamma)) deallocate (gamma)
      call error%set(ERROR_GENERIC, "gamma_block called on an effective density with "// &
                     "no non-separable Gamma; check gamma_ao_map first")
   end subroutine no_gamma_block

end module mqc_czt_effective_density
