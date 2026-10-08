!! The abstract effective-density provider, through a toy method
module test_mqc_effective_density
   !! A provider is only useful if a method can extend it and the contraction
   !! can drive it knowing nothing but the base class. `toy_density_t` is the
   !! smallest such method: two AOs, two columns -- an energy gradient and a
   !! coupling -- and a Gamma on one AO only. Every check here goes through
   !! `class(effective_density_t)`, as the contraction will.
   !!
   !! The coupling column is where the two non-deferred bindings matter: it
   !! overrides `nuclear_weight` to 0 and supplies an antisymmetric overlap
   !! term, while the gradient column keeps both defaults.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_effective_density, only: effective_density_t, separable_pair_t
   implicit none
   private

   public :: collect_mqc_effective_density_tests

   integer, parameter :: N_AO = 2
   integer, parameter :: COLUMN_GRADIENT = 1
   integer, parameter :: COLUMN_COUPLING = 2

   type, extends(effective_density_t) :: toy_density_t
      !! Two columns over two AOs, with fixed matrices
   contains
      procedure :: n_columns => toy_n_columns
      procedure :: energy_weighted => toy_energy_weighted
      procedure :: one_particle => toy_one_particle
      procedure :: density_pool => toy_density_pool
      procedure :: separable_pairs => toy_separable_pairs
      procedure :: gamma_ao_map => toy_gamma_ao_map
      procedure :: gamma_block => toy_gamma_block
      procedure :: overlap_antisymmetric => toy_overlap_antisymmetric
      procedure :: nuclear_weight => toy_nuclear_weight
   end type toy_density_t

contains

   subroutine collect_mqc_effective_density_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("a_method_extends_it_and_dispatches", test_dispatch), &
                  new_unittest("separable_pairs_index_one_shared_pool", test_pool), &
                  new_unittest("gamma_comes_in_compressed_blocks", test_gamma), &
                  new_unittest("a_coupling_overrides_the_defaults", test_coupling), &
                  new_unittest("a_column_out_of_range_is_refused", test_refuse) &
                  ]
   end subroutine collect_mqc_effective_density_tests

   ! ---- the toy method -----------------------------------------------------

   function toy_n_columns(self) result(n)
      class(toy_density_t), intent(in) :: self
      integer :: n
      associate (unused_self => self)
      end associate
      n = 2
   end function toy_n_columns

   subroutine check_column(ic, error)
      integer, intent(in) :: ic
      type(error_t), intent(inout) :: error
      if (ic < 1 .or. ic > 2) call error%set(ERROR_VALIDATION, "toy: no such column")
   end subroutine check_column

   subroutine toy_energy_weighted(self, ic, matrix, error)
      class(toy_density_t), intent(in) :: self
      integer, intent(in) :: ic
      real(dp), allocatable, intent(out) :: matrix(:, :)
      type(error_t), intent(inout) :: error
      associate (unused_self => self)
      end associate
      call check_column(ic, error)
      if (error%has_error()) return
      allocate (matrix(N_AO, N_AO))
      matrix = real(ic, dp)
   end subroutine toy_energy_weighted

   subroutine toy_one_particle(self, ic, matrix, error)
      class(toy_density_t), intent(in) :: self
      integer, intent(in) :: ic
      real(dp), allocatable, intent(out) :: matrix(:, :)
      type(error_t), intent(inout) :: error
      associate (unused_self => self)
      end associate
      call check_column(ic, error)
      if (error%has_error()) return
      allocate (matrix(N_AO, N_AO))
      matrix = 10.0_dp*real(ic, dp)
   end subroutine toy_one_particle

   subroutine toy_density_pool(self, pool, error)
      !! Two densities: one shared by both columns, one the coupling's own
      class(toy_density_t), intent(in) :: self
      real(dp), allocatable, intent(out) :: pool(:, :, :)
      type(error_t), intent(inout) :: error
      associate (unused_self => self, unused_error => error)
      end associate
      allocate (pool(N_AO, N_AO, 2))
      pool(:, :, 1) = 1.0_dp
      pool(:, :, 2) = 2.0_dp
   end subroutine toy_density_pool

   subroutine toy_separable_pairs(self, ic, pairs, error)
      class(toy_density_t), intent(in) :: self
      integer, intent(in) :: ic
      type(separable_pair_t), allocatable, intent(out) :: pairs(:)
      type(error_t), intent(inout) :: error
      associate (unused_self => self)
      end associate
      call check_column(ic, error)
      if (error%has_error()) return
      if (ic == COLUMN_GRADIENT) then
         pairs = [separable_pair_t(left=1, right=1, coulomb=0.5_dp, exchange=0.25_dp)]
      else
         pairs = [separable_pair_t(left=1, right=2, coulomb=1.0_dp, exchange=0.5_dp)]
      end if
   end subroutine toy_separable_pairs

   subroutine toy_gamma_ao_map(self, ao_map)
      !! Only the second AO carries Gamma
      class(toy_density_t), intent(in) :: self
      integer, allocatable, intent(out) :: ao_map(:)
      associate (unused_self => self)
      end associate
      ao_map = [0, 1]
   end subroutine toy_gamma_ao_map

   subroutine toy_gamma_block(self, p_lo, p_hi, gamma, error)
      class(toy_density_t), intent(in) :: self
      integer, intent(in) :: p_lo, p_hi
      real(dp), allocatable, intent(inout) :: gamma(:, :, :, :, :)
      type(error_t), intent(inout) :: error
      integer :: ic
      associate (unused_self => self, unused_error => error)
      end associate
      if (allocated(gamma)) deallocate (gamma)
      allocate (gamma(p_hi - p_lo + 1, 1, 1, 1, 2))
      do ic = 1, 2
         gamma(:, :, :, :, ic) = 100.0_dp*real(ic, dp)
      end do
   end subroutine toy_gamma_block

   subroutine toy_overlap_antisymmetric(self, ic, matrix, error)
      !! The coupling column's antisymmetric term; the gradient has none
      class(toy_density_t), intent(in) :: self
      integer, intent(in) :: ic
      real(dp), allocatable, intent(out) :: matrix(:, :)
      type(error_t), intent(inout) :: error
      associate (unused_self => self)
      end associate
      call check_column(ic, error)
      if (error%has_error()) return
      if (ic /= COLUMN_COUPLING) return
      allocate (matrix(N_AO, N_AO))
      matrix = reshape([0.0_dp, -1.0_dp, 1.0_dp, 0.0_dp], [N_AO, N_AO])
   end subroutine toy_overlap_antisymmetric

   function toy_nuclear_weight(self, ic) result(weight)
      class(toy_density_t), intent(in) :: self
      integer, intent(in) :: ic
      real(dp) :: weight
      associate (unused_self => self)
      end associate
      weight = 1.0_dp
      if (ic == COLUMN_COUPLING) weight = 0.0_dp
   end function toy_nuclear_weight

   ! ---- the tests: everything through the base class -----------------------

   subroutine test_dispatch(error)
      type(error_type), allocatable, intent(out) :: error
      class(effective_density_t), allocatable :: provider
      type(error_t) :: err
      real(dp), allocatable :: w(:, :), p(:, :)

      allocate (toy_density_t :: provider)
      call check(error, provider%n_columns() == 2, "two columns")
      if (allocated(error)) return

      call provider%energy_weighted(2, w, err)
      call provider%one_particle(2, p, err)
      call check(error,.not. err%has_error(), "no error for a valid column")
      if (allocated(error)) return
      call check(error, all(shape(w) == [N_AO, N_AO]) .and. all(w == 2.0_dp), &
                 "W of column 2 through the base class")
      if (allocated(error)) return
      call check(error, all(p == 20.0_dp), "P of column 2 through the base class")
   end subroutine test_dispatch

   subroutine test_pool(error)
      type(error_type), allocatable, intent(out) :: error
      class(effective_density_t), allocatable :: provider
      type(error_t) :: err
      real(dp), allocatable :: pool(:, :, :)
      type(separable_pair_t), allocatable :: pairs(:)
      integer :: ic, k

      allocate (toy_density_t :: provider)
      call provider%density_pool(pool, err)
      call check(error, size(pool, 3) == 2, "two pooled densities")
      if (allocated(error)) return
      do ic = 1, provider%n_columns()
         call provider%separable_pairs(ic, pairs, err)
         do k = 1, size(pairs)
            call check(error, pairs(k)%left >= 1 .and. pairs(k)%left <= size(pool, 3) .and. &
                       pairs(k)%right >= 1 .and. pairs(k)%right <= size(pool, 3), &
                       "every pair indexes into the shared pool")
            if (allocated(error)) return
         end do
      end do
      ! Density 1 is used by both columns and held once.
      call provider%separable_pairs(COLUMN_COUPLING, pairs, err)
      call check(error, pairs(1)%left == 1, "the coupling reuses the gradient's density")
   end subroutine test_pool

   subroutine test_gamma(error)
      type(error_type), allocatable, intent(out) :: error
      class(effective_density_t), allocatable :: provider
      type(error_t) :: err
      integer, allocatable :: ao_map(:)
      real(dp), allocatable :: gamma(:, :, :, :, :)
      integer :: n_sig

      allocate (toy_density_t :: provider)
      call provider%gamma_ao_map(ao_map)
      n_sig = count(ao_map > 0)
      call check(error, size(ao_map) == N_AO .and. n_sig == 1, "one AO carries Gamma")
      if (allocated(error)) return
      call provider%gamma_block(1, n_sig, gamma, err)
      call check(error, all(shape(gamma) == [1, 1, 1, 1, 2]), &
                 "one compressed block, every column at once")
      if (allocated(error)) return
      call check(error, gamma(1, 1, 1, 1, 2) == 200.0_dp, "column 2's Gamma")
   end subroutine test_gamma

   subroutine test_coupling(error)
      type(error_type), allocatable, intent(out) :: error
      class(effective_density_t), allocatable :: provider
      type(error_t) :: err
      real(dp), allocatable :: anti(:, :)

      allocate (toy_density_t :: provider)
      call check(error, provider%nuclear_weight(COLUMN_GRADIENT) == 1.0_dp, &
                 "a gradient column carries the nuclear repulsion")
      if (allocated(error)) return
      call check(error, provider%nuclear_weight(COLUMN_COUPLING) == 0.0_dp, &
                 "a coupling column does not")
      if (allocated(error)) return

      call provider%overlap_antisymmetric(COLUMN_GRADIENT, anti, err)
      call check(error,.not. allocated(anti), "a gradient column has no antisymmetric term")
      if (allocated(error)) return
      call provider%overlap_antisymmetric(COLUMN_COUPLING, anti, err)
      call check(error, allocated(anti), "a coupling column supplies one")
      if (allocated(error)) return
      call check(error, all(anti == -transpose(anti)), "and it is antisymmetric")
   end subroutine test_coupling

   subroutine test_refuse(error)
      type(error_type), allocatable, intent(out) :: error
      class(effective_density_t), allocatable :: provider
      type(error_t) :: err
      real(dp), allocatable :: w(:, :)

      allocate (toy_density_t :: provider)
      call provider%energy_weighted(3, w, err)
      call check(error, err%has_error(), "a column the provider does not have is an error")
      if (allocated(error)) return
      call check(error,.not. allocated(w), "and returns no matrix")
   end subroutine test_refuse

end module test_mqc_effective_density

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_effective_density, only: collect_mqc_effective_density_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_effective_density", collect_mqc_effective_density_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
