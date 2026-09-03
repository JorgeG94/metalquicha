!! One place that decides whether an SCF has converged
module mqc_scf_convergence
   !! What stops an SCF, as a value rather than a test spelled out at each
   !! loop. `scf_convergence_t` carries which metric decides and the threshold
   !! it is measured against, `is_converged` is the only place the comparison
   !! happens, and `describe` is what the iteration table prints.
   !!
   !! ### The metrics are not interchangeable
   !!
   !! Near convergence the energy is variational, so its error falls as the
   !! *square* of the commutator while the density's falls linearly. That makes
   !! the commutator the strictly stronger measure and the energy the weakest:
   !!
   !! * `CONV_METRIC_COMMUTATOR` is what anything reading the density wants --
   !!   multipoles, atomic charges, a fragment potential, any finite difference
   !!   of one of those.
   !! * `CONV_METRIC_ENERGY` is unsafe with an *interpolating* accelerator:
   !!   EDIIS on water/6-31G holds `dE` at 1e-13 while the commutator sits at
   !!   1.1e-2, and an energy-only test stops there. Selecting it is warned
   !!   about.
   !! * `CONV_METRIC_STANDARD` is both, and the default.
   use pic_types, only: dp
   implicit none
   private

   public :: scf_convergence_t
   public :: parse_convergence_metric, convergence_metric_name
   public :: CONV_METRIC_UNKNOWN, CONV_METRIC_STANDARD
   public :: CONV_METRIC_ENERGY, CONV_METRIC_COMMUTATOR, CONV_METRIC_DENSITY

   integer, parameter :: CONV_METRIC_UNKNOWN = 0
      !! Returned by `parse_convergence_metric` for a spelling it does not know.
   ! TODO(mqc): the default still derives its commutator bound as
   ! `sqrt(tolerance)`, which is the derivation this type exists to remove. A
   ! deck asking for `tolerance: 1e-6` is stopped by a commutator gate at 1e-3
   ! that it never wrote. Retiring it means migrating the validation decks'
   ! `tolerance` by `sqrt`.
   integer, parameter :: CONV_METRIC_STANDARD = 1
      !! Energy and commutator together, the commutator bound derived as
      !! `sqrt(tolerance)`. The default.
   integer, parameter :: CONV_METRIC_ENERGY = 2
      !! The change in energy between iterations, alone.
   integer, parameter :: CONV_METRIC_COMMUTATOR = 3
      !! `FDS - SDF` in the orthogonal basis, as a max element. The same
      !! quantity DIIS extrapolates against and the same one pyscf calls
      !! `norm_gorb`, which is why `diis` and `gradient` are accepted spellings.
   integer, parameter :: CONV_METRIC_DENSITY = 4
      !! RMS change in the density matrix between iterations.

   type :: scf_convergence_t
      !! Which measure decides an SCF has converged, and at what threshold
      integer :: metric = CONV_METRIC_STANDARD
      real(dp) :: tolerance = 1.0e-8_dp
         !! In the units of `metric`. Under `CONV_METRIC_STANDARD` that is an
         !! energy, and the commutator bound comes from it by `sqrt`.
      real(dp) :: gradient_tolerance = 0.0_dp
         !! `keywords.scf.gradient_tolerance`. Zero means derive it. Only
         !! consulted under `CONV_METRIC_STANDARD` -- naming any other metric
         !! already states which quantity is being bounded.
   contains
      procedure :: is_converged
      procedure :: describe
      procedure :: commutator_bound
   end type scf_convergence_t

contains

   pure function commutator_bound(this) result(bound)
      !! The commutator threshold in force, derived or stated
      class(scf_convergence_t), intent(in) :: this
      real(dp) :: bound

      if (this%metric == CONV_METRIC_COMMUTATOR) then
         bound = this%tolerance
      else if (this%gradient_tolerance > 0.0_dp) then
         bound = this%gradient_tolerance
      else
         bound = sqrt(this%tolerance)
      end if
   end function commutator_bound

   pure function is_converged(this, de, drms, gnorm) result(done)
      !! The one comparison
      !!
      !! All three quantities are passed whether the metric reads them or not,
      !! so no caller has to branch on the metric.
      class(scf_convergence_t), intent(in) :: this
      real(dp), intent(in) :: de
         !! |E(iter) - E(iter-1)|
      real(dp), intent(in) :: drms
         !! RMS change in the density matrix
      real(dp), intent(in) :: gnorm
         !! max |FDS - SDF| in the orthogonal basis
      logical :: done

      select case (this%metric)
      case (CONV_METRIC_ENERGY)
         done = de < this%tolerance
      case (CONV_METRIC_COMMUTATOR)
         done = gnorm < this%tolerance
      case (CONV_METRIC_DENSITY)
         done = drms < this%tolerance
      case default
         ! CONV_METRIC_STANDARD, and anything unset -- both gates.
         done = de < this%tolerance .and. gnorm < this%commutator_bound()
      end select
   end function is_converged

   pure function describe(this) result(text)
      !! One line naming the rule in force, for the iteration table
      !!
      !! Under `CONV_METRIC_STANDARD` it names both gates, including the
      !! derived commutator bound the deck never wrote.
      class(scf_convergence_t), intent(in) :: this
      character(len=:), allocatable :: text

      character(len=32) :: a, b

      write (a, "(es9.2)") this%tolerance
      select case (this%metric)
      case (CONV_METRIC_ENERGY)
         text = "dE < "//trim(adjustl(a))
      case (CONV_METRIC_COMMUTATOR)
         text = "|FDS-SDF| < "//trim(adjustl(a))
      case (CONV_METRIC_DENSITY)
         text = "dD(rms) < "//trim(adjustl(a))
      case default
         write (b, "(es9.2)") this%commutator_bound()
         text = "dE < "//trim(adjustl(a))//" and |FDS-SDF| < "//trim(adjustl(b))
      end select
   end function describe

   pure subroutine parse_convergence_metric(name, metric, ok)
      !! `keywords.scf.convergence_metric` to a `CONV_METRIC_*` kind
      !!
      !! `diis` and `gradient` are accepted for the commutator: `FDS - SDF` in
      !! the orthogonal basis is both what DIIS extrapolates against and the
      !! orbital gradient.
      character(len=*), intent(in) :: name
      integer, intent(out) :: metric
      logical, intent(out) :: ok

      character(len=len(name)) :: lower
      integer :: i, c

      lower = name
      do i = 1, len_trim(lower)
         c = iachar(lower(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) lower(i:i) = achar(c + 32)
      end do

      ok = .true.
      select case (trim(lower))
      case ("standard", "energy_and_commutator")
         metric = CONV_METRIC_STANDARD
      case ("energy")
         metric = CONV_METRIC_ENERGY
      case ("commutator", "diis", "gradient")
         metric = CONV_METRIC_COMMUTATOR
      case ("density")
         metric = CONV_METRIC_DENSITY
      case default
         metric = CONV_METRIC_UNKNOWN
         ok = .false.
      end select
   end subroutine parse_convergence_metric

   pure function convergence_metric_name(metric) result(name)
      !! The canonical spelling, for a message or the config echo
      integer, intent(in) :: metric
      character(len=:), allocatable :: name

      select case (metric)
      case (CONV_METRIC_STANDARD)
         name = "standard"
      case (CONV_METRIC_ENERGY)
         name = "energy"
      case (CONV_METRIC_COMMUTATOR)
         name = "commutator"
      case (CONV_METRIC_DENSITY)
         name = "density"
      case default
         name = "unknown"
      end select
   end function convergence_metric_name

end module mqc_scf_convergence
