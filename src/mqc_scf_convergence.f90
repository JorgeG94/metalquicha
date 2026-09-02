!! One place that decides whether an SCF has converged
module mqc_scf_convergence
   !! **What stops an SCF, decided in one place instead of at each loop.**
   !!
   !! The test used to be written out at every SCF loop, as
   !! `de < energy_tol .and. gnorm < gtol`, with `gtol` derived from
   !! `sqrt(energy_tol)` a few dozen lines earlier and overridden by an
   !! optional argument. Written twice in `mqc_libcint_rhf` alone, once more
   !! for each cuEST loop, and with a fifth variant inside the fragment
   !! potential builder that quietly replaced the derived threshold with the
   !! density one. Changing the rule meant finding all of them, and the last
   !! attempt did not: the CPU and GPU paths ended up stopping on different
   !! criteria.
   !!
   !! It also made a deck unpredictable. A deck asking for `tolerance: 1e-6`
   !! got an energy gate at 1e-6 and a commutator gate at 1e-3, and where the
   !! second was the binding one the run stopped against a number the deck
   !! never wrote and the iteration table never printed.
   !!
   !! So the rule is a value now. `scf_convergence_t` carries which metric
   !! decides and the threshold it is measured against, `is_converged` is the
   !! only place the comparison happens, and `describe` is what the iteration
   !! table prints so the active rule is never invisible again.
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
   !! * `CONV_METRIC_ENERGY` is the weakest, and unsafe with an *interpolating*
   !!   accelerator. EDIIS on water/6-31G holds `dE` at 1e-13 while the
   !!   commutator sits at 1.1e-2, and an energy-only test stops there --
   !!   5.8e-5 hartree from the answer, measured. Selecting it is a deliberate
   !!   choice and warned about; it is not a default.
   !! * `CONV_METRIC_STANDARD` is both, and is today's behaviour exactly. It is
   !!   the default here so that introducing this type changes no result.
   !!
   !! ### The transitional default
   !!
   !! `CONV_METRIC_STANDARD` derives its commutator bound as `sqrt(tolerance)`,
   !! which is the derivation this type exists to remove. It is kept for one
   !! release so the machinery can land without moving any number, and because
   !! 344 of the 355 validation decks name `tolerance` explicitly -- 288 of them
   !! at 1e-12, which read as a commutator bound would be below the OpenMP
   !! reduction noise floor and unreachable. Retiring it means migrating those
   !! decks by `sqrt`, which is exact rather than approximate, and belongs in
   !! its own change.
   use pic_types, only: dp
   implicit none
   private

   public :: scf_convergence_t
   public :: parse_convergence_metric, convergence_metric_name
   public :: CONV_METRIC_UNKNOWN, CONV_METRIC_STANDARD
   public :: CONV_METRIC_ENERGY, CONV_METRIC_COMMUTATOR, CONV_METRIC_DENSITY

   integer, parameter :: CONV_METRIC_UNKNOWN = 0
      !! Returned by `parse_convergence_metric` for a spelling it does not know.
   integer, parameter :: CONV_METRIC_STANDARD = 1
      !! Energy and commutator together, the commutator bound derived as
      !! `sqrt(tolerance)`. Today's behaviour, and the default.
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
         !! consulted under `CONV_METRIC_STANDARD`: naming a metric already
         !! states which quantity is being bounded, so a second threshold for a
         !! different one would be two rules again.
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
      !! Every SCF loop calls this and none of them decides anything. All three
      !! quantities are passed whether the metric reads them or not: an
      !! argument list that changed with the metric would put a branch back in
      !! the caller, which is the thing being removed.
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
         ! CONV_METRIC_STANDARD, and anything unset -- both gates, which is
         ! what every loop in this program tested before this type existed.
         done = de < this%tolerance .and. gnorm < this%commutator_bound()
      end select
   end function is_converged

   pure function describe(this) result(text)
      !! One line naming the rule in force, for the iteration table
      !!
      !! The complaint that produced this type was not that the threshold was
      !! wrong but that it was invisible: the run stopped against a number the
      !! deck never wrote. Printing it costs one line and removes the whole
      !! class of confusion.
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
      !! `diis` and `gradient` are accepted for the commutator because all
      !! three name the same quantity -- `FDS - SDF` in the orthogonal basis is
      !! both what DIIS extrapolates against and the orbital gradient. A user
      !! reaching for any of the three means the same thing, and refusing two
      !! of them would be pedantry rather than precision.
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
