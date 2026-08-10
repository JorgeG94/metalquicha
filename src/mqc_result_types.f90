!! Quantum chemistry calculation result containers
module mqc_result_types
   !! Defines data structures for storing and managing results from
   !! quantum chemistry calculations including energies, gradients, and properties.
   use pic_types, only: dp, int32
   use pic_mpi_lib, only: comm_t, isend, irecv, send, recv, wait, request_t, MPI_Status
   use mqc_error, only: error_t
   implicit none
   private

   public :: mp2_energy_t          !! MP2 energy components type
   public :: cc_energy_t           !! Coupled cluster energy components type
   public :: energy_t              !! Energy components type
   public :: calculation_result_t  !! Main result container type
   public :: SCF_UNKNOWN, SCF_CONVERGED, SCF_NOT_CONVERGED
   public :: frontier_orbitals
   public :: scf_status_label

   integer, parameter :: SCF_UNKNOWN = 0
      !! The method did not say. Not a claim that it converged.
   integer, parameter :: SCF_CONVERGED = 1
      !! Reached its energy and density thresholds
   integer, parameter :: SCF_NOT_CONVERGED = 2
      !! Ran out of iterations. The energy is whatever the last cycle held,
      !! and putting it into an expansion poisons the total silently -- it is
      !! a number of the right magnitude, so nothing downstream notices.
   public :: mbe_result_t          !! MBE aggregated result container type
   public :: result_send, result_isend  !! Send result over MPI
   public :: result_recv, result_irecv  !! Receive result over MPI

   ! SCS-MP2 scaling parameters
   real(dp), parameter :: SCS_SS_SCALE = 1.0_dp/3.0_dp  !! SCS same-spin scaling factor
   real(dp), parameter :: SCS_OS_SCALE = 1.2_dp         !! SCS opposite-spin scaling factor

   type :: mp2_energy_t
      !! Container for MP2 energy components (SS/OS)
      real(dp) :: ss = 0.0_dp  !! Same-spin correlation energy (Hartree)
      real(dp) :: os = 0.0_dp  !! Opposite-spin correlation energy (Hartree)
   contains
      procedure :: total => mp2_total           !! Compute total MP2 correlation
      procedure :: scs => mp2_scs               !! Compute SCS-MP2 correlation
      procedure :: reset => mp2_reset           !! Reset both components to zero
      procedure :: check_stability => mp2_check_stability  !! Check for positive energies (instability)
   end type mp2_energy_t

   type :: cc_energy_t
      !! Container for coupled cluster energy components
      real(dp) :: singles = 0.0_dp   !! Singles contribution (Hartree)
      real(dp) :: doubles = 0.0_dp   !! Doubles contribution (Hartree)
      real(dp) :: triples = 0.0_dp   !! Triples contribution (Hartree)
   contains
      procedure :: total => cc_total            !! Compute total CC correlation
      procedure :: reset => cc_reset            !! Reset all components to zero
      procedure :: check_stability => cc_check_stability  !! Check for positive energies (instability)
   end type cc_energy_t

   type :: energy_t
      !! Container for quantum chemistry energy components
      !!
      !! Stores energy contributions from different levels of theory.
      !! Total energy is computed as: scf + mp2%total() + cc%total()
      real(dp) :: scf = 0.0_dp           !! SCF/HF reference energy (Hartree)
      type(mp2_energy_t) :: mp2          !! MP2 correlation components
      type(cc_energy_t) :: cc            !! Coupled cluster correlation components
      ! add more as needed, also need to modify the total energy function
   contains
      procedure :: total => energy_total     !! Compute total energy from components
      procedure :: reset => energy_reset     !! Reset all components to zero
   end type energy_t

   type :: calculation_result_t
      !! Container for quantum chemistry calculation results
      !!
      !! Stores computed quantities from QC calculations with flags
      !! indicating which properties have been computed.
      type(energy_t) :: energy                  !! Energy components (Hartree)
      real(dp), allocatable :: gradient(:, :)   !! Energy gradient (3, natoms) (Hartree/Bohr)
      real(dp), allocatable :: sigma(:, :)     !! Stress tensor (3,3) (Hartree/Bohr^3)
      real(dp), allocatable :: hessian(:, :)    !! Energy hessian (future implementation)
      real(dp), allocatable :: dipole(:)        !! Dipole moment vector (3) (Debye)
      real(dp), allocatable :: dipole_derivatives(:, :)  !! Dipole derivatives (3, 3N) in a.u. for IR intensities

      ! Fragment metadata
      real(dp) :: distance = 0.0_dp      !! Minimal atomic distance between monomers (Angstrom, 0 for monomers)

      ! Computation status flags
      logical :: has_energy = .false.    !! Energy has been computed
      logical :: has_gradient = .false.  !! Gradient has been computed
      logical :: has_sigma = .false.     !! Stress tensor has been computed
      logical :: has_hessian = .false.   !! Hessian has been computed
      logical :: has_dipole = .false.    !! Dipole moment has been computed
      logical :: has_dipole_derivatives = .false.  !! Dipole derivatives have been computed

      ! Frontier orbitals
      !
      ! HOMO and LUMO rather than the whole spectrum: the gap is what anyone
      ! asks for, and a per-fragment orbital array would be the largest thing
      ! in the result by a wide margin.
      !
      ! **These do not add up.** Energies and dipoles are additive and the
      ! expansion sums them; a gap is not, and there is no combination of
      ! fragment gaps that is the system's. So these stay per fragment and no
      ! total is formed anywhere -- a summed gap would look like a number and
      ! be nothing at all. For the gap of a molecule, run it unfragmented.
      real(dp) :: homo = 0.0_dp          !! Highest occupied orbital energy (Hartree)
      real(dp) :: lumo = 0.0_dp          !! Lowest unoccupied orbital energy (Hartree)
      logical :: has_orbitals = .false.  !! Whether homo/lumo were reported

      ! SCF convergence
      !
      ! Tri-state rather than a logical, because either default would be a
      ! lie: `.true.` makes a method that never reports look converged, and
      ! `.false.` makes it look failed. A method that does not say leaves this
      ! at SCF_UNKNOWN and the output prints "?" -- which is the truth, and
      ! which shows up as something to fix rather than as a clean bill.
      integer :: scf_status = SCF_UNKNOWN  !! Whether the SCF reached its thresholds
      integer :: scf_iterations = 0        !! Cycles used, 0 if not reported

      ! Error handling
      type(error_t) :: error             !! Calculation error (if any)
      logical :: has_error = .false.     !! True if calculation failed
   contains
      procedure :: destroy => result_destroy  !! Clean up allocated memory
      procedure :: reset => result_reset      !! Reset all values and flags
   end type calculation_result_t

   type :: mbe_result_t
      !! Container for Many-Body Expansion aggregated results
      !!
      !! Stores total properties computed via MBE: energy, gradient, hessian, dipole.
      !! Caller allocates desired components before calling compute_mbe; the function
      !! uses allocated() to determine what to compute and sets has_* flags on success.

      real(dp) :: total_energy = 0.0_dp              !! Total MBE energy (Hartree)
      real(dp), allocatable :: gradient(:, :)        !! Total gradient (3, total_atoms) (Hartree/Bohr)
      real(dp), allocatable :: hessian(:, :)         !! Total Hessian (3*natoms, 3*natoms)
      real(dp), allocatable :: dipole(:)             !! Total dipole moment (3) (e*Bohr)
      real(dp), allocatable :: dipole_derivatives(:, :)  !! Dipole derivatives (3, 3*natoms) for IR intensities

      ! Computation status flags
      logical :: has_energy = .false.                !! Energy has been computed
      logical :: has_gradient = .false.              !! Gradient has been computed
      logical :: has_hessian = .false.               !! Hessian has been computed
      logical :: has_dipole = .false.                !! Dipole has been computed
      logical :: has_dipole_derivatives = .false.    !! Dipole derivatives have been computed
   contains
      procedure :: destroy => mbe_result_destroy            !! Clean up allocated memory
      procedure :: reset => mbe_result_reset                !! Reset all values and flags
      procedure :: allocate_gradient => mbe_result_allocate_gradient
      procedure :: allocate_hessian => mbe_result_allocate_hessian
      procedure :: allocate_dipole => mbe_result_allocate_dipole
   end type mbe_result_t

contains

   pure function mp2_total(this) result(total)
      !! Compute total MP2 correlation energy
      class(mp2_energy_t), intent(in) :: this
      real(dp) :: total

      total = this%ss + this%os
   end function mp2_total

   pure function mp2_scs(this) result(scs_energy)
      !! Compute SCS-MP2 (Spin-Component Scaled MP2) correlation energy
      !! SCS-MP2 uses: E_SCS = (1/3)*E_SS + 1.2*E_OS
      class(mp2_energy_t), intent(in) :: this
      real(dp) :: scs_energy

      scs_energy = SCS_SS_SCALE*this%ss + SCS_OS_SCALE*this%os
   end function mp2_scs

   subroutine mp2_reset(this)
      !! Reset both MP2 components to zero
      class(mp2_energy_t), intent(inout) :: this
      this%ss = 0.0_dp
      this%os = 0.0_dp
   end subroutine mp2_reset

   subroutine mp2_check_stability(this)
      !! Check for positive MP2 correlation energies (instability warning)
      !! Correlation energies should be negative; positive values indicate instability
      use pic_logger, only: logger => global_logger
      class(mp2_energy_t), intent(in) :: this

      if (this%ss > 0.0_dp) then
         call logger%warning("MP2 same-spin correlation energy is positive - possible instability!")
      end if

      if (this%os > 0.0_dp) then
         call logger%warning("MP2 opposite-spin correlation energy is positive - possible instability!")
      end if
   end subroutine mp2_check_stability

   pure function cc_total(this) result(total)
      !! Compute total CC correlation energy
      class(cc_energy_t), intent(in) :: this
      real(dp) :: total

      total = this%singles + this%doubles + this%triples
   end function cc_total

   subroutine cc_reset(this)
      !! Reset all CC components to zero
      class(cc_energy_t), intent(inout) :: this
      this%singles = 0.0_dp
      this%doubles = 0.0_dp
      this%triples = 0.0_dp
   end subroutine cc_reset

   subroutine cc_check_stability(this)
      !! Check for positive CC correlation energies (instability warning)
      !! Correlation energies should be negative; positive values indicate instability
      use pic_logger, only: logger => global_logger
      class(cc_energy_t), intent(in) :: this

      if (this%singles > 0.0_dp) then
         call logger%warning("CC singles correlation energy is positive - possible instability!")
      end if

      if (this%doubles > 0.0_dp) then
         call logger%warning("CC doubles correlation energy is positive - possible instability!")
      end if

      if (this%triples > 0.0_dp) then
         call logger%warning("CC triples correlation energy is positive - possible instability!")
      end if
   end subroutine cc_check_stability

   pure function energy_total(this) result(total)
      !! Compute total energy from all components
      class(energy_t), intent(in) :: this
      real(dp) :: total

      ! this line needs to me modified if more components are added
      total = this%scf + this%mp2%total() + this%cc%total()
   end function energy_total

   subroutine energy_reset(this)
      !! Reset all energy components to zero
      class(energy_t), intent(inout) :: this
      this%scf = 0.0_dp
      call this%mp2%reset()
      call this%cc%reset()
   end subroutine energy_reset

   subroutine result_destroy(this)
      !! Clean up allocated memory in calculation_result_t
      class(calculation_result_t), intent(inout) :: this
      if (allocated(this%gradient)) deallocate (this%gradient)
      if (allocated(this%sigma)) deallocate (this%sigma)
      if (allocated(this%hessian)) deallocate (this%hessian)
      if (allocated(this%dipole)) deallocate (this%dipole)
      if (allocated(this%dipole_derivatives)) deallocate (this%dipole_derivatives)
      call this%reset()
   end subroutine result_destroy

   subroutine result_reset(this)
      !! Reset all values and flags in calculation_result_t
      class(calculation_result_t), intent(inout) :: this
      call this%energy%reset()
      call this%error%clear()
      this%has_energy = .false.
      this%has_gradient = .false.
      this%has_sigma = .false.
      this%has_hessian = .false.
      this%has_dipole = .false.
      this%has_dipole_derivatives = .false.
      this%has_error = .false.
      this%scf_status = SCF_UNKNOWN
      this%scf_iterations = 0
      this%homo = 0.0_dp
      this%lumo = 0.0_dp
      this%has_orbitals = .false.
   end subroutine result_reset

   !---------------------------------------------------------------------------
   ! mbe_result_t type-bound procedures
   !---------------------------------------------------------------------------

   subroutine mbe_result_destroy(this)
      !! Clean up allocated memory in mbe_result_t
      class(mbe_result_t), intent(inout) :: this
      if (allocated(this%gradient)) deallocate (this%gradient)
      if (allocated(this%hessian)) deallocate (this%hessian)
      if (allocated(this%dipole)) deallocate (this%dipole)
      if (allocated(this%dipole_derivatives)) deallocate (this%dipole_derivatives)
      call this%reset()
   end subroutine mbe_result_destroy

   subroutine mbe_result_reset(this)
      !! Reset all values and flags in mbe_result_t
      class(mbe_result_t), intent(inout) :: this
      this%total_energy = 0.0_dp
      this%has_energy = .false.
      this%has_gradient = .false.
      this%has_hessian = .false.
      this%has_dipole = .false.
      this%has_dipole_derivatives = .false.
   end subroutine mbe_result_reset

   subroutine mbe_result_allocate_gradient(this, total_atoms)
      !! Allocate gradient array for total_atoms
      class(mbe_result_t), intent(inout) :: this
      integer, intent(in) :: total_atoms
      if (allocated(this%gradient)) deallocate (this%gradient)
      allocate (this%gradient(3, total_atoms))
      this%gradient = 0.0_dp
   end subroutine mbe_result_allocate_gradient

   subroutine mbe_result_allocate_hessian(this, total_atoms)
      !! Allocate hessian array for total_atoms
      class(mbe_result_t), intent(inout) :: this
      integer, intent(in) :: total_atoms
      integer :: hess_dim
      hess_dim = 3*total_atoms
      if (allocated(this%hessian)) deallocate (this%hessian)
      allocate (this%hessian(hess_dim, hess_dim))
      this%hessian = 0.0_dp
   end subroutine mbe_result_allocate_hessian

   subroutine mbe_result_allocate_dipole(this)
      !! Allocate dipole array (always 3 components)
      class(mbe_result_t), intent(inout) :: this
      if (allocated(this%dipole)) deallocate (this%dipole)
      allocate (this%dipole(3))
      this%dipole = 0.0_dp
   end subroutine mbe_result_allocate_dipole

   subroutine send_error_state(result, comm, dest, tag)
      !! Send the failure flag, code and trace that travel with a result
      !!
      !! This is not optional bookkeeping. Without it a fragment that failed on
      !! a worker arrives at the coordinator looking healthy: `has_error` stays
      !! at its default `.false.`, the untouched SCF energy of 0 is taken at
      !! face value, and the many-body expansion quietly sums zeros. Every
      !! `has_error` check on the receiving side depends on this being sent.
      type(calculation_result_t), intent(in) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: dest, tag

      character(len=:), allocatable :: trace
      integer(int32), allocatable :: encoded(:)
      integer :: i, n

      call send(comm, result%has_error, dest, tag)
      if (.not. result%has_error) return

      call send(comm, result%error%get_code(), dest, tag)

      ! pic-mpi carries no character type, so the trace goes as one character
      ! code per element. It is never empty: a zero-length message would leave
      ! the coordinator with a failure it cannot name.
      trace = result%error%get_full_trace()
      if (len_trim(trace) == 0) trace = "Calculation failed without a diagnostic"
      n = len_trim(trace)
      allocate (encoded(n))
      do i = 1, n
         encoded(i) = iachar(trace(i:i))
      end do
      call send(comm, encoded, dest, tag)
   end subroutine send_error_state

   subroutine recv_error_state(result, comm, source, tag)
      !! Receive the failure flag, code and trace that travel with a result
      type(calculation_result_t), intent(inout) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: source, tag

      type(MPI_Status) :: status
      integer(int32), allocatable :: encoded(:)
      character(len=:), allocatable :: trace
      integer(int32) :: code
      integer :: i

      call recv(comm, result%has_error, source, tag, status)
      if (.not. result%has_error) then
         call result%error%clear()
         return
      end if

      call recv(comm, code, source, tag, status)
      call recv(comm, encoded, source, tag, status)

      allocate (character(len=size(encoded)) :: trace)
      do i = 1, size(encoded)
         trace(i:i) = achar(encoded(i))
      end do
      call result%error%set(int(code), trace)

      ! A failed fragment has no usable energy, whatever the flags said.
      result%has_energy = .false.
   end subroutine recv_error_state

   subroutine result_send(result, comm, dest, tag)
      !! Send calculation result over MPI (blocking)
      !! Sends energy components and conditionally sends gradient based on has_gradient flag
      type(calculation_result_t), intent(in) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: dest, tag

      ! Send energy components
      call send(comm, result%energy%scf, dest, tag)
      call send(comm, result%energy%mp2%ss, dest, tag)
      call send(comm, result%energy%mp2%os, dest, tag)
      call send(comm, result%energy%cc%singles, dest, tag)
      call send(comm, result%energy%cc%doubles, dest, tag)
      call send(comm, result%energy%cc%triples, dest, tag)

      ! Send fragment metadata
      call send(comm, result%distance, dest, tag)
      ! Convergence travels with the energy or it is not knowable at
      ! the coordinator, which is the only rank that writes anything.
      ! A field added to the type and not to these four routines
      ! leaves every parallel run reporting the default.
      call send(comm, result%scf_status, dest, tag)
      call send(comm, result%scf_iterations, dest, tag)
      call send(comm, result%has_orbitals, dest, tag)
      call send(comm, result%homo, dest, tag)
      call send(comm, result%lumo, dest, tag)

      ! Send gradient flag and data if present
      call send(comm, result%has_gradient, dest, tag)
      if (result%has_gradient) then
         call send(comm, result%gradient, dest, tag)
      end if

      ! Send Hessian flag and data if present
      call send(comm, result%has_hessian, dest, tag)
      if (result%has_hessian) then
         call send(comm, result%hessian, dest, tag)
      end if

      ! Send dipole flag and data if present
      call send(comm, result%has_dipole, dest, tag)
      if (result%has_dipole) then
         call send(comm, result%dipole, dest, tag)
      end if

      ! Send dipole derivatives flag and data if present
      call send(comm, result%has_dipole_derivatives, dest, tag)
      if (result%has_dipole_derivatives) then
         call send(comm, result%dipole_derivatives, dest, tag)
      end if

      ! Failure state last, so the receiver has the whole payload drained
      ! before it decides whether to trust any of it.
      call send_error_state(result, comm, dest, tag)
   end subroutine result_send

   subroutine result_isend(result, comm, dest, tag, req)
      !! Send calculation result over MPI (non-blocking)
      !! Sends SCF energy (non-blocking) and other components (blocking)
      type(calculation_result_t), intent(in) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: dest, tag
      type(request_t), intent(out) :: req

      ! Send SCF energy (non-blocking)
      call isend(comm, result%energy%scf, dest, tag, req)

      ! Send other energy components (blocking to avoid needing multiple request handles)
      call send(comm, result%energy%mp2%ss, dest, tag)
      call send(comm, result%energy%mp2%os, dest, tag)
      call send(comm, result%energy%cc%singles, dest, tag)
      call send(comm, result%energy%cc%doubles, dest, tag)
      call send(comm, result%energy%cc%triples, dest, tag)

      ! Send fragment metadata
      call send(comm, result%distance, dest, tag)
      ! Convergence travels with the energy or it is not knowable at
      ! the coordinator, which is the only rank that writes anything.
      ! A field added to the type and not to these four routines
      ! leaves every parallel run reporting the default.
      call send(comm, result%scf_status, dest, tag)
      call send(comm, result%scf_iterations, dest, tag)
      call send(comm, result%has_orbitals, dest, tag)
      call send(comm, result%homo, dest, tag)
      call send(comm, result%lumo, dest, tag)

      ! Send gradient flag and data (blocking to avoid needing multiple request handles)
      call send(comm, result%has_gradient, dest, tag)
      if (result%has_gradient) then
         call send(comm, result%gradient, dest, tag)
      end if

      ! Send Hessian flag and data (blocking to avoid needing multiple request handles)
      call send(comm, result%has_hessian, dest, tag)
      if (result%has_hessian) then
         call send(comm, result%hessian, dest, tag)
      end if

      ! Send dipole flag and data (blocking to avoid needing multiple request handles)
      call send(comm, result%has_dipole, dest, tag)
      if (result%has_dipole) then
         call send(comm, result%dipole, dest, tag)
      end if

      ! Send dipole derivatives flag and data (blocking to avoid needing multiple request handles)
      call send(comm, result%has_dipole_derivatives, dest, tag)
      if (result%has_dipole_derivatives) then
         call send(comm, result%dipole_derivatives, dest, tag)
      end if

      ! Failure state last, so the receiver has the whole payload drained
      ! before it decides whether to trust any of it.
      call send_error_state(result, comm, dest, tag)
   end subroutine result_isend

   subroutine result_recv(result, comm, source, tag, status)
      !! Receive calculation result over MPI (blocking)
      !! Receives energy components and conditionally receives gradient based on flag
      type(calculation_result_t), intent(inout) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: source, tag
      type(MPI_Status), intent(out) :: status

      ! Same reuse hazard as result_irecv -- see the note there.
      call result%destroy()

      ! Receive energy components
      call recv(comm, result%energy%scf, source, tag, status)
      call recv(comm, result%energy%mp2%ss, source, tag, status)
      call recv(comm, result%energy%mp2%os, source, tag, status)
      call recv(comm, result%energy%cc%singles, source, tag, status)
      call recv(comm, result%energy%cc%doubles, source, tag, status)
      call recv(comm, result%energy%cc%triples, source, tag, status)
      result%has_energy = .true.

      ! Receive fragment metadata
      call recv(comm, result%distance, source, tag, status)
      call recv(comm, result%scf_status, source, tag, status)
      call recv(comm, result%scf_iterations, source, tag, status)
      call recv(comm, result%has_orbitals, source, tag, status)
      call recv(comm, result%homo, source, tag, status)
      call recv(comm, result%lumo, source, tag, status)

      ! Receive gradient flag and data if present
      call recv(comm, result%has_gradient, source, tag, status)
      if (result%has_gradient) then
         ! Receive allocatable gradient array (MPI lib handles allocation)
         call recv(comm, result%gradient, source, tag, status)
      end if

      ! Receive Hessian flag and data if present
      call recv(comm, result%has_hessian, source, tag, status)
      if (result%has_hessian) then
         ! Receive allocatable Hessian array (MPI lib handles allocation)
         call recv(comm, result%hessian, source, tag, status)
      end if

      ! Receive dipole flag and data if present
      call recv(comm, result%has_dipole, source, tag, status)
      if (result%has_dipole) then
         ! Receive allocatable dipole array (MPI lib handles allocation)
         call recv(comm, result%dipole, source, tag, status)
      end if

      ! Receive dipole derivatives flag and data if present
      call recv(comm, result%has_dipole_derivatives, source, tag, status)
      if (result%has_dipole_derivatives) then
         ! Receive allocatable dipole derivatives array (MPI lib handles allocation)
         call recv(comm, result%dipole_derivatives, source, tag, status)
      end if

      call recv_error_state(result, comm, source, tag)
   end subroutine result_recv

   subroutine result_irecv(result, comm, source, tag, req)
      !! Receive calculation result over MPI (non-blocking)
      !! Receives SCF energy (non-blocking) and other components (blocking)
      type(calculation_result_t), intent(inout) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: source, tag
      type(request_t), intent(out) :: req
      type(MPI_Status) :: status

      ! Start from an empty result. pic-mpi's array receives allocate only when
      ! the target is unallocated -- they will not resize a buffer that is
      ! already there, yet still receive the sender's full element count into
      ! it. A caller reusing one result across messages of different shapes
      ! (a relay loop forwarding fragments of unequal size, say) would overflow
      ! the smaller allocation and corrupt the heap, surfacing much later as an
      ! abort inside an unrelated deallocate. Clearing here makes reuse safe for
      ! every caller rather than each having to remember.
      call result%destroy()

      ! Receive SCF energy (non-blocking)
      call irecv(comm, result%energy%scf, source, tag, req)

      ! Receive other energy components (blocking to avoid needing multiple request handles)
      call recv(comm, result%energy%mp2%ss, source, tag, status)
      call recv(comm, result%energy%mp2%os, source, tag, status)
      call recv(comm, result%energy%cc%singles, source, tag, status)
      call recv(comm, result%energy%cc%doubles, source, tag, status)
      call recv(comm, result%energy%cc%triples, source, tag, status)
      result%has_energy = .true.

      ! Receive fragment metadata
      call recv(comm, result%distance, source, tag, status)
      call recv(comm, result%scf_status, source, tag, status)
      call recv(comm, result%scf_iterations, source, tag, status)
      call recv(comm, result%has_orbitals, source, tag, status)
      call recv(comm, result%homo, source, tag, status)
      call recv(comm, result%lumo, source, tag, status)

      ! Receive gradient flag and data (blocking to avoid needing multiple request handles)
      call recv(comm, result%has_gradient, source, tag, status)
      if (result%has_gradient) then
         ! Receive allocatable gradient array (MPI lib handles allocation)
         call recv(comm, result%gradient, source, tag, status)
      end if

      ! Receive Hessian flag and data (blocking to avoid needing multiple request handles)
      call recv(comm, result%has_hessian, source, tag, status)
      if (result%has_hessian) then
         ! Receive allocatable Hessian array (MPI lib handles allocation)
         call recv(comm, result%hessian, source, tag, status)
      end if

      ! Receive dipole flag and data (blocking to avoid needing multiple request handles)
      call recv(comm, result%has_dipole, source, tag, status)
      if (result%has_dipole) then
         ! Receive allocatable dipole array (MPI lib handles allocation)
         call recv(comm, result%dipole, source, tag, status)
      end if

      ! Receive dipole derivatives flag and data (blocking to avoid needing multiple request handles)
      call recv(comm, result%has_dipole_derivatives, source, tag, status)
      if (result%has_dipole_derivatives) then
         ! Receive allocatable dipole derivatives array (MPI lib handles allocation)
         call recv(comm, result%dipole_derivatives, source, tag, status)
      end if

      call recv_error_state(result, comm, source, tag)
   end subroutine result_irecv

   pure function scf_status_label(status) result(text)
      !! One word for a status, for a table column or a log line
      integer, intent(in) :: status
      character(len=3) :: text

      select case (status)
      case (SCF_CONVERGED)
         text = "yes"
      case (SCF_NOT_CONVERGED)
         text = "NO "
      case default
         text = "?  "
      end select
   end function scf_status_label

   pure subroutine frontier_orbitals(energies, occupations, homo, lumo, found)
      !! Pick HOMO and LUMO out of an ascending orbital spectrum
      !!
      !! Occupation-driven rather than counting electrons, so the same routine
      !! serves restricted and unrestricted spectra. The threshold is half an
      !! electron: with Fermi smearing the frontier occupations are fractional
      !! rather than 0 or 2, and any cut inside that range names the same pair
      !! for a molecule with a real gap. It is arbitrary for a system without
      !! one, which is a system whose gap was not going to mean much anyway.
      real(dp), intent(in) :: energies(:)
      real(dp), intent(in) :: occupations(:)
      real(dp), intent(out) :: homo, lumo
      logical, intent(out) :: found

      integer :: i, i_homo

      homo = 0.0_dp
      lumo = 0.0_dp
      found = .false.
      i_homo = 0
      do i = 1, min(size(energies), size(occupations))
         if (occupations(i) > 0.5_dp) i_homo = i
      end do

      ! Nothing occupied, or nothing empty: no frontier pair to report, and
      ! saying so beats inventing one from an end of the spectrum.
      if (i_homo < 1 .or. i_homo >= size(energies)) return

      homo = energies(i_homo)
      lumo = energies(i_homo + 1)
      found = .true.
   end subroutine frontier_orbitals

end module mqc_result_types
