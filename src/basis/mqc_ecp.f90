!! Effective core potential data structures
module mqc_ecp
   !! Holds an effective core potential in the semi-local form every code and
   !! every tabulation uses:
   !!
   !!     U(r) = U_L(r) + sum_{l=0}^{L-1} [ U_l(r) - U_L(r) ] P_l
   !!
   !! `U_L` is the **local** channel, the one with the highest angular
   !! momentum; the lower `U_l` act only on their own angular momentum through
   !! the projector `P_l`. Each channel is a sum of primitives
   !!
   !!     U_l(r) = sum_k d_k r^(n_k - 2) exp(-z_k r^2)
   !!
   !! Radial powers are stored exactly as the source file gives them; the `-2`
   !! is a convention difference between tabulations, and the consumer has to
   !! know which one its library wants.
   !!
   !! `core_electrons` is how many electrons the potential stands in for.
   !! Everything downstream has to agree on it: the atom contributes
   !! `Z - core_electrons` electrons to the SCF, screens the nucleus to the
   !! same effective charge, and repels other nuclei with it.
   use pic_types, only: dp
   implicit none
   private

   public :: ecp_shell_type      !! One angular-momentum channel
   public :: atomic_ecp_type     !! The whole potential for one atom
   public :: molecular_ecp_type  !! One entry per atom in a molecule

   type :: ecp_shell_type
      !! One channel of a semi-local ECP
      integer :: ang_mom = -1
         !! The l this channel projects onto; the local channel carries the
         !! highest l in the potential and is not projected at all
      integer :: nprim = 0
      integer, allocatable :: radial_powers(:)  !! r exponents, as tabulated
      real(dp), allocatable :: exponents(:)     !! Gaussian exponents
      real(dp), allocatable :: coefficients(:)
   contains
      procedure :: allocate_arrays => ecp_shell_allocate
      procedure :: destroy => ecp_shell_destroy
   end type ecp_shell_type

   type :: atomic_ecp_type
      !! The effective core potential for one atom
      character(len=:), allocatable :: element
      integer :: core_electrons = 0
         !! Electrons this potential replaces. Zero, with `has_ecp` false, for
         !! a light atom carrying no ECP.
      logical :: has_ecp = .false.
      type(ecp_shell_type) :: local
         !! U_L, the highest angular momentum channel
      type(ecp_shell_type), allocatable :: projected(:)
         !! U_0 .. U_(L-1), each acting through its own projector
      integer :: n_projected = 0
   contains
      procedure :: allocate_projected => atomic_ecp_allocate
      procedure :: destroy => atomic_ecp_destroy
   end type atomic_ecp_type

   type :: molecular_ecp_type
      !! One entry per atom, parallel to `molecular_basis_type%elements`
      !!
      !! Indexed by atom, not by element, so an entry lines up with the atom
      !! index everything else uses. Atoms without an ECP get an entry with
      !! `has_ecp` false rather than being skipped.
      type(atomic_ecp_type), allocatable :: atoms(:)
      integer :: natoms = 0
      integer :: n_ecp_atoms = 0
         !! How many of them actually carry a potential
      ! TODO(mqc): maintained only by `build_molecular_ecp_json`, and `any_ecp`
      ! reads nothing else. A `molecular_ecp_type` filled any other way reports
      ! no ECPs from `any_ecp` while `core_electrons`, which scans `atoms`,
      ! reports them.
   contains
      procedure :: allocate_atoms => molecular_ecp_allocate
      procedure :: destroy => molecular_ecp_destroy
      procedure :: core_electrons => molecular_ecp_core_electrons
      procedure :: any_ecp => molecular_ecp_any
   end type molecular_ecp_type

contains

   subroutine ecp_shell_allocate(this, nprim)
      !! Size one channel for `nprim` primitives
      class(ecp_shell_type), intent(inout) :: this
      integer, intent(in) :: nprim

      call this%destroy()
      this%nprim = nprim
      allocate (this%radial_powers(nprim))
      allocate (this%exponents(nprim))
      allocate (this%coefficients(nprim))
      this%radial_powers = 0
      this%exponents = 0.0_dp
      this%coefficients = 0.0_dp
   end subroutine ecp_shell_allocate

   subroutine ecp_shell_destroy(this)
      class(ecp_shell_type), intent(inout) :: this

      if (allocated(this%radial_powers)) deallocate (this%radial_powers)
      if (allocated(this%exponents)) deallocate (this%exponents)
      if (allocated(this%coefficients)) deallocate (this%coefficients)
      this%nprim = 0
      this%ang_mom = -1
   end subroutine ecp_shell_destroy

   subroutine atomic_ecp_allocate(this, n_projected)
      !! Size the projected channels
      class(atomic_ecp_type), intent(inout) :: this
      integer, intent(in) :: n_projected

      if (allocated(this%projected)) deallocate (this%projected)
      this%n_projected = n_projected
      if (n_projected > 0) allocate (this%projected(n_projected))
   end subroutine atomic_ecp_allocate

   subroutine atomic_ecp_destroy(this)
      class(atomic_ecp_type), intent(inout) :: this
      integer :: i

      if (allocated(this%element)) deallocate (this%element)
      call this%local%destroy()
      if (allocated(this%projected)) then
         do i = 1, size(this%projected)
            call this%projected(i)%destroy()
         end do
         deallocate (this%projected)
      end if
      this%n_projected = 0
      this%core_electrons = 0
      this%has_ecp = .false.
   end subroutine atomic_ecp_destroy

   subroutine molecular_ecp_allocate(this, natoms)
      !! One entry per atom, all of them empty until filled
      class(molecular_ecp_type), intent(inout) :: this
      integer, intent(in) :: natoms

      call this%destroy()
      this%natoms = natoms
      if (natoms > 0) allocate (this%atoms(natoms))
   end subroutine molecular_ecp_allocate

   subroutine molecular_ecp_destroy(this)
      class(molecular_ecp_type), intent(inout) :: this
      integer :: i

      if (allocated(this%atoms)) then
         do i = 1, size(this%atoms)
            call this%atoms(i)%destroy()
         end do
         deallocate (this%atoms)
      end if
      this%natoms = 0
      this%n_ecp_atoms = 0
   end subroutine molecular_ecp_destroy

   pure function molecular_ecp_core_electrons(this) result(total)
      !! Electrons replaced across the whole molecule
      !!
      !! What the SCF's electron count has to drop by. A sum over atoms, not
      !! over elements.
      class(molecular_ecp_type), intent(in) :: this
      integer :: total
      integer :: iatom

      total = 0
      if (.not. allocated(this%atoms)) return
      do iatom = 1, this%natoms
         if (this%atoms(iatom)%has_ecp) total = total + this%atoms(iatom)%core_electrons
      end do
   end function molecular_ecp_core_electrons

   pure function molecular_ecp_any(this) result(any_present)
      !! Whether any atom carries a potential at all
      class(molecular_ecp_type), intent(in) :: this
      logical :: any_present

      any_present = (this%n_ecp_atoms > 0)
   end function molecular_ecp_any

end module mqc_ecp
