!! Constraining a Fock matrix to a block structure the caller chooses
module mqc_fock_projector
   !! A fixed orbital basis, a partition of it into frozen and unfrozen blocks,
   !! and one operation: force the Fock matrix block diagonal in that partition.
   !!
   !! This is GAFO -- the generalized adjusted frozen orbital scheme -- reduced
   !! to the part that has nothing to do with why the orbitals were frozen. What
   !! makes an orbital frozen, and where the frozen orbitals came from, is the
   !! caller's business; what is here is the constraint itself, which is dense
   !! linear algebra over `F`, `C` and `S` and knows nothing about fragments,
   !! bonds or integrals.
   !!
   !! **Why zeroing rather than a level shift.** The obvious way to keep an
   !! orbital out of the occupied space is to add a penalty, `F + B S v v^T S`
   !! with `B` large. That raises the orbital's energy and does not decouple it:
   !! the off-diagonal Fock elements between the frozen and variational
   !! subspaces survive untouched, the two mix, and the matrix is not block
   !! diagonal in the partition the method assumes it is. Raising `B` makes the
   !! diagonal larger and leaves the coupling exactly where it was. Zeroing the
   !! couplings is the fix, and it is why this scheme exists at all.
   !!
   !! **The transform.** With `C` an orthonormal basis in the metric `S` --
   !! `C^T S C = I`, frozen orbitals in the leading columns --
   !!
   !!     F_mo = C^T F C          to the frozen-orbital basis
   !!     F_mo -> block diagonal  the constraint
   !!     F    = (SC) F_mo (SC)^T back again
   !!
   !! The back transform is the inverse of the forward one because `C^T S C = I`
   !! makes `C^-1 = C^T S`, so `(C^T)^-1 = S C`. Not an approximation and not a
   !! symmetrisation: with nothing frozen and `C` spanning the whole AO space
   !! the round trip returns `F` to roundoff, which is what the test asserts.
   !!
   !! `S C` is formed once in `init` rather than each application, since neither
   !! factor moves once the frozen orbitals are chosen.
   !!
   !! **Choosing the shift.** Smaller than instinct suggests. The back
   !! transform `(SC) F_mo (SC)^T` spreads the shifted block over every element,
   !! so the *unfrozen* block -- the one carrying the physics -- is clean only
   !! to about `shift * epsilon`. At the 1e6 the reference implementations use,
   !! that is around 1e-10 Hartree, against a validation suite that works at
   !! 1e-9.
   !!
   !! Nothing here needs a shift that large. A penalty scheme does, because it
   !! is trying to overwhelm a coupling it never removed; here the coupling is
   !! gone, and the shift has only to lift the frozen virtuals clear of the
   !! occupied manifold so that filling by aufbau does not reach them. A few
   !! hundred Hartree is decisive for that and costs three or four digits less.
   !! `projector_shift_sets_the_precision_floor` measures it.
   !!
   !! **Where to apply it.** Immediately after the Fock build and *before* the
   !! DIIS error and extrapolation. The constraint is a linear map and DIIS
   !! extrapolation is a linear combination whose coefficients sum to one, so a
   !! combination of matrices sharing this block structure has it too -- the
   !! zeros stay zero and the shifted diagonal comes through unchanged. Applied
   !! there, every consumer downstream sees the same constrained object; applied
   !! after DIIS, the error vector measures a matrix nobody diagonalizes.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: fock_projector_t

   !> A frozen-orbital partition, and the constraint it implies
   type :: fock_projector_t
      real(dp), allocatable :: basis(:, :)
         !! `C`, `n_ao` by `n_mo`, orthonormal in the metric `S`, with the
         !! frozen orbitals in the leading columns. Their being leading is what
         !! lets a block be named by an index range.
      real(dp), allocatable :: sc(:, :)
         !! `S C`, the back transform, formed once
      integer :: n_frozen_occ = 0
         !! Columns `1 : n_frozen_occ` -- frozen and occupied. Their block is
         !! left alone; only its couplings to everything else are cut.
      integer :: n_frozen = 0
         !! Columns `n_frozen_occ + 1 : n_frozen` are frozen and virtual: made
         !! degenerate at `shift` so nothing occupies them. Zero here means the
         !! projector is inactive and `apply` returns at once.
      real(dp) :: shift = 0.0_dp
         !! The energy the frozen virtual block is held at
   contains
      procedure :: init
      procedure :: apply
   end type fock_projector_t

contains

   subroutine init(this, basis, overlap, n_frozen_occ, n_frozen, shift, error)
      !! Take the basis and partition, and form the back transform
      class(fock_projector_t), intent(out) :: this
      real(dp), intent(in) :: basis(:, :)
      real(dp), intent(in) :: overlap(:, :)
      integer, intent(in) :: n_frozen_occ, n_frozen
      real(dp), intent(in) :: shift
      type(error_t), intent(inout) :: error

      integer :: n_ao, n_mo

      if (error%has_error()) return
      n_ao = size(basis, 1)
      n_mo = size(basis, 2)

      if (size(overlap, 1) /= n_ao .or. size(overlap, 2) /= n_ao) then
         call error%set(ERROR_VALIDATION, "fock projector: the overlap is not "// &
                        "square in the basis the orbitals are expressed in")
         return
      end if
      if (n_frozen_occ < 0 .or. n_frozen < n_frozen_occ .or. n_frozen > n_mo) then
         call error%set(ERROR_VALIDATION, "fock projector: the frozen blocks must "// &
                        "nest as 0 <= n_frozen_occ <= n_frozen <= n_mo")
         return
      end if

      this%n_frozen_occ = n_frozen_occ
      this%n_frozen = n_frozen
      this%shift = shift
      this%basis = basis

      allocate (this%sc(n_ao, n_mo))
      call pic_gemm(overlap, basis, this%sc)
   end subroutine init

   subroutine apply(this, fock, error)
      !! Force `fock` block diagonal in the frozen partition
      class(fock_projector_t), intent(in) :: this
      real(dp), intent(inout) :: fock(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: work(:, :), f_mo(:, :)
      integer :: n_ao, n_mo, nfo, nfr, i

      if (error%has_error()) return

      ! Nothing frozen is not "the identity to roundoff", it is the identity:
      ! returning here rather than transforming there and back is what makes a
      ! calculation with no frozen orbitals bit-identical to one that never had
      ! a projector, which is the regression guarantee this is worth having.
      if (this%n_frozen == 0) return

      n_ao = size(fock, 1)
      n_mo = size(this%basis, 2)
      if (size(this%basis, 1) /= n_ao .or. size(fock, 2) /= n_ao) then
         call error%set(ERROR_VALIDATION, "fock projector: the Fock matrix is not "// &
                        "square in the basis the frozen orbitals were built in")
         return
      end if

      allocate (work(n_ao, n_mo), f_mo(n_mo, n_mo))

      call pic_gemm(fock, this%basis, work)              ! F C
      call pic_gemm(this%basis, work, f_mo, transa="T")  ! C^T F C

      nfo = this%n_frozen_occ
      nfr = this%n_frozen

      ! Every coupling out of the frozen-occupied block. The block itself is
      ! left as it is -- those orbitals stay occupied and keep their own
      ! structure; what must go is their mixing with anything else.
      if (nfo > 0) then
         f_mo(nfo + 1:, :nfo) = 0.0_dp
         f_mo(:nfo, nfo + 1:) = 0.0_dp
      end if

      ! The frozen-virtual block, held degenerate at `shift` and cut from the
      ! unfrozen space. Degenerate rather than merely shifted: an off-diagonal
      ! left inside this block would let its orbitals rotate among themselves,
      ! and they are meant to be the orbitals that were chosen, not a basis for
      ! them.
      if (nfr > nfo) then
         f_mo(nfo + 1:nfr, nfo + 1:nfr) = 0.0_dp
         do i = nfo + 1, nfr
            f_mo(i, i) = this%shift
         end do
         f_mo(nfr + 1:, nfo + 1:nfr) = 0.0_dp
         f_mo(nfo + 1:nfr, nfr + 1:) = 0.0_dp
      end if

      call pic_gemm(this%sc, f_mo, work)                 ! (SC) F_mo
      call pic_gemm(work, this%sc, fock, transb="T")     ! (SC) F_mo (SC)^T

      deallocate (work, f_mo)
   end subroutine apply

end module mqc_fock_projector
