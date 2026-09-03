!! Constraining a Fock matrix to a block structure the caller chooses
module mqc_fock_projector
   !! A fixed orbital basis, a partition of it into frozen and unfrozen blocks,
   !! and one operation: force the Fock matrix block diagonal in that partition.
   !!
   !! This is GAFO -- the generalized adjusted frozen orbital scheme -- reduced
   !! to dense linear algebra over `F`, `C` and `S`. What makes an orbital
   !! frozen, and where the frozen orbitals came from, is the caller's business.
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
   !! the round trip returns `F` to roundoff.
   !!
   !! The couplings are zeroed and not penalised. `F + B S v v^T S` raises the
   !! frozen orbital's energy and leaves the off-diagonal Fock elements between
   !! the frozen and variational subspaces exactly where they were, however
   !! large `B` is.
   !!
   !! **Choosing the shift.** Smaller than instinct suggests. The back
   !! transform `(SC) F_mo (SC)^T` spreads the shifted block over every element,
   !! so the *unfrozen* block is clean only to about `shift * epsilon`. At the
   !! 1e6 the reference implementations use, that is around 1e-10 Hartree,
   !! against a validation suite that works at 1e-9. The shift has only to lift
   !! the frozen virtuals clear of the occupied manifold, for which a few
   !! hundred Hartree is decisive.
   !!
   !! **Where to apply it.** Immediately after the Fock build and *before* the
   !! DIIS error and extrapolation. The constraint is a linear map and DIIS
   !! extrapolation is a linear combination whose coefficients sum to one, so a
   !! combination of matrices sharing this block structure has it too. Applied
   !! after DIIS, the error vector measures a matrix nobody diagonalizes.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use mqc_scf_common, only: build_orthogonalizer, LINEAR_DEPENDENCE_TOL
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: fock_projector_t
   public :: build_frozen_basis

   type :: fock_projector_t
      !! A frozen-orbital partition, and the constraint it implies
      real(dp), allocatable :: basis(:, :)
         !! `C`, `n_ao` by `n_mo`, orthonormal in the metric `S`, with the
         !! frozen orbitals in the leading columns. Their being leading is what
         !! lets a block be named by an index range.
      real(dp), allocatable :: sc(:, :)
         !! `S C`, the back transform, formed once in `init`
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

   subroutine build_frozen_basis(frozen, n_frozen_occ, overlap, basis, n_mo, error)
      !! Turn a set of chosen orbitals into a basis `init` will accept
      !!
      !! The constraint names its blocks by index range, so the frozen orbitals
      !! have to sit in the leading columns of an orthonormal basis, in block
      !! order, with the rest of the space after them. What arrives from a
      !! localization is neither orthogonal nor normalised.
      !!
      !! **Each block is orthonormalised separately.** Doing the frozen set in
      !! one pass mixes the frozen occupied orbitals with the frozen virtual
      !! ones, and then the block boundary means nothing: `apply` would leave
      !! part of a virtual untouched and hold part of an occupied at the level
      !! shift. Mixing *within* a block is harmless, since both blocks are
      !! invariant under a rotation among their own members.
      !!
      !! Canonical rather than symmetric orthonormalisation: the completion
      !! step is deliberately rank deficient, so null directions have to be
      !! dropped rather than amplified.
      !!
      !! Near-dependence among the frozen orbitals themselves is refused --
      !! two localized orbitals that are the same orbital means the bond was
      !! assigned twice, and freezing it twice would be silent.
      real(dp), intent(in) :: frozen(:, :)
         !! `n_ao` by `n_frozen`, occupied first, in the order the blocks want
      integer, intent(in) :: n_frozen_occ
      real(dp), intent(in) :: overlap(:, :)
      real(dp), allocatable, intent(out) :: basis(:, :)
      integer, intent(out) :: n_mo
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: x(:, :), piece(:, :)
      integer :: n_ao, n_frozen, n_virt

      if (error%has_error()) return
      n_ao = size(frozen, 1)
      n_frozen = size(frozen, 2)
      n_virt = n_frozen - n_frozen_occ

      if (n_frozen_occ < 0 .or. n_virt < 0) then
         call error%set(ERROR_VALIDATION, "frozen basis: the occupied block cannot be "// &
                        "negative or larger than the frozen set")
         return
      end if
      if (size(overlap, 1) /= n_ao .or. size(overlap, 2) /= n_ao) then
         call error%set(ERROR_VALIDATION, "frozen basis: the overlap is not square in "// &
                        "the basis the orbitals are expressed in")
         return
      end if

      call build_orthogonalizer(overlap, x, n_mo, error)
      if (error%has_error()) return
      if (n_frozen > n_mo) then
         call error%set(ERROR_VALIDATION, "frozen basis: more orbitals were frozen than "// &
                        "the basis can hold once near-dependent combinations are dropped")
         return
      end if

      allocate (basis(n_ao, n_mo), source=0.0_dp)

      if (n_frozen_occ > 0) then
         piece = frozen(:, :n_frozen_occ)
         call orthonormalize_against(overlap, basis, 0, piece, error)
         if (error%has_error()) return
         if (size(piece, 2) /= n_frozen_occ) then
            call error%set(ERROR_VALIDATION, "frozen basis: the frozen occupied orbitals "// &
                           "are linearly dependent, so one of them is already frozen")
            return
         end if
         basis(:, :n_frozen_occ) = piece
      end if

      if (n_virt > 0) then
         piece = frozen(:, n_frozen_occ + 1:)
         call orthonormalize_against(overlap, basis, n_frozen_occ, piece, error)
         if (error%has_error()) return
         if (size(piece, 2) /= n_virt) then
            call error%set(ERROR_VALIDATION, "frozen basis: a frozen virtual orbital is "// &
                           "already spanned by the frozen occupied ones")
            return
         end if
         basis(:, n_frozen_occ + 1:n_frozen) = piece
      end if

      ! Whatever is left of the space: a basis for all of it, with what is
      ! already placed projected out, leaves exactly `n_mo - n_frozen`
      ! directions.
      piece = x
      call orthonormalize_against(overlap, basis, n_frozen, piece, error)
      if (error%has_error()) return
      if (size(piece, 2) /= n_mo - n_frozen) then
         call error%set(ERROR_VALIDATION, "frozen basis: the complement of the frozen "// &
                        "orbitals did not come back with the dimension it must have")
         return
      end if
      if (n_mo > n_frozen) basis(:, n_frozen + 1:) = piece
   end subroutine build_frozen_basis

   subroutine orthonormalize_against(overlap, prior, n_prior, vecs, error)
      !! Project `vecs` off the first `n_prior` columns of `prior`, then
      !! orthonormalise what is left among itself in the metric `overlap`
      !!
      !! `vecs` comes back with however many independent directions survived,
      !! which may be fewer columns than it went in with. The caller decides
      !! whether losing one was expected.
      real(dp), intent(in) :: overlap(:, :)
      real(dp), intent(in) :: prior(:, :)
      integer, intent(in) :: n_prior
      real(dp), allocatable, intent(inout) :: vecs(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: sv(:, :), gram(:, :), coeff(:, :), kept(:, :), evals(:)
      integer :: n_ao, n_v, n_keep, i, first, info

      if (error%has_error()) return
      n_ao = size(vecs, 1)
      n_v = size(vecs, 2)

      if (n_prior > 0) then
         allocate (sv(n_ao, n_v), coeff(n_prior, n_v))
         call pic_gemm(overlap, vecs, sv)
         call pic_gemm(prior(:, :n_prior), sv, coeff, transa="T")
         call pic_gemm(prior(:, :n_prior), coeff, vecs, alpha=-1.0_dp, beta=1.0_dp)
         deallocate (sv, coeff)
      end if

      allocate (sv(n_ao, n_v), gram(n_v, n_v), evals(n_v))
      call pic_gemm(overlap, vecs, sv)
      call pic_gemm(vecs, sv, gram, transa="T")
      call pic_syev(gram, evals, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "frozen basis: orthonormalisation failed")
         return
      end if

      ! Ascending, so the directions already spanned lead and are dropped.
      n_keep = count(evals > LINEAR_DEPENDENCE_TOL)
      first = n_v - n_keep + 1

      allocate (kept(n_ao, max(n_keep, 1)))
      if (n_keep > 0) then
         call pic_gemm(vecs, gram(:, first:), kept)
         do i = 1, n_keep
            kept(:, i) = kept(:, i)/sqrt(evals(first + i - 1))
         end do
      end if

      deallocate (vecs)
      allocate (vecs(n_ao, n_keep))
      if (n_keep > 0) vecs = kept(:, :n_keep)
   end subroutine orthonormalize_against

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

      ! Returning here rather than transforming there and back is what keeps a
      ! calculation with no frozen orbitals bit-identical to one that never had
      ! a projector.
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
      ! left as it is; what must go is its mixing with anything else.
      if (nfo > 0) then
         f_mo(nfo + 1:, :nfo) = 0.0_dp
         f_mo(:nfo, nfo + 1:) = 0.0_dp
      end if

      ! The frozen-virtual block, held degenerate at `shift` and cut from the
      ! unfrozen space. Degenerate rather than merely shifted: an off-diagonal
      ! left inside this block would let its orbitals rotate among themselves.
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
