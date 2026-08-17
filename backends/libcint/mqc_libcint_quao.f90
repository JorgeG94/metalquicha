!! Quasi-atomic orbitals: the valence-virtual space
module mqc_libcint_quao
   !! The first step of the Ruedenberg construction that needs integrals.
   !!
   !! The virtual space of a large-basis SCF is mostly polarization and diffuse
   !! functions with no chemical content -- you cannot point at an antibonding
   !! orbital in a cc-pVTZ virtual space, and even the *symmetry* of the lowest
   !! virtual changes with the basis. The valence-virtual orbitals are the part
   !! of it that does have chemical content: the antibonding counterparts of the
   !! occupied valence orbitals, recovered by asking which combinations of
   !! virtuals look most like free-atom minimal-basis orbitals.
   !!
   !! What makes them worth having is that the answer barely depends on the
   !! basis. Paper I reports valence-virtual orbital energies that move less
   !! between 6-31G and aug-cc-pVQZ than the occupied ones do, while the
   !! ordinary virtual spectrum changes beyond recognition.
   !!
   !! Reference: West, Schmidt, Gordon, Ruedenberg, J. Chem. Phys. 139, 234107
   !! (2013), section V.A.1 and Appendix A.1.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cgto, only: molecular_basis_type
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use mqc_aambs, only: aambs_file, aambs_dimensions, aambs_dimensions_t
   use mqc_libcint_integrals, only: libcint_molecule_t, mixed_basis_overlap
   implicit none
   private

   public :: build_aambs_molecule    !! The free-atom minimal basis on this geometry
   public :: mo_aambs_overlap        !! < MO | orthogonalized AAMBS >
   public :: valence_virtual_orbitals
   public :: vvo_result_t

   type :: vvo_result_t
      !! What the valence-virtual extraction leaves behind
      real(dp), allocatable :: orbitals(:, :)
         !! (n_ao, n_vvo), in the AO basis, orthonormal and orthogonal to the
         !! occupied space because they are combinations of virtuals only
      real(dp), allocatable :: singular_values(:)
         !! All of them, descending. The first `n_vvo` are retained; the rest
         !! span the valence-external space.
      integer :: n_vvo = 0
      real(dp) :: smallest_retained = 0.0_dp
      real(dp) :: largest_rejected = 0.0_dp
         !! The gap between these two is the diagnostic. Paper I Table I reports
         !! 0.99999 against 0.105-0.272 across eight molecules; anything else
         !! means the projection is not finding a clean valence space.
   end type vvo_result_t

contains

   subroutine build_aambs_molecule(atomic_numbers, element_symbols, coordinates, mol, error)
      !! The accurate atomic minimal basis, on this molecule's nuclei
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)   !! (3, natm), Bohr
      type(libcint_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: error

      type(molecular_basis_type) :: basis
      character(len=:), allocatable :: path

      if (error%has_error()) return
      call aambs_file(path, error)
      if (error%has_error()) return
      call build_molecular_basis_json(path, element_symbols, basis, error)
      if (error%has_error()) return
      call mol%build(atomic_numbers, coordinates, basis, error)
      call basis%destroy()
   end subroutine build_aambs_molecule

   subroutine inverse_sqrt(matrix, result, error, floor)
      !! `A^(-1/2)` for a symmetric positive definite matrix
      real(dp), intent(in) :: matrix(:, :)
      real(dp), allocatable, intent(out) :: result(:, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: floor

      real(dp), allocatable :: vectors(:, :), values(:), scaled(:, :)
      real(dp) :: small
      integer :: n, i, info

      if (error%has_error()) return
      small = 1.0e-10_dp
      if (present(floor)) small = floor

      n = size(matrix, 1)
      allocate (vectors(n, n), values(n), scaled(n, n), result(n, n))
      vectors = matrix
      call pic_syev(vectors, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "the minimal-basis overlap could not be "// &
                        "diagonalized (info = "//to_char(info)//")")
         return
      end if

      ! A minimal basis on well-separated atoms is far from linearly dependent,
      ! so a small eigenvalue here is a statement about the geometry -- atoms on
      ! top of one another -- rather than a threshold to be tuned.
      if (values(1) < small) then
         call error%set(ERROR_VALIDATION, "the free-atom minimal basis is linearly "// &
                        "dependent on this geometry (smallest overlap eigenvalue "// &
                        to_char(values(1))//"). Two atoms are on top of each other, or "// &
                        "close enough that their free-atom orbitals cannot be told apart.")
         return
      end if

      do i = 1, n
         scaled(:, i) = vectors(:, i)/sqrt(values(i))
      end do
      call pic_gemm(scaled, vectors, result, transb="T")
   end subroutine inverse_sqrt

   subroutine mo_aambs_overlap(orbitals, mixed, aambs_overlap, projection, error)
      !! `< MO_p | A#a >`, the molecular orbitals against the orthogonalized AAMBS
      !!
      !! **Both bases have to be orthonormal**, and that is not a detail of
      !! taste. Paper I's Appendix A.1 proves the singular vectors maximize the
      !! projection of one space into the other only under that condition; with
      !! a non-orthogonal basis the largest singular value picks out the
      !! direction with the largest *coefficient*, which is not the same
      !! statement and not the one the method needs. The molecular orbitals are
      !! orthonormal already; the free-atom orbitals are orthonormal within one
      !! atom but not between atoms, so they are symmetrically orthogonalized
      !! across the molecule first -- Paper I's `|A#a>`.
      !!
      !! Paper I notes an earlier formulation of this step "failed to account
      !! for the necessary orthogonalities". This is that account.
      real(dp), intent(in) :: orbitals(:, :)        !! C, (n_ao, n_mo)
      real(dp), intent(in) :: mixed(:, :)           !! < AO | AAMBS >, (n_ao, n_mbs)
      real(dp), intent(in) :: aambs_overlap(:, :)   !! < AAMBS | AAMBS >, (n_mbs, n_mbs)
      real(dp), allocatable, intent(out) :: projection(:, :)   !! (n_mo, n_mbs)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: half(:, :), work(:, :)
      integer :: n_mo, n_mbs

      if (error%has_error()) return
      n_mo = size(orbitals, 2)
      n_mbs = size(mixed, 2)

      call inverse_sqrt(aambs_overlap, half, error)
      if (error%has_error()) return

      allocate (work(n_mo, n_mbs), projection(n_mo, n_mbs))
      call pic_gemm(orbitals, mixed, work, transa="T")
      call pic_gemm(work, half, projection)
      deallocate (work, half)
   end subroutine mo_aambs_overlap

   subroutine valence_virtual_orbitals(orbitals, projection, dims, result, error)
      !! Recover the chemically meaningful part of the virtual space
      !!
      !! The number kept is fixed by counting, not by a threshold:
      !! `n_vvo = n_mbs - n_occ` exactly. The singular values are reported so a
      !! caller can check there is a gap where the cut falls, but they do not
      !! decide where it falls -- a construction that selected on magnitude
      !! would return a different number of orbitals for the same molecule in a
      !! different basis, which is the property the method exists to avoid.
      !!
      !! Solved as the eigenproblem of `Sigma Sigma^T` rather than by an SVD.
      !! Appendix A.1 shows these have the same left singular vectors, and the
      !! external block of the rectangular problem is full of degenerate zero
      !! singular values whose vectors are arbitrary -- GAMESS takes the
      !! eigenvalue route for exactly that reason, and this follows it.
      real(dp), intent(in) :: orbitals(:, :)     !! C, (n_ao, n_mo)
      real(dp), intent(in) :: projection(:, :)   !! < MO | A#a >, (n_mo, n_mbs)
      type(aambs_dimensions_t), intent(in) :: dims
      type(vvo_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: sigma(:, :), b(:, :), values(:), rotated(:, :)
      real(dp), allocatable :: c_virt(:, :)
      integer :: n_ao, n_mo, n_virt, n_occ, i, info

      if (error%has_error()) return
      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_occ = dims%n_occupied
      n_virt = n_mo - n_occ

      if (dims%n_vvo > n_virt) then
         call error%set(ERROR_VALIDATION, "the minimal basis asks for "// &
                        to_char(dims%n_vvo)//" valence-virtual orbitals but there are "// &
                        "only "//to_char(n_virt)//" virtuals. The orbital basis is "// &
                        "smaller than the minimal basis, which cannot describe the "// &
                        "valence space at all.")
         return
      end if

      ! The virtual block of the projection, and its Gram matrix.
      allocate (sigma(n_virt, size(projection, 2)))
      sigma = projection(n_occ + 1:n_mo, :)
      allocate (b(n_virt, n_virt), values(n_virt))
      call pic_gemm(sigma, sigma, b, transb="T")

      ! Negated so the eigenvalues come back ascending in magnitude-reversed
      ! order, i.e. the largest projections first once the sign is undone. The
      ! same device GAMESS uses, and it keeps the retained block at the front
      ! where every downstream index expects it.
      b = -b
      call pic_syev(b, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "the valence-virtual eigenproblem failed "// &
                        "(info = "//to_char(info)//")")
         return
      end if
      values = -values

      allocate (result%singular_values(n_virt))
      do i = 1, n_virt
         ! Negative only by rounding, at eigenvalues that are zero to begin with.
         result%singular_values(i) = sqrt(max(values(i), 0.0_dp))
      end do

      result%n_vvo = dims%n_vvo
      result%smallest_retained = result%singular_values(max(dims%n_vvo, 1))
      if (dims%n_vvo < n_virt) then
         result%largest_rejected = result%singular_values(dims%n_vvo + 1)
      else
         result%largest_rejected = 0.0_dp
      end if

      allocate (c_virt(n_ao, n_virt), rotated(n_ao, dims%n_vvo))
      c_virt = orbitals(:, n_occ + 1:n_mo)
      call pic_gemm(c_virt, b(:, 1:dims%n_vvo), rotated)
      call move_alloc(rotated, result%orbitals)

      deallocate (sigma, b, values, c_virt)
   end subroutine valence_virtual_orbitals

end module mqc_libcint_quao
