!! AVAS: choosing an active space by projecting onto atomic valence orbitals
module mqc_libcint_avas
   !! Which orbitals go in an active space, chosen mechanically: say which
   !! *atomic* orbitals the chemistry lives in -- "the nitrogen 2p", "the
   !! chromium 3d" -- and let the projection decide which molecular orbitals
   !! carry that character. `valence_select` answers the same question by taking
   !! the whole valence shell instead, which is a larger space and no judgement
   !! about the molecule.
   !!
   !! Sayfutyarova, Sun, Chan and Knizia, JCTC 13, 4063 (2017).
   !!
   !! The construction is a projector. With `P` the space spanned by the chosen
   !! atomic orbitals, form
   !!
   !!     sigma = C^T S P S C
   !!
   !! block by block over the occupied and virtual molecular orbitals, and
   !! diagonalize each block. An eigenvalue near one means that combination of
   !! molecular orbitals is almost entirely the atomic character asked for; near
   !! zero means it has none. The active space is everything above a threshold,
   !! and because the spectrum of a sensible request is bimodal, the threshold
   !! does not have to be chosen carefully -- 0.2 is the published default and
   !! usually sits in a gap rather than through a cluster.
   !!
   !! **The reference basis is this code's own accurate atomic minimal basis,
   !! not the one the paper used.** AVAS as published projects onto MINAO, where
   !! `aambs.json` is the free-atom set transcribed from GAMESS, so the
   !! eigenvalues differ. The *selection* should not: a well-posed request
   !! separates near zero and near one, and a cut at 0.2 lands in the same place
   !! either way. The tests therefore check the active space and the converged
   !! CASSCF energy against PySCF, not the intermediate weights.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_io, only: to_char
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_aambs, only: aambs_shell_labels, aambs_dimensions, aambs_dimensions_t
   use mqc_libcint_integrals, only: libcint_molecule_t, mixed_basis_overlap
   use mqc_libcint_quao, only: build_aambs_molecule, mo_aambs_overlap, &
                               valence_virtual_orbitals, vvo_result_t
   implicit none
   private

   public :: avas_select
   public :: valence_select
   public :: avas_result_t
   public :: parse_orbital_label

   real(dp), parameter :: DEFAULT_THRESHOLD = 0.2_dp
      !! The published default. Robust because the eigenvalue spectrum of a
      !! sensible request is bimodal, not because 0.2 is special.

   type :: avas_result_t
      !! The active space AVAS chose, and the orbitals to run it in
      integer :: n_inactive = 0
      integer :: n_active = 0
      integer :: n_active_electrons = 0
      integer :: n_active_occupied = 0     !! Of the active orbitals, how many came from the occupied space
      real(dp), allocatable :: orbitals(:, :)
         !! (n_ao, n_mo), reordered so the active orbitals are contiguous and
         !! begin at `n_inactive + 1`, which is what the CASSCF driver assumes.
      real(dp), allocatable :: occupied_weights(:), virtual_weights(:)
         !! The projector eigenvalues, ascending within each block. They say
         !! whether the cut fell in a gap or through a cluster.
   end type avas_result_t

contains

   subroutine parse_orbital_label(label, symbol, principal, angular, error)
      !! `"N 2p"` into a symbol, a principal quantum number and an l
      !!
      !! Strict about the shape: a label that cannot be parsed is refused rather
      !! than ignored.
      character(len=*), intent(in) :: label
      character(len=:), allocatable, intent(out) :: symbol
      integer, intent(out) :: principal, angular
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: text, shell
      integer :: space, digits, i

      if (error%has_error()) return
      principal = 0
      angular = -1
      text = trim(adjustl(label))
      space = index(text, " ")
      if (space <= 1 .or. space >= len(text)) then
         call error%set(ERROR_VALIDATION, "'"//trim(label)//"' is not an orbital "// &
                        "label. They read like 'N 2p' or 'Cr 3d': an element "// &
                        "symbol, a space, a principal quantum number and a "// &
                        "subshell letter.")
         return
      end if
      symbol = text(1:space - 1)
      shell = trim(adjustl(text(space + 1:)))

      digits = 0
      do i = 1, len(shell)
         if (shell(i:i) < "0" .or. shell(i:i) > "9") exit
         digits = i
      end do
      if (digits == 0 .or. digits >= len(shell)) then
         call error%set(ERROR_VALIDATION, "'"//trim(label)//"' has no principal "// &
                        "quantum number and subshell letter after the element "// &
                        "symbol; it should read like '2p'.")
         return
      end if
      read (shell(1:digits), *) principal

      select case (shell(digits + 1:digits + 1))
      case ("s", "S")
         angular = 0
      case ("p", "P")
         angular = 1
      case ("d", "D")
         angular = 2
      case ("f", "F")
         angular = 3
      case default
         call error%set(ERROR_VALIDATION, "'"//trim(label)//"' names the subshell '"// &
                        shell(digits + 1:digits + 1)//"'. Known subshells are "// &
                        "s, p, d and f.")
      end select
   end subroutine parse_orbital_label

   subroutine valence_select(mol, atomic_numbers, element_symbols, coordinates, &
                             orbitals, n_electrons, result, error, verbose)
      !! The whole valence shell as an active space
      !!
      !! Every occupied orbital that is not core, plus the valence-virtual
      !! orbitals that complete the minimal basis. No threshold and no labels:
      !! the size is fixed by counting the free-atom minimal basis, so the same
      !! molecule gives the same space in any basis set.
      !!
      !! The valence-virtual orbitals come from `mqc_libcint_quao`, which
      !! extracts them for the bonding analysis. No threshold decides how many:
      !! `n_vvo = n_mbs - n_occupied`, exactly.
      !!
      !! **It gets big quickly.** Nitrogen is CAS(10,8) and water is CAS(8,6),
      !! but the count grows with the molecule and not with what is interesting
      !! in it, so a full valence space on anything of size is past what a
      !! complete expansion can hold. `keywords.mcscf.ormas` restricts the
      !! occupations of a space chosen this way.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      real(dp), intent(in) :: orbitals(:, :)   !! Converged reference orbitals
      integer, intent(in) :: n_electrons
      type(avas_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose

      type(libcint_molecule_t) :: aambs
      type(aambs_dimensions_t) :: dims
      type(vvo_result_t) :: vvo
      real(dp), allocatable :: s_mbs(:, :), mixed(:, :), projection(:, :)
      character(len=160) :: line
      integer :: n_ao, n_mo, filled
      logical :: loud

      if (error%has_error()) return
      loud = .false.
      if (present(verbose)) loud = verbose

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)

      ! Counting first: it refuses an element past xenon or an odd electron
      ! count before any integral is built.
      call aambs_dimensions(atomic_numbers, n_electrons, dims, error)
      if (error%has_error()) return

      call build_aambs_molecule(atomic_numbers, element_symbols, coordinates, aambs, error)
      if (error%has_error()) return
      call aambs%overlap(s_mbs)
      call mixed_basis_overlap(mol, aambs, mixed, error)
      if (error%has_error()) return
      call mo_aambs_overlap(orbitals, mixed, s_mbs, projection, error)
      call valence_virtual_orbitals(orbitals, projection, dims, vvo, error)
      if (error%has_error()) return

      result%n_inactive = dims%n_core
      result%n_active = dims%n_valocc + vvo%n_vvo
      result%n_active_occupied = dims%n_valocc
      result%n_active_electrons = 2*dims%n_valocc

      ! Core, then the valence occupied, then the valence virtuals, then
      ! whatever is left. The active block comes out contiguous from
      ! `n_inactive + 1` without any reordering, because that is already the
      ! order the three pieces are in.
      allocate (result%orbitals(n_ao, n_mo))
      result%orbitals(:, 1:dims%n_occupied) = orbitals(:, 1:dims%n_occupied)
      filled = dims%n_occupied
      if (vvo%n_vvo > 0) then
         result%orbitals(:, filled + 1:filled + vvo%n_vvo) = vvo%orbitals
         filled = filled + vvo%n_vvo
      end if
      if (filled < n_mo) then
         result%orbitals(:, filled + 1:n_mo) = vvo%external_orbitals
      end if

      ! The gap between the smallest kept and largest rejected singular value is
      ! the diagnostic the paper reports: a clean valence space keeps values at
      ! essentially one and rejects ones far below. No gap means the minimal
      ! basis is not finding a valence space.
      if (loud) then
         call logger%info("")
         call logger%info("  full valence active space")
         write (line, "(a,i0,a,i0,a)") "    active space                CAS(", &
            result%n_active_electrons, ",", result%n_active, ")"
         call logger%info(trim(line))
         write (line, "(a,i0)") "    inactive (core) orbitals    ", result%n_inactive
         call logger%info(trim(line))
         write (line, "(a,i0,a,i0)") "    valence occupied            ", &
            dims%n_valocc, "   valence virtual  ", vvo%n_vvo
         call logger%info(trim(line))
         write (line, "(a,f10.6,a,f10.6)") "    smallest kept               ", &
            vvo%smallest_retained, "   largest rejected ", vvo%largest_rejected
         call logger%info(trim(line))
      end if

      call aambs%destroy()
   end subroutine valence_select

   subroutine avas_select(mol, atomic_numbers, element_symbols, coordinates, orbitals, &
                          n_occupied, labels, result, error, threshold, verbose)
      !! Choose an active space from a list of atomic orbital labels
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      real(dp), intent(in) :: orbitals(:, :)     !! Converged reference orbitals
      integer, intent(in) :: n_occupied          !! Doubly occupied count
      character(len=*), intent(in) :: labels(:)  !! e.g. ["N 2s", "N 2p"]
      type(avas_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: threshold
      logical, intent(in), optional :: verbose

      type(libcint_molecule_t) :: aambs
      real(dp), allocatable :: s_mbs(:, :), mixed(:, :)
      real(dp), allocatable :: reference(:, :), metric(:, :), inverse(:, :)
      real(dp), allocatable :: projected(:, :), sigma(:, :), work(:, :)
      real(dp), allocatable :: block_matrix(:, :), values(:)
      real(dp), allocatable :: chosen(:, :)
      integer, allocatable :: atom_of(:), principal(:), angular(:), selected(:)
      character(len=:), allocatable :: symbol
      character(len=160) :: line
      real(dp) :: cut
      integer :: n_ao, n_mo, n_mbs, n_ref, n_virtual, i, j, k, info
      integer :: want_n, want_l, n_keep_occ, n_keep_vir, cursor
      logical :: loud
      ! TODO(mqc): `info` is declared and never used here; `block_eigen` and
      ! `symmetric_inverse` each have their own.

      if (error%has_error()) return
      cut = DEFAULT_THRESHOLD
      if (present(threshold)) cut = threshold
      loud = .false.
      if (present(verbose)) loud = verbose

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_virtual = n_mo - n_occupied
      if (n_occupied <= 0 .or. n_virtual < 0) then
         call error%set(ERROR_VALIDATION, "AVAS needs an occupied and a virtual "// &
                        "space to divide; this reference has "//to_char(n_occupied)// &
                        " occupied of "//to_char(n_mo)//" orbitals.")
         return
      end if

      ! ---- which minimal-basis functions were asked for ----------------------
      call build_aambs_molecule(atomic_numbers, element_symbols, coordinates, aambs, error)
      if (error%has_error()) return
      call aambs%overlap(s_mbs)
      call mixed_basis_overlap(mol, aambs, mixed, error)
      call aambs_shell_labels(atomic_numbers, atom_of, principal, angular, error)
      if (error%has_error()) return

      n_mbs = size(s_mbs, 1)
      if (size(atom_of) /= n_mbs) then
         call error%set(ERROR_VALIDATION, "the minimal basis has "//to_char(n_mbs)// &
                        " functions but the shell table describes "// &
                        to_char(size(atom_of))//". The basis file and its shell "// &
                        "list disagree.")
         return
      end if

      allocate (selected(n_mbs))
      n_ref = 0
      do k = 1, size(labels)
         if (len_trim(labels(k)) == 0) cycle
         call parse_orbital_label(labels(k), symbol, want_n, want_l, error)
         if (error%has_error()) return
         do i = 1, n_mbs
            if (principal(i) /= want_n .or. angular(i) /= want_l) cycle
            if (trim(adjustl(element_symbols(atom_of(i)))) /= trim(symbol)) cycle
            if (any(selected(1:n_ref) == i)) cycle
            n_ref = n_ref + 1
            selected(n_ref) = i
         end do
      end do

      if (n_ref == 0) then
         call error%set(ERROR_VALIDATION, "none of the requested orbital labels "// &
                        "match anything in this molecule's minimal basis. Check "// &
                        "the element symbols and that the shells exist -- a free "// &
                        "atom's minimal basis has only its occupied shells, so "// &
                        "there is no '2d' on nitrogen to ask for.")
         return
      end if

      ! ---- the projector, in the molecular orbital basis ---------------------
      !
      !     sigma = C^T S P S C  with  P = R (R^T S R)^-1 R^T
      !
      ! built as `s21^T s2^-1 s21`, which is the same thing with the reference
      ! overlaps done once.
      allocate (reference(n_ref, n_mo), metric(n_ref, n_ref))
      do j = 1, n_ref
         do i = 1, n_ref
            metric(i, j) = s_mbs(selected(i), selected(j))
         end do
      end do
      allocate (work(n_ao, n_ref))
      do j = 1, n_ref
         work(:, j) = mixed(:, selected(j))
      end do
      call pic_gemm(work, orbitals, reference, transa="T")
      deallocate (work)

      call symmetric_inverse(metric, inverse, error)
      if (error%has_error()) return

      allocate (projected(n_ref, n_mo), sigma(n_mo, n_mo))
      call pic_gemm(inverse, reference, projected)
      call pic_gemm(reference, projected, sigma, transa="T")

      ! ---- diagonalize each block, keep what is above the threshold ----------
      allocate (chosen(n_ao, n_mo))
      cursor = 0

      call block_eigen(sigma, 1, n_occupied, block_matrix, values, error)
      if (error%has_error()) return
      result%occupied_weights = values
      n_keep_occ = count(values >= cut)
      ! Ascending, so the rejected ones come first and become inactive.
      call rotate(orbitals(:, 1:n_occupied), block_matrix(:, 1:n_occupied - n_keep_occ), &
                  chosen, cursor)
      result%n_inactive = cursor
      call rotate(orbitals(:, 1:n_occupied), &
                  block_matrix(:, n_occupied - n_keep_occ + 1:n_occupied), chosen, cursor)

      call block_eigen(sigma, n_occupied + 1, n_mo, block_matrix, values, error)
      if (error%has_error()) return
      result%virtual_weights = values
      n_keep_vir = count(values >= cut)
      ! Descending would be tidier here, but `pic_syev` gives ascending and the
      ! kept ones are the tail, so the two slices below are the other way round.
      call rotate(orbitals(:, n_occupied + 1:n_mo), &
                  block_matrix(:, n_virtual - n_keep_vir + 1:n_virtual), chosen, cursor)
      call rotate(orbitals(:, n_occupied + 1:n_mo), &
                  block_matrix(:, 1:n_virtual - n_keep_vir), chosen, cursor)

      result%n_active = n_keep_occ + n_keep_vir
      result%n_active_occupied = n_keep_occ
      result%n_active_electrons = 2*n_keep_occ
      result%orbitals = chosen

      if (result%n_active == 0) then
         call error%set(ERROR_VALIDATION, "AVAS selected no orbitals at a threshold "// &
                        "of "//to_char(cut)//". Either the labels name orbitals this "// &
                        "molecule's occupied space does not use, or the threshold is "// &
                        "too high.")
         return
      end if

      if (loud) then
         call logger%info("")
         call logger%info("  AVAS active space selection")
         write (line, "(a,i0,a)") "    reference orbitals   ", n_ref, &
            " minimal-basis functions"
         call logger%info(trim(line))
         write (line, "(a,f8.3)") "    threshold            ", cut
         call logger%info(trim(line))
         write (line, "(a,i0,a,i0,a)") "    active space         CAS(", &
            result%n_active_electrons, ",", result%n_active, ")"
         call logger%info(trim(line))
         write (line, "(a,i0,a,i0,a)") "    from                 ", n_keep_occ, &
            " occupied and ", n_keep_vir, " virtual orbitals"
         call logger%info(trim(line))
         write (line, "(a,i0)") "    inactive             ", result%n_inactive
         call logger%info(trim(line))
         ! The gap the cut fell in, which says whether the request was well posed.
         call report_gap(result%occupied_weights, cut, "occupied")
         call report_gap(result%virtual_weights, cut, "virtual")
      end if

      call aambs%destroy()
      deallocate (s_mbs, mixed, reference, metric, inverse, projected, sigma)
      deallocate (atom_of, principal, angular, selected, chosen)
   end subroutine avas_select

   subroutine rotate(source, transform, target_block, cursor)
      !! Append `source * transform` to the running orbital set
      real(dp), intent(in) :: source(:, :)
      real(dp), intent(in) :: transform(:, :)
      real(dp), intent(inout) :: target_block(:, :)
      integer, intent(inout) :: cursor

      integer :: n

      n = size(transform, 2)
      if (n == 0) return
      call pic_gemm(source, transform, target_block(:, cursor + 1:cursor + n))
      cursor = cursor + n
   end subroutine rotate

   subroutine block_eigen(matrix, first, last, vectors, values, error)
      !! Eigenpairs of one diagonal block, ascending
      real(dp), intent(in) :: matrix(:, :)
      integer, intent(in) :: first, last
      real(dp), allocatable, intent(out) :: vectors(:, :)
      real(dp), allocatable, intent(out) :: values(:)
      type(error_t), intent(inout) :: error

      integer :: n, info

      if (error%has_error()) return
      n = last - first + 1
      if (allocated(vectors)) deallocate (vectors)
      if (allocated(values)) deallocate (values)
      allocate (vectors(n, n), values(n))
      vectors = matrix(first:last, first:last)
      call pic_syev(vectors, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "an AVAS projector block could not be "// &
                        "diagonalized (info = "//to_char(info)//")")
      end if
   end subroutine block_eigen

   subroutine symmetric_inverse(matrix, inverse, error)
      !! The inverse of a symmetric positive definite matrix
      real(dp), intent(in) :: matrix(:, :)
      real(dp), allocatable, intent(out) :: inverse(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: vectors(:, :), values(:), scaled(:, :)
      integer :: n, i, info

      if (error%has_error()) return
      n = size(matrix, 1)
      allocate (vectors(n, n), values(n), scaled(n, n), inverse(n, n))
      vectors = matrix
      call pic_syev(vectors, values, jobz="V", uplo="U", info=info)
      if (info /= 0 .or. values(1) <= 1.0e-12_dp) then
         call error%set(ERROR_VALIDATION, "the requested reference orbitals are "// &
                        "linearly dependent, so the projection onto them is not "// &
                        "defined. Two atoms may be on top of one another.")
         return
      end if
      do i = 1, n
         scaled(:, i) = vectors(:, i)/values(i)
      end do
      call pic_gemm(scaled, vectors, inverse, transb="T")
      deallocate (vectors, values, scaled)
   end subroutine symmetric_inverse

   subroutine report_gap(weights, cut, name)
      !! How cleanly the threshold separated the two groups
      real(dp), intent(in) :: weights(:)
      real(dp), intent(in) :: cut
      character(len=*), intent(in) :: name

      character(len=160) :: line
      real(dp) :: below, above
      integer :: i

      below = -1.0_dp
      above = 2.0_dp
      do i = 1, size(weights)
         if (weights(i) < cut) below = max(below, weights(i))
         if (weights(i) >= cut) above = min(above, weights(i))
      end do
      if (below < 0.0_dp .or. above > 1.5_dp) return
      write (line, "(a,a,a,f7.4,a,f7.4)") "    ", name, " gap          ", below, &
         " to ", above
      call logger%info(trim(line))
   end subroutine report_gap

end module mqc_libcint_avas
