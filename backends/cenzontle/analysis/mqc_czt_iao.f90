!! Intrinsic atomic orbitals and intrinsic bond orbitals
module mqc_czt_iao
   !! Knizia's intrinsic atomic orbitals (IAOs) and the bond orbitals localised
   !! in them (IBOs). J. Chem. Theory Comput. 9, 4834 (2013).
   !!
   !! The IAOs are a minimal set of atom-centred orbitals, one per function of
   !! a free-atom minimal basis, that span the occupied space of the molecule
   !! exactly. They are polarised by the molecule, unlike the free-atom
   !! functions they start from, and they are orthonormal. A population
   !! analysis in them is therefore free of the two things that trouble
   !! Mulliken's: there is no overlap to divide, and nothing depends on how
   !! many diffuse functions the calculation's basis happens to carry.
   !!
   !! The IBOs are the occupied orbitals rotated among themselves to maximise
   !! `sum_i sum_A (Q_A^i)^p`, with `Q_A^i` the IAO population of orbital `i`
   !! on atom `A` -- Pipek-Mezey's functional on IAO charges, at the fourth
   !! power by default. Each comes out on one atom (a core or a lone pair) or
   !! two (a bond), and the populations say which.
   !!
   !! The minimal basis is `basis_sets/minao/iao.json`, Knizia's MINAO as PySCF
   !! carries it, so the construction is comparable with `pyscf.lo.iao` and
   !! `pyscf.lo.ibo` at their defaults. H through Kr, without potassium.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_number_to_symbol
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, atom_ao_blocks
   use mqc_czt_projection, only: cross_overlap
   implicit none
   private

   public :: build_iao
   public :: iao_charges
   public :: ibo_localize
   public :: ibo_populations
   public :: print_ibo_report
   public :: IBO_DEFAULT_EXPONENT

   character(len=*), parameter :: MINAO_BASIS = "minao/iao"
      !! Knizia's MINAO, the leading cc-pVTZ contractions. **Not**
      !! `minao/minao`, the ANO-RCC set the `minao` initial guess uses: two
      !! minimal bases share that name upstream, and intrinsic charges differ
      !! by a few hundredths between them.
   integer, parameter :: MINAO_Z_MAX = 36
   integer, parameter :: MINAO_Z_MISSING = 19
      !! Potassium: the one element below krypton the upstream table lacks
   integer, parameter :: IBO_DEFAULT_EXPONENT = 4
      !! The power in the localisation functional. Four is Knizia's choice and
      !! gives the more local orbitals; two is Pipek-Mezey's original.
   real(dp), parameter :: OVERLAP_FLOOR = 1.0e-10_dp
      !! Relative eigenvalue below which an overlap direction is dropped
   integer, parameter :: MAX_SWEEPS = 500
   real(dp), parameter :: GRADIENT_TOL = 1.0e-9_dp
      !! `sqrt(sum_ij B_ij^2)` over a sweep below which the localisation stops
   real(dp), parameter :: PAIR_FLOOR = 1.0e-24_dp
      !! `A_ij^2 + B_ij^2` below which a pair is left alone. Two orbitals on
      !! one atom are indistinguishable to the functional, so both vanish to
      !! rounding and the angle between them would be noise.
   integer, parameter :: LINE_LEN = 160

contains

   subroutine build_iao(mol, c_occ, iao, owner, error)
      !! Orthonormal intrinsic atomic orbitals for a set of occupied orbitals
      !!
      !! With `S1` the overlap of the calculation's basis, `S2` that of the
      !! minimal one and `S12` between them, and `P12 = S1^-1 S12`:
      !!
      !!     C~ = orth( P12 S2^-1 S21 C )          the occupied space, depolarised
      !!     A  = P12 + 2 O O~ P12 - O P12 - O~ P12,   O = C C^T S1, O~ = C~ C~^T S1
      !!
      !! and the IAOs are `A` symmetrically orthonormalised under `S1`. The
      !! occupied orbitals lie entirely inside their span.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: c_occ(:, :)
         !! (n_ao, n_occ), orthonormal under the overlap
      real(dp), allocatable, intent(out) :: iao(:, :)
         !! (n_ao, n_iao), orthonormal under the overlap
      integer, allocatable, intent(out) :: owner(:)
         !! (n_iao) the atom each IAO sits on, 1-based
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: minimal
      character(len=2), allocatable :: symbols(:)
      logical, allocatable :: ghost(:)
      integer, allocatable :: offsets(:), counts(:)
      real(dp), allocatable :: s1(:, :), s12(:, :), s2(:, :), s1_inv(:, :), s2_inv(:, :)
      real(dp), allocatable :: p12(:, :), ctild(:, :), work(:, :), o(:, :), ot(:, :)
      real(dp), allocatable :: op12(:, :), otp12(:, :), raw(:, :)
      integer :: n, m, n_occ, iatom, z

      if (error%has_error()) return
      if (.not. allocated(mol%atomic_numbers)) then
         call error%set(ERROR_VALIDATION, "intrinsic orbitals: the molecule does not "// &
                        "record its elements")
         return
      end if
      if (allocated(mol%core_electrons)) then
         if (any(mol%core_electrons > 0)) then
            call error%set(ERROR_VALIDATION, "intrinsic orbitals: not available with "// &
                           "an effective core potential, whose core electrons the "// &
                           "minimal basis would still describe")
            return
         end if
      end if
      allocate (symbols(mol%natm), ghost(mol%natm))
      do iatom = 1, mol%natm
         z = mol%atomic_numbers(iatom)
         if (z < 1 .or. z > MINAO_Z_MAX .or. z == MINAO_Z_MISSING) then
            call error%set(ERROR_VALIDATION, "intrinsic orbitals: the minimal basis "// &
                           "covers H-Kr without potassium, and this molecule has "// &
                           "an element outside it")
            return
         end if
         symbols(iatom) = element_number_to_symbol(z)
         ghost(iatom) = mol%charges(iatom) <= 0.0_dp
      end do
      call build_czt_molecule(mol%atomic_numbers, symbols, mol%coords, MINAO_BASIS, &
                              minimal, error, ghost=ghost, force_cartesian=mol%cartesian)
      if (error%has_error()) then
         call error%add_context("intrinsic orbitals: building the minimal basis")
         return
      end if

      n = mol%nao
      m = minimal%nao
      n_occ = size(c_occ, 2)
      call cross_overlap(mol, minimal, s1, s12, error)
      if (error%has_error()) return
      call minimal%overlap(s2)
      call symmetric_inverse(s1, s1_inv, error)
      call symmetric_inverse(s2, s2_inv, error)
      if (error%has_error()) return

      allocate (p12(n, m))
      call pic_gemm(s1_inv, s12, p12, beta=0.0_dp)

      ! C~ = P12 S2^-1 S21 C, then orthonormalised under S1.
      allocate (work(m, n_occ), ctild(n, n_occ))
      work = matmul(s2_inv, matmul(transpose(s12), c_occ))
      call pic_gemm(p12, work, ctild, beta=0.0_dp)
      deallocate (work)
      call lowdin_orthonormalize(ctild, s1, error)
      if (error%has_error()) return

      ! O = C C^T S1 and O~ = C~ C~^T S1, the two projectors.
      allocate (o(n, n), ot(n, n), op12(n, m), otp12(n, m), raw(n, m))
      o = matmul(c_occ, matmul(transpose(c_occ), s1))
      ot = matmul(ctild, matmul(transpose(ctild), s1))
      call pic_gemm(o, p12, op12, beta=0.0_dp)
      call pic_gemm(ot, p12, otp12, beta=0.0_dp)
      raw = p12 + 2.0_dp*matmul(o, otp12) - op12 - otp12
      call lowdin_orthonormalize(raw, s1, error)
      if (error%has_error()) return
      call move_alloc(raw, iao)

      allocate (owner(m), offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(minimal, offsets, counts)
      do iatom = 1, mol%natm
         owner(offsets(iatom) + 1:offsets(iatom) + counts(iatom)) = iatom
      end do
      call minimal%destroy()
   end subroutine build_iao

   subroutine iao_charges(mol, density, iao, owner, charges, error)
      !! q_A = Z_A minus the population of atom A's intrinsic orbitals
      !!
      !! The density is taken into the IAO basis, `A^T S D S A`, and its
      !! diagonal summed per atom. For the density the IAOs were built from
      !! the populations add up to the electron count exactly, because the
      !! occupied space lies inside their span.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)   !! (n_ao, n_ao), total
      real(dp), intent(in) :: iao(:, :)       !! From `build_iao`
      integer, intent(in) :: owner(:)
      real(dp), allocatable, intent(out) :: charges(:)   !! (natoms)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: s(:, :), sa(:, :), dsa(:, :)
      integer :: i, m

      if (error%has_error()) return
      m = size(iao, 2)
      call mol%overlap(s)
      allocate (sa(size(iao, 1), m), dsa(size(iao, 1), m))
      call pic_gemm(s, iao, sa, beta=0.0_dp)
      call pic_gemm(density, sa, dsa, beta=0.0_dp)
      allocate (charges(mol%natm))
      charges = mol%charges
      do i = 1, m
         charges(owner(i)) = charges(owner(i)) - dot_product(sa(:, i), dsa(:, i))
      end do
   end subroutine iao_charges

   subroutine ibo_localize(mol, c_occ, iao, owner, localized, error, exponent, sweeps)
      !! Rotate occupied orbitals among themselves into intrinsic bond orbitals
      !!
      !! Maximises `sum_i sum_A (Q_A^ii)^p` by Jacobi rotations of orbital
      !! pairs, `Q_A^ij` being the IAO population matrix of the pair on atom
      !! `A`. The density is unchanged: this is a rotation within the occupied
      !! space.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: c_occ(:, :)   !! (n_ao, n_occ)
      real(dp), intent(in) :: iao(:, :)
      integer, intent(in) :: owner(:)
      real(dp), allocatable, intent(out) :: localized(:, :)   !! (n_ao, n_occ)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: exponent   !! 2 or 4; default four
      integer, intent(out), optional :: sweeps    !! Jacobi sweeps taken

      real(dp), allocatable :: s(:, :), sa(:, :), c(:, :), q_ii(:), q_jj(:), q_ij(:)
      real(dp), allocatable :: ci(:), cj(:), li(:), lj(:)
      real(dp) :: aij, bij, phi, cs, sn, gradient
      integer :: n_occ, m, natm, i, j, mu, sweep, p

      if (error%has_error()) return
      p = IBO_DEFAULT_EXPONENT
      if (present(exponent)) p = exponent
      if (p /= 2 .and. p /= 4) then
         call error%set(ERROR_VALIDATION, "intrinsic bond orbitals: the localisation "// &
                        "exponent must be 2 or 4")
         return
      end if
      n_occ = size(c_occ, 2)
      m = size(iao, 2)
      natm = mol%natm
      localized = c_occ
      if (present(sweeps)) sweeps = 0
      if (n_occ < 2) return

      ! The orbitals in the IAO basis: C' = A^T S C, (n_iao, n_occ).
      call mol%overlap(s)
      allocate (sa(size(iao, 1), m), c(m, n_occ))
      call pic_gemm(s, iao, sa, beta=0.0_dp)
      call pic_gemm(sa, c_occ, c, transa="T", beta=0.0_dp)
      deallocate (s, sa)

      allocate (q_ii(natm), q_jj(natm), q_ij(natm), ci(m), cj(m))
      allocate (li(size(c_occ, 1)), lj(size(c_occ, 1)))
      do sweep = 1, MAX_SWEEPS
         gradient = 0.0_dp
         do i = 1, n_occ
            do j = 1, i - 1
               q_ii = 0.0_dp
               q_jj = 0.0_dp
               q_ij = 0.0_dp
               do mu = 1, m
                  q_ii(owner(mu)) = q_ii(owner(mu)) + c(mu, i)*c(mu, i)
                  q_jj(owner(mu)) = q_jj(owner(mu)) + c(mu, j)*c(mu, j)
                  q_ij(owner(mu)) = q_ij(owner(mu)) + c(mu, i)*c(mu, j)
               end do
               if (p == 2) then
                  aij = sum(4.0_dp*q_ij**2 - (q_ii - q_jj)**2)
                  bij = sum(4.0_dp*q_ij*(q_ii - q_jj))
               else
                  aij = sum(-q_ii**4 - q_jj**4 + 6.0_dp*(q_ii**2 + q_jj**2)*q_ij**2 &
                            + q_ii**3*q_jj + q_ii*q_jj**3)
                  bij = sum(4.0_dp*q_ij*(q_ii**3 - q_jj**3))
               end if
               if (aij*aij + bij*bij < PAIR_FLOOR) cycle
               gradient = gradient + bij*bij
               phi = 0.25_dp*atan2(bij, -aij)
               cs = cos(phi)
               sn = sin(phi)
               ci = c(:, i)
               cj = c(:, j)
               c(:, i) = cs*ci + sn*cj
               c(:, j) = -sn*ci + cs*cj
               li = localized(:, i)
               lj = localized(:, j)
               localized(:, i) = cs*li + sn*lj
               localized(:, j) = -sn*li + cs*lj
            end do
         end do
         if (present(sweeps)) sweeps = sweep
         if (sqrt(gradient) < GRADIENT_TOL) return
      end do
      call error%set(ERROR_VALIDATION, "intrinsic bond orbitals: the localisation did "// &
                     "not converge")
   end subroutine ibo_localize

   subroutine ibo_populations(mol, orbitals, iao, owner, populations, error)
      !! The IAO population of each orbital on each atom
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)   !! (n_ao, n_orb)
      real(dp), intent(in) :: iao(:, :)
      integer, intent(in) :: owner(:)
      real(dp), allocatable, intent(out) :: populations(:, :)
         !! (natoms, n_orb). Each column sums to one for an orbital inside the
         !! IAO span, which every occupied orbital is.
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: s(:, :), sa(:, :), c(:, :)
      integer :: m, n_orb, i, mu

      if (error%has_error()) return
      m = size(iao, 2)
      n_orb = size(orbitals, 2)
      call mol%overlap(s)
      allocate (sa(size(iao, 1), m), c(m, n_orb))
      call pic_gemm(s, iao, sa, beta=0.0_dp)
      call pic_gemm(sa, orbitals, c, transa="T", beta=0.0_dp)
      allocate (populations(mol%natm, n_orb))
      populations = 0.0_dp
      do i = 1, n_orb
         do mu = 1, m
            populations(owner(mu), i) = populations(owner(mu), i) + c(mu, i)**2
         end do
      end do
   end subroutine ibo_populations

   subroutine print_ibo_report(populations, charges, element_symbols, exponent, sweeps)
      !! The IAO charges, then each bond orbital and the atoms it sits on
      real(dp), intent(in) :: populations(:, :)   !! (natoms, n_orb)
      real(dp), intent(in) :: charges(:)
      character(len=*), intent(in) :: element_symbols(:)
      integer, intent(in) :: exponent, sweeps

      real(dp), parameter :: SHOWN = 0.05_dp
         !! Smallest atomic share of an orbital that is printed
      character(len=LINE_LEN) :: line
      character(len=16) :: label
      integer :: natm, n_orb, i, iatom, centres

      natm = size(populations, 1)
      n_orb = size(populations, 2)
      call logger%info("")
      call logger%info("  ==== intrinsic atomic and bond orbitals ====")
      call logger%info("")
      call logger%info("     atom      IAO charge")
      do iatom = 1, natm
         write (label, "(a,i0)") trim(adjustl(element_symbols(iatom)))//" ", iatom
         write (line, "(4x,a8,f12.6)") label, charges(iatom)
         call logger%info(trim(line))
      end do
      call logger%info("")
      write (line, "(a,i0,a,i0,a)") "    bond orbitals, exponent ", exponent, ", ", &
         sweeps, " sweeps; electrons per atom, where above 0.10:"
      call logger%info(trim(line))
      do i = 1, n_orb
         centres = count(populations(:, i) > SHOWN)
         write (line, "(i8,a,i0,a)") i, "   ", centres, "-centre:"
         do iatom = 1, natm
            if (populations(iatom, i) <= SHOWN) cycle
            write (label, "(a,i0)") trim(adjustl(element_symbols(iatom))), iatom
            write (line, "(a,2x,a,f6.3)") trim(line), trim(label), 2.0_dp*populations(iatom, i)
         end do
         call logger%info(trim(line))
      end do
      call logger%info("")
   end subroutine print_ibo_report

   subroutine symmetric_inverse(matrix, inverse, error)
      !! The inverse of a symmetric positive matrix, near-null directions dropped
      real(dp), intent(in) :: matrix(:, :)
      real(dp), allocatable, intent(out) :: inverse(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: vec(:, :), val(:), scaled(:, :)
      integer :: n, i, info

      if (error%has_error()) return
      n = size(matrix, 1)
      allocate (vec(n, n), val(n), scaled(n, n), inverse(n, n))
      vec = matrix
      call pic_syev(vec, val, jobz="V", uplo="U", info=info)
      if (info /= 0 .or. maxval(val) <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "intrinsic orbitals: an overlap matrix would "// &
                        "not diagonalise")
         return
      end if
      do i = 1, n
         if (val(i) > OVERLAP_FLOOR*maxval(val)) then
            scaled(:, i) = vec(:, i)/val(i)
         else
            scaled(:, i) = 0.0_dp
         end if
      end do
      call pic_gemm(scaled, vec, inverse, transb="T", beta=0.0_dp)
   end subroutine symmetric_inverse

   subroutine lowdin_orthonormalize(vectors, overlap, error)
      !! Replace `vectors` by `V (V^T S V)^-1/2`, orthonormal under `overlap`
      real(dp), intent(inout) :: vectors(:, :)
      real(dp), intent(in) :: overlap(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: gram(:, :), val(:), half(:, :), root(:, :), sv(:, :)
      real(dp), allocatable :: out(:, :)
      integer :: k, i, info

      if (error%has_error()) return
      k = size(vectors, 2)
      allocate (sv(size(vectors, 1), k), gram(k, k), val(k))
      call pic_gemm(overlap, vectors, sv, beta=0.0_dp)
      call pic_gemm(vectors, sv, gram, transa="T", beta=0.0_dp)
      gram = 0.5_dp*(gram + transpose(gram))
      call pic_syev(gram, val, jobz="V", uplo="U", info=info)
      if (info /= 0 .or. minval(val) <= OVERLAP_FLOOR*maxval(val)) then
         call error%set(ERROR_VALIDATION, "intrinsic orbitals: a set of vectors is "// &
                        "linearly dependent and cannot be orthonormalised")
         return
      end if
      allocate (half(k, k), root(k, k), out(size(vectors, 1), k))
      do i = 1, k
         half(:, i) = gram(:, i)/sqrt(val(i))
      end do
      call pic_gemm(half, gram, root, transb="T", beta=0.0_dp)
      call pic_gemm(vectors, root, out, beta=0.0_dp)
      vectors = out
   end subroutine lowdin_orthonormalize

end module mqc_czt_iao
