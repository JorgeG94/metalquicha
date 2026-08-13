!! Distributed multipole analysis, in GAMESS's MAKEFP convention
module mqc_libcint_dma
   !! The electrostatics half of an effective fragment potential: charge through
   !! octopole at every expansion point, rather than one multipole expansion at
   !! one origin. A single-centre expansion is hopeless a few Angstrom from a
   !! molecule, which is exactly where a fragment's neighbours are.
   !!
   !! **The partition is Stone's, and simpler than its reputation.** For each pair
   !! of primitive Gaussians the product is itself a Gaussian, centred at
   !!
   !!     P = (a A + b B) / (a + b)
   !!
   !! and the whole of that pair's density goes to the nearest expansion point.
   !! Winner takes all: no overlap weighting, no Becke partition, no fitting.
   !! Ties -- within `TIE_TOLERANCE` in squared distance -- are split equally
   !! between the points that tie. GAMESS has a Becke-grid alternative behind
   !! `BIGEXP`, which defaults to zero, so a default deck never reaches it.
   !!
   !! **Why this needs primitive resolution, and what that costs.** The product
   !! centre depends on the two exponents, so different primitives of the *same*
   !! contracted shell pair can land on different expansion points. Assigning at
   !! the shell-pair level would therefore give a different answer, not a rounder
   !! one. So the density is transformed into an uncontracted primitive basis
   !! first, the assignment is made there, and the moments are accumulated over
   !! primitive pairs. `uncontract` builds that basis; it is the only real
   !! machinery here.
   !!
   !! **Moments are raw Cartesian and electronic only.** Not traceless: GAMESS's
   !! own quadrupole traces are nowhere near zero (water's oxygen comes to -12.55),
   !! because the Buckingham conversion happens in whatever consumes the potential
   !! and not in what writes it. And the nucleus enters the monopole alone -- the
   !! dipole and above are the electrons' alone -- which is what makes the sum rule
   !! below come out.
   !!
   !! **The check that matters here is the sum rule.** Translate every distributed
   !! moment to a common origin, add, and the total molecular moments come back.
   !! It is an exact identity, not an approximation, and it catches a partition
   !! that loses density. What it does *not* catch is a wrong distribution or a
   !! transposed component: a transposed quadrupole sums to a transposed total. So
   !! it is necessary and not sufficient, and the comparison against a GAMESS
   !! potential is what closes the gap -- the same lesson the distributed
   !! polarizabilities taught, where the sum rule was blind to a transpose that a
   !! per-point comparison found immediately.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim
   use mqc_libcint_multipole, only: multipole_matrices
   use libcint_fortran, only: LIBCINT_ATOM_OF, LIBCINT_ANG_OF, LIBCINT_NPRIM_OF, &
                              LIBCINT_NCTR_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_COEFF, &
                              LIBCINT_BAS_SLOTS
   implicit none
   private

   public :: dma_result_t
   public :: distributed_multipoles
   public :: expansion_points
   public :: N_QUAD, N_OCT

   !> Two expansion points count as tied when their squared distances to a
   !> product centre agree to this, in Bohr^2. GAMESS's own tolerance.
   real(dp), parameter :: TIE_TOLERANCE = 1.0e-6_dp

   !> Highest element with a covalent radius here, so a bond can be decided.
   integer, parameter :: MAX_ELEMENT = 36

   !> Unique components of a Cartesian quadrupole and octopole.
   integer, parameter :: N_QUAD = 6
   integer, parameter :: N_OCT = 10

   !> Covalent radii in Angstrom, Emsley's table as GAMESS uses it, indexed by Z.
   !>
   !> A bond is `|Ri - Rj| <= rad_i + rad_j`, evaluated in Angstrom, which is what
   !> decides where the bond midpoints go and therefore how many expansion points
   !> there are. Only through the first two rows plus the common heavier atoms --
   !> anything beyond returns a default and is flagged, because silently inventing
   !> a radius would silently invent or drop a bond.
   real(dp), parameter :: COVALENT_RADII(MAX_ELEMENT) = [ &
                          0.32_dp, 0.93_dp, &                                     ! H  He
                          1.23_dp, 0.90_dp, 0.82_dp, 0.77_dp, 0.75_dp, 0.73_dp, 0.72_dp, 0.71_dp, &  ! Li..Ne
                          1.54_dp, 1.36_dp, 1.18_dp, 1.11_dp, 1.06_dp, 1.02_dp, 0.99_dp, 0.98_dp, &  ! Na..Ar
                          2.03_dp, 1.74_dp, 1.44_dp, 1.32_dp, 1.22_dp, 1.18_dp, 1.17_dp, 1.17_dp, &  ! K..Fe
                          1.16_dp, 1.15_dp, 1.17_dp, 1.25_dp, 1.26_dp, 1.22_dp, 1.20_dp, 1.16_dp, &  ! Co..Se
                          1.14_dp, 1.12_dp]                                       ! Br Kr

   real(dp), parameter :: BOHR_PER_ANGSTROM = 1.8897261254578281_dp

   !> Component order for the packed quadrupole, as offsets into libcint's full
   !> 3x3: `XX YY ZZ XY XZ YZ`, which is what a `.efp` carries.
   integer, parameter :: QUAD_PACK(N_QUAD) = [1, 5, 9, 2, 3, 6]

   !> And the octopole, `XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ`, as offsets into
   !> the full 27 with z running fastest.
   integer, parameter :: OCT_PACK(N_OCT) = [1, 14, 27, 2, 3, 5, 15, 9, 18, 6]

   type :: dma_result_t
      !! One expansion point per column, laid out the way a `.efp` wants it
      real(dp), allocatable :: points(:, :)        !! (3, n_points), Bohr
      character(len=8), allocatable :: labels(:)   !! `A01O`, `BO21`, ...
      real(dp), allocatable :: electronic(:)       !! (n_points), negative
      real(dp), allocatable :: nuclear(:)          !! (n_points), Z on atoms else 0
      real(dp), allocatable :: dipole(:, :)        !! (3, n_points)
      real(dp), allocatable :: quadrupole(:, :)    !! (6, n_points), QUAD_PACK order
      real(dp), allocatable :: octopole(:, :)      !! (10, n_points), OCT_PACK order
   end type dma_result_t

contains

   subroutine expansion_points(atomic_numbers, coords, points, labels, nuclear, error)
      !! Every atom, then every bond midpoint
      !!
      !! Labels follow GAMESS: `A<nn><symbol>` for atoms in input order and
      !! `BO<hi><lo>` for the midpoint of atoms hi and lo, the higher index first.
      !! A midpoint carries no charge and no mass.
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coords(:, :)         !! (3, natm), Bohr
      real(dp), allocatable, intent(out) :: points(:, :)
      character(len=8), allocatable, intent(out) :: labels(:)
      real(dp), allocatable, intent(out) :: nuclear(:)
      type(error_t), intent(inout) :: error

      integer :: natm, i, j, n_bond, k
      real(dp) :: r, limit
      integer, allocatable :: bond_i(:), bond_j(:)
      character(len=8) :: text

      natm = size(atomic_numbers)
      do i = 1, natm
         if (atomic_numbers(i) < 1 .or. atomic_numbers(i) > size(COVALENT_RADII)) then
            write (text, "(i0)") atomic_numbers(i)
            call error%set(ERROR_VALIDATION, "distributed multipoles: no covalent "// &
                           "radius for element "//trim(text)//", so its bonds -- and "// &
                           "therefore its expansion points -- cannot be decided")
            return
         end if
      end do

      ! Bonds first, so the total count is known before anything is allocated.
      allocate (bond_i(natm*(natm - 1)/2), bond_j(natm*(natm - 1)/2))
      n_bond = 0
      do i = 2, natm
         do j = 1, i - 1
            r = norm2(coords(:, i) - coords(:, j))/BOHR_PER_ANGSTROM
            limit = COVALENT_RADII(atomic_numbers(i)) + COVALENT_RADII(atomic_numbers(j))
            if (r <= limit) then
               n_bond = n_bond + 1
               bond_i(n_bond) = i
               bond_j(n_bond) = j
            end if
         end do
      end do

      allocate (points(3, natm + n_bond), labels(natm + n_bond))
      allocate (nuclear(natm + n_bond))
      do i = 1, natm
         points(:, i) = coords(:, i)
         nuclear(i) = real(atomic_numbers(i), dp)
         write (labels(i), "(a,i2.2,a)") "A", i, trim(element_symbol(atomic_numbers(i)))
      end do
      do k = 1, n_bond
         i = bond_i(k)
         j = bond_j(k)
         ! The plain arithmetic mean, confirmed against GAMESS's own output to
         ! 1e-9 -- not a weighted or electronegativity-shifted midpoint.
         points(:, natm + k) = 0.5_dp*(coords(:, i) + coords(:, j))
         nuclear(natm + k) = 0.0_dp
         write (labels(natm + k), "(a,i0,i0)") "BO", i, j
      end do

      deallocate (bond_i, bond_j)
   end subroutine expansion_points

   function element_symbol(z) result(symbol)
      !! Enough of a symbol table to label an expansion point
      integer, intent(in) :: z
      character(len=2) :: symbol
      character(len=2), parameter :: table(MAX_ELEMENT) = [ &
                                     "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", &
                                     "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", &
                                     "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
                                     "Ga", "Ge", "As", "Se", "Br", "Kr"]
      if (z >= 1 .and. z <= size(table)) then
         symbol = table(z)
      else
         symbol = "X "
      end if
   end function element_symbol

   subroutine uncontract(mol, unc, transform, shell_atom, shell_exponent, error)
      !! An uncontracted copy of the basis, and the map onto the contracted one
      !!
      !! One shell per (original shell, primitive), each with a single primitive
      !! whose stored coefficient is 1 -- so libcint returns the bare primitive
      !! integral and the contraction coefficients are applied here instead, where
      !! the assignment can see them.
      !!
      !! `transform` is `(n_ao, n_unc_ao)`, carrying the contraction coefficient
      !! that links each contracted basis function to each primitive one. A shell
      !! with `nctr` contraction columns contributes `nctr` blocks of AOs and every
      !! one of them draws on the same primitives, which is exactly why the map is
      !! a matrix rather than an index list.
      type(libcint_molecule_t), intent(in) :: mol
      type(libcint_molecule_t), intent(out) :: unc
      real(dp), allocatable, intent(out) :: transform(:, :)
      integer, allocatable, intent(out) :: shell_atom(:)
      real(dp), allocatable, intent(out) :: shell_exponent(:)
      type(error_t), intent(inout) :: error

      integer :: ish, iprim, ictr, nprim, nctr, ang, ncomp, n_unc, off_env
      integer :: unc_shell, m, ao_c, ao_u, exp_ptr, coeff_ptr
      real(dp) :: alpha, coefficient

      ! How many uncontracted shells, and how much env they need.
      n_unc = 0
      do ish = 1, mol%nbas
         n_unc = n_unc + mol%bas(LIBCINT_NPRIM_OF, ish)
      end do
      if (n_unc < 1) then
         call error%set(ERROR_VALIDATION, "distributed multipoles: the basis has no shells")
         return
      end if

      unc%natm = mol%natm
      unc%cartesian = mol%cartesian
      unc%atm = mol%atm
      unc%charges = mol%charges
      unc%coords = mol%coords
      unc%nbas = n_unc
      allocate (unc%bas(LIBCINT_BAS_SLOTS, n_unc))
      unc%bas = 0
      ! The original env holds the coordinates the atm records point at, so it is
      ! kept whole and the primitive data appended after it.
      allocate (unc%env(size(mol%env) + 2*n_unc))
      unc%env = 0.0_dp
      unc%env(1:size(mol%env)) = mol%env

      allocate (shell_atom(n_unc), shell_exponent(n_unc))

      off_env = size(mol%env)
      unc_shell = 0
      do ish = 1, mol%nbas
         nprim = mol%bas(LIBCINT_NPRIM_OF, ish)
         ang = mol%bas(LIBCINT_ANG_OF, ish)
         exp_ptr = mol%bas(LIBCINT_PTR_EXP, ish)
         do iprim = 1, nprim
            unc_shell = unc_shell + 1
            alpha = mol%env(exp_ptr + iprim)

            unc%bas(LIBCINT_ATOM_OF, unc_shell) = mol%bas(LIBCINT_ATOM_OF, ish)
            unc%bas(LIBCINT_ANG_OF, unc_shell) = ang
            unc%bas(LIBCINT_NPRIM_OF, unc_shell) = 1
            unc%bas(LIBCINT_NCTR_OF, unc_shell) = 1
            unc%bas(LIBCINT_PTR_EXP, unc_shell) = off_env
            unc%env(off_env + 1) = alpha
            off_env = off_env + 1
            unc%bas(LIBCINT_PTR_COEFF, unc_shell) = off_env
            ! One, not `gto_norm`: the normalisation is already folded into the
            ! contracted basis's stored coefficients, and those are what
            ! `transform` carries. Putting it here as well would apply it twice.
            unc%env(off_env + 1) = 1.0_dp
            off_env = off_env + 1

            shell_atom(unc_shell) = mol%bas(LIBCINT_ATOM_OF, ish) + 1
            shell_exponent(unc_shell) = alpha
         end do
      end do

      ! Offsets and the AO count, in the same convention the original uses.
      allocate (unc%shell_offset(n_unc))
      unc%nao = 0
      do ish = 1, n_unc
         unc%shell_offset(ish) = unc%nao
         unc%nao = unc%nao + shell_dim(unc%cartesian, ish - 1, unc%bas)
      end do

      ! The map. Contracted AO (shell, column, component) draws on primitive AO
      ! (shell-primitive, same component) with the stored contraction coefficient.
      allocate (transform(mol%nao, unc%nao))
      transform = 0.0_dp
      unc_shell = 0
      do ish = 1, mol%nbas
         nprim = mol%bas(LIBCINT_NPRIM_OF, ish)
         nctr = mol%bas(LIBCINT_NCTR_OF, ish)
         coeff_ptr = mol%bas(LIBCINT_PTR_COEFF, ish)
         ! Components per contraction column: the shell's total divided by nctr,
         ! since libcint counts a general contraction's columns into its shell
         ! dimension and lays them out contraction-outermost.
         ncomp = shell_dim(mol%cartesian, ish - 1, mol%bas)/nctr
         do iprim = 1, nprim
            unc_shell = unc_shell + 1
            do ictr = 1, nctr
               ! (nprim, nctr) with stride nprim, as libcint reads it.
               coefficient = mol%env(coeff_ptr + (ictr - 1)*nprim + iprim)
               do m = 1, ncomp
                  ao_c = mol%shell_offset(ish) + (ictr - 1)*ncomp + m
                  ao_u = unc%shell_offset(unc_shell) + m
                  transform(ao_c, ao_u) = coefficient
               end do
            end do
         end do
      end do
   end subroutine uncontract

   subroutine distributed_multipoles(mol, density, atomic_numbers, result, error)
      !! Charge through octopole at every expansion point
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)        !! Total AO density, 2 C C^T
      integer, intent(in) :: atomic_numbers(:)
      type(dma_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: unc
      real(dp), allocatable :: transform(:, :), shell_exponent(:)
      real(dp), allocatable :: d_prim(:, :), work(:, :), ovl(:, :)
      real(dp), allocatable :: dip(:, :, :), quad(:, :, :), oct(:, :, :)
      integer, allocatable :: shell_atom(:), owner(:, :), n_owner(:, :)
      integer :: n_points, n_unc, ip, iq, k, kk, c, mu, nu, du, dq
      integer :: ou, oq, i, j
      real(dp) :: centre(3)
      real(dp) :: best, d2, weight, pop
      real(dp), allocatable :: mono(:)

      call expansion_points(atomic_numbers, mol%coords, result%points, result%labels, &
                            result%nuclear, error)
      if (error%has_error()) return
      n_points = size(result%points, 2)

      call uncontract(mol, unc, transform, shell_atom, shell_exponent, error)
      if (error%has_error()) return
      n_unc = unc%nbas

      ! The density in the primitive basis: D' = T^T D T. Every later sum runs
      ! over primitive pairs, because that is the level the partition works at.
      allocate (work(mol%nao, unc%nao), d_prim(unc%nao, unc%nao))
      call pic_gemm(density, transform, work)
      call pic_gemm(transform, work, d_prim, transa="T")

      ! Assignment. `owner` holds the winning points for each primitive shell
      ! pair, `n_owner` how many tied, so the weight is 1/n_owner.
      allocate (owner(n_points, 1), n_owner(n_unc, n_unc))
      deallocate (owner)
      allocate (owner(n_unc*n_unc, n_points))
      owner = 0
      n_owner = 0
      do iq = 1, n_unc
         do ip = 1, n_unc
            centre = (shell_exponent(ip)*mol%coords(:, shell_atom(ip)) &
                      + shell_exponent(iq)*mol%coords(:, shell_atom(iq))) &
                     /(shell_exponent(ip) + shell_exponent(iq))
            best = huge(1.0_dp)
            do k = 1, n_points
               d2 = sum((centre - result%points(:, k))**2)
               if (d2 < best) best = d2
            end do
            kk = 0
            do k = 1, n_points
               d2 = sum((centre - result%points(:, k))**2)
               if (d2 - best <= TIE_TOLERANCE) then
                  kk = kk + 1
                  owner((iq - 1)*n_unc + ip, kk) = k
               end if
            end do
            n_owner(ip, iq) = kk
         end do
      end do

      allocate (mono(n_points))
      allocate (result%electronic(n_points), result%dipole(3, n_points))
      allocate (result%quadrupole(6, n_points), result%octopole(10, n_points))
      mono = 0.0_dp
      result%dipole = 0.0_dp
      result%quadrupole = 0.0_dp
      result%octopole = 0.0_dp

      call unc%overlap(ovl)

      ! The monopole needs no origin, so it is done once outside the point loop.
      do iq = 1, n_unc
         dq = shell_dim(unc%cartesian, iq - 1, unc%bas)
         oq = unc%shell_offset(iq)
         do ip = 1, n_unc
            du = shell_dim(unc%cartesian, ip - 1, unc%bas)
            ou = unc%shell_offset(ip)
            if (n_owner(ip, iq) == 0) cycle
            weight = 1.0_dp/real(n_owner(ip, iq), dp)
            pop = 0.0_dp
            do j = 1, dq
               do i = 1, du
                  pop = pop + d_prim(ou + i, oq + j)*ovl(ou + i, oq + j)
               end do
            end do
            do kk = 1, n_owner(ip, iq)
               k = owner((iq - 1)*n_unc + ip, kk)
               mono(k) = mono(k) + weight*pop
            end do
         end do
      end do
      ! Electrons are negative; the nuclear column is separate, as a `.efp` has it.
      result%electronic = -mono

      ! Dipole and above are expanded about each point in turn, so the integrals
      ! are recomputed per point rather than translated. n_points is a handful and
      ! a translation is one more place to get a binomial wrong.
      do k = 1, n_points
         call multipole_matrices(unc, result%points(:, k), 1, dip, error)
         if (.not. error%has_error()) then
            call multipole_matrices(unc, result%points(:, k), 2, quad, error)
         end if
         if (.not. error%has_error()) then
            call multipole_matrices(unc, result%points(:, k), 3, oct, error)
         end if
         if (error%has_error()) then
            call unc%destroy()
            return
         end if

         do iq = 1, n_unc
            dq = shell_dim(unc%cartesian, iq - 1, unc%bas)
            oq = unc%shell_offset(iq)
            do ip = 1, n_unc
               du = shell_dim(unc%cartesian, ip - 1, unc%bas)
               ou = unc%shell_offset(ip)
               if (n_owner(ip, iq) == 0) cycle
               weight = 0.0_dp
               do kk = 1, n_owner(ip, iq)
                  if (owner((iq - 1)*n_unc + ip, kk) == k) then
                     weight = 1.0_dp/real(n_owner(ip, iq), dp)
                  end if
               end do
               if (weight == 0.0_dp) cycle

               do j = 1, dq
                  do i = 1, du
                     mu = ou + i
                     nu = oq + j
                     do c = 1, 3
                        result%dipole(c, k) = result%dipole(c, k) &
                                              - weight*d_prim(mu, nu)*dip(mu, nu, c)
                     end do
                     do c = 1, 6
                        result%quadrupole(c, k) = result%quadrupole(c, k) &
                                                  - weight*d_prim(mu, nu)*quad(mu, nu, QUAD_PACK(c))
                     end do
                     do c = 1, 10
                        result%octopole(c, k) = result%octopole(c, k) &
                                                - weight*d_prim(mu, nu)*oct(mu, nu, OCT_PACK(c))
                     end do
                  end do
               end do
            end do
         end do

         deallocate (dip, quad, oct)
      end do

      call unc%destroy()
      deallocate (transform, shell_atom, shell_exponent, d_prim, work, ovl)
      deallocate (owner, n_owner, mono)
   end subroutine distributed_multipoles

end module mqc_libcint_dma
