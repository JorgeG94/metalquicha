!! Atomic partial charges, by population analysis and by potential fitting
module mqc_libcint_charges
   !! Two schemes that answer the same question differently, and disagree for a
   !! reason worth knowing.
   !!
   !! **Mulliken** splits the density by which basis functions carry it, halving
   !! every overlap between the two atoms it spans. It costs one trace and is
   !! notoriously basis-set dependent -- add diffuse functions centred on one
   !! atom and its share of a neighbour's density grows, without the molecule
   !! having changed. Good for a trend, bad for a number.
   !!
   !! **CHELPG** asks what point charges would reproduce the molecule's own
   !! electrostatic potential on a grid of points outside it, and solves for
   !! them. That is a physically meaningful question -- the potential is an
   !! observable of the density -- and the answer is far less basis-set
   !! sensitive. It costs an ESP evaluation on a few thousand points, which in
   !! practice is not the expensive part: both schemes need the same SCF, and on
   !! 24 to 32 atoms in 6-31G the whole CHELPG step adds about two seconds to an
   !! eighteen-second calculation. Choosing the cheaper scheme buys ~10%, so
   !! choose on which answer you want rather than on cost.
   !!
   !! For embedding one fragment in the field of another, CHELPG is the one that
   !! matters: the whole point is to reproduce a field, and that is what it is
   !! fitted to do.
   use pic_types, only: dp
   use pic_lapack_interfaces, only: pic_getrf, pic_getrs
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_vdw_radius
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim
   use mqc_libcint_esp, only: esp_contract
   use libcint_fortran, only: LIBCINT_BAS_SLOTS, LIBCINT_ATOM_OF
   implicit none
   private

   public :: mulliken_charges
   public :: chelpg_charges
   public :: chelpg_grid

   !> Grid spacing, Angstrom. Breneman and Wiberg's value.
   real(dp), parameter :: CHELPG_SPACING = 0.3_dp
   !> How far past the van der Waals surface points are kept, Angstrom.
   !>
   !> Both a box margin and a cutoff: the lattice extends this far beyond the
   !> molecule, and a point further than this from every atom is dropped. The
   !> shell that survives is where a neighbouring molecule's electrons would
   !> actually sit, which is the region the charges have to get right.
   real(dp), parameter :: CHELPG_HEAD_SPACE = 2.8_dp
   real(dp), parameter :: ANGSTROM_TO_BOHR = 1.8897261254578281_dp

contains

   subroutine ao_to_atom(mol, owner)
      !! Which atom each basis function belongs to
      type(libcint_molecule_t), intent(in) :: mol
      integer, allocatable, intent(out) :: owner(:)

      integer :: ish, k, first, dim

      allocate (owner(mol%nao))
      do ish = 1, mol%nbas
         first = mol%shell_offset(ish)
         dim = shell_dim(mol%cartesian, ish - 1, mol%bas)
         do k = 1, dim
            ! libcint's ATOM_OF is 0-based, as is shell_offset.
            owner(first + k) = mol%bas(LIBCINT_ATOM_OF, ish) + 1
         end do
      end do
   end subroutine ao_to_atom

   subroutine mulliken_charges(mol, density, overlap, charges, error)
      !! q_A = Z_A - sum_{mu in A} (D S)_mu,mu
      !!
      !! The diagonal of `D S` is the gross population of each basis function,
      !! and summing it over an atom's functions charges that atom with every
      !! overlap it takes part in, half of which belongs to its neighbour. That
      !! halving is the whole of Mulliken's arbitrariness and the whole of its
      !! cheapness.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :), overlap(:, :)
      real(dp), allocatable, intent(out) :: charges(:)
      type(error_t), intent(inout) :: error

      integer, allocatable :: owner(:)
      real(dp) :: population
      integer :: mu, nu

      if (size(density, 1) /= mol%nao .or. size(overlap, 1) /= mol%nao) then
         call error%set(ERROR_VALIDATION, "mulliken charges: density and overlap must "// &
                        "be the size of the basis")
         return
      end if

      call ao_to_atom(mol, owner)
      allocate (charges(mol%natm))
      charges = mol%charges     !! nuclear charge to start from

      do mu = 1, mol%nao
         population = 0.0_dp
         do nu = 1, mol%nao
            population = population + density(mu, nu)*overlap(nu, mu)
         end do
         charges(owner(mu)) = charges(owner(mu)) - population
      end do
   end subroutine mulliken_charges

   subroutine chelpg_grid(mol, points, error, spacing, head_space)
      !! The shell of points a CHELPG fit is made on
      !!
      !! A cubic lattice over the molecule's bounding box, grown by
      !! `head_space`, keeping points that are outside every atom's van der
      !! Waals sphere and within `head_space` of at least one atom.
      !!
      !! Both tests matter and for different reasons. Inside a sphere the
      !! classical potential of a point charge diverges while the real one does
      !! not, so fitting there would fit nonsense. Far outside, every charge
      !! distribution with the right multipoles looks alike, so points there
      !! carry no information about where the charge sits and would dilute the
      !! points that do.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), allocatable, intent(out) :: points(:, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: spacing, head_space

      real(dp), allocatable :: keep(:, :), radii(:)
      real(dp) :: lo(3), hi(3), r(3)
      real(dp) :: h, margin, d, dmin
      integer :: nx(3)
      integer :: i, j, k, iatom, n_keep, pass
      logical :: inside

      h = CHELPG_SPACING
      if (present(spacing)) h = spacing
      margin = CHELPG_HEAD_SPACE
      if (present(head_space)) margin = head_space
      h = h*ANGSTROM_TO_BOHR
      margin = margin*ANGSTROM_TO_BOHR

      allocate (radii(mol%natm))
      do iatom = 1, mol%natm
         radii(iatom) = element_vdw_radius(nint(mol%charges(iatom)))*ANGSTROM_TO_BOHR
         if (radii(iatom) <= 0.0_dp) then
            call error%set(ERROR_VALIDATION, "chelpg: no van der Waals radius for element "// &
                           to_char(nint(mol%charges(iatom)))//", so the excluded region "// &
                           "around it is undefined")
            return
         end if
      end do

      do i = 1, 3
         lo(i) = minval(mol%coords(i, :)) - margin
         hi(i) = maxval(mol%coords(i, :)) + margin
         nx(i) = int((hi(i) - lo(i))/h) + 1
      end do

      ! Counted, then filled. The alternative is growing an array point by point
      ! over a lattice that is tens of thousands of candidates wide.
      n_keep = 0
      do pass = 1, 2
         if (pass == 2) then
            allocate (keep(3, max(n_keep, 1)))
            n_keep = 0
         end if
         do k = 0, nx(3) - 1
            do j = 0, nx(2) - 1
               do i = 0, nx(1) - 1
                  r = [lo(1) + i*h, lo(2) + j*h, lo(3) + k*h]
                  inside = .false.
                  dmin = huge(1.0_dp)
                  do iatom = 1, mol%natm
                     d = norm2(r - mol%coords(:, iatom))
                     if (d < radii(iatom)) then
                        inside = .true.
                        exit
                     end if
                     dmin = min(dmin, d)
                  end do
                  if (inside) cycle
                  if (dmin > margin) cycle
                  n_keep = n_keep + 1
                  if (pass == 2) keep(:, n_keep) = r
               end do
            end do
         end do
      end do

      if (n_keep < mol%natm + 1) then
         call error%set(ERROR_VALIDATION, "chelpg: the grid kept "//to_char(n_keep)// &
                        " points, fewer than the "//to_char(mol%natm + 1)//" unknowns "// &
                        "the fit has; use a finer spacing or a larger head space")
         return
      end if

      allocate (points(3, n_keep), source=keep(:, 1:n_keep))
   end subroutine chelpg_grid

   subroutine chelpg_charges(mol, density, charges, error, total_charge, spacing, head_space)
      !! Point charges that best reproduce the molecule's own potential
      !!
      !! Least squares over the grid, constrained to the right total charge:
      !!
      !!     minimise  sum_g ( V(g) - sum_A q_A / |r_g - R_A| )^2
      !!     subject to  sum_A q_A = Q
      !!
      !! which is a linear system in `q` and one Lagrange multiplier. The
      !! multiplier is why the matrix is `natm+1` square and why it is symmetric
      !! but not positive definite -- so it is factorised generally rather than
      !! by Cholesky.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: charges(:)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: total_charge
      real(dp), intent(in), optional :: spacing, head_space

      real(dp), allocatable :: points(:, :), v(:), inv_d(:, :), a(:, :), rhs(:, :)
      integer, allocatable :: ipiv(:)
      real(dp) :: q_total, d
      integer :: n, npt, iatom, jatom, g, info

      n = mol%natm
      q_total = 0.0_dp
      if (present(total_charge)) q_total = total_charge

      call chelpg_grid(mol, points, error, spacing, head_space)
      if (error%has_error()) return
      npt = size(points, 2)

      ! The potential the molecule actually makes: electrons from the density,
      ! nuclei by Coulomb's law. `esp_contract` returns the electronic part
      ! already negative, so the two are added rather than subtracted.
      call esp_contract(mol, points, density, v, error)
      if (error%has_error()) return

      allocate (inv_d(n, npt))
      do g = 1, npt
         do iatom = 1, n
            d = norm2(points(:, g) - mol%coords(:, iatom))
            inv_d(iatom, g) = 1.0_dp/d
            v(g) = v(g) + mol%charges(iatom)/d
         end do
      end do

      ! `rhs` is a one-column matrix because the solver takes a rank-two
      ! right-hand side; there is only ever one system here.
      allocate (a(n + 1, n + 1), rhs(n + 1, 1), ipiv(n + 1))
      a = 0.0_dp
      rhs = 0.0_dp
      do iatom = 1, n
         do jatom = 1, n
            a(iatom, jatom) = sum(inv_d(iatom, :)*inv_d(jatom, :))
         end do
         rhs(iatom, 1) = sum(v*inv_d(iatom, :))
         a(iatom, n + 1) = 1.0_dp
         a(n + 1, iatom) = 1.0_dp
      end do
      rhs(n + 1, 1) = q_total

      call pic_getrf(a, ipiv, info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "chelpg: the fit matrix is singular, which "// &
                        "means two atoms are indistinguishable to every grid point")
         return
      end if
      call pic_getrs(a, ipiv, rhs, info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "chelpg: the constrained fit would not solve")
         return
      end if

      allocate (charges(n), source=rhs(1:n, 1))
      call logger%verbose("  chelpg: fitted "//to_char(n)//" charges to "// &
                          to_char(npt)//" grid points")
   end subroutine chelpg_charges

end module mqc_libcint_charges
