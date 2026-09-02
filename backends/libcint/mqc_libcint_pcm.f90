!! A polarizable continuum for the CPU backend
module mqc_libcint_pcm
   !! The cavity, the surface charges, and the one-electron operator they add.
   !!
   !! This is the smooth switching/Gaussian (SWIG) discretization of Lange and
   !! Herbert [J. Chem. Phys. 133, 244111 (2010)], solved either as C-PCM or as
   !! IEF-PCM, and it follows PySCF's `pyscf.solvent.pcm` term for term -- the
   !! same Lebedev points per sphere, the same switching function, the same
   !! per-point Gaussian exponents from the fitted xi of Table II in
   !! [J. Chem. Phys. 122, 194110 (2005)], the same S and D matrices. That is a
   !! deliberate choice of reference: every piece of this file can be diffed
   !! against a working Python implementation, and the validation suite does so.
   !!
   !! **What each surface point is.** Not a point charge but a normalized
   !! spherical Gaussian of exponent `zeta**2`, which is what keeps every matrix
   !! element finite as points approach each other and the cavity smooth as
   !! points cross the switching region. Its electrostatic potential at distance
   !! d is `erf(zeta*d)/d` -- exactly the attenuated `1/r` that libcint's
   !! `int1e_grids` produces when `env(PTR_RANGE_OMEGA)` holds `zeta`. So the
   !! solute-surface integrals need no fake basis functions: the points are
   !! grouped by their zeta (there are few distinct values -- one per cavity
   !! radius per Lebedev weight orbit) and each group is one batched grids call
   !! per shell pair, attenuated by its own omega.
   !!
   !! **The charge solve is direct.** K is formed once per geometry -- it depends
   !! on the cavity and the dielectric, neither of which moves during an SCF --
   !! and LU-factored; each iteration is two triangular solves. So
   !! `keywords.pcm.tolerance` and `max_iter`, which bound cuEST's conjugate
   !! gradient solve, have nothing to bound here: the direct solve meets any
   !! tolerance they could ask for.
   !!
   !! The energy is `E_diel = q . v / 2` and the Fock matrix takes the *full*
   !! potential of the charges, because the charges are determined variationally
   !! and the half does not appear in the derivative -- the same split the cuEST
   !! path documents in `system_pcm_device`.
   use pic_types, only: dp
   use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_loc, c_null_ptr
   use pic_blas_interfaces, only: pic_gemm, pic_gemv
   use pic_lapack_interfaces, only: pic_getrf, pic_getrs
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_physical_constants, only: PI
   use mqc_lebedev, only: lebedev_grid
   use mqc_pcm_radii, only: cavity_radius
   use mqc_method_config, only: pcm_config_t
   use mqc_calculation_defaults, only: DEFAULT_PCM_ZETA
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim
   use libcint_fortran, only: LIBCINT_PTR_RANGE_OMEGA
   use pic_io, only: to_char
   implicit none
   private

   public :: pcm_context_t
   public :: pcm_surface_t
   public :: build_pcm_surface
   public :: PCM_MODEL_CPCM
   public :: PCM_MODEL_IEFPCM

   integer, parameter :: PCM_MODEL_CPCM = 1
      !! Conductor-like: K = S, scaled by (eps-1)/eps
   integer, parameter :: PCM_MODEL_IEFPCM = 2
      !! Integral equation formalism: the D-matrix terms, scaled by (eps-1)/(eps+1)

   !! `PTR_GRIDS` from libcint's `cint.h`. Not exported by the Fortran
   !! interface, and 0-based like every other `PTR_*`, so it is used as `+ 1`.
   integer, parameter :: LIBCINT_PTR_GRIDS = 12

   !! The Lebedev orders the SWIG exponents are fitted for, and the fitted
   !! values. Table II of [J. Chem. Phys. 122, 194110 (2005)], copied from
   !! PySCF's `pcm.XI` so the two discretizations are the same numbers.
   !! `mqc_lebedev` also carries orders 74, 230 and 266, which the paper does
   !! not tabulate; those are refused rather than interpolated.
   integer, parameter :: N_XI_ORDERS = 29
   integer, parameter :: XI_ORDERS(N_XI_ORDERS) = [ &
                         6, 14, 26, 38, 50, 86, 110, 146, 170, 194, 302, 350, 434, 590, &
                         770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, &
                         4334, 4802, 5294, 5810]
   real(dp), parameter :: XI_VALUES(N_XI_ORDERS) = [ &
                          4.84566077868_dp, 4.86458714334_dp, 4.85478226219_dp, 4.90105812685_dp, &
                          4.89250673295_dp, 4.89741372580_dp, 4.90101060987_dp, 4.89825187392_dp, &
                          4.90685517725_dp, 4.90337644248_dp, 4.90498088169_dp, 4.86879474832_dp, &
                          4.90567349080_dp, 4.90624071359_dp, 4.90656435779_dp, 4.90685167998_dp, &
                          4.90704098216_dp, 4.90721023869_dp, 4.90733270691_dp, 4.90744499142_dp, &
                          4.90753082825_dp, 4.90760972766_dp, 4.90767282394_dp, 4.90773141371_dp, &
                          4.90777965981_dp, 4.90782469526_dp, 4.90749125553_dp, 4.90762073452_dp, &
                          4.90792902522_dp]

   !! A point whose quadrature weight times switching value falls below this is
   !! buried inside a neighbouring sphere and dropped, as PySCF drops it.
   real(dp), parameter :: BURIED_CUTOFF = 1.0e-16_dp

   type :: pcm_surface_t
      !! The discretized cavity: one smooth Gaussian per surviving Lebedev point
      integer :: n_points = 0
      integer :: n_discarded = 0            !! Buried points that were dropped
      integer, allocatable :: parent(:)     !! Owning atom, 1-based
      real(dp), allocatable :: coords(:, :)  !! (3, n), Bohr
      real(dp), allocatable :: normals(:, :)  !! (3, n), outward unit normals
      real(dp), allocatable :: zeta(:)      !! Gaussian exponent is zeta**2
      real(dp), allocatable :: weights(:)   !! Lebedev weight times 4 pi
      real(dp), allocatable :: switch_f(:)  !! SWIG switching value, in (0, 1]
      real(dp), allocatable :: area(:)      !! weight * R**2 * switch, Bohr**2
      real(dp), allocatable :: r_parent(:)  !! Cavity radius of the owning atom
   end type pcm_surface_t

   type :: pcm_group_t
      !! Surface points sharing one zeta, batched into one grids call
      !!
      !! `env` is the molecule's env with this group's point coordinates
      !! appended and `PTR_RANGE_OMEGA` set to the group's zeta, prebuilt once
      !! because nothing in it moves during the SCF.
      real(dp) :: zeta = 0.0_dp
      integer, allocatable :: idx(:)        !! Global point indices
      real(dp), allocatable :: env(:)
   end type pcm_group_t

   type :: pcm_context_t
      !! Everything the SCF iteration needs from the continuum
      logical :: enabled = .false.
      integer :: model = PCM_MODEL_CPCM
      real(dp) :: dielectric = 0.0_dp
      real(dp) :: f_eps = 0.0_dp            !! The dielectric scale factor of the model
      type(pcm_surface_t) :: surface
      real(dp), allocatable :: kmat(:, :)   !! K, LU-factored in place
      integer, allocatable :: ipiv(:)
      real(dp), allocatable :: rmat(:, :)   !! R, IEF-PCM only; C-PCM's R is -f_eps * I
      real(dp), allocatable :: v_nuc(:)     !! Nuclear potential at the points
      integer :: n_groups = 0
      type(pcm_group_t), allocatable :: groups(:)
      ! What the last solve did, for reporting beside the SCF table.
      real(dp) :: e_diel = 0.0_dp
      real(dp) :: q_total = 0.0_dp          !! Total apparent surface charge
   contains
      procedure :: build => pcm_build
      procedure :: init_model => pcm_init_model
      procedure :: build_operators => pcm_build_operators
      procedure :: solve => pcm_solve
      procedure :: operator_matrix => pcm_operator_matrix
      procedure :: destroy => pcm_destroy
   end type pcm_context_t

   interface
      function cint1e_grids_sph(buf, dims, shls, atm, natm, bas, nbas, env, opt, &
                                cache) result(ret) bind(C, name="int1e_grids_sph")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value, intent(in) :: buf, dims, shls, atm
         integer(c_int), value, intent(in) :: natm
         type(c_ptr), value, intent(in) :: bas
         integer(c_int), value, intent(in) :: nbas
         type(c_ptr), value, intent(in) :: env, opt, cache
         integer(c_int) :: ret
      end function cint1e_grids_sph

      function cint1e_grids_cart(buf, dims, shls, atm, natm, bas, nbas, env, opt, &
                                 cache) result(ret) bind(C, name="int1e_grids_cart")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value, intent(in) :: buf, dims, shls, atm
         integer(c_int), value, intent(in) :: natm
         type(c_ptr), value, intent(in) :: bas
         integer(c_int), value, intent(in) :: nbas
         type(c_ptr), value, intent(in) :: env, opt, cache
         integer(c_int) :: ret
      end function cint1e_grids_cart
   end interface

contains

   function switch_h(x) result(h)
      !! The elementary switching function, eq. 3.19 of the SWIG paper
      !!
      !! A quintic that rises smoothly from 0 to 1 across [0, 1] with two
      !! vanishing derivatives at each end. (The paper prints the polynomial
      !! with a typo; this is the corrected form PySCF uses.)
      real(dp), intent(in) :: x
      real(dp) :: h

      if (x <= 0.0_dp) then
         h = 0.0_dp
      else if (x >= 1.0_dp) then
         h = 1.0_dp
      else
         h = x**3*(10.0_dp - 15.0_dp*x + 6.0_dp*x**2)
      end if
   end function switch_h

   subroutine build_pcm_surface(atomic_numbers, coords, angular_points, radii_scale, &
                                surface, error)
      !! Tessellate the cavity: Lebedev points per atom, switched and pruned
      !!
      !! Free of the molecule type on purpose: the cavity depends on nuclei and
      !! radii, not on the basis, and the unit tests exercise it against
      !! analytic spheres with no basis in sight.
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coords(:, :)      !! (3, natm), Bohr
      integer, intent(in) :: angular_points     !! Lebedev points per sphere
      real(dp), intent(in) :: radii_scale       !! Scaling from van der Waals radii
      type(pcm_surface_t), intent(out) :: surface
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: unit_points(:, :), unit_weights(:)
      real(dp), allocatable :: radii(:), r_switch(:), r_inner(:)
      real(dp), allocatable :: pt(:, :), nrm(:, :), zet(:), wgt(:), swf(:)
      integer, allocatable :: own(:)
      real(dp) :: point(3)
      real(dp) :: xi, ratio, alpha, w, f, d, dist
      integer :: natm, iatom, jatom, ipt, n_kept, xi_index, i

      if (error%has_error()) return

      natm = size(atomic_numbers)
      if (size(coords, 1) /= 3 .or. size(coords, 2) /= natm) then
         call error%set(ERROR_VALIDATION, "PCM surface: coordinates must be (3, natm)")
         return
      end if

      ! The exponent table is fitted per Lebedev order, so an order it does not
      ! cover has no exponents to use -- refused rather than interpolated, since
      ! a wrong exponent smooths the cavity by the wrong amount and changes the
      ! solvation energy without failing.
      xi_index = 0
      do i = 1, size(XI_ORDERS)
         if (XI_ORDERS(i) == angular_points) then
            xi_index = i
            exit
         end if
      end do
      if (xi_index == 0) then
         call error%set(ERROR_VALIDATION, "keywords.pcm.angular_points = "// &
                        to_char(angular_points)//" has no fitted SWIG exponent. The "// &
                        "continuum's Gaussian widths are tabulated per Lebedev order "// &
                        "[J. Chem. Phys. 122, 194110 (2005), Table II]; use one of "// &
                        "6, 14, 26, 38, 50, 86, 110, 146, 170, 194, 302, 350, 434, 590 "// &
                        "or denser.")
         return
      end if
      xi = XI_VALUES(xi_index)

      call lebedev_grid(angular_points, unit_points, unit_weights, error)
      if (error%has_error()) return

      ! The switching region of each sphere, eq. 3.11-3.13 of the SWIG paper:
      ! a shell of width R_sw just inside the cavity radius, sized so that the
      ! region a point sweeps as it fades matches the point spacing.
      allocate (radii(natm), r_switch(natm), r_inner(natm))
      do iatom = 1, natm
         call cavity_radius(atomic_numbers(iatom), radii_scale, radii(iatom), error)
         if (error%has_error()) return
         r_switch(iatom) = radii(iatom)*sqrt(14.0_dp/real(angular_points, dp))
         ratio = radii(iatom)/r_switch(iatom)
         alpha = 0.5_dp + ratio - sqrt(ratio**2 - 1.0_dp/28.0_dp)
         r_inner(iatom) = radii(iatom) - alpha*r_switch(iatom)
      end do

      allocate (pt(3, natm*angular_points), nrm(3, natm*angular_points))
      allocate (zet(natm*angular_points), wgt(natm*angular_points))
      allocate (swf(natm*angular_points), own(natm*angular_points))

      n_kept = 0
      surface%n_discarded = 0
      do iatom = 1, natm
         do ipt = 1, angular_points
            point = radii(iatom)*unit_points(:, ipt) + coords(:, iatom)
            w = 4.0_dp*PI*unit_weights(ipt)

            ! The product of every other sphere's elementary switch. One factor
            ! at zero -- the point is inside that sphere's inner radius -- kills
            ! the product, which is what "buried" means here.
            f = 1.0_dp
            do jatom = 1, natm
               if (jatom == iatom) cycle
               dist = norm2(point - coords(:, jatom))
               d = (dist - r_inner(jatom))/r_switch(jatom)
               f = f*switch_h(d)
               if (f <= 0.0_dp) exit
            end do

            if (w*f <= BURIED_CUTOFF) then
               surface%n_discarded = surface%n_discarded + 1
               cycle
            end if

            n_kept = n_kept + 1
            pt(:, n_kept) = point
            nrm(:, n_kept) = unit_points(:, ipt)
            wgt(n_kept) = w
            swf(n_kept) = f
            own(n_kept) = iatom
            zet(n_kept) = xi/(radii(iatom)*sqrt(w))
         end do
      end do

      if (n_kept == 0) then
         call error%set(ERROR_VALIDATION, "PCM surface: every point was buried, "// &
                        "which no molecular cavity produces -- check the geometry's units")
         return
      end if

      surface%n_points = n_kept
      surface%parent = own(1:n_kept)
      surface%coords = pt(:, 1:n_kept)
      surface%normals = nrm(:, 1:n_kept)
      surface%zeta = zet(1:n_kept)
      surface%weights = wgt(1:n_kept)
      surface%switch_f = swf(1:n_kept)
      allocate (surface%area(n_kept), surface%r_parent(n_kept))
      do ipt = 1, n_kept
         surface%r_parent(ipt) = radii(own(ipt))
         surface%area(ipt) = wgt(ipt)*radii(own(ipt))**2*swf(ipt)
      end do
   end subroutine build_pcm_surface

   subroutine pcm_init_model(this, method, dielectric, error)
      !! Settle which continuum model this is and its dielectric scale factor
      class(pcm_context_t), intent(inout) :: this
      character(len=*), intent(in) :: method
      real(dp), intent(in) :: dielectric
      type(error_t), intent(inout) :: error

      if (error%has_error()) return

      ! A dielectric is required rather than defaulted: every solvent has a
      ! different one, and a default would silently pick a solvent. Same
      ! refusal the cuEST plan makes.
      if (dielectric <= 1.0_dp) then
         call error%set(ERROR_VALIDATION, "the polarizable continuum needs a solvent "// &
                        "dielectric greater than one; set keywords.pcm.dielectric. "// &
                        "There is no solvent-name table on this path, so nothing can "// &
                        "be assumed from a name.")
         return
      end if
      this%dielectric = dielectric

      select case (trim(adjustl(method)))
      case ("cpcm", "c-pcm", "CPCM", "C-PCM")
         this%model = PCM_MODEL_CPCM
         this%f_eps = (dielectric - 1.0_dp)/dielectric
      case ("iefpcm", "ief-pcm", "IEFPCM", "IEF-PCM")
         this%model = PCM_MODEL_IEFPCM
         this%f_eps = (dielectric - 1.0_dp)/(dielectric + 1.0_dp)
      case default
         call error%set(ERROR_VALIDATION, "keywords.pcm.method = '"//trim(method)// &
                        "' is not a continuum model here. The CPU path solves "// &
                        "'cpcm' (conductor-like) or 'iefpcm' (integral equation "// &
                        "formalism); they are different models with different "// &
                        "energies, so neither substitutes for the other.")
         return
      end select
   end subroutine pcm_init_model

   subroutine pcm_build_operators(this, error)
      !! Form K from the surface, factor it, and keep R where the model has one
      !!
      !! The matrices are eqs. 3.35-3.43 of the SWIG paper as PySCF writes
      !! them: S carries `erf(zeta_ij r)/r` off the diagonal and
      !! `zeta sqrt(2/pi)/F` on it; D is its normal derivative with the
      !! Gaussian correction term and `-zeta sqrt(2/pi)/(2R)` on the diagonal.
      !!
      !!     C-PCM:   K = S                          R = -f_eps
      !!     IEF-PCM: K = S - f_eps/(2 pi) D A S     R = -f_eps (1 - D A/(2 pi))
      !!
      !! Once per geometry, hence dense LAPACK rather than anything iterative:
      !! the LU is paid once and every SCF iteration is two triangular solves.
      class(pcm_context_t), intent(inout) :: this
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: smat(:, :), dmat(:, :)
      real(dp) :: dr(3)
      real(dp) :: r, nr, zi, zj, zij, xr
      integer :: n, i, j, info

      if (error%has_error()) return
      n = this%surface%n_points

      allocate (smat(n, n))
      if (this%model == PCM_MODEL_IEFPCM) allocate (dmat(n, n))

      !$omp parallel do default(none) schedule(static) &
      !$omp    shared(this, smat, dmat, n) private(i, j, dr, r, nr, zi, zj, zij, xr)
      do j = 1, n
         zj = this%surface%zeta(j)
         do i = 1, n
            if (i == j) then
               smat(j, j) = zj*sqrt(2.0_dp/PI)/this%surface%switch_f(j)
               if (allocated(dmat)) then
                  dmat(j, j) = -zj*sqrt(2.0_dp/PI)/(2.0_dp*this%surface%r_parent(j))
               end if
            else
               zi = this%surface%zeta(i)
               dr = this%surface%coords(:, i) - this%surface%coords(:, j)
               r = norm2(dr)
               zij = zi*zj/sqrt(zi**2 + zj**2)
               xr = zij*r
               smat(i, j) = erf(xr)/r
               if (allocated(dmat)) then
                  nr = dot_product(dr, this%surface%normals(:, j))
                  dmat(i, j) = smat(i, j)*nr/r**2 &
                               - 2.0_dp*xr/sqrt(PI)*exp(-xr**2)*nr/r**3
               end if
            end if
         end do
      end do
      !$omp end parallel do

      if (this%model == PCM_MODEL_IEFPCM) then
         ! D becomes D*A in place -- a column scaling, since A weights the
         ! charge each point stands for -- and then feeds both K and R.
         do j = 1, n
            dmat(:, j) = dmat(:, j)*this%surface%area(j)
         end do
         allocate (this%kmat(n, n))
         this%kmat = smat
         call pic_gemm(dmat, smat, this%kmat, alpha=-this%f_eps/(2.0_dp*PI), beta=1.0_dp)
         deallocate (smat)
         ! R = -f_eps (I - DA/(2 pi)), kept dense: IEF-PCM's right-hand side
         ! mixes the potentials, so the solve needs the whole matrix.
         call move_alloc(dmat, this%rmat)
         this%rmat = this%rmat*(this%f_eps/(2.0_dp*PI))
         do j = 1, n
            this%rmat(j, j) = this%rmat(j, j) - this%f_eps
         end do
      else
         call move_alloc(smat, this%kmat)
      end if

      allocate (this%ipiv(n))
      call pic_getrf(this%kmat, this%ipiv, info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "the PCM K matrix is singular (dgetrf info = "// &
                        to_char(info)//"), which a physical cavity does not produce -- "// &
                        "two surface points may coincide; check the geometry")
         return
      end if
   end subroutine pcm_build_operators

   subroutine pcm_solve(this, v_grids, q_sym, error)
      !! Apparent surface charges for one solute potential
      !!
      !! `K q = R v`, and for the nonsymmetric K of IEF-PCM also the transposed
      !! system, averaging the two as PySCF does -- the symmetrized charge is
      !! what makes `E = q . v / 2` a variational functional of the density.
      !! C-PCM's K is symmetric and R is a scalar, so there one solve is both.
      class(pcm_context_t), intent(inout) :: this
      real(dp), intent(in) :: v_grids(:)
      real(dp), allocatable, intent(out) :: q_sym(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: rhs(:, :), adj(:, :), qt(:)
      integer :: n, info

      if (error%has_error()) return
      n = this%surface%n_points

      allocate (rhs(n, 1))
      if (this%model == PCM_MODEL_IEFPCM) then
         rhs(:, 1) = 0.0_dp
         call pic_gemv(this%rmat, v_grids, rhs(:, 1))
      else
         rhs(:, 1) = -this%f_eps*v_grids
      end if
      call pic_getrs(this%kmat, this%ipiv, rhs, info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "PCM charge solve failed (dgetrs info = "// &
                        to_char(info)//")")
         return
      end if

      if (this%model == PCM_MODEL_IEFPCM) then
         allocate (adj(n, 1), qt(n))
         adj(:, 1) = v_grids
         call pic_getrs(this%kmat, this%ipiv, adj, trans="T", info=info)
         if (info /= 0) then
            call error%set(ERROR_VALIDATION, "PCM adjoint charge solve failed "// &
                           "(dgetrs info = "//to_char(info)//")")
            return
         end if
         qt = 0.0_dp
         call pic_gemv(this%rmat, adj(:, 1), qt, trans_a="T")
         q_sym = 0.5_dp*(rhs(:, 1) + qt)
      else
         q_sym = rhs(:, 1)
      end if
   end subroutine pcm_solve

   subroutine pcm_build(this, mol, atomic_numbers, config, error)
      !! Cavity, operators, nuclear potential and integral batches, once per geometry
      class(pcm_context_t), intent(inout) :: this
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: atomic_numbers(:)
      type(pcm_config_t), intent(in) :: config
      type(error_t), intent(inout) :: error

      integer :: ipt, iatom
      real(dp) :: d

      if (error%has_error()) return

      ! `zeta` tunes the exponent convention cuEST's plan takes, and this path
      ! does not have that freedom: its exponents are the fitted xi values the
      ! discretization was parameterized with, one per Lebedev order. A deck
      ! that moves the knob would silently get the same surface, which is the
      ! kind of ignored keyword this codebase refuses.
      if (abs(config%zeta - DEFAULT_PCM_ZETA) > 0.0_dp) then
         call error%set(ERROR_VALIDATION, "keywords.pcm.zeta tunes the Gaussian "// &
                        "switching exponents of the cuEST plan. The CPU continuum "// &
                        "uses the fitted SWIG exponents for its Lebedev order "// &
                        "[J. Chem. Phys. 122, 194110 (2005)] and has no free "// &
                        "prefactor; remove the keyword or run the cuEST backend.")
         return
      end if

      call this%init_model(config%method, config%dielectric, error)
      if (error%has_error()) return

      call build_pcm_surface(atomic_numbers, mol%coords, config%angular_points, &
                             config%radii_scale, this%surface, error)
      if (error%has_error()) return

      call pcm_build_operators(this, error)
      if (error%has_error()) return

      ! The nuclear potential each smooth point charge feels. A normalized
      ! Gaussian of exponent zeta**2 against a point nucleus integrates to
      ! erf(zeta*d)/d exactly; ghosts carry zero charge in `mol%charges` and so
      ! contribute nothing, though their spheres still shape the cavity above.
      allocate (this%v_nuc(this%surface%n_points))
      this%v_nuc = 0.0_dp
      do ipt = 1, this%surface%n_points
         do iatom = 1, mol%natm
            d = norm2(this%surface%coords(:, ipt) - mol%coords(:, iatom))
            this%v_nuc(ipt) = this%v_nuc(ipt) &
                              + mol%charges(iatom)*erf(this%surface%zeta(ipt)*d)/d
         end do
      end do

      call build_zeta_groups(this, mol, error)
      if (error%has_error()) return

      this%enabled = .true.
   end subroutine pcm_build

   subroutine build_zeta_groups(this, mol, error)
      !! Batch the surface points by zeta for the attenuated grids integrals
      !!
      !! `int1e_grids` takes one `PTR_RANGE_OMEGA` per call, so points sharing
      !! an exponent -- same cavity radius, same Lebedev weight orbit -- go in
      !! one call. A 302-point sphere has about a dozen distinct weights, so a
      !! molecule yields tens of groups, not thousands of single-point calls.
      class(pcm_context_t), intent(inout) :: this
      type(libcint_molecule_t), intent(in) :: mol
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: unique_zeta(:)
      integer, allocatable :: group_of(:), counts(:)
      integer :: n, ipt, ig, n_unique, grid_offset, m
      real(dp) :: z

      if (error%has_error()) return
      n = this%surface%n_points

      allocate (unique_zeta(n), group_of(n))
      n_unique = 0
      do ipt = 1, n
         z = this%surface%zeta(ipt)
         group_of(ipt) = 0
         do ig = 1, n_unique
            if (abs(unique_zeta(ig) - z) <= 1.0e-12_dp*z) then
               group_of(ipt) = ig
               exit
            end if
         end do
         if (group_of(ipt) == 0) then
            n_unique = n_unique + 1
            unique_zeta(n_unique) = z
            group_of(ipt) = n_unique
         end if
      end do

      this%n_groups = n_unique
      allocate (this%groups(n_unique), counts(n_unique))
      counts = 0
      do ipt = 1, n
         counts(group_of(ipt)) = counts(group_of(ipt)) + 1
      end do

      grid_offset = size(mol%env)
      do ig = 1, n_unique
         this%groups(ig)%zeta = unique_zeta(ig)
         allocate (this%groups(ig)%idx(counts(ig)))
         allocate (this%groups(ig)%env(grid_offset + 3*counts(ig)))
         this%groups(ig)%env(1:grid_offset) = mol%env
         this%groups(ig)%env(LIBCINT_PTR_RANGE_OMEGA + 1) = unique_zeta(ig)
         this%groups(ig)%env(LIBCINT_PTR_GRIDS + 1) = real(grid_offset, dp)
      end do

      counts = 0
      do ipt = 1, n
         ig = group_of(ipt)
         counts(ig) = counts(ig) + 1
         m = counts(ig)
         this%groups(ig)%idx(m) = ipt
         this%groups(ig)%env(grid_offset + 3*(m - 1) + 1:grid_offset + 3*m) = &
            this%surface%coords(:, ipt)
      end do
   end subroutine build_zeta_groups

   subroutine surface_potential(groups, n_groups, mol, n_points, density, values, error)
      !! `values(k) = sum_uv D_uv <u| erf(zeta_k |r - r_k|)/|r - r_k| |v>`
      !!
      !! The electronic part of the solute potential on the surface, contracted
      !! inside the integral loop -- the tensor this avoids storing is
      !! `(n_ao, n_ao, n_points)`, the same reasoning as `esp_contract`, whose
      !! loop this is with the group batching and the omega attenuation added.
      type(pcm_group_t), intent(in), target :: groups(:)
      integer, intent(in) :: n_groups
      type(libcint_molecule_t), intent(in), target :: mol
      integer, intent(in) :: n_points
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(out) :: values(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable, target :: buf(:)
      real(dp), allocatable :: local(:)
      integer, allocatable :: pair_i(:), pair_j(:)
      integer, target :: shls(4)
      integer :: block_max, max_group, n_pair, p, ish, jsh, di, dj, io, jo
      integer :: ig, ng, i, j, g, ret
      real(dp) :: dij, wt

      if (error%has_error()) return

      values(1:n_points) = 0.0_dp
      block_max = 1
      do ish = 1, mol%nbas
         block_max = max(block_max, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do
      max_group = 0
      do ig = 1, n_groups
         max_group = max(max_group, size(groups(ig)%idx))
      end do

      n_pair = mol%nbas*(mol%nbas + 1)/2
      allocate (pair_i(n_pair), pair_j(n_pair))
      p = 0
      do ish = 1, mol%nbas
         do jsh = 1, ish
            p = p + 1
            pair_i(p) = ish
            pair_j(p) = jsh
         end do
      end do

      !$omp parallel default(none) &
      !$omp    shared(groups, n_groups, mol, density, values, n_points, &
      !$omp           block_max, max_group, n_pair, pair_i, pair_j) &
      !$omp    private(p, ish, jsh, di, dj, io, jo, ig, ng, i, j, g, ret, &
      !$omp            shls, buf, local, dij, wt)
      allocate (buf(block_max*block_max*max_group), local(n_points))
      local = 0.0_dp
      !$omp do schedule(dynamic)
      do p = 1, n_pair
         ish = pair_i(p)
         jsh = pair_j(p)
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
         jo = mol%shell_offset(jsh)
         ! An off-diagonal block stands for itself and its transpose, counted
         ! once and doubled, exactly as `esp_contract` counts them.
         wt = 2.0_dp
         if (ish == jsh) wt = 1.0_dp

         do ig = 1, n_groups
            ng = size(groups(ig)%idx)
            shls = [ish - 1, jsh - 1, 0, ng]
            if (mol%cartesian) then
               ret = cint1e_grids_cart(c_loc(buf), c_null_ptr, c_loc(shls), &
                                       c_loc(mol%atm), mol%natm, c_loc(mol%bas), &
                                       mol%nbas, c_loc(groups(ig)%env), c_null_ptr, &
                                       c_null_ptr)
            else
               ret = cint1e_grids_sph(c_loc(buf), c_null_ptr, c_loc(shls), &
                                      c_loc(mol%atm), mol%natm, c_loc(mol%bas), &
                                      mol%nbas, c_loc(groups(ig)%env), c_null_ptr, &
                                      c_null_ptr)
            end if
            if (ret == 0) cycle
            do j = 1, dj
               do i = 1, di
                  dij = wt*density(io + i, jo + j)
                  if (dij == 0.0_dp) cycle
                  do g = 1, ng
                     local(groups(ig)%idx(g)) = local(groups(ig)%idx(g)) &
                                                + dij*buf(g + (i - 1)*ng + (j - 1)*ng*di)
                  end do
               end do
            end do
         end do
      end do
      !$omp end do
      !$omp critical(mqc_pcm_potential_accumulate)
      values(1:n_points) = values(1:n_points) + local
      !$omp end critical(mqc_pcm_potential_accumulate)
      deallocate (buf, local)
      !$omp end parallel
   end subroutine surface_potential

   subroutine surface_fock(groups, n_groups, mol, q, vmat, error)
      !! `vmat(u,v) = -sum_k q_k <u| erf(zeta_k |r - r_k|)/|r - r_k| |v>`
      !!
      !! The one-electron operator the surface charges add to the Fock matrix.
      !! The minus sign is the electron's charge against the apparent charges'
      !! potential. Same loop as `surface_potential` contracted the other way;
      !! a shell pair owns its block of `vmat` outright, so the threads share
      !! nothing but the read-only tables.
      type(pcm_group_t), intent(in), target :: groups(:)
      integer, intent(in) :: n_groups
      type(libcint_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: vmat(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable, target :: buf(:)
      real(dp), allocatable :: blk(:, :)
      integer, allocatable :: pair_i(:), pair_j(:)
      integer, target :: shls(4)
      integer :: block_max, max_group, n_pair, p, ish, jsh, di, dj, io, jo
      integer :: ig, ng, i, j, g, ret

      if (error%has_error()) return

      vmat = 0.0_dp
      block_max = 1
      do ish = 1, mol%nbas
         block_max = max(block_max, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do
      max_group = 0
      do ig = 1, n_groups
         max_group = max(max_group, size(groups(ig)%idx))
      end do

      n_pair = mol%nbas*(mol%nbas + 1)/2
      allocate (pair_i(n_pair), pair_j(n_pair))
      p = 0
      do ish = 1, mol%nbas
         do jsh = 1, ish
            p = p + 1
            pair_i(p) = ish
            pair_j(p) = jsh
         end do
      end do

      !$omp parallel default(none) &
      !$omp    shared(groups, n_groups, mol, q, vmat, block_max, max_group, &
      !$omp           n_pair, pair_i, pair_j) &
      !$omp    private(p, ish, jsh, di, dj, io, jo, ig, ng, i, j, g, ret, &
      !$omp            shls, buf, blk)
      allocate (buf(block_max*block_max*max_group), blk(block_max, block_max))
      !$omp do schedule(dynamic)
      do p = 1, n_pair
         ish = pair_i(p)
         jsh = pair_j(p)
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
         jo = mol%shell_offset(jsh)

         blk(1:di, 1:dj) = 0.0_dp
         do ig = 1, n_groups
            ng = size(groups(ig)%idx)
            shls = [ish - 1, jsh - 1, 0, ng]
            if (mol%cartesian) then
               ret = cint1e_grids_cart(c_loc(buf), c_null_ptr, c_loc(shls), &
                                       c_loc(mol%atm), mol%natm, c_loc(mol%bas), &
                                       mol%nbas, c_loc(groups(ig)%env), c_null_ptr, &
                                       c_null_ptr)
            else
               ret = cint1e_grids_sph(c_loc(buf), c_null_ptr, c_loc(shls), &
                                      c_loc(mol%atm), mol%natm, c_loc(mol%bas), &
                                      mol%nbas, c_loc(groups(ig)%env), c_null_ptr, &
                                      c_null_ptr)
            end if
            if (ret == 0) cycle
            do j = 1, dj
               do i = 1, di
                  do g = 1, ng
                     blk(i, j) = blk(i, j) &
                                 + q(groups(ig)%idx(g))*buf(g + (i - 1)*ng + (j - 1)*ng*di)
                  end do
               end do
            end do
         end do

         ! A pair owns both its block and the transposed one, so this write is
         ! race-free across the pair loop. On the diagonal they coincide.
         do j = 1, dj
            do i = 1, di
               vmat(io + i, jo + j) = -blk(i, j)
            end do
         end do
         if (ish /= jsh) then
            do j = 1, dj
               do i = 1, di
                  vmat(jo + j, io + i) = -blk(i, j)
               end do
            end do
         end if
      end do
      !$omp end do
      deallocate (buf, blk)
      !$omp end parallel
   end subroutine surface_fock

   subroutine pcm_operator_matrix(this, mol, density, vmat, energy, error)
      !! One SCF iteration's continuum: charges from this density, then their operator
      !!
      !! The caller adds `vmat` -- the full potential -- to its Fock matrix and
      !! `energy` -- which already carries the half -- to its electronic energy,
      !! the same contract `system_pcm_device` keeps on the GPU path.
      class(pcm_context_t), intent(inout) :: this
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)   !! Total density (alpha plus beta)
      real(dp), intent(out) :: vmat(:, :)
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: v_grids(:), q_sym(:)

      energy = 0.0_dp
      if (error%has_error()) return

      allocate (v_grids(this%surface%n_points))
      call surface_potential(this%groups, this%n_groups, mol, this%surface%n_points, &
                             density, v_grids, error)
      if (error%has_error()) return
      ! Nuclei push the potential up, electrons pull it down.
      v_grids = this%v_nuc - v_grids

      call this%solve(v_grids, q_sym, error)
      if (error%has_error()) return

      this%e_diel = 0.5_dp*dot_product(q_sym, v_grids)
      this%q_total = sum(q_sym)
      energy = this%e_diel

      call surface_fock(this%groups, this%n_groups, mol, q_sym, vmat, error)
   end subroutine pcm_operator_matrix

   subroutine pcm_destroy(this)
      !! Release everything; the context is reusable after another `build`
      class(pcm_context_t), intent(inout) :: this

      this%enabled = .false.
      this%n_groups = 0
      this%e_diel = 0.0_dp
      this%q_total = 0.0_dp
      if (allocated(this%kmat)) deallocate (this%kmat)
      if (allocated(this%ipiv)) deallocate (this%ipiv)
      if (allocated(this%rmat)) deallocate (this%rmat)
      if (allocated(this%v_nuc)) deallocate (this%v_nuc)
      if (allocated(this%groups)) deallocate (this%groups)
      this%surface = pcm_surface_t()
   end subroutine pcm_destroy

end module mqc_libcint_pcm
