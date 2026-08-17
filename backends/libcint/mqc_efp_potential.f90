!! A complete effective fragment potential, computed and written here
module mqc_efp_potential
   !! `RUNTYP=MAKEFP`: a geometry and a basis name in, a `.efp` file out.
   !!
   !! Each parameter block already had a `validation/check_*` program computing it
   !! and dumping it for comparison against GAMESS's printed numbers. This is the
   !! assembly -- the SCF, the localization, the multipoles, the static and dynamic
   !! polarizabilities, the projection data and both screening fits, written as one
   !! file -- and it is where a run turns into a fragment someone else can use.
   !!
   !! **All seventeen sections GAMESS's reader recognises are written**, so the
   !! potential is complete: electrostatics to octupole with charge-penetration
   !! screening, polarization, exchange repulsion, dispersion through `E6`, `E7` and
   !! `E8`, and charge transfer. Seventeen and not eighteen because `CTFOK` is a
   !! *subsection* of `CTVEC` rather than a section: GAMESS accepts it only directly
   !! behind one and aborts on a standalone one, so the two share a `STOP`.
   !!
   !! **The formats target GAMESS's reader, not byte-identity with its writer.** The
   !! test that matters is whether GAMESS accepts the file and agrees with the
   !! energies it computes from it, which `tools/efp_validation/dimer_energy.py`
   !! asks: every term agrees, exchange repulsion and charge transfer exactly and the
   !! rest between 1e-10 and 2e-07.
   !!
   !! **`LMOQQPOL` is the one block validated by its energy rather than its values.**
   !! Our per-orbital tensors differ from GAMESS's written ones by up to 0.162, and
   !! deliberately so: our response summed over the orbitals reproduces GAMESS's own
   !! *molecular* quadrupole-quadrupole polarizability -- which
   !! `$MAKEFP MOLPOL=.TRUE.` writes with no translation applied -- to 1.6e-05, while
   !! its written per-orbital block misses that same tensor by 1.98e-02. The
   !! difference is confined to the part antisymmetric under exchanging the two index
   !! pairs, which cancels on summing, and it moves the interaction energy by 2e-07.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_mass
   use mqc_cgto, only: molecular_basis_type
   use mqc_basis_utils, only: find_basis_file
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule, shell_dim
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_atomic_guess, only: build_restricted_guess, guess_display_name
   use mqc_libcint_localize, only: boys_localize
   use mqc_libcint_dma, only: dma_result_t, distributed_multipoles
   use mqc_libcint_cphf, only: response_hessian_t, distributed_polarizability, &
                               distributed_dynamic_polarizability, &
                               distributed_dynamic_cross, &
                               casimir_polder_frequencies, N_CASIMIR_POLDER
   use mqc_libcint_multipole, only: multipole_matrices
   use mqc_libcint_screening, only: fit_screening, screening_target_t, &
                                    SCREEN_EXPONENTIAL, SCREEN_GAUSSIAN
   use pic_timer, only: timer_type
   use libcint_fortran, only: LIBCINT_ANG_OF
   use pic_logger, only: logger => global_logger
   use mqc_program_limits, only: MAX_LINE_LENGTH
   implicit none
   private

   public :: efp_potential_t
   public :: make_efp_potential
   public :: write_efp_potential
   !> Exposed so a reader can invert it. The printed contraction coefficients carry
   !> this factor, and recovering the raw ones is an exact division by it.
   public :: gamess_primitive_norm
   public :: from_gamess_ao_order
   !> Its inverse, for anything that changes an orbital and must hand the result back
   !> to a fragment still storing GAMESS's order. Rotation is the first such caller.
   public :: to_gamess_ao_order
   public :: frozen_core

   !> Longest line any section emits, with room to spare.
   integer, parameter :: MAX_LINE = 160

   !> Flattened extents of the Cartesian tensors the polarizability blocks carry:
   !> a pair of directions, a triple, and a quadruple. Nine is also the number of
   !> slots GAMESS writes a polarizability in, those being the same thing.
   integer, parameter :: N_CART_PAIR = 9
   integer, parameter :: N_CART_TRIPLE = 27
   integer, parameter :: N_CART_QUAD = 81

   !> Components in a Cartesian d and f shell.
   integer, parameter :: N_CART_D = 6
   integer, parameter :: N_CART_F = 10

   !> Row and column of each of GAMESS's nine polarizability slots. The
   !> off-diagonal triples are the transpose of what its labels suggest; that was
   !> measured in `validation/check_distributed_polarizability.py` rather than
   !> assumed, and it is the one convention here that a symmetric test tensor
   !> would not have caught.
   !>
   !> It is also the transpose of what `efinp.src:7552-7561` writes, and that is not a
   !> conflict: GAMESS indexes the tensor `(field, dipole)` where this code indexes it
   !> `(dipole, field)`, so both put the same number in slot 4 and the files agree. See
   !> the longer note on `POL_ROW` in `mqc_efp_read`, which carries the measurement.
   integer, parameter :: POL_ROW(N_CART_PAIR) = [1, 2, 3, 2, 3, 3, 1, 1, 2]
   integer, parameter :: POL_COL(N_CART_PAIR) = [1, 2, 3, 1, 1, 2, 2, 3, 3]

   !> libcint's index for each of GAMESS's six Cartesian d slots, and the
   !> normalization between the two codes' d functions. Both established in
   !> `validation/check_projection.py` against GAMESS's own coefficients.
   !> libcint's full-Cartesian quadrupole slots, which run xx,xy,xz,yx,...,zz.
   integer, parameter :: QXX = 1, QXY = 2, QXZ = 3, QYX = 4, QYY = 5
   integer, parameter :: QYZ = 6, QZX = 7, QZY = 8, QZZ = 9

   integer, parameter :: D_FROM_LIBCINT(N_CART_D) = [1, 4, 6, 2, 3, 5]
   real(dp), parameter :: D_NORMALIZATION = 1.585330892_dp

   !> The same for the ten Cartesian f slots, read off GAMESS's own coefficients in
   !> `validation/check_projection` and solved for in its Python half. Derived rather
   !> than reasoned about, exactly as the d map was.
   !>
   !> The molecule has to be in a frame with no zero coordinate for this to be
   !> solvable at all: planar water puts an exact zero in every function with an odd
   !> power of y, and a slot that is zero on both sides admits any scale factor. A
   !> first attempt in the planar frame returned a map that looked plausible, sent
   !> four different GAMESS slots to the same one of ours, and was wrong.
   integer, parameter :: F_FROM_LIBCINT(N_CART_F) = [1, 7, 10, 2, 3, 4, 8, 6, 9, 5]
   real(dp), parameter :: F_NORMALIZATION = 1.339849174_dp

   !> Which of the three normalization classes each GAMESS f slot belongs to. The
   !> measured scales come out as `F_NORMALIZATION` divided by one of 1, sqrt(5) and
   !> sqrt(15) -- the ratios are exact to eight figures, so they are written as the
   !> square roots rather than as the fitted decimals.
   integer, parameter :: F_CLASS(N_CART_F) = [1, 1, 1, 2, 2, 2, 2, 2, 2, 3]

   type :: efp_potential_t
      !! Every parameter a `.efp` carries that we can compute
      character(len=:), allocatable :: name        !! `$FRAGNAME`, without the `$`
      character(len=:), allocatable :: basis_name
      integer :: n_points = 0     !! Expansion points: atoms then bond midpoints
      integer :: n_atoms = 0
      integer :: nao = 0
      integer :: n_occ = 0        !! Including the core, which `CTFOK` needs
      integer :: n_lmo = 0        !! Valence localized orbitals
      integer :: multiplicity = 1
      real(dp) :: vdwscl = 0.7_dp  !! The screening grid's van der Waals scale
      character(len=8), allocatable :: labels(:)      !! `A01O`, `BO21`, ...
      real(dp), allocatable :: points(:, :)           !! (3, n_points), Bohr
      real(dp), allocatable :: mass(:)                !! amu, zero at a midpoint
      real(dp), allocatable :: charge(:)              !! Z, zero at a midpoint
      real(dp), allocatable :: q_elec(:), q_nuc(:)    !! MONOPOLES
      real(dp), allocatable :: dipole(:, :)           !! (3, n_points)
      real(dp), allocatable :: quadrupole(:, :)       !! (6, n_points)
      real(dp), allocatable :: octopole(:, :)         !! (10, n_points)
      real(dp), allocatable :: centroids(:, :)        !! (3, n_lmo)
      real(dp), allocatable :: static_pol(:, :, :)    !! (3, 3, n_lmo)
      real(dp), allocatable :: dynamic_pol(:, :, :, :)  !! (3, 3, n_lmo, n_freq)
      real(dp), allocatable :: dipquad_pre(:, :, :, :, :)
         !! The dipole-quadrupole tensor *before* the translation, which is what
         !! `QQSHIFT` takes as an input -- so it has to be kept, not just the
         !! shifted form the file carries.
      real(dp), allocatable :: quadquad(:, :, :, :, :, :)
         !! `(3, 3, 3, 3, n_lmo, n_freq)`, after the write-time translation.
      real(dp), allocatable :: dipquad(:, :, :, :, :)
         !! `(3, 3, 3, n_lmo, n_freq)` as `A'(a,b,c)`, **after** the write-time
         !! translation to each centroid. Post-shift because that is the form the
         !! file carries and the shift needs the dipole-dipole tensors, so keeping
         !! the pre-shift tensor would mean carrying its inputs too.
      real(dp), allocatable :: frequencies(:)         !! Imaginary, a.u.
      real(dp), allocatable :: fock_lmo(:, :)         !! (n_lmo, n_lmo)
      real(dp), allocatable :: orbitals(:, :)         !! LMOs in GAMESS's AO order
      real(dp), allocatable :: canonical(:, :)
         !! All the canonical MOs, in GAMESS's AO order. `CTVEC` in its
         !! canonical-orbital form is exactly this matrix.
      real(dp), allocatable :: eps_occ(:)             !! CTFOK
      real(dp), allocatable :: screen2(:)             !! Exponential alpha per point
      real(dp), allocatable :: screen(:)              !! Gaussian alpha per point
      character(len=MAX_LINE), allocatable :: basis_lines(:)
   contains
      procedure :: destroy => potential_destroy
   end type efp_potential_t

contains

   subroutine potential_destroy(self)
      class(efp_potential_t), intent(inout) :: self

      if (allocated(self%name)) deallocate (self%name)
      if (allocated(self%basis_name)) deallocate (self%basis_name)
      if (allocated(self%labels)) deallocate (self%labels)
      if (allocated(self%points)) deallocate (self%points)
      if (allocated(self%mass)) deallocate (self%mass)
      if (allocated(self%charge)) deallocate (self%charge)
      if (allocated(self%q_elec)) deallocate (self%q_elec)
      if (allocated(self%q_nuc)) deallocate (self%q_nuc)
      if (allocated(self%dipole)) deallocate (self%dipole)
      if (allocated(self%quadrupole)) deallocate (self%quadrupole)
      if (allocated(self%octopole)) deallocate (self%octopole)
      if (allocated(self%centroids)) deallocate (self%centroids)
      if (allocated(self%static_pol)) deallocate (self%static_pol)
      if (allocated(self%dynamic_pol)) deallocate (self%dynamic_pol)
      if (allocated(self%dipquad)) deallocate (self%dipquad)
      if (allocated(self%quadquad)) deallocate (self%quadquad)
      if (allocated(self%dipquad_pre)) deallocate (self%dipquad_pre)
      if (allocated(self%frequencies)) deallocate (self%frequencies)
      if (allocated(self%fock_lmo)) deallocate (self%fock_lmo)
      if (allocated(self%orbitals)) deallocate (self%orbitals)
      if (allocated(self%canonical)) deallocate (self%canonical)
      if (allocated(self%eps_occ)) deallocate (self%eps_occ)
      if (allocated(self%screen2)) deallocate (self%screen2)
      if (allocated(self%screen)) deallocate (self%screen)
      if (allocated(self%basis_lines)) deallocate (self%basis_lines)
      self%n_points = 0
      self%n_atoms = 0
      self%nao = 0
      self%n_occ = 0
      self%n_lmo = 0
   end subroutine potential_destroy

   subroutine make_efp_potential(atomic_numbers, element_symbols, coordinates, &
                                 basis_name, name, pot, error, charge, n_core, &
                                 vdwscl, verbose, aux_basis, guess)
      !! The whole pipeline: SCF, localization, and every parameter block
      !!
      !! The order is forced by what depends on what. The SCF gives the density
      !! the multipoles are taken from and the orbitals everything else needs;
      !! localization gives the centres the polarizabilities and the exchange
      !! repulsion data are expressed on; the multipoles have to exist before the
      !! screening can be fitted, because what screening fits is the error the
      !! *damped multipole* potential makes against the quantum one.
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)    !! (3, natm), Bohr
      character(len=*), intent(in) :: basis_name
      character(len=*), intent(in) :: name         !! Fragment name, e.g. `WATER`
      type(efp_potential_t), intent(out) :: pot
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: charge
         !! Net charge. Ignoring it is not a small error: the electron count comes out
         !! wrong, and for a cation it comes out odd, so a closed-shell reference is
         !! refused outright. Ionic-liquid fragments are the obvious case.
      integer, intent(in), optional :: n_core
         !! Orbitals excluded from the localized set. Default is the standard
         !! frozen core, which is what MAKEFP uses: its polarizable points and
         !! its exchange-repulsion orbitals are valence only.
      real(dp), intent(in), optional :: vdwscl
      logical, intent(in), optional :: verbose
      character(len=*), intent(in), optional :: guess
         !! Initial-guess name from the deck (`keywords.guess.type`). Default is
         !! "auto", which resolves to SAD on this backend -- the same guess the
         !! Energy driver uses. Without it MAKEFP ran the core guess: a system the
         !! Energy driver converged in ten iterations did not converge in two
         !! hundred here, and the guess was the whole of the difference.
      character(len=*), intent(in), optional :: aux_basis
         !! Fit the dynamic response rather than building its Hessian exactly. That
         !! build is `n_ov` Fock builds and is most of what a potential costs; fitted,
         !! it is two matrix products. The approximation is real and is measured by
         !! `validation/check_df_hessian`, so this is asked for rather than inferred.
         !!
         !! The auxiliary basis must match the orbital basis in angular form, which
         !! `build_df_mo_block` checks: libcint builds all three centres of a fitting
         !! integral in one form. In practice that means a Cartesian Pople orbital
         !! basis has no usable fitting set here, since the ones on hand are spherical.

      type(libcint_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(dma_result_t) :: dma
      type(response_hessian_t) :: shared_hessian
      type(screening_target_t) :: screen_target
      real(dp), allocatable :: loc(:, :), ovl(:, :), sc(:, :), w(:, :), scaled(:, :)
      real(dp), allocatable :: alpha(:)
      real(dp) :: rms_exp, rms_gauss
      integer :: natm, core, i, j, k, n_valence, n_electrons
      integer :: guess_kind
      real(dp), allocatable :: guess_total(:, :)
      character(len=:), allocatable :: guess_name
      type(timer_type) :: stage
      logical :: talk

      character(len=MAX_LINE_LENGTH) :: line

      talk = .false.
      if (present(verbose)) talk = verbose
      if (talk) call stage%start()
      natm = size(atomic_numbers)
      n_electrons = sum(atomic_numbers)
      if (present(charge)) n_electrons = n_electrons - charge
      pot%name = trim(name)
      pot%basis_name = trim(basis_name)
      pot%n_atoms = natm
      if (present(vdwscl)) pot%vdwscl = vdwscl

      if (size(coordinates, 1) /= 3 .or. size(coordinates, 2) /= natm) then
         call error%set(ERROR_VALIDATION, "makefp: coordinates must be (3, natm)")
         return
      end if

      call build_libcint_molecule(atomic_numbers, element_symbols, coordinates, &
                                  basis_name, mol, error)
      if (error%has_error()) return
      pot%nao = mol%nao

      if (talk) then
         write (line, "(A,A,A,I0,A)") "  basis ", trim(basis_name), ", ", mol%nao, " functions"
         call logger%info(trim(line))
      end if

      if (present(aux_basis)) then
         ! Read the fitting set in whatever angular form the orbital basis is in.
         ! libcint builds all three centres of a fitting integral in one form, so the
         ! two have to agree -- and the writer needs the orbital basis Cartesian,
         ! while every fitting set on hand is declared spherical. Without this the
         ! fitted path and a written potential are mutually exclusive.
         !
         ! Legitimate because an auxiliary basis is a fitting space, not a
         ! wavefunction: taking the Cartesian components of the same primitives
         ! enlarges that space rather than changing what is being fitted. It does
         ! make the space redundant -- the Cartesian d shells carry an s contaminant
         ! that duplicates the aux s functions -- so the metric is more nearly
         ! singular than it would otherwise be, and what that costs is measured in
         ! `validation/check_df_hessian` rather than assumed.
         call build_libcint_molecule(atomic_numbers, element_symbols, coordinates, &
                                     aux_basis, aux, error, &
                                     force_cartesian=mol%cartesian)
         if (error%has_error()) then
            call mol%destroy()
            return
         end if
         if (talk) then
            write (line, "(A,A,A,I0,A)") "  density fitting: ", trim(aux_basis), ", ", aux%nao, &
               " functions (the SCF and the dynamic response)"
            call logger%info(trim(line))
         end if
      else if (talk) then
         call logger%info("  density fitting: off (the SCF and the dynamic response are exact)")
      end if

      ! Checked here, not where the ordering map needs it. The map runs after the
      ! SCF, the localization and every response solve, so a basis this cannot emit
      ! would otherwise be refused several minutes in.
      call check_angular_form(mol, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if

      ! Honour the deck's guess, the same as the Energy path: MAKEFP passing no
      ! guess meant the core guess, which on anything larger than a few atoms does
      ! not converge in 200 iterations where SAD converges in ten. `guess_total`
      ! is left unallocated for core/gwh and the SCF reads it only for the atomic
      ! guesses, so passing it unconditionally is safe.
      guess_name = "auto"
      if (present(guess)) guess_name = guess
      call build_restricted_guess(mol, guess_name, guess_kind, guess_total, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      if (talk) call logger%info("  initial guess: "//guess_display_name(guess_kind))

      ! `talk` (the makefp verbosity, true for a real run from the driver) rather
      ! than a hardcoded .false.: the SCF that builds the potential is the main
      ! computation here, so its iteration table belongs alongside the stage
      ! timings this routine already prints, not hidden a level down.
      !
      ! 1e-10 energy / 1e-8 density: the density is what the multipoles and the
      ! response are taken from, and 1e-8 is tight enough that the potential does
      ! not move in any digit it reports. The energy keeps the 100:1 ratio the
      ! rest of the code uses (SCF_ENERGY_TOL / SCF_DENSITY_TOL), so the density is
      ! the binding criterion -- by the time it reaches 1e-8 the energy is already
      ! past 1e-12. The old 1e-12 / 1e-10 chased digits below what a fragment
      ! potential can carry and cost several iterations doing it.
      !
      ! `aux` present means density-fit the SCF too, not just the response: the
      ! `scf.density_fitting` flag now fits both stages against the one auxiliary
      ! basis, the same as the key does on the Energy path. Absent, the SCF is
      ! exact.
      if (present(aux_basis)) then
         call run_libcint_rhf(mol, n_electrons, 200, 1.0e-10_dp, 1.0e-8_dp, &
                              talk, scf, error, guess=guess_kind, guess_density=guess_total, aux=aux)
      else
         call run_libcint_rhf(mol, n_electrons, 200, 1.0e-10_dp, 1.0e-8_dp, &
                              talk, scf, error, guess=guess_kind, guess_density=guess_total)
      end if
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      if (.not. scf%converged) then
         call error%set(ERROR_VALIDATION, "makefp: the SCF did not converge")
         call mol%destroy()
         return
      end if
      pot%n_occ = scf%n_occupied
      if (talk) call report(stage, "SCF", talk)
      if (talk) then
         write (line, "(A,F18.10)") "  RHF energy ", scf%energy
         call logger%info(trim(line))
      end if

      core = frozen_core(atomic_numbers)
      if (present(n_core)) core = n_core
      n_valence = pot%n_occ - core
      if (n_valence < 1) then
         call error%set(ERROR_VALIDATION, "makefp: no valence orbitals to localize")
         call mol%destroy()
         return
      end if
      pot%n_lmo = n_valence

      allocate (pot%eps_occ(pot%n_occ))
      pot%eps_occ = scf%orbital_energies(1:pot%n_occ)

      ! --- localization: the centres everything below is expressed on ------------
      call boys_localize(mol, scf%orbitals(:, core + 1:pot%n_occ), n_valence, &
                         loc, pot%centroids, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      if (talk) then
         write (line, "(A,I0,A,I0,A)") "  localized ", n_valence, " of ", pot%n_occ, " occupied orbitals"
         call logger%info(trim(line))
      end if
      if (talk) call report(stage, "localization", talk)

      ! --- electrostatics -------------------------------------------------------
      call distributed_multipoles(mol, scf%density, atomic_numbers, dma, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      pot%n_points = size(dma%labels)
      allocate (pot%labels(pot%n_points), pot%points(3, pot%n_points))
      pot%labels = dma%labels
      pot%points = dma%points
      allocate (pot%q_elec(pot%n_points), pot%q_nuc(pot%n_points))
      pot%q_elec = dma%electronic
      pot%q_nuc = dma%nuclear
      allocate (pot%dipole(3, pot%n_points), pot%quadrupole(6, pot%n_points))
      allocate (pot%octopole(10, pot%n_points))
      pot%dipole = dma%dipole
      pot%quadrupole = dma%quadrupole
      pot%octopole = dma%octopole

      ! Mass and charge sit on atoms only; a bond midpoint carries neither, which
      ! is how GAMESS's reader tells the two kinds of point apart.
      allocate (pot%mass(pot%n_points), pot%charge(pot%n_points))
      pot%mass = 0.0_dp
      pot%charge = 0.0_dp
      do i = 1, natm
         pot%mass(i) = element_mass(atomic_numbers(i))
         pot%charge(i) = real(atomic_numbers(i), dp)
      end do

      ! --- polarization ---------------------------------------------------------
      if (talk) call report(stage, "distributed multipoles", talk)

      allocate (pot%frequencies(N_CASIMIR_POLDER))
      pot%frequencies = casimir_polder_frequencies()
      ! One Hessian for all three dynamic blocks. It depends on the reference alone,
      ! so rebuilding it per block is three identical builds -- 55 seconds each at 115
      ! orbitals, which was most of the potential's cost.
      ! Only this call can build the Hessian: it is the first of the three, and the
      ! two blocks after it are handed the built one and reuse it. So this is the one
      ! place the auxiliary basis has to reach, and an optional dummy cannot be passed
      ! conditionally from a local -- hence the branch rather than one call.
      !
      ! Ahead of the static block, which used to come first. The static response is
      ! `(A+B) U = -h`, the same operator these blocks build, so once it exists the
      ! static solve is a factorization rather than a conjugate-gradient run --
      ! forty seconds of the hundred a tripeptide took, for a result that agrees to
      ! thirteen digits.
      if (present(aux_basis)) then
         call dipole_quadrupole_block(mol, scf, coordinates, atomic_numbers, core, pot, &
                                      shared_hessian, error, progress=talk, aux=aux)
      else
         call dipole_quadrupole_block(mol, scf, coordinates, atomic_numbers, core, pot, &
                                      shared_hessian, error, progress=talk)
      end if
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      if (talk) then
         write (line, "(A,I0,A)") "  polarizabilities at ", N_CASIMIR_POLDER, " imaginary frequencies"
         call logger%info(trim(line))
         write (line, "(A)") "  dipole-quadrupole dispersion tensors"
         call logger%info(trim(line))
      end if
      if (talk) call report(stage, "all three dynamic blocks, with the Hessian build", talk)

      ! The static block, solved against the Hessian the dynamic blocks just
      ! built rather than by iterating for it again.
      call distributed_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                      pot%n_occ, pot%static_pol, pot%centroids, &
                                      error, n_core=core, hessian=shared_hessian)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      if (talk) call report(stage, "static polarizability", talk)

      ! --- exchange repulsion: the LMO Fock matrix and the orbitals themselves --
      ! F in the LMO basis is W^T diag(eps) W with W = C_occ^T S C_loc, so no AO
      ! Fock matrix is needed and nothing here depends on basis function ordering.
      call mol%overlap(ovl)
      allocate (sc(mol%nao, n_valence), w(pot%n_occ, n_valence))
      call pic_gemm(ovl, loc, sc)
      call pic_gemm(scf%orbitals(:, 1:pot%n_occ), sc, w, transa="T")
      allocate (scaled(pot%n_occ, n_valence), pot%fock_lmo(n_valence, n_valence))
      do j = 1, n_valence
         do k = 1, pot%n_occ
            scaled(k, j) = scf%orbital_energies(k)*w(k, j)
         end do
      end do
      call pic_gemm(w, scaled, pot%fock_lmo, transa="T")

      call to_gamess_ao_order(mol, loc, pot%orbitals, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if

      ! The charge-transfer basis. `CTVEC` has two forms in GAMESS and this is the
      ! second: with `$MAKEFP CTVVO=.FALSE.` it writes the whole canonical MO
      ! matrix and the header `CTVEC NA NUM` -- the branch GAMESS labels CMO -- where
      ! its default path writes `NOCC` occupied orbitals plus a set of
      ! quasi-atomic valence virtuals built by `VVOS`. The canonical form needs no
      ! extra machinery and is what GAMESS itself recommends when the valence
      ! virtuals cannot be formed -- "PLEASE TRY IT AGAIN WITH CANONICAL ORBITALS".
      call to_gamess_ao_order(mol, scf%orbitals, pot%canonical, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if

      call projection_basis_lines(atomic_numbers, pot%labels, pot%points, &
                                  element_symbols, basis_name, pot%basis_lines, &
                                  error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if

      ! --- charge penetration screening ----------------------------------------
      ! Last, because it is fitted to the error the multipoles above make.
      !
      ! One grid and one quantum potential for both damping forms. They are fitted
      ! to the same target and differ only in the damping term of the objective, so
      ! the first fit hands its target to the second rather than the second building
      ! an identical one. Building it is nearly the whole cost: the grid runs to
      ! tens of thousands of points and each needs an integral over every shell pair.
      call fit_screening(mol, scf%density, dma, atomic_numbers, SCREEN_EXPONENTIAL, &
                         alpha, error, target=screen_target, residual=rms_exp)
      if (error%has_error()) then
         call screen_target%destroy()
         call mol%destroy()
         return
      end if
      allocate (pot%screen2(pot%n_points))
      pot%screen2 = alpha
      deallocate (alpha)
      call fit_screening(mol, scf%density, dma, atomic_numbers, SCREEN_GAUSSIAN, &
                         alpha, error, target=screen_target, residual=rms_gauss)
      call screen_target%destroy()
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      allocate (pot%screen(pot%n_points))
      pot%screen = alpha
      deallocate (alpha)
      if (talk) then
         write (line, "(A,F0.4,A,F0.4,A)") "  screening fitted: exponential misses by ", &
            rms_exp, " kcal/mol, Gaussian by ", rms_gauss, " kcal/mol"
         call logger%info(trim(line))
      end if
      if (talk) call report(stage, "charge-penetration screening", talk)

      call shared_hessian%destroy()
      call mol%destroy()
      deallocate (loc, ovl, sc, w, scaled)
   end subroutine make_efp_potential

   subroutine dipole_quadrupole_block(mol, scf, coordinates, atomic_numbers, core, &
                                      pot, hessian, error, progress, aux)
      !! `DIPOLE-QUADRUPOLE DYNAMIC POLARIZABLE POINTS`, ready to write
      !!
      !! Three conventions here were established by
      !! `validation/check_dipquad_sumrule`, which pins them by structure rather
      !! than by fitting -- see its docstring and the record in
      !! `mqc_libcint_cphf`. Getting any one of them wrong leaves a tensor that
      !! passes every internal check and disagrees with GAMESS, which is what
      !! happened for a long time.
      !!
      !!   * **The quadrupole measures and the dipole drives**, both expanded about
      !!     the centre of mass, and the quadrupole is the traceless Buckingham
      !!     form. Per orbital that is not the same as the reverse assignment,
      !!     because the projector onto the localized set does not commute with the
      !!     response operator; summed over orbitals they agree.
      !!   * **The write-time translation to each centroid**, `DQSHIFT`, whose
      !!     `delta_bc` term takes the dipole-dipole tensor **transposed**:
      !!     `alpha(a,d)`, not `alpha(d,a)` as the rest of the formula reads. That
      !!     transpose is the whole of what used to be a 14% discrepancy; with it
      !!     the block agrees to 8.9e-05, which is GAMESS's own precision.
      !!   * **The dipole-dipole tensor in the shift is the dynamic one at the same
      !!     frequency**, not the static one.
      type(libcint_molecule_t), intent(in) :: mol
      type(rhf_result_t), intent(in) :: scf
      real(dp), intent(in) :: coordinates(:, :)
      integer, intent(in) :: atomic_numbers(:)
      integer, intent(in) :: core
      type(efp_potential_t), intent(inout) :: pot
      type(response_hessian_t), intent(inout) :: hessian
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: progress
      type(libcint_molecule_t), intent(in), optional :: aux
         !! Fit the Hessian rather than build it exactly; passed straight through.

      real(dp), allocatable :: dip(:, :, :), quad(:, :, :), buck(:, :, :)
      real(dp), allocatable :: both(:, :, :), all_blocks(:, :, :, :)
      real(dp), allocatable :: raw(:, :, :, :), centroids(:, :), qq(:, :, :, :)
      real(dp) :: com(3), r(3), alpha(3, 3)
      real(dp) :: mass_total, isotropic
      integer :: i, a, b, c, d, k, f, n_freq, n_both

      com = 0.0_dp
      mass_total = 0.0_dp
      do i = 1, size(atomic_numbers)
         com = com + element_mass(atomic_numbers(i))*coordinates(:, i)
         mass_total = mass_total + element_mass(atomic_numbers(i))
      end do
      com = com/mass_total

      call multipole_matrices(mol, com, 1, dip, error)
      if (error%has_error()) return
      call multipole_matrices(mol, com, 2, quad, error)
      if (error%has_error()) return

      ! The traceless Buckingham quadrupole, as GAMESS builds it
      ! -- the traceless Buckingham form -- kept as all nine Cartesian slots so that no expansion
      ! of six unique values into nine has to be guessed at.
      allocate (buck(mol%nao, mol%nao, 9))
      buck(:, :, QXX) = 0.5_dp*(2.0_dp*quad(:, :, QXX) - quad(:, :, QYY) &
                                - quad(:, :, QZZ))
      buck(:, :, QYY) = 0.5_dp*(2.0_dp*quad(:, :, QYY) - quad(:, :, QXX) &
                                - quad(:, :, QZZ))
      buck(:, :, QZZ) = 0.5_dp*(2.0_dp*quad(:, :, QZZ) - quad(:, :, QXX) &
                                - quad(:, :, QYY))
      buck(:, :, QXY) = 1.5_dp*quad(:, :, QXY)
      buck(:, :, QXZ) = 1.5_dp*quad(:, :, QXZ)
      buck(:, :, QYZ) = 1.5_dp*quad(:, :, QYZ)
      buck(:, :, QYX) = buck(:, :, QXY)
      buck(:, :, QZX) = buck(:, :, QXZ)
      buck(:, :, QZY) = buck(:, :, QYZ)

      ! All three dynamic blocks from one solve.
      !
      ! They are three contractions of one response, not three responses: the
      ! dipole-dipole block measures and drives with `dip`, the mixed one measures
      ! with `buck` and drives with `dip`, and the quadrupole one uses `buck` both
      ! sides. So the *driving* operators are `dip` and `buck` together -- twelve of
      ! them -- and every block is a slice of the twelve-by-twelve result.
      !
      ! Run as three calls this cost three sets of twelve factorizations of the same
      ! shifted matrix, differing only in which right-hand sides came along. Two
      ! thirds of that was redundant, and at 8350 pairs a factorization is eleven
      ! seconds on one core.
      !
      ! The dipole-dipole slice is what `distributed_dynamic_polarizability` used to
      ! return, which is the same projection with `measure` and `respond` both
      ! `dip` -- so it is computed here now and the separate call is gone.
      n_both = 3 + size(buck, 3)
      allocate (both(mol%nao, mol%nao, n_both))
      both(:, :, 1:3) = dip
      both(:, :, 4:n_both) = buck

      if (present(aux)) then
         call distributed_dynamic_cross(mol, scf%orbitals, scf%orbital_energies, &
                                        pot%n_occ, pot%frequencies, both, both, &
                                        all_blocks, centroids, error, n_core=core, &
                                        hessian=hessian, progress=progress, aux=aux)
      else
         call distributed_dynamic_cross(mol, scf%orbitals, scf%orbital_energies, &
                                        pot%n_occ, pot%frequencies, both, both, &
                                        all_blocks, centroids, error, n_core=core, &
                                        hessian=hessian, progress=progress)
      end if
      if (error%has_error()) return
      deallocate (both)

      n_freq = size(pot%frequencies)
      if (allocated(pot%centroids)) deallocate (pot%centroids)
      allocate (pot%centroids, source=centroids)
      allocate (pot%dynamic_pol(3, 3, pot%n_lmo, n_freq))
      pot%dynamic_pol = all_blocks(1:3, 1:3, :, :)
      allocate (raw(size(buck, 3), 3, pot%n_lmo, n_freq))
      raw = all_blocks(4:n_both, 1:3, :, :)
      allocate (pot%dipquad(3, 3, 3, pot%n_lmo, n_freq))
      allocate (pot%dipquad_pre(3, 3, 3, pot%n_lmo, n_freq))
      do f = 1, n_freq
         do k = 1, pot%n_lmo
            r = pot%centroids(:, k) - com
            alpha = pot%dynamic_pol(:, :, k, f)
            do a = 1, 3
               ! The delta_bc term's transpose: alpha(a,d), not alpha(d,a).
               isotropic = 0.0_dp
               do d = 1, 3
                  isotropic = isotropic + r(d)*alpha(a, d)
               end do
               do b = 1, 3
                  do c = 1, 3
                     ! raw is (quadrupole slot, dipole, orbital, frequency); the
                     ! nine quadrupole slots run with the second index fastest.
                     pot%dipquad_pre(a, b, c, k, f) = raw((b - 1)*3 + c, a, k, f)
                     pot%dipquad(a, b, c, k, f) = raw((b - 1)*3 + c, a, k, f) &
                                                  - 1.5_dp*(r(b)*alpha(c, a) &
                                                            + r(c)*alpha(a, b))
                     if (b == c) then
                        pot%dipquad(a, b, c, k, f) = pot%dipquad(a, b, c, k, f) &
                                                     + isotropic
                     end if
                  end do
               end do
            end do
         end do
      end do

      ! --- the quadrupole-quadrupole block ------------------------------------
      ! Same operator on both sides, and the factor is 1/3 rather than 1: `LQQPOL`
      ! carries a factor of a third that the dipole-quadrupole one does not. Confirmed
      ! against GAMESS's own molecular `QUAD-QUAD POLARIZABILITY`, which
      ! `$MAKEFP MOLPOL=.TRUE.` writes with no translation applied: our response
      ! summed over the orbitals and divided by three reproduces it to 1.5e-05.
      allocate (qq(size(buck, 3), size(buck, 3), pot%n_lmo, n_freq))
      qq = all_blocks(4:n_both, 4:n_both, :, :)
      deallocate (all_blocks)

      allocate (pot%quadquad(3, 3, 3, 3, pot%n_lmo, n_freq))
      do f = 1, n_freq
         do k = 1, pot%n_lmo
            r = pot%centroids(:, k) - com
            alpha = pot%dynamic_pol(:, :, k, f)
            call qq_shift(qq(:, :, k, f)/3.0_dp, pot%dipquad_pre(:, :, :, k, f), &
                          alpha, r, pot%quadquad(:, :, :, :, k, f))
         end do
      end do

      deallocate (dip, quad, buck, raw, qq, centroids)
   end subroutine dipole_quadrupole_block

   subroutine qq_shift(qq, dq, alpha, r, shifted)
      !! `QQSHIFT`, the quadrupole-quadrupole translation
      !!
      !! Transcribed term for term. It mixes in both the dipole-dipole and the
      !! *pre-shift* dipole-quadrupole tensors, which is why the two blocks are
      !! built together.
      real(dp), intent(in) :: qq(N_CART_PAIR, N_CART_PAIR)       !! pre-shift, already scaled by 1/3
      real(dp), intent(in) :: dq(3, 3, 3)    !! pre-shift dipole-quadrupole
      real(dp), intent(in) :: alpha(3, 3)
      real(dp), intent(in) :: r(3)
      real(dp), intent(out) :: shifted(3, 3, 3, 3)

      real(dp) :: a1, a2, rralph, adelt1, adelt2, rrad1, rrad2, rrdalph
      integer :: a, b, c, e, i, j

      do a = 1, 3
         do b = 1, 3
            do c = 1, 3
               do e = 1, 3
                  a1 = r(a)*dq(b, c, e) + r(b)*dq(a, c, e)
                  a2 = r(c)*dq(e, a, b) + r(e)*dq(c, a, b)
                  rralph = r(a)*r(c)*alpha(b, e) + r(a)*r(e)*alpha(b, c) &
                           + r(b)*r(c)*alpha(a, e) + r(b)*r(e)*alpha(a, c)
                  adelt1 = 0.0_dp
                  adelt2 = 0.0_dp
                  rrad1 = 0.0_dp
                  rrad2 = 0.0_dp
                  rrdalph = 0.0_dp
                  do i = 1, 3
                     if (a == b) adelt1 = adelt1 + r(i)*dq(i, c, e)
                     if (c == e) adelt2 = adelt2 + r(i)*dq(i, a, b)
                     if (a == b) rrad1 = rrad1 + r(c)*r(i)*alpha(i, e) &
                                         + r(e)*r(i)*alpha(i, c)
                     if (c == e) rrad2 = rrad2 + r(a)*r(i)*alpha(b, i) &
                                         + r(b)*r(i)*alpha(a, i)
                     if (a == b .and. c == e) then
                        do j = 1, 3
                           rrdalph = rrdalph + r(i)*r(j)*alpha(i, j)
                        end do
                     end if
                  end do
                  shifted(a, b, c, e) = qq((a - 1)*3 + b, (c - 1)*3 + e) &
                                        - 0.5_dp*(a1 + a2 + rrad1 + rrad2) &
                                        + (adelt1 + adelt2 + rrdalph)/3.0_dp &
                                        + 0.75_dp*rralph
               end do
            end do
         end do
      end do
   end subroutine qq_shift

   pure function frozen_core(atomic_numbers) result(n)
      !! The standard frozen core, which is the set MAKEFP excludes
      integer, intent(in) :: atomic_numbers(:)
      integer :: n
      integer :: i, z

      n = 0
      do i = 1, size(atomic_numbers)
         z = atomic_numbers(i)
         if (z > 2 .and. z <= 10) then
            n = n + 1
         else if (z > 10 .and. z <= 18) then
            n = n + 5
         else if (z > 18 .and. z <= 36) then
            n = n + 9
         else if (z > 36) then
            n = n + 18
         end if
      end do
   end function frozen_core

   subroutine report(stage, what, talk)
      !! Seconds for the stage just finished, then restart the clock
      !!
      !! Printed by the emitter rather than only by a test harness: a fragment of any
      !! size takes long enough that "where did that go" is a question the program
      !! should answer without being rebuilt. At 115 orbitals the dynamic response is
      !! 88% of the run and everything else is under four seconds, but that split
      !! moves with the fragment -- the localization is quadratic in the occupied
      !! count and the multipoles work over primitive pairs, so neither can be
      !! assumed small for something larger.
      !!
      !! `pic_timer` rather than `system_clock` directly: it uses the OpenMP clock
      !! when PIC is built with OpenMP, which is what a threaded stage wants, and it
      !! is what the rest of this codebase already times with.
      type(timer_type), intent(inout) :: stage
      character(len=*), intent(in) :: what
      logical, intent(in) :: talk

      character(len=MAX_LINE_LENGTH) :: line

      if (.not. talk) return
      call stage%stop()
      write (line, "(A,F9.1,A,A)") "      ", stage%get_elapsed_time(), " s  ", what
      call logger%info(trim(line))
      flush (6)
      call stage%start()
   end subroutine report

   subroutine check_angular_form(mol, error)
      !! Refuse a basis whose angular form this cannot write
      !!
      !! **The ordering map is Cartesian s, p and d, and both codes have to agree on
      !! that.** GAMESS defaults ISPHER = -1, Cartesian, which is what a Pople set
      !! wants -- the Basis Set Exchange declares 6-31G*'s d Cartesian too, so the
      !! validated path lines up without either side being asked.
      !!
      !! A Dunning or def2 set is not simply unsupported: BSE declares those
      !! spherical and GAMESS will not run them Cartesian, so both sides would agree
      !! on spherical. What is missing is the spherical ordering map and `ISPHER=1`
      !! in the deck, not the possibility. Emitting 58 orbitals where the reader
      !! expands 65 would be rejected or read as nonsense, so it is refused.
      type(libcint_molecule_t), intent(in) :: mol
      type(error_t), intent(inout) :: error

      integer :: ish
      logical :: has_high_l

      ! Only an issue if there is actually a shell it applies to: s and p are the
      ! same either way, so 6-31G* on hydrogen alone is not a spherical basis in any
      ! sense that matters, and refusing it was wrong.
      has_high_l = .false.
      do ish = 1, mol%nbas
         if (mol%bas(LIBCINT_ANG_OF, ish) >= 2) has_high_l = .true.
      end do
      if (.not. mol%cartesian .and. has_high_l) then
         call error%set(ERROR_VALIDATION, &
                        "makefp: this basis is spherical and only Cartesian s, p "// &
                        "and d are mapped to GAMESS's ordering. GAMESS reads a "// &
                        "spherical potential with ISPHER=1, so this needs a "// &
                        "spherical ordering map rather than a Cartesian basis.")
      end if
   end subroutine check_angular_form

   subroutine from_gamess_ao_order(mol, mapped, coefficients, error)
      !! The inverse of `to_gamess_ao_order`
      !!
      !! A potential stores its orbitals in GAMESS's AO order, so anything that reads
      !! one back and wants to use it with our own integrals has to undo the
      !! permutation and the normalization. Written here, beside the forward map, so
      !! the two cannot drift: a permutation inverted in another file is a permutation
      !! that gets out of step the next time this one changes.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: mapped(:, :)
      real(dp), allocatable, intent(out) :: coefficients(:, :)
      type(error_t), intent(inout) :: error

      integer :: ish, l, off, dim, slot
      real(dp) :: scale

      allocate (coefficients(size(mapped, 1), size(mapped, 2)))
      coefficients = mapped
      do ish = 1, mol%nbas
         l = mol%bas(LIBCINT_ANG_OF, ish)
         off = mol%shell_offset(ish)
         dim = shell_dim(mol%cartesian, ish - 1, mol%bas)
         if (l == 2 .and. dim == 6) then
            do slot = 1, 6
               scale = D_NORMALIZATION
               if (slot > 3) scale = D_NORMALIZATION/sqrt(3.0_dp)
               coefficients(off + D_FROM_LIBCINT(slot), :) = mapped(off + slot, :)/scale
            end do
         else if (l == 3 .and. dim == 10) then
            do slot = 1, 10
               select case (F_CLASS(slot))
               case (1)
                  scale = F_NORMALIZATION
               case (2)
                  scale = F_NORMALIZATION/sqrt(5.0_dp)
               case default
                  scale = F_NORMALIZATION/sqrt(15.0_dp)
               end select
               coefficients(off + F_FROM_LIBCINT(slot), :) = mapped(off + slot, :)/scale
            end do
         else if (l >= 2) then
            call error%set(ERROR_VALIDATION, "efp: this basis has g functions or "// &
                           "higher, whose ordering against GAMESS is not mapped")
            return
         end if
      end do
   end subroutine from_gamess_ao_order

   subroutine to_gamess_ao_order(mol, coefficients, mapped, error)
      !! Orbital coefficients in the AO order and normalization GAMESS reads
      !!
      !! Only the Cartesian d shells move. Our s and p already agree, including
      !! the interleaving inside a shared-exponent `L` shell, which is why the
      !! projection data was accepted as-is once the d block was fixed.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coefficients(:, :)
      real(dp), allocatable, intent(out) :: mapped(:, :)
      type(error_t), intent(inout) :: error

      integer :: ish, l, off, dim, slot
      real(dp) :: scale

      allocate (mapped(size(coefficients, 1), size(coefficients, 2)))
      mapped = coefficients
      do ish = 1, mol%nbas
         l = mol%bas(LIBCINT_ANG_OF, ish)
         off = mol%shell_offset(ish)
         dim = shell_dim(mol%cartesian, ish - 1, mol%bas)
         if (l == 2 .and. dim == 6) then
            do slot = 1, 6
               scale = D_NORMALIZATION
               if (slot > 3) scale = D_NORMALIZATION/sqrt(3.0_dp)
               mapped(off + slot, :) = coefficients(off + D_FROM_LIBCINT(slot), :)*scale
            end do
         else if (l == 3 .and. dim == 10) then
            do slot = 1, 10
               select case (F_CLASS(slot))
               case (1)
                  scale = F_NORMALIZATION
               case (2)
                  scale = F_NORMALIZATION/sqrt(5.0_dp)
               case default
                  scale = F_NORMALIZATION/sqrt(15.0_dp)
               end select
               mapped(off + slot, :) = coefficients(off + F_FROM_LIBCINT(slot), :)*scale
            end do
         else if (l >= 2) then
            ! g and up need their own permutation and normalizations against
            ! GAMESS's ordering, derived the way the d and f ones were: read off its
            ! own printed coefficients for a basis that has them. Fifteen Cartesian
            ! components with several normalization classes, so guessing is worse
            ! than usual here.
            call error%set(ERROR_VALIDATION, &
                           "makefp: this basis has g functions or higher, and only "// &
                           "Cartesian s, p, d and f are mapped to GAMESS's ordering "// &
                           "so far")
            return
         end if
      end do
   end subroutine to_gamess_ao_order

   subroutine projection_basis_lines(atomic_numbers, labels, points, symbols, &
                                     basis_name, lines, error)
      !! `PROJECTION BASIS SET`, in GAMESS's columns and its normalization
      !!
      !! Two things here are GAMESS's conventions rather than ours. Shells are
      !! named the way it names them, including `L` for a shared-exponent sp pair,
      !! because that is what its reader expects -- our own reader emits one shell
      !! per coefficient column, so an `L` has to be recognised by finding an s
      !! and a p over identical exponents. And the printed coefficient has the
      !! primitive normalization folded in, which is the factor
      !! `gamess_primitive_norm` supplies.
      !!
      !! The primitive counter runs across the whole file rather than restarting
      !! per atom.
      integer, intent(in) :: atomic_numbers(:)
      character(len=8), intent(in) :: labels(:)
      real(dp), intent(in) :: points(:, :)
      character(len=*), intent(in) :: symbols(:)
      character(len=*), intent(in) :: basis_name
      character(len=MAX_LINE), allocatable, intent(out) :: lines(:)
      type(error_t), intent(inout) :: error

      type(molecular_basis_type) :: basis
      character(len=:), allocatable :: path
      character(len=MAX_LINE), allocatable :: buffer(:)
      character(len=MAX_LINE) :: text
      integer :: n, natm, iatom, ish, nsh, k, primitive, valence, ncol, col
      integer :: l_of(2)
      logical :: is_l
      real(dp) :: expo, coefficient

      natm = size(atomic_numbers)
      call find_basis_file(basis_name, path, error)
      if (error%has_error()) return
      call build_molecular_basis_json(path, symbols, basis, error)
      if (error%has_error()) return

      allocate (buffer(4096))
      n = 0
      primitive = 0
      do iatom = 1, natm
         valence = atomic_numbers(iatom) - frozen_core([atomic_numbers(iatom)])*2
         ! Ten wide here, where COORDINATES uses eight -- what MAKEFP writes, and
         ! the two sections genuinely differ. Written as A8 then two blanks
         ! rather than A10, because A10 on an eight-character string pads on the
         ! left and would indent the label instead of the number.
         write (text, "(A8,A2,3F15.10,F7.1)") labels(iatom), "  ", &
            points(:, iatom), real(valence, dp)
         n = n + 1
         buffer(n) = text

         nsh = basis%elements(iatom)%nshells
         ish = 1
         do while (ish <= nsh)
            ! An sp pair arrives as two consecutive shells over the same
            ! exponents; GAMESS wants them as one "L" with two coefficient
            ! columns.
            is_l = .false.
            if (ish < nsh) then
               if (basis%elements(iatom)%shells(ish)%ang_mom == 0 .and. &
                   basis%elements(iatom)%shells(ish + 1)%ang_mom == 1) then
                  if (size(basis%elements(iatom)%shells(ish)%exponents) == &
                      size(basis%elements(iatom)%shells(ish + 1)%exponents)) then
                     is_l = maxval(abs(basis%elements(iatom)%shells(ish)%exponents - &
                                       basis%elements(iatom)%shells(ish + 1)%exponents)) &
                            < 1.0e-12_dp
                  end if
               end if
            end if

            if (is_l) then
               ncol = 2
               l_of = [0, 1]
               write (text, "(A,I11)") "   L", &
                  size(basis%elements(iatom)%shells(ish)%exponents)
            else
               ncol = 1
               l_of(1) = basis%elements(iatom)%shells(ish)%ang_mom
               write (text, "(A,A,I11)") "   ", &
                  shell_letter(l_of(1)), size(basis%elements(iatom)%shells(ish)%exponents)
            end if
            n = n + 1
            buffer(n) = text

            do k = 1, size(basis%elements(iatom)%shells(ish)%exponents)
               primitive = primitive + 1
               expo = basis%elements(iatom)%shells(ish)%exponents(k)
               write (text, "(I6,F21.10)") primitive, expo
               do col = 1, ncol
                  coefficient = basis%elements(iatom)%shells(ish + col - 1) &
                                %coefficients(k)
                  write (text(len_trim(text) + 1:), "(F15.8)") &
                     coefficient*gamess_primitive_norm(l_of(col), expo)
               end do
               n = n + 1
               buffer(n) = text
            end do
            ish = ish + ncol
         end do
         n = n + 1
         buffer(n) = "  "
      end do

      allocate (lines(n))
      lines = buffer(1:n)
      deallocate (buffer)
      call basis%destroy()
   end subroutine projection_basis_lines

   pure function shell_letter(l) result(c)
      integer, intent(in) :: l
      character(len=1) :: c
      character(len=*), parameter :: LETTERS = "SPDFGH"

      if (l >= 0 .and. l < len(LETTERS)) then
         c = LETTERS(l + 1:l + 1)
      else
         c = "?"
      end if
   end function shell_letter

   pure function gamess_primitive_norm(l, exponent) result(factor)
      !! The normalization GAMESS folds into a printed contraction coefficient
      integer, intent(in) :: l
      real(dp), intent(in) :: exponent
      real(dp) :: factor
      real(dp), parameter :: PI = 3.141592653589793_dp
      real(dp) :: double_factorial
      integer :: n

      double_factorial = 1.0_dp
      n = 2*l - 1
      do while (n > 1)
         double_factorial = double_factorial*real(n, dp)
         n = n - 2
      end do
      factor = (2.0_dp*exponent/PI)**0.75_dp*(4.0_dp*exponent)**(0.5_dp*real(l, dp)) &
               /sqrt(double_factorial)
   end function gamess_primitive_norm

   subroutine write_efp_potential(pot, path, error, omitted)
      !! The potential as a `.efp` file
      !!
      !! `omitted` comes back with the sections a complete MAKEFP would have
      !! written and this did not, so a caller can say so rather than let their
      !! absence pass unremarked.
      type(efp_potential_t), intent(in) :: pot
      character(len=*), intent(in) :: path
      type(error_t), intent(inout) :: error
      character(len=:), allocatable, intent(out), optional :: omitted

      integer :: unit, i, k, f, stat, a, b, c, e
      character(len=8) :: label
      real(dp) :: tensor(N_CART_PAIR)
      real(dp) :: wide(N_CART_TRIPLE)
      real(dp) :: broad(N_CART_QUAD)

      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call error%set(ERROR_VALIDATION, "makefp: cannot write "//trim(path))
         return
      end if

      write (unit, "(A)") &
         "          RUNTYP=MAKEFP EFFECTIVE FRAGMENT POTENTIAL DATA FOLLOWS..."
      write (unit, "(A)") "          "//pot%name//" GENERATED BY METALQUICHA"
      write (unit, "(A)") " $"//pot%name
      write (unit, "(A)") "EFP DATA FOR "//pot%name//" SCFTYP=RHF     "// &
         "... GENERATED WITH BASIS SET="//pot%basis_name

      write (unit, "(A)") " COORDINATES (BOHR)"
      do i = 1, pot%n_points
         write (unit, "(A8,3F15.10,F12.7,F5.1)") pot%labels(i), pot%points(:, i), &
            pot%mass(i), pot%charge(i)
      end do
      write (unit, "(A)") " STOP"

      write (unit, "(A)") " MONOPOLES "
      do i = 1, pot%n_points
         write (unit, "(A8,F15.10,F10.5)") pot%labels(i), pot%q_elec(i), pot%q_nuc(i)
      end do
      write (unit, "(A)") " STOP"

      write (unit, "(A)") " DIPOLES "
      do i = 1, pot%n_points
         call write_record(unit, pot%labels(i), pot%dipole(:, i), 16, 10, 3)
      end do
      write (unit, "(A)") " STOP"

      write (unit, "(A)") " QUADRUPOLES "
      do i = 1, pot%n_points
         call write_record(unit, pot%labels(i), pot%quadrupole(:, i), 16, 10, 4)
      end do
      write (unit, "(A)") " STOP"

      write (unit, "(A)") " OCTUPOLES  "
      do i = 1, pot%n_points
         call write_record(unit, pot%labels(i), pot%octopole(:, i), 17, 9, 4)
      end do
      write (unit, "(A)") " STOP"

      write (unit, "(A)") " POLARIZABLE POINTS"
      do k = 1, pot%n_lmo
         write (label, "(A,I0)") "CT", k
         do i = 1, 9
            tensor(i) = pot%static_pol(POL_ROW(i), POL_COL(i), k)
         end do
         call write_tensor_point(unit, label, pot%centroids(:, k), tensor)
      end do
      write (unit, "(A)") " STOP"

      ! The dynamic blocks run frequency-outermost, each stamped with its own
      ! frequency, and the label carries a gap -- "CT  2" where every other
      ! section writes "CT2".
      write (unit, "(A)") " DYNAMIC POLARIZABLE POINTS"
      do f = 1, size(pot%frequencies)
         do k = 1, pot%n_lmo
            write (label, "(A,I3)") "CT", k
            do i = 1, 9
               tensor(i) = pot%dynamic_pol(POL_ROW(i), POL_COL(i), k, f)
            end do
            ! Only the first point of a block carries the frequency, which is how
            ! MAKEFP writes it: the stamp opens a block rather than labelling a
            ! point. GAMESS's reader tolerates one per point, but a parser keying
            ! on the stamp -- ours included -- then sees every point as a new
            ! block.
            if (k == 1) then
               call write_tensor_point(unit, label, pot%centroids(:, k), tensor, &
                                       frequency=pot%frequencies(f))
            else
               call write_tensor_point(unit, label, pot%centroids(:, k), tensor)
            end if
         end do
      end do
      write (unit, "(A)") " STOP"

      ! The dipole-quadrupole block, 27 values a point. The slot order is
      ! `(a-1)*9 + (c-1)*3 + b` -- the *first* quadrupole index runs fastest, which
      ! is transposed from how the `DQSHIFT` source reads, and was pinned by
      ! requiring the pre-shift tensor's symmetry in `bc` to come back.
      write (unit, "(A)") " DIPOLE-QUADRUPOLE DYNAMIC POLARIZABLE POINTS"
      do f = 1, size(pot%frequencies)
         do k = 1, pot%n_lmo
            write (label, "(A,I3)") "CT", k
            do a = 1, 3
               do b = 1, 3
                  do c = 1, 3
                     wide((a - 1)*9 + (c - 1)*3 + b) = pot%dipquad(a, b, c, k, f)
                  end do
               end do
            end do
            if (k == 1) then
               write (unit, "(A,3F15.10,A,F9.6,A)") trim(label), &
                  pot%centroids(:, k), " -- FOR W=", pot%frequencies(f), "I A.U."
            else
               write (unit, "(A,3F15.10)") trim(label), pot%centroids(:, k)
            end if
            call write_values(unit, wide, 16, 10, 4)
         end do
      end do
      write (unit, "(A)") " STOP"

      ! The quadrupole-quadrupole block, 81 values a point, written with the last
      ! index fastest, the last of the four varying first. No
      ! transposition here, unlike the dipole-quadrupole slots: every `QQSHIFT`
      ! term is symmetric within each index pair, so the written values are too.
      write (unit, "(A)") " LMOQQPOL DYNAMIC POLARIZABLE POINTS"
      do f = 1, size(pot%frequencies)
         do k = 1, pot%n_lmo
            write (label, "(A,I3)") "CT", k
            i = 0
            do a = 1, 3
               do b = 1, 3
                  do c = 1, 3
                     do e = 1, 3
                        i = i + 1
                        broad(i) = pot%quadquad(a, b, c, e, k, f)
                     end do
                  end do
               end do
            end do
            if (k == 1) then
               write (unit, "(A,3F15.10,A,F9.6,A)") trim(label), &
                  pot%centroids(:, k), " -- FOR W=", pot%frequencies(f), "I A.U."
            else
               write (unit, "(A,3F15.10)") trim(label), pot%centroids(:, k)
            end if
            call write_values(unit, broad, 16, 10, 4)
         end do
      end do
      write (unit, "(A)") " STOP"

      write (unit, "(A)") " PROJECTION BASIS SET"
      do i = 1, size(pot%basis_lines)
         write (unit, "(A)") trim(pot%basis_lines(i))
      end do
      write (unit, "(A)") " STOP"

      write (unit, "(A,I5)") " MULTIPLICITY", pot%multiplicity
      write (unit, "(A)") " STOP"

      write (unit, "(A,2I7)") " PROJECTION WAVEFUNCTION", pot%n_lmo, pot%nao
      call write_wavefunction(unit, pot%orbitals)

      write (unit, "(A)") " FOCK MATRIX ELEMENTS"
      call write_lower_triangle(unit, pot%fock_lmo)

      write (unit, "(A)") " LMO CENTROIDS"
      do k = 1, pot%n_lmo
         write (label, "(A,I0)") "CT", k
         write (unit, "(A3,3F15.10)") adjustl(label), pot%centroids(:, k)
      end do
      write (unit, "(A)") " STOP"

      ! CTFOK is deliberately not written, and cannot be. It is not a top-level
      ! section but a *subsection of* CTVEC: GAMESS's reader looks for it only
      ! immediately after reading a CTVEC block, and a standalone one falls through
      ! its keyword dispatcher to an abort. So writing the occupied orbital energies without CTVEC does
      ! not produce a file with one more section in it -- it produces a file
      ! GAMESS refuses to read at all. They are kept in the potential type
      ! against CTVEC being solved, and `pot%eps_occ` is what to emit then.

      ! Beta is frozen at one, as MAKEFP freezes it: its ICFIX flag fixes the
      ! prefactor and fits the exponent alone.
      ! Charge transfer, in the canonical-orbital form: the header carries the
      ! occupied count and the number of vectors, then the whole MO matrix in the
      ! same five-to-a-line layout the projection wavefunction uses. `CTFOK` is a
      ! *subsection* of this one, not a section of its own -- GAMESS's reader looks
      ! for it only directly behind a `CTVEC` block and aborts on a standalone one --
      ! so the two are written together and share one `STOP`.
      write (unit, "(A,I8,A,I8)") " CTVEC   ", pot%n_occ, "  ", pot%nao
      call write_wavefunction(unit, pot%canonical)
      write (unit, "(A)") " CTFOK   "
      call write_values(unit, pot%eps_occ, 16, 10, 4)
      write (unit, "(A)") " STOP"

      write (unit, "(A,F8.3,A)") "SCREEN2      (FROM VDWSCL=", pot%vdwscl, ")"
      do i = 1, pot%n_points
         write (unit, "(1X,A8,2F14.9)") pot%labels(i), 1.0_dp, pot%screen2(i)
      end do
      write (unit, "(A)") "STOP"

      write (unit, "(A,F8.3,A)") "SCREEN       (FROM VDWSCL=", pot%vdwscl, ")"
      do i = 1, pot%n_points
         write (unit, "(1X,A8,2F14.9)") pot%labels(i), 1.0_dp, pot%screen(i)
      end do
      write (unit, "(A)") "STOP"

      write (unit, "(A)") " $END"
      close (unit)

      if (present(omitted)) then
         omitted = "nothing"
      end if
   end subroutine write_efp_potential

   subroutine write_record(unit, label, values, width, decimals, per_line)
      !! A labelled record, continued with `>` and indented under its label
      integer, intent(in) :: unit
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: values(:)
      integer, intent(in) :: width, decimals, per_line

      character(len=MAX_LINE) :: line
      character(len=16) :: value_format
      integer :: i, first, last, at

      write (value_format, "(A,I0,A,I0,A)") "(F", width, ".", decimals, ")"
      first = 1
      do while (first <= size(values))
         last = min(first + per_line - 1, size(values))
         line = ""
         if (first == 1) then
            line(1:8) = label
         end if
         at = 9
         do i = first, last
            write (line(at:at + width - 1), value_format) values(i)
            at = at + width
         end do
         if (last < size(values)) line(at:at + 1) = " >"
         write (unit, "(A)") trim(line)
         first = last + 1
      end do
   end subroutine write_record

   subroutine write_tensor_point(unit, label, xyz, tensor, frequency)
      !! A point's coordinates on its label line, its nine components beneath
      integer, intent(in) :: unit
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: xyz(3)
      real(dp), intent(in) :: tensor(N_CART_PAIR)
      real(dp), intent(in), optional :: frequency

      if (present(frequency)) then
         ! "CT  1", five characters, then the coordinates -- not the label padded
         ! to eight the way the labelled record sections pad theirs.
         write (unit, "(A,3F15.10,A,F9.6,A)") trim(label), xyz, " -- FOR W=", &
            frequency, "I A.U."
      else
         ! trim, not a fixed width: this path serves both "CT1" from the static
         ! section and "CT  2" from a dynamic block's continuation points, and a
         ! fixed A3 would truncate the second.
         write (unit, "(A,3F15.10)") trim(label), xyz
      end if
      call write_values(unit, tensor, 16, 10, 4)
   end subroutine write_tensor_point

   subroutine write_values(unit, values, width, decimals, per_line)
      !! Unlabelled values, continued with `>`
      integer, intent(in) :: unit
      real(dp), intent(in) :: values(:)
      integer, intent(in) :: width, decimals, per_line

      character(len=MAX_LINE) :: line
      character(len=16) :: value_format
      integer :: i, first, last, at

      write (value_format, "(A,I0,A,I0,A)") "(F", width, ".", decimals, ")"
      first = 1
      do while (first <= size(values))
         last = min(first + per_line - 1, size(values))
         line = ""
         at = 1
         do i = first, last
            write (line(at:at + width - 1), value_format) values(i)
            at = at + width
         end do
         if (last < size(values)) line(at:at + 1) = " >"
         write (unit, "(A)") trim(line)
         first = last + 1
      end do
   end subroutine write_values

   subroutine write_lower_triangle(unit, matrix)
      !! The lower triangle, row by row, four values to a line
      integer, intent(in) :: unit
      real(dp), intent(in) :: matrix(:, :)

      real(dp), allocatable :: packed(:)
      integer :: n, i, j, at

      n = size(matrix, 1)
      allocate (packed(n*(n + 1)/2))
      at = 0
      do i = 1, n
         do j = 1, i
            at = at + 1
            packed(at) = matrix(i, j)
         end do
      end do
      call write_values(unit, packed, 16, 10, 4)
      deallocate (packed)
   end subroutine write_lower_triangle

   subroutine write_wavefunction(unit, orbitals)
      !! `PROJECTION WAVEFUNCTION`: five coefficients a line, orbital by orbital
      !!
      !! The line carries the orbital index and a chunk counter, which is what
      !! GAMESS's reader keys on, so the chunk counter restarts with each orbital.
      integer, intent(in) :: unit
      real(dp), intent(in) :: orbitals(:, :)

      integer :: nao, n_lmo, i, start, last, chunk, k
      character(len=MAX_LINE) :: line
      integer :: at

      nao = size(orbitals, 1)
      n_lmo = size(orbitals, 2)
      do i = 1, n_lmo
         chunk = 0
         start = 1
         do while (start <= nao)
            last = min(start + 4, nao)
            chunk = chunk + 1
            line = ""
            write (line(1:5), "(I2,I3)") i, chunk
            at = 6
            do k = start, last
               write (line(at:at + 14), "(ES15.8E2)") orbitals(k, i)
               at = at + 15
            end do
            write (unit, "(A)") trim(line)
            start = last + 1
         end do
      end do
   end subroutine write_wavefunction

end module mqc_efp_potential
