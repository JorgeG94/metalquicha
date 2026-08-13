!! A complete effective fragment potential, computed and written here
module mqc_efp_potential
   !! `RUNTYP=MAKEFP`: everything from a geometry and a basis name to a `.efp`
   !! file on disk, in Fortran, with no reference potential involved.
   !!
   !! **Why this module exists separately from the pieces it calls.** Each
   !! parameter block already had a `validation/check_*` program that computed it
   !! and dumped it for a Python comparison against GAMESS's printed numbers.
   !! That establishes that we compute the same quantities. It does not produce a
   !! potential: the only thing that had ever *written* one was a Python script
   !! that took GAMESS's own file and substituted our numbers into it, section by
   !! section. That was the right instrument for the question it answered -- does
   !! GAMESS agree with our parameters -- and it is not a makefp. The sections
   !! were computed in six unrelated programs, the file format lived in Python,
   !! and nothing assembled the two.
   !!
   !! So this is the assembly, and it is where a run turns into a fragment
   !! someone else can use.
   !!
   !! **Fifteen of the seventeen sections GAMESS's reader recognises.**
   !! `LMOQQPOL DYNAMIC
   !! POLARIZABLE POINTS` and `CTVEC` are
   !! not written, because we cannot reproduce GAMESS's values for them -- see the
   !! record in `mqc_libcint_cphf`, which ends in two span arguments showing the
   !! discrepancy is neither a rearrangement nor a translation of what we compute.
   !! `write_efp_potential` therefore reports what it omitted rather than emitting
   !! a plausible guess: a wrong dispersion tensor in a file someone runs is worse
   !! than an absent one, and their absence is legible to a reader of the file.
   !!
   !! **Where `LMOQQPOL` stands**, since it is now the only dispersion block
   !! missing and it is much closer than it was. Its 81 values are
   !! `QQL_SFT(3,3,3,3)` written with the last index fastest, and its write-time
   !! translation `QQSHIFT` (`efinp.src:13275`) takes the dipole-dipole *and*
   !! pre-shift dipole-quadrupole tensors as inputs -- both of which this module
   !! computes. Subtracting that shift from GAMESS's values recovers a tensor
   !! symmetric in both index pairs to 4.4e-16, which validates the formula and the
   !! index reading, and the recovered tensor is our quadrupole-quadrupole response
   !! times roughly 1/3.
   !!
   !! **That factor is not exact, and an earlier version of this comment said it
   !! was.** The median ratio is 0.3333 in every one of the four index-pattern
   !! groups, which is what made it look like a clean constant, but the best single
   !! scale is 0.3312 and the ratio scatters over [0.312, 0.354] between the 16th
   !! and 84th percentiles -- 5.8% residual, far above the 1.7e-04 this block's
   !! numbers are printed to. So 1/3 is the leading behaviour and something
   !! structural sits on top of it. Quoting the median as the relation was the same
   !! error as reading a global scale fit off a tensor whose components follow
   !! different conventions, which is what a scale fit did on the dipole-quadrupole
   !! block before the index order was pinned.
   !!
   !! **Summed over the orbitals the factor is clean and the gap shrinks.** The same
   !! move that cracked the dipole-quadrupole block -- summing to remove the
   !! projection -- gives a best scale of 0.333336, 1/3 to five digits, with the
   !! residual down to 2.0e-02 from 5.8e-02 per orbital. So 1/3 is real, and part of
   !! what is left is the per-orbital decomposition rather than the quantity.
   !!
   !! **And there is one concrete inconsistency to chase.** A summed response tensor
   !! must be symmetric under exchanging its two index pairs, because summed over
   !! every orbital the projector is the identity and the response function is
   !! symmetric in measure and respond. Ours is, to 1.1e-04. The tensor recovered
   !! from GAMESS by subtracting `QQSHIFT` is *not*: 2.2e-01. So the shift being
   !! subtracted carries a pair-asymmetric error, and the formula shows where it can
   !! come from -- `RRalphdel1` contracts the dipole-dipole tensor on its first index
   !! (`DD(i,e)`) while its mirror `RRalphdel2` contracts on the second (`DD(b,i)`),
   !! so the two are mirror images only if `DD` is symmetric, and per orbital it is
   !! not. That is the thread, and it is the same class of thing as the transposed
   !! `alpha` that closed `DQSHIFT` -- but transposing these two terms was tried and
   !! makes it worse, so it is not the same fix.
   !!
   !! With the factor the forward recipe reaches 1.3e-01 relative, and the residual
   !! sits on the components the `delta` terms touch: 1.6e-01 where both index pairs
   !! are diagonal against 2.6e-02 where neither is. Transposing the dipole-dipole
   !! tensor in any combination of those terms makes it worse, which is the fix that
   !! worked for `DQSHIFT`. And a least-squares fit of all six term coefficients over
   !! all 3888 values returns no clean fraction and only improves to 1.28e-01, so the
   !! reference is not in the span of our tensor plus the five shift terms -- the
   !! quadrupole-quadrupole response itself still differs, as the dipole-quadrupole
   !! one appeared to before its operator convention was found. Not written until
   !! that closes.
   !!
   !! Seventeen and not eighteen because `CTFOK` is not a section. It is a
   !! subsection of `CTVEC`, accepted only directly behind one, so it goes when
   !! `CTVEC` goes -- see the comment where it would have been written. That was
   !! found by handing a file to GAMESS, not by reading its output: every
   !! parameter in it agreed with GAMESS's own and the file was still unreadable.
   !!
   !! What the fifteen support, in EFP terms: electrostatics to octupole with
   !! charge-penetration screening, polarization, exchange repulsion, and
   !! dispersion through its `E6` and `E7` terms. What they do not is `E8`, which
   !! needs the quadrupole-quadrupole block, and charge transfer, which needs
   !! `CTVEC`.
   !!
   !! **The formats are GAMESS's reader's, not its writer's.** Free-form values in
   !! measured columns with the continuation markers each section expects. Byte
   !! identity with MAKEFP's output would be a much larger job for no extra
   !! confidence, since the test that matters is whether GAMESS accepts the file
   !! and agrees with the energies it computes from it.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_mass
   use mqc_cgto, only: molecular_basis_type
   use mqc_basis_utils, only: find_basis_file
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule, shell_dim
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_localize, only: boys_localize
   use mqc_libcint_dma, only: dma_result_t, distributed_multipoles
   use mqc_libcint_cphf, only: distributed_polarizability, &
                               distributed_dynamic_polarizability, &
                               distributed_dynamic_cross, &
                               casimir_polder_frequencies, N_CASIMIR_POLDER
   use mqc_libcint_multipole, only: multipole_matrices
   use mqc_libcint_screening, only: fit_screening, SCREEN_EXPONENTIAL, SCREEN_GAUSSIAN
   use libcint_fortran, only: LIBCINT_ANG_OF
   implicit none
   private

   public :: efp_potential_t
   public :: make_efp_potential
   public :: write_efp_potential

   !> Longest line any section emits, with room to spare.
   integer, parameter :: MAX_LINE = 160

   !> Row and column of each of GAMESS's nine polarizability slots. The
   !> off-diagonal triples are the transpose of what its labels suggest; that was
   !> measured in `validation/check_distributed_polarizability.py` rather than
   !> assumed, and it is the one convention here that a symmetric test tensor
   !> would not have caught.
   integer, parameter :: POL_ROW(9) = [1, 2, 3, 2, 3, 3, 1, 1, 2]
   integer, parameter :: POL_COL(9) = [1, 2, 3, 1, 1, 2, 2, 3, 3]

   !> libcint's index for each of GAMESS's six Cartesian d slots, and the
   !> normalization between the two codes' d functions. Both established in
   !> `validation/check_projection.py` against GAMESS's own coefficients.
   !> libcint's full-Cartesian quadrupole slots, which run xx,xy,xz,yx,...,zz.
   integer, parameter :: QXX = 1, QXY = 2, QXZ = 3, QYX = 4, QYY = 5
   integer, parameter :: QYZ = 6, QZX = 7, QZY = 8, QZZ = 9

   integer, parameter :: D_FROM_LIBCINT(6) = [1, 4, 6, 2, 3, 5]
   real(dp), parameter :: D_NORMALIZATION = 1.585330892_dp

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
      real(dp) :: vdwscl = 0.7_dp !! The screening grid's van der Waals scale
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
      real(dp), allocatable :: dynamic_pol(:, :, :, :)!! (3, 3, n_lmo, n_freq)
      real(dp), allocatable :: dipquad(:, :, :, :, :)
         !! `(3, 3, 3, n_lmo, n_freq)` as `A'(a,b,c)`, **after** the write-time
         !! translation to each centroid. Post-shift because that is the form the
         !! file carries and the shift needs the dipole-dipole tensors, so keeping
         !! the pre-shift tensor would mean carrying its inputs too.
      real(dp), allocatable :: frequencies(:)         !! Imaginary, a.u.
      real(dp), allocatable :: fock_lmo(:, :)         !! (n_lmo, n_lmo)
      real(dp), allocatable :: orbitals(:, :)         !! LMOs in GAMESS's AO order
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
      if (allocated(self%frequencies)) deallocate (self%frequencies)
      if (allocated(self%fock_lmo)) deallocate (self%fock_lmo)
      if (allocated(self%orbitals)) deallocate (self%orbitals)
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
                                 basis_name, name, pot, error, n_core, vdwscl, &
                                 verbose)
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
      integer, intent(in), optional :: n_core
         !! Orbitals excluded from the localized set. Default is the standard
         !! frozen core, which is what MAKEFP uses: its polarizable points and
         !! its exchange-repulsion orbitals are valence only.
      real(dp), intent(in), optional :: vdwscl
      logical, intent(in), optional :: verbose

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(dma_result_t) :: dma
      real(dp), allocatable :: loc(:, :), ovl(:, :), sc(:, :), w(:, :), scaled(:, :)
      real(dp), allocatable :: alpha(:)
      integer :: natm, core, i, j, k, n_valence
      logical :: talk

      talk = .false.
      if (present(verbose)) talk = verbose
      natm = size(atomic_numbers)
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

      if (talk) write (*, "(A,A,A,I0,A)") "  basis ", trim(basis_name), ", ", &
         mol%nao, " functions"

      call run_libcint_rhf(mol, sum(atomic_numbers), 200, 1.0e-12_dp, 1.0e-10_dp, &
                           .false., scf, error)
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
      if (talk) write (*, "(A,F18.10)") "  RHF energy ", scf%energy

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
      if (talk) write (*, "(A,I0,A,I0,A)") "  localized ", n_valence, " of ", &
         pot%n_occ, " occupied orbitals"

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
      call distributed_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                      pot%n_occ, pot%static_pol, pot%centroids, &
                                      error, n_core=core)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if

      allocate (pot%frequencies(N_CASIMIR_POLDER))
      pot%frequencies = casimir_polder_frequencies()
      call distributed_dynamic_polarizability(mol, scf%orbitals, &
                                              scf%orbital_energies, pot%n_occ, &
                                              pot%frequencies, pot%dynamic_pol, &
                                              pot%centroids, error, n_core=core)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      if (talk) write (*, "(A,I0,A)") "  polarizabilities at ", N_CASIMIR_POLDER, &
         " imaginary frequencies"

      call dipole_quadrupole_block(mol, scf, coordinates, atomic_numbers, core, pot, &
                                   error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      if (talk) write (*, "(A)") "  dipole-quadrupole dispersion tensors"

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

      call projection_basis_lines(atomic_numbers, pot%labels, pot%points, &
                                  element_symbols, basis_name, pot%basis_lines, &
                                  error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if

      ! --- charge penetration screening ----------------------------------------
      ! Last, because it is fitted to the error the multipoles above make.
      call fit_screening(mol, scf%density, dma, atomic_numbers, SCREEN_EXPONENTIAL, &
                         alpha, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      allocate (pot%screen2(pot%n_points))
      pot%screen2 = alpha
      deallocate (alpha)
      call fit_screening(mol, scf%density, dma, atomic_numbers, SCREEN_GAUSSIAN, &
                         alpha, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      allocate (pot%screen(pot%n_points))
      pot%screen = alpha
      deallocate (alpha)
      if (talk) write (*, "(A)") "  screening fitted for both damping forms"

      call mol%destroy()
      deallocate (loc, ovl, sc, w, scaled)
   end subroutine make_efp_potential

   subroutine dipole_quadrupole_block(mol, scf, coordinates, atomic_numbers, core, &
                                      pot, error)
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
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dip(:, :, :), quad(:, :, :), buck(:, :, :)
      real(dp), allocatable :: raw(:, :, :, :), centroids(:, :)
      real(dp) :: com(3), mass_total, r(3), alpha(3, 3), isotropic
      integer :: i, a, b, c, d, k, f, n_freq

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
      ! (`prpel.src:5625`), kept as all nine Cartesian slots so that no expansion
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

      call distributed_dynamic_cross(mol, scf%orbitals, scf%orbital_energies, &
                                     pot%n_occ, pot%frequencies, buck, dip, raw, &
                                     centroids, error, n_core=core)
      if (error%has_error()) return

      n_freq = size(pot%frequencies)
      allocate (pot%dipquad(3, 3, 3, pot%n_lmo, n_freq))
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

      deallocate (dip, quad, buck, raw, centroids)
   end subroutine dipole_quadrupole_block

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
         else if (l >= 2) then
            call error%set(ERROR_VALIDATION, &
                           "makefp: no basis function map for l > 1 beyond Cartesian d")
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

      integer :: unit, i, k, f, stat, a, b, c
      character(len=8) :: label
      real(dp) :: tensor(9)
      real(dp) :: wide(27)

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
      ! immediately after reading a CTVEC block (`efinp.src`, RDCANV, which
      ! aborts with "SUBSECTION 'CTFOK' ... MUST FOLLOW RIGHT BEHIND ANY 'CTVEC'
      ! SUBSECTION"), and a standalone one falls through its keyword dispatcher
      ! to an abort. So writing the occupied orbital energies without CTVEC does
      ! not produce a file with one more section in it -- it produces a file
      ! GAMESS refuses to read at all. They are kept in the potential type
      ! against CTVEC being solved, and `pot%eps_occ` is what to emit then.

      ! Beta is frozen at one, as MAKEFP freezes it: its ICFIX flag fixes the
      ! prefactor and fits the exponent alone.
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
         omitted = "LMOQQPOL DYNAMIC POLARIZABLE POINTS, CTVEC (with CTFOK)"
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
      real(dp), intent(in) :: tensor(9)
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
