!! Gaussian integrals on the CPU, through libcint
module mqc_libcint_integrals
   !! Packs a molecule and its basis into libcint's arrays, and hands back the
   !! matrices an SCF needs: overlap, core Hamiltonian, and the two-electron
   !! integrals.
   !!
   !! **This exists to be checked against, not to be fast.** Everything ab
   !! initio here runs on cuEST, which needs an A100, so the HF and DFT paths
   !! are compile-checked and never executed where anyone is working. A CPU
   !! path makes them runnable on a laptop and gives the GPU one a second
   !! implementation to disagree with.
   !!
   !! Two things about libcint's conventions are worth knowing before reading
   !! the packing, because both fail quietly:
   !!
   !!   * **Contraction coefficients must arrive pre-multiplied by each
   !!     primitive's normalisation.** libcint does not apply it. Without it
   !!     every integral is wrong by a constant per shell, nothing errors, and
   !!     the symptom is an overlap diagonal that is not 1.
   !!   * **The slot constants from `libcint_fortran` are already 1-based.**
   !!     The C header's are 0-based and the Fortran interface converts them;
   !!     adding one again writes into the neighbouring slot and leaves the
   !!     real one uninitialised, which surfaces as a crash inside libcint.
   !!
   !! The ERIs are held in core, all n^4 of them. That is the right choice for
   !! what this is for -- a hundred basis functions is 800 MB and a reference
   !! calculation is small by construction -- and the wrong one for anything
   !! else. Direct or density-fitted assembly is a different backend.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_timing, only: timing_report_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cgto, only: molecular_basis_type, atomic_basis_type
   use mqc_basis_utils, only: find_basis_file
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use pic_lapack_interfaces, only: pic_syev
   use libcint_fortran, only: libcint_1e_ovlp_sph, libcint_1e_kin_sph, &
                              libcint_3c2e_sph, libcint_2c2e_sph, &
                              libcint_1e_nuc_sph, libcint_2e_sph, &
                              libcint_cgto_sph, libcint_tot_cgto_sph, &
                              libcint_1e_ovlp_cart, libcint_1e_kin_cart, &
                              libcint_3c2e_cart, libcint_2c2e_cart, &
                              libcint_1e_nuc_cart, libcint_2e_cart, &
                              libcint_cgto_cart, libcint_tot_cgto_cart, &
                              libcint_2e_sph_optimizer, libcint_2e_cart_optimizer, &
                              libcint_3c2e_sph_optimizer, libcint_3c2e_cart_optimizer, &
                              libcint_2c2e_sph_optimizer, libcint_2c2e_cart_optimizer, &
                              libcint_gto_norm, libcint_del_optimizer, &
                              LIBCINT_PTR_RANGE_OMEGA, &
                              LIBCINT_ATM_SLOTS, LIBCINT_BAS_SLOTS, &
                              LIBCINT_CHARGE_OF, LIBCINT_PTR_COORD, &
                              LIBCINT_ATOM_OF, LIBCINT_ANG_OF, LIBCINT_NPRIM_OF, &
                              LIBCINT_NCTR_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_COEFF, &
                              LIBCINT_PTR_ENV_START, LIBCINT_KAPPA_OF
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
   implicit none
   private

   integer, parameter :: KAPPA_SP_SHELL = 64
      !! What libfint reads in KAPPA_OF to mean "this shell holds s and p".
      !!
      !! A wire-format constant, like the slot indices this file already
      !! imports, and here for the same reason: it is part of the layout of
      !! `bas`, not something the module that fills `bas` can ask for. Defined
      !! at libfint's `src/cint_bas.f90:61`; if that number moves, this is the
      !! other end of the agreement.
      !!
      !! KAPPA_OF is the relativistic kappa, which a non-relativistic shell
      !! leaves at zero -- which is what frees the slot to carry this, and also
      !! why it is written only under MQC_WITH_SP_SHELLS. libcint would read a
      !! 64 there as a real kappa and act on it.

   public :: libcint_molecule_t
   public :: build_libcint_molecule
   public :: build_df_tensor
   public :: build_df_mo_tensor
   public :: build_df_mo_block
   ! The one place that maps "is this basis Cartesian?" onto libcint's two sets
   ! of entry points. Exported because the direct Fock build needs the same
   ! mapping, and two copies of it is two chances to route half the calls.
   public :: shell_dim
   public :: pair_index
   ! The shell set a four-centre integral loop should run over, and the
   ! Schwarz bounds re-blocked to match it. Exported because the direct Fock
   ! build runs the same loops over the same choice, and two copies of the
   ! choice is two chances to hand a fused shell to a driver that cannot
   ! read one.
   public :: eri_shell_table_t
   public :: eri_shell_table
   public :: eri_schwarz_collapse
   public :: two_electron_block
   public :: two_electron_optimizer
   public :: ket_transformed_pairs
   ! Where an atom's functions live, and what angular momentum each carries.
   ! Both follow from the packing order in `molecule_build`, so they belong
   ! beside it rather than being re-derived by whatever needs them.
   public :: atom_ao_blocks
   public :: build_df_shell_table
   public :: mixed_basis_overlap
   ! The undifferentiated three- and two-centre integrals, and the metric's
   ! inverse square root. Public because the gradient needs the same three
   ! quantities the fitted Fock build does -- and, for the metric, needs them
   ! from the *same* eigenvalue threshold: a gradient built on a pseudo-inverse
   ! that kept different modes than the SCF's would not differentiate the
   ! energy the SCF actually converged.
   public :: three_centre
   public :: two_centre
   public :: metric_inverse_sqrt
   public :: subshell_layout

   type :: libcint_molecule_t
      !! One molecule, packed the way libcint wants it
      integer :: natm = 0
      integer :: nbas = 0
      integer :: nao = 0
      logical :: cartesian = .false.
         !! Whether every integral over this molecule must use libcint's
         !! Cartesian entry points rather than its spherical ones. Set from the
         !! basis file's `function_type`, and one value for the whole molecule:
         !! libcint chooses the form per call, so a d shell cannot be Cartesian
         !! while an f shell beside it is spherical. It also changes `nao` --
         !! six functions per d shell instead of five -- so anything sized from
         !! `nao` follows it automatically, and anything that calls `_sph`
         !! directly does not.
      integer, allocatable :: atm(:, :)
      integer, allocatable :: bas(:, :)
      real(dp), allocatable :: env(:)
      integer, allocatable :: shell_offset(:)  !! First AO of each shell, 0-based
      ! The fused-sp view of the same shells: every consecutive s/p pair on
      ! shared exponents -- an L shell, which is what every Pople basis is
      ! written in -- collapsed to one `bas` entry that libfint's `int2e`
      ! evaluates in a single pass over the shared pair data. Built by
      ! `build_sp_view`, and empty (`nbas_sp` zero) whenever the molecule has
      ! no L shell or the build is against libcint, which cannot read the
      ! marker. The AO order of the view is identical to the split order --
      ! one s function then three p -- so anything indexed by AO needs no map.
      integer :: nbas_sp = 0
      integer, allocatable :: bas_sp(:, :)
      real(dp), allocatable :: env_sp(:)
         !! The split `env` with the fused coefficient blocks appended, so the
         !! atom coordinate pointers in `atm` stay valid against either table
      integer, allocatable :: shell_offset_sp(:)  !! First AO per view shell, 0-based
      integer, allocatable :: sp_split_first(:)
         !! First split shell of each view shell, 1-based, with an (nbas_sp+1)
         !! sentinel -- the map `eri_schwarz_collapse` re-blocks bounds through
      real(dp), allocatable :: charges(:)      !! Nuclear charges, for repulsion
      real(dp), allocatable :: coords(:, :)    !! (3, natm), Bohr
   contains
      procedure :: build => molecule_build
      procedure :: overlap => molecule_overlap
      procedure :: core_hamiltonian => molecule_core_hamiltonian
      procedure :: kinetic => molecule_kinetic
      procedure :: eris => molecule_eris
      procedure :: eris_packed => molecule_eris_packed
      procedure :: nuclear_repulsion => molecule_nuclear_repulsion
      procedure :: atom_subset => molecule_atom_subset
      procedure :: destroy => molecule_destroy
   end type libcint_molecule_t

   type :: eri_shell_table_t
      !! The shell table a four-centre integral loop runs over
      !!
      !! Either the molecule's own shells or its fused-sp view, chosen once in
      !! `eri_shell_table` so a loop body cannot mix the two. Widths and
      !! offsets are carried explicitly: inside the loops a shell's dimension
      !! is a table lookup, not a `cgto` call per quartet.
      integer :: nbas = 0
      integer :: block_max = 1                 !! Largest shell width, for buffers
      integer, allocatable :: bas(:, :)
      real(dp), allocatable :: env(:)
      integer, allocatable :: offs(:)          !! (nbas+1) first AO per shell, 0-based
      integer, allocatable :: dims(:)          !! (nbas) functions per shell
   end type eri_shell_table_t

contains

   ! ---------------------------------------------------------------------------
   ! Convention routing
   !
   ! libcint spells the two angular conventions as two families of entry points
   ! -- cint2e_sph and cint2e_cart, and so on -- with identical signatures. The
   ! choice is per call, which is why it has to be one choice for a whole
   ! molecule, and these four functions are the only places the choice is made.
   ! Everything else asks for "this molecule's" integrals and gets the right
   ! family.
   ! ---------------------------------------------------------------------------

   function shell_dim(cartesian, bas_id, bas) result(dim)
      !! How many functions one shell contributes, in the given convention
      !!
      !! Six for a Cartesian d shell against five for a spherical one, and the
      !! gap widens with l. This is where the basis function count comes from,
      !! so getting it wrong does not produce a slightly wrong energy -- it
      !! produces a consistent SCF over the wrong basis.
      logical, intent(in) :: cartesian
      integer, intent(in) :: bas_id      !! 0-based, as libcint counts
      integer, intent(in) :: bas(:, :)
      integer :: dim

      if (cartesian) then
         dim = libcint_cgto_cart(bas_id, bas)
      else
         dim = libcint_cgto_sph(bas_id, bas)
      end if
   end function shell_dim

   function total_shell_dim(cartesian, bas, nbas) result(ntot)
      !! Total basis functions over all shells, in the given convention
      logical, intent(in) :: cartesian
      integer, intent(in) :: bas(:, :)
      integer, intent(in) :: nbas
      integer :: ntot

      if (cartesian) then
         ntot = libcint_tot_cgto_cart(bas, nbas)
      else
         ntot = libcint_tot_cgto_sph(bas, nbas)
      end if
   end function total_shell_dim

   function two_electron_block(cartesian, buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
      !! One (ij|kl) shell quartet, in the given convention
      !!
      !! `buf` is assumed-shape rather than the assumed-size the libcint
      !! interface declares. Callers pass a whole contiguous scratch array
      !! sized for the largest quartet, so the two are interchangeable here,
      !! and an explicit shape is what the style guide asks for.
      logical, intent(in) :: cartesian
      real(dp), intent(out) :: buf(:)
      integer, intent(in) :: shls(4)
      integer, intent(in) :: atm(:, :)
      integer, intent(in) :: natm
      integer, intent(in) :: bas(:, :)
      integer, intent(in) :: nbas
      real(dp), intent(in) :: env(:)
      type(c_ptr), intent(in), optional :: opt
      integer :: ret

      if (cartesian) then
         if (present(opt)) then
            ret = libcint_2e_cart(buf, shls, atm, natm, bas, nbas, env, opt)
         else
            ret = libcint_2e_cart(buf, shls, atm, natm, bas, nbas, env)
         end if
      else
         if (present(opt)) then
            ret = libcint_2e_sph(buf, shls, atm, natm, bas, nbas, env, opt)
         else
            ret = libcint_2e_sph(buf, shls, atm, natm, bas, nbas, env)
         end if
      end if
   end function two_electron_block

   subroutine two_electron_optimizer(cartesian, opt, atm, natm, bas, nbas, env)
      !! The per-shell-pair precomputation, for the matching convention
      !!
      !! The optimizer is not interchangeable between the two families: it
      !! caches data laid out for the convention it was built for.
      logical, intent(in) :: cartesian
      type(c_ptr), intent(inout) :: opt
      integer, intent(in) :: atm(:, :)
      integer, intent(in) :: natm
      integer, intent(in) :: bas(:, :)
      integer, intent(in) :: nbas
      real(dp), intent(in) :: env(:)

      if (cartesian) then
         call libcint_2e_cart_optimizer(opt, atm, natm, bas, nbas, env)
      else
         call libcint_2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
      end if
   end subroutine two_electron_optimizer

   subroutine build_libcint_molecule(atomic_numbers, element_symbols, coordinates, &
                                     basis_name, mol, error, normalize_contractions, &
                                     force_cartesian, ghost)
      !! A molecule from a basis set *name*, through the ordinary reader
      !!
      !! This is what makes the backend general rather than a demonstration:
      !! any of the basis sets in `basis_sets/` rather than whatever was
      !! typed into a test.
      !!
      !! **Two normalisations, and they are not the same thing.** The BSE files
      !! give contraction coefficients against normalised *primitives*, and
      !! `molecule_build` folds in `libcint_gto_norm`, which is the primitive
      !! convention libcint asks for. That leaves the *contracted* function
      !! normalised only to the precision the coefficients were published at --
      !! about 1e-6 for cc-pVDZ -- so `<chi|chi>` is 1.000001, not 1.
      !!
      !! This comment used to say applying a second normalisation "would count
      !! it twice", which conflated the two. It does not: the contraction norm
      !! is a separate sum over primitive pairs, and PySCF applies it.
      !!
      !! Why it went unnoticed for so long is worth stating, because it says
      !! which tests can and cannot see it: an SCF energy is **invariant** to the
      !! normalisation of a basis function. Scaling one does not change the space
      !! the basis spans, so the orbital coefficients absorb it and the energy is
      !! identical. Every energy validated against PySCF to 1e-9 was therefore
      !! correct and stayed correct. It shows up only in per-AO quantities -- an
      !! overlap diagonal, a multipole matrix element, anything an effective
      !! fragment potential is made of -- and nothing looked at one until the
      !! multipole integrals arrived.
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)   !! (3, natm), Bohr
      character(len=*), intent(in) :: basis_name
      type(libcint_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: ghost(:)
         !! Atoms keeping their basis and losing their nuclear charge; see
         !! `molecule_build`.
      logical, intent(in), optional :: force_cartesian
         !! Read the basis in Cartesian form whatever the file declares. See the note
         !! in `build`: BSE is inconsistent about Pople sets and GAMESS assumes
         !! Cartesian for them.
      logical, intent(in), optional :: normalize_contractions
         !! Scale each contracted shell so `<chi|chi>` is exactly one. Default
         !! true, which is both physically right and what PySCF does. False
         !! reproduces the behaviour from before this was fixed, which is worth
         !! being able to ask for: it is how the effect on any given quantity can
         !! be measured rather than argued about.

      type(molecular_basis_type) :: basis
      character(len=:), allocatable :: path
      type(error_t) :: read_error

      call find_basis_file(basis_name, path, read_error)
      if (read_error%has_error()) then
         call error%set(ERROR_VALIDATION, "no basis set file for '"//trim(basis_name)// &
                        "': "//read_error%get_message())
         return
      end if

      call build_molecular_basis_json(path, element_symbols, basis, read_error)
      if (read_error%has_error()) then
         call error%set(ERROR_VALIDATION, "could not read "//trim(basis_name)//": "// &
                        read_error%get_message())
         return
      end if

      if (present(ghost)) then
         call mol%build(atomic_numbers, coordinates, basis, error, &
                        normalize_contractions=normalize_contractions, &
                        force_cartesian=force_cartesian, ghost=ghost)
      else
         call mol%build(atomic_numbers, coordinates, basis, error, &
                        normalize_contractions=normalize_contractions, &
                        force_cartesian=force_cartesian)
      end if
      call basis%destroy()
   end subroutine build_libcint_molecule

   function contraction_group_size(atom_basis, first) result(nctr)
      !! How many consecutive shells from `first` are one general contraction
      !!
      !! The reader emits one shell per coefficient column, because that is what
      !! a Basis Set Exchange entry lists. For a general contraction -- cc-pVDZ
      !! oxygen is nine s primitives carrying three columns -- those columns
      !! share one set of exponents, and libcint can take them as a single shell
      !! with `NCTR_OF` set to the column count. It then evaluates the primitives
      !! once and contracts them into every column, instead of repeating the
      !! primitive work per column.
      !!
      !! Only *consecutive* shells merge, and that is the whole safety argument.
      !! libcint lays a shell's functions out with the contraction index
      !! outermost, so merging columns that were already adjacent reproduces the
      !! basis function order exactly -- the same functions, the same sequence,
      !! the same matrix. Merging across a gap would silently permute the AO
      !! basis, and every matrix built on it.
      !!
      !! An SP shell never merges: its columns carry different angular momenta,
      !! which is the first thing tested here.
      type(atomic_basis_type), intent(in) :: atom_basis
      integer, intent(in) :: first
      integer :: nctr

      integer :: k

      nctr = 1
      do k = first + 1, atom_basis%nshells
         if (atom_basis%shells(k)%ang_mom /= atom_basis%shells(first)%ang_mom) exit
         if (atom_basis%shells(k)%nfunc /= atom_basis%shells(first)%nfunc) exit
         ! Bit-for-bit: both columns were parsed from the same exponent strings
         ! in the same file, so anything short of equality means they are not
         ! the same primitives and must not share a shell.
         if (any(atom_basis%shells(k)%exponents /= atom_basis%shells(first)%exponents)) exit
         nctr = nctr + 1
      end do
   end function contraction_group_size

   subroutine molecule_build(this, atomic_numbers, coordinates, basis, error, &
                             normalize_contractions, force_cartesian, ghost)
      !! Pack atoms and shells into libcint's atm/bas/env
      class(libcint_molecule_t), intent(inout) :: this
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)   !! (3, natm), Bohr
      type(molecular_basis_type), intent(in) :: basis
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: force_cartesian  !! See the note below
      logical, intent(in), optional :: normalize_contractions
         !! See `build_libcint_molecule`. Default true.
      logical, intent(in), optional :: ghost(:)
         !! Which atoms carry their basis functions but no nuclear charge.
         !!
         !! A counterpoise-corrected monomer is exactly this: every basis
         !! function of the dimer, and only its own nuclei. The basis is taken
         !! from `basis` and is untouched here, so ghosting changes the nuclear
         !! attraction and the nuclear repulsion and nothing else -- in
         !! particular the AO count and ordering are identical to the unghosted
         !! molecule's, which is what lets matrices from two of them be
         !! contracted against each other.

      logical :: do_normalize
      integer :: jprim
      real(dp) :: norm2, scale, ai, aj

      integer :: iatom, ishell, iprim, nprim, off, env_size, ang
      integer :: shell_index, ictr, nctr, first, z_eff

      do_normalize = .true.
      if (present(normalize_contractions)) do_normalize = normalize_contractions

      this%natm = size(atomic_numbers)
      ! What the basis file said, carried through the reader. Every integral
      ! call below routes on it, and so does nao.
      !
      ! `force_cartesian` overrides it, and the reason it exists is that the Basis Set
      ! Exchange is not consistent about Pople sets: 6-31G* declares its d Cartesian,
      ! while 6-31G(2df,p) and 6-311++G** declare theirs spherical. The Pople
      ! convention is Cartesian, and it is what GAMESS assumes -- its ISPHER default
      ! is -1 -- so reading one of those as spherical is a different calculation from
      ! the one anybody comparing against GAMESS means. Not automatic: it changes the
      ! variational space, and for a Dunning set spherical is right and GAMESS uses it
      ! too, so which one is correct depends on the basis and the caller has to say.
      this%cartesian = basis%is_cartesian()
      if (present(force_cartesian)) then
         if (force_cartesian) this%cartesian = .true.
      end if
      if (basis%nelements /= this%natm) then
         call error%set(ERROR_VALIDATION, "libcint: the basis covers a different "// &
                        "number of atoms than the geometry has")
         return
      end if

      ! Count the shells libcint will see, which is fewer than the reader
      ! produced wherever a general contraction was split into one shell per
      ! coefficient column. See `contraction_group_size` for why they merge.
      this%nbas = 0
      env_size = LIBCINT_PTR_ENV_START + 3*this%natm
      do iatom = 1, this%natm
         ishell = 1
         do while (ishell <= basis%elements(iatom)%nshells)
            nctr = contraction_group_size(basis%elements(iatom), ishell)
            nprim = basis%elements(iatom)%shells(ishell)%nfunc
            this%nbas = this%nbas + 1
            ! One set of exponents for the group, one coefficient column each.
            env_size = env_size + nprim + nprim*nctr
            ishell = ishell + nctr
         end do
      end do

      allocate (this%atm(LIBCINT_ATM_SLOTS, this%natm))
      allocate (this%bas(LIBCINT_BAS_SLOTS, this%nbas))
      allocate (this%env(env_size))
      allocate (this%shell_offset(this%nbas + 1))
      allocate (this%charges(this%natm))
      allocate (this%coords(3, this%natm))

      ! Zeroed rather than merely filled: libcint reads slots nothing here
      ! sets -- NUC_MOD_OF, PTR_ZETA, KAPPA_OF -- and stack garbage in a
      ! pointer slot crashes inside the library with nothing to say why.
      this%atm = 0
      this%bas = 0
      this%env = 0.0_dp

      off = LIBCINT_PTR_ENV_START
      do iatom = 1, this%natm
         z_eff = atomic_numbers(iatom)
         if (present(ghost)) then
            if (ghost(iatom)) z_eff = 0
         end if
         this%atm(LIBCINT_CHARGE_OF, iatom) = z_eff
         this%atm(LIBCINT_PTR_COORD, iatom) = off
         this%env(off + 1:off + 3) = coordinates(1:3, iatom)
         this%charges(iatom) = real(z_eff, dp)
         this%coords(:, iatom) = coordinates(1:3, iatom)
         off = off + 3
      end do

      shell_index = 0
      do iatom = 1, this%natm
         first = 1
         do while (first <= basis%elements(iatom)%nshells)
            nctr = contraction_group_size(basis%elements(iatom), first)
            shell_index = shell_index + 1
            ang = basis%elements(iatom)%shells(first)%ang_mom
            nprim = basis%elements(iatom)%shells(first)%nfunc

            this%bas(LIBCINT_ATOM_OF, shell_index) = iatom - 1   ! libcint counts from 0
            this%bas(LIBCINT_ANG_OF, shell_index) = ang
            this%bas(LIBCINT_NPRIM_OF, shell_index) = nprim
            this%bas(LIBCINT_NCTR_OF, shell_index) = nctr

            this%bas(LIBCINT_PTR_EXP, shell_index) = off
            this%env(off + 1:off + nprim) = basis%elements(iatom)%shells(first)%exponents
            off = off + nprim

            ! libcint reads the coefficients as an (nprim, nctr) matrix with
            ! stride nprim -- `coeff[nprim*ictr + iprim]` in g1e.c -- so the
            ! columns go down one after another, in the order the reader
            ! produced them. That order is what keeps the basis functions where
            ! they were: a shell with nctr contractions emits its contractions
            ! outermost, which is the same sequence the split shells gave.
            this%bas(LIBCINT_PTR_COEFF, shell_index) = off
            do ictr = 1, nctr
               do iprim = 1, nprim
                  ! The normalisation libcint expects to have been applied. It
                  ! depends on l and the exponent alone, so every column of a
                  ! group takes the same factor.
                  this%env(off + iprim) = &
                     basis%elements(iatom)%shells(first + ictr - 1)%coefficients(iprim) &
                     *libcint_gto_norm(ang, basis%elements(iatom)%shells(first)%exponents(iprim))
               end do

               ! And now the *contraction* normalisation, which is a different
               ! sum: <chi|chi> over primitive pairs, where two normalised
               ! primitives of the same l on the same centre overlap as
               !
               !     S_ij = (2 sqrt(a_i a_j) / (a_i + a_j))^(l + 3/2)
               !
               ! A single primitive gives S = 1 and a norm of exactly one, which
               ! is why uncontracted shells were already right and only the
               ! contracted ones were off. Done per column, since each
               ! contraction has its own coefficients and its own norm.
               if (do_normalize) then
                  norm2 = 0.0_dp
                  do iprim = 1, nprim
                     do jprim = 1, nprim
                        ai = basis%elements(iatom)%shells(first)%exponents(iprim)
                        aj = basis%elements(iatom)%shells(first)%exponents(jprim)
                        norm2 = norm2 + this%env(off + iprim)*this%env(off + jprim) &
                                /(libcint_gto_norm(ang, ai)*libcint_gto_norm(ang, aj)) &
                                *(2.0_dp*sqrt(ai*aj)/(ai + aj))**(real(ang, dp) + 1.5_dp)
                     end do
                  end do
                  if (norm2 > 0.0_dp) then
                     scale = 1.0_dp/sqrt(norm2)
                     do iprim = 1, nprim
                        this%env(off + iprim) = this%env(off + iprim)*scale
                     end do
                  end if
               end if

               off = off + nprim
            end do

            first = first + nctr
         end do
      end do

      ! Where each shell's functions start, so a shell-pair block knows where
      ! it lands in the matrix.
      this%shell_offset(1) = 0
      do shell_index = 1, this%nbas
         this%shell_offset(shell_index + 1) = this%shell_offset(shell_index) &
                                              + shell_dim(this%cartesian, shell_index - 1, this%bas)
      end do
      this%nao = this%shell_offset(this%nbas + 1)

      if (this%nao /= total_shell_dim(this%cartesian, this%bas, this%nbas)) then
         call error%set(ERROR_VALIDATION, "libcint: shell offsets disagree with the "// &
                        "basis function count")
         return
      end if

      call build_sp_view(this)
   end subroutine molecule_build

#ifdef MQC_WITH_SP_SHELLS
   subroutine build_sp_view(this)
      !! Fuse each s/p pair on shared exponents into one L shell, as a view
      !!
      !! A 6-31G Fock build spends its time recomputing, for the p half of
      !! every valence shell, the exponential prefactors, pair data and Rys
      !! roots it just computed for the s half -- the halves share exponents,
      !! and splitting them is purely an artefact of `bas` having one ANG_OF
      !! slot. libfint reads a marker in KAPPA_OF (see `KAPPA_SP_SHELL`) and
      !! evaluates both halves in one pass, measured at 2-3x on a Pople basis.
      !!
      !! A VIEW, not the representation: the split `bas` stays the molecule's
      !! own, because libfint carries an L shell through `int2e_sph` and
      !! `int2e_cart` only. Every other driver -- one-electron, three-centre,
      !! every derivative -- deliberately behaves as if L shells did not
      !! exist, and would read the fused entry as a plain p shell over the s
      !! coefficients: not an error, a silently wrong overlap. So the fused
      !! table is a companion that exactly one constructor
      !! (`eri_shell_table`) ever hands out, and only the plain four-centre
      !! loops consume it.
      !!
      !! Fusing was first tried in the packer itself, making the fused form
      !! primary and teaching every reader of `bas` about it. That founders on
      !! the paragraph above -- the readers include libfint drivers that
      !! cannot be taught -- and is why this is a post-pass here instead.
      type(libcint_molecule_t), intent(inout) :: this

      integer :: ish, iview, nview, nprim, extra, off, ptr_s, ptr_p

      ! Count first. The two passes ask the same question in the same order,
      ! so they cannot disagree about what a shell is.
      nview = 0
      extra = 0
      ish = 1
      do while (ish <= this%nbas)
         nview = nview + 1
         if (fused_sp_pair(this, ish)) then
            extra = extra + 2*this%bas(LIBCINT_NPRIM_OF, ish)
            ish = ish + 2
         else
            ish = ish + 1
         end if
      end do

      ! No L shell anywhere -- cc-pVDZ, say -- keeps no view: `nbas_sp` stays
      ! zero and every consumer falls back to the split table it always used.
      if (nview == this%nbas) return

      this%nbas_sp = nview
      allocate (this%bas_sp(LIBCINT_BAS_SLOTS, nview))
      allocate (this%shell_offset_sp(nview + 1))
      allocate (this%sp_split_first(nview + 1))
      allocate (this%env_sp(size(this%env) + extra))
      this%bas_sp = 0
      this%env_sp(1:size(this%env)) = this%env
      this%env_sp(size(this%env) + 1:) = 0.0_dp

      off = size(this%env)
      iview = 0
      ish = 1
      do while (ish <= this%nbas)
         iview = iview + 1
         this%sp_split_first(iview) = ish
         ! One s then three p is exactly the split AO order, so the view
         ! shell starts where its first split shell did.
         this%shell_offset_sp(iview) = this%shell_offset(ish)
         if (fused_sp_pair(this, ish)) then
            nprim = this%bas(LIBCINT_NPRIM_OF, ish)
            this%bas_sp(LIBCINT_ATOM_OF, iview) = this%bas(LIBCINT_ATOM_OF, ish)
            ! The higher l present, so every ceiling, stride and Rys order
            ! libfint derives is the one the p sub-block needs; the s
            ! sub-block then runs with more roots than it needs, which is
            ! exact -- Rys is exact to degree 2n-1.
            this%bas_sp(LIBCINT_ANG_OF, iview) = 1
            this%bas_sp(LIBCINT_NPRIM_OF, iview) = nprim
            ! NOT doubled. libfint's convention (src/cint_bas.f90:40) is one
            ! contraction whose coefficient block is nprim s coefficients then
            ! nprim p ones -- the block a plain shell with two contractions
            ! would have, but still one contraction. Writing 2 here gives a
            ! shell of two p contractions: it runs, and produces 139 basis
            ! functions where the molecule has 83.
            this%bas_sp(LIBCINT_NCTR_OF, iview) = 1
            this%bas_sp(LIBCINT_KAPPA_OF, iview) = KAPPA_SP_SHELL
            ! One set of exponents: the s shell's, which `fused_sp_pair` has
            ! checked is bit-for-bit the p shell's too.
            this%bas_sp(LIBCINT_PTR_EXP, iview) = this%bas(LIBCINT_PTR_EXP, ish)
            ! Copied, not renormalised. Each split column already carries its
            ! own l's normalisation -- primitive norm, overlap exponent, and
            ! the norm divided back out all depend on l -- and that per-column
            ! form is exactly what the fused convention wants. Recomputing it
            ! here with the fused shell's ANG_OF would put the s column on a p
            ! normalisation: a wrong basis that still converges.
            ptr_s = this%bas(LIBCINT_PTR_COEFF, ish)
            ptr_p = this%bas(LIBCINT_PTR_COEFF, ish + 1)
            this%bas_sp(LIBCINT_PTR_COEFF, iview) = off
            this%env_sp(off + 1:off + nprim) = this%env(ptr_s + 1:ptr_s + nprim)
            this%env_sp(off + nprim + 1:off + 2*nprim) = this%env(ptr_p + 1:ptr_p + nprim)
            off = off + 2*nprim
            ish = ish + 2
         else
            this%bas_sp(:, iview) = this%bas(:, ish)
            ish = ish + 1
         end if
      end do
      this%sp_split_first(nview + 1) = this%nbas + 1
      this%shell_offset_sp(nview + 1) = this%nao
   end subroutine build_sp_view

   function fused_sp_pair(this, ish) result(is_sp)
      !! Whether split shells `ish` and `ish+1` are the s and p of one L shell
      !!
      !! The reader emits one shell per angular momentum, so a basis entry
      !! with `"angular_momentum": [0, 1]` -- every Pople set's valence shells
      !! -- arrives as an s shell and a p shell, consecutive and sharing one
      !! set of exponents. That adjacency is what is detected here, on the
      !! packed table rather than the reader's, so a molecule assembled any
      !! other way (an atom subset, a test fixture) gets the same answer from
      !! the same bytes.
      !!
      !! Only the k=1 case, one s column and one p column, which is what every
      !! Pople set is. A generally contracted sp block arrives as an s and a p
      !! shell of several columns each -- `contraction_group_size` has already
      !! merged same-l columns -- and is left to the split path rather than
      !! guessed at.
      type(libcint_molecule_t), intent(in) :: this
      integer, intent(in) :: ish
      logical :: is_sp

      integer :: nprim, ptr_s, ptr_p

      is_sp = .false.
      if (ish + 1 > this%nbas) return
      if (this%bas(LIBCINT_ATOM_OF, ish) /= this%bas(LIBCINT_ATOM_OF, ish + 1)) return
      if (this%bas(LIBCINT_ANG_OF, ish) /= 0) return
      if (this%bas(LIBCINT_ANG_OF, ish + 1) /= 1) return
      if (this%bas(LIBCINT_NCTR_OF, ish) /= 1) return
      if (this%bas(LIBCINT_NCTR_OF, ish + 1) /= 1) return
      nprim = this%bas(LIBCINT_NPRIM_OF, ish)
      if (this%bas(LIBCINT_NPRIM_OF, ish + 1) /= nprim) return
      ptr_s = this%bas(LIBCINT_PTR_EXP, ish)
      ptr_p = this%bas(LIBCINT_PTR_EXP, ish + 1)
      ! Bit-for-bit, like the contraction grouping in the packer: both blocks
      ! were copied from the same parsed strings, so anything short of
      ! equality means different primitives that must not share a shell.
      if (any(this%env(ptr_s + 1:ptr_s + nprim) /= this%env(ptr_p + 1:ptr_p + nprim))) return
      is_sp = .true.
   end function fused_sp_pair
#else
   subroutine build_sp_view(this)
      !! Without libfint there is nothing to fuse: libcint would read
      !! `KAPPA_SP_SHELL` as a real relativistic kappa and act on it. The view
      !! stays empty and every consumer keeps the split table.
      type(libcint_molecule_t), intent(inout) :: this

      this%nbas_sp = 0
   end subroutine build_sp_view
#endif

   subroutine eri_shell_table(mol, tab)
      !! The shell set to hand libfint's `int2e`: fused where the view exists
      !!
      !! Everything else -- one-electron integrals, three-centre fitting
      !! integrals, every derivative, AO evaluation on a grid -- keeps reading
      !! `mol%bas`, because `int2e_sph`/`int2e_cart` are the only drivers that
      !! understand an L shell. Routing the choice through one constructor is
      !! what keeps that boundary in one place.
      type(libcint_molecule_t), intent(in) :: mol
      type(eri_shell_table_t), intent(out) :: tab

      integer :: ish

      if (mol%nbas_sp > 0) then
         tab%nbas = mol%nbas_sp
         tab%bas = mol%bas_sp
         tab%env = mol%env_sp
         tab%offs = mol%shell_offset_sp
      else
         tab%nbas = mol%nbas
         tab%bas = mol%bas
         tab%env = mol%env
         tab%offs = mol%shell_offset
      end if
      allocate (tab%dims(tab%nbas))
      tab%block_max = 1
      do ish = 1, tab%nbas
         ! Offset differences, not `cgto` calls: the offsets were built from
         ! the dimensions, and libfint already answers four for a fused shell,
         ! so the two agree by construction.
         tab%dims(ish) = tab%offs(ish + 1) - tab%offs(ish)
         tab%block_max = max(tab%block_max, tab%dims(ish))
      end do
   end subroutine eri_shell_table

   subroutine eri_schwarz_collapse(mol, bounds, collapsed)
      !! Schwarz bounds re-blocked onto the eri view's shells
      !!
      !! `schwarz_bounds` works per split shell and stays that way: the
      !! gradient's quartet loop indexes its result directly and runs over
      !! split shells, because its integrals are `int2e_ip1`, which the view
      !! must never reach. A fused shell's bound is the largest of its
      !! sub-shells': the Schwarz inequality only needs the diagonal elements
      !! (mn|mn), every one of which lives in some split sub-block, so the
      !! maximum over sub-blocks bounds the fused quartet at least as tightly
      !! as the split bounds bounded the split ones.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: bounds(:, :)
      real(dp), allocatable, intent(out) :: collapsed(:, :)

      integer :: i, j, a0, a1, b0, b1

      if (mol%nbas_sp <= 0) then
         collapsed = bounds
         return
      end if

      allocate (collapsed(mol%nbas_sp, mol%nbas_sp))
      do j = 1, mol%nbas_sp
         b0 = mol%sp_split_first(j)
         b1 = mol%sp_split_first(j + 1) - 1
         do i = 1, mol%nbas_sp
            a0 = mol%sp_split_first(i)
            a1 = mol%sp_split_first(i + 1) - 1
            collapsed(i, j) = maxval(bounds(a0:a1, b0:b1))
         end do
      end do
   end subroutine eri_schwarz_collapse

   subroutine molecule_overlap(this, s)
      !! S, shell pair by shell pair
      class(libcint_molecule_t), intent(in) :: this
      real(dp), allocatable, intent(out) :: s(:, :)

      call one_electron(this, s, 1)
   end subroutine molecule_overlap

   subroutine molecule_kinetic(this, t)
      !! T on its own
      !!
      !! Exchange repulsion between two fragments needs the kinetic energy separately;
      !! the core Hamiltonian folds the nuclear attraction in with it and cannot serve.
      class(libcint_molecule_t), intent(in) :: this
      real(dp), allocatable, intent(out) :: t(:, :)

      call one_electron(this, t, 2)
   end subroutine molecule_kinetic

   subroutine molecule_core_hamiltonian(this, h)
      !! H = T + V, the one-electron part of the Fock matrix
      class(libcint_molecule_t), intent(in) :: this
      real(dp), allocatable, intent(out) :: h(:, :)

      real(dp), allocatable :: v(:, :)

      call one_electron(this, h, 2)   ! kinetic
      call one_electron(this, v, 3)   ! nuclear attraction
      h = h + v
   end subroutine molecule_core_hamiltonian

   subroutine one_electron(this, matrix, which)
      !! Any of the three one-electron matrices, by the same loop
      !!
      !! One routine rather than three copies: the only thing that differs is
      !! which libcint entry point is called, and three copies of a shell-pair
      !! loop is three places for an offset to be wrong in.
      type(libcint_molecule_t), intent(in) :: this
      real(dp), allocatable, intent(out) :: matrix(:, :)
      integer, intent(in) :: which   !! 1 overlap, 2 kinetic, 3 nuclear

      real(dp), allocatable :: buf(:)
      integer :: shls(2)
      integer :: ish, jsh, di, dj, i, j, ret, io, jo

      allocate (matrix(this%nao, this%nao))
      matrix = 0.0_dp
      ! A shell *pair* block is di*dj, so the bound is the square of the
      ! largest shell -- not the largest shell. With s functions only,
      ! di*dj is 1 and the wrong bound is indistinguishable from the right
      ! one; the first p shell writes nine doubles into three and corrupts
      ! the heap, which surfaces as a free() failure somewhere else later.
      allocate (buf(max_block(this)**2))

      do ish = 1, this%nbas
         di = shell_dim(this%cartesian, ish - 1, this%bas)
         io = this%shell_offset(ish)
         do jsh = 1, this%nbas
            dj = shell_dim(this%cartesian, jsh - 1, this%bas)
            jo = this%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]

            if (this%cartesian) then
               select case (which)
               case (1)
                  ret = libcint_1e_ovlp_cart(buf, shls, this%atm, this%natm, &
                                             this%bas, this%nbas, this%env)
               case (2)
                  ret = libcint_1e_kin_cart(buf, shls, this%atm, this%natm, &
                                            this%bas, this%nbas, this%env)
               case default
                  ret = libcint_1e_nuc_cart(buf, shls, this%atm, this%natm, &
                                            this%bas, this%nbas, this%env)
               end select
            else
               select case (which)
               case (1)
                  ret = libcint_1e_ovlp_sph(buf, shls, this%atm, this%natm, &
                                            this%bas, this%nbas, this%env)
               case (2)
                  ret = libcint_1e_kin_sph(buf, shls, this%atm, this%natm, &
                                           this%bas, this%nbas, this%env)
               case default
                  ret = libcint_1e_nuc_sph(buf, shls, this%atm, this%natm, &
                                           this%bas, this%nbas, this%env)
               end select
            end if
            if (ret == 0) cycle   ! screened away, leave the block zero

            ! libcint fills the block in column-major order, which is what
            ! Fortran wants -- element (i,j) of the block is buf(i + (j-1)*di).
            do j = 1, dj
               do i = 1, di
                  matrix(io + i, jo + j) = buf(i + (j - 1)*di)
               end do
            end do
         end do
      end do
      deallocate (buf)
   end subroutine one_electron

   subroutine build_df_tensor(orb, aux, b, error, omega)
      !! B(mu nu, P) = sum_Q (mu nu | Q) [(P|Q)^(-1/2)]_QP
      !!
      !! The fitted J and K are contractions of this one tensor, which is why
      !! it is formed once and kept rather than recomputed per iteration. It
      !! is (nao^2, naux) -- large, but n^2 rather than the n^4 of the exact
      !! integrals, and that is the whole point of fitting.
      !!
      !! The metric is inverted through its eigendecomposition rather than a
      !! Cholesky, and modes below a threshold are dropped. A JKFIT auxiliary
      !! basis is close to linearly dependent by construction, and a Cholesky
      !! of a near-singular metric fails outright where this degrades.
      type(libcint_molecule_t), intent(in) :: orb   !! Orbital basis
      type(libcint_molecule_t), intent(in) :: aux   !! Auxiliary basis, same atoms
      real(dp), allocatable, intent(out) :: b(:, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: omega
         !! Fit the *attenuated* kernel `erf(omega r)/r` instead of `1/r`, which
         !! is what a range-separated functional's long-range exchange needs.
         !!
         !! Both halves are attenuated, not just the three-centre one. The
         !! fitted approximation is
         !!
         !!     (mu nu|lam sig)_omega ~ sum_PQ (mu nu|P)_omega [J_omega]^-1_PQ (Q|lam sig)_omega
         !!
         !! so the metric has to be `(P|Q)_omega` as well. Mixing an attenuated
         !! three-centre tensor with a full-range metric is not a worse fit of
         !! the same thing -- it fits nothing, and the error does not shrink as
         !! the auxiliary basis grows, which is how it would be noticed.

      real(dp), parameter :: NULL_THRESHOLD = 1.0e-10_dp
      real(dp), allocatable :: metric(:, :), vectors(:, :), values(:), half(:, :)
      real(dp), allocatable :: three(:, :)
      integer :: naux, nao, i, j, kept, info

      ! (mu nu | P) has orbital shells on two centres and auxiliary shells on
      ! the third, in one libcint call -- and that call is spherical or
      ! Cartesian for all three at once. A Cartesian orbital basis fitted with
      ! a spherical auxiliary one cannot be expressed, and guessing which side
      ! to honour would silently fit in a basis neither deck asked for.
      type(timing_report_t) :: clk

      if (orb%cartesian .neqv. aux%cartesian) then
         call error%set(ERROR_VALIDATION, "density fitting: the orbital basis is "// &
                        angular_form_name(orb%cartesian)//" and the auxiliary basis is "// &
                        angular_form_name(aux%cartesian)//". libcint builds all three "// &
                        "centres of a fitting integral in one form, so the two must "// &
                        "agree; choose an auxiliary basis of the same kind.")
         return
      end if

      nao = orb%nao
      naux = aux%nao

      call clk%start()
      if (present(omega)) then
         call two_centre(aux, metric, omega=omega)
      else
         call two_centre(aux, metric)
      end if
      call clk%lap("2-centre metric")
      if (present(omega)) then
         call three_centre(orb, aux, three, omega=omega)
      else
         call three_centre(orb, aux, three)
      end if
      call clk%lap("3-centre integrals")

      call metric_inverse_sqrt(metric, half, error)
      call clk%lap("metric^-1/2")
      if (error%has_error()) return

      allocate (b(nao*nao, naux))
      call pic_gemm(three, half, b)
      call clk%lap("B = (mn|Q) J^-1/2")
      call clk%finish()
      call clk%report("density fitting")
   end subroutine build_df_tensor

   subroutine metric_inverse_sqrt(metric, half, error)
      !! J^(-1/2) = U s^(-1/2) U^T over the modes that survive the threshold
      !!
      !! Through the eigendecomposition rather than a Cholesky. A JKFIT or
      !! RIFIT set is close to linearly dependent by construction, and a
      !! Cholesky of a near-singular metric fails outright where this degrades.
      real(dp), intent(in) :: metric(:, :)
      real(dp), allocatable, intent(out) :: half(:, :)
      type(error_t), intent(inout) :: error

      real(dp), parameter :: NULL_THRESHOLD = 1.0e-10_dp
      real(dp), allocatable :: vectors(:, :), values(:), scaled(:, :)
      integer :: naux, i, kept, info

      naux = size(metric, 1)
      allocate (vectors(naux, naux), values(naux))
      vectors = metric
      call pic_syev(vectors, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "density fitting: the metric would not diagonalise")
         return
      end if

      kept = count(values > NULL_THRESHOLD)
      if (kept == 0) then
         call error%set(ERROR_VALIDATION, "density fitting: the auxiliary metric is singular")
         return
      end if

      ! U s^(-1/2) U^T as one gemm rather than an outer-product loop. The loop
      ! this replaces was naux^2 strided vector updates -- the same arithmetic,
      ! but memory-bound and invisible to BLAS, and it showed up as the largest
      ! non-library cost in a profile of RI-MP2. Dropped modes are zeroed in the
      ! scaled copy, which is how they stay out of the product.
      allocate (scaled(naux, naux))
      do i = 1, naux
         if (values(i) > NULL_THRESHOLD) then
            scaled(:, i) = vectors(:, i)/sqrt(values(i))
         else
            scaled(:, i) = 0.0_dp
         end if
      end do

      allocate (half(naux, naux))
      call pic_gemm(scaled, vectors, half, transb="T")
      deallocate (scaled)
   end subroutine metric_inverse_sqrt

   subroutine build_df_mo_tensor(orb, aux, c_occ, c_vir, bia, error)
      !! B^P_ia directly, transforming to MO before fitting rather than after
      !!
      !! `build_df_tensor` fits first and hands back B(mu nu, P), which is what
      !! a Fock build wants because it contracts against a density in the AO
      !! basis. A correlation treatment does not: it only ever reads the
      !! occupied-virtual block, and fitting the whole AO pair space to get
      !! there is most of the work thrown away.
      !!
      !! The metric contraction is where it shows. Applied to the AO tensor it
      !! is nao^2 by naux by naux; applied after the transform it is
      !! n_occ*n_vir by naux by naux. On a water hexamer in cc-pVDZ that is
      !! 20736 pair functions against 2736, so the same gemm costs 7.6 times
      !! less for the same answer.
      !!
      !! The result is laid out (n_vir, naux, n_occ) so that one occupied
      !! orbital's block is contiguous, which is the shape the energy step
      !! hands to BLAS.
      type(libcint_molecule_t), intent(in) :: orb, aux
      real(dp), intent(in) :: c_occ(:, :)   !! (nao, n_occ), correlated occupied only
      real(dp), intent(in) :: c_vir(:, :)   !! (nao, n_vir)
      real(dp), allocatable, intent(out) :: bia(:, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: fitted(:, :)
      integer :: naux, n_o, n_v, p_index, i

      n_o = size(c_occ, 2)
      n_v = size(c_vir, 2)

      call build_df_mo_block(orb, aux, c_occ, c_vir, fitted, error)
      if (error%has_error()) return
      naux = size(fitted, 2)

      allocate (bia(n_v, naux, n_o))
      do p_index = 1, naux
         do i = 1, n_o
            ! `fitted` is flattened (i, a), so one occupied orbital's virtuals
            ! sit at stride n_o.
            bia(:, p_index, i) = fitted(i::n_o, p_index)
         end do
      end do
      deallocate (fitted)
   end subroutine build_df_mo_tensor

   subroutine build_df_mo_block(orb, aux, c_left, c_right, b, error)
      !! B^P_pq for any two coefficient blocks, laid out (pq, P)
      !!
      !! The general form of `build_df_mo_tensor`, which assumed
      !! occupied-virtual because MP2 wants nothing else. Coupled cluster wants
      !! every block, so the choice moves to the caller -- the same move
      !! `transform_ovov` made when it became `transform_block`, and for the same
      !! reason: one routine cannot transpose a half differently in six places.
      !!
      !! The compound index runs left-fastest, `pq = (q-1) n_left + p`, which is
      !! the layout that makes `B B^T` the fitted `(pq|rs)` with no repacking.
      type(libcint_molecule_t), intent(in) :: orb, aux
      real(dp), intent(in) :: c_left(:, :), c_right(:, :)   !! (nao, n_left), (nao, n_right)
      real(dp), allocatable, intent(out) :: b(:, :)         !! (n_left*n_right, naux)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: metric(:, :), half(:, :), three(:, :)
      real(dp), allocatable :: ovl(:, :), bp(:, :), tmp(:, :), full(:, :)
      integer :: nao, naux, n_l, n_r, p_index

      if (orb%cartesian .neqv. aux%cartesian) then
         call error%set(ERROR_VALIDATION, "density fitting: the orbital basis is "// &
                        angular_form_name(orb%cartesian)//" and the auxiliary basis is "// &
                        angular_form_name(aux%cartesian)//". libcint builds all three "// &
                        "centres of a fitting integral in one form, so the two must "// &
                        "agree; choose an auxiliary basis of the same kind.")
         return
      end if

      nao = orb%nao
      naux = aux%nao
      n_l = size(c_left, 2)
      n_r = size(c_right, 2)

      call two_centre(aux, metric)
      call metric_inverse_sqrt(metric, half, error)
      if (error%has_error()) return
      deallocate (metric)

      call three_centre(orb, aux, three)

      ! (mu nu|P) -> (pq|P), one auxiliary function at a time. Transforming
      ! before fitting rather than after, which is the ordering PySCF uses and
      ! the one that keeps the metric contraction off the AO pair space.
      allocate (ovl(n_l*n_r, naux))
      allocate (bp(nao, nao), tmp(n_l, nao), full(n_l, n_r))
      do p_index = 1, naux
         bp = reshape(three(:, p_index), [nao, nao])
         call pic_gemm(c_left, bp, tmp, transa="T")
         call pic_gemm(tmp, c_right, full)
         ovl(:, p_index) = reshape(full, [n_l*n_r])
      end do
      deallocate (bp, tmp, full, three)

      allocate (b(n_l*n_r, naux))
      call pic_gemm(ovl, half, b)
      deallocate (ovl, half)
   end subroutine build_df_mo_block

   subroutine two_centre(aux, metric, omega)
      !! (P|Q) over the auxiliary basis
      !!
      !! `omega` switches the kernel to `erf(omega r)/r`, which is what a
      !! range-separated functional's long-range exchange is fitted against.
      !! libcint takes it through `env` rather than through a different entry
      !! point, so the attenuated metric is the same loop over a modified copy.
      !!
      !! **The slot is one-based here and zero-based in libcint's headers**, so
      !! it is `LIBCINT_PTR_RANGE_OMEGA + 1`. Getting it wrong is silent: the
      !! neighbouring slot is ignored by a two-centre integral, so the
      !! "attenuated" metric would come back full-range and the long- and
      !! short-range pieces would sum to unscaled exchange.
      type(libcint_molecule_t), intent(in) :: aux
      real(dp), allocatable, intent(out) :: metric(:, :)
      real(dp), intent(in), optional :: omega

      real(dp), allocatable :: buf(:), env_local(:)
      integer :: shls(2)
      integer :: ish, jsh, di, dj, i, j, io, jo, ret

      allocate (metric(aux%nao, aux%nao))
      metric = 0.0_dp
      allocate (buf(max_block(aux)**2))

      ! A copy, because `aux%env` is shared and an attenuated build must not
      ! leave the omega set behind for the next caller.
      env_local = aux%env
      if (present(omega)) env_local(LIBCINT_PTR_RANGE_OMEGA + 1) = omega

      do ish = 1, aux%nbas
         di = shell_dim(aux%cartesian, ish - 1, aux%bas)
         io = aux%shell_offset(ish)
         do jsh = 1, aux%nbas
            dj = shell_dim(aux%cartesian, jsh - 1, aux%bas)
            jo = aux%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]
            if (aux%cartesian) then
               ret = libcint_2c2e_cart(buf, shls, aux%atm, aux%natm, aux%bas, aux%nbas, &
                                       env_local)
            else
               ret = libcint_2c2e_sph(buf, shls, aux%atm, aux%natm, aux%bas, aux%nbas, &
                                      env_local)
            end if
            if (ret == 0) cycle
            do j = 1, dj
               do i = 1, di
                  metric(io + i, jo + j) = buf(i + (j - 1)*di)
               end do
            end do
         end do
      end do
   end subroutine two_centre

   subroutine mixed_basis_overlap(orb, aux, overlap, error)
      !! The overlap between two different bases on the same atoms
      !!
      !!     S_mixed(mu, k) = < chi_mu | phi_k >
      !!
      !! with `chi` the orbital basis the SCF ran in and `phi` a second basis --
      !! for the quasi-atomic construction, the free-atom minimal basis the
      !! molecular orbitals are projected onto. Nothing in the quasi-atomic
      !! papers is computable without it: every space in the construction is
      !! defined by an SVD of exactly this matrix.
      !!
      !! It reuses `build_df_shell_table` rather than concatenating the two
      !! bases again. That routine exists because libcint addresses every centre
      !! by index into one shell table, and a second construction of the same
      !! table would agree with the first until one of them was edited.
      !!
      !! The dummy shell that table appends for three-centre calls is simply not
      !! referenced here; a one-electron call names two centres and stops.
      type(libcint_molecule_t), intent(in) :: orb   !! Orbital basis
      type(libcint_molecule_t), intent(in) :: aux   !! Second basis, same atoms
      real(dp), allocatable, intent(out) :: overlap(:, :)  !! (orb%nao, aux%nao)
      type(error_t), intent(inout) :: error

      integer, allocatable :: bas(:, :)
      real(dp), allocatable :: env(:), buf(:)
      integer :: shls(2)
      integer :: nbas_orb, dummy, ish, jsh, di, dj, i, j, io, jo, ret

      if (error%has_error()) return

      ! One angular convention runs the whole call, exactly as it does for the
      ! three-centre fit. Mixing them is not a thing libcint can be asked for:
      ! the entry point is chosen once, so a spherical orbital basis and a
      ! Cartesian auxiliary would silently be integrated as though both were
      ! whichever the branch below picked.
      if (orb%cartesian .neqv. aux%cartesian) then
         call error%set(ERROR_VALIDATION, "mixed-basis overlap: the two bases must "// &
                        "use the same angular form. libcint takes one convention "// &
                        "for the whole call, so a mismatch would integrate one of "// &
                        "them in a form it was not built in.")
         return
      end if

      if (orb%natm /= aux%natm) then
         call error%set(ERROR_VALIDATION, "mixed-basis overlap: the two bases are "// &
                        "on different numbers of atoms, so they do not describe "// &
                        "the same molecule.")
         return
      end if

      call build_df_shell_table(orb, aux, bas, env, dummy)
      nbas_orb = orb%nbas

      allocate (overlap(orb%nao, aux%nao))
      overlap = 0.0_dp
      allocate (buf(max(max_block(orb), max_block(aux))**2))

      do ish = 1, nbas_orb
         di = shell_dim(orb%cartesian, ish - 1, bas)
         io = orb%shell_offset(ish)
         do jsh = 1, aux%nbas
            dj = shell_dim(orb%cartesian, nbas_orb + jsh - 1, bas)
            jo = aux%shell_offset(jsh)
            shls = [ish - 1, nbas_orb + jsh - 1]

            if (orb%cartesian) then
               ret = libcint_1e_ovlp_cart(buf, shls, orb%atm, orb%natm, bas, &
                                          nbas_orb + aux%nbas, env)
            else
               ret = libcint_1e_ovlp_sph(buf, shls, orb%atm, orb%natm, bas, &
                                         nbas_orb + aux%nbas, env)
            end if
            if (ret == 0) cycle   ! screened away, leave the block zero

            do j = 1, dj
               do i = 1, di
                  overlap(io + i, jo + j) = buf(i + (j - 1)*di)
               end do
            end do
         end do
      end do

      deallocate (bas, env, buf)
   end subroutine mixed_basis_overlap

   subroutine build_df_shell_table(orb, aux, bas, env, dummy)
      !! Orbital and auxiliary shells in one table, as libcint needs them
      !!
      !! libcint addresses all four centres of a three-centre call by index
      !! into a single shell table, so the two bases have to be concatenated
      !! and the auxiliary env offsets shifted past the orbital env. The fourth
      !! index names a dummy s shell of zero exponent, which is how libcint
      !! spells "only three centres".
      !!
      !! Extracted so that the energy's `three_centre` and the gradient's
      !! derivative equivalents build the same table. Two constructions of this
      !! would agree until one of them was edited, and the failure would be a
      !! fitted quantity differentiated in a basis it was not built in.
      type(libcint_molecule_t), intent(in) :: orb, aux
      integer, allocatable, intent(out) :: bas(:, :)
      real(dp), allocatable, intent(out) :: env(:)
      integer, intent(out) :: dummy   !! 1-based index of the dummy shell

      integer :: nbas_orb, nbas_aux, n_env_orb, ish

      nbas_orb = orb%nbas
      nbas_aux = aux%nbas
      n_env_orb = size(orb%env)

      allocate (bas(LIBCINT_BAS_SLOTS, nbas_orb + nbas_aux + 1))
      allocate (env(n_env_orb + size(aux%env) + 1))
      bas = 0
      env = 0.0_dp
      env(1:n_env_orb) = orb%env
      env(n_env_orb + 1:n_env_orb + size(aux%env)) = aux%env

      bas(:, 1:nbas_orb) = orb%bas
      bas(:, nbas_orb + 1:nbas_orb + nbas_aux) = aux%bas
      do ish = 1, nbas_aux
         bas(LIBCINT_PTR_EXP, nbas_orb + ish) = aux%bas(LIBCINT_PTR_EXP, ish) + n_env_orb
         bas(LIBCINT_PTR_COEFF, nbas_orb + ish) = aux%bas(LIBCINT_PTR_COEFF, ish) + n_env_orb
      end do

      dummy = nbas_orb + nbas_aux + 1
      bas(LIBCINT_ATOM_OF, dummy) = 0
      bas(LIBCINT_ANG_OF, dummy) = 0
      bas(LIBCINT_NPRIM_OF, dummy) = 1
      bas(LIBCINT_NCTR_OF, dummy) = 1
      bas(LIBCINT_PTR_EXP, dummy) = size(env) - 1
      bas(LIBCINT_PTR_COEFF, dummy) = size(env) - 1
      env(size(env)) = 0.0_dp
   end subroutine build_df_shell_table

   subroutine three_centre(orb, aux, three, omega)
      !! (mu nu | P), flattened to (nao*nao, naux)
      !!
      !! The orbital and auxiliary shells are concatenated into one bas array,
      !! because libcint addresses all four centres of a 3c2e call by index
      !! into a single table. The fourth index names a dummy s shell with a
      !! zero exponent, which is how libcint spells "only three centres".
      !!
      !! One convention runs the whole call, `orb%cartesian`, auxiliary shells
      !! included. `build_df_tensor` has already refused the case where the two
      !! bases disagree, so that is not a choice being made here.
      type(libcint_molecule_t), intent(in) :: orb, aux
      real(dp), allocatable, intent(out) :: three(:, :)
      real(dp), intent(in), optional :: omega
         !! Attenuate the kernel to `erf(omega r)/r`, for the long-range
         !! exchange of a range-separated functional. See `two_centre` for why
         !! the slot index is what it is.

      integer, allocatable :: bas(:, :)
      real(dp), allocatable :: env(:), buf(:)
      integer :: nbas_orb, nbas_aux, dummy
      integer :: shls(4)
      integer :: ish, jsh, ksh, di, dj, dk, i, j, k, io, jo, ko, ret, idx
      integer :: npair, ipair
      integer, allocatable :: pair_i(:), pair_j(:)
      type(c_ptr) :: opt

      nbas_orb = orb%nbas
      nbas_aux = aux%nbas
      call build_df_shell_table(orb, aux, bas, env, dummy)
      ! `env` is built fresh above, so this needs no copy of its own.
      if (present(omega)) env(LIBCINT_PTR_RANGE_OMEGA + 1) = omega

      allocate (three(orb%nao*orb%nao, aux%nao))
      three = 0.0_dp

      ! (ij|K) is symmetric in i and j, so only the lower triangle of the
      ! orbital shell pairs is computed and each block is written into both
      ! halves. The loop ran the full square before, which evaluated every
      ! off-diagonal block twice for the same numbers.
      npair = nbas_orb*(nbas_orb + 1)/2
      allocate (pair_i(npair), pair_j(npair))
      ipair = 0
      do ish = 1, nbas_orb
         do jsh = 1, ish
            ipair = ipair + 1
            pair_i(ipair) = ish
            pair_j(ipair) = jsh
         end do
      end do

      ! Built once for the combined orbital+auxiliary shell set and read by
      ! every thread. Without it libcint redoes the per-shell-pair setup on each
      ! call -- `CINTinit_int3c2e_EnvVars` and `CINTset_pairdata` show up in a
      ! profile of this loop, which is that work repeating.
      opt = c_null_ptr
      if (orb%cartesian) then
         call libcint_3c2e_cart_optimizer(opt, orb%atm, orb%natm, bas, &
                                          nbas_orb + nbas_aux + 1, env)
      else
         call libcint_3c2e_sph_optimizer(opt, orb%atm, orb%natm, bas, &
                                         nbas_orb + nbas_aux + 1, env)
      end if

      ! No accumulator to reduce, unlike the Fock build: a shell pair owns its
      ! block of `three` outright, so no two iterations write the same element
      ! and the threads need nothing shared but a private buffer each. `opt` is
      ! shared and read-only, the same arrangement the Fock build uses.
      !$omp parallel default(none) &
      !$omp    shared(orb, aux, bas, env, three, nbas_orb, nbas_aux, dummy, npair, pair_i, pair_j, opt) &
      !$omp    private(ipair, ish, jsh, ksh, di, dj, dk, io, jo, ko, i, j, k, shls, ret, idx, buf)
      allocate (buf(max_block(orb)**2*max_block(aux)))

      !$omp do schedule(dynamic)
      do ipair = 1, npair
         ish = pair_i(ipair)
         jsh = pair_j(ipair)
         di = shell_dim(orb%cartesian, ish - 1, bas)
         io = orb%shell_offset(ish)
         dj = shell_dim(orb%cartesian, jsh - 1, bas)
         jo = orb%shell_offset(jsh)

         do ksh = 1, nbas_aux
            dk = shell_dim(orb%cartesian, nbas_orb + ksh - 1, bas)
            ko = aux%shell_offset(ksh)
            shls = [ish - 1, jsh - 1, nbas_orb + ksh - 1, dummy - 1]
            if (orb%cartesian) then
               ret = libcint_3c2e_cart(buf, shls, orb%atm, orb%natm, bas, &
                                       nbas_orb + nbas_aux + 1, env, opt)
            else
               ret = libcint_3c2e_sph(buf, shls, orb%atm, orb%natm, bas, &
                                      nbas_orb + nbas_aux + 1, env, opt)
            end if
            if (ret == 0) cycle
            do k = 1, dk
               do j = 1, dj
                  do i = 1, di
                     idx = i + (j - 1)*di + (k - 1)*di*dj
                     three((jo + j - 1)*orb%nao + io + i, ko + k) = buf(idx)
                     ! The mirrored element. On a diagonal shell pair this
                     ! writes the partner entry of the same block, which the
                     ! loop also visits with i and j exchanged and the same
                     ! value, so it is a repeat rather than a conflict.
                     three((io + i - 1)*orb%nao + jo + j, ko + k) = buf(idx)
                  end do
               end do
            end do
         end do
      end do
      !$omp end do

      deallocate (buf)
      !$omp end parallel

      call libcint_del_optimizer(opt)
      deallocate (pair_i, pair_j)
   end subroutine three_centre

   subroutine molecule_eris(this, eri)
      !! Every two-electron integral, in core, as (ij|kl)
      !!
      !! Stored as a full n^4 array. The eightfold symmetry decides which
      !! integrals to *evaluate*, not how they are stored -- the array stays
      !! dense, so every consumer indexes it directly and needs no index map.
      !! That keeps what made this readable as a reference while removing the
      !! eight times too many calls into libcint that filling it naively cost.
      !!
      !! Scattering all eight positions is safe here precisely because this
      !! assigns rather than accumulates: when shells coincide some permutations
      !! name the same element, and writing the same value twice is harmless. An
      !! accumulating version would have to know which had collapsed, which is
      !! the bookkeeping this routine was written to avoid.
      class(libcint_molecule_t), intent(in) :: this
      real(dp), allocatable, intent(out) :: eri(:, :, :, :)

      real(dp), allocatable :: buf(:)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl, lsh_max
      ! `ksh_local` was a BLOCK-scoped declaration inside the parallel region.
      ! nvfortran cannot compile a BLOCK construct in the scope of a parallel
      ! directive (NVFORTRAN-S-1219), and the construct existed only to
      ! declare this one integer, so it is declared here and privatised with
      ! the rest instead.
      integer :: ksh_local
      integer :: shls(4)
      integer :: i, j, k, l, io, jo, ko, lo, ret, idx
      integer :: p, q, r, t
      integer :: npair, ipair
      integer, allocatable :: pair_i(:), pair_j(:)
      real(dp) :: value
      type(c_ptr) :: opt
      type(eri_shell_table_t) :: tab

      allocate (eri(this%nao, this%nao, this%nao, this%nao))
      eri = 0.0_dp

      ! The fused-sp view where it exists, the split shells otherwise. The
      ! tensor is AO-indexed and the view's AO order is the split order, so
      ! nothing past the shell loops changes.
      call eri_shell_table(this, tab)

      opt = c_null_ptr
      call two_electron_optimizer(this%cartesian, opt, this%atm, this%natm, tab%bas, &
                                  tab%nbas, tab%env)

      ! Threaded over bra pairs, exactly as the direct Fock build is, and for
      ! once with nothing to reduce: a canonical quartet owns the eight
      ! positions it scatters into and no other canonical quartet names any of
      ! them, so the threads write disjoint elements of a shared array.
      npair = tab%nbas*(tab%nbas + 1)/2
      allocate (pair_i(npair), pair_j(npair))
      ipair = 0
      do ish = 1, tab%nbas
         do jsh = 1, ish
            ipair = ipair + 1
            pair_i(ipair) = ish
            pair_j(ipair) = jsh
         end do
      end do

      !$omp parallel default(none) &
      !$omp    shared(this, tab, eri, opt, npair, pair_i, pair_j) &
      !$omp    private(ipair, ish, jsh, ksh, lsh, lsh_max, di, dj, dk, dl, io, jo, ko, lo, &
      !$omp            i, j, k, l, p, q, r, t, shls, ret, idx, value, buf, ksh_local)
      allocate (buf(tab%block_max**4))
      !$omp do schedule(dynamic)
      do ipair = 1, npair
         ish = pair_i(ipair)
         jsh = pair_j(ipair)
         di = tab%dims(ish)
         io = tab%offs(ish)
         dj = tab%dims(jsh)
         jo = tab%offs(jsh)
         do ksh_local = 1, ish
            ksh = ksh_local
            dk = tab%dims(ksh)
            ko = tab%offs(ksh)
            lsh_max = ksh
            if (ksh == ish) lsh_max = jsh
            do lsh = 1, lsh_max
               dl = tab%dims(lsh)
               lo = tab%offs(lsh)
               shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]
               ret = two_electron_block(this%cartesian, buf, shls, this%atm, &
                                        this%natm, tab%bas, tab%nbas, tab%env, opt)
               if (ret == 0) cycle

               do l = 1, dl
                  t = lo + l
                  do k = 1, dk
                     r = ko + k
                     do j = 1, dj
                        q = jo + j
                        do i = 1, di
                           p = io + i
                           idx = i + (j - 1)*di + (k - 1)*di*dj + (l - 1)*di*dj*dk
                           value = buf(idx)
                           eri(p, q, r, t) = value
                           eri(q, p, r, t) = value
                           eri(p, q, t, r) = value
                           eri(q, p, t, r) = value
                           eri(r, t, p, q) = value
                           eri(t, r, p, q) = value
                           eri(r, t, q, p) = value
                           eri(t, r, q, p) = value
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end do
      deallocate (buf)
      !$omp end parallel

      call libcint_del_optimizer(opt)
      deallocate (pair_i, pair_j)
   end subroutine molecule_eris

   subroutine ket_transformed_pairs(this, c_vir, c_occ, pair_ov, pair_oo)
      !! `(mu nu | b j)` and `(mu nu | i j)` in one pass, without the AO tensor
      !!
      !! The four-index transform behind an exact response Hessian starts by
      !! contracting the ket pair. Doing that from a stored `eri` costs `n_ao^4`
      !! to hold -- 16.5 GB at 217 orbitals, which is what stops a tripeptide
      !! using the transform at all -- and the tensor is read once and dropped.
      !! Contracting as the quartets come out needs only what is being built.
      !!
      !! **Threaded over bra pairs so the accumulation needs no reduction.**
      !! `molecule_eris` can scatter all eight permutations of a canonical quartet
      !! because it *assigns*: no two canonical quartets name the same element.
      !! A contraction accumulates, so that argument does not carry over. Giving
      !! each thread a bra pair makes it the sole writer of those rows, and the
      !! price is the bra-ket permutation: every `(mu nu | lambda sigma)` is
      !! computed once here where the stored build computes it for eight
      !! positions, so roughly twice the integral work of a symmetric pass. The
      !! `(lambda sigma)` symmetry is still used.
      !!
      !! Both blocks come out of one pass because they share every integral and
      !! differ only in what the ket is contracted with.
      class(libcint_molecule_t), intent(in) :: this
      real(dp), intent(in) :: c_vir(:, :), c_occ(:, :)
      real(dp), intent(out) :: pair_ov(:, :)   !! `(n_ao^2, n_vir*n_occ)`, `b` fastest
      real(dp), intent(out) :: pair_oo(:, :)   !! `(n_ao^2, n_occ*n_occ)`, `i` fastest

      real(dp), allocatable :: buf(:), ket(:, :, :), tmp(:, :), slab(:, :)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo
      integer :: i, j, k, l, ret, idx, m, n_ao, n_vir, n_occ, row
      integer :: npair, ipair, nket, iket
      integer, allocatable :: pair_i(:), pair_j(:)
      real(dp) :: value
      type(c_ptr) :: opt
      type(eri_shell_table_t) :: tab

      n_ao = this%nao
      n_vir = size(c_vir, 2)
      n_occ = size(c_occ, 2)
      pair_ov = 0.0_dp
      pair_oo = 0.0_dp

      ! The fused-sp view where it exists; AO-indexed output, so only the
      ! shell loops feel it.
      call eri_shell_table(this, tab)

      opt = c_null_ptr
      call two_electron_optimizer(this%cartesian, opt, this%atm, this%natm, tab%bas, &
                                  tab%nbas, tab%env)

      npair = tab%nbas*tab%nbas
      nket = tab%nbas*(tab%nbas + 1)/2
      allocate (pair_i(npair), pair_j(npair))
      ipair = 0
      do ish = 1, tab%nbas
         do jsh = 1, tab%nbas
            ipair = ipair + 1
            pair_i(ipair) = ish
            pair_j(ipair) = jsh
         end do
      end do

      !$omp parallel default(none) &
      !$omp    shared(this, tab, c_vir, c_occ, pair_ov, pair_oo, opt, npair, pair_i, pair_j, &
      !$omp           n_ao, n_vir, n_occ) &
      !$omp    private(ipair, ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo, &
      !$omp            i, j, k, l, m, ret, idx, value, buf, ket, tmp, slab, row, iket)
      allocate (buf(tab%block_max**4))
      allocate (ket(n_ao, n_ao, tab%block_max**2))
      allocate (tmp(n_ao, max(n_vir, n_occ)), slab(max(n_vir, n_occ), max(n_vir, n_occ)))
      !$omp do schedule(dynamic)
      do ipair = 1, npair
         ish = pair_i(ipair)
         jsh = pair_j(ipair)
         di = tab%dims(ish)
         io = tab%offs(ish)
         dj = tab%dims(jsh)
         jo = tab%offs(jsh)

         ! Every ket pair against this bra pair, filled both ways round so the
         ! contraction below can be a plain matrix product.
         ket(:, :, 1:di*dj) = 0.0_dp
         do ksh = 1, tab%nbas
            dk = tab%dims(ksh)
            ko = tab%offs(ksh)
            do lsh = 1, ksh
               dl = tab%dims(lsh)
               lo = tab%offs(lsh)
               ret = two_electron_block(this%cartesian, buf, &
                                        [ish - 1, jsh - 1, ksh - 1, lsh - 1], this%atm, &
                                        this%natm, tab%bas, tab%nbas, tab%env, opt)
               if (ret == 0) cycle
               do l = 1, dl
                  do k = 1, dk
                     do j = 1, dj
                        do i = 1, di
                           idx = i + (j - 1)*di + (k - 1)*di*dj + (l - 1)*di*dj*dk
                           value = buf(idx)
                           m = i + (j - 1)*di
                           ket(ko + k, lo + l, m) = value
                           ket(lo + l, ko + k, m) = value
                        end do
                     end do
                  end do
               end do
            end do
         end do

         ! Contract the ket pair. The rows written are this bra pair's own, so
         ! no two threads touch the same one.
         do j = 1, dj
            do i = 1, di
               m = i + (j - 1)*di
               row = (io + i) + n_ao*(jo + j - 1)

               call pic_gemm(ket(:, :, m), c_occ, tmp(:, 1:n_occ))
               call pic_gemm(c_vir, tmp(:, 1:n_occ), slab(1:n_vir, 1:n_occ), transa="T")
               pair_ov(row, :) = reshape(slab(1:n_vir, 1:n_occ), [n_vir*n_occ])

               call pic_gemm(c_occ, tmp(:, 1:n_occ), slab(1:n_occ, 1:n_occ), transa="T")
               pair_oo(row, :) = reshape(slab(1:n_occ, 1:n_occ), [n_occ*n_occ])
            end do
         end do
      end do
      !$omp end do
      deallocate (buf, ket, tmp, slab)
      !$omp end parallel

      call libcint_del_optimizer(opt)
      deallocate (pair_i, pair_j)
   end subroutine ket_transformed_pairs

   subroutine molecule_eris_packed(this, eri)
      !! Every unique two-electron integral, as (pq|rs) over packed AO pairs
      !!
      !! The dense n^4 form that `eris` returns costs more to fill than to
      !! compute: eight positions written per integral, across eight different
      !! stride patterns. That is bandwidth rather than arithmetic, and it does
      !! not improve with more cores -- the threaded dense build gets 4.3x on
      !! forty of them.
      !!
      !! Packing the two AO pairs collapses (mu nu) with (nu mu), so of the
      !! eightfold symmetry only the bra-ket swap is left to write out: two
      !! stores per integral instead of eight, into an array a quarter the size.
      !! A water hexamer in cc-pVDZ goes from 3.4 GB to 872 MB.
      !!
      !! `eris` stays. The in-core SCF and `check_direct` read the dense form,
      !! and a reference implementation is worth more readable than packed.
      class(libcint_molecule_t), intent(in) :: this
      real(dp), allocatable, intent(out) :: eri(:, :)

      real(dp), allocatable :: buf(:)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl, lsh_max
      integer :: shls(4)
      integer :: i, j, k, l, io, jo, ko, lo, ret, idx
      integer :: p, q, r, t, pq, rt, n_pair
      integer :: npair_sh, ipair
      integer, allocatable :: pair_i(:), pair_j(:)
      real(dp) :: value
      type(c_ptr) :: opt
      type(eri_shell_table_t) :: tab

      n_pair = this%nao*(this%nao + 1)/2
      allocate (eri(n_pair, n_pair))
      eri = 0.0_dp

      ! The fused-sp view where it exists; AO-indexed output, so only the
      ! shell loops feel it.
      call eri_shell_table(this, tab)

      opt = c_null_ptr
      call two_electron_optimizer(this%cartesian, opt, this%atm, this%natm, tab%bas, &
                                  tab%nbas, tab%env)

      npair_sh = tab%nbas*(tab%nbas + 1)/2
      allocate (pair_i(npair_sh), pair_j(npair_sh))
      ipair = 0
      do ish = 1, tab%nbas
         do jsh = 1, ish
            ipair = ipair + 1
            pair_i(ipair) = ish
            pair_j(ipair) = jsh
         end do
      end do

      !$omp parallel default(none) &
      !$omp    shared(this, tab, eri, opt, npair_sh, pair_i, pair_j) &
      !$omp    private(ipair, ish, jsh, ksh, lsh, lsh_max, di, dj, dk, dl, io, jo, ko, lo, &
      !$omp            i, j, k, l, p, q, r, t, pq, rt, shls, ret, idx, value, buf)
      allocate (buf(tab%block_max**4))
      !$omp do schedule(dynamic)
      do ipair = 1, npair_sh
         ish = pair_i(ipair)
         jsh = pair_j(ipair)
         di = tab%dims(ish)
         io = tab%offs(ish)
         dj = tab%dims(jsh)
         jo = tab%offs(jsh)
         do ksh = 1, ish
            dk = tab%dims(ksh)
            ko = tab%offs(ksh)
            lsh_max = ksh
            if (ksh == ish) lsh_max = jsh
            do lsh = 1, lsh_max
               dl = tab%dims(lsh)
               lo = tab%offs(lsh)
               shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]
               ret = two_electron_block(this%cartesian, buf, shls, this%atm, &
                                        this%natm, tab%bas, tab%nbas, tab%env, opt)
               if (ret == 0) cycle

               do l = 1, dl
                  t = lo + l
                  do k = 1, dk
                     r = ko + k
                     rt = pair_index(r, t)
                     do j = 1, dj
                        q = jo + j
                        do i = 1, di
                           p = io + i
                           idx = i + (j - 1)*di + (k - 1)*di*dj + (l - 1)*di*dj*dk
                           value = buf(idx)
                           pq = pair_index(p, q)
                           eri(pq, rt) = value
                           eri(rt, pq) = value
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end do
      deallocate (buf)
      !$omp end parallel

      call libcint_del_optimizer(opt)
      deallocate (pair_i, pair_j)
   end subroutine molecule_eris_packed

   pure function pair_index(p, q) result(pq)
      !! Where the AO pair (p,q) sits once the two orderings are collapsed
      !!
      !! Lower triangle, one-based: (1,1) is 1, (2,1) is 2, (2,2) is 3. The
      !! caller need not order its arguments.
      integer, intent(in) :: p, q
      integer :: pq

      if (p >= q) then
         pq = p*(p - 1)/2 + q
      else
         pq = q*(q - 1)/2 + p
      end if
   end function pair_index
   pure function molecule_nuclear_repulsion(this) result(energy)
      !! Sum of Z_a Z_b / R_ab over pairs
      class(libcint_molecule_t), intent(in) :: this
      real(dp) :: energy

      integer :: a, b
      real(dp) :: r

      energy = 0.0_dp
      do a = 1, this%natm
         do b = a + 1, this%natm
            r = norm2(this%coords(:, a) - this%coords(:, b))
            if (r > 0.0_dp) energy = energy + this%charges(a)*this%charges(b)/r
         end do
      end do
   end function molecule_nuclear_repulsion

   subroutine molecule_atom_subset(this, iatom, atom_mol, error)
      !! One atom of this molecule, on its own, in the same basis
      !!
      !! Built by copying the parent's shell data rather than re-reading the
      !! basis file. That is the whole point: the free atom then carries
      !! bit-identical exponents and contraction coefficients to the block it
      !! will be dropped back into, so its AO count cannot disagree with that
      !! block's and its orbitals are expressed in exactly those functions. A
      !! re-read would usually agree and would occasionally not -- a general
      !! contraction merged differently, a Cartesian/spherical flag taken from a
      !! different file -- and the symptom would be a guess quietly built in the
      !! wrong basis.
      !!
      !! The coefficients already carry `libcint_gto_norm` from `molecule_build`,
      !! so they are copied untouched. Re-normalising here would apply it twice.
      !!
      !! Placed at the origin. A free atom has no field to be oriented in, and
      !! nothing downstream reads the coordinate except the nuclear repulsion,
      !! which is zero for one atom.
      class(libcint_molecule_t), intent(in) :: this
      integer, intent(in) :: iatom                          !! 1-based
      type(libcint_molecule_t), intent(out) :: atom_mol
      type(error_t), intent(inout) :: error

      integer :: ish, k, nprim, nctr, off, env_size, target_atom, ncoeff

      if (iatom < 1 .or. iatom > this%natm) then
         call error%set(ERROR_VALIDATION, "libcint: atom index outside the molecule")
         return
      end if

      target_atom = iatom - 1   ! libcint counts atoms from zero

      atom_mol%nbas = 0
      env_size = LIBCINT_PTR_ENV_START + 3
      do ish = 1, this%nbas
         if (this%bas(LIBCINT_ATOM_OF, ish) /= target_atom) cycle
         nprim = this%bas(LIBCINT_NPRIM_OF, ish)
         nctr = this%bas(LIBCINT_NCTR_OF, ish)
         atom_mol%nbas = atom_mol%nbas + 1
         env_size = env_size + nprim + nprim*nctr
      end do

      if (atom_mol%nbas == 0) then
         call error%set(ERROR_VALIDATION, "libcint: this atom carries no basis functions")
         return
      end if

      atom_mol%natm = 1
      atom_mol%cartesian = this%cartesian

      allocate (atom_mol%atm(LIBCINT_ATM_SLOTS, 1))
      allocate (atom_mol%bas(LIBCINT_BAS_SLOTS, atom_mol%nbas))
      allocate (atom_mol%env(env_size))
      allocate (atom_mol%shell_offset(atom_mol%nbas + 1))
      allocate (atom_mol%charges(1))
      allocate (atom_mol%coords(3, 1))

      ! Zeroed for the same reason `molecule_build` zeroes: libcint reads slots
      ! nothing here sets, and garbage in a pointer slot crashes inside the
      ! library with nothing to say why.
      atom_mol%atm = 0
      atom_mol%bas = 0
      atom_mol%env = 0.0_dp

      atom_mol%atm(LIBCINT_CHARGE_OF, 1) = this%atm(LIBCINT_CHARGE_OF, iatom)
      atom_mol%atm(LIBCINT_PTR_COORD, 1) = LIBCINT_PTR_ENV_START
      atom_mol%charges(1) = this%charges(iatom)
      atom_mol%coords = 0.0_dp

      off = LIBCINT_PTR_ENV_START + 3
      k = 0
      do ish = 1, this%nbas
         if (this%bas(LIBCINT_ATOM_OF, ish) /= target_atom) cycle
         k = k + 1
         nprim = this%bas(LIBCINT_NPRIM_OF, ish)
         nctr = this%bas(LIBCINT_NCTR_OF, ish)
         ncoeff = nprim*nctr

         atom_mol%bas(LIBCINT_ATOM_OF, k) = 0
         atom_mol%bas(LIBCINT_ANG_OF, k) = this%bas(LIBCINT_ANG_OF, ish)
         atom_mol%bas(LIBCINT_NPRIM_OF, k) = nprim
         atom_mol%bas(LIBCINT_NCTR_OF, k) = nctr

         atom_mol%bas(LIBCINT_PTR_EXP, k) = off
         atom_mol%env(off + 1:off + nprim) = &
            this%env(this%bas(LIBCINT_PTR_EXP, ish) + 1:this%bas(LIBCINT_PTR_EXP, ish) + nprim)
         off = off + nprim

         atom_mol%bas(LIBCINT_PTR_COEFF, k) = off
         atom_mol%env(off + 1:off + ncoeff) = &
            this%env(this%bas(LIBCINT_PTR_COEFF, ish) + 1:this%bas(LIBCINT_PTR_COEFF, ish) + ncoeff)
         off = off + ncoeff
      end do

      atom_mol%shell_offset(1) = 0
      do k = 1, atom_mol%nbas
         atom_mol%shell_offset(k + 1) = atom_mol%shell_offset(k) &
                                        + shell_dim(atom_mol%cartesian, k - 1, atom_mol%bas)
      end do
      atom_mol%nao = atom_mol%shell_offset(atom_mol%nbas + 1)
   end subroutine molecule_atom_subset

   pure subroutine atom_ao_blocks(mol, offsets, counts)
      !! Where each atom's basis functions sit in the molecular matrices
      !!
      !! `molecule_build` walks atoms outermost when it packs shells, so every
      !! atom's functions are one contiguous run and an atomic block can be
      !! copied into a row range rather than scattered. This is the routine that
      !! records that fact; if the packing order ever changed, this is what
      !! would have to change with it.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(out) :: offsets(:)   !! First AO of each atom, 0-based
      integer, intent(out) :: counts(:)    !! Functions on each atom

      integer :: ish, iatom

      counts = 0
      do ish = 1, mol%nbas
         iatom = mol%bas(LIBCINT_ATOM_OF, ish) + 1
         counts(iatom) = counts(iatom) + mol%shell_offset(ish + 1) - mol%shell_offset(ish)
      end do

      offsets(1) = 0
      do iatom = 2, mol%natm
         offsets(iatom) = offsets(iatom - 1) + counts(iatom - 1)
      end do
   end subroutine atom_ao_blocks

   subroutine subshell_layout(mol, ang, first, ncomp, n_sub)
      !! One entry per contraction column: its angular momentum and AO range
      !!
      !! A libcint shell with `NCTR_OF` columns lays its functions out with the
      !! contraction index outermost -- all components of column one, then all
      !! of column two -- so each column is a contiguous run of `2l+1` (or
      !! `(l+1)(l+2)/2`) functions. Those runs are the sets a spherical average
      !! has to act within, which is the only reason this is worth extracting.
      type(libcint_molecule_t), intent(in) :: mol
      integer, allocatable, intent(out) :: ang(:)     !! Angular momentum per column
      integer, allocatable, intent(out) :: first(:)   !! First AO of the column, 0-based
      integer, allocatable, intent(out) :: ncomp(:)   !! Functions in the column
      integer, intent(out) :: n_sub

      integer :: ish, ictr, nctr, per, total

      n_sub = 0
      do ish = 1, mol%nbas
         n_sub = n_sub + mol%bas(LIBCINT_NCTR_OF, ish)
      end do

      allocate (ang(n_sub), first(n_sub), ncomp(n_sub))

      n_sub = 0
      do ish = 1, mol%nbas
         nctr = mol%bas(LIBCINT_NCTR_OF, ish)
         total = mol%shell_offset(ish + 1) - mol%shell_offset(ish)
         per = total/nctr
         do ictr = 1, nctr
            n_sub = n_sub + 1
            ang(n_sub) = mol%bas(LIBCINT_ANG_OF, ish)
            first(n_sub) = mol%shell_offset(ish) + (ictr - 1)*per
            ncomp(n_sub) = per
         end do
      end do
   end subroutine subshell_layout

   pure function angular_form_name(cartesian) result(name)
      !! "Cartesian" or "spherical", for an error message to say which is which
      logical, intent(in) :: cartesian
      character(len=:), allocatable :: name

      if (cartesian) then
         name = "Cartesian"
      else
         name = "spherical"
      end if
   end function angular_form_name

   pure function max_block(this) result(n)
      !! Largest number of functions any one shell contributes
      type(libcint_molecule_t), intent(in) :: this
      integer :: n

      integer :: ish

      n = 1
      do ish = 1, this%nbas
         n = max(n, this%shell_offset(ish + 1) - this%shell_offset(ish))
      end do
   end function max_block

   subroutine molecule_destroy(this)
      class(libcint_molecule_t), intent(inout) :: this

      if (allocated(this%atm)) deallocate (this%atm)
      if (allocated(this%bas)) deallocate (this%bas)
      if (allocated(this%env)) deallocate (this%env)
      if (allocated(this%shell_offset)) deallocate (this%shell_offset)
      if (allocated(this%bas_sp)) deallocate (this%bas_sp)
      if (allocated(this%env_sp)) deallocate (this%env_sp)
      if (allocated(this%shell_offset_sp)) deallocate (this%shell_offset_sp)
      if (allocated(this%sp_split_first)) deallocate (this%sp_split_first)
      if (allocated(this%charges)) deallocate (this%charges)
      if (allocated(this%coords)) deallocate (this%coords)
      this%natm = 0
      this%nbas = 0
      this%nbas_sp = 0
      this%nao = 0
      this%cartesian = .false.
   end subroutine molecule_destroy

end module mqc_libcint_integrals
