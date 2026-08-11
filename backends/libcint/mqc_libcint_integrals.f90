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
                              LIBCINT_ATM_SLOTS, LIBCINT_BAS_SLOTS, &
                              LIBCINT_CHARGE_OF, LIBCINT_PTR_COORD, &
                              LIBCINT_ATOM_OF, LIBCINT_ANG_OF, LIBCINT_NPRIM_OF, &
                              LIBCINT_NCTR_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_COEFF, &
                              LIBCINT_PTR_ENV_START
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
   implicit none
   private

   public :: libcint_molecule_t
   public :: build_libcint_molecule
   public :: build_df_tensor
   public :: build_df_mo_tensor
   ! The one place that maps "is this basis Cartesian?" onto libcint's two sets
   ! of entry points. Exported because the direct Fock build needs the same
   ! mapping, and two copies of it is two chances to route half the calls.
   public :: shell_dim
   public :: two_electron_block
   public :: two_electron_optimizer

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
      real(dp), allocatable :: charges(:)      !! Nuclear charges, for repulsion
      real(dp), allocatable :: coords(:, :)    !! (3, natm), Bohr
   contains
      procedure :: build => molecule_build
      procedure :: overlap => molecule_overlap
      procedure :: core_hamiltonian => molecule_core_hamiltonian
      procedure :: eris => molecule_eris
      procedure :: nuclear_repulsion => molecule_nuclear_repulsion
      procedure :: destroy => molecule_destroy
   end type libcint_molecule_t

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
                                     basis_name, mol, error)
      !! A molecule from a basis set *name*, through the ordinary reader
      !!
      !! This is what makes the backend general rather than a demonstration:
      !! any of the basis sets in `basis_sets/` rather than whatever was
      !! typed into a test.
      !!
      !! No normalisation is applied here beyond libcint's own. The BSE files
      !! give contraction coefficients against normalised primitives, the
      !! reader passes them through untouched, and `molecule_build` folds in
      !! `libcint_gto_norm` -- which is the one convention libcint asks for.
      !! Applying `normalize_shell_coefficients` as well would count it twice,
      !! and the symptom is an overlap diagonal that is not 1.
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)   !! (3, natm), Bohr
      character(len=*), intent(in) :: basis_name
      type(libcint_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: error

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

      call mol%build(atomic_numbers, coordinates, basis, error)
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

   subroutine molecule_build(this, atomic_numbers, coordinates, basis, error)
      !! Pack atoms and shells into libcint's atm/bas/env
      class(libcint_molecule_t), intent(inout) :: this
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)   !! (3, natm), Bohr
      type(molecular_basis_type), intent(in) :: basis
      type(error_t), intent(inout) :: error

      integer :: iatom, ishell, iprim, nprim, off, env_size, ang
      integer :: shell_index, ictr, nctr, first

      this%natm = size(atomic_numbers)
      ! What the basis file said, carried through the reader. Every integral
      ! call below routes on it, and so does nao.
      this%cartesian = basis%is_cartesian()
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
         this%atm(LIBCINT_CHARGE_OF, iatom) = atomic_numbers(iatom)
         this%atm(LIBCINT_PTR_COORD, iatom) = off
         this%env(off + 1:off + 3) = coordinates(1:3, iatom)
         this%charges(iatom) = real(atomic_numbers(iatom), dp)
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
      end if
   end subroutine molecule_build

   subroutine molecule_overlap(this, s)
      !! S, shell pair by shell pair
      class(libcint_molecule_t), intent(in) :: this
      real(dp), allocatable, intent(out) :: s(:, :)

      call one_electron(this, s, 1)
   end subroutine molecule_overlap

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

   subroutine build_df_tensor(orb, aux, b, error)
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

      real(dp), parameter :: NULL_THRESHOLD = 1.0e-10_dp
      real(dp), allocatable :: metric(:, :), vectors(:, :), values(:), half(:, :)
      real(dp), allocatable :: three(:, :)
      integer :: naux, nao, i, j, kept, info

      ! (mu nu | P) has orbital shells on two centres and auxiliary shells on
      ! the third, in one libcint call -- and that call is spherical or
      ! Cartesian for all three at once. A Cartesian orbital basis fitted with
      ! a spherical auxiliary one cannot be expressed, and guessing which side
      ! to honour would silently fit in a basis neither deck asked for.
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

      call two_centre(aux, metric)
      call three_centre(orb, aux, three)

      call metric_inverse_sqrt(metric, half, error)
      if (error%has_error()) return

      allocate (b(nao*nao, naux))
      call pic_gemm(three, half, b)
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

      real(dp), allocatable :: metric(:, :), half(:, :), three(:, :)
      real(dp), allocatable :: ovl(:, :), bp(:, :), tmp(:, :), full(:, :)
      integer :: nao, naux, n_o, n_v, p_index, i

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
      n_o = size(c_occ, 2)
      n_v = size(c_vir, 2)

      call two_centre(aux, metric)
      call metric_inverse_sqrt(metric, half, error)
      if (error%has_error()) return
      deallocate (metric)

      call three_centre(orb, aux, three)

      ! (mu nu|P) -> (ia|P), one auxiliary function at a time.
      allocate (ovl(n_o*n_v, naux))
      allocate (bp(nao, nao), tmp(n_o, nao), full(n_o, n_v))
      do p_index = 1, naux
         bp = reshape(three(:, p_index), [nao, nao])
         call pic_gemm(c_occ, bp, tmp, transa="T")
         call pic_gemm(tmp, c_vir, full)
         ovl(:, p_index) = reshape(full, [n_o*n_v])
      end do
      deallocate (bp, tmp, full, three)

      ! The fit, now over the occupied-virtual block instead of every AO pair.
      block
         real(dp), allocatable :: fitted(:, :)
         allocate (fitted(n_o*n_v, naux))
         call pic_gemm(ovl, half, fitted)
         deallocate (ovl, half)
         allocate (bia(n_v, naux, n_o))
         do p_index = 1, naux
            do i = 1, n_o
               ! `fitted` is flattened (i, a), so one occupied orbital's
               ! virtuals sit at stride n_o.
               bia(:, p_index, i) = fitted(i::n_o, p_index)
            end do
         end do
         deallocate (fitted)
      end block
   end subroutine build_df_mo_tensor

   subroutine two_centre(aux, metric)
      !! (P|Q) over the auxiliary basis
      type(libcint_molecule_t), intent(in) :: aux
      real(dp), allocatable, intent(out) :: metric(:, :)

      real(dp), allocatable :: buf(:)
      integer :: shls(2)
      integer :: ish, jsh, di, dj, i, j, io, jo, ret

      allocate (metric(aux%nao, aux%nao))
      metric = 0.0_dp
      allocate (buf(max_block(aux)**2))

      do ish = 1, aux%nbas
         di = shell_dim(aux%cartesian, ish - 1, aux%bas)
         io = aux%shell_offset(ish)
         do jsh = 1, aux%nbas
            dj = shell_dim(aux%cartesian, jsh - 1, aux%bas)
            jo = aux%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]
            if (aux%cartesian) then
               ret = libcint_2c2e_cart(buf, shls, aux%atm, aux%natm, aux%bas, aux%nbas, aux%env)
            else
               ret = libcint_2c2e_sph(buf, shls, aux%atm, aux%natm, aux%bas, aux%nbas, aux%env)
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

   subroutine three_centre(orb, aux, three)
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

      integer, allocatable :: bas(:, :)
      real(dp), allocatable :: env(:), buf(:)
      integer :: nbas_orb, nbas_aux, dummy, n_env_orb
      integer :: shls(4)
      integer :: ish, jsh, ksh, di, dj, dk, i, j, k, io, jo, ko, ret, idx
      integer :: npair, ipair
      integer, allocatable :: pair_i(:), pair_j(:)
      type(c_ptr) :: opt

      nbas_orb = orb%nbas
      nbas_aux = aux%nbas

      ! Orbital shells, then auxiliary shells, then one dummy. The auxiliary
      ! env offsets are shifted by the orbital env length because both live in
      ! one array now.
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

      ! The dummy: one s shell, one primitive, exponent and coefficient zero.
      dummy = nbas_orb + nbas_aux + 1
      bas(LIBCINT_ATOM_OF, dummy) = 0
      bas(LIBCINT_ANG_OF, dummy) = 0
      bas(LIBCINT_NPRIM_OF, dummy) = 1
      bas(LIBCINT_NCTR_OF, dummy) = 1
      bas(LIBCINT_PTR_EXP, dummy) = size(env) - 1
      bas(LIBCINT_PTR_COEFF, dummy) = size(env) - 1
      env(size(env)) = 0.0_dp

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
      integer :: shls(4)
      integer :: i, j, k, l, io, jo, ko, lo, ret, idx
      integer :: p, q, r, t
      real(dp) :: value
      type(c_ptr) :: opt

      allocate (eri(this%nao, this%nao, this%nao, this%nao))
      eri = 0.0_dp
      allocate (buf(max_block(this)**4))

      opt = c_null_ptr
      call two_electron_optimizer(this%cartesian, opt, this%atm, this%natm, this%bas, &
                                  this%nbas, this%env)

      do ish = 1, this%nbas
         di = shell_dim(this%cartesian, ish - 1, this%bas)
         io = this%shell_offset(ish)
         do jsh = 1, ish
            dj = shell_dim(this%cartesian, jsh - 1, this%bas)
            jo = this%shell_offset(jsh)
            do ksh = 1, ish
               dk = shell_dim(this%cartesian, ksh - 1, this%bas)
               ko = this%shell_offset(ksh)
               lsh_max = ksh
               if (ksh == ish) lsh_max = jsh
               do lsh = 1, lsh_max
                  dl = shell_dim(this%cartesian, lsh - 1, this%bas)
                  lo = this%shell_offset(lsh)
                  shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]
                  ret = two_electron_block(this%cartesian, buf, shls, this%atm, &
                                           this%natm, this%bas, this%nbas, this%env, opt)
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
      end do

      call libcint_del_optimizer(opt)
      deallocate (buf)
   end subroutine molecule_eris

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
      if (allocated(this%charges)) deallocate (this%charges)
      if (allocated(this%coords)) deallocate (this%coords)
      this%natm = 0
      this%nbas = 0
      this%nao = 0
      this%cartesian = .false.
   end subroutine molecule_destroy

end module mqc_libcint_integrals
