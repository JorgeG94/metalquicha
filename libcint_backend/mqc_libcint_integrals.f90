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
   use mqc_cgto, only: molecular_basis_type
   use mqc_basis_utils, only: find_basis_file
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use pic_lapack_interfaces, only: pic_syev
   use libcint_fortran, only: libcint_1e_ovlp_sph, libcint_1e_kin_sph, &
                              libcint_3c2e_sph, libcint_2c2e_sph, &
                              libcint_1e_nuc_sph, libcint_2e_sph, &
                              libcint_cgto_sph, libcint_tot_cgto_sph, &
                              libcint_gto_norm, &
                              LIBCINT_ATM_SLOTS, LIBCINT_BAS_SLOTS, &
                              LIBCINT_CHARGE_OF, LIBCINT_PTR_COORD, &
                              LIBCINT_ATOM_OF, LIBCINT_ANG_OF, LIBCINT_NPRIM_OF, &
                              LIBCINT_NCTR_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_COEFF, &
                              LIBCINT_PTR_ENV_START
   implicit none
   private

   public :: libcint_molecule_t
   public :: build_libcint_molecule
   public :: build_df_tensor

   type :: libcint_molecule_t
      !! One molecule, packed the way libcint wants it
      integer :: natm = 0
      integer :: nbas = 0
      integer :: nao = 0
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

   subroutine molecule_build(this, atomic_numbers, coordinates, basis, error)
      !! Pack atoms and shells into libcint's atm/bas/env
      class(libcint_molecule_t), intent(inout) :: this
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coordinates(:, :)   !! (3, natm), Bohr
      type(molecular_basis_type), intent(in) :: basis
      type(error_t), intent(inout) :: error

      integer :: iatom, ishell, iprim, nprim, off, env_size, ang
      integer :: shell_index

      this%natm = size(atomic_numbers)
      if (basis%nelements /= this%natm) then
         call error%set(ERROR_VALIDATION, "libcint: the basis covers a different "// &
                        "number of atoms than the geometry has")
         return
      end if

      this%nbas = 0
      env_size = LIBCINT_PTR_ENV_START + 3*this%natm
      do iatom = 1, this%natm
         this%nbas = this%nbas + basis%elements(iatom)%nshells
         do ishell = 1, basis%elements(iatom)%nshells
            ! Exponents and coefficients, stored back to back.
            env_size = env_size + 2*basis%elements(iatom)%shells(ishell)%nfunc
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
         do ishell = 1, basis%elements(iatom)%nshells
            shell_index = shell_index + 1
            ang = basis%elements(iatom)%shells(ishell)%ang_mom
            nprim = basis%elements(iatom)%shells(ishell)%nfunc

            this%bas(LIBCINT_ATOM_OF, shell_index) = iatom - 1   ! libcint counts from 0
            this%bas(LIBCINT_ANG_OF, shell_index) = ang
            this%bas(LIBCINT_NPRIM_OF, shell_index) = nprim
            this%bas(LIBCINT_NCTR_OF, shell_index) = 1

            this%bas(LIBCINT_PTR_EXP, shell_index) = off
            this%env(off + 1:off + nprim) = basis%elements(iatom)%shells(ishell)%exponents
            off = off + nprim

            this%bas(LIBCINT_PTR_COEFF, shell_index) = off
            do iprim = 1, nprim
               ! The normalisation libcint expects to have been applied.
               this%env(off + iprim) = &
                  basis%elements(iatom)%shells(ishell)%coefficients(iprim) &
                  *libcint_gto_norm(ang, basis%elements(iatom)%shells(ishell)%exponents(iprim))
            end do
            off = off + nprim
         end do
      end do

      ! Where each shell's functions start, so a shell-pair block knows where
      ! it lands in the matrix.
      this%shell_offset(1) = 0
      do shell_index = 1, this%nbas
         this%shell_offset(shell_index + 1) = this%shell_offset(shell_index) &
                                              + libcint_cgto_sph(shell_index - 1, this%bas)
      end do
      this%nao = this%shell_offset(this%nbas + 1)

      if (this%nao /= libcint_tot_cgto_sph(this%bas, this%nbas)) then
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
         di = libcint_cgto_sph(ish - 1, this%bas)
         io = this%shell_offset(ish)
         do jsh = 1, this%nbas
            dj = libcint_cgto_sph(jsh - 1, this%bas)
            jo = this%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]

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

      nao = orb%nao
      naux = aux%nao

      call two_centre(aux, metric)
      call three_centre(orb, aux, three)

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

      ! J^(-1/2) = U s^(-1/2) U^T over the surviving modes.
      allocate (half(naux, naux))
      half = 0.0_dp
      do i = 1, naux
         if (values(i) <= NULL_THRESHOLD) cycle
         do j = 1, naux
            half(:, j) = half(:, j) + vectors(:, i)*vectors(j, i)/sqrt(values(i))
         end do
      end do

      allocate (b(nao*nao, naux))
      call pic_gemm(three, half, b)
   end subroutine build_df_tensor

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
         di = libcint_cgto_sph(ish - 1, aux%bas)
         io = aux%shell_offset(ish)
         do jsh = 1, aux%nbas
            dj = libcint_cgto_sph(jsh - 1, aux%bas)
            jo = aux%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]
            ret = libcint_2c2e_sph(buf, shls, aux%atm, aux%natm, aux%bas, aux%nbas, aux%env)
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
      type(libcint_molecule_t), intent(in) :: orb, aux
      real(dp), allocatable, intent(out) :: three(:, :)

      integer, allocatable :: bas(:, :)
      real(dp), allocatable :: env(:), buf(:)
      integer :: nbas_orb, nbas_aux, dummy, n_env_orb
      integer :: shls(4)
      integer :: ish, jsh, ksh, di, dj, dk, i, j, k, io, jo, ko, ret, idx

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
      allocate (buf(max_block(orb)**2*max_block(aux)))

      do ish = 1, nbas_orb
         di = libcint_cgto_sph(ish - 1, bas)
         io = orb%shell_offset(ish)
         do jsh = 1, nbas_orb
            dj = libcint_cgto_sph(jsh - 1, bas)
            jo = orb%shell_offset(jsh)
            do ksh = 1, nbas_aux
               dk = libcint_cgto_sph(nbas_orb + ksh - 1, bas)
               ko = aux%shell_offset(ksh)
               shls = [ish - 1, jsh - 1, nbas_orb + ksh - 1, dummy - 1]
               ret = libcint_3c2e_sph(buf, shls, orb%atm, orb%natm, bas, &
                                      nbas_orb + nbas_aux + 1, env)
               if (ret == 0) cycle
               do k = 1, dk
                  do j = 1, dj
                     do i = 1, di
                        idx = i + (j - 1)*di + (k - 1)*di*dj
                        three((jo + j - 1)*orb%nao + io + i, ko + k) = buf(idx)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine three_centre

   subroutine molecule_eris(this, eri)
      !! Every two-electron integral, in core, as (ij|kl)
      !!
      !! Stored as a full n^4 array with no permutational symmetry exploited.
      !! For a reference implementation that is the right trade: the eightfold
      !! symmetry saves memory and costs an index map, and an index map is
      !! exactly the sort of thing that would need its own validation.
      class(libcint_molecule_t), intent(in) :: this
      real(dp), allocatable, intent(out) :: eri(:, :, :, :)

      real(dp), allocatable :: buf(:)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: shls(4)
      integer :: i, j, k, l, io, jo, ko, lo, ret, idx

      allocate (eri(this%nao, this%nao, this%nao, this%nao))
      eri = 0.0_dp
      ! Four shells, so the fourth power. Same reasoning as the pair above.
      allocate (buf(max_block(this)**4))

      do ish = 1, this%nbas
         di = libcint_cgto_sph(ish - 1, this%bas)
         io = this%shell_offset(ish)
         do jsh = 1, this%nbas
            dj = libcint_cgto_sph(jsh - 1, this%bas)
            jo = this%shell_offset(jsh)
            do ksh = 1, this%nbas
               dk = libcint_cgto_sph(ksh - 1, this%bas)
               ko = this%shell_offset(ksh)
               do lsh = 1, this%nbas
                  dl = libcint_cgto_sph(lsh - 1, this%bas)
                  lo = this%shell_offset(lsh)
                  shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]
                  ret = libcint_2e_sph(buf, shls, this%atm, this%natm, &
                                       this%bas, this%nbas, this%env)
                  if (ret == 0) cycle

                  do l = 1, dl
                     do k = 1, dk
                        do j = 1, dj
                           do i = 1, di
                              idx = i + (j - 1)*di + (k - 1)*di*dj + (l - 1)*di*dj*dk
                              eri(io + i, jo + j, ko + k, lo + l) = buf(idx)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
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
   end subroutine molecule_destroy

end module mqc_libcint_integrals
