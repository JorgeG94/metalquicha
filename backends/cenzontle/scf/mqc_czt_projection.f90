!! A converged small-basis density, carried into a larger basis as a starting point
module mqc_czt_projection
   !! Converge the cheap basis first, then start the expensive one from it.
   !!
   !! For a system that converges badly the cost is rarely the iterations in the
   !! target basis, but that a core or GWH guess is far enough from the answer
   !! that the iteration wanders. A small-basis calculation on the same geometry
   !! costs almost nothing and lands close to the right density.
   !!
   !! **What is projected is the occupied orbitals, not the density.** The obvious
   !! thing is
   !!
   !!     D_B = S_BB^-1 S_BA D_A S_AB S_BB^-1
   !!
   !! which is one expression and wrong in a way that matters: the result is not
   !! idempotent and its trace against S is not the electron count, so the first
   !! Fock matrix is built from something that is not a density. Projecting the
   !! occupied orbitals and re-orthonormalising them gives a D that is idempotent
   !! and has exactly the right trace by construction, which is worth the extra
   !! eigendecomposition of an n_occ-square matrix.
   !!
   !!     C~ = S_BB^-1 S_BA C_A(occ)      project
   !!     M  = C~^T S_BB C~               overlap of the projected orbitals
   !!     C  = C~ M^-1/2                  re-orthonormalise
   !!     D  = 2 C C^T
   !!
   !! **The cross-basis overlap needs both bases in one molecule.** libcint
   !! computes integrals over shells of the molecule it is given, so S_BA -- rows
   !! in the target basis, columns in the small one -- does not exist for either
   !! molecule alone. `merge_basis_sets` builds a third whose shell list is the
   !! target's followed by the small one's over the same atoms, and the wanted
   !! block is the off-diagonal corner of its overlap, as PySCF's `intor_cross`
   !! does.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_scf_types, only: guess_step_t, scf_numerics_t
   use mqc_czt_integrals, only: czt_molecule_t, shell_dim, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, SCF_GUESS_SAD, SCF_GUESS_GWH
   use mqc_diis, only: parse_accelerator_name, ACCEL_DIIS
   use libcint_fortran, only: LIBCINT_BAS_SLOTS, LIBCINT_PTR_EXP, &
                              LIBCINT_PTR_COEFF, LIBCINT_PTR_ENV_START
   implicit none
   private

   public :: merge_basis_sets
   public :: cross_overlap
   public :: project_occupied
   public :: climb_basis_ladder

   ! Overlap eigenvalues below this are dropped when inverting. The same cutoff
   ! `build_orthogonalizer` uses: a near-dependent target basis has directions the
   ! small basis cannot inform, and amplifying them by dividing through a tiny
   ! eigenvalue produces a guess worse than the one it replaces.
   real(dp), parameter :: OVERLAP_FLOOR = 1.0e-7_dp

contains

   subroutine climb_basis_ladder(steps, atomic_numbers, symbols, coords, nelec, &
                                 mol_target, density, error, verbose, scf_in)
      !! Converge each basis in turn, projecting forward, and return a density
      !! in the target basis
      !!
      !! Each rung starts from the previous rung's density rather than from a
      !! core guess, so the first SCF is the only one starting from nothing.
      !!
      !! The target basis is not a rung. `mol_target` is already built by the
      !! caller and the last projection lands in it, so a two-step ladder plus a
      !! model basis is three SCFs, the last of which is the caller's.
      type(guess_step_t), intent(in) :: steps(:)
         !! Only `maxiter` and `tolerance` are per-rung. How the iteration is
         !! *driven* is a property of the calculation and lives in `scf_in`.
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: nelec
      type(czt_molecule_t), intent(in) :: mol_target
      real(dp), allocatable, intent(out) :: density(:, :)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose
      type(scf_numerics_t), intent(in), optional :: scf_in
         !! How to drive each rung's SCF: the accelerator, the DIIS subspace,
         !! the level shift, the linear-dependence threshold and incremental
         !! Fock building. Absent leaves every rung on the defaults.
         !!
         !! `linear_dependence` is the one that matters: a ladder exists to climb
         !! into a large diffuse basis, which is exactly where the overlap goes
         !! near-singular.

      type(czt_molecule_t), allocatable :: rungs(:)
      type(rhf_result_t) :: scf
      real(dp) :: shift, lindep
      integer :: accel, ndiis
      logical :: incr, accel_ok
      real(dp), allocatable :: carried(:, :)
      integer :: i, n, n_occ
      logical :: talk

      talk = .false.
      if (present(verbose)) talk = verbose
      n = size(steps)
      n_occ = nelec/2

      if (mod(nelec, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "basis projection guess: the ladder converges "// &
                        "closed-shell SCFs and this fragment has an odd electron count")
         return
      end if

      ! The SCF settings are resolved once for the whole ladder.
      shift = 0.0_dp
      lindep = 0.0_dp
      accel = ACCEL_DIIS
      ndiis = 8
      incr = .true.
      if (present(scf_in)) then
         shift = scf_in%level_shift
         lindep = scf_in%linear_dependence
         ndiis = scf_in%diis_size
         if (.not. scf_in%use_diis) ndiis = 0
         incr = scf_in%incremental_fock
         call parse_accelerator_name(scf_in%accelerator, accel, accel_ok)
         ! A spelling this routine cannot parse has already been refused by the
         ! caller, which runs the same parse before any of this. Falling back
         ! rather than erroring keeps a guess from failing a calculation.
         if (.not. accel_ok) accel = ACCEL_DIIS
      end if

      ! Every rung is built up front. The last one is needed after its SCF to
      ! project into the next, and a molecule that has been destroyed cannot be
      ! asked for its overlap.
      allocate (rungs(n))
      do i = 1, n
         call build_czt_molecule(atomic_numbers, symbols, coords, steps(i)%basis, &
                                 rungs(i), error)
         if (error%has_error()) then
            call error%set(ERROR_VALIDATION, "basis projection guess: rung "//to_char(i)// &
                           " ("//trim(steps(i)%basis)//"): "//error%get_message())
            return
         end if
      end do

      do i = 1, n
         if (i == 1) then
            ! GWH rather than SAD for the bottom rung: SAD needs a density handed
            ! to it and there is none yet, and this basis is small enough that
            ! the difference between the two starting points costs an iteration
            ! rather than a minute.
            call run_czt_rhf(rungs(i), nelec, steps(i)%maxiter, steps(i)%tolerance, &
                             steps(i)%tolerance, .false., scf, error, &
                             guess=SCF_GUESS_GWH, &
                             level_shift=shift, linear_dependence=lindep, &
                             accelerator=accel, diis_vectors=ndiis, &
                             incremental_fock=incr)
         else
            call run_czt_rhf(rungs(i), nelec, steps(i)%maxiter, steps(i)%tolerance, &
                             steps(i)%tolerance, .false., scf, error, &
                             guess=SCF_GUESS_SAD, guess_density=carried, &
                             level_shift=shift, linear_dependence=lindep, &
                             accelerator=accel, diis_vectors=ndiis, &
                             incremental_fock=incr)
         end if
         if (error%has_error()) return
         ! A rung that does not converge is not fatal: a half-converged density
         ! from a smaller basis is still closer to the answer than a core guess.
         ! Reported, because a ladder that never converges anywhere is worth
         ! knowing about.
         if (talk) then
            call logger%info("  guess rung "//to_char(i)//": "//trim(steps(i)%basis)// &
                             ", "//to_char(rungs(i)%nao)//" functions, "// &
                             to_char(scf%iterations)//" iterations"// &
                             merge("                ", " (not converged)", scf%converged))
         end if

         if (i < n) then
            if (allocated(carried)) deallocate (carried)
            call project_occupied(rungs(i + 1), rungs(i), scf%orbitals, n_occ, carried, error)
         else
            call project_occupied(mol_target, rungs(i), scf%orbitals, n_occ, density, error)
         end if
         if (error%has_error()) return
      end do
   end subroutine climb_basis_ladder

   subroutine merge_basis_sets(mol_target, mol_small, merged, error)
      !! One molecule carrying both shell lists, target first
      !!
      !! Only the shells are concatenated. Both molecules describe the same atoms
      !! at the same coordinates, and `molecule_build` writes those into `env`
      !! before any basis data, at offsets that depend on the atom count alone, so
      !! the target's `atm` already points at coordinates that are correct for the
      !! merged molecule and the small basis contributes exponents and
      !! coefficients only.
      type(czt_molecule_t), intent(in) :: mol_target, mol_small
      type(czt_molecule_t), intent(out) :: merged
      type(error_t), intent(inout) :: error

      integer :: coord_end, tail, shift, ish, ns, nt

      if (mol_target%natm /= mol_small%natm) then
         call error%set(ERROR_VALIDATION, "basis projection: the two bases describe "// &
                        "different numbers of atoms")
         return
      end if
      if (maxval(abs(mol_target%coords - mol_small%coords)) > 1.0e-10_dp) then
         call error%set(ERROR_VALIDATION, "basis projection: the two bases are on "// &
                        "different geometries")
         return
      end if
      ! One angular form for the whole merged molecule, because libcint chooses it
      ! per call rather than per shell -- the same constraint `cartesian` carries
      ! for a single molecule, which here means the two bases have to agree.
      if (mol_target%cartesian .neqv. mol_small%cartesian) then
         call error%set(ERROR_VALIDATION, "basis projection: one basis is Cartesian "// &
                        "and the other spherical; a merged molecule can only be one "// &
                        "or the other")
         return
      end if

      nt = mol_target%nbas
      ns = mol_small%nbas
      coord_end = LIBCINT_PTR_ENV_START + 3*mol_target%natm
      tail = size(mol_small%env) - coord_end          !! small basis data only
      shift = size(mol_target%env) - coord_end        !! where it lands in `merged`

      merged%natm = mol_target%natm
      merged%nbas = nt + ns
      merged%cartesian = mol_target%cartesian
      merged%atm = mol_target%atm
      merged%charges = mol_target%charges
      merged%coords = mol_target%coords

      allocate (merged%env(size(mol_target%env) + tail))
      merged%env(1:size(mol_target%env)) = mol_target%env
      merged%env(size(mol_target%env) + 1:) = mol_small%env(coord_end + 1:)

      allocate (merged%bas(LIBCINT_BAS_SLOTS, merged%nbas))
      merged%bas(:, 1:nt) = mol_target%bas
      merged%bas(:, nt + 1:) = mol_small%bas
      do ish = nt + 1, merged%nbas
         merged%bas(LIBCINT_PTR_EXP, ish) = merged%bas(LIBCINT_PTR_EXP, ish) + shift
         merged%bas(LIBCINT_PTR_COEFF, ish) = merged%bas(LIBCINT_PTR_COEFF, ish) + shift
      end do

      allocate (merged%shell_offset(merged%nbas + 1))
      merged%shell_offset(1) = 0
      do ish = 1, merged%nbas
         merged%shell_offset(ish + 1) = merged%shell_offset(ish) + &
                                        shell_dim(merged%cartesian, ish - 1, merged%bas)
      end do
      merged%nao = merged%shell_offset(merged%nbas + 1)

      if (merged%nao /= mol_target%nao + mol_small%nao) then
         call error%set(ERROR_VALIDATION, "basis projection: the merged molecule has "// &
                        to_char(merged%nao)//" functions where the two bases have "// &
                        to_char(mol_target%nao + mol_small%nao)//" between them")
         return
      end if
   end subroutine merge_basis_sets

   subroutine cross_overlap(mol_target, mol_small, s_target, s_cross, error)
      !! S_BB and S_BA in one pass over the merged molecule
      !!
      !! Both are blocks of the merged molecule's overlap, so they are returned
      !! together rather than computing it twice.
      type(czt_molecule_t), intent(in) :: mol_target, mol_small
      real(dp), allocatable, intent(out) :: s_target(:, :)   !! (n_B, n_B)
      real(dp), allocatable, intent(out) :: s_cross(:, :)    !! (n_B, n_A)
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: merged
      real(dp), allocatable :: s_all(:, :)
      integer :: nb, na

      call merge_basis_sets(mol_target, mol_small, merged, error)
      if (error%has_error()) return

      call merged%overlap(s_all)
      nb = mol_target%nao
      na = mol_small%nao
      s_target = s_all(1:nb, 1:nb)
      s_cross = s_all(1:nb, nb + 1:nb + na)
   end subroutine cross_overlap

   subroutine project_occupied(mol_target, mol_small, coeff_small, n_occ, density, error)
      !! A target-basis density from the small basis's occupied orbitals
      !!
      !! The result is idempotent against S_BB and its trace against S_BB is
      !! exactly 2*n_occ, both by construction.
      type(czt_molecule_t), intent(in) :: mol_target, mol_small
      real(dp), intent(in) :: coeff_small(:, :)   !! (n_A, n_mo_A), occupied first
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: density(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: s_target(:, :), s_cross(:, :)
      real(dp), allocatable :: rhs(:, :), proj(:, :), m(:, :), work(:, :)
      real(dp), allocatable :: vec(:, :), val(:), half(:, :)
      integer :: nb, i, kept

      if (n_occ < 1) then
         call error%set(ERROR_VALIDATION, "basis projection: no occupied orbitals to project")
         return
      end if
      if (size(coeff_small, 2) < n_occ) then
         call error%set(ERROR_VALIDATION, "basis projection: fewer orbitals than occupied ones")
         return
      end if

      call cross_overlap(mol_target, mol_small, s_target, s_cross, error)
      if (error%has_error()) return
      nb = mol_target%nao

      ! S_BA C_A(occ)
      allocate (rhs(nb, n_occ))
      call pic_gemm(s_cross, coeff_small(:, 1:n_occ), rhs)

      ! S_BB^-1 applied through its eigendecomposition rather than a solve, so the
      ! near-null directions of a nearly dependent target basis can be dropped
      ! instead of divided by.
      call invert_overlap(s_target, work, kept, error)
      if (error%has_error()) return
      if (kept < n_occ) then
         call error%set(ERROR_VALIDATION, "basis projection: the target basis keeps only "// &
                        to_char(kept)//" independent directions, fewer than the "// &
                        to_char(n_occ)//" occupied orbitals being projected")
         return
      end if
      allocate (proj(nb, n_occ))
      call pic_gemm(work, rhs, proj)

      ! M = C~^T S_BB C~, the overlap the projected orbitals have among themselves.
      ! It is the identity only if the projection happened to be exact.
      allocate (m(n_occ, n_occ))
      deallocate (rhs)
      allocate (rhs(nb, n_occ))
      call pic_gemm(s_target, proj, rhs)
      call pic_gemm(proj, rhs, m, transa="T")

      ! C = C~ M^-1/2, by eigendecomposition. M is small -- occupied square -- and
      ! positive definite unless the projection collapsed two orbitals onto one,
      ! which the floor below catches.
      allocate (vec(n_occ, n_occ), val(n_occ))
      vec = m
      call pic_syev(vec, val, jobz="V", uplo="U")
      if (minval(val) <= OVERLAP_FLOOR) then
         call error%set(ERROR_VALIDATION, "basis projection: the projected orbitals are "// &
                        "linearly dependent; the small basis cannot represent this "// &
                        "occupied space")
         return
      end if
      ! M^-1/2 = V d^-1/2 V^T, so exactly one of the two factors carries the
      ! eigenvalue scaling. Scaling `vec` in place and multiplying it by itself
      ! would put d^-1/2 in twice and produce M^-1, which the self-projection
      ! check cannot see: there M is the identity and every power of it is the
      ! same matrix.
      allocate (half(n_occ, n_occ))
      do i = 1, n_occ
         half(:, i) = vec(:, i)/sqrt(val(i))
      end do
      deallocate (m)
      allocate (m(n_occ, n_occ))
      call pic_gemm(half, vec, m, transb="T")    !! M^-1/2

      deallocate (rhs)
      allocate (rhs(nb, n_occ))
      call pic_gemm(proj, m, rhs)                !! C = C~ M^-1/2

      allocate (density(nb, nb))
      call pic_gemm(rhs, rhs, density, transb="T")
      density = 2.0_dp*density
   end subroutine project_occupied

   subroutine invert_overlap(overlap, inverse, kept, error)
      !! S^-1 with near-null directions dropped rather than amplified
      real(dp), intent(in) :: overlap(:, :)
      real(dp), allocatable, intent(out) :: inverse(:, :)
      integer, intent(out) :: kept
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: vec(:, :), val(:), scaled(:, :)
      integer :: n, i

      n = size(overlap, 1)
      allocate (vec(n, n), val(n), scaled(n, n), inverse(n, n))
      vec = overlap
      call pic_syev(vec, val, jobz="V", uplo="U")

      kept = 0
      scaled = 0.0_dp
      do i = 1, n
         if (val(i) > OVERLAP_FLOOR) then
            kept = kept + 1
            scaled(:, i) = vec(:, i)/val(i)
         end if
      end do
      if (kept == 0) then
         call error%set(ERROR_VALIDATION, "basis projection: the target overlap has no "// &
                        "eigenvalue above the linear-dependence floor")
         return
      end if
      call pic_gemm(scaled, vec, inverse, transb="T")
   end subroutine invert_overlap

end module mqc_czt_projection
