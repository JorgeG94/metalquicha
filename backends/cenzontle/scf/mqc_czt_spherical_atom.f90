!! Spherically averaged free atoms, by an SCF with fractional occupations
module mqc_czt_spherical_atom
   !! One isolated atom solved so that its density is spherically symmetric by
   !! construction rather than by averaging afterwards.
   !!
   !! Within each angular momentum `l` the Fock and overlap matrices are
   !! averaged over the `2l+1` components and diagonalised once, and every
   !! orbital of that `l` is occupied equally across its components. The
   !! electrons per `l` are the spherically averaged ground-state configuration
   !! of `minao_occupation` -- PySCF's `NRSRHF_CONFIGURATION`, from Lehtola,
   !! Phys. Rev. A 101, 012516 (2020) -- so the open shell holds a fraction of
   !! an electron in each of its components and the density has one form only.
   !!
   !! Two spin treatments:
   !!
   !! * **polarized** (the default) splits each `l` by Hund's rule: its closed
   !!   shells to both spins, its open shell to alpha first. Exact for the
   !!   hydrogen atom.
   !! * **restricted** gives each spin half of every `l`, which is PySCF's
   !!   `atom_hf.AtomSphericAverageRHF` and `atom_ks.AtomSphericAverageRKS`, and
   !!   so the free atom of the Hirshfeld extension at
   !!   github.com/frobnitzem/hirshfeld. A lone electron then sees half of its
   !!   own Coulomb repulsion: hydrogen in cc-pVDZ comes out at -0.354 hartree
   !!   rather than -0.499, with a density to match.
   !!
   !! The Fock matrices are `assemble_fock_uhf`'s, so Hartree-Fock and every
   !! functional the molecular SCF accepts give the same operator here.
   !!
   !! Refused: more than one centre, a ghost, an effective core potential, an
   !! element outside H-Kr, a Cartesian shell of `l >= 2`, and a basis with too
   !! few functions of some `l` to hold its configuration.
   use pic_types, only: dp
   use pic_lapack_interfaces, only: pic_syev
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_string_utils, only: int_to_text
   use mqc_diis, only: diis_state_t
   use mqc_scf_common, only: build_orthogonalizer
   use mqc_czt_integrals, only: czt_molecule_t, subshell_layout
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_xc, only: xc_context_t, xc_context_create
   use mqc_czt_rhf, only: assemble_fock_uhf
   use mqc_czt_atomic_guess, only: spherical_average
   use mqc_czt_minao, only: minao_occupation
   implicit none
   private

   public :: spherical_atom_settings_t   !! How a free atom is solved
   public :: spherical_atom_t            !! One converged free atom
   public :: solve_spherical_atom        !! Converge one atom, spherically averaged
   public :: spherical_atom_electrons    !! Electrons per l and spin for an element

   integer, parameter :: L_OCCUPIED_MAX = 3
      !! The configurations stop at f.

   real(dp), parameter :: FRACTION_FLOOR = 1.0e-12_dp
      !! An occupation below this is an empty orbital, not a fractional one.

   type :: spherical_atom_settings_t
      !! How a free atom is solved
      character(len=:), allocatable :: functional
         !! Unallocated, empty or "hf" is Hartree-Fock; anything else is a
         !! functional name `xc_context_create` accepts.
      integer :: grid_level = -1
         !! Quadrature level for a functional. Negative keeps the default.
      logical :: polarized = .true.
         !! Hund's-rule spin polarisation within each l. False is restricted.
      integer :: max_iter = 200
      real(dp) :: energy_tol = 1.0e-10_dp
         !! Change in the energy between iterations, hartree.
      real(dp) :: commutator_tol = 1.0e-7_dp
         !! Largest element of `F D S - S D F`, either spin.
      integer :: diis_vectors = 8
   end type spherical_atom_settings_t

   type :: spherical_atom_t
      !! One converged free atom
      integer :: atomic_number = 0
      logical :: polarized = .true.
      real(dp) :: energy = 0.0_dp            !! Total energy, hartree
      real(dp), allocatable :: d_alpha(:, :)  !! (n_ao, n_ao), the atom's own basis
      real(dp), allocatable :: d_beta(:, :)   !! (n_ao, n_ao); equal to alpha when restricted
      integer :: iterations = 0
      logical :: converged = .false.
   end type spherical_atom_t

   type :: l_block_t
      ! The contraction columns of one angular momentum, and the canonical
      ! orthogonaliser of their averaged overlap. Fixed for the whole SCF.
      integer :: n_col = 0
      integer :: n_comp = 0                  ! 2l+1
      integer, allocatable :: column(:)      ! (n_col) index into the subshell layout
      real(dp), allocatable :: x(:, :)       ! (n_col, n_mo)
      integer :: n_mo = 0
   end type l_block_t

contains

   subroutine spherical_atom_electrons(atomic_number, polarized, electrons, error)
      !! Electrons in each angular momentum and spin, `electrons(l, spin)`
      !!
      !! `l` runs 0..3 and `spin` is 1 for alpha, 2 for beta. Integers when
      !! polarized; half-integers can appear when restricted, which halves every
      !! `l`.
      integer, intent(in) :: atomic_number
      logical, intent(in) :: polarized
      real(dp), intent(out) :: electrons(0:L_OCCUPIED_MAX, 2)
      type(error_t), intent(inout) :: error

      integer :: l, n_double, per_spin_shell, n_l, n_open
      real(dp) :: open_fraction

      electrons = 0.0_dp
      do l = 0, L_OCCUPIED_MAX
         call minao_occupation(atomic_number, l, n_double, open_fraction)
         per_spin_shell = 2*l + 1
         n_open = nint(open_fraction*real(per_spin_shell, dp))
         n_l = 2*per_spin_shell*n_double + n_open
         if (polarized) then
            electrons(l, 1) = real(per_spin_shell*n_double + min(n_open, per_spin_shell), dp)
            electrons(l, 2) = real(per_spin_shell*n_double + max(0, n_open - per_spin_shell), dp)
         else
            electrons(l, :) = 0.5_dp*real(n_l, dp)
         end if
      end do

      if (nint(sum(electrons)) /= atomic_number .or. atomic_number < 1) then
         call error%set(ERROR_VALIDATION, "spherical atom: no ground-state configuration "// &
                        "for Z="//int_to_text(atomic_number)//"; the table covers H-Kr")
      end if
   end subroutine spherical_atom_electrons

   subroutine solve_spherical_atom(atom_mol, settings, atom, error)
      !! Converge one atom with its density spherically symmetric throughout
      !!
      !! Starts from the core Hamiltonian, so it never needs a guess that itself
      !! needs an atom. `atom%converged` false is reported, not raised: the
      !! caller decides whether an unconverged free atom is good enough.
      type(czt_molecule_t), intent(in) :: atom_mol
      type(spherical_atom_settings_t), intent(in) :: settings
      type(spherical_atom_t), intent(out) :: atom
      type(error_t), intent(inout) :: error

      type(l_block_t) :: blocks(0:L_OCCUPIED_MAX)
      type(xc_context_t) :: xc
      type(diis_state_t) :: diis
      real(dp), allocatable :: h(:, :), s(:, :), bounds(:, :), eri(:, :, :, :)
      real(dp), allocatable :: fock(:, :, :), averaged(:, :), residual(:, :, :)
      real(dp), allocatable :: density(:, :, :)
      integer, allocatable :: ang(:), first(:), ncomp(:)
      real(dp) :: electrons(0:L_OCCUPIED_MAX, 2)
      real(dp) :: energy, previous, commutator
      integer :: n, n_sub, z, iter, spin
      logical :: exact, ok

      if (error%has_error()) return
      call check_atom(atom_mol, z, error)
      if (error%has_error()) return
      call spherical_atom_electrons(z, settings%polarized, electrons, error)
      if (error%has_error()) return

      n = atom_mol%nao
      call subshell_layout(atom_mol, ang, first, ncomp, n_sub)
      call atom_mol%overlap(s)
      call atom_mol%core_hamiltonian(h)
      call build_l_blocks(s, ang, first, ncomp, n_sub, electrons, blocks, error)
      if (error%has_error()) return
      call schwarz_bounds(atom_mol, bounds, error)
      if (error%has_error()) return

      if (wants_functional(settings)) then
         ! Spin-polarised whichever treatment was asked for: the restricted one
         ! is the polarised functional at equal spin densities.
         if (settings%grid_level >= 0) then
            call xc_context_create(atom_mol, settings%functional, xc, error, &
                                   level=settings%grid_level, polarized=.true.)
         else
            call xc_context_create(atom_mol, settings%functional, xc, error, &
                                   polarized=.true.)
         end if
         if (error%has_error()) then
            call error%add_context("spherical atom Z="//int_to_text(z))
            return
         end if
      end if

      allocate (fock(n, n, 2), density(n, n, 2), residual(n, n, 2), averaged(n, n))
      do spin = 1, 2
         call occupy(h, first, ncomp, blocks, electrons(:, spin), density(:, :, spin), &
                     error)
      end do
      if (error%has_error()) return

      call diis%init(settings%diis_vectors, 2*n*n, 2*n*n)
      atom%atomic_number = z
      atom%polarized = settings%polarized
      previous = huge(1.0_dp)

      do iter = 1, settings%max_iter
         call assemble_fock_uhf(atom_mol, h, density(:, :, 1), density(:, :, 2), eri, bounds, &
                                xc, fock(:, :, 1), fock(:, :, 2), energy, error)
         if (error%has_error()) exit

         ! The quadrature is not exactly spherical, so the operator is averaged
         ! every iteration rather than trusted to come out averaged.
         do spin = 1, 2
            call spherical_average(atom_mol, fock(:, :, spin), averaged, exact)
            fock(:, :, spin) = averaged
            call commutator_residual(fock(:, :, spin), density(:, :, spin), s, &
                                     residual(:, :, spin))
         end do
         commutator = maxval(abs(residual))

         atom%iterations = iter
         atom%energy = energy
         if (abs(energy - previous) < settings%energy_tol .and. &
             commutator < settings%commutator_tol) then
            atom%converged = .true.
            exit
         end if
         previous = energy

         call diis%push(fock, residual)
         call diis%extrapolate(fock, ok)
         if (.not. ok .and. diis%count() >= 2) then
            ! A singular subspace stays singular: its vectors have become
            ! nearly parallel, which is what a slowly converging mode makes of
            ! them, and every later extrapolation would be skipped and leave a
            ! plain fixed-point iteration. Restarting from the newest entry
            ! lets the next two extrapolate again.
            call diis%destroy()
            call diis%init(settings%diis_vectors, 2*n*n, 2*n*n)
            call diis%push(fock, residual)
         end if

         do spin = 1, 2
            call occupy(fock(:, :, spin), first, ncomp, blocks, electrons(:, spin), &
                        density(:, :, spin), error)
         end do
         if (error%has_error()) exit
      end do

      call diis%destroy()
      call xc%destroy()
      if (error%has_error()) then
         call error%add_context("spherical atom Z="//int_to_text(z))
         return
      end if

      atom%d_alpha = density(:, :, 1)
      atom%d_beta = density(:, :, 2)
   end subroutine solve_spherical_atom

   subroutine check_atom(atom_mol, z, error)
      ! What this solver will not take, refused before any integral is built.
      type(czt_molecule_t), intent(in) :: atom_mol
      integer, intent(out) :: z
      type(error_t), intent(inout) :: error

      integer, allocatable :: ang(:), first(:), ncomp(:)
      integer :: n_sub

      z = 0
      if (atom_mol%natm /= 1) then
         call error%set(ERROR_VALIDATION, "spherical atom: expected one centre, got "// &
                        int_to_text(atom_mol%natm))
         return
      end if
      if (.not. allocated(atom_mol%atomic_numbers)) then
         call error%set(ERROR_VALIDATION, "spherical atom: the molecule does not record "// &
                        "its element")
         return
      end if
      z = atom_mol%atomic_numbers(1)
      if (atom_mol%charges(1) <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "spherical atom: a ghost centre has no electrons")
         return
      end if
      if (allocated(atom_mol%core_electrons)) then
         if (any(atom_mol%core_electrons > 0)) then
            call error%set(ERROR_VALIDATION, "spherical atom: not available with an "// &
                           "effective core potential")
            return
         end if
      end if
      if (atom_mol%cartesian) then
         call subshell_layout(atom_mol, ang, first, ncomp, n_sub)
         if (any(ang >= 2)) then
            ! Six Cartesian d functions hold an s-like combination, so averaging
            ! over components does not leave a spherical density.
            call error%set(ERROR_VALIDATION, "spherical atom: a Cartesian basis with "// &
                           "l >= 2 cannot be averaged over its components")
            return
         end if
      end if
   end subroutine check_atom

   pure function wants_functional(settings) result(wants)
      ! Whether the settings name anything other than Hartree-Fock.
      type(spherical_atom_settings_t), intent(in) :: settings
      logical :: wants

      wants = .false.
      if (.not. allocated(settings%functional)) return
      select case (trim(adjustl(settings%functional)))
      case ("", "hf", "HF", "scf", "SCF")
         wants = .false.
      case default
         wants = .true.
      end select
   end function wants_functional

   subroutine build_l_blocks(s, ang, first, ncomp, n_sub, electrons, blocks, error)
      ! Group the contraction columns by l and orthogonalise each group's
      ! averaged overlap, refusing an l whose configuration needs more
      ! orbitals than the basis gives it.
      real(dp), intent(in) :: s(:, :)
      integer, intent(in) :: ang(:), first(:), ncomp(:)
      integer, intent(in) :: n_sub
      real(dp), intent(in) :: electrons(0:L_OCCUPIED_MAX, 2)
      type(l_block_t), intent(inout) :: blocks(0:L_OCCUPIED_MAX)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: s_l(:, :)
      integer :: l, a, p, q, needed

      do l = 0, L_OCCUPIED_MAX
         blocks(l)%n_comp = 2*l + 1
         blocks(l)%column = pack([(a, a=1, n_sub)], ang(1:n_sub) == l)
         blocks(l)%n_col = size(blocks(l)%column)

         needed = orbitals_needed(maxval(electrons(l, :)), blocks(l)%n_comp)
         if (needed == 0) cycle
         if (blocks(l)%n_col == 0) then
            call error%set(ERROR_VALIDATION, "spherical atom: the basis has no l="// &
                           int_to_text(l)//" functions and the configuration occupies "// &
                           int_to_text(needed))
            return
         end if

         allocate (s_l(blocks(l)%n_col, blocks(l)%n_col))
         do q = 1, blocks(l)%n_col
            do p = 1, blocks(l)%n_col
               s_l(p, q) = component_mean(s, first(blocks(l)%column(p)), &
                                          first(blocks(l)%column(q)), ncomp(blocks(l)%column(p)))
            end do
         end do
         call build_orthogonalizer(s_l, blocks(l)%x, blocks(l)%n_mo, error)
         deallocate (s_l)
         if (error%has_error()) return

         if (blocks(l)%n_mo < needed) then
            call error%set(ERROR_VALIDATION, "spherical atom: l="//int_to_text(l)// &
                           " keeps "//int_to_text(blocks(l)%n_mo)//" orbitals and the "// &
                           "configuration occupies "//int_to_text(needed))
            return
         end if
      end do
   end subroutine build_l_blocks

   subroutine occupy(fock, first, ncomp, blocks, electrons, density, error)
      ! Diagonalise each l's averaged Fock block and fill its lowest orbitals
      ! with `electrons(l)` spread over the 2l+1 components: one electron per
      ! component in each full orbital, the remainder evenly in the next.
      real(dp), intent(in) :: fock(:, :)
      integer, intent(in) :: first(:), ncomp(:)
      type(l_block_t), intent(in) :: blocks(0:L_OCCUPIED_MAX)
      real(dp), intent(in) :: electrons(0:L_OCCUPIED_MAX)
      real(dp), intent(out) :: density(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: f_l(:, :), work(:, :), f_ortho(:, :), coeff(:, :)
      real(dp), allocatable :: eigenvalues(:), occupation(:), d_l(:, :)
      integer :: l, p, q, k, m, n_col, n_mo, info, cp, cq

      density = 0.0_dp
      do l = 0, L_OCCUPIED_MAX
         if (electrons(l) <= FRACTION_FLOOR) cycle
         n_col = blocks(l)%n_col
         n_mo = blocks(l)%n_mo

         allocate (f_l(n_col, n_col), work(n_col, n_mo), f_ortho(n_mo, n_mo), &
                   coeff(n_col, n_mo), eigenvalues(n_mo), occupation(n_mo), d_l(n_col, n_col))
         do q = 1, n_col
            do p = 1, n_col
               f_l(p, q) = component_mean(fock, first(blocks(l)%column(p)), &
                                          first(blocks(l)%column(q)), ncomp(blocks(l)%column(p)))
            end do
         end do
         call pic_gemm(f_l, blocks(l)%x, work)                      ! f X
         call pic_gemm(blocks(l)%x, work, f_ortho, transa="T")      ! X^T f X
         call pic_syev(f_ortho, eigenvalues, jobz="V", uplo="U", info=info)
         if (info /= 0) then
            call error%set(ERROR_VALIDATION, "spherical atom: the l="//int_to_text(l)// &
                           " Fock block did not diagonalise")
            return
         end if
         call pic_gemm(blocks(l)%x, f_ortho, coeff)                 ! C = X C'

         call orbital_occupations(electrons(l), blocks(l)%n_comp, occupation)
         d_l = 0.0_dp
         do k = 1, n_mo
            if (occupation(k) <= 0.0_dp) exit
            do q = 1, n_col
               do p = 1, n_col
                  d_l(p, q) = d_l(p, q) + occupation(k)*coeff(p, k)*coeff(q, k)
               end do
            end do
         end do

         ! The same block in every component, and nothing between components.
         do q = 1, n_col
            cq = blocks(l)%column(q)
            do p = 1, n_col
               cp = blocks(l)%column(p)
               do m = 1, ncomp(cp)
                  density(first(cp) + m, first(cq) + m) = d_l(p, q)
               end do
            end do
         end do
         deallocate (f_l, work, f_ortho, coeff, eigenvalues, occupation, d_l)
      end do
   end subroutine occupy

   pure subroutine orbital_occupations(electrons, n_comp, occupation)
      ! One spin's electrons in one l, as an occupation per component of each
      ! orbital in energy order: whole ones, then the fraction, then zeros.
      real(dp), intent(in) :: electrons
      integer, intent(in) :: n_comp
      real(dp), intent(out) :: occupation(:)

      real(dp) :: per_component
      integer :: k

      occupation = 0.0_dp
      per_component = electrons/real(n_comp, dp)
      do k = 1, size(occupation)
         if (per_component <= FRACTION_FLOOR) exit
         occupation(k) = min(1.0_dp, per_component)
         per_component = per_component - occupation(k)
      end do
   end subroutine orbital_occupations

   pure function orbitals_needed(electrons, n_comp) result(needed)
      ! How many orbitals of one l a spin's electrons reach into.
      real(dp), intent(in) :: electrons
      integer, intent(in) :: n_comp
      integer :: needed

      needed = ceiling(electrons/real(n_comp, dp) - FRACTION_FLOOR)
      needed = max(needed, 0)
   end function orbitals_needed

   pure function component_mean(x, first_p, first_q, n_comp) result(mean)
      ! The mean over components of one column pair's diagonal, (1/(2l+1))
      ! sum_m x(p_m, q_m). `first_*` are 0-based, as the layout gives them.
      real(dp), intent(in) :: x(:, :)
      integer, intent(in) :: first_p, first_q, n_comp
      real(dp) :: mean

      integer :: m

      mean = 0.0_dp
      do m = 1, n_comp
         mean = mean + x(first_p + m, first_q + m)
      end do
      mean = mean/real(n_comp, dp)
   end function component_mean

   subroutine commutator_residual(fock, density, s, residual)
      ! F D S - S D F, zero exactly when D is built from F's own eigenvectors.
      real(dp), intent(in) :: fock(:, :), density(:, :), s(:, :)
      real(dp), intent(out) :: residual(:, :)

      real(dp), allocatable :: fd(:, :), fds(:, :)

      allocate (fd(size(fock, 1), size(fock, 2)), fds(size(fock, 1), size(fock, 2)))
      call pic_gemm(fock, density, fd)
      call pic_gemm(fd, s, fds)
      residual = fds - transpose(fds)
   end subroutine commutator_residual

end module mqc_czt_spherical_atom
