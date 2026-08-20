!! CASCI: the CI problem for a chosen active space, in a molecule's orbitals
module mqc_libcint_casci
   !! Everything between a converged set of molecular orbitals and a CI energy.
   !! The CI machinery itself -- `mqc_determinants`, `mqc_ci`, `mqc_davidson` --
   !! knows nothing about molecules and takes a one- and two-electron integral
   !! set over some orbitals. This is what produces that set.
   !!
   !! The orbitals are partitioned into three groups. **Inactive** orbitals are
   !! doubly occupied in every determinant, so they never appear in the CI at
   !! all; what they contribute is a fixed energy and a mean field that the
   !! active electrons move in. **Active** orbitals are the ones the CI
   !! distributes electrons over. **Virtual** orbitals are empty in every
   !! determinant and contribute nothing whatsoever -- a CASCI energy does not
   !! depend on them, which is worth knowing because it means the size of the
   !! basis set past the active space affects this calculation only through the
   !! orbitals themselves.
   !!
   !! So the active-space Hamiltonian is
   !!
   !!     h_tu^eff = h_tu + sum_i [2 (tu|ii) - (ti|iu)]
   !!
   !! -- the bare one-electron integral plus the inactive mean field, which is
   !! exactly the closed-shell Fock matrix built from the inactive density and
   !! then projected onto the active orbitals. That is why no separate
   !! Coulomb-exchange assembly appears below: `build_fock_direct` already
   !! computes `H + J - K/2`, and fed an inactive density it is the inactive
   !! Fock.
   !!
   !! The constant is the energy of the inactive electrons in their own field
   !! plus nuclear repulsion, and it is not a detail -- it is most of the total
   !! energy. Only the difference from it is variational here.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_io, only: to_char
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_libcint_direct, only: schwarz_bounds, build_fock_direct, direct_stats_t
   use mqc_libcint_mp2, only: transform_block
   use mqc_determinants, only: link_table_t, build_link_table, n_strings
   use mqc_ci, only: absorb_one_electron, ci_diagonal
   use mqc_rdm, only: active_space_rdms
   use mqc_ormas_space, only: ormas_space_t, build_ormas_space
   use mqc_ormas_ci, only: ormas_solve, ormas_density_matrices
   use mqc_davidson, only: davidson_lowest, davidson_result_t
   implicit none
   private

   public :: active_space_integrals
   public :: run_libcint_casci
   public :: run_libcint_ormas_ci
   public :: casci_result_t

   type :: casci_result_t
      !! What a CASCI leaves behind
      real(dp) :: energy = 0.0_dp          !! Total, including the inactive constant
      real(dp) :: core_energy = 0.0_dp     !! Inactive electrons plus nuclear repulsion
      real(dp) :: active_energy = 0.0_dp   !! The CI eigenvalue alone
      real(dp), allocatable :: energies(:)     !! (n_roots) totals, ascending
      real(dp), allocatable :: ci_vector(:, :)
         !! (n_alpha_strings, n_beta_strings), the ground root. Left unallocated
         !! by a restricted space, whose determinants are not a rectangle --
         !! `ci_flat` carries it there instead.
      real(dp), allocatable :: ci_flat(:)
         !! (n_determinants), the ground root of a restricted space
      type(ormas_space_t) :: ormas
         !! The partition `ci_flat` is addressed by, kept because a flat list of
         !! coefficients means nothing without it. Populated only by a
         !! restricted solve; `ormas%n_determinants` is zero otherwise, which is
         !! how a caller tells.
      real(dp), allocatable :: vectors(:, :, :)   !! All roots
      real(dp), allocatable :: dm1(:, :)
         !! (n_active, n_active) one-particle density of the ground root, in the
         !! active orbital basis. Built here rather than left to the caller
         !! because doing it outside means rebuilding the excitation tables,
         !! which this routine has already built and thrown away.
      real(dp), allocatable :: dm2(:, :, :, :)
         !! Active two-particle density, spin-traced, in the convention of
         !! `mqc_rdm`. Kept for the same reason and by the same argument, and
         !! because an energy decomposition of a correlated wave function cannot
         !! be rebuilt from `dm1`: the single-determinant
         !! `Gamma = gamma gamma - (1/2) gamma gamma` is exactly what
         !! correlation makes false, and the difference between the two is the
         !! whole correlation contribution.
      integer :: n_determinants = 0
      integer :: iterations = 0
      integer :: sigma_products = 0
      logical :: converged = .false.
   end type casci_result_t

contains

   subroutine active_space_integrals(mol, orbitals, n_inactive, n_active, &
                                     h_eff, eri_act, core_energy, error)
      !! The Hamiltonian the active electrons see
      !!
      !! Returns the effective one-electron integrals over active orbitals, the
      !! two-electron integrals over active orbitals, and the constant.
      !!
      !! Costs one packed AO integral build, which is `n_ao^4 / 4` and is the
      !! memory ceiling on this routine rather than anything about the active
      !! space. An active-space transformation driven directly from AO integrals
      !! would not need it and is the obvious thing to do when it starts to
      !! matter; nothing above here would change.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)      !! (n_ao, n_mo), the SCF's
      integer, intent(in) :: n_inactive           !! Doubly occupied, frozen
      integer, intent(in) :: n_active
      real(dp), allocatable, intent(out) :: h_eff(:, :)        !! (n_active, n_active)
      real(dp), allocatable, intent(out) :: eri_act(:, :, :, :)
      real(dp), intent(out) :: core_energy
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: h_ao(:, :), density(:, :), fock(:, :), bounds(:, :)
      real(dp), allocatable :: eri_packed(:, :), work(:, :)
      real(dp), allocatable :: c_inactive(:, :), c_active(:, :)
      type(direct_stats_t) :: stats
      integer :: n_ao, n_mo, p, q

      if (error%has_error()) return
      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      core_energy = 0.0_dp

      if (n_inactive < 0 .or. n_active < 0) then
         call error%set(ERROR_VALIDATION, "an active space cannot have "// &
                        to_char(n_inactive)//" inactive and "//to_char(n_active)// &
                        " active orbitals.")
         return
      end if
      if (n_inactive + n_active > n_mo) then
         call error%set(ERROR_VALIDATION, to_char(n_inactive)//" inactive plus "// &
                        to_char(n_active)//" active orbitals is more than the "// &
                        to_char(n_mo)//" the basis provides.")
         return
      end if

      allocate (c_active(n_ao, n_active))
      c_active = orbitals(:, n_inactive + 1:n_inactive + n_active)

      call mol%core_hamiltonian(h_ao)

      ! The inactive mean field. With no inactive orbitals the density is zero,
      ! the Fock build returns the bare core Hamiltonian, and the constant is
      ! nuclear repulsion alone -- which is the full-CI case and needs no
      ! special path.
      allocate (density(n_ao, n_ao), fock(n_ao, n_ao))
      density = 0.0_dp
      if (n_inactive > 0) then
         allocate (c_inactive(n_ao, n_inactive))
         c_inactive = orbitals(:, 1:n_inactive)
         call pic_gemm(c_inactive, c_inactive, density, transb="T", &
                       alpha=2.0_dp, beta=0.0_dp)
         deallocate (c_inactive)
      end if

      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) return
      call build_fock_direct(mol, h_ao, density, bounds, fock, stats, error)
      if (error%has_error()) return

      ! E_core = E_nuc + (1/2) sum_pq D_pq (h_pq + F_pq). The half is because
      ! the two-electron part of F double counts the interaction it describes.
      core_energy = mol%nuclear_repulsion()
      do q = 1, n_ao
         do p = 1, n_ao
            core_energy = core_energy + 0.5_dp*density(p, q)*(h_ao(p, q) + fock(p, q))
         end do
      end do

      allocate (h_eff(n_active, n_active), work(n_ao, n_active))
      call pic_gemm(fock, c_active, work)
      call pic_gemm(c_active, work, h_eff, transa="T")
      deallocate (work)

      call mol%eris_packed(eri_packed)
      call transform_block(eri_packed, c_active, c_active, c_active, c_active, eri_act)

      deallocate (h_ao, density, fock, bounds, eri_packed, c_active)
   end subroutine active_space_integrals

   subroutine run_libcint_casci(mol, orbitals, n_inactive, n_active, &
                                n_alpha, n_beta, result, error, n_roots, verbose, &
                                tolerance, guess)
      !! A complete-active-space CI on converged orbitals
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active
      integer, intent(in) :: n_alpha, n_beta      !! Active electrons, per spin
      type(casci_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_roots
      logical, intent(in), optional :: verbose
      real(dp), intent(in), optional :: tolerance
      real(dp), intent(in), optional :: guess(:, :, :)

      real(dp), allocatable :: h_eff(:, :), eri_act(:, :, :, :)
      real(dp), allocatable :: folded(:, :), diagonal(:, :)
      type(link_table_t) :: alpha, beta
      type(davidson_result_t) :: davidson
      character(len=128) :: line
      integer :: roots, i
      logical :: loud

      if (error%has_error()) return
      roots = 1
      if (present(n_roots)) roots = n_roots
      loud = .false.
      if (present(verbose)) loud = verbose

      if (n_alpha < 0 .or. n_beta < 0 .or. n_alpha > n_active .or. n_beta > n_active) then
         call error%set(ERROR_VALIDATION, "an active space of "//to_char(n_active)// &
                        " orbitals cannot hold "//to_char(n_alpha)//" alpha and "// &
                        to_char(n_beta)//" beta electrons.")
         return
      end if

      call active_space_integrals(mol, orbitals, n_inactive, n_active, h_eff, &
                                  eri_act, result%core_energy, error)
      if (error%has_error()) return

      call build_link_table(n_active, n_alpha, alpha, error)
      call build_link_table(n_active, n_beta, beta, error)
      if (error%has_error()) return
      result%n_determinants = alpha%n_strings*beta%n_strings

      if (n_alpha + n_beta == 0) then
         ! No active electrons: the CI is one determinant with no electrons in
         ! it and the energy is the constant. Worth handling rather than
         ! dividing by zero in the folding.
         result%active_energy = 0.0_dp
         result%energy = result%core_energy
         result%converged = .true.
         allocate (result%energies(1))
         result%energies(1) = result%energy
         call alpha%destroy()
         call beta%destroy()
         return
      end if

      call absorb_one_electron(h_eff, eri_act, n_alpha + n_beta, folded, error)
      call ci_diagonal(h_eff, eri_act, alpha, beta, diagonal, error)
      if (error%has_error()) return

      call davidson_lowest(folded, diagonal, alpha, beta, roots, davidson, error, &
                           tolerance=tolerance, guess=guess)
      if (error%has_error()) return

      result%converged = davidson%converged
      result%iterations = davidson%iterations
      result%sigma_products = davidson%sigma_products
      result%active_energy = davidson%values(1)
      result%energy = result%core_energy + davidson%values(1)
      allocate (result%energies(roots))
      do i = 1, roots
         result%energies(i) = result%core_energy + davidson%values(i)
      end do
      call move_alloc(davidson%vectors, result%vectors)
      allocate (result%ci_vector(alpha%n_strings, beta%n_strings))
      result%ci_vector = result%vectors(:, :, 1)
      call active_space_rdms(result%ci_vector, alpha, beta, result%dm1, result%dm2, error)
      if (error%has_error()) return

      if (loud) then
         call logger%info("")
         call logger%info("  complete active space CI")
         write (line, "(a,i0,a,i0,a)") "    active space                CAS(", &
            n_alpha + n_beta, ",", n_active, ")"
         call logger%info(trim(line))
         write (line, "(a,i0)") "    inactive orbitals           ", n_inactive
         call logger%info(trim(line))
         write (line, "(a,i0)") "    determinants                ", result%n_determinants
         call logger%info(trim(line))
         write (line, "(a,i0,a,i0,a)") "    iterations                  ", &
            davidson%iterations, " (", davidson%sigma_products, " sigma products)"
         call logger%info(trim(line))
         write (line, "(a,f22.10)") "    inactive plus nuclear       ", result%core_energy
         call logger%info(trim(line))
         write (line, "(a,f22.10)") "    active                      ", &
            result%active_energy
         call logger%info(trim(line))
         write (line, "(a,f22.10)") "    total                       ", result%energy
         call logger%info(trim(line))
         if (.not. result%converged) then
            call logger%warning("    the CI did not converge to the requested residual")
         end if
      end if

      call alpha%destroy()
      call beta%destroy()
      deallocate (h_eff, eri_act, folded, diagonal)
   end subroutine run_libcint_casci

   subroutine run_libcint_ormas_ci(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                   subspaces, min_electrons, max_electrons, result, &
                                   error, n_roots, verbose, tolerance, guess)
      !! A CI over an occupation-restricted active space, on converged orbitals
      !!
      !! The same integrals as a CASCI -- a restricted space changes which
      !! determinants are kept, not what the Hamiltonian over them is -- and
      !! then `ormas_solve` instead of a Davidson over the rectangle.
      !!
      !! Fixed orbitals only. Optimising them for a restricted space is not the
      !! same problem as for a complete one: rotating one active orbital into
      !! another stops being redundant the moment the space is not complete, so
      !! the orbital gradient acquires a block a CASSCF is entitled to ignore.
      !! Pretending otherwise would converge to the wrong answer quietly, so the
      !! caller is refused instead.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active
      integer, intent(in) :: n_alpha, n_beta
      integer, intent(in) :: subspaces(:), min_electrons(:), max_electrons(:)
      type(casci_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_roots
      logical, intent(in), optional :: verbose
      real(dp), intent(in), optional :: tolerance
      real(dp), intent(in), optional :: guess(:, :)   !! (n_determinants, n_roots)

      real(dp), allocatable :: h_eff(:, :), eri_act(:, :, :, :)
      real(dp), allocatable :: energies(:), vectors(:, :)
      type(ormas_space_t) :: space
      character(len=128) :: line
      integer :: roots, i
      logical :: loud

      if (error%has_error()) return
      roots = 1
      if (present(n_roots)) roots = n_roots
      loud = .false.
      if (present(verbose)) loud = verbose

      call build_ormas_space(subspaces, n_active, n_alpha, n_beta, min_electrons, &
                             max_electrons, space, error)
      if (error%has_error()) return

      call active_space_integrals(mol, orbitals, n_inactive, n_active, h_eff, &
                                  eri_act, result%core_energy, error)
      if (error%has_error()) return

      result%n_determinants = int(space%n_determinants)
      if (n_alpha + n_beta == 0) then
         result%active_energy = 0.0_dp
         result%energy = result%core_energy
         result%converged = .true.
         allocate (result%energies(1))
         result%energies(1) = result%energy
         call space%destroy()
         return
      end if

      if (present(guess)) then
         call ormas_solve(space, h_eff, eri_act, roots, energies, vectors, error, &
                          tolerance=tolerance, guess=guess)
      else
         call ormas_solve(space, h_eff, eri_act, roots, energies, vectors, error, &
                          tolerance=tolerance)
      end if
      if (error%has_error()) return

      result%converged = .true.
      result%active_energy = energies(1)
      result%energy = result%core_energy + energies(1)
      allocate (result%energies(roots))
      do i = 1, roots
         result%energies(i) = result%core_energy + energies(i)
      end do
      result%ci_flat = vectors(:, 1)
      ! Deep-copied rather than referenced, then destroyed below as before: the
      ! flat coefficients are unreadable without the addressing that produced
      ! them, and anything wanting to re-express this wave function in another
      ! orbital basis needs both.
      result%ormas = space

      call ormas_density_matrices(space, result%ci_flat, result%dm1, result%dm2, error)
      if (error%has_error()) return

      if (loud) then
         call logger%info("")
         call logger%info("  occupation-restricted active space CI")
         write (line, "(a,i0,a,i0,a)") "    active space                ORMAS(", &
            n_alpha + n_beta, ",", n_active, ")"
         call logger%info(trim(line))
         write (line, "(a,i0)") "    inactive orbitals           ", n_inactive
         call logger%info(trim(line))
         write (line, "(a,i0)") "    subspaces                   ", size(subspaces)
         call logger%info(trim(line))
         do i = 1, size(subspaces)
            write (line, "(a,i0,a,i0,a,i0,a,i0)") "      space ", i, ": orbitals from ", &
               subspaces(i), ", electrons ", space%min_electrons(i), " to ", &
               space%max_electrons(i)
            call logger%info(trim(line))
         end do
         write (line, "(a,i0)") "    determinants                ", result%n_determinants
         call logger%info(trim(line))
         write (line, "(a,f22.10)") "    inactive plus nuclear       ", result%core_energy
         call logger%info(trim(line))
         write (line, "(a,f22.10)") "    active                      ", &
            result%active_energy
         call logger%info(trim(line))
         write (line, "(a,f22.10)") "    total                       ", result%energy
         call logger%info(trim(line))
      end if

      call space%destroy()
      deallocate (h_eff, eri_act)
   end subroutine run_libcint_ormas_ci

end module mqc_libcint_casci
